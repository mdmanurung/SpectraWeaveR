#!/usr/bin/env Rscript
# Purpose: correct the CytoBatchFlagR cells with one method for one leave-one-donor-out fold.
# Run:     pixi run Rscript experiments/benchmarks/06_cytether-cytobatchflagr/02_correct.R <prep_dir> <out_dir> <method> <held_out>
#          method: none, paired_shift, heb, cycombine, cytonorm, cytofin; held_out: tcu1, tcu2, tcs3, tcs4
# Inputs:  <prep_dir>/cells.rds, markers.csv, populations.rds (from 01_prep.R)
# Outputs (in <out_dir>/<method>/fold_<held_out>/, which must not exist):
#   corrected.rds   numeric matrix, same rows and marker columns as cells.rds
#   effects.csv.gz  cyTether engines only: per batch, population and marker effect, se_donor, df
#   status.csv      one row: method, fold, supported, abstained, reason, fits_on_heldout, runtime_s

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
source("experiments/benchmarks/05_cytether-synthetic/adapters.R")
source("experiments/benchmarks/_common/adapter_cytofin.R")

# Parameters (approved, plan T6) -------------------------------------------------
ANCHOR <- "batch01" # Q8
MIN_CELLS <- 10 # G1
VARIANCE_FLOOR <- 1e-8 # heb numerical floor, as in the synthetic benchmark
HEB_POOL <- "population" # N4: lead variant
CYTONORM_SEED <- 1 # N6 for T6; overrides the synthetic adapter's 431
CYTOFIN_MODE <- "meanshift" # N7
MIN_BATCHES <- 2
MIN_BATCHES_ORIGIN <- "two batches are the least that define a shift"

args <- commandArgs(trailingOnly = TRUE)
prep_dir <- args[[1]]
method <- args[[3]]
held_out <- args[[4]]
stopifnot("unknown method" = method %in% c("none", "paired_shift", "heb", "cycombine", "cytonorm", "cytofin"))
out_dir <- file.path(args[[2]], method, paste0("fold_", held_out))
if (dir.exists(out_dir)) stop(sprintf("%s exists; use a new out_dir.", out_dir), call. = FALSE)

# Inputs ---------------------------------------------------------------------------
cells <- readRDS(file.path(prep_dir, "cells.rds"))
markers <- read_csv(file.path(prep_dir, "markers.csv"), show_col_types = FALSE)$channel
population <- readRDS(file.path(prep_dir, "populations.rds"))[[paste0("fold_", held_out)]]
stopifnot("populations align with cells" = length(population) == nrow(cells))
x <- as.matrix(cells[, markers])

# Roles: the held-out control never calibrates; the other three are reference donors in every batch.
meta <- cells |>
  select(sample_id, batch, control) |>
  mutate(
    donor_id = control,
    is_reference = control != held_out,
    is_holdout_reference = control == held_out,
    arrival = "initial"
  )

# Correct ----------------------------------------------------------------------------
t0 <- Sys.time()
effects <- NULL
res <- switch(method,
  none = list(x = x, supported = TRUE, abstained = 0L, reason = "raw", fits_on_heldout = FALSE),
  cycombine = run_cycombine(x, meta, markers, ANCHOR),
  cytonorm = run_cytonorm(x, meta, markers),
  cytofin = run_cytofin(x, meta, markers, ANCHOR, mode = CYTOFIN_MODE),
  {
    # cyTether engines: fit on reference cells, integrate every cell with abstention.
    samples <- meta |>
      distinct(sample_id, batch, donor_id, is_reference) |>
      mutate(
        reference_role = if_else(is_reference, "technical_replicate", "none"),
        reference_id = if_else(is_reference, donor_id, NA_character_),
        aliquot = if_else(is_reference, sample_id, NA_character_)
      )
    manifest <- sw_sample_manifest(
      samples, sample_id = "sample_id", batch = "batch", donor = "donor_id", aliquot = "aliquot",
      reference_role = "reference_role", reference_id = "reference_id"
    )
    audit <- sw_design_audit(
      samples, minimum_batches = MIN_BATCHES, minimum_batches_origin = MIN_BATCHES_ORIGIN,
      require_shared_reference = TRUE, anchor = ANCHOR
    )
    # N4: flat hierarchy, one random effect per batch.
    hierarchy <- if (method == "heb") distinct(meta, batch) |> mutate(study = batch) else NULL
    xt <- as_tibble(x) |> mutate(sample_id = meta$sample_id, .before = 1)
    fit_rows <- meta$is_reference
    model <- sw_corrector_fit(
      xt[fit_rows, ], manifest, audit, population = population[fit_rows], markers = markers,
      transform = list(name = "asinh", cofactor = 6000), anchor = ANCHOR, engine = method,
      min_cells = MIN_CELLS, seed = 1, hierarchy = hierarchy, variance_floor = VARIANCE_FLOOR,
      pool = if (method == "heb") HEB_POOL else NULL,
      anchor_effect = if (method == "heb") "random" else NULL
    )
    out <- sw_corrector_integrate(model, xt, manifest, population = population,
                                  transform = list(name = "asinh", cofactor = 6000), unsupported = "abstain")
    effects <- out$model$effects
    list(x = as.matrix(out$x[, markers]), supported = TRUE,
         abstained = sum(out$actions$action == "abstain"), reason = paste("engine", method),
         fits_on_heldout = FALSE)
  }
)
runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
stopifnot("corrected matrix keeps the shape" = identical(dim(res$x), dim(x)))

# Write --------------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE)
saveRDS(res$x, file.path(out_dir, "corrected.rds"))
if (!is.null(effects)) write_csv(effects, file.path(out_dir, "effects.csv.gz"))
write_csv(
  tibble(method = method, fold = held_out, supported = res$supported, abstained = res$abstained,
         reason = res$reason, fits_on_heldout = isTRUE(res$fits_on_heldout), runtime_s = runtime),
  file.path(out_dir, "status.csv")
)
