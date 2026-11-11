#!/usr/bin/env Rscript
# Purpose: anchor plan A8: in 3 seeded batches no control is a reference; compare heb_kv_bap (abstains there)
#          with heb_kv_bap_anc (internal anchors) for one leave-one-donor-out fold.
# Run:     Rscript experiments/benchmarks/06_cytether-cytobatchflagr/02c_correct_strip.R <prep_dir> <out_dir> <method> <held_out>
#          method: heb_kv_bap, heb_kv_bap_anc or heb_kv_bap_anc_u (A14)
# Inputs:  <prep_dir>/cells.rds, markers.csv, populations.rds (01_prep.R)
# Outputs (in <out_dir>/<method>/fold_<held_out>/, which must not exist):
#   corrected.rds   numeric matrix, same rows and marker columns as cells.rds
#   effects.csv.gz  per batch, population and marker effect, se_donor, df, support
#   anchors.csv.gz  anchor modes (anchor arms only) and anchor_markers.csv (separation, eligible)
#   status.csv      method, fold, stripped batches, abstained, runtime_s

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
source("experiments/benchmarks/06_cytether-cytobatchflagr/correct_lib.R")

# Parameters (approved, anchor plan) ------------------------------------------------
N_STRIP <- 3 # A8
STRIP_SEED <- 1 # A8
ANCHOR_POSITIVE <- "BUV563-A" # A1: CD3 channel (markers.csv)
ANCHOR_MIN_SEPARATION <- 3 # A2
ANCHOR_SEED <- 1 # A7

args <- commandArgs(trailingOnly = TRUE)
prep_dir <- args[[1]]
method <- args[[3]]
held_out <- args[[4]]
stopifnot("method is heb_kv_bap, heb_kv_bap_anc or heb_kv_bap_anc_u" = method %in% c("heb_kv_bap", "heb_kv_bap_anc", "heb_kv_bap_anc_u"))
out_dir <- file.path(args[[2]], method, paste0("fold_", held_out))
if (dir.exists(out_dir)) stop(sprintf("%s exists; use a new out_dir.", out_dir), call. = FALSE)

# Inputs ---------------------------------------------------------------------------
cells <- readRDS(file.path(prep_dir, "cells.rds"))
markers <- read_csv(file.path(prep_dir, "markers.csv"), show_col_types = FALSE)$channel
population <- readRDS(file.path(prep_dir, "populations.rds"))[[paste0("fold_", held_out)]]
stopifnot("populations align with cells" = length(population) == nrow(cells))
x <- as.matrix(cells[, markers])

# Stripped batches: same draw for every fold and method --------------------------------
set.seed(STRIP_SEED)
strip <- sort(sample(setdiff(sort(unique(cells$batch)), ANCHOR), N_STRIP))
meta <- cells |>
  select(sample_id, batch, control) |>
  mutate(
    donor_id = control,
    is_holdout_reference = control == held_out,
    is_reference = control != held_out & !batch %in% strip,
    arrival = "initial"
  )

# Anchors from every non-held-out sample, stripped batches included --------------------
anchors <- NULL
if (method %in% c("heb_kv_bap_anc", "heb_kv_bap_anc_u")) {
  rows <- !meta$is_holdout_reference
  anchors <- sw_corrector_anchors(x[rows, ], meta$sample_id[rows], meta$batch[rows], markers,
                                  positive_markers = ANCHOR_POSITIVE, min_separation = ANCHOR_MIN_SEPARATION,
                                  seed = ANCHOR_SEED)
}

# Correct ----------------------------------------------------------------------------
t0 <- Sys.time()
res <- correct_t6(x, meta, population, method, markers, anchors = anchors)
runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
stopifnot("corrected matrix keeps the shape" = identical(dim(res$x), dim(x)))

# Write --------------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE)
saveRDS(res$x, file.path(out_dir, "corrected.rds"))
write_csv(res$effects, file.path(out_dir, "effects.csv.gz"))
if (!is.null(anchors)) {
  write_csv(anchors$modes, file.path(out_dir, "anchors.csv.gz"))
  write_csv(anchors$markers, file.path(out_dir, "anchor_markers.csv"))
}
write_csv(
  tibble(method = method, fold = held_out, supported = res$supported, abstained = res$abstained,
         reason = res$reason, fits_on_heldout = FALSE, stripped = paste(strip, collapse = ";"), runtime_s = runtime),
  file.path(out_dir, "status.csv")
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
