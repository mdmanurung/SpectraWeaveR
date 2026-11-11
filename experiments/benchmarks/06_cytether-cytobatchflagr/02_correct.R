#!/usr/bin/env Rscript
# Purpose: correct the CytoBatchFlagR cells with one method for one leave-one-donor-out fold.
# Run:     pixi run Rscript experiments/benchmarks/06_cytether-cytobatchflagr/02_correct.R <prep_dir> <out_dir> <method> <held_out>
#          method: one of T6_METHODS in correct_lib.R; held_out: tcu1, tcu2, tcs3, tcs4
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
source("experiments/benchmarks/06_cytether-cytobatchflagr/correct_lib.R") # parameters and methods

args <- commandArgs(trailingOnly = TRUE)
prep_dir <- args[[1]]
method <- args[[3]]
held_out <- args[[4]]
stopifnot("unknown method" = method %in% T6_METHODS)
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
res <- correct_t6(x, meta, population, method, markers)
effects <- res$effects
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
