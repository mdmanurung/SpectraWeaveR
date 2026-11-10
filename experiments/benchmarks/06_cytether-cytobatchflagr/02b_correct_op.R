#!/usr/bin/env Rscript
# Purpose: OpenProblems track of T6: correct one replicate split of one 2-donor fold with one method.
# Run:     pixi run Rscript experiments/benchmarks/06_cytether-cytobatchflagr/02b_correct_op.R <prep_dir> <out_dir> <method> <fold> <split>
#          fold: op_tcu1_tcs3, op_tcu1_tcs4, op_tcu2_tcs3, op_tcu2_tcs4; split: 1 or 2
# Inputs:  <prep_dir>/cells.rds, markers.csv, populations_op.rds, split_op.rds (01_prep.R, 01b_prep_op.R)
# Outputs (in <out_dir>/<method>/<fold>/split_<k>/, which must not exist):
#   corrected.f64   float64, column-major (n_rows x 20) for the rows in rows.i32
#   rows.i32        int32, 1-based row indices into cells.rds: the two reference donors plus split k of the evaluated pair
#   status.csv      method, fold, split, supported, abstained, reason, fits_on_heldout, runtime_s

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
source("experiments/benchmarks/06_cytether-cytobatchflagr/correct_lib.R")

args <- commandArgs(trailingOnly = TRUE)
prep_dir <- args[[1]]
method <- args[[3]]
fold <- args[[4]]
split_k <- as.integer(args[[5]])
stopifnot("unknown method" = method %in% c("none", "paired_shift", "heb", "cycombine", "cytonorm", "cytofin"))
out_dir <- file.path(args[[2]], method, fold, paste0("split_", split_k))
if (dir.exists(out_dir)) stop(sprintf("%s exists; use a new out_dir.", out_dir), call. = FALSE)

# Inputs ---------------------------------------------------------------------------
cells <- readRDS(file.path(prep_dir, "cells.rds"))
markers <- read_csv(file.path(prep_dir, "markers.csv"), show_col_types = FALSE)$channel
population <- readRDS(file.path(prep_dir, "populations_op.rds"))[[fold]]
split <- readRDS(file.path(prep_dir, "split_op.rds"))
stopifnot("populations align with cells" = length(population) == nrow(cells))
evaluated <- strsplit(sub("^op_", "", fold), "_")[[1]]

# The task's integrated split k: control (reference) samples plus split k of the evaluated samples.
is_ref <- !cells$control %in% evaluated
rows <- which(is_ref | split == split_k)
meta <- cells[rows, ] |>
  select(sample_id, batch, control) |>
  mutate(donor_id = control, is_reference = is_ref[rows], is_holdout_reference = !is_ref[rows], arrival = "initial")

# Correct ----------------------------------------------------------------------------
t0 <- Sys.time()
res <- correct_t6(as.matrix(cells[rows, markers]), meta, population[rows], method, markers)
runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
stopifnot("corrected matrix keeps the shape" = identical(dim(res$x), c(length(rows), length(markers))))

# Write --------------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE)
writeBin(as.vector(res$x), file.path(out_dir, "corrected.f64"), size = 8)
writeBin(as.integer(rows), file.path(out_dir, "rows.i32"), size = 4)
write_csv(
  tibble(method = method, fold = fold, split = split_k, supported = res$supported, abstained = res$abstained,
         reason = res$reason, fits_on_heldout = isTRUE(res$fits_on_heldout), runtime_s = runtime),
  file.path(out_dir, "status.csv")
)
