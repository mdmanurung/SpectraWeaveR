#!/usr/bin/env Rscript
# Purpose: OpenProblems track of T6 (O1, O3, O5): frozen FlowSOM populations per 2-donor fold,
#          the seeded replicate split, and a binary export of the raw cells for the Python scorer.
# Run:     pixi run Rscript experiments/benchmarks/06_cytether-cytobatchflagr/01b_prep_op.R <prep_dir>
# Inputs:  <prep_dir>/cells.rds, markers.csv (from 01_prep.R)
# Outputs (new files in <prep_dir>; the script stops if any exists):
#   populations_op.rds  one integer column per OpenProblems fold (op_tcu1_tcs3 ...): FlowSOM metacluster
#   split_op.rds        integer per cell, 1 or 2: replicate split of the task (used for evaluated cells)
#   cells_raw.f64       the 20 marker columns as float64, column-major (n_cells x 20)
#   obs.csv             sample_id, batch, control, group (stimulated / unstimulated) per cell
#   populations_op.i32, split_op.i32   int32 copies of the two .rds files for the Python scorer

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(FlowSOM)
})

# Parameters (approved, plan T6) -------------------------------------------------
ANCHOR <- "batch01" # Q8
FLOWSOM_XDIM <- 10 # N5
FLOWSOM_YDIM <- 10 # N5
FLOWSOM_NCLUS <- 20 # N5
FLOWSOM_SEED <- 1 # N6
SPLIT_SEED <- 1 # O3
UNSTIM <- c("tcu1", "tcu2")
STIM <- c("tcs3", "tcs4")

prep_dir <- commandArgs(trailingOnly = TRUE)[[1]]
outs <- file.path(prep_dir, c("populations_op.rds", "split_op.rds", "cells_raw.f64", "obs.csv"))
if (any(file.exists(outs))) stop("OpenProblems prep outputs exist; use a new prep dir.", call. = FALSE)
cells <- readRDS(file.path(prep_dir, "cells.rds"))
markers <- read_csv(file.path(prep_dir, "markers.csv"), show_col_types = FALSE)$channel
xm <- as.matrix(cells[, markers])

# Populations per 2-donor fold --------------------------------------------------------
# O1: hold out one unstimulated and one stimulated control; FlowSOM on the anchor's other two.
folds <- tidyr::expand_grid(u = UNSTIM, s = STIM) |> mutate(name = paste0("op_", u, "_", s))
fit_fold <- function(u, s) {
  anchor_ref <- cells$batch == ANCHOR & !cells$control %in% c(u, s)
  fsom <- FlowSOM(xm[anchor_ref, ], colsToUse = markers, xdim = FLOWSOM_XDIM, ydim = FLOWSOM_YDIM,
                  nClus = FLOWSOM_NCLUS, seed = FLOWSOM_SEED, silent = TRUE)
  as.integer(GetMetaclusters(NewData(fsom, xm)))
}
populations_op <- map2(folds$u, folds$s, fit_fold) |>
  set_names(folds$name) |>
  as_tibble()
saveRDS(populations_op, outs[[1]])

# Replicate split ------------------------------------------------------------------------
# O3: seeded random 50/50 split of every cell; only evaluated (held-out) cells use it.
set.seed(SPLIT_SEED)
saveRDS(sample(rep(1:2, length.out = nrow(cells))), outs[[2]])

# Export for Python ------------------------------------------------------------------------
writeBin(as.vector(xm), outs[[3]], size = 8)
cells |>
  select(sample_id, batch, control) |>
  mutate(group = if_else(control %in% STIM, "stimulated", "unstimulated")) |>
  write_csv(outs[[4]])
# int32 copies of the populations (fold-major, in column order) and the split, for 04_score_op.py.
writeBin(unlist(populations_op, use.names = FALSE), file.path(prep_dir, "populations_op.i32"), size = 4)
writeBin(readRDS(outs[[2]]), file.path(prep_dir, "split_op.i32"), size = 4)
