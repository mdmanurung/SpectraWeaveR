#!/usr/bin/env Rscript
# Purpose: score the held-out donor with the CytoNorm 2020 EMD (Q1, D1) for every method and fold.
# Run:     pixi run Rscript experiments/benchmarks/06_cytether-cytobatchflagr/03_score_emd.R <prep_dir> <results_dir> <out_dir>
# Inputs:  <prep_dir>/cells.rds, markers.csv, populations.rds; <results_dir>/<method>/fold_<c>/corrected.rds
# Outputs (in <out_dir>, which must not exist):
#   emd_long.csv     method, fold, population, marker, emd (max over batch pairs of the held-out donor's files)
#   emd_summary.csv  method, fold, median_emd over population x marker, n_pairs scored
#   status_all.csv   every status.csv stacked
#   sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
})

# Parameters (CytoNorm 2020 emdEvaluation defaults; D1) -------------------------------
BIN_SIZE <- 0.1
MIN_RANGE <- -100
MAX_RANGE <- 100
MANUAL_THRESHOLD <- 2 # a file needs this many cells of a population to enter the comparison

args <- commandArgs(trailingOnly = TRUE)
prep_dir <- args[[1]]
results_dir <- args[[2]]
out_dir <- args[[3]]
if (dir.exists(out_dir)) stop(sprintf("%s exists; use a new out_dir.", out_dir), call. = FALSE)

cells <- readRDS(file.path(prep_dir, "cells.rds"))
markers <- read_csv(file.path(prep_dir, "markers.csv"), show_col_types = FALSE)
populations <- readRDS(file.path(prep_dir, "populations.rds"))
breaks <- seq(MIN_RANGE, MAX_RANGE, by = BIN_SIZE)

# EMD between 1D histograms in bin units equals sum |CDF1 - CDF2|; this is what emdist::emd2d
# returns for one-column histograms. Max over all batch pairs, NA files skipped.
max_pair_emd <- function(values, batch) {
  cdf <- map(split(values, batch), function(v) {
    if (length(v) < MANUAL_THRESHOLD) return(NULL)
    cumsum(hist(v, breaks = breaks, plot = FALSE)$counts) / length(v)
  }) |>
    compact()
  if (length(cdf) < 2) return(NA_real_)
  m <- do.call(cbind, cdf)
  pairs <- combn(ncol(m), 2)
  max(vapply(seq_len(ncol(pairs)), \(k) sum(abs(m[, pairs[1, k]] - m[, pairs[2, k]])), numeric(1)))
}

# Check the shortcut against emdist once.
set.seed(1)
a <- hist(rnorm(500), breaks = breaks, plot = FALSE)$counts / 500
b <- hist(rnorm(500, 0.7), breaks = breaks, plot = FALSE)$counts / 500
stopifnot("CDF shortcut matches emdist::emd2d" = abs(sum(abs(cumsum(a) - cumsum(b))) - emdist::emd2d(matrix(a), matrix(b))) < 1e-6)

# Score -----------------------------------------------------------------------------------
runs <- list.files(results_dir, pattern = "corrected.rds", recursive = TRUE, full.names = TRUE)
score_run <- function(path) {
  fold <- sub("fold_", "", basename(dirname(path)))
  method <- basename(dirname(dirname(path)))
  x <- readRDS(path)
  held <- cells$control == fold
  pop <- populations[[paste0("fold_", fold)]][held]
  batch <- cells$batch[held]
  grid <- tidyr::expand_grid(population = sort(unique(pop)), marker = markers$channel)
  grid$emd <- map2_dbl(grid$population, grid$marker, function(p, m) {
    rows <- pop == p
    max_pair_emd(x[held, m][rows], batch[rows])
  })
  grid |> mutate(method = method, fold = fold, .before = 1)
}
emd_long <- map(runs, score_run) |> list_rbind()

dir.create(out_dir, recursive = TRUE)
write_csv(emd_long, file.path(out_dir, "emd_long.csv"))
emd_long |>
  summarise(median_emd = median(emd, na.rm = TRUE), n_scored = sum(!is.na(emd)), .by = c(method, fold)) |>
  write_csv(file.path(out_dir, "emd_summary.csv"))
list.files(results_dir, pattern = "status.csv", recursive = TRUE, full.names = TRUE) |>
  map(\(f) read_csv(f, show_col_types = FALSE, col_types = cols(.default = "c"))) |>
  list_rbind() |>
  write_csv(file.path(out_dir, "status_all.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
