#!/usr/bin/env Rscript
# Purpose: split the held-out donor's EMD into location and shape for one fold, and score a prototype
#          heb_kv_bap_scale arm (HEB shift plus a reference-cell MAD scale), to test whether scale explains
#          cycombine_ref's lead over heb_kv_bap (plan T6 S1-S4, user 2026-11-10).
# Run:     pixi run Rscript experiments/benchmarks/06_cytether-cytobatchflagr/08_emd_decompose.R <prep_dir> <results_dir> <out_dir> <fold>
# Inputs:  <prep_dir>/cells.rds, markers.csv, populations.rds (01_prep.R)
#          <results_dir>/<method>/fold_<fold>/corrected.rds for METHODS (02_correct.R)
# Outputs (in <out_dir>, created if absent; the fold file must not exist):
#   fold_<fold>.csv   method, fold, population, marker, emd_total, emd_shape (max over batch pairs)
#   scale_<fold>.csv  batch, population, marker, scale factor s, n calibrating cells, whether s was applied
#   sessionInfo_<fold>.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
})

# Parameters -------------------------------------------------------------------------------
METHODS <- c("none", "heb_kv_bap", "cycombine_ref", "cytonorm") # Q4, inherited from emd_scores_v5
BIN_SIZE <- 0.1 # 03_score_emd.R (CytoNorm 2020 defaults, D1)
MIN_RANGE <- -100
MAX_RANGE <- 100
MANUAL_THRESHOLD <- 2
ANCHOR <- "batch01" # Q8
MIN_CELLS <- 10 # G1: calibrating cells per batch x population needed for a scale factor (Q2)

args <- commandArgs(trailingOnly = TRUE)
prep_dir <- args[[1]]
results_dir <- args[[2]]
out_dir <- args[[3]]
fold <- args[[4]]
out_file <- file.path(out_dir, paste0("fold_", fold, ".csv"))
if (file.exists(out_file)) {
  stop(sprintf("%s exists; use a new out_dir.", out_file), call. = FALSE)
}

cells <- readRDS(file.path(prep_dir, "cells.rds"))
markers <- read_csv(
  file.path(prep_dir, "markers.csv"),
  show_col_types = FALSE
)$channel
pop <- readRDS(file.path(prep_dir, "populations.rds"))[[paste0("fold_", fold)]]
stopifnot("populations align with cells" = length(pop) == nrow(cells))
held <- cells$control == fold
calib <- !held # the other three controls calibrate (02_correct.R)
breaks <- seq(MIN_RANGE, MAX_RANGE, by = BIN_SIZE)

# EMD as in 03_score_emd.R; centre = TRUE moves each file to median 0 first (shape only, Q1) -----
max_pair_emd <- function(values, batch, centre) {
  cdf <- map(split(values, batch), function(v) {
    if (length(v) < MANUAL_THRESHOLD) {
      return(NULL)
    }
    if (centre) {
      v <- v - median(v)
    }
    cumsum(hist(v, breaks = breaks, plot = FALSE)$counts) / length(v)
  }) |>
    compact()
  if (length(cdf) < 2) {
    return(NA_real_)
  }
  m <- do.call(cbind, cdf)
  pairs <- combn(ncol(m), 2)
  max(vapply(
    seq_len(ncol(pairs)),
    \(k) sum(abs(m[, pairs[1, k]] - m[, pairs[2, k]])),
    numeric(1)
  ))
}

# Prototype scale arm (Q2): on the HEB-corrected values, per batch x population x marker,
# x' = m_b + s * (x - m_b), m_b = median of the batch's calibrating cells,
# s = MAD(anchor-batch calibrating cells) / MAD(batch calibrating cells); donors pooled; s = 1 below MIN_CELLS.
heb <- readRDS(file.path(
  results_dir,
  "heb_kv_bap",
  paste0("fold_", fold),
  "corrected.rds"
))
scaled <- heb
key <- paste(cells$batch, pop, sep = "\r")
anchor_rows <- split(
  which(calib & cells$batch == ANCHOR),
  pop[calib & cells$batch == ANCHOR]
)
scale_log <- list()
for (k in unique(key[held])) {
  rows <- which(key == k)
  b <- cells$batch[rows[[1]]]
  p <- pop[rows[[1]]]
  cal <- rows[calib[rows]]
  arows <- anchor_rows[[p]]
  ok <- b != ANCHOR && length(cal) >= MIN_CELLS && length(arows) >= MIN_CELLS
  for (m in markers) {
    s <- if (ok) mad(heb[arows, m]) / mad(heb[cal, m]) else NA_real_
    use <- is.finite(s) && s > 0
    if (use) {
      mb <- median(heb[cal, m])
      scaled[rows, m] <- mb + s * (heb[rows, m] - mb)
    }
    scale_log[[length(scale_log) + 1]] <- tibble(
      batch = b,
      population = p,
      marker = m,
      s = s,
      n_calib = length(cal),
      applied = use
    )
  }
}

# Score every method on the held-out donor --------------------------------------------------
corrected <- c(
  set_names(
    map(METHODS, \(mt) {
      readRDS(file.path(
        results_dir,
        mt,
        paste0("fold_", fold),
        "corrected.rds"
      ))
    }),
    METHODS
  ),
  list(heb_kv_bap_scale = scaled)
)
stopifnot(
  "corrected matrices keep the cell rows" = all(
    map_int(corrected, nrow) == nrow(cells)
  )
)
pop_h <- pop[held]
batch_h <- cells$batch[held]
grid <- tidyr::expand_grid(population = sort(unique(pop_h)), marker = markers)
scores <- imap(corrected, function(x, mt) {
  xh <- x[held, , drop = FALSE]
  grid |>
    mutate(
      method = mt,
      fold = fold,
      emd_total = map2_dbl(population, marker, \(p, m) {
        max_pair_emd(xh[pop_h == p, m], batch_h[pop_h == p], FALSE)
      }),
      emd_shape = map2_dbl(population, marker, \(p, m) {
        max_pair_emd(xh[pop_h == p, m], batch_h[pop_h == p], TRUE)
      }),
      .before = 1
    )
}) |>
  list_rbind()

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(scores, out_file)
write_csv(
  list_rbind(scale_log),
  file.path(out_dir, paste0("scale_", fold, ".csv"))
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, paste0("sessionInfo_", fold, ".txt"))
)
