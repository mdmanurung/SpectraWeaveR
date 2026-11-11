#!/usr/bin/env Rscript
# Purpose: anchor plan A8/A14 EMD check: in the stripped batches, is anchor-corrected EMD below abstaining?
#          Abstaining (A14 Q1) = the arm's own corrected matrix with raw values put back in the stripped batches.
# Run:     Rscript experiments/benchmarks/06_cytether-cytobatchflagr/03c_score_strip.R <prep_dir> <strip_results> <out_dir>
# Inputs:  <prep_dir>/cells.rds, markers.csv, populations.rds (01_prep.R)
#          <strip_results>/<arm>/fold_<c>/{corrected.rds,status.csv} (02c_correct_strip.R), arm in ARMS
# Outputs (in <out_dir>, which must not exist):
#   emd_strip_long.csv     arm, version (corrected/abstain), fold, population, marker, emd
#   emd_strip_summary.csv  arm, version, fold, median_emd, n_scored
#   emd_strip_gap.csv      arm, fold, abstain and corrected median, gap = corrected - abstain, below (gap < 0)
#   fig_emd_strip.png/.pdf median EMD per fold, corrected vs abstain, one panel per arm
#   sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
})

# Parameters ---------------------------------------------------------------------------
ARMS <- c("heb_kv_bap_anc", "heb_kv_bap_anc_u") # A14 Q4
LINEAGE <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD16", "CD45RA", "CD103", "CCR6") # A14 Q3 = O4 lineage set
BIN_SIZE <- 0.1 # D1, as 03_score_emd.R
MIN_RANGE <- -100
MAX_RANGE <- 100
MANUAL_THRESHOLD <- 2 # a file needs this many cells of a population, as 03_score_emd.R

args <- commandArgs(trailingOnly = TRUE)
prep_dir <- args[[1]]
strip_dir <- args[[2]]
out_dir <- args[[3]]
if (dir.exists(out_dir)) stop(sprintf("%s exists; use a new out_dir.", out_dir), call. = FALSE)

# Inputs -------------------------------------------------------------------------------
cells <- readRDS(file.path(prep_dir, "cells.rds"))
markers <- read_csv(file.path(prep_dir, "markers.csv"), show_col_types = FALSE)
populations <- readRDS(file.path(prep_dir, "populations.rds"))
lineage <- markers$channel[match(LINEAGE, markers$label)]
stopifnot("every lineage marker is in markers.csv" = !anyNA(lineage))
breaks <- seq(MIN_RANGE, MAX_RANGE, by = BIN_SIZE)

# EMD over batch pairs with >= 1 stripped batch (A14 Q2) ------------------------------
# 1D EMD in bin units = sum |CDF1 - CDF2| (checked against emdist in 03_score_emd.R).
max_strip_pair_emd <- function(values, batch, stripped) {
  cdf <- map(split(values, batch), function(v) {
    if (length(v) < MANUAL_THRESHOLD) return(NULL)
    cumsum(hist(v, breaks = breaks, plot = FALSE)$counts) / length(v)
  }) |>
    compact()
  if (length(cdf) < 2 || !any(names(cdf) %in% stripped)) return(NA_real_)
  pairs <- combn(names(cdf), 2)
  pairs <- pairs[, pairs[1, ] %in% stripped | pairs[2, ] %in% stripped, drop = FALSE]
  max(apply(pairs, 2, \(pr) sum(abs(cdf[[pr[1]]] - cdf[[pr[2]]]))))
}

# Score each arm, corrected and abstain, per fold ---------------------------------------
score_fold <- function(arm, fold) {
  run <- file.path(strip_dir, arm, paste0("fold_", fold))
  stripped <- strsplit(read_csv(file.path(run, "status.csv"), show_col_types = FALSE)$stripped, ";")[[1]]
  corrected <- readRDS(file.path(run, "corrected.rds"))
  held <- cells$control == fold
  in_strip <- cells$batch %in% stripped
  abstain <- corrected
  abstain[in_strip, lineage] <- as.matrix(cells[in_strip, lineage])
  pop <- populations[[paste0("fold_", fold)]][held]
  batch <- cells$batch[held]
  versions <- list(corrected = corrected, abstain = abstain)
  expand_grid(version = names(versions), population = sort(unique(pop)), marker = lineage) |>
    mutate(emd = pmap_dbl(list(version, population, marker), function(v, p, m) {
      rows <- pop == p
      max_strip_pair_emd(versions[[v]][held, m][rows], batch[rows], stripped)
    })) |>
    mutate(arm = arm, fold = fold, .before = 1)
}
emd_long <- expand_grid(arm = ARMS, fold = c("tcu1", "tcu2", "tcs3", "tcs4")) |>
  pmap(score_fold) |>
  list_rbind() |>
  mutate(marker = markers$label[match(marker, markers$channel)])

summary <- emd_long |>
  summarise(median_emd = median(emd, na.rm = TRUE), n_scored = sum(!is.na(emd)), .by = c(arm, version, fold))
gap <- summary |>
  select(arm, version, fold, median_emd) |>
  pivot_wider(names_from = version, values_from = median_emd) |>
  mutate(gap = corrected - abstain, below = gap < 0)

# Write --------------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE)
write_csv(emd_long, file.path(out_dir, "emd_strip_long.csv"))
write_csv(summary, file.path(out_dir, "emd_strip_summary.csv"))
write_csv(gap, file.path(out_dir, "emd_strip_gap.csv"))

# Figure: per fold, abstain to corrected, one line per fold; one panel per arm.
fig_frame <- mutate(summary, version = factor(version, levels = c("abstain", "corrected")))
p <- ggplot(fig_frame, aes(version, median_emd, group = fold, colour = fold)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.5) +
  facet_wrap(~arm, nrow = 1) +
  expand_limits(y = 0) +
  labs(x = NULL, y = "median EMD, stripped-batch pairs, held-out donor\n(lineage markers, bin units)",
       colour = "held-out", title = "A8/A14: stripped batches, anchor-corrected vs abstain") +
  theme_bw(base_size = 8)
ggsave(file.path(out_dir, "fig_emd_strip.png"), p, width = 120, height = 70, units = "mm", dpi = 200)
ggsave(file.path(out_dir, "fig_emd_strip.pdf"), p, width = 120, height = 70, units = "mm")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
