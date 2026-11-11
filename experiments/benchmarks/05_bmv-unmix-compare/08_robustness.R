# Purpose: check whether the arm ordering survives leaving out each donor and each channel.
#          Deterministic: no resampling, no new parameters; the readouts are the ones in 02 and 02b.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/08_robustness.R <run_dir>
#          from the repo root.
# Inputs:  <run_dir>/negative_spread.csv, stain_index.csv, negative_cells.csv [02_metrics.R, 02b_shape.R]
# Outputs: <run_dir>/robustness_leave_out.csv  per readout, left-out unit, arm: score and rank
#          <run_dir>/robustness_summary.csv    per readout, arm: rank range across all leave-outs

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

REFERENCE_ARM <- "as188_joint" # same reference as 02 (user, 2026-11-09)

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_files <- file.path(run_dir, c("robustness_leave_out.csv", "robustness_summary.csv"))
stopifnot("08 outputs must not exist" = !any(file.exists(out_files)))

spread <- read_csv(file.path(run_dir, "negative_spread.csv"), show_col_types = FALSE)
stain <- read_csv(file.path(run_dir, "stain_index.csv"), show_col_types = FALSE)
negc <- read_csv(file.path(run_dir, "negative_cells.csv"), show_col_types = FALSE)

# Per arm, sample, channel score against the reference arm; lower is better for all three.
scores <- bind_rows(
  spread |>
    mutate(score = rsd_neg, readout = "negative_rsd") |>
    select(readout, arm, sample, unit = channel, score),
  stain |>
    mutate(score = -stain_index, readout = "neg_stain_index") |>
    select(readout, arm, sample, unit = marker, score),
  negc |>
    mutate(score = frac_below, readout = "negative_cell_frac") |>
    select(readout, arm, sample, unit = channel, score)
) |>
  mutate(
    ref = score[match(paste(readout, REFERENCE_ARM, sample, unit), paste(readout, arm, sample, unit))],
    .by = readout
  )

# Arm score: ratio readouts use the median log2 ratio to the reference, stain index the median difference.
arm_score <- function(d) {
  d |>
    summarise(
      value = if (first(readout) == "neg_stain_index") median(score - ref, na.rm = TRUE) else median(log2(score / ref), na.rm = TRUE),
      .by = c(readout, arm)
    ) |>
    mutate(rank = rank(value), .by = readout)
}

samples <- sort(unique(scores$sample))
units <- scores |> distinct(readout, unit)
loo <- bind_rows(
  tibble(left_out = "none", kind = "none") |> mutate(res = list(arm_score(scores))),
  tibble(left_out = samples, kind = "donor") |>
    mutate(res = map(left_out, \(s) arm_score(filter(scores, sample != s)))),
  units |>
    mutate(left_out = unit, kind = "channel") |>
    distinct(left_out, kind) |>
    mutate(res = map(left_out, \(u) arm_score(filter(scores, unit != u))))
) |>
  unnest(res)

summ <- loo |>
  summarise(
    rank_full = rank[kind == "none"],
    rank_min = min(rank),
    rank_max = max(rank),
    .by = c(readout, arm)
  )

write_csv(loo, out_files[[1]])
write_csv(summ, out_files[[2]])
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo_08.txt"))
