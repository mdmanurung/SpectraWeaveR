# Purpose: redraw the swap-test (37b) and variant-count (39b) figures from their saved tables, with
#          readable fonts and no captions. No FCS file is read and nothing is recomputed.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/41_replot_scores.R <source_run_dir> <out_dir>
#          from the repo root; out_dir must not exist.
# Inputs:  <source_run_dir>/scores.csv           one row per arm, donor, population, channel [37b or 39b]
#          <source_run_dir>/bv510_events.csv.gz  BV510 arcsinh per event [37b or 39b]
#          experiments/benchmarks/05_bmv-unmix-compare/helper41_plot-scores.R  plot_bv510_density(), plot_scores()
# Outputs: <out_dir>/fig_bv510_density.png  BV510 density per arm, donor x population
#          <out_dir>/fig_scores.png         three scores per arm, donors joined by lines
#          <out_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/helper41_plot-scores.R") # plot_bv510_density(), plot_scores()
theme_set(ggpubr::theme_pubr(base_size = 11))

# Parameters ------------------------------------------------------------------------------
NEG_CUT <- 1 # dotted line, as in 37b and 39b (user, 2026-11-10)
# colours and order as in 37b and 39b; one palette for both runs
ARM_COLORS <- c(
  as188_joint = "grey30",
  s102_c188 = "#D55E00",
  s188_c102 = "#0072B2",
  swap_af102 = "#E69F00",
  swap_var102 = "#009E73",
  var102_n26 = "#CC79A7",
  var188_n100 = "#56B4E9",
  var102_bv510only = "#F0E442",
  var102_except_bv510 = "#999933"
)

cli_args <- commandArgs(trailingOnly = TRUE)
source_dir <- cli_args[[1]]
out_dir <- cli_args[[2]]
stopifnot("out_dir must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)

# Read the saved tables -------------------------------------------------------------------
# an arm missing from ARM_COLORS becomes NA and fails the factor check below
scores <- read_csv(file.path(source_dir, "scores.csv"), show_col_types = FALSE)
dens <- read_csv(
  file.path(source_dir, "bv510_events.csv.gz"),
  show_col_types = FALSE
)
arms <- names(ARM_COLORS)[names(ARM_COLORS) %in% scores$arm]
stopifnot(
  "every arm has a colour" = setequal(arms, unique(c(scores$arm, dens$arm)))
)
scores <- scores |> mutate(arm = factor(arm, levels = arms))
dens <- dens |> mutate(arm = factor(arm, levels = arms))

# Figures ---------------------------------------------------------------------------------
ggsave(
  file.path(out_dir, "fig_bv510_density.png"),
  plot_bv510_density(dens, ARM_COLORS[arms], NEG_CUT),
  width = 220,
  height = 120,
  units = "mm",
  dpi = 200,
  bg = "white"
)
ggsave(
  file.path(out_dir, "fig_scores.png"),
  plot_scores(scores, ARM_COLORS[arms]),
  width = 200,
  height = 200,
  units = "mm",
  dpi = 200,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
