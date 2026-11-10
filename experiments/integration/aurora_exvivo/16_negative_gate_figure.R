# experiments/integration/aurora_exvivo/16_negative_gate_figure.R
# Purpose: show which channels the separation gate (G24-G27) keeps for the negative checks:
#          median GMM separation over files per channel against min_separation.
# Run:     Rscript experiments/integration/aurora_exvivo/16_negative_gate_figure.R <run_dir>
# Inputs:  <run_dir>/negative_separation.csv   from 06_negative_real.R (gated run)
# Outputs: <run_dir>/fig_negative_gate.png, figure_source_data_negative_gate.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

MIN_SEPARATION <- 3 # G24, as in 06

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]

# Gate per channel ------------------------------------------------------------------------------
gate <- read.csv(file.path(run_dir, "negative_separation.csv")) |>
  mutate(
    scored = if_else(separation >= MIN_SEPARATION, "scored", "gated (NA)"),
    channel = reorder(channel, separation)
  )
write.csv(
  gate,
  file.path(run_dir, "figure_source_data_negative_gate.csv"),
  row.names = FALSE
)

p <- ggplot(gate, aes(separation, channel, colour = scored)) +
  geom_vline(xintercept = MIN_SEPARATION, linetype = 2, colour = "grey50") +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c(scored = "black", `gated (NA)` = "grey65")) +
  labs(
    x = "Median GMM separation over files",
    y = NULL,
    colour = NULL,
    caption = sprintf(
      "Aurora Batch1 donor files; dashed line: min_separation = %g; %d of %d channels scored",
      MIN_SEPARATION,
      sum(gate$scored == "scored"),
      nrow(gate)
    )
  )
ggsave(
  file.path(run_dir, "fig_negative_gate.png"),
  p,
  width = 110,
  height = 130,
  units = "mm",
  dpi = 200
)
