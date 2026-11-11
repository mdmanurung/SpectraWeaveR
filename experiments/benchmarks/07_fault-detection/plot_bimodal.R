# experiments/benchmarks/07_fault-detection/plot_bimodal.R
# Purpose: figure for the file-level two-mode split in benchmark 07: AUROC per tool over seeds,
#          with the new bimodality_z (G20) and dip_z (G35) arms beside the existing tools and the
#          shuffled-label null.
# Run:     Rscript experiments/benchmarks/07_fault-detection/plot_bimodal.R <run_dir>
# Inputs:  <run_dir>/detection_long.csv   from run.R (a run with the dip_z arm)
# Outputs: <run_dir>/fig_bimodal_tools.png, figure_source_data_bimodal_tools.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

BAR <- 0.9 # D1 mean AUROC bar

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]

# Split fault and null, per seed --------------------------------------------------------------
long <- read.csv(file.path(run_dir, "detection_long.csv")) |>
  filter(level == "file", fault %in% c("bimodal", "none (shuffled labels)")) |>
  mutate(
    fault = if_else(
      fault == "bimodal",
      "two-mode split (+/-0.5)",
      "null (shuffled labels)"
    ),
    new_arm = tool %in% c("bimodality_z", "dip_z")
  )
order <- long |>
  filter(fault == "two-mode split (+/-0.5)") |>
  summarise(m = mean(auroc, na.rm = TRUE), .by = tool) |>
  arrange(m)
long <- mutate(long, tool = factor(tool, levels = order$tool))
write.csv(
  long,
  file.path(run_dir, "figure_source_data_bimodal_tools.csv"),
  row.names = FALSE
)

p <- ggplot(long, aes(auroc, tool, colour = new_arm)) +
  geom_vline(xintercept = c(0.5, BAR), linetype = 2, colour = "grey60") +
  geom_point(
    position = position_jitter(height = 0.15, width = 0, seed = 1),
    size = 0.6,
    alpha = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 3,
    size = 2.5,
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(`TRUE` = "#D55E00", `FALSE` = "grey40"),
    labels = c(`TRUE` = "new arm", `FALSE` = "existing")
  ) +
  facet_wrap(vars(fault)) +
  labs(
    x = "AUROC (dots: 10 seeds; cross: mean)",
    y = NULL,
    colour = NULL,
    caption = "Benchmark 07, simulated files; dashed lines: chance 0.5 and the bar 0.9"
  )
ggsave(
  file.path(run_dir, "fig_bimodal_tools.png"),
  p,
  width = 170,
  height = 90,
  units = "mm",
  dpi = 200
)
