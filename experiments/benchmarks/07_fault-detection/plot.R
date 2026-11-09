# experiments/benchmarks/07_fault-detection/plot.R
# Purpose: draw the detection matrix (mean AUROC per fault x tool) with seeds as points
#          in a second panel, so a reader sees the spread and not only the mean.
# Run:     Rscript experiments/benchmarks/07_fault-detection/plot.R <run_dir>
# Inputs:  <run_dir>/detection_long.csv, <run_dir>/detection_matrix.csv (from run.R)
# Outputs: <run_dir>/detection_matrix.png    heatmap, mean AUROC, text = value
#          <run_dir>/detection_seeds.png     per-seed AUROC for the best tool of each fault
#          <run_dir>/figure_source_data.csv  the numbers plotted, for both panels

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
mat <- read.csv(file.path(run_dir, "detection_matrix.csv"))
long <- read.csv(file.path(run_dir, "detection_long.csv"))

# Plot frames ---------------------------------------------------------------
fault_order <- c(
  "scale_all",
  "shift",
  "stretch",
  "bimodal",
  "debris",
  "saturation",
  "dead_channel",
  "nonfinite",
  "clog",
  "spillover",
  "batch_shift",
  "batch_stretch",
  "batch_bimodal",
  "batch_drift",
  "none (shuffled labels)"
)
heat <- mat |>
  filter(!is.na(auroc)) |>
  mutate(
    fault = factor(fault, levels = rev(fault_order)),
    label = sprintf("%.2f", auroc)
  )
write.csv(heat, file.path(run_dir, "figure_source_data.csv"), row.names = FALSE)

p_heat <- ggplot(heat, aes(tool, fault, fill = auroc)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label), size = 2.2) +
  scale_fill_gradient2(
    low = "#b2182b",
    mid = "white",
    high = "#2166ac",
    midpoint = 0.5,
    limits = c(0, 1)
  ) +
  labs(
    x = "Tool (score used)",
    y = "Injected fault",
    fill = "Mean AUROC\n(10 seeds)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(run_dir, "detection_matrix.png"),
  p_heat,
  width = 140,
  height = 110,
  units = "mm",
  dpi = 200
)

# Best tool per fault, with each seed as a point.
best <- mat |>
  filter(fault != "none (shuffled labels)") |>
  slice_max(auroc, by = fault, n = 1, with_ties = FALSE) |>
  select(fault, tool)
seeds <- long |>
  inner_join(best, by = c("fault", "tool")) |>
  mutate(fault = factor(fault, levels = rev(fault_order)))
p_seed <- ggplot(seeds, aes(auroc, fault)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_point(
    position = position_jitter(height = 0.15, seed = 1),
    size = 1,
    alpha = 0.7
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    fill = "#d95f02",
    size = 2
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "AUROC of the best tool per fault (points: seeds; diamond: mean)",
    y = NULL
  )
ggsave(
  file.path(run_dir, "detection_seeds.png"),
  p_seed,
  width = 120,
  height = 90,
  units = "mm",
  dpi = 200
)
cat("saved figures in", run_dir, "\n")
