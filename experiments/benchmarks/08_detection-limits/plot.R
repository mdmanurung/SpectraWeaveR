# experiments/benchmarks/08_detection-limits/plot.R
# Purpose: AUROC against effect size for each tool, one panel per sweep.
# Run:     Rscript experiments/benchmarks/08_detection-limits/plot.R <out_dir>
# Inputs:  <out_dir>/limits_matrix.csv   from run.R
# Outputs: <out_dir>/fig_detection_limits.png
#          <out_dir>/figure_source_data_limits.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))
out_dir <- commandArgs(trailingOnly = TRUE)[[1]]
mat <- read.csv(file.path(out_dir, "limits_matrix.csv"))
write.csv(mat, file.path(out_dir, "figure_source_data_limits.csv"), row.names = FALSE)
p <- ggplot(mat, aes(size, auroc, colour = tool)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = auroc - ci95_half, ymax = auroc + ci95_half), size = 0.25) +
  geom_line() +
  facet_wrap(vars(sweep), scales = "free_x", nrow = 1) +
  scale_y_continuous(limits = c(0.3, 1)) +
  labs(x = "Effect size (units differ by panel)", y = "AUROC (mean, 95% CI over seeds)", colour = NULL)
ggsave(file.path(out_dir, "fig_detection_limits.png"), p, width = 180, height = 70, units = "mm", dpi = 200)
