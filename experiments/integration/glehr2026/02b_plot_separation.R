# experiments/integration/glehr2026/02b_plot_separation.R
# Purpose: show whether the control's split separation predicts how far our M is from the authors'.
# Run:     Rscript experiments/integration/glehr2026/02b_plot_separation.R <reproduce_run_dir>
# Inputs:  <reproduce_run_dir>/compare_m_M.csv   from 02_reproduce_cd3.R
# Outputs: <reproduce_run_dir>/fig_separation_vs_error.png, figure_source_data_separation_vs_error.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
plot_frame <- read.csv(file.path(run_dir, "compare_m_M.csv")) |>
  mutate(abs_log_ratio_M = abs(log_ratio_M))
write.csv(plot_frame, file.path(run_dir, "figure_source_data_separation_vs_error.csv"), row.names = FALSE)
print(as.data.frame(summarise(plot_frame, min_sep = min(separation), max_sep = max(separation), .by = channel)))

p <- ggplot(plot_frame, aes(separation, abs_log_ratio_M, colour = channel == "KrO-A")) +
  geom_hline(yintercept = c(0.05, 0.15), linetype = 2, colour = "grey60") +
  geom_point(size = 1, alpha = 0.8) +
  scale_x_log10() +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D"), labels = c("other", "KrO-A")) +
  labs(
    x = "Control split separation, gain / (SD neg + SD pos), log scale",
    y = "|log(M ours / M authors)|",
    colour = NULL
  )
ggsave(file.path(run_dir, "fig_separation_vs_error.png"), p, width = 110, height = 80, units = "mm", dpi = 200)
