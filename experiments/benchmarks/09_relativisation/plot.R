# experiments/benchmarks/09_relativisation/plot.R
# Purpose: show batch EMD per channel before and after relativisation, by CD3 donor variation.
# Run:     Rscript experiments/benchmarks/09_relativisation/plot.R <out_dir>
# Inputs:  <out_dir>/relativisation_long.csv, gain_recovery.csv   from run.R
# Outputs: <out_dir>/fig_relativisation.png
#          <out_dir>/figure_source_data_relativisation.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

out_dir <- commandArgs(trailingOnly = TRUE)[[1]]
long <- read.csv(file.path(out_dir, "relativisation_long.csv"))
recovery <- read.csv(file.path(out_dir, "gain_recovery.csv"))

# One point per seed and channel; the panel title gives the gain recovery.
plot_df <- long |>
  mutate(
    version = factor(version, levels = c("observed", "relativised", "truth (no instrument)")),
    panel = sprintf("CD3 donor SD (log) = %.1f", cd3_donor_sd)
  )
write.csv(plot_df, file.path(out_dir, "figure_source_data_relativisation.csv"), row.names = FALSE)
cat("gain recovery (cor of log batch gain ratio with truth), median per CD3 SD:\n")
print(summarise(recovery, cor = median(cor_log_gain), .by = cd3_donor_sd))

p <- ggplot(plot_df, aes(version, max_pair_emd)) +
  geom_boxplot(outlier.shape = NA, colour = "grey50") +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 0.6) +
  facet_wrap(vars(panel), nrow = 1) +
  labs(
    x = NULL,
    y = "Largest pairwise batch EMD per channel\n(arcsinh scale; points = seed x channel)"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(
  file.path(out_dir, "fig_relativisation.png"),
  p,
  width = 150,
  height = 80,
  units = "mm",
  dpi = 200
)
