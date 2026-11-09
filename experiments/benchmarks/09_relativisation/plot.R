# experiments/benchmarks/09_relativisation/plot.R
# Purpose: show batch EMD per channel and OTD per file before and after relativisation, by CD3
#          donor variation.
# Run:     Rscript experiments/benchmarks/09_relativisation/plot.R <out_dir>
# Inputs:  <out_dir>/relativisation_long.csv, otd_long.csv, gain_recovery.csv   from run.R
# Outputs: <out_dir>/fig_relativisation.png   EMD panel (top) and OTD panel (bottom)
#          <out_dir>/figure_source_data_relativisation_emd.csv
#          <out_dir>/figure_source_data_relativisation_otd.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

out_dir <- commandArgs(trailingOnly = TRUE)[[1]]
VERSIONS <- c(
  "observed",
  "relativised per sample",
  "relativised per batch (E1)",
  "truth (no instrument)"
)
label_sd <- \(v) sprintf("CD3 donor SD (log) = %.1f", v)

emd_df <- read.csv(file.path(out_dir, "relativisation_long.csv")) |>
  mutate(version = factor(version, levels = VERSIONS), panel = label_sd(cd3_donor_sd))
# OTD: one value per file, summarised to the median per seed so points are seeds.
otd_df <- read.csv(file.path(out_dir, "otd_long.csv")) |>
  summarise(otd = median(otd), all_converged = all(converged), .by = c(seed, cd3_donor_sd, version)) |>
  mutate(version = factor(version, levels = VERSIONS), panel = label_sd(cd3_donor_sd))
write.csv(emd_df, file.path(out_dir, "figure_source_data_relativisation_emd.csv"), row.names = FALSE)
write.csv(otd_df, file.path(out_dir, "figure_source_data_relativisation_otd.csv"), row.names = FALSE)

cat("median largest batch EMD per channel:\n")
print(as.data.frame(summarise(emd_df, median = round(median(max_pair_emd), 3), .by = c(cd3_donor_sd, version))))
cat("median per-seed median OTD:\n")
print(as.data.frame(summarise(otd_df, median = signif(median(otd), 3), .by = c(cd3_donor_sd, version))))
cat("all OTD solves converged:", all(otd_df$all_converged), "\n")
recovery <- read.csv(file.path(out_dir, "gain_recovery.csv"))
print(summarise(recovery, cor = median(cor_log_gain), .by = cd3_donor_sd))

p_emd <- ggplot(emd_df, aes(version, max_pair_emd)) +
  geom_boxplot(outlier.shape = NA, colour = "grey50") +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 0.5) +
  facet_wrap(vars(panel), nrow = 1) +
  labs(x = NULL, y = "Largest batch-pair EMD\nper channel (arcsinh)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_otd <- ggplot(otd_df, aes(version, otd)) +
  geom_boxplot(outlier.shape = NA, colour = "grey50") +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 0.8) +
  facet_wrap(vars(panel), nrow = 1) +
  labs(x = NULL, y = "Median OTD of a file to\nthe pooled reference (per seed)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p <- patchwork::wrap_plots(p_emd, p_otd, ncol = 1)
ggsave(file.path(out_dir, "fig_relativisation.png"), p, width = 170, height = 130, units = "mm", dpi = 200)
