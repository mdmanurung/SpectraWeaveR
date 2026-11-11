# experiments/benchmarks/10_negative-faults/plot.R
# Purpose: AUROC per check for each negative-population fault and size.
# Run:     Rscript experiments/benchmarks/10_negative-faults/plot.R <out_dir>
# Inputs:  <out_dir>/negative_long.csv   from run.R
# Outputs: <out_dir>/fig_negative_faults.png
#          <out_dir>/figure_source_data_negative.csv  mean AUROC and 95% half-width

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

out_dir <- commandArgs(trailingOnly = TRUE)[[1]]
long <- read.csv(file.path(out_dir, "negative_long.csv"))
mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(fault, size, score)
  ) |>
  mutate(case = paste0(fault, " ", size))
write.csv(mat, file.path(out_dir, "figure_source_data_negative.csv"), row.names = FALSE)
print(
  mat |>
    mutate(auroc = round(auroc, 2)) |>
    select(case, score, auroc) |>
    tidyr::pivot_wider(names_from = case, values_from = auroc) |>
    as.data.frame()
)

p <- ggplot(mat, aes(auroc, score)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)), size = 0.2) +
  facet_wrap(vars(case), nrow = 2) +
  labs(x = "AUROC (mean, 95% CI over 5 seeds)", y = NULL)
ggsave(file.path(out_dir, "fig_negative_faults.png"), p, width = 220, height = 110, units = "mm", dpi = 200)
