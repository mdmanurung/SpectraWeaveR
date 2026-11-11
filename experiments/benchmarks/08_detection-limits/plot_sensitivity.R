# experiments/benchmarks/08_detection-limits/plot_sensitivity.R
# Purpose: AUROC at three fixed fault sizes when k (clusters) or the number of faulty files changes.
# Run:     Rscript experiments/benchmarks/08_detection-limits/plot_sensitivity.R <runs_dir> <out_dir>
# Inputs:  <runs_dir>/2026-11-09_limits-k*-bad*/limits_matrix.csv   from run.R with SW_K and SW_N_BAD
#          <runs_dir>/2026-11-09_detection-limits/limits_matrix.csv  baseline (k 6, 12 faulty files)
# Outputs: <out_dir>/fig_sensitivity.png, <out_dir>/figure_source_data_sensitivity.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))
args <- commandArgs(trailingOnly = TRUE)
runs_dir <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
settings <- tibble::tribble(
  ~dir, ~k, ~n_bad,
  "2026-11-09_limits-k3-bad12", 3, 12,
  "2026-11-09_limits-k12-bad12", 12, 12,
  "2026-11-09_limits-k6-bad4", 6, 4,
  "2026-11-09_limits-k6-bad24", 6, 24,
  "2026-11-09_detection-limits", 6, 12
)
keep <- c("shift:0.3", "stretch:1.1", "rare_pop_fraction:0.1")
res <- settings |>
  rowwise() |>
  reframe(
    k = k,
    n_bad = n_bad,
    read.csv(file.path(runs_dir, dir, "limits_matrix.csv"))
  ) |>
  mutate(case = paste0(sweep, ":", size)) |>
  filter(case %in% keep) |>
  mutate(setting = paste0("k=", k, ", ", n_bad, " faulty"))
write.csv(res, file.path(out_dir, "figure_source_data_sensitivity.csv"), row.names = FALSE)
p <- ggplot(res, aes(setting, auroc, colour = tool)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = auroc - ci95_half, ymax = auroc + ci95_half), size = 0.25,
                  position = position_dodge(width = 0.5)) +
  facet_wrap(vars(case), nrow = 1) +
  labs(x = NULL, y = "AUROC (mean, 95% CI over seeds)", colour = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "fig_sensitivity.png"), p, width = 180, height = 80, units = "mm", dpi = 200)
