# experiments/integration/aurora_exvivo/07_cluster_stability_real.R
# Purpose: how stable is a k-means clustering of the real files across seeds and event counts,
#          for several k? Numbers only; tells whether cluster-level checks (sw_cluster_scores)
#          rest on a stable partition.
# Run:     Rscript experiments/integration/aurora_exvivo/07_cluster_stability_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/cluster_stability.csv   per k: check, item, value
#          <run_dir>/fig_cluster_stability.png, figure_source_data_cluster_stability.csv
#
# Every value below is illustrative; none is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
K_VALUES <- c(6, 12, 20)
SEEDS <- 1:5
SIZES <- c(5000, 20000)
N_EVENTS <- 500 # per file into the pool
NSTART <- 3
ITER_MAX <- 100

ev <- readRDS(file.path(run_dir, "events.rds"))
ch <- setdiff(colnames(ev[[1]]), "Time")
names(ev) <- paste0("f", seq_along(ev)) # file names stay out of outputs

res <- lapply(K_VALUES, function(k) {
  st <- sw_cluster_stability(ev, ch, k, SEEDS, SIZES, N_EVENTS, draw_seed = 1, NSTART, ITER_MAX)
  st$k <- k
  st$n_kmeans_warnings <- length(attr(st, "kmeans_warnings"))
  st
}) |>
  bind_rows()
write.csv(res, file.path(run_dir, "cluster_stability.csv"), row.names = FALSE)
write.csv(res, file.path(run_dir, "figure_source_data_cluster_stability.csv"), row.names = FALSE)
print(summarise(res, median = round(median(value), 3), min = round(min(value), 3), .by = c(k, check)))

p <- ggplot(res, aes(factor(k), value)) +
  geom_point(position = position_jitter(width = 0.1, height = 0, seed = 1), size = 1) +
  facet_wrap(vars(check), scales = "free_y", nrow = 1) +
  labs(x = "k (clusters)", y = "Value (ARI: 1 = same partition; moved share: 0 = same)")
ggsave(file.path(run_dir, "fig_cluster_stability.png"), p, width = 180, height = 65, units = "mm", dpi = 200)
