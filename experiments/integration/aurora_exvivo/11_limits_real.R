# experiments/integration/aurora_exvivo/11_limits_real.R
# Purpose: detection limits on real files for the two fault classes that missed the AUROC bar at
#          one size in 10_file_faults_real.R: a marker split into two modes, and a shift inside
#          the smallest cluster. Sweeps the fault size and scores sd_z, a per-channel EMD to the
#          pooled cohort (robust z over files) and, for the rare population, cluster_median_z.
# Run:     Rscript experiments/integration/aurora_exvivo/11_limits_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/limits_real_long.csv
#          <run_dir>/fig_limits_real.png, figure_source_data_limits_real.csv
#
# Every value below is illustrative; none is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
BATCH <- "Batch1"
SEEDS <- 1:5
FRACTION_BAD <- 0.25
SPLIT_SD <- c(0.5, 1, 1.5) # half-distance between the two modes, in marker SD
RARE_SD <- c(1, 2, 3) # shift inside the smallest cluster, in marker SD
CLUSTER_K <- 6 # as in 10_file_faults_real.R
CLUSTER_EVENTS <- 2000
EMD_EVENTS <- 2000
EMD_POOL <- 20000

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
base <- ev[info$file[info$batch == BATCH & grepl("BMV_[0-9]+[A-Z]", info$file)]]
names(base) <- paste0("f", seq_along(base))
ch <- setdiff(colnames(base[[1]]), "Time")
pooled <- do.call(rbind, base)
marker_sd <- apply(pooled[, ch], 2, sd)
set.seed(1)
km <- stats::kmeans(pooled[sample.int(nrow(pooled), 20000), ch], CLUSTER_K, nstart = 3, iter.max = 100)
smallest <- which.min(km$size)
nearest <- function(x) {
  d2 <- outer(rowSums(x^2), rowSums(km$centers^2), "+") - 2 * tcrossprod(x, km$centers)
  max.col(-d2, ties.method = "first")
}

auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}
robust_z <- function(x) (x - median(x)) / (1.4826 * mad(x, constant = 1))

grid <- bind_rows(
  tibble(fault = "bimodal", size = SPLIT_SD),
  tibble(fault = "rare_pop", size = RARE_SD)
)
long <- lapply(SEEDS, function(seed) {
  lapply(seq_len(nrow(grid)), function(g) {
    fault <- grid$fault[g]
    size <- grid$size[g]
    set.seed(seed * 1000 + g)
    bad <- seq_along(base) %in% sample(seq_along(base), round(FRACTION_BAD * length(base)))
    b <- sample(ch, 1)
    files <- base
    files[bad] <- lapply(files[bad], function(x) {
      if (fault == "bimodal") {
        x[, b] <- x[, b] + sample(c(-size, size), nrow(x), replace = TRUE) * marker_sd[[b]]
      } else {
        hit <- nearest(x[, ch]) == smallest
        x[hit, b] <- x[hit, b] + size * marker_sd[[b]]
      }
      x
    })
    sa <- sw_fcs_signal_audit(files, ch, NULL, NULL, time_channel = NULL, time_bin = NULL)
    sd_z <- sa[sa$check == "sd_z" & sa$channel == b, ]
    af <- sw_anomaly_features(
      lapply(files, \(x) x[, b, drop = FALSE]),
      channels = b,
      features = "emd",
      quantile_step = NULL,
      n_events = EMD_EVENTS,
      seed = seed,
      scale = "OMIQ-transformed",
      aggregate_from = "all",
      agg_size = EMD_POOL,
      reference = NULL
    )
    scores <- list(
      sd_z = abs(sd_z$value[match(names(files), sd_z$file)]),
      emd_z = abs(robust_z(af[[paste0(b, "_emd")]][match(names(files), af$file)]))
    )
    if (fault == "rare_pop") {
      cs <- sw_cluster_scores(files, ch, CLUSTER_K, CLUSTER_EVENTS, seed, 3, 30)
      d <- cs[cs$check == "median_z" & cs$marker == b, ]
      top <- tapply(abs(d$value), d$file, \(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
      scores$cluster_median_z <- as.numeric(top[names(files)])
    }
    tibble(seed = seed, fault = fault, size = size, check = names(scores),
      auroc = vapply(scores, \(v) auroc(v, bad), numeric(1)))
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(long, file.path(run_dir, "limits_real_long.csv"), row.names = FALSE)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(fault, size, check)
  )
write.csv(mat, file.path(run_dir, "figure_source_data_limits_real.csv"), row.names = FALSE)
print(as.data.frame(mutate(mat, across(c(auroc, ci95_half), \(x) round(x, 2)))))
p <- ggplot(mat, aes(size, auroc, colour = check)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = auroc - ci95_half, ymax = pmin(auroc + ci95_half, 1)), size = 0.2) +
  geom_line() +
  facet_wrap(vars(fault), scales = "free_x") +
  labs(x = "Fault size (marker SD)", y = "AUROC on real files (mean, 95% CI over 5 seeds)", colour = NULL)
ggsave(file.path(run_dir, "fig_limits_real.png"), p, width = 150, height = 70, units = "mm", dpi = 200)
