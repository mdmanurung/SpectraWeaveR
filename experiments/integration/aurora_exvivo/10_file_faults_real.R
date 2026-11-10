# experiments/integration/aurora_exvivo/10_file_faults_real.R
# Purpose: inject the file-level fault classes of benchmark 07 into 25% of the real donor files of
#          one batch, sized in units of each marker's real SD, and score the checks that target
#          them. Fills the "not tested" rows of the coverage table in reference/diagnostics.qmd.
# Run:     Rscript experiments/integration/aurora_exvivo/10_file_faults_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/file_faults_long.csv
#          <run_dir>/fig_file_faults_real.png, figure_source_data_file_faults_real.csv
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
TIME_BINS <- 50
CLUSTER_K <- 6 # stable on these files (fig_cluster_stability.png: seed-pair ARI 1.00)
CLUSTER_EVENTS <- 2000
FAULT_NAMES <- c(
  "saturation", # 5% of events in one channel set to the cohort maximum of that channel
  "dead_channel", # one channel constant
  "nonfinite", # 2% of events NA in one channel
  "clog", # 80% of events dropped in the middle 30% of the run
  "time_drift", # one channel rises by 1 SD over the run
  "spillover", # channel B += 0.3 x centred channel A, all events
  "debris", # 10% of events replaced by low values in every channel
  "bimodal", # one channel split +-0.5 SD at random
  "rare_pop" # one channel +1 SD in the events of the smallest k-means cluster
)

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
base <- ev[info$file[info$batch == BATCH & grepl("BMV_[0-9]+[A-Z]", info$file)]]
names(base) <- paste0("f", seq_along(base))
ch <- setdiff(colnames(base[[1]]), "Time")
pooled <- do.call(rbind, base)
marker_sd <- apply(pooled[, ch], 2, sd)
cohort_max <- apply(pooled[, ch], 2, max)
cohort_low <- apply(pooled[, ch], 2, stats::quantile, probs = 0.01)
time_bin <- as.numeric(median(vapply(base, \(m) diff(range(m[, "Time"])), numeric(1))) / TIME_BINS)
# Clusters for the rare-population fault, from the clean pool (seeded).
set.seed(1)
km <- stats::kmeans(pooled[sample.int(nrow(pooled), 20000), ch], CLUSTER_K, nstart = 3, iter.max = 100)
smallest <- which.min(km$size)
nearest <- function(x) {
  d2 <- outer(rowSums(x^2), rowSums(km$centers^2), "+") - 2 * tcrossprod(x, km$centers)
  max.col(-d2, ties.method = "first")
}

inject <- function(x, fault, a, b) {
  t_frac <- (x[, "Time"] - min(x[, "Time"])) / diff(range(x[, "Time"]))
  switch(
    fault,
    saturation = {
      idx <- sample.int(nrow(x), round(0.05 * nrow(x)))
      x[idx, b] <- cohort_max[[b]]
    },
    dead_channel = {
      x[, b] <- stats::median(x[, b])
    },
    nonfinite = {
      x[sample.int(nrow(x), round(0.02 * nrow(x))), b] <- NA
    },
    clog = {
      drop <- t_frac > 0.35 & t_frac < 0.65 & stats::runif(nrow(x)) < 0.8
      x <- x[!drop, ]
    },
    time_drift = {
      x[, b] <- x[, b] + marker_sd[[b]] * t_frac
    },
    spillover = {
      x[, b] <- x[, b] + 0.3 * (x[, a] - stats::median(x[, a]))
    },
    debris = {
      idx <- sample.int(nrow(x), round(0.1 * nrow(x)))
      x[idx, ch] <- matrix(cohort_low[ch], length(idx), length(ch), byrow = TRUE) +
        stats::rnorm(length(idx) * length(ch), 0, 0.05 * marker_sd[ch])
    },
    bimodal = {
      x[, b] <- x[, b] + sample(c(-0.5, 0.5), nrow(x), replace = TRUE) * marker_sd[[b]]
    },
    rare_pop = {
      hit <- nearest(x[, ch]) == smallest
      x[hit, b] <- x[hit, b] + marker_sd[[b]]
    }
  )
  x
}

auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}
top <- \(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)

long <- lapply(SEEDS, function(seed) {
  lapply(FAULT_NAMES, function(fault) {
    set.seed(seed * 1000 + match(fault, FAULT_NAMES))
    bad <- seq_along(base) %in% sample(seq_along(base), round(FRACTION_BAD * length(base)))
    pick <- sample(ch, 2)
    files <- base
    files[bad] <- lapply(files[bad], inject, fault = fault, a = pick[1], b = pick[2])
    sa <- sw_fcs_signal_audit(
      files,
      ch,
      range_max = cohort_max,
      range_min = NULL,
      time_channel = "Time",
      time_bin = time_bin
    )
    # Hit channel's row for channel checks; the file's single row for file-level checks.
    row_of <- function(check, item) {
      d <- sa[sa$check == check & (is.na(sa$channel) | sa$channel %in% item), ]
      as.numeric(tapply(abs(d$value), d$file, top)[names(files)])
    }
    scores <- list(
      saturation_share = row_of("saturation_share", pick[2]),
      zero_variance = row_of("zero_variance", pick[2]),
      nonfinite_share = row_of("nonfinite_share", pick[2]),
      rate_cv = row_of("rate_cv", NA),
      max_gap = row_of("max_gap", NA),
      time_trend = row_of("time_trend", pick[2]),
      corr_shift = row_of("corr_shift", NA),
      median_z = row_of("median_z", pick[2]),
      sd_z = row_of("sd_z", pick[2])
    )
    if (fault == "rare_pop") {
      # Cluster checks need complete data; only this fault is scored with them.
      cs <- sw_cluster_scores(files, ch, CLUSTER_K, CLUSTER_EVENTS, seed, 3, 30)
      d <- cs[cs$check == "median_z" & cs$marker == pick[2], ]
      scores$cluster_median_z <- as.numeric(tapply(abs(d$value), d$file, top)[names(files)])
    }
    tibble(
      seed = seed,
      fault = fault,
      check = names(scores),
      auroc = vapply(scores, \(v) auroc(v, bad), numeric(1))
    )
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(long, file.path(run_dir, "file_faults_long.csv"), row.names = FALSE)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(fault, check)
  )
write.csv(mat, file.path(run_dir, "figure_source_data_file_faults_real.csv"), row.names = FALSE)
print(as.data.frame(slice_max(mutate(mat, auroc = round(auroc, 2)), auroc, n = 2, by = fault)))
p <- ggplot(mat, aes(auroc, check)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)), size = 0.15) +
  facet_wrap(vars(fault), nrow = 2) +
  labs(x = "AUROC on real files (mean, 95% CI over 5 seeds)", y = NULL)
ggsave(file.path(run_dir, "fig_file_faults_real.png"), p, width = 220, height = 110, units = "mm", dpi = 200)
