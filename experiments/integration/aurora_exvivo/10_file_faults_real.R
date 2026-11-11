# experiments/integration/aurora_exvivo/10_file_faults_real.R
# Purpose: inject the file-level fault classes of benchmark 07 into 25% of the real donor files of
#          one batch, sized in units of each marker's real SD, and score the checks that target
#          them. Fills the "not tested" rows of the coverage table in reference/diagnostics.qmd.
# G40:     null arm (no injection, same draws)
#          arm = "null" rows use the same bad files and hit channels as the fault arm, score the
#          clean base files, and give the AUROC of every check against those same bad labels.
# G41:     SW_SEEDS, SW_NULL_ONLY, SW_ONLY_CASE (defaults keep G40 behaviour)
# Run:     Rscript experiments/integration/aurora_exvivo/10_file_faults_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/file_faults_long.csv                  seed, fault, arm (fault/null), check, auroc
#          <run_dir>/figure_source_data_file_faults_real.csv  per fault, arm, check: mean AUROC, 95% CI
#          <run_dir>/fig_file_faults_real.png               fault arm only
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
# G41 switches (defaults keep G40 behaviour): seeds, null arm only, restricted grid.
SEEDS <- eval(parse(text = Sys.getenv("SW_SEEDS", "1:5")))
NULL_ONLY <- Sys.getenv("SW_NULL_ONLY") == "1" # injection draws still run; fault scoring is skipped
ONLY_CASE <- Sys.getenv("SW_ONLY_CASE") # comma list of fault names; empty = all
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

FAULTS_RUN <- if (nzchar(ONLY_CASE)) {
  strsplit(ONLY_CASE, ",")[[1]]
} else {
  FAULT_NAMES
}
stopifnot("SW_ONLY_CASE names known faults" = all(FAULTS_RUN %in% FAULT_NAMES))

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
time_bin <- as.numeric(
  median(vapply(base, \(m) diff(range(m[, "Time"])), numeric(1))) / TIME_BINS
)
# Clusters for the rare-population fault, from the clean pool (seeded).
set.seed(1)
km <- stats::kmeans(
  pooled[sample.int(nrow(pooled), 20000), ch],
  CLUSTER_K,
  nstart = 3,
  iter.max = 100
)
smallest <- which.min(km$size)
nearest <- function(x) {
  d2 <- outer(rowSums(x^2), rowSums(km$centers^2), "+") -
    2 * tcrossprod(x, km$centers)
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
      x[idx, ch] <- matrix(
        cohort_low[ch],
        length(idx),
        length(ch),
        byrow = TRUE
      ) +
        stats::rnorm(length(idx) * length(ch), 0, 0.05 * marker_sd[ch])
    },
    bimodal = {
      x[, b] <- x[, b] +
        sample(c(-0.5, 0.5), nrow(x), replace = TRUE) * marker_sd[[b]]
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
  lapply(FAULTS_RUN, function(fault) {
    set.seed(seed * 1000 + match(fault, FAULT_NAMES))
    bad <- seq_along(base) %in%
      sample(seq_along(base), round(FRACTION_BAD * length(base)))
    pick <- sample(ch, 2)
    files <- base
    files[bad] <- lapply(
      files[bad],
      inject,
      fault = fault,
      a = pick[1],
      b = pick[2]
    )
    # Scores every check on one set of files against the bad labels of this iteration.
    score_files <- function(files) {
      sa <- sw_fcs_signal_audit(
        files,
        ch,
        range_max = cohort_max,
        range_min = NULL,
        time_channel = "Time",
        time_bin = time_bin,
        area_height = NULL,
        area_height_quantile = NULL,
        reference_group = NULL,
        reference = NULL
      )
      # Hit channel's row for channel checks; the file's single row for file-level checks.
      row_of <- function(check, item) {
        d <- sa[
          sa$check == check & (is.na(sa$channel) | sa$channel %in% item),
        ]
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
        # The injected pair's row (oracle, like the channel rows), and the largest over all pairs.
        # Pair names follow the channel order of `ch`.
        corr_pair_z = row_of(
          "corr_pair_z",
          paste(pick[order(match(pick, ch))], collapse = " ~ ")
        ),
        corr_pair_z_max = row_of(
          "corr_pair_z",
          unique(sa$channel[sa$check == "corr_pair_z"])
        ),
        median_z = row_of("median_z", pick[2]),
        sd_z = row_of("sd_z", pick[2])
      )
      if (fault == "rare_pop") {
        # Cluster checks need complete data; only this fault is scored with them.
        cs <- sw_cluster_scores(
          files,
          ch,
          CLUSTER_K,
          CLUSTER_EVENTS,
          seed,
          3,
          30,
          reference = NULL,
          reference_group = NULL
        )
        d <- cs[cs$check == "median_z" & cs$marker == pick[2], ]
        scores$cluster_median_z <- as.numeric(tapply(
          abs(d$value),
          d$file,
          top
        )[names(files)])
      }
      tibble(
        seed = seed,
        fault = fault,
        check = names(scores),
        auroc = vapply(scores, \(v) auroc(v, bad), numeric(1))
      )
    }
    # Fault arm first: it uses the random stream exactly as before. The null arm comes last, so
    # it cannot change any draw of this iteration; the next iteration starts with its own set.seed().
    # The injection above and its draws always run; only the fault scoring is skipped (the null
    # arm needs nothing from it: scoring reads no global RNG and sets no state).
    fault_rows <- if (NULL_ONLY) {
      NULL
    } else {
      mutate(score_files(files), arm = "fault", .after = fault)
    }
    null_rows <- mutate(score_files(base), arm = "null", .after = fault)
    bind_rows(fault_rows, null_rows)
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
    .by = c(fault, arm, check)
  )
write.csv(
  mat,
  file.path(run_dir, "figure_source_data_file_faults_real.csv"),
  row.names = FALSE
)
if (!NULL_ONLY) {
  print(as.data.frame(slice_max(
    mutate(filter(mat, arm == "fault"), auroc = round(auroc, 2)),
    auroc,
    n = 2,
    by = fault
  )))
}
if (!NULL_ONLY) {
  p <- ggplot(filter(mat, arm == "fault"), aes(auroc, check)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
    geom_pointrange(
      aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)),
      size = 0.15
    ) +
    facet_wrap(vars(fault), nrow = 2) +
    labs(x = "AUROC on real files (mean, 95% CI over 5 seeds)", y = NULL)
  ggsave(
    file.path(run_dir, "fig_file_faults_real.png"),
    p,
    width = 220,
    height = 110,
    units = "mm",
    dpi = 200
  )
}
