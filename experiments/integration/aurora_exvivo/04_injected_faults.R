# experiments/integration/aurora_exvivo/04_injected_faults.R
# Purpose: score the diagnostics on real files with injected faults. The background is real
#          (one batch of donor files, so no batch effect is in the cohort); the fault is
#          added to a random subset. Ground truth is known, so AUROC can be computed.
# Run:     Rscript experiments/integration/aurora_exvivo/04_injected_faults.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/injected_long.csv    AUROC per seed, fault, size and tool
#          <run_dir>/injected_matrix.csv  mean AUROC and 95% half-width
#          <run_dir>/fig_injected_faults.png
#          <run_dir>/figure_source_data_injected.csv
#
# Every value below is illustrative; none is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
BATCH <- "Batch1" # cohort: donor files of this batch
SEEDS <- 1:5
FRACTION_BAD <- 0.25 # share of files that get the fault
SHIFT_SD <- c(0.3, 0.5, 1) # added to one marker, in units of that marker's SD in the cohort
STRETCH <- c(1.1, 1.3, 1.5) # width factor on one marker
QSTEP <- 0.1
NEV <- 2000
AGG <- 20000
NTREES <- 300
FOREST_SAMPLE <- 32
CLUSTER_K <- 12
TIME_BINS <- 50

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
keep <- info$file[info$batch == BATCH & grepl("BMV_[0-9]+[A-Z]", info$file)]
base <- ev[keep]
ch <- setdiff(colnames(base[[1]]), "Time")
marker_sd <- apply(do.call(rbind, base)[, ch], 2, sd)
time_bin <- as.numeric(median(vapply(base, \(m) diff(range(m[, "Time"])), numeric(1))) / TIME_BINS)
cat("cohort files:", length(base), "\n")

auroc <- function(score, pos) {
  r <- rank(score)
  n1 <- sum(pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * sum(!pos))
}

tool_scores <- function(files, seed) {
  ids <- names(files)
  feat <- sw_anomaly_features(files, ch, c("quantile", "emd"), QSTEP, NEV, seed, "real, OMIQ-transformed", "all", AGG)
  fl <- sw_anomaly_flag_outlier(feat, NTREES, FOREST_SAMPLE, 1, 0.5, seed)
  aud <- sw_fcs_signal_audit(files, ch, range_max = NULL, range_min = NULL, time_channel = "Time", time_bin = time_bin)
  by_abs <- function(check) {
    a <- aud[aud$check == check, ]
    as.numeric(tapply(abs(a$value), a$file, max, na.rm = TRUE)[ids])
  }
  cs <- sw_cluster_scores(files, ch, CLUSTER_K, NEV, seed, 3, 30, reference = NULL)
  cm <- cs[cs$check == "median_z", ]
  cm <- vapply(split(abs(cm$value), cm$file)[ids], max, numeric(1), na.rm = TRUE)
  tibble(
    forest = fl$score[match(ids, fl$file)],
    median_z = by_abs("median_z"),
    sd_z = by_abs("sd_z"),
    cluster_median_z = as.numeric(cm)
  )
}

jobs <- c(
  lapply(SHIFT_SD, \(s) list(fault = "shift", size = s, fn = \(x, m) {
    x[, m] <- x[, m] + s * marker_sd[[m]]
    x
  })),
  lapply(STRETCH, \(s) list(fault = "stretch", size = s, fn = \(x, m) {
    x[, m] <- mean(x[, m]) + (x[, m] - mean(x[, m])) * s
    x
  }))
)

long <- lapply(SEEDS, \(sd) {
  lapply(jobs, \(j) {
    set.seed(sd * 1000 + round(j$size * 100))
    n_bad <- round(FRACTION_BAD * length(base))
    bad <- seq_along(base) %in% sample(seq_along(base), n_bad)
    marker <- sample(ch, 1)
    files <- base
    files[bad] <- lapply(files[bad], j$fn, m = marker)
    sc <- tool_scores(files, sd)
    tibble(
      seed = sd,
      fault = j$fault,
      size = j$size,
      tool = names(sc),
      auroc = vapply(sc, \(v) auroc(v, bad), numeric(1))
    )
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(long, file.path(run_dir, "injected_long.csv"), row.names = FALSE)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc) / sqrt(n()),
    auroc = mean(auroc),
    .by = c(fault, size, tool)
  )
write.csv(mat, file.path(run_dir, "injected_matrix.csv"), row.names = FALSE)
write.csv(mat, file.path(run_dir, "figure_source_data_injected.csv"), row.names = FALSE)

p <- ggplot(mat, aes(size, auroc, colour = tool)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = auroc - ci95_half, ymax = auroc + ci95_half), size = 0.25) +
  geom_line() +
  facet_wrap(vars(fault), scales = "free_x", nrow = 1) +
  labs(x = "Effect size (SD added, or width factor)", y = "AUROC (mean, 95% CI over seeds)", colour = NULL)
ggsave(file.path(run_dir, "fig_injected_faults.png"), p, width = 140, height = 70, units = "mm", dpi = 200)
cat("wrote", nrow(long), "rows\n")
