# experiments/benchmarks/08_detection-limits/run.R
# Purpose: find how small a fault can be before each tool stops ranking the faulty files
#          above the clean ones. One fault family at a time, effect size on a grid.
# Run:     Rscript experiments/benchmarks/08_detection-limits/run.R <out_dir>
# Inputs:  none (simulated here)
# Outputs: <out_dir>/limits_long.csv   per seed, sweep, level, tool AUROC
#          <out_dir>/limits_matrix.csv mean AUROC and 95% half-width
#
# Every value below is benchmark-only and illustrative. None is a package default.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
for (f in c("corrector", "anomaly_features", "anomaly_flag", "fcs_signal_audit", "cluster_scores")) {
  source(file.path("R", paste0(f, ".R")))
}

out_dir <- commandArgs(trailingOnly = TRUE)[[1]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Benchmark settings ---------------------------------------------------------
SEEDS <- 1:6
N_CLEAN <- 60
N_BAD <- 12
N_CH <- 8
N_EV <- 1000
FILE_JITTER <- 0.1
QSTEP <- 0.1
NEV <- 500
AGG <- 5000
NTREES <- 300
FOREST_SAMPLE <- 32
CLUSTER_K <- 6 # clusters for the cluster-level scores
SHIFT_SIZES <- c(0.05, 0.1, 0.2, 0.3, 0.5) # SD units added to one marker
STRETCH_SIZES <- c(1.05, 1.1, 1.2, 1.3, 1.5) # width factor on one marker
RARE_FRACTIONS <- c(0.02, 0.05, 0.1, 0.2) # share of events in the shifted population
RARE_SHIFT <- 1 # SD units added in that population
CHANNELS <- paste0("C", seq_len(N_CH))

make_file <- function(rare_fraction = 0.2) {
  p <- c(0.5, 0.5 - rare_fraction, rare_fraction)
  pop <- sample(1:3, N_EV, replace = TRUE, prob = p)
  mu <- matrix(rep(c(-2, 0, 2), N_CH), nrow = 3)
  jit <- matrix(stats::rnorm(3 * N_CH, 0, FILE_JITTER), 3)
  x <- mu[pop, ] + jit[pop, ] + matrix(stats::rnorm(N_EV * N_CH), N_EV)
  colnames(x) <- CHANNELS
  out <- cbind(Time = sort(stats::runif(N_EV, 0, 100)), x)
  attr(out, "pop") <- pop
  out
}

auroc <- function(score, pos) {
  r <- rank(score)
  n1 <- sum(pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * sum(!pos))
}

# Scores per file for four tools: forest, median_z, sd_z, cluster median_z.
tool_scores <- function(files, seed) {
  ids <- names(files)
  feat <- sw_anomaly_features(files, CHANNELS, c("quantile", "emd"), QSTEP, NEV, seed, "simulated", "all", AGG)
  fl <- sw_anomaly_flag_outlier(feat, NTREES, FOREST_SAMPLE, 1, 0.5, seed)
  aud <- sw_fcs_signal_audit(files, CHANNELS, range_max = NULL, range_min = NULL, time_channel = "Time", time_bin = 5)
  by_abs <- function(check) {
    a <- aud[aud$check == check, ]
    as.numeric(tapply(abs(a$value), a$file, max, na.rm = TRUE)[ids])
  }
  cs <- sw_cluster_scores(files, CHANNELS, CLUSTER_K, NEV, seed, 3, 30)
  cm <- cs[cs$check == "median_z", ]
  cm <- vapply(split(abs(cm$value), cm$file)[ids], max, numeric(1), na.rm = TRUE)
  tibble(
    forest = fl$score[match(ids, fl$file)],
    median_z = by_abs("median_z"),
    sd_z = by_abs("sd_z"),
    cluster_median_z = as.numeric(cm)
  )
}

run_one <- function(seed, sweep, size, fault, rare_fraction = 0.2) {
  set.seed(seed * 1000 + round(size * 100))
  files <- lapply(seq_len(N_CLEAN + N_BAD), \(i) make_file(rare_fraction))
  bad <- c(rep(FALSE, N_CLEAN), rep(TRUE, N_BAD))
  files[bad] <- lapply(files[bad], fault)
  ord <- sample(length(files))
  files <- files[ord]
  bad <- bad[ord]
  names(files) <- paste0("f", seq_along(files))
  sc <- tool_scores(files, seed)
  tibble(
    seed = seed,
    sweep = sweep,
    size = size,
    tool = names(sc),
    auroc = vapply(sc, \(v) auroc(v, bad), numeric(1))
  )
}

jobs <- c(
  lapply(SHIFT_SIZES, \(s) list(sweep = "shift", size = s, rf = 0.2, fault = \(x) {
    ch <- sample(CHANNELS, 1)
    x[, ch] <- x[, ch] + s
    x
  })),
  lapply(STRETCH_SIZES, \(s) list(sweep = "stretch", size = s, rf = 0.2, fault = \(x) {
    ch <- sample(CHANNELS, 1)
    x[, ch] <- mean(x[, ch]) + (x[, ch] - mean(x[, ch])) * s
    x
  })),
  lapply(RARE_FRACTIONS, \(f) list(sweep = "rare_pop_fraction", size = f, rf = f, fault = \(x) {
    ch <- sample(CHANNELS, 1)
    x[attr(x, "pop") == 3, ch] <- x[attr(x, "pop") == 3, ch] + RARE_SHIFT
    x
  }))
)
long <- lapply(SEEDS, \(sd) {
  lapply(jobs, \(j) run_one(sd, j$sweep, j$size, j$fault, j$rf))
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(long, file.path(out_dir, "limits_long.csv"), row.names = FALSE)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc) / sqrt(n()),
    auroc = mean(auroc),
    .by = c(sweep, size, tool)
  )
write.csv(mat, file.path(out_dir, "limits_matrix.csv"), row.names = FALSE)
cat("wrote", nrow(long), "rows to", out_dir, "\n")
