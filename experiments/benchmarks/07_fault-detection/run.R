# experiments/benchmarks/07_fault-detection/run.R
# Purpose: detection matrix. Inject one fault class at a time into simulated
#          cohorts and score how well each diagnostic tool ranks the faulty
#          files (or markers) above the clean ones. AUROC per fault x tool.
# Run:     Rscript experiments/benchmarks/07_fault-detection/run.R <out_dir>
# Inputs:  none (simulated here)
# Outputs: <out_dir>/detection_long.csv   per seed, fault, tool AUROC
#          <out_dir>/detection_matrix.csv mean AUROC and 95% half-width
#
# Every value below is benchmark-only and illustrative. None is a package default.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
for (f in c(
  "corrector",
  "anomaly_features",
  "anomaly_flag",
  "fcs_signal_audit",
  "batch_scores"
)) {
  source(file.path("R", paste0(f, ".R")))
}

out_dir <- commandArgs(trailingOnly = TRUE)[[1]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Benchmark settings ---------------------------------------------------------
SEEDS <- 1:10
N_CLEAN <- 60 # clean files in a file-level cohort
N_BAD <- 12 # faulty files in a file-level cohort
N_CH <- 8 # markers
N_EV <- 1000 # events per file
FILE_JITTER <- 0.1 # SD of file-to-file population-mean jitter
# Instrument maximum on the simulated scale. It sits above every clean, stretched
# or debris value (debris is 8), so only the saturation fault reaches it. At 6 the
# scale, stretch and debris faults crossed it and scored 1.00 by construction.
RANGE_MAX <- 9
N_BATCH <- 4 # batches in a batch-level cohort
FILES_PER_BATCH <- 12
# Feature and flag settings (as in benchmark v1, not tuned):
QSTEP <- 0.1
NEV <- 500
AGG <- 5000
NTREES <- 300
FOREST_SAMPLE <- 32
CHANNELS <- paste0("C", seq_len(N_CH))

# Simulation -----------------------------------------------------------------
make_file <- function(mu) {
  pop <- sample(1:3, N_EV, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  jit <- matrix(stats::rnorm(3 * N_CH, 0, FILE_JITTER), 3)
  x <- mu[pop, ] + jit[pop, ] + matrix(stats::rnorm(N_EV * N_CH), N_EV)
  colnames(x) <- CHANNELS
  cbind(Time = sort(stats::runif(N_EV, 0, 100)), x)
}
MU <- matrix(rep(c(-2, 0, 2), N_CH), nrow = 3)

# File-level faults: each takes a file matrix and returns the faulty one.
FILE_FAULTS <- list(
  scale_all = function(x) {
    ch <- CHANNELS
    mu <- colMeans(x[, ch])
    x[, ch] <- sweep(sweep(x[, ch], 2, mu), 2, 1.3, "*") +
      rep(mu, each = nrow(x))
    x
  },
  shift = function(x) {
    ch <- sample(CHANNELS, 1)
    x[, ch] <- x[, ch] + 0.5
    x
  },
  stretch = function(x) {
    ch <- sample(CHANNELS, 1)
    x[, ch] <- mean(x[, ch]) + (x[, ch] - mean(x[, ch])) * 1.5
    x
  },
  bimodal = function(x) {
    ch <- sample(CHANNELS, 1)
    x[, ch] <- x[, ch] + sample(c(-0.5, 0.5), nrow(x), replace = TRUE)
    x
  },
  debris = function(x) {
    ch <- sample(CHANNELS, 1)
    x[sample(nrow(x), 0.1 * nrow(x)), ch] <- 8
    x
  },
  saturation = function(x) {
    ch <- sample(CHANNELS, 1)
    x[sample(nrow(x), 0.15 * nrow(x)), ch] <- RANGE_MAX
    x
  },
  dead_channel = function(x) {
    x[, sample(CHANNELS, 1)] <- 0
    x
  },
  nonfinite = function(x) {
    ch <- sample(CHANNELS, 1)
    x[sample(nrow(x), 0.05 * nrow(x)), ch] <- NA
    x
  },
  clog = function(x) {
    # Events in a time window are thinned by 80%: the event rate drops.
    win <- x[, "Time"] > 40 & x[, "Time"] < 70
    drop <- win & stats::runif(nrow(x)) < 0.8
    x[!drop, ]
  },
  spillover = function(x) {
    ch <- sample(CHANNELS, 2)
    x[, ch[2]] <- 0.6 * x[, ch[1]] + 0.8 * x[, ch[2]]
    x
  },
  # Added after the first run so earlier faults keep their seeds.
  time_drift = function(x) {
    ch <- sample(CHANNELS, 1)
    x[, ch] <- x[, ch] + 0.03 * x[, "Time"] # about 3 SD over the run
    x
  }
)

auroc <- function(score, pos) {
  ok <- !is.na(score)
  score <- score[ok]
  pos <- pos[ok]
  if (sum(pos) == 0 || sum(!pos) == 0) {
    return(NA_real_)
  }
  r <- rank(score)
  n1 <- sum(pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * sum(!pos))
}

# Per-file scores from every file-level tool.
file_tool_scores <- function(files, seed) {
  ids <- names(files)
  with_na <- vapply(files, \(m) anyNA(m[, CHANNELS]), logical(1))
  # Forest score needs complete data; files with NA get no forest score.
  forest <- rep(NA_real_, length(ids))
  if (!any(with_na)) {
    feat <- sw_anomaly_features(
      files,
      CHANNELS,
      c("quantile", "emd"),
      QSTEP,
      NEV,
      seed,
      "simulated",
      "all",
      AGG
    )
    fl <- sw_anomaly_flag_outlier(feat, NTREES, FOREST_SAMPLE, 1, 0.5, seed)
    forest <- fl$score[match(ids, fl$file)]
  }
  aud <- sw_fcs_signal_audit(
    files,
    CHANNELS,
    range_max = RANGE_MAX,
    range_min = NULL,
    time_channel = "Time",
    time_bin = 5
  )
  by <- function(check) {
    a <- aud[aud$check == check, ]
    tapply(a$value, a$file, max, na.rm = TRUE)[ids]
  }
  # Signed checks (time trend, z-scores): the size of the deviation is what counts.
  by_abs <- function(check) {
    a <- aud[aud$check == check, ]
    tapply(abs(a$value), a$file, max, na.rm = TRUE)[ids]
  }
  tibble(
    file = ids,
    forest = forest,
    saturation_share = as.numeric(by("saturation_share")),
    zero_variance = as.numeric(by("zero_variance")),
    nonfinite_share = as.numeric(by("nonfinite_share")),
    rate_cv = as.numeric(by("rate_cv")),
    max_gap = as.numeric(by("max_gap")),
    corr_shift = as.numeric(by("corr_shift")),
    time_trend = as.numeric(by_abs("time_trend")),
    median_z = as.numeric(by_abs("median_z")),
    sd_z = as.numeric(by_abs("sd_z"))
  )
}

file_level <- function(seed) {
  lapply(names(FILE_FAULTS), function(fault) {
    set.seed(seed * 100 + match(fault, names(FILE_FAULTS)))
    files <- lapply(seq_len(N_CLEAN + N_BAD), \(i) make_file(MU))
    bad <- c(rep(FALSE, N_CLEAN), rep(TRUE, N_BAD))
    files[bad] <- lapply(files[bad], FILE_FAULTS[[fault]])
    ord <- sample(length(files))
    files <- files[ord]
    bad <- bad[ord]
    names(files) <- paste0("f", seq_along(files))
    sc <- file_tool_scores(files, seed)
    tool_cols <- setdiff(names(sc), "file")
    tibble(
      seed = seed,
      fault = fault,
      level = "file",
      tool = tool_cols,
      auroc = vapply(tool_cols, \(t) auroc(sc[[t]], bad), numeric(1))
    )
  }) |>
    bind_rows()
}

# Batch-level faults: one marker in one batch is altered; score the markers.
BATCH_FAULTS <- c(
  "batch_shift",
  "batch_stretch",
  "batch_bimodal",
  "batch_drift"
)
batch_level <- function(seed) {
  lapply(BATCH_FAULTS, function(fault) {
    set.seed(seed * 100 + 50 + match(fault, BATCH_FAULTS))
    n <- N_BATCH * FILES_PER_BATCH
    batch <- rep(paste0("b", seq_len(N_BATCH)), each = FILES_PER_BATCH)
    files <- lapply(seq_len(n), \(i) make_file(MU)[, CHANNELS])
    names(files) <- paste0("f", seq_len(n))
    names(batch) <- names(files)
    ord <- stats::setNames(seq_len(n), names(files)) # acquisition order
    mk <- sample(CHANNELS, 1) # the faulty marker
    hit <- which(batch == "b2")
    for (k in seq_along(hit)) {
      x <- files[[hit[k]]]
      x[, mk] <- switch(
        fault,
        batch_shift = x[, mk] + 0.5,
        batch_stretch = mean(x[, mk]) + (x[, mk] - mean(x[, mk])) * 1.5,
        batch_bimodal = x[, mk] + sample(c(-0.5, 0.5), nrow(x), replace = TRUE),
        batch_drift = x[, mk] + 0.1 * k # steady drift over the batch
      )
      files[[hit[k]]] <- x
    }
    sc <- sw_batch_scores(files, batch, CHANNELS, NEV, seed, order = ord)
    pe <- sc[sc$check == "pair_emd", ] |>
      summarise(v = max(value), .by = marker)
    mr <- sc[sc$check == "median_range", ]
    e2 <- sc[sc$check == "eta2", ]
    rho <- sc[sc$check == "order_rho", ] |>
      summarise(v = max(abs(value)), .by = marker)
    per_marker <- list(
      pair_emd = setNames(pe$v, pe$marker),
      median_range = setNames(mr$value, mr$marker),
      eta2 = setNames(e2$value, e2$marker),
      order_rho = setNames(rho$v, rho$marker)
    )
    tibble(
      seed = seed,
      fault = fault,
      level = "batch",
      tool = names(per_marker),
      auroc = vapply(
        per_marker,
        \(v) auroc(v[CHANNELS], CHANNELS == mk),
        numeric(1)
      )
    )
  }) |>
    bind_rows()
}

# Null check: no fault injected, tools scored against shuffled labels.
null_check <- function(seed) {
  set.seed(seed * 100 + 99)
  files <- lapply(seq_len(N_CLEAN + N_BAD), \(i) make_file(MU))
  names(files) <- paste0("f", seq_along(files))
  sc <- file_tool_scores(files, seed)
  bad <- sample(rep(c(FALSE, TRUE), c(N_CLEAN, N_BAD)))
  tool_cols <- setdiff(names(sc), "file")
  tibble(
    seed = seed,
    fault = "none (shuffled labels)",
    level = "file",
    tool = tool_cols,
    auroc = vapply(tool_cols, \(t) auroc(sc[[t]], bad), numeric(1))
  )
}

long <- lapply(SEEDS, \(s) {
  bind_rows(file_level(s), batch_level(s), null_check(s))
}) |>
  bind_rows()
write.csv(long, file.path(out_dir, "detection_long.csv"), row.names = FALSE)

mat <- long |>
  summarise(
    auroc = mean(auroc, na.rm = TRUE),
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    n_seeds = sum(!is.na(auroc)),
    .by = c(level, fault, tool)
  )
write.csv(mat, file.path(out_dir, "detection_matrix.csv"), row.names = FALSE)
cat(
  "wrote",
  nrow(long),
  "rows and",
  nrow(mat),
  "matrix cells to",
  out_dir,
  "\n"
)
