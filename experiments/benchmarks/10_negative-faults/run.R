# experiments/benchmarks/10_negative-faults/run.R
# Purpose: do the negative-population checks (sw_negative_audit) find faults that live only in
#          the negative population, and do they beat the whole-file checks (median_z, sd_z,
#          corr_shift) on those faults? Fault classes from the Colibri blog: shifted negative
#          (unmixing skew), widened negative (spread / AF mismatch), hyper-negative events
#          (AF over-correction), angled negative (spillover between two channels), and a leak
#          from a bright channel into a negative one (tandem breakdown, dye interaction).
# Run:     Rscript experiments/benchmarks/10_negative-faults/run.R <out_dir>
# Inputs:  none (simulated here, analysis scale)
# Outputs: <out_dir>/negative_long.csv    per seed, fault, size, check, view: AUROC
#
# Every value below is benchmark-only and illustrative. None is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
})

out_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run dir must be new" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)

# Benchmark settings ---------------------------------------------------------
SEEDS <- 1:5
N_CLEAN <- 60
N_BAD <- 12
N_EV <- 2000
CHANNELS <- paste0("C", 1:6)
NEG_SD <- 0.5 # SD of the negative population
POS_MEAN <- 3
NEG_EVENTS <- 1000 # events per file for the audit
HYPER_K <- 4 # hyper-negative cut: negative median minus HYPER_K robust SDs
FAULTS <- tibble::tribble(
  ~fault, ~size,
  "neg_shift", 0.1, # added to the negative population, in units of NEG_SD
  "neg_shift", 0.3,
  "neg_widen", 1.2, # factor on the negative SD
  "neg_widen", 1.5,
  "hyper_negative", 0.01, # share of negative events moved to -8 negative SDs
  "hyper_negative", 0.03,
  "angled_negative", 0.1, # slope added from channel 1 into channel 2, negatives only
  "angled_negative", 0.3,
  "spill_leak", 0.05, # share of the C1 signal leaking into C2 where C1 is positive and C2 negative
  "spill_leak", 0.15
)

# Simulation -----------------------------------------------------------------
make_file <- function() {
  x <- sapply(CHANNELS, function(ch) {
    neg <- stats::runif(N_EV) < stats::runif(1, 0.5, 0.7)
    ifelse(neg, stats::rnorm(N_EV, 0, NEG_SD), stats::rnorm(N_EV, POS_MEAN, 0.8))
  })
  colnames(x) <- CHANNELS
  x
}
# Faults act on the negative population of the hit channel (C2; C1 also for the angle).
inject <- function(x, fault, size) {
  neg2 <- x[, "C2"] < POS_MEAN / 2
  neg1 <- x[, "C1"] < POS_MEAN / 2
  switch(
    fault,
    neg_shift = {
      x[neg2, "C2"] <- x[neg2, "C2"] + size * NEG_SD
    },
    neg_widen = {
      x[neg2, "C2"] <- x[neg2, "C2"] * size
    },
    hyper_negative = {
      idx <- which(neg2)[seq_len(round(size * sum(neg2)))]
      x[idx, "C2"] <- -8 * NEG_SD
    },
    angled_negative = {
      both <- neg1 & neg2
      x[both, "C2"] <- x[both, "C2"] + size * x[both, "C1"]
    },
    spill_leak = {
      leak <- !neg1 & neg2
      x[leak, "C2"] <- x[leak, "C2"] + size * x[leak, "C1"]
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

one_run <- function(seed, fault, size) {
  set.seed(seed * 1000 + match(fault, unique(FAULTS$fault)) * 10 + round(size * 100))
  files <- lapply(seq_len(N_CLEAN + N_BAD), \(i) make_file())
  bad <- c(rep(FALSE, N_CLEAN), rep(TRUE, N_BAD))
  files[bad] <- lapply(files[bad], inject, fault = fault, size = size)
  ord <- sample(length(files))
  files <- files[ord]
  bad <- bad[ord]
  ids <- paste0("f", seq_along(files))
  names(files) <- ids

  na <- sw_negative_audit(files, CHANNELS, NEG_EVENTS, seed, HYPER_K, pairs = "all")
  sa <- sw_fcs_signal_audit(files, CHANNELS, NULL, NULL, time_channel = NULL, time_bin = NULL)
  hit_row <- function(tab, check, item, file_col, item_col) {
    d <- tab[tab$check == check & tab[[item_col]] %in% item, ]
    as.numeric(tapply(abs(d$value), d[[file_col]], max, na.rm = TRUE)[ids])
  }
  max_row <- function(tab, check, file_col) {
    d <- tab[tab$check == check, ]
    as.numeric(tapply(abs(d$value), d[[file_col]], max, na.rm = TRUE)[ids])
  }
  hit <- if (fault == "angled_negative") "C1 ~ C2" else "C2"
  scores <- list(
    "neg_median_z | hit row" = hit_row(na, "neg_median_z", "C2", "file", "channel"),
    "neg_sd_z | hit row" = hit_row(na, "neg_sd_z", "C2", "file", "channel"),
    "hyper_negative_share | hit row" = hit_row(na, "hyper_negative_share", "C2", "file", "channel"),
    "neg_corr_z | hit row" = hit_row(na, "neg_corr_z", "C1 ~ C2", "file", "channel"),
    "neg_corr_z | max over pairs" = max_row(na, "neg_corr_z", "file"),
    "spill_corr_z | hit row" = hit_row(na, "spill_corr_z", "C1 > C2", "file", "channel"),
    "spill_corr_z | max over pairs" = max_row(na, "spill_corr_z", "file"),
    "median_z (whole file) | hit row" = hit_row(sa, "median_z", "C2", "file", "channel"),
    "sd_z (whole file) | hit row" = hit_row(sa, "sd_z", "C2", "file", "channel"),
    "corr_shift (whole file) | max over pairs" = max_row(sa, "corr_shift", "file")
  )
  tibble(
    seed = seed,
    fault = fault,
    size = size,
    score = names(scores),
    auroc = vapply(scores, \(v) auroc(v, bad), numeric(1))
  )
}

null_run <- function(seed) {
  # No fault; labels shuffled. Every check should sit near 0.5.
  r <- one_run(seed, "neg_shift", 0)
  r$fault <- "none (shuffled)"
  r
}

long <- bind_rows(
  lapply(SEEDS, \(s) bind_rows(lapply(seq_len(nrow(FAULTS)), \(i) one_run(s, FAULTS$fault[[i]], FAULTS$size[[i]])))),
  lapply(SEEDS, null_run)
)
write.csv(long, file.path(out_dir, "negative_long.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("wrote", nrow(long), "rows\n")
