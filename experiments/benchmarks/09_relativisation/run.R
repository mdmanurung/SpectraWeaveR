# experiments/benchmarks/09_relativisation/run.R
# Purpose: check that CD3 relativisation (Glehr et al. 2026, Box 1) removes a per-instrument
#          offset and gain, and measure how much it degrades when the CD3 level differs between
#          donors (the method's assumption A1 broken on purpose).
# Run:     Rscript experiments/benchmarks/09_relativisation/run.R <out_dir>
# Inputs:  none (simulated here, linear scale)
# Outputs: <out_dir>/relativisation_long.csv  per seed, CD3 donor SD, data version, channel:
#                                             largest pairwise batch EMD (arcsinh scale)
#          <out_dir>/otd_long.csv             per seed, CD3 donor SD, version, file: OTD to
#                                             the pooled reference (Glehr et al. score)
#          <out_dir>/gain_recovery.csv        per seed and CD3 donor SD: correlation of the
#                                             fitted batch gain ratio with the true gain
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
N_BATCH <- 4 # instruments
DONORS_PER_BATCH <- 8
CHANNELS <- paste0("ch", 1:6)
N_EVENTS <- 4000 # events per sample file
N_CONTROL <- 2000 # events per CD3 control
CD3_LEVEL <- 1e4 # true CD3+ median on the linear scale, so truth and relativised share a scale
NEG_LEVEL <- 150 # true median of reference-negative (autofluorescence) events
GAIN_SD <- 0.4 # SD of log instrument gain per batch and channel
OFFSET_MAX <- 300 # instrument offsets drawn from U(0, OFFSET_MAX)
CD3_DONOR_SD <- c(0, 0.1, 0.3) # SD of log CD3 level between donors (0 = assumption A1 holds)
COFACTOR <- 6000 # arcsinh cofactor for scoring and for the k-means split
KM_NSTART <- 5
MIN_EVENTS <- 100
SCORE_EVENTS <- 1000 # events per file for the EMD scores
OTD_EVENTS <- 300 # events per file for the OTD (smaller than the paper's 10,000 for speed)
OTD_REFERENCE <- 1280 # pooled reference events (40 per file)
OTD_BLUR <- 0.75 # paper's blur, arcsinh units
OTD_MAX_ITER <- 500
OTD_TOL <- 1e-5

# Simulation -----------------------------------------------------------------
true_sample <- function() {
  # Each channel: half the events negative (median 200), half positive (median 3000), with a
  # donor-specific positive share.
  x <- vapply(CHANNELS, \(ch) {
    share <- stats::runif(1, 0.3, 0.7)
    pos <- stats::runif(N_EVENTS) < share
    ifelse(pos, 3000, 200) * exp(stats::rnorm(N_EVENTS, 0, 0.4))
  }, numeric(N_EVENTS))
  colnames(x) <- CHANNELS
  x
}
true_control <- function(channel, cd3_factor) {
  # CD3 single stain in this channel: half CD3- (autofluorescence), half CD3+.
  v <- c(
    NEG_LEVEL * exp(stats::rnorm(N_CONTROL / 2, 0, 0.3)),
    CD3_LEVEL * cd3_factor * exp(stats::rnorm(N_CONTROL / 2, 0, 0.3))
  )
  m <- matrix(0, length(v), length(CHANNELS), dimnames = list(NULL, CHANNELS))
  m[, channel] <- v
  m
}
instrument <- function(x, offset, gain) {
  # Linear detector response, per channel; small additive noise.
  sweep(sweep(x, 2, gain, "*"), 2, offset, "+") + stats::rnorm(length(x), 0, 20)
}
arcsinh <- function(files) lapply(files, \(m) asinh(m / COFACTOR))
max_pair_emd <- function(files, batch, seed) {
  sc <- sw_batch_scores(files, batch, CHANNELS, SCORE_EVENTS, seed, order = NULL, reference_group = NULL, population = NULL)
  sc[sc$check == "pair_emd", ] |>
    summarise(value = max(value), .by = marker)
}

one_run <- function(seed, cd3_sd) {
  set.seed(seed * 100 + round(cd3_sd * 100))
  batches <- paste0("b", seq_len(N_BATCH))
  gain <- matrix(exp(stats::rnorm(N_BATCH * length(CHANNELS), 0, GAIN_SD)), N_BATCH)
  offset <- matrix(stats::runif(N_BATCH * length(CHANNELS), 0, OFFSET_MAX), N_BATCH)
  dimnames(gain) <- dimnames(offset) <- list(batches, CHANNELS)

  ids <- paste0("d", seq_len(N_BATCH * DONORS_PER_BATCH))
  batch <- stats::setNames(rep(batches, each = DONORS_PER_BATCH), ids)
  truth <- stats::setNames(lapply(ids, \(i) true_sample()), ids)
  observed <- lapply(ids, \(i) instrument(truth[[i]], offset[batch[[i]], ], gain[batch[[i]], ]))
  names(observed) <- ids

  # One CD3 control per donor and channel, through the donor's instrument.
  cd3_factor <- stats::setNames(exp(stats::rnorm(length(ids), 0, cd3_sd)), ids)
  control_map <- expand.grid(unit = ids, channel = CHANNELS, stringsAsFactors = FALSE)
  control_map$control <- paste(control_map$unit, control_map$channel, sep = "_")
  controls <- lapply(seq_len(nrow(control_map)), \(i) {
    u <- control_map$unit[[i]]
    raw <- true_control(control_map$channel[[i]], cd3_factor[[u]])
    instrument(raw, offset[batch[[u]], ], gain[batch[[u]], ])
  })
  names(controls) <- control_map$control

  positive <- sw_relativise_split_kmeans(
    controls,
    control_map,
    \(v) asinh(v / COFACTOR),
    seed = seed,
    nstart = KM_NSTART
  )
  model <- sw_relativise_fit(controls, control_map, positive, MIN_EVENTS)
  relativised <- sw_relativise_apply(
    observed,
    model,
    unit = stats::setNames(ids, ids),
    channels = CHANNELS,
    low = 0,
    high = CD3_LEVEL,
    input_scale = "linear"
  )

  # E1: one control per batch and channel from one fixed reference material (CD3 factor 1 in
  # every batch), acquired on that batch's instrument; unit = batch.
  batch_info <- expand.grid(batch = batches, channel = CHANNELS, stringsAsFactors = FALSE)
  batch_info$control <- paste(batch_info$batch, batch_info$channel, sep = "_")
  batch_info$reference_marker <- "CD3"
  batch_controls <- lapply(seq_len(nrow(batch_info)), \(i) {
    b <- batch_info$batch[[i]]
    instrument(true_control(batch_info$channel[[i]], 1), offset[b, ], gain[b, ])
  })
  names(batch_controls) <- batch_info$control
  batch_map <- sw_relativise_batch_map(batch_info)
  batch_positive <- sw_relativise_split_kmeans(
    batch_controls,
    batch_map,
    \(v) asinh(v / COFACTOR),
    seed = seed,
    nstart = KM_NSTART
  )
  batch_model <- sw_relativise_fit(batch_controls, batch_map, batch_positive, MIN_EVENTS)
  per_batch <- sw_relativise_apply(
    observed,
    batch_model,
    unit = batch,
    channels = CHANNELS,
    low = 0,
    high = CD3_LEVEL,
    input_scale = "linear"
  )

  # Batch EMD per channel, and the paper's OTD per file, for every version.
  versions <- list(
    "truth (no instrument)" = truth,
    "observed" = observed,
    "relativised per sample" = relativised,
    "relativised per batch (E1)" = per_batch
  )
  emd <- bind_rows(lapply(names(versions), \(v) {
    max_pair_emd(arcsinh(versions[[v]]), batch, seed) |> mutate(version = v)
  })) |>
    mutate(seed = seed, cd3_donor_sd = cd3_sd)
  otd <- bind_rows(lapply(names(versions), \(v) {
    sw_otd_scores(
      arcsinh(versions[[v]]),
      CHANNELS,
      OTD_EVENTS,
      OTD_REFERENCE,
      OTD_BLUR,
      seed,
      OTD_MAX_ITER,
      OTD_TOL
    ) |>
      mutate(version = v, batch = batch[file])
  })) |>
    mutate(seed = seed, cd3_donor_sd = cd3_sd)

  # Does the diagnose step recover the true batch gain?
  dx <- sw_relativise_diagnose(model, batch = batch, order = NULL)
  ratio <- dx[dx$check == "gain_batch_ratio", ]
  true_rel <- vapply(seq_len(nrow(ratio)), \(i) {
    g <- gain[, ratio$channel[[i]]]
    g[[ratio$unit[[i]]]] / stats::median(g[batch])
  }, numeric(1))
  recovery <- tibble(
    seed = seed,
    cd3_donor_sd = cd3_sd,
    cor_log_gain = stats::cor(log(ratio$value), log(true_rel)),
    max_abs_gain_z = max(abs(dx$value[dx$check == "gain_z"]), na.rm = TRUE)
  )
  list(emd = emd, otd = otd, recovery = recovery)
}

runs <- lapply(SEEDS, \(s) lapply(CD3_DONOR_SD, \(sd) one_run(s, sd))) |>
  unlist(recursive = FALSE)
long <- bind_rows(lapply(runs, `[[`, "emd")) |>
  rename(channel = marker, max_pair_emd = value)
recovery <- bind_rows(lapply(runs, `[[`, "recovery"))
otd <- bind_rows(lapply(runs, `[[`, "otd"))
write.csv(long, file.path(out_dir, "relativisation_long.csv"), row.names = FALSE)
write.csv(otd, file.path(out_dir, "otd_long.csv"), row.names = FALSE)
write.csv(recovery, file.path(out_dir, "gain_recovery.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("wrote", nrow(long), "rows to", out_dir, "\n")
