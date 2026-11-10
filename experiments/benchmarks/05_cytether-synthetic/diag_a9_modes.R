#!/usr/bin/env Rscript
# Purpose: record the A9 anchor mode heb assigns to each population x marker, by refitting heb_pool_pop_kv_bap_anc
#          on one stage of each seed with a recorder on .sw_heb_design. Never reads the truth file.
# Run:     pixi run Rscript experiments/benchmarks/05_cytether-synthetic/diag_a9_modes.R <sim_dir> <stage> <out_csv> <seed> [<seed> ...]
# Inputs:  <sim_dir>/seed_<n>/<stage>/observed_X.csv.gz and cell_metadata.csv.gz (from export_observed.py)
# Outputs: <out_csv>   seed, anchor, population, marker, mode (low, high or NA when no anchor row joins),
#                      n_anchor_rows (anchor rows in that design)
#          <out_csv without .csv>_theta.csv   fitted variance components per seed and pooled set (marker)
# SW_DIAG_ANCHOR_DISCREPANCY = "none" (default) or "batch" sets plan A11 for the refit.
# Settings copy run_method.R for the heb_pool_pop_kv_bap_anc arm (plan A0-A9, F6, F8); keep them in step.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters (inherited from run_method.R, confirmed C5 2026-11-10) --------------
MIN_CELLS <- 10L # G1
SEED <- 1L
TRANSFORM <- list(name = "asinh", cofactor = 5)
ANCHOR_POSITIVE <- "CD3" # A1
ANCHOR_MIN_SEPARATION <- 3 # A2
ANCHOR_SEED <- 1 # A7
VARIANCE_FLOOR <- 1e-8
MIN_BATCHES <- 2
MIN_BATCHES_ORIGIN <- "two batches are the least that define a shift (synthetic harness)"

args <- commandArgs(trailingOnly = TRUE)
sim_dir <- args[[1]]
stage <- args[[2]]
out_csv <- args[[3]]
seeds <- args[-(1:3)]
stopifnot("output must not exist yet" = !file.exists(out_csv))

# Recorder on the design builder ----------------------------------------------------
# Each heb design is one population x marker; its anchor rows carry the A9 mode.
ns <- asNamespace("SpectraWeaveR")
design_orig <- get(".sw_heb_design", envir = ns)
log_env <- new.env()
log_env$rows <- list()
recorder <- function(
  d,
  hierarchy,
  anchor,
  anchor_effect,
  residual_model,
  anchor_rows = NULL,
  ...
) {
  log_env$rows[[length(log_env$rows) + 1]] <- tibble(
    population = as.character(d$population[[1]]),
    marker = as.character(d$marker[[1]]),
    mode = if (is.null(anchor_rows)) {
      NA_character_
    } else {
      paste(unique(anchor_rows$mode), collapse = "+")
    },
    n_anchor_rows = if (is.null(anchor_rows)) 0L else nrow(anchor_rows)
  )
  design_orig(
    d,
    hierarchy,
    anchor,
    anchor_effect,
    residual_model,
    anchor_rows = anchor_rows,
    ...
  )
}
unlockBinding(".sw_heb_design", ns)
assign(".sw_heb_design", recorder, envir = ns)

# Recorder on the variance components (A11 diagnostic) ---------------------------------
# One row per pooled set: study, run, residual, anchor-row extra and anchor discrepancy variances.
theta_orig <- get(".sw_heb_theta", envir = ns)
log_env$theta <- list()
theta_recorder <- function(designs, variance_floor, variance_prior) {
  res <- theta_orig(designs, variance_floor, variance_prior)
  th <- c(res$theta, rep(NA_real_, 5 - length(res$theta)))
  log_env$theta[[length(log_env$theta) + 1]] <- tibble(
    marker = as.character(designs[[1]]$d$marker[[1]]),
    n_designs = length(designs),
    n_designs_with_anchors = sum(vapply(
      designs,
      function(z) any(z$ra > 0),
      logical(1)
    )),
    study = th[[1]],
    run = th[[2]],
    resid_donor = th[[3]],
    resid_anchor_extra = th[[4]],
    tau_d = th[[5]]
  )
  res
}
unlockBinding(".sw_heb_theta", ns)
assign(".sw_heb_theta", theta_recorder, envir = ns)
# A11 setting for this diagnostic; A9 modes do not depend on it.
ANCHOR_DISCREPANCY <- Sys.getenv("SW_DIAG_ANCHOR_DISCREPANCY", "none")
# SW_DIAG_ANCHORS = "0" refits the lead arm (no anchors) to record its variance components (plan A13 step 1).
USE_ANCHORS <- Sys.getenv("SW_DIAG_ANCHORS", "1") == "1"
# SW_DIAG_ANCHOR_BATCHES = "all" (default) or "unreferenced" sets plan A13.
ANCHOR_BATCHES <- Sys.getenv("SW_DIAG_ANCHOR_BATCHES", "all")

# One seed: the run_method.R setup for the anchor arm, fit only ----------------------
modes_one_seed <- function(seed) {
  stage_dir <- file.path(sim_dir, paste0("seed_", seed), stage)
  x <- as.matrix(read_csv(
    file.path(stage_dir, "observed_X.csv.gz"),
    show_col_types = FALSE
  ))
  markers <- colnames(x)
  meta <- read_csv(
    file.path(stage_dir, "cell_metadata.csv.gz"),
    show_col_types = FALSE
  ) |>
    mutate(
      role = if_else(
        is_reference & !is_holdout_reference,
        "technical_replicate",
        "none"
      )
    )
  samples <- meta |>
    distinct(
      sample_id,
      batch,
      donor_id,
      role,
      is_reference,
      is_holdout_reference,
      arrival
    ) |>
    mutate(
      reference_role = role,
      reference_id = if_else(
        role == "technical_replicate",
        donor_id,
        NA_character_
      ),
      aliquot = if_else(role == "technical_replicate", sample_id, NA_character_)
    )
  fit_samples <- meta |>
    filter(arrival == "initial") |>
    distinct(sample_id)
  fit_rows <- meta$arrival == "initial" & !meta$is_holdout_reference

  # Anchor batch (Q8 rule): most reference donors, ties to the first sorted batch.
  anchor <- samples |>
    filter(role == "technical_replicate", arrival == "initial") |>
    summarise(n = n_distinct(donor_id), .by = batch) |>
    arrange(desc(n), batch) |>
    pull(batch) |>
    first()
  manifest <- sw_sample_manifest(
    samples,
    sample_id = "sample_id",
    batch = "batch",
    donor = "donor_id",
    aliquot = "aliquot",
    reference_role = "reference_role",
    reference_id = "reference_id"
  )
  audit <- sw_design_audit(
    samples |>
      filter(
        sample_id %in% fit_samples$sample_id,
        role == "technical_replicate" | !is_reference
      ),
    minimum_batches = MIN_BATCHES,
    minimum_batches_origin = MIN_BATCHES_ORIGIN,
    require_shared_reference = TRUE,
    anchor = anchor
  )
  anchors <- if (USE_ANCHORS) {
    sw_corrector_anchors(
      x[fit_rows, , drop = FALSE],
      meta$sample_id[fit_rows],
      meta$batch[fit_rows],
      markers,
      positive_markers = ANCHOR_POSITIVE,
      min_separation = ANCHOR_MIN_SEPARATION,
      seed = ANCHOR_SEED
    )
  } else {
    NULL
  }
  xt <- as_tibble(x) |>
    mutate(sample_id = meta$sample_id, .before = 1)

  log_env$rows <- list()
  log_env$theta <- list()
  sw_corrector_fit(
    xt[fit_rows, ],
    manifest,
    audit,
    population = meta$cell_type[fit_rows],
    markers = markers,
    transform = TRANSFORM,
    anchor = anchor,
    engine = "heb",
    min_cells = MIN_CELLS,
    seed = SEED,
    hierarchy = distinct(meta, batch, study),
    variance_floor = VARIANCE_FLOOR,
    pool = "population",
    anchor_effect = "random",
    variance_prior = "boundary_avoiding",
    residual_model = "known_plus_extra",
    anchors = anchors,
    anchor_discrepancy = ANCHOR_DISCREPANCY,
    anchor_batches = ANCHOR_BATCHES,
    scale_model = "none"
  )
  list(
    modes = list_rbind(log_env$rows) |>
      mutate(seed = as.integer(seed), anchor = anchor, .before = 1),
    theta = list_rbind(log_env$theta) |>
      mutate(
        seed = as.integer(seed),
        anchors = USE_ANCHORS,
        anchor_discrepancy = ANCHOR_DISCREPANCY,
        .before = 1
      )
  )
}

# All seeds --------------------------------------------------------------------------
fits <- map(seeds, modes_one_seed)
modes <- map(fits, "modes") |>
  list_rbind()
theta <- map(fits, "theta") |>
  list_rbind()
write_csv(modes, out_csv)
write_csv(theta, sub("\\.csv$", "_theta.csv", out_csv))
writeLines(
  capture.output(sessionInfo()),
  file.path(dirname(out_csv), "sessionInfo_modes.txt")
)
