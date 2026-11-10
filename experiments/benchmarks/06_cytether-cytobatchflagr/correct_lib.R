# Purpose: one correction call for T6, shared by 02_correct.R (EMD folds) and 02b_correct_op.R
#          (OpenProblems folds). Sourced; nothing runs here.
# Inputs:  x (cells x markers, arcsinh), meta (sample_id, batch, donor_id, is_reference,
#          is_holdout_reference, arrival), population (integer per cell), method, markers
# Outputs: list(x, supported, abstained, reason, fits_on_heldout, effects)

source("experiments/benchmarks/05_cytether-synthetic/adapters.R")
source("experiments/benchmarks/_common/adapter_cytofin.R")

# Parameters (approved, plan T6) -------------------------------------------------
ANCHOR <- "batch01" # Q8
MIN_CELLS <- 10 # G1
VARIANCE_FLOOR <- 1e-8 # heb numerical floor, as in the synthetic benchmark
HEB_POOL <- "population" # N4: lead variant
CYTONORM_SEED <- 1 # N6 for T6; overrides the synthetic adapter's 431
CYTOFIN_MODE <- "meanshift" # N7
MIN_BATCHES <- 2
MIN_BATCHES_ORIGIN <- "two batches are the least that define a shift"
T6_TRANSFORM <- list(name = "asinh", cofactor = 6000) # Q4

correct_t6 <- function(x, meta, population, method, markers) {
  if (method == "none") {
    return(list(x = x, supported = TRUE, abstained = 0L, reason = "raw", fits_on_heldout = FALSE))
  }
  if (method == "cycombine") return(run_cycombine(x, meta, markers, ANCHOR))
  if (method == "cytonorm") return(run_cytonorm(x, meta, markers))
  if (method == "cytofin") return(run_cytofin(x, meta, markers, ANCHOR, mode = CYTOFIN_MODE))
  # cyTether engines: fit on reference cells, integrate every cell with abstention.
  samples <- meta |>
    dplyr::distinct(sample_id, batch, donor_id, is_reference) |>
    dplyr::mutate(
      reference_role = dplyr::if_else(is_reference, "technical_replicate", "none"),
      reference_id = dplyr::if_else(is_reference, donor_id, NA_character_),
      aliquot = dplyr::if_else(is_reference, sample_id, NA_character_)
    )
  manifest <- sw_sample_manifest(
    samples, sample_id = "sample_id", batch = "batch", donor = "donor_id", aliquot = "aliquot",
    reference_role = "reference_role", reference_id = "reference_id"
  )
  audit <- sw_design_audit(
    samples, minimum_batches = MIN_BATCHES, minimum_batches_origin = MIN_BATCHES_ORIGIN,
    require_shared_reference = TRUE, anchor = ANCHOR
  )
  # N4: flat hierarchy, one random effect per batch.
  hierarchy <- if (method == "heb") dplyr::mutate(dplyr::distinct(meta, batch), study = batch) else NULL
  xt <- tibble::as_tibble(x) |> dplyr::mutate(sample_id = meta$sample_id, .before = 1)
  fit_rows <- meta$is_reference
  model <- sw_corrector_fit(
    xt[fit_rows, ], manifest, audit, population = population[fit_rows], markers = markers,
    transform = T6_TRANSFORM, anchor = ANCHOR, engine = method,
    min_cells = MIN_CELLS, seed = 1, hierarchy = hierarchy, variance_floor = VARIANCE_FLOOR,
    pool = if (method == "heb") HEB_POOL else NULL,
    anchor_effect = if (method == "heb") "random" else NULL
  )
  out <- sw_corrector_integrate(model, xt, manifest, population = population,
                                transform = T6_TRANSFORM, unsupported = "abstain")
  list(x = as.matrix(out$x[, markers]), supported = TRUE,
       abstained = sum(out$actions$action == "abstain"), reason = paste("engine", method),
       fits_on_heldout = FALSE, effects = out$model$effects)
}
