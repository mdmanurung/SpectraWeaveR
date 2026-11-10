#!/usr/bin/env Rscript
# Purpose: correct the simulated stages S0 to S12 of one seed with one method and write the corrected matrices.
# Run:     pixi run Rscript experiments/benchmarks/05_cytether-synthetic/run_method.R <seed_dir> <out_dir> <method>
#          method is one of: none, naive_batch_center (unsafe stress control), paired_shift, heb (pool none), heb_flat (no study level), heb_pool_pop, heb_pool_all, cycombine, cytonorm (adapters.R), reversed_shift (sanity: paired_shift with the correction applied in the wrong direction)
# Inputs:  <seed_dir>/<stage>/observed_X.csv.gz and cell_metadata.csv.gz (from export_observed.py; never the truth file)
# Outputs (in <out_dir>/<method>/<stage>/):
#   corrected_X.csv.gz   corrected matrix, same shape and row order as observed_X.csv.gz
#   status.csv      one row: supported, abstained_pairs, reason

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
source("experiments/benchmarks/05_cytether-synthetic/adapters.R")

# Parameters (approved by the user unless marked) ---------------------------------
MIN_CELLS <- as.integer(Sys.getenv("SW_MIN_CELLS", "10")) # G1: approved value 10; the env override is for labelled diagnostics only
SEED <- 1L # paired_shift draws no random numbers; recorded in the model only
TRANSFORM <- list(name = "asinh", cofactor = 5) # nominal label: simulated data are already on one scale
ANCHOR_POSITIVE <- "CD3" # plan A1: CD3+ is the only positive-mode anchor
ANCHOR_MIN_SEPARATION <- 3 # plan A2: pooled SDs
ANCHOR_SEED <- 1 # plan A7: k-means seed
VARIANCE_FLOOR <- 1e-8 # numerical floor on variance components (heb only); sensitivity to be ablated
MIN_BATCHES <- 2 # logical minimum to estimate a shift; gates eligibility only
MIN_BATCHES_ORIGIN <- "two batches are the least that define a shift (synthetic harness)"

args <- commandArgs(trailingOnly = TRUE)
seed_dir <- args[[1]]
method <- args[[3]]
stopifnot(
  "unknown method" = method %in%
    c(
      "none",
      "naive_batch_center",
      "paired_shift",
      "heb",
      "heb_afix",
      "heb_shared",
      "heb_flat",
      "heb_pool_pop",
      "heb_pool_all",
      "heb_bap",
      "heb_pool_pop_bap",
      "heb_pool_pop_kv",
      "heb_pool_pop_kv_bap",
      "heb_pool_all_kv_bap",
      "heb_afix_pop_kv_bap",
      "heb_pool_pop_kv_bap_anc",
      "heb_pool_pop_kv_bap_anc_d",
      "heb_pool_pop_kv_bap_anc_u",
      "heb_pool_pop_kv_bap_sc",
      "heb_pool_pop_kv_bap_anc_u_sc",
      "reversed_shift",
      "cycombine",
      "cytonorm"
    )
)

# One stage ------------------------------------------------------------------------
correct_stage <- function(stage_dir) {
  x <- as.matrix(read_csv(
    file.path(stage_dir, "observed_X.csv.gz"),
    show_col_types = FALSE
  ))
  markers <- colnames(x)
  meta <- read_csv(
    file.path(stage_dir, "cell_metadata.csv.gz"),
    show_col_types = FALSE
  )
  if (method == "none") {
    return(list(x = x, supported = TRUE, abstained = 0L, reason = "raw"))
  }

  if (method == "naive_batch_center") {
    # Bundle baseline: per population, move each batch mean onto the first sorted batch; uses all cells.
    out <- x
    first <- sort(unique(meta$batch))[[1]]
    for (pop in unique(meta$cell_type)) {
      target <- colMeans(x[
        meta$batch == first & meta$cell_type == pop,
        ,
        drop = FALSE
      ])
      for (b in unique(meta$batch)) {
        rows <- meta$batch == b & meta$cell_type == pop
        if (
          sum(rows) >= 3 &&
            sum(meta$batch == first & meta$cell_type == pop) >= 3
        ) {
          out[rows, ] <- sweep(
            x[rows, , drop = FALSE],
            2,
            colMeans(x[rows, , drop = FALSE]) - target
          )
        }
      }
    }
    return(list(
      x = out,
      supported = TRUE,
      abstained = 0L,
      reason = "first sorted batch",
      fits_on_heldout = TRUE
    ))
  }

  # Which cells may calibrate: initial, non-holdout reference cells only.
  meta <- meta |>
    mutate(
      calibrates = is_reference & !is_holdout_reference & arrival == "initial",
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
      )
    )
  fit_samples <- meta |>
    filter(arrival == "initial") |>
    distinct(sample_id)
  fit_rows <- meta$arrival == "initial" & !meta$is_holdout_reference
  initial_refs <- samples |>
    filter(role == "technical_replicate", arrival == "initial")
  if (nrow(initial_refs) == 0) {
    return(list(
      x = NULL,
      supported = FALSE,
      abstained = NA_integer_,
      reason = "no reference samples"
    ))
  }

  # Anchor (Q8 rule): batch with the most reference donors, ties to the first sorted batch.
  anchor <- initial_refs |>
    summarise(n = n_distinct(donor_id), .by = batch) |>
    arrange(desc(n), batch) |>
    pull(batch) |>
    first()

  if (method == "cycombine") {
    return(run_cycombine(x, meta, markers, anchor))
  }
  if (method == "cytonorm") {
    return(run_cytonorm(x, meta, markers))
  }

  manifest_df <- samples |>
    mutate(
      aliquot = if_else(role == "technical_replicate", sample_id, NA_character_)
    )
  manifest <- sw_sample_manifest(
    manifest_df,
    sample_id = "sample_id",
    batch = "batch",
    donor = "donor_id",
    aliquot = "aliquot",
    reference_role = "reference_role",
    reference_id = "reference_id"
  )
  audit <- sw_design_audit(
    manifest_df |>
      filter(
        sample_id %in% fit_samples$sample_id,
        role == "technical_replicate" | !is_reference
      ),
    minimum_batches = MIN_BATCHES,
    minimum_batches_origin = MIN_BATCHES_ORIGIN,
    require_shared_reference = TRUE,
    anchor = anchor
  )

  xt <- as_tibble(x) |> mutate(sample_id = meta$sample_id, .before = 1)
  pop <- meta$cell_type
  engine <- if (method == "reversed_shift") {
    "paired_shift"
  } else if (method == "heb_shared") {
    "heb_shared"
  } else if (startsWith(method, "heb")) {
    "heb"
  } else {
    method
  }
  heb_pool <- c(
    heb = "none",
    heb_afix = "none",
    heb_flat = "none",
    heb_pool_pop = "population",
    heb_pool_all = "all",
    heb_shared = NA,
    heb_bap = "none",
    heb_pool_pop_bap = "population",
    heb_pool_pop_kv = "population",
    heb_pool_pop_kv_bap = "population",
    heb_pool_all_kv_bap = "all",
    heb_afix_pop_kv_bap = "population",
    heb_pool_pop_kv_bap_anc = "population",
    heb_pool_pop_kv_bap_anc_d = "population",
    heb_pool_pop_kv_bap_anc_u = "population",
    heb_pool_pop_kv_bap_sc = "population", # plan F12
    heb_pool_pop_kv_bap_anc_u_sc = "population" # F12 on S18, where the lead abstains
  )[method] # *_bap: F6 prior; *_kv: F8 known variance
  hierarchy <- if (engine %in% c("heb", "heb_shared")) {
    distinct(meta, batch, study)
  } else {
    NULL
  }
  # Plan A0-A9 (*_anc): internal anchors from every sample that may calibrate an anchor: initial and not a
  # held-out reference. Participant samples count, so batches without reference donors get anchor rows.
  anc_rows <- fit_rows
  heb_anchors <- if (grepl("_anc", method, fixed = TRUE)) {
    # *_anc, *_anc_d (plan A11), *_anc_u (plan A13)
    sw_corrector_anchors(
      x[anc_rows, , drop = FALSE],
      meta$sample_id[anc_rows],
      meta$batch[anc_rows],
      markers,
      positive_markers = ANCHOR_POSITIVE,
      min_separation = ANCHOR_MIN_SEPARATION,
      seed = ANCHOR_SEED
    )
  } else {
    NULL
  }
  if (method == "heb_flat") {
    hierarchy$study <- hierarchy$batch # flat arm: one random effect per batch, no study level
  }
  model <- sw_corrector_fit(
    xt[fit_rows, ],
    manifest,
    audit,
    population = pop[fit_rows],
    markers = markers,
    transform = TRANSFORM,
    anchor = anchor,
    engine = engine,
    min_cells = MIN_CELLS,
    seed = SEED,
    hierarchy = hierarchy,
    variance_floor = VARIANCE_FLOOR,
    pool = if (engine == "heb") unname(heb_pool) else NULL,
    anchor_effect = if (engine == "heb") {
      (if (startsWith(method, "heb_afix")) "fixed" else "random")
    } else {
      NULL
    },
    variance_prior = if (engine == "heb") {
      (if (grepl("_bap", method, fixed = TRUE)) "boundary_avoiding" else "none")
    } else {
      NULL
    },
    residual_model = if (engine == "heb") {
      (if (grepl("_kv", method, fixed = TRUE)) {
        "known_plus_extra"
      } else {
        "estimated"
      })
    } else {
      NULL
    },
    anchors = if (engine == "heb") heb_anchors else NULL,
    anchor_discrepancy = if (engine == "heb") {
      (if (endsWith(method, "_anc_d")) "batch" else "none")
    } else {
      NULL
    }, # plan A11
    anchor_batches = if (engine == "heb") {
      (if (grepl("_anc_u", method, fixed = TRUE)) "unreferenced" else "all")
    } else {
      NULL
    }, # plan A13
    scale_model = if (engine == "heb") {
      (if (endsWith(method, "_sc")) "log_sd" else "none")
    } else {
      NULL
    } # plan F12
  )

  # New batches first, without holdout cells, so holdout cells never calibrate.
  non_holdout <- !meta$is_holdout_reference
  res1 <- sw_corrector_integrate(
    model,
    xt[non_holdout, ],
    manifest,
    population = pop[non_holdout],
    transform = TRANSFORM,
    unsupported = "abstain"
  )
  # Historical invariance: integrating old rows alone must give the same values as integrating them with the new data.
  old <- non_holdout & meta$arrival == "initial"
  res_old <- sw_corrector_integrate(
    model,
    xt[old, ],
    manifest,
    population = pop[old],
    transform = TRANSFORM,
    unsupported = "abstain"
  )
  historical_max_delta <- max(abs(
    as.matrix(res_old$x[, markers]) -
      as.matrix(res1$x[old[non_holdout], markers])
  ))
  res2 <- sw_corrector_integrate(
    res1$model,
    xt[!non_holdout, ],
    manifest,
    population = pop[!non_holdout],
    transform = TRANSFORM,
    unsupported = "abstain"
  )
  out <- x
  out[non_holdout, ] <- as.matrix(res1$x[, markers])
  out[!non_holdout, ] <- as.matrix(res2$x[, markers])
  abstained <- sum(res1$actions$action == "abstain") +
    sum(res2$actions$action == "abstain")
  if (method == "reversed_shift") {
    out <- 2 * x - out # undo the correction, then apply it with the opposite sign
  }
  list(
    x = out,
    supported = TRUE,
    abstained = abstained,
    reason = paste0("anchor=", anchor),
    hist_delta = historical_max_delta,
    effects = res2$model$effects
  )
}

# All stages ----------------------------------------------------------------------
out_root <- file.path(args[[2]], method)
stages <- list.files(seed_dir, pattern = "^S[0-9]+[bw]?$") # S13b and S18w are stages too
# SW_STAGES="S18w" runs only the listed stages, so a missed stage can run beside finished ones.
if (nzchar(Sys.getenv("SW_STAGES"))) {
  stages <- intersect(stages, strsplit(Sys.getenv("SW_STAGES"), ",")[[1]])
}
for (stage in stages) {
  out_dir <- file.path(out_root, stage)
  stopifnot("output dir must not exist yet" = !dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)
  res <- tryCatch(
    correct_stage(file.path(seed_dir, stage)),
    error = function(e) {
      list(
        x = NULL,
        supported = FALSE,
        abstained = NA_integer_,
        reason = conditionMessage(e)
      )
    }
  )
  if (!is.null(res$x)) {
    write_csv(as_tibble(res$x), file.path(out_dir, "corrected_X.csv.gz"))
  }
  if (!is.null(res$effects)) {
    # Estimates and standard errors per (batch, population, marker); evaluate.py scores their coverage.
    write_csv(res$effects, file.path(out_dir, "effects.csv.gz"))
  }
  write_csv(
    tibble(
      stage = stage,
      method = method,
      supported = res$supported,
      abstained_pairs = res$abstained,
      reason = res$reason,
      historical_max_delta = res$hist_delta %||% NA_real_,
      min_cells = MIN_CELLS, # recorded so a diagnostic with another value is visible in every status file
      fits_on_heldout = isTRUE(res$fits_on_heldout), # TRUE: the method saw held-out reference or future cells while fitting
      ref_mode_ok = res$ref_mode_ok %||% NA
    ),
    file.path(out_dir, "status.csv")
  )
}
writeLines(
  capture.output(sessionInfo()),
  file.path(out_root, "sessionInfo.txt")
)
