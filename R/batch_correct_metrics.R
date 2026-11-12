# R/batch_correct_metrics.R
# Metric records for batch-correction candidates, in the long shape that
# sw_batch_candidates_from_metrics() reads (ROADMAP A2c). The three decision
# metrics reuse the helpers behind sw_corrector_compare(); LISI and kBET come
# from sw_batch_mixing() and are report only.

#' Metric records for batch-correction candidates
#'
#' @description
#' Scores each corrected version of the data against the uncorrected data and
#' returns one long table. Pass it to [sw_batch_candidates_from_metrics()] and
#' then [sw_batch_decide()]. The uncorrected data enter as the identity
#' candidate, named `identity_id`.
#'
#' Metrics, per candidate:
#'
#' * `emd_within_population`: for every pair of batches, population and
#'   marker, the mean absolute difference between the two batches' quantile
#'   functions at `n_quantiles` probability midpoints (a 1-Wasserstein
#'   distance on pooled cells). Cells of all samples in a batch are pooled.
#' * `mad_compression_fraction`: per marker, the median over samples of
#'   `max(0, 1 - MAD after / MAD before)`.
#' * `joint_structure_distance`: the mean over samples of the mean absolute
#'   difference between the marker correlation matrices before and after.
#' * `lisi_batch` and `kbet_rejection`: batch mixing from [sw_batch_mixing()]
#'   (Korsunsky 2019; Buttner 2019), overall (`comparison_id` =
#'   `"all_batches"`) and per batch (`"batch:<label>"`). These are report
#'   only. [sw_batch_candidates_from_metrics()] does not read them, and no
#'   limit on them exists unless you pass one to the decision layer.
#'
#' @details
#' **No donor or reference structure.** The within-population EMD here compares
#' pooled cells of two batches, so it contains biological differences between
#' the samples in each batch. For an EMD on repeated reference samples only,
#' use [sw_corrector_compare()] with a design audit.
#'
#' **Units.** Samples are the replicates for the MAD and structure metrics.
#' Cells are not replicates. A sample with fewer than `min_cells` cells in a
#' population is left out of that population's MAD and structure records.
#'
#' **Not assessed.** An EMD, MAD or structure record that lacks cells has a
#' status other than `"ok"` and `NA` value, as in [sw_corrector_compare()].
#' Mixing units under `min_events` or with fewer than two batches get status
#' `"not_assessed"`. Samples with fewer than `n_events` cells are left out of
#' the mixing draw with a warning that names them.
#'
#' **Weights.** `weights = "samples"` sets each record's weight to the number
#' of samples behind it (for EMD, the samples with at least `min_cells` cells
#' in the pair and population). `weights = "equal"` sets every weight to 1.
#'
#' **Seed.** The mixing draw uses `seed` and the caller's RNG state is
#' restored. Every candidate gets the same draw.
#'
#' No pass or fail, and no roll-up over markers.
#'
#' @param candidates A named list of corrected data, each with the same rows
#'   and the `markers` columns of `uncorrected`. Names must be unique and
#'   differ from `identity_id`.
#' @param uncorrected A numeric matrix or data frame, one row per cell, with
#'   columns named by `markers`. Already transformed by the caller.
#' @param sample_id One sample label per row of `uncorrected`.
#' @param batch A named vector of batch labels, one per sample, named by
#'   sample label.
#' @param population `NULL` (one unit, reported with an `NA` population), or
#'   one label per row of `uncorrected`. No `NA`.
#' @param markers Character vector of columns to score.
#' @param identity_id String, the candidate id given to `uncorrected`.
#' @param n_quantiles Whole number of quantile probabilities for the EMD
#'   (uncalibrated).
#' @param min_cells Whole number of at least 2, the fewest cells a sample or
#'   batch needs in a population (uncalibrated).
#' @param weights `"samples"` or `"equal"`, how records are weighted
#'   (uncalibrated).
#' @param k Whole number, neighbours per event for LISI and kBET
#'   (uncalibrated).
#' @param perplexity Number above 1 and below `k`, the LISI kernel
#'   perplexity (uncalibrated).
#' @param alpha Number in (0, 1), the kBET test level (uncalibrated).
#' @param n_events Whole number, events drawn per sample for mixing
#'   (uncalibrated).
#' @param seed Whole number for the mixing draw (uncalibrated).
#' @param min_events Whole number, fewest drawn events a mixing unit needs
#'   (uncalibrated).
#'
#' @return A tibble with columns `candidate_id`, `metric`, `population`,
#'   `marker`, `comparison_id`, `value`, `status`, `weight` and
#'   `calibration` (always `"uncalibrated"`).
#'
#' @seealso [sw_batch_candidates_from_metrics()], [sw_batch_decide()],
#'   [sw_corrector_compare()], [sw_batch_mixing()]
#'
#' @examples
#' set.seed(1)
#' ids <- rep(paste0("s", 1:4), each = 250)
#' batch <- c(s1 = "b1", s2 = "b1", s3 = "b2", s4 = "b2")
#' x <- cbind(CD3 = rnorm(1000), CD4 = rnorm(1000))
#' x[ids %in% c("s3", "s4"), "CD3"] <- x[ids %in% c("s3", "s4"), "CD3"] + 1
#' fixed <- x
#' fixed[ids %in% c("s3", "s4"), "CD3"] <- fixed[ids %in% c("s3", "s4"), "CD3"] - 1
#' records <- sw_correct_metrics(
#'   list(shift = fixed),
#'   uncorrected = x,
#'   sample_id = ids,
#'   batch = batch,
#'   population = NULL,
#'   markers = c("CD3", "CD4"),
#'   identity_id = "uncorrected",
#'   n_quantiles = 50, # illustrative
#'   min_cells = 50, # illustrative
#'   weights = "samples",
#'   k = 30, # illustrative
#'   perplexity = 10, # illustrative
#'   alpha = 0.05, # illustrative
#'   n_events = 200, # illustrative
#'   seed = 1,
#'   min_events = 100 # illustrative
#' )
#' sw_batch_candidates_from_metrics(records, identity_id = "uncorrected")
#'
#' @export
sw_correct_metrics <- function(
  candidates,
  uncorrected,
  sample_id,
  batch,
  population,
  markers,
  identity_id,
  n_quantiles,
  min_cells,
  weights,
  k,
  perplexity,
  alpha,
  n_events,
  seed,
  min_events
) {
  .sw_stop_missing(c(
    candidates = missing(candidates),
    uncorrected = missing(uncorrected),
    sample_id = missing(sample_id),
    batch = missing(batch),
    population = missing(population),
    markers = missing(markers),
    identity_id = missing(identity_id),
    n_quantiles = missing(n_quantiles),
    min_cells = missing(min_cells),
    weights = missing(weights),
    k = missing(k),
    perplexity = missing(perplexity),
    alpha = missing(alpha),
    n_events = missing(n_events),
    seed = missing(seed),
    min_events = missing(min_events)
  ))
  chk <- .sw_check_metrics_args(
    candidates,
    uncorrected,
    sample_id,
    batch,
    population,
    markers,
    identity_id,
    n_quantiles,
    min_cells,
    weights
  )
  mats <- chk$mats
  cells <- chk$cells
  pops <- sort(unique(cells$population))
  samples <- sort(unique(cells$sample_id))
  batches <- sort(unique(cells$batch))
  if (length(batches) < 2) {
    stop("`batch` needs at least two batches.", call. = FALSE)
  }
  probs <- (seq_len(n_quantiles) - 0.5) / n_quantiles
  pairs <- utils::combn(batches, 2, simplify = FALSE)

  # EMD, MAD compression and joint structure per pair and population, with
  # the helpers behind sw_corrector_compare().
  records <- list()
  for (pair in pairs) {
    cid <- paste(pair, collapse = "__vs__")
    in_pair <- cells$batch %in% pair
    for (pop in pops) {
      in_pop <- in_pair & cells$population == pop
      s_pop <- sort(unique(cells$sample_id[in_pop]))
      sample_rows <- lapply(s_pop, \(s) which(in_pop & cells$sample_id == s))
      sample_rows <- sample_rows[lengths(sample_rows) >= min_cells]
      pop_label <- if (chk$no_population) NA_character_ else pop
      for (cand in names(mats)) {
        emd <- .sw_compare_emd(
          mats[[cand]],
          cells,
          in_pop,
          pair,
          "all",
          markers,
          min_cells,
          probs,
          cand,
          pop_label,
          cid
        )
        emd$weight <- length(sample_rows)
        records[[length(records) + 1]] <- emd
        records[[length(records) + 1]] <- .sw_compare_mad(
          mats[[cand]],
          mats[[identity_id]],
          sample_rows,
          markers,
          cand,
          pop_label,
          cid
        )
        records[[length(records) + 1]] <- .sw_compare_structure(
          mats[[cand]],
          mats[[identity_id]],
          sample_rows,
          markers,
          cand,
          pop_label,
          cid
        )
      }
    }
  }

  # LISI and kBET: one seeded draw per sample, the same rows for every candidate.
  per_sample <- table(cells$sample_id)
  thin <- names(per_sample)[per_sample < n_events]
  if (length(thin) > 0) {
    warning(
      "Left out of the mixing draw (fewer than `n_events` cells): ",
      paste(thin, collapse = ", "),
      call. = FALSE
    )
  }
  mix_samples <- setdiff(samples, thin)
  for (cand in names(mats)) {
    records[[length(records) + 1]] <- .sw_metrics_mixing(
      mats[[cand]],
      cells,
      mix_samples,
      markers,
      cand,
      chk$no_population,
      k,
      perplexity,
      alpha,
      n_events,
      seed,
      min_events
    )
  }

  out <- dplyr::bind_rows(records)
  out <- out[
    order(
      match(out$candidate_id, names(mats)),
      out$metric,
      out$population,
      out$comparison_id,
      out$marker
    ),
  ]
  if (weights == "equal") {
    out$weight <- rep(1, nrow(out))
  }
  out$value <- unname(out$value)
  out$calibration <- rep("uncalibrated", nrow(out))
  tibble::as_tibble(out)
}

# Validate the arguments and build the cell table and the marker matrices.
.sw_check_metrics_args <- function(
  candidates,
  uncorrected,
  sample_id,
  batch,
  population,
  markers,
  identity_id,
  n_quantiles,
  min_cells,
  weights
) {
  if (
    !is.character(identity_id) ||
      length(identity_id) != 1 ||
      is.na(identity_id) ||
      !nzchar(identity_id)
  ) {
    stop("`identity_id` must be one non-empty string.", call. = FALSE)
  }
  if (
    !is.list(candidates) ||
      is.data.frame(candidates) ||
      is.null(names(candidates)) ||
      any(!nzchar(names(candidates))) ||
      anyDuplicated(names(candidates)) > 0
  ) {
    stop(
      "`candidates` must be a list with unique, non-empty names.",
      call. = FALSE
    )
  }
  if (identity_id %in% names(candidates)) {
    stop(
      "The name '",
      identity_id,
      "' is reserved for the uncorrected data (`identity_id`).",
      call. = FALSE
    )
  }
  if (!is.character(markers) || length(markers) == 0 || anyNA(markers)) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  .sw_check_count(n_quantiles, "n_quantiles")
  if (!is.numeric(min_cells) || length(min_cells) != 1 || min_cells < 2) {
    stop("`min_cells` must be one number of at least 2.", call. = FALSE)
  }
  if (!identical(weights, "samples") && !identical(weights, "equal")) {
    stop("`weights` must be \"samples\" or \"equal\".", call. = FALSE)
  }
  as_marker_matrix <- function(x, what) {
    absent <- setdiff(markers, colnames(x))
    if (length(absent) > 0) {
      stop(
        what,
        " lacks marker column(s): ",
        paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
    m <- as.matrix(x[, markers, drop = FALSE])
    if (!is.numeric(m) || any(!is.finite(m))) {
      stop(
        what,
        " must be numeric with no NA or infinite values.",
        call. = FALSE
      )
    }
    m
  }
  base <- as_marker_matrix(uncorrected, "`uncorrected`")
  n <- nrow(base)
  mats <- c(
    stats::setNames(list(base), identity_id),
    lapply(
      stats::setNames(names(candidates), names(candidates)),
      \(nm) {
        cand <- candidates[[nm]]
        if (nrow(cand) != n) {
          stop(
            "Candidate '",
            nm,
            "' must have the same rows as `uncorrected`.",
            call. = FALSE
          )
        }
        as_marker_matrix(cand, paste0("Candidate '", nm, "'"))
      }
    )
  )
  if (length(sample_id) != n || anyNA(sample_id)) {
    stop(
      "`sample_id` needs one non-missing label per row of `uncorrected`.",
      call. = FALSE
    )
  }
  sample_id <- as.character(sample_id)
  if (
    is.null(names(batch)) ||
      !all(unique(sample_id) %in% names(batch)) ||
      anyNA(batch[unique(sample_id)])
  ) {
    stop(
      "`batch` must be a named vector with a label for every sample.",
      call. = FALSE
    )
  }
  no_population <- is.null(population)
  if (no_population) {
    population <- rep("all", n)
  } else if (length(population) != n || anyNA(population)) {
    stop(
      "`population` must be NULL or have one non-missing label per row.",
      call. = FALSE
    )
  }
  cells <- data.frame(
    sample_id = sample_id,
    batch = as.character(batch[sample_id]),
    population = as.character(population),
    is_reference = TRUE,
    reference_id = "all",
    stringsAsFactors = FALSE
  )
  list(mats = mats, cells = cells, no_population = no_population)
}

# LISI and kBET records for one candidate, in the long metric shape.
.sw_metrics_mixing <- function(
  xm,
  cells,
  mix_samples,
  markers,
  cand,
  no_population,
  k,
  perplexity,
  alpha,
  n_events,
  seed,
  min_events
) {
  empty <- tibble::tibble(
    candidate_id = character(),
    metric = character(),
    population = character(),
    marker = character(),
    comparison_id = character(),
    value = double(),
    status = character(),
    weight = double()
  )
  keep <- cells$sample_id %in% mix_samples
  batch_by_sample <- tapply(cells$batch[keep], cells$sample_id[keep], `[[`, 1)
  if (length(unique(batch_by_sample)) < 2) {
    return(empty)
  }
  data <- lapply(
    stats::setNames(mix_samples, mix_samples),
    \(s) xm[cells$sample_id == s, , drop = FALSE]
  )
  pop <- if (no_population) {
    NULL
  } else {
    lapply(
      stats::setNames(mix_samples, mix_samples),
      \(s) cells$population[cells$sample_id == s]
    )
  }
  mix <- sw_batch_mixing(
    data,
    batch = batch_by_sample,
    markers = markers,
    k = k,
    perplexity = perplexity,
    alpha = alpha,
    n_events = n_events,
    seed = seed,
    min_events = min_events,
    population = pop
  )
  tibble::tibble(
    candidate_id = cand,
    metric = mix$check,
    population = mix$population,
    marker = NA_character_,
    comparison_id = .sw_mixing_comparison(mix$batch),
    value = mix$value,
    status = ifelse(mix$status == "assessed", "ok", "not_assessed"),
    weight = as.double(mix$n_events)
  )
}

# "all_batches" for the overall row, "batch:<label>" for a per-batch row.
.sw_mixing_comparison <- function(batch) {
  ifelse(is.na(batch), "all_batches", paste0("batch:", batch))
}
