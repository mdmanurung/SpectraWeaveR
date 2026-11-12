# R/ref_coverage.R
# Does the reference span the study range? CytoNorm asks that references cover
# the expression range of the study (Van Gassen 2020, idea only; GPL code not
# read). This function counts study events outside chosen reference quantiles.

#' Share of study events outside the reference range
#'
#' @description
#' For each batch, marker and population, take the reference quantiles at
#' `probs[1]` and `probs[2]` and report the share of study events below the
#' lower one and above the upper one. If the reference spans the study, the
#' shares sit near `probs[1]` and `1 - probs[2]`. Much larger shares mean the
#' study has events in a range where the reference has few or none, so a
#' correction fitted on the reference extrapolates there.
#'
#' @details
#' **Scale is yours.** `reference` and `study` must be on the same, already
#' transformed scale (DR37). Nothing is transformed here. `transform` and
#' `cofactor` are recorded only.
#'
#' **Quantiles** use `stats::quantile()` type 7 on all reference events of the
#' population. All events are used; there is no draw, so no seed.
#'
#' **Units.** Shares are over events and describe coverage. They are not
#' replicate-level evidence: donors and aliquots are the replicates.
#'
#' **Not assessed.** A row with fewer than `min_events` reference events or
#' fewer than `min_events` study events in the batch and population has `NA`
#' shares and `status = "not_assessed"`.
#'
#' **No roll-up.** Rows are per batch, marker and population (Q14).
#'
#' @param reference A numeric matrix or data frame of reference events (pool
#'   the reference files first), one column per marker.
#' @param study A named list of numeric matrices or data frames, one per study
#'   file, with the same marker columns.
#' @param batch A named vector of batch labels covering every file in `study`.
#' @param markers Character vector of columns to score.
#' @param probs Numeric vector `c(lower, upper)` with
#'   `0 < lower < upper < 1`, the reference quantile probabilities
#'   (uncalibrated: the plan proposes 0.01 and 0.99 from a weak source).
#' @param reference_labels,study_labels Both `NULL` (one population `"all"`),
#'   or population labels: a vector with one label per reference event, and a
#'   named list with one label vector per study file.
#' @param min_events Whole number, fewest events a population needs in the
#'   reference and in a study batch (uncalibrated: no sourced value).
#' @param transform,cofactor Text note and number for the scale. Recorded only.
#'
#' @return A tibble of class `sw_ref_coverage` with columns `batch`, `marker`,
#'   `population`, `share_below`, `share_above`, `expected_below`,
#'   `expected_above`, `n_reference_cells`, `n_study_cells`, `status` and
#'   `reason`.
#'
#' @seealso [sw_ref_vs_study()], [sw_ref_composition()]
#'
#' @examples
#' set.seed(1)
#' ref <- cbind(A = stats::rnorm(2000))
#' study <- list(s1 = cbind(A = stats::rnorm(500)), s2 = cbind(A = stats::rnorm(500, 1)))
#' sw_ref_coverage(
#'   ref, study,
#'   batch = c(s1 = "b1", s2 = "b2"),
#'   markers = "A",
#'   probs = c(0.01, 0.99), # illustrative
#'   reference_labels = NULL, study_labels = NULL,
#'   min_events = 100, # illustrative
#'   transform = "none", cofactor = NA_real_
#' )
#'
#' @export
sw_ref_coverage <- function(
  reference,
  study,
  batch,
  markers,
  probs,
  reference_labels,
  study_labels,
  min_events,
  transform,
  cofactor
) {
  .sw_stop_missing(c(
    reference = missing(reference),
    study = missing(study),
    batch = missing(batch),
    markers = missing(markers),
    probs = missing(probs),
    reference_labels = missing(reference_labels),
    study_labels = missing(study_labels),
    min_events = missing(min_events),
    transform = missing(transform),
    cofactor = missing(cofactor)
  ))
  .sw_check_count(min_events, "min_events")
  .sw_check_scale_note(transform, cofactor)
  if (
    !is.numeric(probs) ||
      length(probs) != 2 ||
      anyNA(probs) ||
      probs[[1]] <= 0 ||
      probs[[2]] >= 1 ||
      probs[[1]] >= probs[[2]]
  ) {
    stop(
      "`probs` must be c(lower, upper) with 0 < lower < upper < 1.",
      call. = FALSE
    )
  }
  if (!is.character(markers) || length(markers) == 0) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.list(study) || is.null(names(study)) || anyDuplicated(names(study))) {
    stop("`study` must be a named list of files.", call. = FALSE)
  }
  if (is.null(names(batch)) || !all(names(study) %in% names(batch))) {
    stop(
      "`batch` must be a named vector covering every file in `study`.",
      call. = FALSE
    )
  }
  if (is.null(reference_labels) != is.null(study_labels)) {
    stop(
      "Give both `reference_labels` and `study_labels`, or neither.",
      call. = FALSE
    )
  }
  ref <- as.matrix(reference)
  stud <- lapply(study, as.matrix)
  has_cols <- \(m) all(markers %in% colnames(m))
  if (!has_cols(ref) || !all(vapply(stud, has_cols, logical(1)))) {
    stop(
      "`reference` and every study file need all `markers` as columns.",
      call. = FALSE
    )
  }
  if (is.null(reference_labels)) {
    ref_pop <- rep("all", nrow(ref))
    stud_pop <- lapply(stud, \(m) rep("all", nrow(m)))
  } else {
    ref_pop <- as.character(reference_labels)
    stud_pop <- lapply(names(stud), \(f) as.character(study_labels[[f]])) |>
      stats::setNames(names(stud))
    ok <- length(ref_pop) == nrow(ref) &&
      all(vapply(
        names(stud),
        \(f) length(stud_pop[[f]]) == nrow(stud[[f]]),
        logical(1)
      ))
    if (!ok) {
      stop("Labels must have one entry per event.", call. = FALSE)
    }
    if (anyNA(ref_pop) || anyNA(unlist(stud_pop))) {
      stop("Labels must not contain NA.", call. = FALSE)
    }
  }
  bat <- as.character(batch[names(stud)])
  if (anyNA(bat)) {
    stop("`batch` has missing labels.", call. = FALSE)
  }

  # Study events of a batch are pooled over its files.
  rows <- list()
  for (b in sort(unique(bat))) {
    fs <- names(stud)[bat == b]
    sm <- do.call(rbind, lapply(stud[fs], \(m) m[, markers, drop = FALSE]))
    sp <- unlist(stud_pop[fs], use.names = FALSE)
    for (pop in sort(unique(c(ref_pop, sp)))) {
      rsel <- ref_pop == pop
      ssel <- sp == pop
      n_ref <- sum(rsel)
      n_stu <- sum(ssel)
      enough <- n_ref >= min_events && n_stu >= min_events
      for (mk in markers) {
        rv <- ref[rsel, mk]
        sv <- sm[ssel, mk]
        if (enough && anyNA(c(rv, sv))) {
          stop(
            "Events must not contain NA in the chosen markers.",
            call. = FALSE
          )
        }
        q <- if (enough) {
          stats::quantile(rv, probs, names = FALSE)
        } else {
          c(NA_real_, NA_real_)
        }
        rows[[length(rows) + 1]] <- tibble::tibble(
          batch = b,
          marker = mk,
          population = pop,
          share_below = if (enough) mean(sv < q[[1]]) else NA_real_,
          share_above = if (enough) mean(sv > q[[2]]) else NA_real_,
          expected_below = probs[[1]],
          expected_above = 1 - probs[[2]],
          n_reference_cells = n_ref,
          n_study_cells = n_stu,
          status = if (enough) "ok" else "not_assessed",
          reason = if (enough) {
            NA_character_
          } else {
            "fewer than min_events events in the reference or the study batch"
          }
        )
      }
    }
  }
  out <- dplyr::bind_rows(rows)
  structure(
    out,
    class = c("sw_ref_coverage", class(out)),
    probs = probs,
    transform = transform,
    cofactor = cofactor
  )
}
