#' Declare a batch-correction decision policy
#'
#' A policy is the set of gates that \code{\link{sw_batch_decide}} applies to
#' each correction candidate. Every threshold changes which candidate can be
#' recommended, so every argument is required and none has a default. There
#' are no preset policies.
#'
#' Ported from the decision policy of batchdx-cyto (MIT). That project ships
#' the values \code{min_coverage = 0.6}, \code{min_emd_reduction = 0.2},
#' \code{max_mad_compression_fraction = 0.25},
#' \code{max_joint_structure_change = 0.1} and \code{no_effect_emd = 0.05}.
#' Those values are **uncalibrated engineering defaults** from batchdx-cyto.
#' They are not validated cutoffs and must not be cited as such.
#'
#' A threshold counts as calibrated only when its origin does not contain the
#' word "uncalibrated" and \code{calibration_artifacts} names an artifact for
#' it. A policy with any uncalibrated threshold can yield
#' \code{eligible_for_review}, never \code{recommended_under_policy}.
#'
#' @param min_coverage Numeric in \[0, 1\]. Minimum fraction of the identity
#'   candidate's usable comparisons that a candidate must also cover.
#' @param min_emd_reduction Numeric. Minimum relative reduction in
#'   within-population EMD against the identity candidate.
#' @param max_mad_compression_fraction Numeric. Largest acceptable fraction of
#'   MAD compression; above it the candidate is rejected for distortion.
#' @param max_joint_structure_change Numeric. Largest acceptable median joint
#'   structure distance; above it the candidate is rejected for distortion.
#' @param no_effect_emd Numeric. If the identity candidate's worst-marker
#'   median EMD is below this value, correction is \code{not_needed}.
#' @param require_reference_supported Logical. If \code{TRUE}, a candidate is
#'   recommended only when every comparison its evidence rests on is
#'   \code{reference_supported} in the design audit.
#' @param threshold_origins Named character vector with exactly the five
#'   numeric threshold names. Each value states where the threshold came from,
#'   for example \code{"uncalibrated engineering default"}.
#' @param calibration_artifacts Named character vector mapping threshold names
#'   to a calibration artifact (path, URI or hash). Pass
#'   \code{character(0)} when no threshold is calibrated. The artifact is not
#'   re-verified here.
#'
#' @return An object of class \code{sw_batch_policy}: a list with the six
#'   settings, \code{threshold_origins}, \code{calibration_artifacts},
#'   \code{any_threshold_uncalibrated} and \code{settings} (a tibble with one
#'   row per threshold).
#'
#' @seealso \code{\link{sw_batch_decide}}
#'
#' @examples
#' origin <- "uncalibrated engineering default (batchdx-cyto)"
#' policy <- sw_batch_policy(
#'   min_coverage = 0.6,
#'   min_emd_reduction = 0.2,
#'   max_mad_compression_fraction = 0.25,
#'   max_joint_structure_change = 0.1,
#'   no_effect_emd = 0.05,
#'   require_reference_supported = TRUE,
#'   threshold_origins = c(
#'     min_coverage = origin, min_emd_reduction = origin,
#'     max_mad_compression_fraction = origin,
#'     max_joint_structure_change = origin, no_effect_emd = origin
#'   ),
#'   calibration_artifacts = character(0)
#' )
#' policy$any_threshold_uncalibrated
#'
#' @export
sw_batch_policy <- function(
  min_coverage,
  min_emd_reduction,
  max_mad_compression_fraction,
  max_joint_structure_change,
  no_effect_emd,
  require_reference_supported,
  threshold_origins,
  calibration_artifacts
) {
  args <- c(
    "min_coverage",
    "min_emd_reduction",
    "max_mad_compression_fraction",
    "max_joint_structure_change",
    "no_effect_emd",
    "require_reference_supported",
    "threshold_origins",
    "calibration_artifacts"
  )
  absent <- setdiff(args, names(match.call())[-1])
  if (length(absent) > 0) {
    stop(
      "Argument(s) ",
      paste0("'", absent, "'", collapse = ", "),
      " have no default; set each and record why.",
      call. = FALSE
    )
  }

  values <- list(
    min_coverage = min_coverage,
    min_emd_reduction = min_emd_reduction,
    max_mad_compression_fraction = max_mad_compression_fraction,
    max_joint_structure_change = max_joint_structure_change,
    no_effect_emd = no_effect_emd
  )
  thresholds <- names(values)
  for (nm in thresholds) {
    v <- values[[nm]]
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v)) {
      stop("'", nm, "' must be a single finite number.", call. = FALSE)
    }
  }
  if (min_coverage < 0 || min_coverage > 1) {
    stop("'min_coverage' must lie in [0, 1].", call. = FALSE)
  }
  if (
    !is.logical(require_reference_supported) ||
      length(require_reference_supported) != 1 ||
      is.na(require_reference_supported)
  ) {
    stop("'require_reference_supported' must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.character(threshold_origins) ||
      is.null(names(threshold_origins)) ||
      !setequal(names(threshold_origins), thresholds) ||
      anyDuplicated(names(threshold_origins)) > 0 ||
      anyNA(threshold_origins) ||
      any(threshold_origins == "")
  ) {
    stop(
      "'threshold_origins' must be a named character vector with one non-empty origin for each of: ",
      paste(thresholds, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (length(calibration_artifacts) == 0) {
    calibration_artifacts <- stats::setNames(character(0), character(0))
  }
  if (
    !is.character(calibration_artifacts) ||
      (length(calibration_artifacts) > 0 &&
        (is.null(names(calibration_artifacts)) ||
          !all(names(calibration_artifacts) %in% thresholds)))
  ) {
    stop(
      "'calibration_artifacts' must be a character vector named by threshold (",
      paste(thresholds, collapse = ", "),
      "), or character(0).",
      call. = FALSE
    )
  }

  # A threshold is calibrated only with a non-"uncalibrated" origin and a named artifact.
  origin <- unname(threshold_origins[thresholds])
  artifact <- unname(calibration_artifacts[thresholds])
  has_artifact <- !is.na(artifact) & artifact != ""
  calibrated <- !grepl("uncalibrated", origin, fixed = TRUE) & has_artifact

  structure(
    c(
      values,
      list(
        require_reference_supported = require_reference_supported,
        threshold_origins = threshold_origins[thresholds],
        calibration_artifacts = calibration_artifacts,
        any_threshold_uncalibrated = !all(calibrated),
        settings = tibble::tibble(
          threshold = thresholds,
          value = unlist(values[thresholds], use.names = FALSE),
          origin = origin,
          calibration_artifact = ifelse(has_artifact, artifact, NA_character_),
          calibrated = calibrated
        )
      )
    ),
    class = "sw_batch_policy"
  )
}

#' @export
print.sw_batch_policy <- function(x, ...) {
  cat(
    "<sw_batch_policy> require_reference_supported =",
    x$require_reference_supported,
    "\n"
  )
  print(x$settings)
  if (x$any_threshold_uncalibrated) {
    cat(
      "At least one threshold is uncalibrated; no candidate can be recommended.\n"
    )
  }
  invisible(x)
}

#' Decide between batch-correction candidates under a declared policy
#'
#' Assigns one state to each correction candidate and names a recommendation
#' only when exactly one candidate clears every gate. There is no weighted
#' score: a candidate must pass coverage, improvement and preservation gates
#' independently. Keeping the uncorrected data is always a valid outcome.
#'
#' Ported from the decision layer of batchdx-cyto (MIT). Gates run in this
#' order, and the first that applies sets the state:
#' \enumerate{
#'   \item status \code{"failed"}: \code{execution_failed}.
#'   \item status \code{"ineligible"} or \code{"not_available"}:
#'     \code{insufficient_evidence}.
#'   \item the design is below its floor: \code{below_design_floor}.
#'   \item every population-free batch contrast is non-identifiable:
#'     \code{non_identifiable}.
#'   \item identity candidate with worst-marker median EMD below
#'     \code{no_effect_emd}: \code{not_needed}; otherwise the identity
#'     candidate is \code{eligible_for_review}.
#'   \item a protected factor is confounded with batch: \code{non_identifiable}.
#'   \item MAD compression or joint structure change above its limit:
#'     \code{rejected_distortion}.
#'   \item no relative reduction, coverage missing or below
#'     \code{min_coverage}, or a guardrail metric missing:
#'     \code{insufficient_evidence}.
#'   \item reduction below \code{min_emd_reduction}, evidence not fully
#'     reference-supported (when required), labels from clustering, or any
#'     threshold uncalibrated: \code{eligible_for_review}.
#'   \item otherwise: \code{recommended_under_policy}.
#' }
#' Labels from clustering of the uncorrected baseline may themselves be
#' batch-driven, so \code{label_source = "clustering"} caps every candidate at
#' \code{eligible_for_review} under any policy.
#'
#' @param design An \code{sw_design_audit} object from
#'   \code{\link{sw_design_audit}}.
#' @param candidates A data frame with one row per candidate, including the
#'   identity candidate. Columns: \code{candidate_id}, \code{status} (one of
#'   \code{"ok"}, \code{"failed"}, \code{"ineligible"},
#'   \code{"not_available"}), \code{emd_median},
#'   \code{emd_worst_marker_median}, \code{emd_relative_reduction},
#'   \code{coverage}, \code{mad_compression_fraction},
#'   \code{joint_structure_distance_median}, and the list column
#'   \code{evidence_comparison_ids} (population-free comparison ids from
#'   \code{design$comparisons}). An optional \code{status_reason} column gives
#'   the reason for a non-\code{"ok"} status. Missing metrics are \code{NA}.
#'   \code{\link{sw_batch_candidates_from_metrics}} builds this table.
#' @param policy An \code{sw_batch_policy} object from
#'   \code{\link{sw_batch_policy}}.
#' @param identity_id The \code{candidate_id} of the uncorrected (identity)
#'   candidate. **No default.**
#' @param label_source Where the evaluation population labels came from: one
#'   of \code{"reference"}, \code{"manual_gating"} or \code{"clustering"}.
#'   **No default.**
#'
#' @return An object of class \code{sw_batch_decision}: a list with
#'   \code{assessments} (tibble: \code{candidate_id}, \code{state},
#'   \code{reason}, \code{guardrails_failed} list, the evidence columns and
#'   \code{reference_supported}), \code{recommended_candidate_id} (\code{NA}
#'   unless exactly one candidate is recommended),
#'   \code{recommendation_withheld_reason}, \code{pareto_front},
#'   \code{pareto_comparable}, \code{warnings}, \code{policy_settings}
#'   (thresholds with origins), \code{require_reference_supported},
#'   \code{identity_id} and \code{label_source}.
#'
#' @seealso \code{\link{sw_batch_policy}}, \code{\link{sw_design_audit}},
#'   \code{\link{sw_batch_candidates_from_metrics}}
#'
#' @examples
#' design <- sw_design_audit(
#'   data.frame(
#'     sample_id = paste0("s", 1:4),
#'     batch = rep(c("b1", "b2"), each = 2),
#'     condition = c("ctrl", "case", "ctrl", "case"),
#'     reference_role = c("anchor", "none", "anchor", "none"),
#'     reference_id = c("r1", NA, "r1", NA)
#'   ),
#'   minimum_batches = 2,
#'   minimum_batches_origin = "analyst choice for this example",
#'   require_shared_reference = TRUE,
#'   protected = "condition"
#' )
#' candidates <- tibble::tibble(
#'   candidate_id = c("uncorrected", "cycombine"),
#'   status = "ok",
#'   emd_median = c(0.8, 0.3),
#'   emd_worst_marker_median = c(1.2, 0.4),
#'   emd_relative_reduction = c(0, 0.6),
#'   coverage = c(1, 1),
#'   mad_compression_fraction = c(0, 0.05),
#'   joint_structure_distance_median = c(0, 0.02),
#'   evidence_comparison_ids = list("b1__vs__b2", "b1__vs__b2")
#' )
#' origin <- "uncalibrated engineering default (batchdx-cyto)"
#' policy <- sw_batch_policy(
#'   min_coverage = 0.6, min_emd_reduction = 0.2,
#'   max_mad_compression_fraction = 0.25, max_joint_structure_change = 0.1,
#'   no_effect_emd = 0.05, require_reference_supported = TRUE,
#'   threshold_origins = c(
#'     min_coverage = origin, min_emd_reduction = origin,
#'     max_mad_compression_fraction = origin,
#'     max_joint_structure_change = origin, no_effect_emd = origin
#'   ),
#'   calibration_artifacts = character(0)
#' )
#' decision <- sw_batch_decide(
#'   design, candidates, policy,
#'   identity_id = "uncorrected", label_source = "manual_gating"
#' )
#' decision$assessments[, c("candidate_id", "state")]
#'
#' @export
sw_batch_decide <- function(
  design,
  candidates,
  policy,
  identity_id,
  label_source
) {
  args <- c("design", "candidates", "policy", "identity_id", "label_source")
  absent <- setdiff(args, names(match.call())[-1])
  if (length(absent) > 0) {
    stop(
      "Argument(s) ",
      paste0("'", absent, "'", collapse = ", "),
      " are required and have no default.",
      call. = FALSE
    )
  }
  if (!inherits(design, "sw_design_audit")) {
    stop("'design' must be an sw_design_audit object.", call. = FALSE)
  }
  if (!inherits(policy, "sw_batch_policy")) {
    stop("'policy' must be an sw_batch_policy object.", call. = FALSE)
  }
  label_tiers <- c("reference", "manual_gating", "clustering")
  if (
    !is.character(label_source) ||
      length(label_source) != 1 ||
      !label_source %in% label_tiers
  ) {
    stop(
      "'label_source' must be one of: ",
      paste(label_tiers, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  candidates <- tibble::as_tibble(candidates)
  needed <- c(
    "candidate_id",
    "status",
    .sw_candidate_metric_cols,
    "evidence_comparison_ids"
  )
  missing_cols <- setdiff(needed, names(candidates))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) missing from 'candidates': ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(candidates$candidate_id) > 0) {
    stop("'candidate_id' must be unique.", call. = FALSE)
  }
  bad_status <- setdiff(unique(candidates$status), .sw_candidate_status)
  if (length(bad_status) > 0) {
    stop(
      "Unknown candidate status: ",
      paste(bad_status, collapse = ", "),
      ". Allowed: ",
      paste(.sw_candidate_status, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.list(candidates$evidence_comparison_ids)) {
    stop("'evidence_comparison_ids' must be a list column.", call. = FALSE)
  }
  if (!identity_id %in% candidates$candidate_id) {
    stop(
      "The identity candidate '",
      identity_id,
      "' is not in 'candidates'; it must always be assessed.",
      call. = FALSE
    )
  }
  if (!"status_reason" %in% names(candidates)) {
    candidates$status_reason <- NA_character_
  }
  candidates <- candidates[order(candidates$candidate_id), ]

  # Support class of each population-free batch contrast.
  comps <- design$comparisons[is.na(design$comparisons$population), ]
  support_by_comparison <- stats::setNames(comps$support, comps$comparison_id)
  all_non_identifiable <- nrow(comps) > 0 &&
    all(comps$support == "non_identifiable")

  # No-effect check on the identity candidate's worst marker.
  baseline_emd <- candidates$emd_worst_marker_median[
    candidates$candidate_id == identity_id
  ]
  no_effect <- !is.na(baseline_emd) && baseline_emd < policy$no_effect_emd

  rows <- lapply(seq_len(nrow(candidates)), function(i) {
    .sw_assess_one(
      candidates[i, ],
      design = design,
      policy = policy,
      identity_id = identity_id,
      label_source = label_source,
      support_by_comparison = support_by_comparison,
      all_non_identifiable = all_non_identifiable,
      no_effect = no_effect,
      baseline_emd = baseline_emd
    )
  })
  assessments <- dplyr::bind_rows(rows)

  recommended <- assessments$candidate_id[
    assessments$state == "recommended_under_policy"
  ]
  withheld <- NA_character_
  if (length(recommended) == 0) {
    withheld <- paste0(
      "no candidate cleared every gate under this policy; review the assessments and ",
      "their reasons. Retaining the uncorrected data is a valid outcome."
    )
  } else if (length(recommended) > 1) {
    # A tie is never silent: the report says why no single candidate is named.
    withheld <- sprintf(
      paste0(
        "%d candidates cleared every gate under this policy (%s); the policy does not ",
        "break ties between fully qualified candidates. Choose on grounds the policy ",
        "does not encode (cost, throughput, laboratory preference) and record the choice."
      ),
      length(recommended),
      paste(sort(recommended), collapse = ", ")
    )
  }

  # Pareto front over the largest group of candidates sharing the same evidence comparisons.
  warnings_out <- character(0)
  ok_rows <- candidates$status == "ok" &
    lengths(candidates$evidence_comparison_ids) > 0
  keys <- vapply(
    candidates$evidence_comparison_ids,
    function(ids) paste(sort(unique(as.character(ids))), collapse = "|"),
    character(1)
  )
  groups <- split(candidates$candidate_id[ok_rows], keys[ok_rows])
  front <- character(0)
  comparable <- FALSE
  if (length(groups) > 0) {
    largest <- groups[[which.max(lengths(groups))]]
    if (length(largest) > 1) {
      comparable <- TRUE
      front <- .sw_pareto_front(assessments[
        assessments$candidate_id %in% largest,
      ])
    } else {
      warnings_out <- c(
        warnings_out,
        paste0(
          "candidates were not evaluated on comparable evidence, so no Pareto comparison ",
          "is shown"
        )
      )
    }
  }

  if (policy$any_threshold_uncalibrated) {
    warnings_out <- c(
      warnings_out,
      paste0(
        "the policy thresholds in use are uncalibrated. Uncalibrated engineering ",
        "defaults are not validated universal cutoffs and must not be cited as such."
      )
    )
  }

  structure(
    list(
      assessments = assessments,
      recommended_candidate_id = if (length(recommended) == 1) {
        recommended
      } else {
        NA_character_
      },
      recommendation_withheld_reason = withheld,
      pareto_front = front,
      pareto_comparable = comparable,
      warnings = warnings_out,
      policy_settings = policy$settings,
      require_reference_supported = policy$require_reference_supported,
      identity_id = identity_id,
      label_source = label_source
    ),
    class = "sw_batch_decision"
  )
}

#' @export
print.sw_batch_decision <- function(x, ...) {
  cat(
    "<sw_batch_decision>",
    nrow(x$assessments),
    "candidate(s); labels:",
    x$label_source,
    "\n"
  )
  print(x$assessments[, c("candidate_id", "state", "reason")])
  if (is.na(x$recommended_candidate_id)) {
    cat("No recommendation:", x$recommendation_withheld_reason, "\n")
  } else {
    cat("Recommended under policy:", x$recommended_candidate_id, "\n")
  }
  if (x$pareto_comparable) {
    cat("Pareto front:", paste(x$pareto_front, collapse = ", "), "\n")
  }
  for (w in x$warnings) {
    cat("Warning:", w, "\n")
  }
  invisible(x)
}

#' Build decision candidates from long per-comparison metric records
#'
#' Summarises a long table of metric records into the one-row-per-candidate
#' table that \code{\link{sw_batch_decide}} takes. It also computes each
#' candidate's relative EMD reduction and coverage against the identity
#' candidate.
#'
#' The return shapes of \code{\link{sw_correct_evaluate_emd}} and
#' \code{\link{sw_correct_evaluate_mad}} come from cyCombine and carry no
#' per-comparison EMD, joint structure distance or compression fraction. This
#' helper therefore reads a documented long table instead. Ported from
#' batchdx-cyto (MIT).
#'
#' Summaries per candidate:
#' \describe{
#'   \item{\code{emd_median}}{Weighted median of the \code{"ok"}
#'     \code{emd_within_population} values.}
#'   \item{\code{emd_worst_marker_median}}{Median within each marker, then the
#'     largest. A batch effect confined to one marker is not averaged away.}
#'   \item{\code{emd_relative_reduction}}{(sum of identity EMD - sum of
#'     candidate EMD) / sum of identity EMD, over the keys (population,
#'     marker, comparison) that are \code{"ok"} for both. \code{NA} when no key
#'     is shared or the identity sum is at most 1e-9.}
#'   \item{\code{coverage}}{Shared keys divided by the identity candidate's
#'     \code{"ok"} keys; 0 when the identity has none, \code{NA} when the
#'     candidate has no EMD records.}
#'   \item{\code{mad_compression_fraction},
#'     \code{joint_structure_distance_median}}{Weighted medians of the
#'     \code{"ok"} values of those metrics.}
#'   \item{\code{evidence_comparison_ids}}{Comparison ids with an \code{"ok"}
#'     EMD record.}
#' }
#'
#' @param records A data frame with columns \code{candidate_id},
#'   \code{metric} (\code{"emd_within_population"},
#'   \code{"mad_compression_fraction"} or \code{"joint_structure_distance"}),
#'   \code{marker}, \code{comparison_id} (a population-free id from
#'   \code{sw_design_audit()$comparisons}), \code{value}, \code{status}
#'   (\code{"ok"} or any other string) and \code{weight} (aggregation weight
#'   for the median). An optional \code{population} column completes the
#'   record key.
#' @param identity_id The \code{candidate_id} of the uncorrected candidate.
#'   **No default.**
#' @param candidate_status Optional named character vector of run status per
#'   candidate (\code{"ok"}, \code{"failed"}, \code{"ineligible"},
#'   \code{"not_available"}). Candidates not named get \code{"ok"}.
#'
#' @return A tibble with one row per candidate, in the shape
#'   \code{\link{sw_batch_decide}} expects.
#'
#' @seealso \code{\link{sw_batch_decide}}
#'
#' @examples
#' records <- tibble::tibble(
#'   candidate_id = rep(c("uncorrected", "cycombine"), each = 3),
#'   metric = rep(c("emd_within_population", "emd_within_population",
#'                  "mad_compression_fraction"), 2),
#'   marker = rep(c("CD3", "CD4", NA), 2),
#'   comparison_id = "b1__vs__b2",
#'   value = c(1.0, 0.8, 0, 0.3, 0.2, 0.05),
#'   status = "ok",
#'   weight = 1
#' )
#' sw_batch_candidates_from_metrics(records, identity_id = "uncorrected")
#'
#' @export
sw_batch_candidates_from_metrics <- function(
  records,
  identity_id,
  candidate_status = NULL
) {
  if (missing(records) || missing(identity_id)) {
    stop(
      "'records' and 'identity_id' are required and have no default.",
      call. = FALSE
    )
  }
  records <- tibble::as_tibble(records)
  needed <- c(
    "candidate_id",
    "metric",
    "marker",
    "comparison_id",
    "value",
    "status",
    "weight"
  )
  missing_cols <- setdiff(needed, names(records))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) missing from 'records': ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!"population" %in% names(records)) {
    records$population <- NA_character_
  }
  records$key <- paste(
    records$population,
    records$marker,
    records$comparison_id,
    sep = "|"
  )
  if (anyDuplicated(records[, c("candidate_id", "metric", "key")]) > 0) {
    stop(
      "'records' has duplicate (candidate_id, metric, population, marker, comparison_id) rows.",
      call. = FALSE
    )
  }
  records$ok <- records$status == "ok" & !is.na(records$value)

  # Identity EMD keys that are usable as the reduction denominator.
  emd <- records[records$metric == "emd_within_population", ]
  identity_ok <- emd[emd$candidate_id == identity_id & emd$ok, ]
  identity_value <- stats::setNames(identity_ok$value, identity_ok$key)

  ids <- sort(unique(records$candidate_id))
  rows <- lapply(ids, function(id) {
    own <- records[records$candidate_id == id, ]
    own_emd <- own[own$metric == "emd_within_population", ]
    ok_emd <- own_emd[own_emd$ok, ]

    # Reduction and coverage on the keys OK for both this candidate and the identity.
    reduction <- NA_real_
    coverage <- NA_real_
    if (nrow(own_emd) > 0 && identity_id %in% records$candidate_id) {
      usable <- intersect(names(identity_value), ok_emd$key)
      coverage <- if (length(identity_value) > 0) {
        length(usable) / length(identity_value)
      } else {
        0
      }
      denominator <- sum(identity_value[usable])
      if (length(usable) > 0 && denominator > 1e-9) {
        candidate_sum <- sum(ok_emd$value[match(usable, ok_emd$key)])
        reduction <- (denominator - candidate_sum) / denominator
      }
    }

    tibble::tibble(
      candidate_id = id,
      status = if (
        !is.null(candidate_status) && id %in% names(candidate_status)
      ) {
        unname(candidate_status[[id]])
      } else {
        "ok"
      },
      emd_median = .sw_weighted_median(own, "emd_within_population"),
      emd_worst_marker_median = .sw_worst_marker_median(
        own,
        "emd_within_population"
      ),
      emd_relative_reduction = reduction,
      coverage = coverage,
      mad_compression_fraction = .sw_weighted_median(
        own,
        "mad_compression_fraction"
      ),
      joint_structure_distance_median = .sw_weighted_median(
        own,
        "joint_structure_distance"
      ),
      evidence_comparison_ids = list(sort(unique(ok_emd$comparison_id)))
    )
  })
  dplyr::bind_rows(rows)
}

.sw_candidate_metric_cols <- c(
  "emd_median",
  "emd_worst_marker_median",
  "emd_relative_reduction",
  "coverage",
  "mad_compression_fraction",
  "joint_structure_distance_median"
)

.sw_candidate_status <- c("ok", "failed", "ineligible", "not_available")

# One candidate through the gates, in the batchdx-cyto order.
.sw_assess_one <- function(
  cand,
  design,
  policy,
  identity_id,
  label_source,
  support_by_comparison,
  all_non_identifiable,
  no_effect,
  baseline_emd
) {
  id <- cand$candidate_id
  evidence_ids <- unique(as.character(cand$evidence_comparison_ids[[1]]))
  reference_supported <- length(evidence_ids) > 0 &&
    all(evidence_ids %in% names(support_by_comparison)) &&
    all(support_by_comparison[evidence_ids] == "reference_supported")
  reduction <- cand$emd_relative_reduction
  coverage <- cand$coverage
  compression <- cand$mad_compression_fraction
  joint <- cand$joint_structure_distance_median
  status_reason <- cand$status_reason
  guardrails <- character(0)

  if (cand$status == "failed") {
    state <- "execution_failed"
    reason <- if (is.na(status_reason)) "the backend failed" else status_reason
  } else if (cand$status %in% c("ineligible", "not_available")) {
    state <- "insufficient_evidence"
    reason <- if (is.na(status_reason)) {
      "the backend did not run"
    } else {
      status_reason
    }
  } else if (isTRUE(design$below_design_floor)) {
    state <- "below_design_floor"
    reason <- if (is.na(design$design_floor_reason)) {
      "the design is below the declared floor"
    } else {
      design$design_floor_reason
    }
  } else if (all_non_identifiable) {
    state <- "non_identifiable"
    reason <- "every batch contrast is non-identifiable given the protected factors"
  } else if (id == identity_id && no_effect) {
    state <- "not_needed"
    reason <- sprintf(
      paste0(
        "the worst marker's median within-population EMD on the uncorrected baseline ",
        "is %.3g, below the declared no-effect threshold %.3g; no marker shows a batch ",
        "difference worth correcting"
      ),
      baseline_emd,
      policy$no_effect_emd
    )
  } else if (id == identity_id) {
    state <- "eligible_for_review"
    reason <- "the identity candidate is always retained as a valid outcome"
  } else {
    # Guardrails are recorded even when a later design gate sets the state.
    if (
      !is.na(compression) && compression > policy$max_mad_compression_fraction
    ) {
      guardrails <- c(
        guardrails,
        sprintf(
          "mad_compression_fraction %.3g exceeds %.3g",
          compression,
          policy$max_mad_compression_fraction
        )
      )
    }
    if (!is.na(joint) && joint > policy$max_joint_structure_change) {
      guardrails <- c(
        guardrails,
        sprintf(
          "joint_structure_distance median %.3g exceeds %.3g",
          joint,
          policy$max_joint_structure_change
        )
      )
    }
    if (length(design$confounded_factors) > 0) {
      state <- "non_identifiable"
      reason <- sprintf(
        paste0(
          "protected factor(s) %s are perfectly confounded with the batch factor '%s'; ",
          "no automatic correction recommendation is possible for the affected contrast, ",
          "and adding the factor as a covariate does not help"
        ),
        paste(design$confounded_factors, collapse = ", "),
        design$batch_col
      )
    } else if (length(guardrails) > 0) {
      state <- "rejected_distortion"
      reason <- paste(guardrails, collapse = "; ")
    } else if (is.na(reduction)) {
      state <- "insufficient_evidence"
      reason <- "no relative reduction could be computed against the identity candidate"
    } else if (is.na(coverage)) {
      state <- "insufficient_evidence"
      reason <- "comparison coverage is unavailable for this candidate"
    } else if (coverage < policy$min_coverage) {
      state <- "insufficient_evidence"
      reason <- sprintf(
        "comparison coverage %.2f is below the declared minimum %.2f",
        coverage,
        policy$min_coverage
      )
    } else if (is.na(compression) || is.na(joint)) {
      # An unguarded candidate is insufficient evidence, not a silent pass.
      absent <- c(
        "mad_compression_fraction",
        "joint_structure_distance"
      )[is.na(c(compression, joint))]
      state <- "insufficient_evidence"
      reason <- sprintf(
        paste0(
          "preservation guardrail(s) %s have no usable record for this candidate, ",
          "so distortion cannot be ruled out"
        ),
        paste(absent, collapse = ", ")
      )
    } else if (reduction < policy$min_emd_reduction) {
      state <- "eligible_for_review"
      reason <- sprintf(
        "technical improvement %.3g is below the declared minimum %.3g",
        reduction,
        policy$min_emd_reduction
      )
    } else if (policy$require_reference_supported && !reference_supported) {
      state <- "eligible_for_review"
      reason <- paste0(
        "not every batch contrast this candidate's evidence rests on is supported by a ",
        "measured reference; the evidence is conditional on comparable biological ",
        "samples behaving identically"
      )
    } else if (label_source == "clustering") {
      # Invariant, not a policy option: clustering labels can never alone support a recommendation.
      state <- "eligible_for_review"
      reason <- paste0(
        "evaluation labels come from clustering of the uncorrected baseline and may ",
        "themselves be batch-driven; they cannot alone support a recommendation"
      )
    } else if (policy$any_threshold_uncalibrated) {
      state <- "eligible_for_review"
      reason <- paste0(
        "every gate passed, but the policy thresholds are uncalibrated, so the outcome ",
        "is review rather than automatic approval"
      )
    } else {
      state <- "recommended_under_policy"
      reason <- sprintf(
        "coverage %.2f, improvement %.3g, no guardrail triggered",
        coverage,
        reduction
      )
    }
  }

  tibble::tibble(
    candidate_id = id,
    state = state,
    reason = reason,
    guardrails_failed = list(guardrails),
    emd_median = cand$emd_median,
    emd_worst_marker_median = cand$emd_worst_marker_median,
    emd_relative_reduction = reduction,
    coverage = coverage,
    mad_compression_fraction = compression,
    joint_structure_distance_median = joint,
    reference_supported = reference_supported
  )
}

# Candidates not dominated on (reduction up, joint change down, compression down).
.sw_pareto_front <- function(assessments) {
  pts <- cbind(
    assessments$emd_relative_reduction,
    -assessments$joint_structure_distance_median,
    -assessments$mad_compression_fraction
  )
  keep <- stats::complete.cases(pts)
  pts <- pts[keep, , drop = FALSE]
  ids <- assessments$candidate_id[keep]
  dominated <- vapply(
    seq_along(ids),
    function(i) {
      any(vapply(
        seq_along(ids)[-i],
        function(j) all(pts[j, ] >= pts[i, ]) && any(pts[j, ] > pts[i, ]),
        logical(1)
      ))
    },
    logical(1)
  )
  sort(ids[!dominated])
}

# Weighted median of the OK values of one metric; NA when none.
.sw_weighted_median <- function(records, metric) {
  ok <- records[records$metric == metric & records$ok, ]
  if (nrow(ok) == 0) {
    return(NA_real_)
  }
  ord <- order(ok$value)
  values <- ok$value[ord]
  cumulative <- cumsum(ok$weight[ord])
  values[which(cumulative >= cumulative[length(cumulative)] / 2)[1]]
}

# Median within each marker, then the worst (largest) marker; NA when none.
.sw_worst_marker_median <- function(records, metric) {
  ok <- records[records$metric == metric & records$ok, ]
  if (nrow(ok) == 0) {
    return(NA_real_)
  }
  marker <- ifelse(is.na(ok$marker), "*", as.character(ok$marker))
  max(vapply(split(ok$value, marker), stats::median, numeric(1)))
}
