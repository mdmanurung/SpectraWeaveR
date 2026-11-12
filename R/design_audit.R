#' Audit a batch-correction design before correcting
#'
#' Classifies every pair of batches by what the sample-level design can
#' support, and reports which correction and detection methods are eligible.
#' The audit works on samples, never on cells, because cells from one donor
#' are not independent replicates.
#'
#' Ported from the design audit of batchdx-cyto (MIT). Each pair of batches
#' gets one support class:
#' \describe{
#'   \item{\code{reference_supported}}{A measured reference (technical
#'     replicate or anchor) is present in both batches and no protected factor
#'     is disjoint between them.}
#'   \item{\code{biologically_stratified_conditional}}{No shared reference, but
#'     a protected factor has overlapping levels in both batches.}
#'   \item{\code{descriptive_only}}{Neither a shared reference nor an
#'     overlapping protected factor exists.}
#'   \item{\code{non_identifiable}}{A protected factor takes disjoint levels in
#'     the two batches, so the batch and biological contrasts are the same
#'     contrast. A shared reference and a covariate adjustment do not change
#'     this.}
#' }
#' A reference with role \code{"biological_pseudo_control"} or
#' \code{"synthetic_reference"} never counts as a measured technical reference
#' and never bridges batches. A \code{synthetic_reference} row needs a
#' \code{source_batch} column value (the batch its cells came from).
#'
#' @param samples A data frame with one row per sample. It needs a
#'   \code{sample_id} column, the \code{batch_col} column, and one column for
#'   each name in \code{protected}. For reference-based support it also needs
#'   \code{reference_role} and \code{reference_id} columns.
#' @param minimum_batches Integer. The design floor: fewer batches than this
#'   sets \code{below_design_floor}. **No default.** This choice changes which
#'   designs pass, so the analyst must set it and record the reason.
#' @param minimum_batches_origin One string: why that floor (a protocol, a
#'   paper, or "analyst choice"). It is stored in the result. **No default.**
#' @param require_shared_reference Logical. If \code{TRUE}, a design with no
#'   reference-supported batch pair sets \code{below_design_floor}. **No
#'   default**, for the same reason as \code{minimum_batches}.
#' @param batch_col Name of the correction factor column (default:
#'   \code{"batch"}).
#' @param protected Character vector of biological factor columns that must be
#'   preserved (for example \code{"condition"}). Default: none.
#' @param nested_within Optional column name (for example \code{"study"}).
#'   The audit warns when a batch spans more than one level. The column is
#'   kept out of the rank check because nested batches make it collinear.
#'   Kept for backward compatibility; \code{study_col} runs the same check.
#' @param donor_col Optional donor column. If absent from \code{samples}, the
#'   audit warns that donor-level uncertainty is unavailable. When present it
#'   fills \code{donor} in \code{reference_graph} and \code{donors_per_batch}.
#' @param study_col Optional study column. Each batch must sit in one study;
#'   the audit warns and names the batches that span more than one. Like
#'   \code{nested_within}, it stays out of the rank check.
#' @param run_col Optional acquisition-run column. Each run must sit in one
#'   batch; the audit warns and names the runs that span more than one.
#' @param aliquot_col Optional aliquot column. It fills \code{aliquot} in
#'   \code{reference_graph} and \code{aliquots_per_donor}.
#' @param anchor Optional anchor batch, one value of \code{batch_col}. It
#'   must be present in the data. **No default**: the anchor sets the target
#'   scale, so an automatic choice would turn its noise into bias in every
#'   other batch. Without it the \code{paired_shift} engine is ineligible.
#' @param crossed Character vector of crossed technical factor columns (for
#'   example \code{"instrument"}). The audit reports a batch-by-level table
#'   and flags factors perfectly confounded with batch. They are not added to
#'   the design matrix. Default: none.
#' @param populations Optional character vector of cell populations. Pairs are
#'   re-classified within each population.
#' @param population_counts Optional long data frame with columns
#'   \code{sample_id}, \code{population} and \code{n}. Required when
#'   \code{populations} is set.
#'
#' If \code{samples} is an \code{sw_manifest}, it is used as is. Its
#' \code{study}, \code{run}, \code{donor} and \code{aliquot} columns fill
#' \code{study_col}, \code{run_col}, \code{donor_col} and \code{aliquot_col}
#' when those arguments are \code{NULL}, because the manifest already maps
#' them. \code{anchor} and \code{crossed} are never filled.
#'
#' @return An object of class \code{sw_design_audit}: a list with
#'   \code{batches}, \code{n_batches}, \code{contingency},
#'   \code{design_rank}, \code{design_full_rank}, \code{confounded_factors},
#'   \code{reference_coverage}, \code{reference_components},
#'   \code{below_design_floor}, \code{design_floor_reason},
#'   \code{comparisons} (tibble), \code{eligibility} (tibble) and
#'   \code{warnings}. It also holds:
#'   \describe{
#'     \item{\code{reference_graph}}{Tibble with one row per measured
#'       reference observation group: \code{batch}, \code{reference_id},
#'       \code{donor}, \code{aliquot}, \code{n_samples}. \code{donor} and
#'       \code{aliquot} are \code{NA} when their column is not declared.}
#'     \item{\code{synthetic_references}}{Tibble of rows with role
#'       \code{synthetic_reference}, kept apart from measured references:
#'       \code{sample_id}, \code{batch}, \code{reference_id},
#'       \code{source_batch}. Zero rows when there are none. These rows
#'       appear in no coverage, component, graph or bridge table.}
#'     \item{\code{anchor}}{The declared anchor batch, or \code{NULL}.}
#'     \item{\code{anchor_bridge}}{Tibble with one row per non-anchor batch:
#'       \code{batch}, \code{n_shared_with_anchor}, \code{shared_ids}
#'       (list), \code{bridged} (at least one reference id shared with the
#'       anchor batch itself) and \code{indirect_only} (connected to the
#'       anchor only through other batches). Indirect connection is not
#'       bridged. Empty when no anchor is declared.}
#'     \item{\code{n_studies}, \code{n_runs}}{Distinct studies and runs, or
#'       \code{NA} when the column is not declared.}
#'     \item{\code{donors_per_batch}}{Tibble: \code{batch},
#'       \code{n_donors}, \code{n_reference_ids}. \code{n_donors} is
#'       \code{NA} without a donor column.}
#'     \item{\code{aliquots_per_donor}}{Tibble: \code{donor},
#'       \code{n_aliquots}. Empty unless both donor and aliquot columns are
#'       declared.}
#'     \item{\code{crossed}, \code{crossed_contingency},
#'       \code{confounded_crossed}}{The declared crossed factors, their
#'       batch-by-level counts (same columns as \code{contingency}), and the
#'       ones perfectly confounded with batch. These stay out of
#'       \code{confounded_factors}, which holds protected factors only.}
#'   }
#'   The \code{eligibility} table has a \code{paired_shift} row. It is
#'   eligible only when an anchor is declared, at least 2 batches exist, and
#'   every non-anchor batch shares a measured reference id with the anchor.
#'
#' @seealso \code{\link{sw_batch_decide}}, \code{\link{sw_correct_run}},
#'   \code{\link{sw_flag_iqr}}
#'
#' @examples
#' samples <- data.frame(
#'   sample_id = paste0("s", 1:8),
#'   batch = rep(c("b1", "b2", "b3", "b4"), each = 2),
#'   condition = rep(c("ctrl", "case"), 4),
#'   reference_role = rep(c("anchor", "none"), 4),
#'   reference_id = rep(c("ref1", NA), 4)
#' )
#' audit <- sw_design_audit(
#'   samples,
#'   minimum_batches = 3,
#'   minimum_batches_origin = "analyst choice for this example",
#'   require_shared_reference = TRUE,
#'   protected = "condition"
#' )
#' audit$comparisons[, c("comparison_id", "support")]
#'
#' # Declare an anchor to check the paired-reference shift engine.
#' audit <- sw_design_audit(
#'   samples,
#'   minimum_batches = 3,
#'   minimum_batches_origin = "analyst choice for this example",
#'   require_shared_reference = TRUE,
#'   protected = "condition",
#'   anchor = "b1"
#' )
#' audit$anchor_bridge
#' audit$eligibility[audit$eligibility$name == "paired_shift", ]
#'
#' @export
sw_design_audit <- function(
  samples,
  minimum_batches,
  minimum_batches_origin,
  require_shared_reference,
  batch_col = "batch",
  protected = character(0),
  nested_within = NULL,
  donor_col = NULL,
  populations = NULL,
  population_counts = NULL,
  study_col = NULL,
  run_col = NULL,
  aliquot_col = NULL,
  anchor = NULL,
  crossed = character(0)
) {
  if (
    missing(minimum_batches) ||
      missing(minimum_batches_origin) ||
      missing(require_shared_reference)
  ) {
    stop(
      "'minimum_batches', 'minimum_batches_origin' and ",
      "'require_shared_reference' have no default; set all three and ",
      "record why.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(minimum_batches) ||
      length(minimum_batches) != 1 ||
      is.na(minimum_batches)
  ) {
    stop("'minimum_batches' must be a single number.", call. = FALSE)
  }
  if (
    !is.character(minimum_batches_origin) ||
      length(minimum_batches_origin) != 1 ||
      !nzchar(minimum_batches_origin)
  ) {
    stop(
      "'minimum_batches_origin' must be one non-empty string saying where ",
      "the floor comes from.",
      call. = FALSE
    )
  }
  if (
    !is.logical(require_shared_reference) ||
      length(require_shared_reference) != 1 ||
      is.na(require_shared_reference)
  ) {
    stop("'require_shared_reference' must be TRUE or FALSE.", call. = FALSE)
  }
  # A manifest already maps its standard columns, so use them when unset.
  if (inherits(samples, "sw_manifest")) {
    if (is.null(study_col) && "study" %in% names(samples)) {
      study_col <- "study"
    }
    if (is.null(run_col) && "run" %in% names(samples)) {
      run_col <- "run"
    }
    if (is.null(donor_col) && "donor" %in% names(samples)) {
      donor_col <- "donor"
    }
    if (is.null(aliquot_col) && "aliquot" %in% names(samples)) {
      aliquot_col <- "aliquot"
    }
  }
  if (
    !is.null(anchor) &&
      (length(anchor) != 1 || is.na(anchor) || !nzchar(as.character(anchor)))
  ) {
    stop("'anchor' must be one batch value.", call. = FALSE)
  }
  samples <- tibble::as_tibble(samples)
  needed <- c(
    "sample_id",
    batch_col,
    protected,
    nested_within,
    study_col,
    run_col,
    aliquot_col,
    crossed
  )
  absent <- setdiff(needed, names(samples))
  if (length(absent) > 0) {
    stop(
      "Column(s) missing from 'samples': ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(samples$sample_id) > 0) {
    stop("'sample_id' must be unique.", call. = FALSE)
  }
  if (!is.null(populations) && is.null(population_counts)) {
    stop(
      "'population_counts' is required when 'populations' is set.",
      call. = FALSE
    )
  }
  has_ref <- c("reference_role", "reference_id") %in% names(samples)
  if (any(has_ref) && !all(has_ref)) {
    stop(
      "'reference_role' and 'reference_id' must be given together.",
      call. = FALSE
    )
  }
  if (!any(has_ref)) {
    samples$reference_role <- "none"
    samples$reference_id <- NA_character_
  }
  bad_role <- setdiff(
    unique(stats::na.omit(as.character(samples$reference_role))),
    c(
      "none", "technical_replicate", "anchor", "biological_pseudo_control",
      "synthetic_reference"
    )
  )
  if (length(bad_role) > 0) {
    stop(
      "Unknown reference_role value(s): ",
      paste(bad_role, collapse = ", "),
      ". Use none, technical_replicate, anchor, biological_pseudo_control or ",
      "synthetic_reference.",
      call. = FALSE
    )
  }
  is_synthetic <- samples$reference_role %in% "synthetic_reference"
  if (
    any(is_synthetic) &&
      (!"source_batch" %in% names(samples) ||
        anyNA(samples$source_batch[is_synthetic]))
  ) {
    stop(
      "Samples with role 'synthetic_reference' need a 'source_batch'.",
      call. = FALSE
    )
  }

  samples$.batch <- as.character(samples[[batch_col]])
  batches <- sort(unique(stats::na.omit(samples$.batch)))
  n_batches <- length(batches)
  warnings_out <- character(0)
  if (!is.null(anchor)) {
    anchor <- as.character(anchor)
    if (!anchor %in% batches) {
      stop(
        sprintf("Anchor batch '%s' is not present in '%s'.", anchor, batch_col),
        call. = FALSE
      )
    }
  }

  n_missing <- sum(is.na(samples$.batch))
  if (n_missing > 0) {
    warnings_out <- c(
      warnings_out,
      sprintf(
        "%d sample(s) have no value for '%s' and are excluded from every batch contrast.",
        n_missing,
        batch_col
      )
    )
  }
  if (
    !is.null(donor_col) &&
      (!donor_col %in% names(samples) || all(is.na(samples[[donor_col]])))
  ) {
    warnings_out <- c(
      warnings_out,
      sprintf(
        "No sample has a value for '%s'; donor-matched metrics and donor-level uncertainty are unavailable.",
        donor_col
      )
    )
  }

  # Batch-by-level counts for each protected factor and each crossed factor.
  contingency <- .sw_batch_contingency(samples, protected)
  crossed_contingency <- .sw_batch_contingency(samples, crossed)

  # Treatment-coded design matrix rank over the correction and protected factors.
  design_cols <- c(".batch", protected)
  design_df <- samples[stats::complete.cases(samples[design_cols]), design_cols]
  design_df[] <- lapply(design_df, factor)
  usable <- vapply(design_df, nlevels, integer(1)) > 1
  rank <- 0L
  full_rank <- FALSE
  if (nrow(design_df) > 0) {
    if (any(usable)) {
      mm <- stats::model.matrix(
        stats::reformulate(names(design_df)[usable]),
        data = design_df
      )
    } else {
      mm <- matrix(1, nrow = nrow(design_df), ncol = 1)
    }
    rank <- qr(mm)$rank
    full_rank <- rank == ncol(mm)
  }

  # Batches must nest in studies; study_col runs the old nested_within check.
  for (nest_col in unique(c(nested_within, study_col))) {
    spans <- tapply(
      as.character(samples[[nest_col]]),
      samples$.batch,
      function(v) length(unique(stats::na.omit(v)))
    )
    multi <- names(spans)[spans > 1]
    if (length(multi) > 0) {
      warnings_out <- c(
        warnings_out,
        sprintf(
          "Batch(es) %s span more than one '%s' level, so batches are not nested in '%s'.",
          paste(multi, collapse = ", "),
          nest_col,
          nest_col
        )
      )
    }
  }

  # Runs must nest in batches.
  if (!is.null(run_col)) {
    d <- samples[!is.na(samples[[run_col]]) & !is.na(samples$.batch), ]
    spans <- tapply(d$.batch, as.character(d[[run_col]]), function(v) {
      length(unique(v))
    })
    multi <- names(spans)[spans > 1]
    if (length(multi) > 0) {
      warnings_out <- c(
        warnings_out,
        sprintf(
          "Run(s) %s span more than one '%s' level, so runs are not nested in batches.",
          paste(multi, collapse = ", "),
          batch_col
        )
      )
    }
  }
  n_distinct_or_na <- function(col) {
    if (is.null(col)) {
      return(NA_integer_)
    }
    length(unique(stats::na.omit(as.character(samples[[col]]))))
  }
  n_studies <- n_distinct_or_na(study_col)
  n_runs <- n_distinct_or_na(run_col)

  confounded <- Filter(
    function(f) .sw_perfectly_confounded(samples, f),
    protected
  )
  batch_identifiable <- length(confounded) == 0
  confounded_crossed <- Filter(
    function(f) .sw_perfectly_confounded(samples, f),
    crossed
  )

  refs <- .sw_reference_components(samples)
  n_with_ref <- sum(vapply(refs$coverage, length, integer(1)) > 0)
  has_donor <- !is.null(donor_col) && donor_col %in% names(samples)
  graph <- .sw_reference_graph(samples, if (has_donor) donor_col, aliquot_col)
  bridge <- .sw_anchor_bridge(anchor, batches, refs)
  counts <- .sw_donor_counts(
    samples,
    batches,
    refs$coverage,
    if (has_donor) donor_col,
    aliquot_col
  )

  # Classify every batch pair, and again within each population when requested.
  pairs <- if (n_batches >= 2) {
    utils::combn(batches, 2, simplify = FALSE)
  } else {
    list()
  }
  comparisons <- lapply(pairs, function(p) {
    rows <- list(.sw_classify_pair(
      p[1],
      p[2],
      samples,
      protected,
      refs$coverage
    ))
    for (pop in populations) {
      rows[[length(rows) + 1]] <- .sw_classify_pair(
        p[1],
        p[2],
        samples,
        protected,
        refs$coverage,
        population = pop,
        population_counts = population_counts
      )
    }
    dplyr::bind_rows(rows)
  }) |>
    dplyr::bind_rows()
  if (nrow(comparisons) == 0) {
    comparisons <- .sw_empty_comparisons()
  }

  floor_reasons <- character(0)
  if (n_batches < minimum_batches) {
    floor_reasons <- c(
      floor_reasons,
      sprintf(
        "%d batch(es) present; the declared minimum is %d.",
        n_batches,
        as.integer(minimum_batches)
      )
    )
  }
  any_supported <- any(
    comparisons$support[is.na(comparisons$population)] == "reference_supported"
  )
  if (isTRUE(require_shared_reference) && !any_supported) {
    floor_reasons <- c(
      floor_reasons,
      "No measured reference sample is shared across any pair of batches."
    )
  }

  eligibility <- lapply(seq_len(nrow(.sw_eligibility_rules)), function(i) {
    .sw_evaluate_rule(
      .sw_eligibility_rules[i, ],
      n_batches,
      n_with_ref,
      batch_identifiable
    )
  }) |>
    dplyr::bind_rows()
  ps_row <- .sw_evaluate_paired_shift(anchor, n_batches, bridge)
  # heb v0 uses the same rule (every batch shares a reference with the anchor);
  # a rule based on a connected reference graph is a later variant.
  heb_row <- ps_row
  heb_row$name <- "heb"
  heb_row$sw_function <- "sw_corrector_fit(engine = 'heb')"
  shared_row <- ps_row
  shared_row$name <- "heb_shared"
  shared_row$sw_function <- "sw_corrector_fit(engine = 'heb_shared')"
  eligibility <- dplyr::bind_rows(eligibility, ps_row, heb_row, shared_row)

  if (length(confounded) > 0) {
    warnings_out <- c(
      warnings_out,
      sprintf(
        paste0(
          "Protected factor(s) %s are perfectly confounded with '%s'. ",
          "No automatic correction recommendation is possible for the ",
          "affected contrast, and adding the factor as a covariate does not help."
        ),
        paste(confounded, collapse = ", "),
        batch_col
      )
    )
  }
  if (length(confounded_crossed) > 0) {
    warnings_out <- c(
      warnings_out,
      sprintf(
        paste0(
          "Crossed factor(s) %s are perfectly confounded with '%s'. ",
          "Their effect cannot be separated from the batch effect, so a ",
          "batch correction absorbs both."
        ),
        paste(confounded_crossed, collapse = ", "),
        batch_col
      )
    )
  }
  if (!is.null(anchor) && length(refs$coverage[[anchor]]) == 0) {
    warnings_out <- c(
      warnings_out,
      sprintf(
        "Anchor batch '%s' holds no measured reference sample, so no batch is bridged to it.",
        anchor
      )
    )
  }
  unbridged <- bridge$batch[!bridge$bridged]
  if (!is.null(anchor) && length(unbridged) > 0) {
    warnings_out <- c(
      warnings_out,
      sprintf(
        "Batch(es) %s share no measured reference with anchor '%s'.",
        paste(unbridged, collapse = ", "),
        anchor
      )
    )
  }
  if (length(refs$components) > 1 && any(lengths(refs$coverage) > 0)) {
    warnings_out <- c(
      warnings_out,
      sprintf(
        paste0(
          "The measured-reference graph has %d disconnected component(s); ",
          "correction cannot be assumed to propagate between them."
        ),
        length(refs$components)
      )
    )
  }
  if (!full_rank) {
    warnings_out <- c(
      warnings_out,
      sprintf(
        "Design matrix is rank-deficient (rank %d).",
        rank
      )
    )
  }

  structure(
    list(
      batch_col = batch_col,
      protected = protected,
      batches = batches,
      n_batches = n_batches,
      contingency = contingency,
      design_rank = rank,
      design_full_rank = full_rank,
      confounded_factors = confounded,
      reference_coverage = refs$coverage,
      reference_components = refs$components,
      reference_graph = graph,
      synthetic_references = .sw_synthetic_references(samples),
      anchor = anchor,
      anchor_bridge = bridge,
      n_studies = n_studies,
      n_runs = n_runs,
      donors_per_batch = counts$donors_per_batch,
      aliquots_per_donor = counts$aliquots_per_donor,
      crossed = crossed,
      crossed_contingency = crossed_contingency,
      confounded_crossed = confounded_crossed,
      minimum_batches = minimum_batches,
      minimum_batches_origin = minimum_batches_origin,
      below_design_floor = length(floor_reasons) > 0,
      design_floor_reason = if (length(floor_reasons) > 0) {
        paste(floor_reasons, collapse = " ")
      } else {
        NA_character_
      },
      comparisons = comparisons,
      eligibility = eligibility,
      warnings = warnings_out
    ),
    class = "sw_design_audit"
  )
}

#' @export
print.sw_design_audit <- function(x, ...) {
  cat(
    "<sw_design_audit>",
    x$n_batches,
    "batches;",
    nrow(x$comparisons),
    "comparison(s)\n"
  )
  cat("Support: ")
  print(table(x$comparisons$support))
  cat(
    "Below design floor:",
    x$below_design_floor,
    if (x$below_design_floor) paste0("(", x$design_floor_reason, ")") else "",
    "\n"
  )
  if (is.null(x$anchor)) {
    cat("Anchor: none declared\n")
  } else {
    cat(sprintf(
      "Anchor: %s; %d of %d other batch(es) bridged directly\n",
      x$anchor,
      sum(x$anchor_bridge$bridged),
      nrow(x$anchor_bridge)
    ))
  }
  for (w in x$warnings) {
    cat("Warning:", w, "\n")
  }
  invisible(x)
}

# Method constraints with a stated origin, so no threshold is a hidden constant.
.sw_eligibility_rules <- tibble::tibble(
  name = c(
    "cycombine_cluster",
    "cycombine_express",
    "iqr_flag",
    "none",
    "cycombine"
  ),
  sw_function = c(
    "sw_correct_detect_batch",
    "sw_correct_detect_batch",
    "sw_flag_iqr",
    "(identity)",
    "sw_correct_run"
  ),
  kind = c("detector", "detector", "detector", "backend", "backend"),
  min_batches = c(3L, 2L, 4L, 1L, 2L),
  max_batches = c(30L, NA, NA, NA, NA),
  needs_reference_every_batch = c(FALSE, FALSE, TRUE, FALSE, FALSE),
  needs_identifiable_batch = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  origin = c(
    paste0(
      "cyCombine::detect_batch_effect computes EMD and cluster-percentage ",
      "outlier analysis only when 3 <= n_batches <= 30; outside that range ",
      "it emits UMAPs alone."
    ),
    paste0(
      "cyCombine::detect_batch_effect_express stops with an error below 2 ",
      "batches."
    ),
    paste0(
      "CytoBatchFlagR flags a batch when a per-control per-marker value ",
      "falls beyond 1.5 IQR of the quartiles taken across batches. The rule ",
      "needs a control sample in every batch, and for n <= 3 values the ",
      "Tukey fence can never be exceeded, so the assessment is degenerate ",
      "below 4 batches."
    ),
    "The identity candidate is always eligible and is always benchmarked.",
    paste0(
      "cyCombine::correct_data runs ComBat per SOM node; the batch term must ",
      "be estimable given the covariates supplied (cyCombine:::check_confound)."
    )
  )
)

.sw_empty_comparisons <- function() {
  tibble::tibble(
    comparison_id = character(),
    batch_a = character(),
    batch_b = character(),
    population = character(),
    n_samples_a = integer(),
    n_samples_b = integer(),
    shared_reference_ids = list(),
    support = character(),
    reason = character(),
    overlapping_biology = list()
  )
}

# TRUE when every batch holds exactly one level of `factor` and the batches
# differ. The batch and factor terms then span the same subspace.
.sw_perfectly_confounded <- function(samples, factor) {
  d <- samples[!is.na(samples$.batch) & !is.na(samples[[factor]]), ]
  if (length(unique(d$.batch)) < 2) {
    return(FALSE)
  }
  levels_by_batch <- split(as.character(d[[factor]]), d$.batch)
  one_each <- all(vapply(
    levels_by_batch,
    function(v) length(unique(v)) == 1L,
    logical(1)
  ))
  one_each && length(unique(vapply(levels_by_batch, `[`, character(1), 1))) > 1
}

# Rows with role synthetic_reference, listed apart from measured references.
.sw_synthetic_references <- function(samples) {
  syn <- samples$reference_role %in% "synthetic_reference"
  tibble::tibble(
    sample_id = as.character(samples$sample_id[syn]),
    batch = as.character(samples$.batch[syn]),
    reference_id = as.character(samples$reference_id[syn]),
    source_batch = as.character(samples[["source_batch"]][syn])
  )
}

# Measured-reference coverage per batch and connected components of batches.
.sw_reference_components <- function(samples) {
  measured <- samples$reference_role %in%
    c("technical_replicate", "anchor") &
    !is.na(samples$reference_id) &
    !is.na(samples$.batch)
  batches <- sort(unique(stats::na.omit(samples$.batch)))
  coverage <- lapply(batches, function(b) {
    sort(unique(samples$reference_id[measured & samples$.batch == b]))
  })
  names(coverage) <- batches

  adjacency <- lapply(batches, function(b) character(0))
  names(adjacency) <- batches
  for (id in unique(samples$reference_id[measured])) {
    members <- sort(unique(samples$.batch[
      measured & samples$reference_id == id
    ]))
    for (a in members) {
      adjacency[[a]] <- union(adjacency[[a]], setdiff(members, a))
    }
  }

  seen <- character(0)
  components <- list()
  for (node in batches) {
    if (node %in% seen) {
      next
    }
    stack <- node
    seen <- c(seen, node)
    comp <- character(0)
    while (length(stack) > 0) {
      cur <- stack[length(stack)]
      stack <- stack[-length(stack)]
      comp <- c(comp, cur)
      new <- setdiff(sort(adjacency[[cur]]), seen)
      seen <- c(seen, new)
      stack <- c(stack, new)
    }
    components[[length(components) + 1]] <- sort(comp)
  }
  list(coverage = coverage, components = components)
}

.sw_classify_pair <- function(
  batch_a,
  batch_b,
  samples,
  protected,
  coverage,
  population = NULL,
  population_counts = NULL
) {
  in_a <- samples[!is.na(samples$.batch) & samples$.batch == batch_a, ]
  in_b <- samples[!is.na(samples$.batch) & samples$.batch == batch_b, ]

  if (!is.null(population) && !is.null(population_counts)) {
    present <- population_counts$sample_id[
      population_counts$population == population & population_counts$n > 0
    ]
    in_a <- in_a[in_a$sample_id %in% present, ]
    in_b <- in_b[in_b$sample_id %in% present, ]
  }

  shared <- sort(intersect(coverage[[batch_a]], coverage[[batch_b]]))
  out <- function(support, reason, overlapping = character(0)) {
    tibble::tibble(
      comparison_id = paste0(
        batch_a,
        "__vs__",
        batch_b,
        if (!is.null(population)) paste0("@", population) else ""
      ),
      batch_a = batch_a,
      batch_b = batch_b,
      population = if (is.null(population)) NA_character_ else population,
      n_samples_a = nrow(in_a),
      n_samples_b = nrow(in_b),
      shared_reference_ids = list(shared),
      support = support,
      reason = reason,
      overlapping_biology = list(overlapping)
    )
  }

  if (nrow(in_a) == 0 || nrow(in_b) == 0) {
    return(out(
      "descriptive_only",
      if (!is.null(population)) {
        sprintf(
          "Population '%s' has no cells in at least one batch of the pair; the contrast is not supported within this population.",
          population
        )
      } else {
        "One batch of the pair contains no samples."
      }
    ))
  }

  overlapping <- character(0)
  disjoint <- character(0)
  for (f in protected) {
    la <- unique(stats::na.omit(as.character(in_a[[f]])))
    lb <- unique(stats::na.omit(as.character(in_b[[f]])))
    both <- intersect(la, lb)
    if (length(both) > 0) {
      overlapping <- c(overlapping, paste0(f, "=", sort(both)))
    } else if (length(la) > 0 && length(lb) > 0) {
      disjoint <- c(disjoint, f)
    }
  }

  # Identifiability comes first: a shared reference never rescues a protected
  # factor that is disjoint between the two batches.
  if (length(disjoint) > 0) {
    reason <- sprintf(
      "Protected factor(s) %s take disjoint levels across the two batches; the batch contrast and the biological contrast are not separately estimable. Adding the same factor as a covariate does not resolve this.",
      paste(disjoint, collapse = ", ")
    )
    if (length(shared) > 0) {
      reason <- sprintf(
        "Measured reference(s) %s are present in both batches, but %s The shared reference does not change that.",
        paste(shared, collapse = ", "),
        sub("^Protected", "protected", reason)
      )
    }
    return(out("non_identifiable", reason, overlapping))
  }
  if (length(shared) > 0) {
    return(out(
      "reference_supported",
      sprintf(
        "Measured reference identity/identities %s present in both batches.",
        paste(shared, collapse = ", ")
      ),
      overlapping
    ))
  }
  if (length(overlapping) > 0) {
    return(out(
      "biologically_stratified_conditional",
      sprintf(
        "No measured reference is shared; evidence is conditional on comparable biological samples (%s) behaving identically.",
        paste(overlapping, collapse = ", ")
      ),
      overlapping
    ))
  }
  out(
    "descriptive_only",
    "No measured reference and no declared protected factor to stratify on."
  )
}

.sw_evaluate_rule <- function(rule, n_batches, n_with_ref, batch_identifiable) {
  failures <- character(0)
  codes <- character(0)
  if (!is.na(rule$min_batches) && n_batches < rule$min_batches) {
    codes <- c(codes, "batch_count_below_minimum")
    failures <- c(
      failures,
      sprintf("%d batch(es) present, %d required", n_batches, rule$min_batches)
    )
  }
  if (!is.na(rule$max_batches) && n_batches > rule$max_batches) {
    codes <- c(codes, "batch_count_above_maximum")
    failures <- c(
      failures,
      sprintf(
        "%d batch(es) present, at most %d supported",
        n_batches,
        rule$max_batches
      )
    )
  }
  if (
    rule$needs_reference_every_batch &&
      !(n_batches > 0 && n_with_ref == n_batches)
  ) {
    codes <- c(codes, "missing_reference_in_some_batch")
    failures <- c(
      failures,
      sprintf(
        "%d of %d batch(es) contain a measured reference sample",
        n_with_ref,
        n_batches
      )
    )
  }
  if (rule$needs_identifiable_batch && !batch_identifiable) {
    codes <- c(codes, "batch_not_identifiable")
    failures <- c(
      failures,
      "at least one protected factor is perfectly confounded with the batch factor"
    )
  }

  eligible <- length(failures) == 0
  tibble::tibble(
    name = rule$name,
    sw_function = rule$sw_function,
    kind = rule$kind,
    eligible = eligible,
    reason_code = if (eligible) "eligible" else paste(codes, collapse = "+"),
    reason = if (eligible) {
      paste0("All declared requirements met. Constraint origin: ", rule$origin)
    } else {
      paste0(
        paste(failures, collapse = "; "),
        ". Constraint origin: ",
        rule$origin
      )
    },
    override_allowed = !"batch_not_identifiable" %in% codes
  )
}

# Batch-by-level counts for each named factor column, in long format.
.sw_batch_contingency <- function(samples, factors) {
  out <- lapply(factors, function(f) {
    d <- samples[!is.na(samples$.batch), ]
    counts <- as.data.frame(
      table(batch = d$.batch, level = as.character(d[[f]]), useNA = "ifany"),
      stringsAsFactors = FALSE
    )
    tibble::tibble(
      factor = f,
      batch = counts$batch,
      level = counts$level,
      n = as.integer(counts$Freq)
    )
  }) |>
    dplyr::bind_rows()
  if (nrow(out) == 0) {
    out <- tibble::tibble(
      factor = character(),
      batch = character(),
      level = character(),
      n = integer()
    )
  }
  out
}

# Measured reference observations grouped by (batch, reference_id, donor,
# aliquot). Missing donor or aliquot columns give NA, not an invented level.
.sw_reference_graph <- function(samples, donor_col, aliquot_col) {
  measured <- samples$reference_role %in%
    c("technical_replicate", "anchor") &
    !is.na(samples$reference_id) &
    !is.na(samples$.batch)
  d <- samples[measured, ]
  col_or_na <- function(col) {
    if (is.null(col)) rep(NA_character_, nrow(d)) else as.character(d[[col]])
  }
  obs <- tibble::tibble(
    batch = d$.batch,
    reference_id = as.character(d$reference_id),
    donor = col_or_na(donor_col),
    aliquot = col_or_na(aliquot_col)
  )
  key <- paste(obs$batch, obs$reference_id, obs$donor, obs$aliquot, sep = "\r")
  first <- !duplicated(key)
  out <- obs[first, ]
  out$n_samples <- as.integer(table(key)[key[first]])
  out[order(out$batch, out$reference_id, out$donor, out$aliquot), ]
}

# Reference ids each non-anchor batch shares with the anchor batch itself.
# A batch reached only through other batches is indirect_only, not bridged,
# because the paired shift averages (batch - anchor) over shared references.
.sw_anchor_bridge <- function(anchor, batches, refs) {
  if (is.null(anchor)) {
    return(tibble::tibble(
      batch = character(),
      n_shared_with_anchor = integer(),
      shared_ids = list(),
      bridged = logical(),
      indirect_only = logical()
    ))
  }
  others <- setdiff(batches, anchor)
  anchor_component <- Filter(function(comp) anchor %in% comp, refs$components)[[
    1
  ]]
  shared <- lapply(others, function(b) {
    sort(intersect(refs$coverage[[b]], refs$coverage[[anchor]]))
  })
  bridged <- lengths(shared) > 0
  tibble::tibble(
    batch = others,
    n_shared_with_anchor = as.integer(lengths(shared)),
    shared_ids = shared,
    bridged = bridged,
    indirect_only = !bridged & others %in% anchor_component
  )
}

# Donors and reference ids per batch; aliquots per donor.
.sw_donor_counts <- function(
  samples,
  batches,
  coverage,
  donor_col,
  aliquot_col
) {
  d <- samples[!is.na(samples$.batch), ]
  n_donors <- vapply(
    batches,
    function(b) {
      if (is.null(donor_col)) {
        return(NA_integer_)
      }
      length(unique(stats::na.omit(as.character(d[[donor_col]][
        d$.batch == b
      ]))))
    },
    integer(1)
  )
  donors_per_batch <- tibble::tibble(
    batch = batches,
    n_donors = unname(n_donors),
    n_reference_ids = as.integer(lengths(coverage[batches]))
  )
  aliquots_per_donor <- tibble::tibble(
    donor = character(),
    n_aliquots = integer()
  )
  if (!is.null(donor_col) && !is.null(aliquot_col)) {
    keep <- !is.na(d[[donor_col]]) & !is.na(d[[aliquot_col]])
    by_donor <- split(
      as.character(d[[aliquot_col]][keep]),
      as.character(d[[donor_col]][keep])
    )
    aliquots_per_donor <- tibble::tibble(
      donor = names(by_donor),
      n_aliquots = vapply(
        by_donor,
        function(v) length(unique(v)),
        integer(1),
        USE.NAMES = FALSE
      )
    )
  }
  list(
    donors_per_batch = donors_per_batch,
    aliquots_per_donor = aliquots_per_donor
  )
}

# Eligibility of the paired-reference shift. It depends on the anchor and the
# bridge table, so it is evaluated here, not by the generic rule columns.
.sw_paired_shift_origin <- paste0(
  "The paired-reference shift estimates each batch effect as the mean of ",
  "(batch - anchor) over reference ids measured in both that batch and the ",
  "anchor batch. A batch without a reference shared with the anchor batch ",
  "has no such difference; a path through other batches is not used. It ",
  "needs at least 2 batches. A missing anchor or bridge cannot be overridden."
)

.sw_evaluate_paired_shift <- function(anchor, n_batches, bridge) {
  codes <- character(0)
  failures <- character(0)
  if (n_batches < 2) {
    codes <- c(codes, "batch_count_below_minimum")
    failures <- c(
      failures,
      sprintf("%d batch(es) present, 2 required", n_batches)
    )
  }
  if (is.null(anchor)) {
    codes <- c(codes, "anchor_not_declared")
    failures <- c(failures, "no anchor batch was declared")
  } else if (any(!bridge$bridged)) {
    codes <- c(codes, "batch_without_anchor_bridge")
    failures <- c(
      failures,
      sprintf(
        "batch(es) %s share no measured reference id with anchor '%s'",
        paste(bridge$batch[!bridge$bridged], collapse = ", "),
        anchor
      )
    )
  }
  eligible <- length(failures) == 0
  tibble::tibble(
    name = "paired_shift",
    sw_function = "sw_corrector_fit(engine = 'paired_shift')",
    kind = "backend",
    eligible = eligible,
    reason_code = if (eligible) "eligible" else paste(codes, collapse = "+"),
    reason = if (eligible) {
      paste0(
        "All declared requirements met. Constraint origin: ",
        .sw_paired_shift_origin
      )
    } else {
      paste0(
        paste(failures, collapse = "; "),
        ". Constraint origin: ",
        .sw_paired_shift_origin
      )
    },
    override_allowed = !any(
      codes %in% c("anchor_not_declared", "batch_without_anchor_bridge")
    )
  )
}
