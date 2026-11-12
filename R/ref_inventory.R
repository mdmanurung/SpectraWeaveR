# R/ref_inventory.R
# Reference checks, part 1: what reference material exists in each batch.
# Also holds the shared selection helpers used by sw_ref_consistency() and
# sw_ref_noise_floor(). Own code; no upstream code is used.

# Shared helpers ------------------------------------------------------------------

# Check `data`, `manifest` and `roles`, then return the manifest rows of the
# measured reference samples (columns sample_id, batch, reference_id and, when
# the manifest has them, aliquot and donor) and the event matrices of those
# samples. Stops when a selected sample has no data.
.sw_ref_select <- function(data, manifest, roles) {
  if (
    !is.list(data) ||
      is.data.frame(data) ||
      length(data) == 0 ||
      is.null(names(data)) ||
      anyNA(names(data)) ||
      anyDuplicated(names(data))
  ) {
    stop("`data` must be a uniquely named list of files.", call. = FALSE)
  }
  need <- c("sample_id", "batch", "reference_role", "reference_id")
  if (!is.data.frame(manifest) || !all(need %in% names(manifest))) {
    stop(
      "`manifest` needs columns ",
      paste0("'", need, "'", collapse = ", "),
      " (see sw_sample_manifest()).",
      call. = FALSE
    )
  }
  if (anyDuplicated(manifest[["sample_id"]])) {
    stop("`manifest$sample_id` must be unique.", call. = FALSE)
  }
  if (!is.character(roles) || length(roles) == 0 || anyNA(roles)) {
    stop(
      "`roles` must be a character vector of reference_role values, ",
      "for example c(\"technical_replicate\", \"anchor\").",
      call. = FALSE
    )
  }
  keep_cols <- intersect(
    c("sample_id", "batch", "reference_id", "aliquot", "donor"),
    names(manifest)
  )
  man <- as.data.frame(manifest)
  is_ref <- as.character(man[["reference_role"]]) %in%
    roles &
    !is.na(man[["reference_id"]])
  sel <- tibble::as_tibble(man[is_ref, keep_cols, drop = FALSE])
  sel[["sample_id"]] <- as.character(sel[["sample_id"]])
  sel[["batch"]] <- as.character(sel[["batch"]])
  sel[["reference_id"]] <- as.character(sel[["reference_id"]])
  if (anyNA(sel[["batch"]])) {
    stop("Reference samples have missing `batch` labels.", call. = FALSE)
  }
  if (nrow(sel) == 0) {
    stop(
      "No manifest row has `reference_role` in `roles` and a `reference_id`.",
      call. = FALSE
    )
  }
  absent <- setdiff(sel[["sample_id"]], names(data))
  if (length(absent) > 0) {
    stop(
      "Reference sample(s) not in `data`: ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  sel <- sel[order(sel[["reference_id"]], sel[["batch"]], sel[["sample_id"]]), ]
  list(sel = sel, data = data[sel[["sample_id"]]])
}

# Event matrices with the chosen markers; stops on a missing marker.
.sw_ref_matrices <- function(data, markers) {
  if (!is.character(markers) || length(markers) == 0 || anyNA(markers)) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  lapply(names(data), function(f) {
    m <- as.matrix(data[[f]])
    if (!is.numeric(m) || !all(markers %in% colnames(m))) {
      stop(
        "File '",
        f,
        "' is not numeric or lacks marker(s): ",
        paste(setdiff(markers, colnames(m)), collapse = ", "),
        call. = FALSE
      )
    }
    m[, markers, drop = FALSE]
  }) |>
    stats::setNames(names(data))
}

.sw_ref_check_seed <- function(seed) {
  if (
    !is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed != round(seed)
  ) {
    stop("`seed` must be one whole number.", call. = FALSE)
  }
}

# Inventory -----------------------------------------------------------------------

#' List the reference samples in each batch
#'
#' @description
#' Count what reference material each batch holds: events per sample, which
#' markers are present, and (with labels) events per population. The table has
#' counts only. It never says a reference is good or bad. `under_min` marks a
#' count below your `min_events`, so thin rows are visible before any
#' score is read.
#'
#' @details
#' `sw_design_audit()` never reads events, so it cannot see a reference file that
#' is nearly empty or lacks a channel. This function reads the event counts.
#' Only samples whose manifest `reference_role` is in `roles` and that carry a
#' `reference_id` are listed.
#'
#' Rows come in two kinds, named in `row_type`. A `"marker"` row has one marker
#' and `present` (the column exists and has at least one finite value). A
#' `"population"` row has one population and its event count; a population a
#' sample does not contain gets a count of 0, so a missing population shows up.
#'
#' @param data A named list of event matrices or data frames (one per file).
#' @param manifest A manifest with `sample_id`, `batch`, `reference_role` and
#'   `reference_id` (see [sw_sample_manifest()]).
#' @param markers Character vector of marker (column) names to check.
#' @param labels `NULL` (no population rows, explicit skip), or a named list
#'   with one vector of population labels per reference file (same length as
#'   its rows).
#' @param min_events Whole number, events below which a sample or population
#'   count is marked `under_min` (uncalibrated; you set it).
#' @param roles Character vector of `reference_role` values that count as
#'   reference material (uncalibrated; you set it).
#'
#' @return A tibble of class `sw_ref_inventory` with columns `reference_id`,
#'   `batch`, `sample_id`, `n_events`, `row_type`, `marker`, `present`,
#'   `population`, `n_population_events`, `under_min`. Attribute `settings`
#'   echoes `min_events` and `roles`.
#'
#' @seealso [sw_ref_consistency()], [sw_design_audit()]
#'
#' @examples
#' set.seed(1)
#' d <- list(
#'   r1 = cbind(A = stats::rnorm(80), B = stats::rnorm(80)),
#'   r2 = cbind(A = stats::rnorm(30), B = stats::rnorm(30))
#' )
#' man <- data.frame(
#'   sample_id = c("r1", "r2"),
#'   batch = c("b1", "b2"),
#'   reference_role = "technical_replicate",
#'   reference_id = "ref1"
#' )
#' sw_ref_inventory(
#'   d,
#'   man,
#'   markers = c("A", "B"),
#'   labels = NULL,
#'   min_events = 50, # illustrative, uncalibrated
#'   roles = "technical_replicate"
#' )
#'
#' @export
sw_ref_inventory <- function(
  data,
  manifest,
  markers,
  labels,
  min_events,
  roles
) {
  .sw_stop_missing(c(
    data = missing(data),
    manifest = missing(manifest),
    markers = missing(markers),
    labels = missing(labels),
    min_events = missing(min_events),
    roles = missing(roles)
  ))
  .sw_check_count(min_events, "min_events")
  picked <- .sw_ref_select(data, manifest, roles)
  sel <- picked[["sel"]]
  mats <- picked[["data"]]
  if (!is.character(markers) || length(markers) == 0 || anyNA(markers)) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.null(labels)) {
    ok <- is.list(labels) &&
      !is.null(names(labels)) &&
      all(sel[["sample_id"]] %in% names(labels)) &&
      all(vapply(
        sel[["sample_id"]],
        \(f) length(labels[[f]]) == NROW(mats[[f]]),
        logical(1)
      ))
    if (!ok) {
      stop(
        "`labels` must be NULL or a named list with one label per event for ",
        "every reference file.",
        call. = FALSE
      )
    }
  }

  # One block of marker rows per sample.
  marker_rows <- lapply(sel[["sample_id"]], \(f) {
    m <- as.matrix(mats[[f]])
    tibble::tibble(
      sample_id = f,
      n_events = nrow(m),
      row_type = "marker",
      marker = markers,
      present = vapply(
        markers,
        \(mk) mk %in% colnames(m) && any(is.finite(m[, mk])),
        logical(1),
        USE.NAMES = FALSE
      ),
      population = NA_character_,
      n_population_events = NA_integer_
    )
  }) |>
    dplyr::bind_rows()

  # Population counts, with 0 for a population this sample lacks.
  pop_rows <- NULL
  if (!is.null(labels)) {
    all_pops <- sel[["sample_id"]] |>
      lapply(\(f) as.character(labels[[f]])) |>
      unlist() |>
      unique() |>
      sort()
    pop_rows <- lapply(sel[["sample_id"]], \(f) {
      lab <- as.character(labels[[f]])
      tibble::tibble(
        sample_id = f,
        n_events = length(lab),
        row_type = "population",
        marker = NA_character_,
        present = NA,
        population = all_pops,
        n_population_events = vapply(
          all_pops,
          \(p) sum(lab == p, na.rm = TRUE),
          integer(1),
          USE.NAMES = FALSE
        )
      )
    }) |>
      dplyr::bind_rows()
  }

  out <- dplyr::bind_rows(marker_rows, pop_rows) |>
    dplyr::left_join(
      sel[c("sample_id", "batch", "reference_id")],
      by = "sample_id",
      relationship = "many-to-one",
      unmatched = "error"
    ) |>
    dplyr::mutate(
      under_min = dplyr::if_else(
        row_type == "marker",
        n_events < min_events,
        n_population_events < min_events
      )
    ) |>
    dplyr::select(
      reference_id,
      batch,
      sample_id,
      n_events,
      row_type,
      marker,
      present,
      population,
      n_population_events,
      under_min
    ) |>
    dplyr::arrange(reference_id, batch, sample_id, row_type)
  structure(
    out,
    settings = list(min_events = min_events, roles = roles),
    class = c("sw_ref_inventory", class(out))
  )
}
