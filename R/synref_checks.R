# Synthetic references: overlap guard ----------------------------------------------

#' Stop when a synthetic reference shares events, files or donors with other data
#'
#' A synthetic reference is built from real events. If it shares events, files
#' or donors with the data it calibrates or is scored on, the score is
#' optimistic. The function reports the shared share per synthetic reference
#' and unit, and stops when any share is above `max_overlap_share`. The donor
#' is the dependence unit: files of one donor are not independent, and cells
#' are not replicates.
#'
#' @param synref An `sw_synref`.
#' @param against Data frame with `sample_id` (the calibration or scoring
#'   files) and optionally `row` (event index within the file), or a character
#'   vector of `sample_id`s. Without `row`, the event unit is `not_assessed`;
#'   the file unit covers it.
#' @param manifest Data frame with `sample_id` and `donor`, or `NULL` for an
#'   explicit "donor not assessed".
#' @param max_overlap_share Largest allowed share of shared units, from 0 to 1
#'   (uncalibrated). Any share above it stops the run.
#'
#' @return A tibble with `synref_id`, `unit` (event, file, donor), `n_units`,
#'   `n_shared`, `share` and `status` (`ok`, `overlap` or `not_assessed`),
#'   returned only when no share is above the limit.
#'
#' @examples
#' set.seed(1)
#' d <- list(a = cbind(M1 = rnorm(100)), b = cbind(M1 = rnorm(100)))
#' man <- data.frame(sample_id = c("a", "b"), batch = "b1", donor = c("d1", "d2"))
#' ref <- sw_synref_downsample(
#'   d, man, groups = list(r = "a"), strata = "file", labels = NULL,
#'   n_events = 20, n_draws = 1, seed = 1
#' )
#' sw_synref_check_overlap(ref, against = "b", manifest = man, max_overlap_share = 0)
#'
#' @export
sw_synref_check_overlap <- function(
  synref,
  against,
  manifest,
  max_overlap_share
) {
  .sw_stop_missing(c(
    synref = missing(synref),
    against = missing(against),
    manifest = missing(manifest),
    max_overlap_share = missing(max_overlap_share)
  ))
  .sw_validate_synref(synref)
  if (
    !is.numeric(max_overlap_share) ||
      length(max_overlap_share) != 1 ||
      is.na(max_overlap_share) ||
      max_overlap_share < 0 ||
      max_overlap_share > 1
  ) {
    stop("'max_overlap_share' must be one number from 0 to 1.", call. = FALSE)
  }
  if (is.character(against)) {
    against <- data.frame(sample_id = against)
  }
  if (!is.data.frame(against) || !"sample_id" %in% names(against)) {
    stop("'against' needs a 'sample_id' column.", call. = FALSE)
  }
  against_ids <- as.character(against[["sample_id"]])
  against_donors <- NULL
  if (!is.null(manifest)) {
    if (!is.data.frame(manifest) || !"sample_id" %in% names(manifest)) {
      stop("'manifest' needs a 'sample_id' column, or use NULL.", call. = FALSE)
    }
    if ("donor" %in% names(manifest)) {
      pos <- match(against_ids, as.character(manifest[["sample_id"]]))
      against_donors <- unique(stats::na.omit(as.character(manifest[["donor"]])[
        pos
      ]))
    }
  }
  against_keys <- if ("row" %in% names(against)) {
    paste(against_ids, against[["row"]], sep = "\r")
  }

  prov <- synref[["provenance"]]
  rows <- lapply(unique(prov[["synref_id"]]), function(id) {
    p <- prov[prov[["synref_id"]] == id, , drop = FALSE]
    files <- unique(p[["source_sample"]])
    donors <- unique(stats::na.omit(p[["donor"]]))
    n_event <- NA_integer_
    s_event <- NA_integer_
    if (!is.null(against_keys)) {
      keys <- paste(p[["source_sample"]], p[["source_row"]], sep = "\r")
      n_event <- nrow(p)
      s_event <- sum(keys %in% against_keys)
    }
    n_donor <- NA_integer_
    s_donor <- NA_integer_
    if (!is.null(against_donors) && length(donors) > 0) {
      n_donor <- length(donors)
      s_donor <- sum(donors %in% against_donors)
    }
    tibble::tibble(
      synref_id = id,
      unit = c("event", "file", "donor"),
      n_units = as.integer(c(n_event, length(files), n_donor)),
      n_shared = as.integer(c(s_event, sum(files %in% against_ids), s_donor))
    )
  })
  out <- dplyr::bind_rows(rows) |>
    dplyr::mutate(
      share = .data[["n_shared"]] / .data[["n_units"]],
      status = dplyr::case_when(
        is.na(.data[["share"]]) ~ "not_assessed",
        .data[["share"]] > max_overlap_share ~ "overlap",
        .default = "ok"
      )
    )
  bad <- out[out[["status"]] == "overlap", , drop = FALSE]
  if (nrow(bad) > 0) {
    stop(
      sprintf(
        "Synthetic reference shares data with 'against' above max_overlap_share = %s: %s.",
        format(max_overlap_share),
        paste(
          sprintf(
            "%s (%s, share %.3f)",
            bad[["synref_id"]],
            bad[["unit"]],
            bad[["share"]]
          ),
          collapse = "; "
        )
      ),
      call. = FALSE
    )
  }
  out
}
