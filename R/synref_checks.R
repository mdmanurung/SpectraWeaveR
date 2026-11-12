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
#'   explicit "donor not assessed". Every `against` sample needs a row with a
#'   donor; otherwise the function stops.
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
      donor_of <- as.character(manifest[["donor"]])[pos]
      if (anyNA(donor_of)) {
        stop(
          "'against' sample(s) missing from 'manifest' or without a donor: ",
          paste(unique(against_ids[is.na(donor_of)]), collapse = ", "),
          ". The donor is the dependence unit; add them, or use manifest = NULL.",
          call. = FALSE
        )
      }
      against_donors <- unique(donor_of)
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

# Spread across draws ---------------------------------------------------------------

#' Spread of marker quantiles across draws of a synthetic reference
#'
#' For each reference group, marker and quantile, reports how far the quantile
#' moves between the draws of an `sw_synref` (standard deviation and range).
#' The draws come from one recipe repeated `n_draws` times (see
#' [sw_synref_downsample()] and [sw_synref_mix()]), so the spread is the
#' event-level sampling noise of that recipe. It shrinks as `n_events` grows.
#' It is not aliquot-to-aliquot or run-to-run noise. A group with fewer than
#' two draws gets `not_assessed`. Rows are not rolled up over markers.
#'
#' @param synref An `sw_synref` with draws named `<group>_dNN` (from
#'   [sw_synref_downsample()] or [sw_synref_mix()]) or `<group>_partNN`.
#' @param probs Quantile probabilities in [0, 1] (uncalibrated).
#'
#' @return A tibble: `group`, `marker`, `prob`, `sd_across_draws`,
#'   `range_across_draws`, `n_draws`, `mean_events` and `status` (`ok` or
#'   `not_assessed`).
#'
#' @examples
#' set.seed(1)
#' d <- list(a = cbind(M1 = rnorm(500), M2 = rnorm(500)))
#' man <- data.frame(sample_id = "a", batch = "b1")
#' ref <- sw_synref_downsample(
#'   d, man, groups = list(r = "a"), strata = "file", labels = NULL,
#'   n_events = 100, n_draws = 5, seed = 1
#' )
#' sw_synref_spread(ref, probs = c(0.5, 0.9))
#'
#' @export
sw_synref_spread <- function(synref, probs) {
  .sw_stop_missing(c(synref = missing(synref), probs = missing(probs)))
  .sw_validate_synref(synref)
  if (
    !is.numeric(probs) ||
      length(probs) == 0 ||
      anyNA(probs) ||
      any(probs < 0 | probs > 1) ||
      anyDuplicated(probs)
  ) {
    stop("'probs' must be unique numbers from 0 to 1.", call. = FALSE)
  }
  ids <- names(synref[["events"]])
  group_of <- sub("_(d|part)[0-9]+$", "", ids)
  rows <- lapply(unique(group_of), function(g) {
    mats <- synref[["events"]][ids[group_of == g]]
    markers <- colnames(mats[[1]])
    n_draws <- length(mats)
    q <- lapply(mats, function(m) {
      apply(m, 2, stats::quantile, probs = probs, names = FALSE)
    })
    # Array: quantile x marker x draw.
    arr <- array(
      unlist(q),
      dim = c(length(probs), length(markers), n_draws)
    )
    enough <- n_draws >= 2
    tibble::tibble(
      group = g,
      marker = rep(markers, each = length(probs)),
      prob = rep(probs, times = length(markers)),
      sd_across_draws = if (enough) {
        as.vector(apply(arr, c(1, 2), stats::sd))
      } else {
        NA_real_
      },
      range_across_draws = if (enough) {
        as.vector(apply(arr, c(1, 2), \(v) diff(range(v))))
      } else {
        NA_real_
      },
      n_draws = as.integer(n_draws),
      mean_events = mean(vapply(mats, nrow, integer(1))),
      status = if (enough) "ok" else "not_assessed"
    )
  })
  dplyr::bind_rows(rows)
}

# Fidelity to the source ----------------------------------------------------------

#' What a synthetic reference changed against its source
#'
#' For each synthetic reference, compares its events with the pooled events of
#' the files it was drawn from. Marker rows give the exact one-dimensional EMD
#' ([sw_emd_1d()]) in the units of the transformed data. Population rows give
#' the share in the source, the share in the reference and their ratio. With
#' `labels = NULL` there are no population rows.
#' Rank rows within check and marker only; nothing is rolled up per reference.
#'
#' @param synref An `sw_synref`.
#' @param data Named list of the numeric event matrices the reference was drawn
#'   from (the same transformed data).
#' @param labels Named list of population labels per file, or `NULL` for
#'   "no population rows".
#' @param n_events Events drawn from the source pool for the EMD
#'   (uncalibrated), or `NULL` to use the whole pool.
#' @param seed Random seed for that draw (uncalibrated). The caller's RNG state
#'   is restored.
#'
#' @return A tibble: `synref_id`, `check` (`emd_to_source` or `share_ratio`),
#'   `item` (marker or population), `emd_to_source`, `share_source`,
#'   `share_synref`, `share_ratio` and `status`.
#'
#' @examples
#' set.seed(1)
#' d <- list(a = cbind(M1 = rnorm(500), M2 = rnorm(500)))
#' man <- data.frame(sample_id = "a", batch = "b1")
#' ref <- sw_synref_downsample(
#'   d, man, groups = list(r = "a"), strata = "file", labels = NULL,
#'   n_events = 100, n_draws = 1, seed = 1
#' )
#' sw_synref_fidelity(ref, d, labels = NULL, n_events = NULL, seed = 1)
#'
#' @export
sw_synref_fidelity <- function(synref, data, labels, n_events, seed) {
  .sw_stop_missing(c(
    synref = missing(synref),
    data = missing(data),
    labels = missing(labels),
    n_events = missing(n_events),
    seed = missing(seed)
  ))
  .sw_validate_synref(synref)
  .sw_synref_check_data(data)
  .sw_synref_check_seed(seed)
  if (!is.null(n_events)) {
    .sw_synref_check_count(n_events, "n_events")
  }
  prov <- synref[["provenance"]]
  absent <- setdiff(unique(prov[["source_sample"]]), names(data))
  if (length(absent) > 0) {
    stop(
      "Source sample(s) not in 'data': ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.null(labels)) {
    ok <- all(vapply(
      unique(prov[["source_sample"]]),
      \(id) length(labels[[id]]) == nrow(data[[id]]),
      logical(1)
    ))
    if (!ok) {
      stop(
        "'labels' needs one vector per source file with one value per event.",
        call. = FALSE
      )
    }
  }
  markers <- colnames(synref[["events"]][[1]])
  if (!identical(markers, colnames(data[[1]]))) {
    stop("'data' and 'synref' need the same markers.", call. = FALSE)
  }
  rows <- .sw_with_seed(seed, {
    lapply(names(synref[["events"]]), function(id) {
      p <- prov[prov[["synref_id"]] == id, , drop = FALSE]
      files <- unique(p[["source_sample"]])
      pool <- do.call(rbind, unname(data[files]))
      pool_lab <- if (!is.null(labels)) {
        as.character(unlist(labels[files], use.names = FALSE))
      }
      keep <- seq_len(nrow(pool))
      if (!is.null(n_events) && n_events < nrow(pool)) {
        keep <- sort(sample.int(nrow(pool), n_events))
      }
      ev <- synref[["events"]][[id]]
      emd <- tibble::tibble(
        synref_id = id,
        check = "emd_to_source",
        item = markers,
        emd_to_source = vapply(
          markers,
          \(m) sw_emd_1d(ev[, m], pool[keep, m]),
          numeric(1)
        ),
        share_source = NA_real_,
        share_synref = NA_real_,
        share_ratio = NA_real_,
        status = "ok"
      )
      if (is.null(labels)) {
        return(emd)
      }
      # Population of each drawn event, looked up from the source labels.
      offs <- c(0, cumsum(vapply(data[files], nrow, integer(1))))[
        match(p[["source_sample"]], files)
      ]
      syn_lab <- pool_lab[offs + p[["source_row"]]]
      pops <- sort(unique(pool_lab))
      sh_src <- vapply(pops, \(q) mean(pool_lab == q), numeric(1))
      sh_syn <- vapply(pops, \(q) mean(syn_lab == q), numeric(1))
      pop <- tibble::tibble(
        synref_id = id,
        check = "share_ratio",
        item = pops,
        emd_to_source = NA_real_,
        share_source = unname(sh_src),
        share_synref = unname(sh_syn),
        share_ratio = unname(sh_syn / sh_src),
        status = "ok"
      )
      dplyr::bind_rows(emd, pop)
    })
  })
  dplyr::bind_rows(rows)
}
