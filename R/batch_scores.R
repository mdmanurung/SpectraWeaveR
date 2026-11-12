# R/batch_scores.R
# Marker-level batch-effect scores from files with a batch label: pairwise EMD,
# spread of batch medians, variance explained by batch, and acquisition-order
# drift within each batch. Written from scratch from standard definitions.

#' Score batch effects per marker
#'
#' @description
#' Give every marker four plain views of a batch effect, from files that carry
#' a batch label. The result is numbers only; no batch or marker is marked bad.
#'
#' * `pair_emd`: exact 1D earth mover's distance ([sw_emd_1d()]) between the
#'   pooled events of two batches, for each batch pair. It is in the units of
#'   the marker scale.
#' * `rest_emd`: per batch, exact 1D EMD between that batch's pooled events and
#'   the pooled events of all other batches. `rest_emd_z` is its robust z over
#'   batches within the marker (`.sw_robust_z()`), so one odd batch stands out
#'   even when every pair involving it looks alike. It needs at least 3
#'   batches (fewer give `NA`) and is `NA` when the batch spread is 0. With
#'   `reference_group`, both are computed within each group across batches,
#'   then combined by the median over groups like every other row.
#' * `median_range`: highest minus lowest batch median of the marker. A batch
#'   median is the median over its files of each file's median.
#' * `eta2`: share of the variance of file medians that the batch label
#'   explains (between-batch sum of squares over total). Near 0 means files
#'   from different batches look alike; near 1 means the batch decides the
#'   value.
#' * `spread_ratio`: per batch, the median over its files of each file's SD,
#'   divided by the same median over all files. 1 means the same width as the
#'   cohort; a batch with a stretched marker or a split into two modes is above
#'   1. Medians and pooled EMD can miss a symmetric widening.
#' * `order_rho`: Spearman correlation between a file's marker median and its
#'   acquisition order, within each batch. A steady drift inside a batch shows
#'   here. Only when `order` is given and a batch has at least 4 files.
#'
#' @details
#' **Batch and biology.** A batch label that is confounded with a condition
#' (all cases in one batch) cannot be separated from biology by any of these
#' scores. Check the design first with [sw_design_audit()].
#'
#' **Scale is yours.** `data` must be on the analysis scale. This function
#' never transforms.
#'
#' **No defaults.** `n_events`, `seed` and `order` change the numbers, so each
#' must be passed. Events are drawn per file without replacement with `seed`;
#' a file with fewer events is an error. The caller's RNG state is restored.
#'
#' @param data A named list. Each element is a numeric matrix or data frame
#'   with one column per marker. Names are the file ids.
#' @param batch A named vector of batch labels, one per file in `data`.
#' @param markers Character vector of marker (column) names to score.
#' @param n_events Whole number, events drawn per file.
#' @param seed Whole number for the draws.
#' @param order `NULL`, or a named numeric vector with the acquisition order
#'   (or time) of each file. Used only for `order_rho`.
#' @param reference_group `NULL`, or a named character vector with one group
#'   label per file (names = file ids, covering every file in `data`), for
#'   example the donor. With groups, every score is computed within each group
#'   from that group's files only, across batches, and the groups are combined
#'   per row (marker and batch, or batch pair) by the median over groups. This
#'   removes donor differences from the batch scores. A group with files from
#'   fewer than two batches is skipped. `n_files` is the sum over groups.
#'   `NULL` scores all files together. A data frame (two-way form) is an
#'   error here, because removing the batch effect would remove what is measured.
#'
#' @param population `NULL` (skip), or a list with two elements. `labels` is a
#'   named list with one vector per file in `data`: the population label of
#'   every event, in the row order of that file (gate or cluster label; `NA`
#'   events are left out). `label_source` is `"gate"` or `"cluster"`. With
#'   `population`, the same scores are also computed inside each population,
#'   and two more checks are added (see Population rows). `NULL` returns the
#'   scores exactly as before.
#' @param min_events `NULL` when `population` is `NULL`; otherwise a whole
#'   number, the fewest events a file needs in a population to be scored there
#'   (uncalibrated: there is no data-derived value). A thinner file x population
#'   cell gives a row with `status = "not_assessed"` and `NA` value. It cannot
#'   exceed `n_events`.
#' @param group_summary `NULL` when `population` is `NULL`; otherwise `"median"`
#'   or `"mean"` (uncalibrated: the choice moves the numbers). It sets how
#'   per-group population rows are combined across the groups of
#'   `reference_group`. Rows without `population` always use the median. It is
#'   recorded but unused when `reference_group` is `NULL`.
#'
#' @section Population rows:
#' With `population`, up to `n_events` events are drawn per file and population
#' (all of them when the file has fewer, but at least `min_events`). Every check
#' above is computed inside each population from the files that pass
#' `min_events`; a population with fewer than two batches among them gets one
#' row with `status = "not_assessed"`. Two checks are added:
#'
#' * `pop_share` and `share_rest_z`: per batch, the median over its files of the
#'   share of the file's labelled events in the population, and its robust z over
#'   batches (needs 3 batches, else `NA`). Files, not cells, are the replicates.
#' * `batch_order_rho`: only when `order` is given. Spearman correlation across
#'   batches between the batch median of the marker and the batch order (the
#'   median `order` of its files). It needs 3 batches, else `NA`. It is also
#'   computed on all events (population `NA`).
#'
#' Cluster labels depend on the batch they were fitted in, so with
#' `label_source = "cluster"` every population row of the result is capped at
#' `eligible_for_review` (column `eligibility_cap`).
#'
#' @return A tibble of class `sw_batch_scores` with columns `check`, `marker`,
#'   `batch` (the batch, or `"b1 ~ b2"` for a pair, or `NA` for the all-batch
#'   scores `median_range` and `eta2`), `value`, and `n_files` (files behind
#'   the value; for a pair, the smaller of the two batches). With `population`
#'   it also has `population` (`NA` for the all-event rows), `file` (set only
#'   on thin-cell rows), `n_events` (events drawn behind a population row),
#'   `status` (`"assessed"` or `"not_assessed"`), `label_source` and
#'   `eligibility_cap`.
#'
#' @seealso [sw_design_audit()], [sw_anomaly_features()]
#'
#' @examples
#' set.seed(1)
#' files <- lapply(1:6, function(i) {
#'   matrix(stats::rnorm(600) + (i > 3), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' })
#' names(files) <- paste0("f", 1:6)
#' sw_batch_scores(
#'   files,
#'   batch = stats::setNames(rep(c("b1", "b2"), each = 3), names(files)),
#'   markers = c("A", "B"),
#'   n_events = 100, # illustrative
#'   seed = 1,
#'   order = NULL,
#'   reference_group = NULL,
#'   population = NULL
#' )
#'
#' @export
sw_batch_scores <- function(
  data,
  batch,
  markers,
  n_events,
  seed,
  order,
  reference_group,
  population,
  min_events,
  group_summary
) {
  .sw_stop_missing(c(
    data = missing(data),
    batch = missing(batch),
    markers = missing(markers),
    n_events = missing(n_events),
    seed = missing(seed),
    order = missing(order),
    reference_group = missing(reference_group),
    population = missing(population)
  ))
  if (!is.null(population)) {
    .sw_stop_missing(c(
      min_events = missing(min_events),
      group_summary = missing(group_summary)
    ))
    .sw_check_count(min_events, "min_events")
    if (min_events > n_events) {
      stop("`min_events` cannot exceed `n_events`.", call. = FALSE)
    }
    if (
      !is.character(group_summary) ||
        length(group_summary) != 1 ||
        !group_summary %in% c("median", "mean")
    ) {
      stop("`group_summary` must be \"median\" or \"mean\".", call. = FALSE)
    }
    .sw_check_population(population, names(data), data)
  }
  if (!is.list(data) || is.null(names(data)) || anyDuplicated(names(data))) {
    stop("`data` must be a named list of files.", call. = FALSE)
  }
  if (!is.character(markers) || length(markers) == 0) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  .sw_check_count(n_events, "n_events")
  if (
    !is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed != round(seed)
  ) {
    stop("`seed` must be one whole number.", call. = FALSE)
  }
  if (is.null(names(batch)) || !all(names(data) %in% names(batch))) {
    stop(
      "`batch` must be a named vector covering every file in `data`.",
      call. = FALSE
    )
  }
  if (
    !is.null(order) &&
      (is.null(names(order)) || !all(names(data) %in% names(order)))
  ) {
    stop(
      "`order` must be NULL or a named vector covering every file.",
      call. = FALSE
    )
  }
  files <- sort(names(data))
  batch <- as.character(batch[files])
  if (anyNA(batch)) {
    stop("`batch` has missing labels.", call. = FALSE)
  }
  grp <- .sw_check_reference_group(reference_group, files)
  if (inherits(grp, "sw_two_way")) {
    stop(
      "A two-way (data-frame) `reference_group` would remove the batch effect ",
      "that sw_batch_scores() measures. Pass a named vector (for example the ",
      "donor) or NULL.",
      call. = FALSE
    )
  }
  batches <- sort(unique(batch))
  if (length(batches) < 2) {
    stop("Need at least two batches.", call. = FALSE)
  }

  mats <- lapply(files, function(f) {
    m <- as.matrix(data[[f]])
    if (!all(markers %in% colnames(m))) {
      stop(
        "File '",
        f,
        "' lacks marker(s): ",
        paste(setdiff(markers, colnames(m)), collapse = ", "),
        call. = FALSE
      )
    }
    if (nrow(m) < n_events) {
      stop(
        "File '",
        f,
        "' has ",
        nrow(m),
        " events, fewer than `n_events` (",
        n_events,
        ").",
        call. = FALSE
      )
    }
    m[, markers, drop = FALSE]
  })
  names(mats) <- files
  drawn <- .sw_with_seed(
    seed,
    lapply(mats, function(m) {
      m[sample.int(nrow(m), n_events), , drop = FALSE]
    })
  )
  if (anyNA(do.call(rbind, drawn))) {
    stop("Events must not contain NA in the chosen markers.", call. = FALSE)
  }

  extra <- !is.null(population) && !is.null(order)
  # One set of scores over all files, or one per reference group (only that
  # group's files, across batches) combined by the median over groups.
  if (is.null(grp)) {
    out <- .sw_batch_scores_core(drawn, files, batch, markers, order, extra)
  } else {
    per_group <- lapply(unique(grp), function(g) {
      keep <- files[grp == g]
      if (length(unique(batch[files %in% keep])) < 2) {
        return(NULL)
      }
      .sw_batch_scores_core(
        drawn[keep],
        keep,
        batch[files %in% keep],
        markers,
        order,
        extra
      )
    })
    per_group <- dplyr::bind_rows(per_group)
    if (nrow(per_group) == 0) {
      stop(
        "No reference group has files from at least two batches.",
        call. = FALSE
      )
    }
    out <- per_group |>
      dplyr::summarise(
        value = stats::median(value),
        n_files = sum(n_files),
        .by = c(check, marker, batch)
      )
    # Keep the row order of the all-files result: check, marker order as first seen.
    out <- out[
      order(
        match(out$check, unique(per_group$check)),
        match(out$marker, markers)
      ),
    ]
    out <- tibble::as_tibble(out)
  }
  if (!is.null(population)) {
    out <- .sw_batch_scores_population(
      out,
      population,
      mats,
      files,
      batch,
      markers,
      order,
      grp,
      n_events,
      min_events,
      group_summary,
      seed
    )
  }
  class(out) <- c("sw_batch_scores", class(out))
  out
}

# Scores for one set of files. `drawn` is a named list of event matrices in the
# order of `files`; `batch` has one label per file, in that order.
.sw_batch_scores_core <- function(
  drawn,
  files,
  batch,
  markers,
  order,
  extra = FALSE
) {
  batches <- sort(unique(batch))
  # File medians per marker (rows = files).
  fmed <- t(vapply(
    drawn,
    function(m) apply(m, 2, stats::median),
    numeric(length(markers))
  ))
  if (length(markers) == 1) {
    fmed <- matrix(fmed, ncol = 1)
  }
  colnames(fmed) <- markers
  rownames(fmed) <- files
  # File SDs per marker, for the spread ratio.
  fsd <- t(vapply(
    drawn,
    function(m) apply(m, 2, stats::sd),
    numeric(length(markers))
  ))
  if (length(markers) == 1) {
    fsd <- matrix(fsd, ncol = 1)
  }
  colnames(fsd) <- markers
  rownames(fsd) <- files
  n_by_batch <- table(batch)[batches]

  rows <- list()
  add <- function(check, marker, b, value, n) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(
      check = check,
      marker = marker,
      batch = b,
      value = value,
      n_files = as.integer(n)
    )
  }

  # Pooled events per batch, for pairwise EMD.
  pooled <- lapply(batches, function(b) {
    do.call(rbind, drawn[files[batch == b]])
  })
  names(pooled) <- batches
  pairs <- utils::combn(batches, 2, simplify = FALSE)
  for (p in pairs) {
    plan <- .sw_emd_plan(nrow(pooled[[p[1]]]), nrow(pooled[[p[2]]]))
    for (mk in markers) {
      v <- .sw_emd_sorted(
        sort(pooled[[p[1]]][, mk]),
        sort(pooled[[p[2]]][, mk]),
        plan
      )
      add("pair_emd", mk, paste(p, collapse = " ~ "), v, min(n_by_batch[p]))
    }
  }

  # One batch against all other batches pooled, then robust z over batches.
  for (mk in markers) {
    rest <- vapply(
      batches,
      function(b) {
        mine <- sort(pooled[[b]][, mk])
        others <- sort(unlist(lapply(pooled[setdiff(batches, b)], \(m) {
          m[, mk]
        })))
        .sw_emd_sorted(mine, others, .sw_emd_plan(length(mine), length(others)))
      },
      numeric(1)
    )
    rest_z <- .sw_robust_z(rest)
    for (i in seq_along(batches)) {
      add("rest_emd", mk, batches[i], rest[[i]], n_by_batch[[batches[i]]])
      add("rest_emd_z", mk, batches[i], rest_z[[i]], n_by_batch[[batches[i]]])
    }
  }

  for (mk in markers) {
    x <- fmed[, mk]
    bmed <- tapply(x, batch, stats::median)
    add("median_range", mk, NA_character_, max(bmed) - min(bmed), length(files))
    sst <- sum((x - mean(x))^2)
    bmean <- as.numeric(tapply(x, batch, mean)[batches])
    ssb <- sum(as.numeric(n_by_batch) * (bmean - mean(x))^2)
    add(
      "eta2",
      mk,
      NA_character_,
      if (sst > 0) ssb / sst else NA_real_,
      length(files)
    )
    # Spread: batch median of file SDs over the median of all file SDs. Catches a
    # change of width or a split into modes that medians and pooled EMD miss.
    all_sd <- stats::median(fsd[, mk])
    for (b in batches) {
      add(
        "spread_ratio",
        mk,
        b,
        if (all_sd > 0) {
          stats::median(fsd[batch == b, mk]) / all_sd
        } else {
          NA_real_
        },
        n_by_batch[[b]]
      )
    }
    if (!is.null(order)) {
      for (b in batches) {
        idx <- batch == b
        if (sum(idx) >= 4) {
          rho <- suppressWarnings(stats::cor(
            x[idx],
            order[files][idx],
            method = "spearman"
          ))
          add("order_rho", mk, b, rho, sum(idx))
        }
      }
      if (extra) {
        # Batch median of the marker against batch order, across batches.
        b_order <- tapply(order[files], batch, stats::median)[batches]
        rho_b <- if (length(batches) >= 3) {
          suppressWarnings(stats::cor(
            as.numeric(bmed[batches]),
            as.numeric(b_order),
            method = "spearman"
          ))
        } else {
          NA_real_
        }
        add("batch_order_rho", mk, NA_character_, rho_b, length(files))
      }
    }
  }

  dplyr::bind_rows(rows)
}
