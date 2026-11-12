# R/ref_noise_floor.R
# How different are two draws of the same reference by chance? A cell split of
# one file gives a sampling-only floor. Two aliquots of one reference in one batch
# give an aliquot floor. Own code; uses sw_synref_split() for the cell split.

# EMD of two samples. `n_quantiles` NULL gives the exact 1D EMD; a number gives
# the mean absolute gap of the two quantile functions on the grid
# (k - 0.5) / K, type-7 quantiles, as sw_corrector_compare() does.
.sw_ref_emd <- function(x, y, n_quantiles) {
  if (is.null(n_quantiles)) {
    return(sw_emd_1d(x, y))
  }
  probs <- (seq_len(n_quantiles) - 0.5) / n_quantiles
  mean(abs(
    stats::quantile(x, probs, type = 7, names = FALSE) -
      stats::quantile(y, probs, type = 7, names = FALSE)
  ))
}

#' Measure the chance difference between two draws of a reference
#'
#' @description
#' Give, per reference and marker, how far apart two draws of the same material
#' land by chance, so the between-batch EMD of [sw_ref_consistency()] can be
#' read against it. Two floors are reported, and each row says which one it is.
#'
#' * `floor_type = "sampling only"`: each reference file is cut into two
#'   disjoint halves of `n_events` events ([sw_synref_split()]) and the EMD of
#'   the halves is taken, `n_splits` times with different seeds. The halves have
#'   identical biology and the same run, so this floor holds event-level noise
#'   only. It is a lower bound. It does not hold aliquot-to-aliquot or
#'   run-to-run noise.
#' * `floor_type = "aliquot"`: added only where a reference has at least 2
#'   aliquots in one batch (manifest column `aliquot`). The events of each
#'   aliquot are pooled, `n_events` are drawn from each, and every aliquot pair
#'   is compared, `n_splits` times. This floor holds sampling noise and the
#'   difference between aliquots of the same material.
#'
#' @details
#' **Rows, not a summary.** Every split is its own row. This function picks no
#' summary statistic and no cut: read the spread of `value` across splits and
#' units yourself. Cells are not replicates, so the splits of one file are not
#' independent runs; more events do not make a floor more independent.
#'
#' **Thin rows.** A file with fewer than `2 * n_events` events cannot be split,
#' and an aliquot with fewer than `n_events` pooled events cannot be drawn. Such
#' units get `status = "not_assessed"`, `value = NA` and a `reason`.
#'
#' **Scale is yours.** `data` must already be transformed; nothing is
#' transformed here.
#'
#' @param data A named list of event matrices or data frames, one per file.
#' @param manifest A manifest with `sample_id`, `batch`, `reference_role` and
#'   `reference_id`, and optionally `aliquot` (see [sw_sample_manifest()]).
#' @param markers Character vector of marker (column) names.
#' @param n_events Whole number, events per half or per aliquot draw
#'   (uncalibrated).
#' @param n_splits Whole number, repeated splits per unit (uncalibrated; the
#'   plan proposes 20).
#' @param n_quantiles `NULL` for the exact 1D EMD ([sw_emd_1d()]), or a whole
#'   number of quantile points K. With K, the EMD is the mean absolute gap of the
#'   two quantile functions at `(k - 0.5) / K`, the form `sw_corrector_compare()`
#'   uses (uncalibrated).
#' @param seed Whole number. Split seeds are drawn from it; the caller's RNG
#'   state is restored (uncalibrated).
#' @param roles Character vector of `reference_role` values that count as
#'   reference material (uncalibrated).
#'
#' @return A tibble of class `sw_ref_noise_floor` with columns `reference_id`,
#'   `batch`, `floor_type`, `unit` (the file for a sampling floor, `"a1 ~ a2"`
#'   for an aliquot pair), `marker`, `split`, `value`, `status` and `reason`.
#'   Attribute `settings` echoes the arguments.
#'
#' @seealso [sw_ref_consistency()], [sw_synref_split()], [sw_emd_1d()]
#'
#' @examples
#' set.seed(1)
#' d <- list(
#'   r1 = cbind(A = stats::rnorm(400)),
#'   r2 = cbind(A = stats::rnorm(400))
#' )
#' man <- data.frame(
#'   sample_id = c("r1", "r2"),
#'   batch = c("b1", "b2"),
#'   reference_role = "technical_replicate",
#'   reference_id = "ref1"
#' )
#' sw_ref_noise_floor(
#'   d,
#'   man,
#'   markers = "A",
#'   n_events = 100, # illustrative, uncalibrated
#'   n_splits = 3, # illustrative
#'   n_quantiles = NULL,
#'   seed = 1,
#'   roles = "technical_replicate"
#' )
#'
#' @export
sw_ref_noise_floor <- function(
  data,
  manifest,
  markers,
  n_events,
  n_splits,
  n_quantiles,
  seed,
  roles
) {
  .sw_stop_missing(c(
    data = missing(data),
    manifest = missing(manifest),
    markers = missing(markers),
    n_events = missing(n_events),
    n_splits = missing(n_splits),
    n_quantiles = missing(n_quantiles),
    seed = missing(seed),
    roles = missing(roles)
  ))
  .sw_check_count(n_events, "n_events")
  .sw_check_count(n_splits, "n_splits")
  if (!is.null(n_quantiles)) {
    .sw_check_count(n_quantiles, "n_quantiles")
  }
  .sw_ref_check_seed(seed)
  picked <- .sw_ref_select(data, manifest, roles)
  sel <- picked[["sel"]]
  mats <- .sw_ref_matrices(picked[["data"]], markers)
  if (anyNA(do.call(rbind, mats))) {
    stop("Events must not contain NA in the chosen markers.", call. = FALSE)
  }
  split_seeds <- .sw_with_seed(seed, sample.int(1e9, n_splits))
  man <- as.data.frame(sel[c("sample_id", "batch")])

  # Sampling-only floor: cut each file into two disjoint halves.
  sampling <- lapply(seq_len(nrow(sel)), \(i) {
    f <- sel[["sample_id"]][i]
    if (nrow(mats[[f]]) < 2 * n_events) {
      return(.sw_ref_floor_rows(
        sel[i, ],
        "sampling only",
        f,
        markers,
        NA_integer_,
        NA_real_,
        sprintf("file has %d events, fewer than 2 * n_events", nrow(mats[[f]]))
      ))
    }
    lapply(seq_len(n_splits), \(k) {
      parts <- sw_synref_split(
        mats[f],
        man,
        source = f,
        n_parts = 2,
        events_per_part = n_events,
        seed = split_seeds[k]
      )[["events"]]
      v <- vapply(
        markers,
        \(mk) .sw_ref_emd(parts[[1]][, mk], parts[[2]][, mk], n_quantiles),
        numeric(1)
      )
      .sw_ref_floor_rows(
        sel[i, ],
        "sampling only",
        f,
        markers,
        k,
        v,
        NA_character_
      )
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_rows()

  # Aliquot floor: pairs of aliquots of one reference inside one batch.
  aliquot_rows <- NULL
  if ("aliquot" %in% names(sel)) {
    groups <- sel |>
      dplyr::filter(!is.na(aliquot)) |>
      dplyr::distinct(reference_id, batch, aliquot) |>
      dplyr::summarise(n_aliquots = dplyr::n(), .by = c(reference_id, batch)) |>
      dplyr::filter(n_aliquots >= 2)
    aliquot_rows <- lapply(seq_len(nrow(groups)), \(g) {
      rid <- groups[["reference_id"]][g]
      bt <- groups[["batch"]][g]
      rows <- sel[sel[["reference_id"]] == rid & sel[["batch"]] == bt, ]
      pooled <- split(rows[["sample_id"]], rows[["aliquot"]]) |>
        lapply(\(ids) do.call(rbind, mats[ids]))
      pooled <- pooled[sort(names(pooled))]
      info <- tibble::tibble(reference_id = rid, batch = bt)
      lapply(utils::combn(names(pooled), 2, simplify = FALSE), \(pr) {
        unit <- paste(pr, collapse = " ~ ")
        short <- pr[vapply(pr, \(a) nrow(pooled[[a]]) < n_events, logical(1))]
        if (length(short) > 0) {
          return(.sw_ref_floor_rows(
            info,
            "aliquot",
            unit,
            markers,
            NA_integer_,
            NA_real_,
            paste0(
              "aliquot ",
              short[1],
              " has fewer than n_events pooled events"
            )
          ))
        }
        lapply(seq_len(n_splits), \(k) {
          draw <- .sw_with_seed(
            split_seeds[k],
            lapply(pooled[pr], \(m) {
              m[sample.int(nrow(m), n_events), , drop = FALSE]
            })
          )
          v <- vapply(
            markers,
            \(mk) .sw_ref_emd(draw[[1]][, mk], draw[[2]][, mk], n_quantiles),
            numeric(1)
          )
          .sw_ref_floor_rows(
            info,
            "aliquot",
            unit,
            markers,
            k,
            v,
            NA_character_
          )
        }) |>
          dplyr::bind_rows()
      }) |>
        dplyr::bind_rows()
    }) |>
      dplyr::bind_rows()
  }

  out <- dplyr::bind_rows(sampling, aliquot_rows) |>
    dplyr::arrange(reference_id, floor_type, batch, unit, marker, split)
  structure(
    out,
    settings = list(
      markers = markers,
      n_events = n_events,
      n_splits = n_splits,
      n_quantiles = n_quantiles,
      seed = seed,
      roles = roles
    ),
    class = c("sw_ref_noise_floor", class(out))
  )
}

# Rows for one unit and one split: one row per marker. `value` is a named
# vector in `markers` order, or a single NA with a reason.
.sw_ref_floor_rows <- function(
  info,
  floor_type,
  unit,
  markers,
  split,
  value,
  reason
) {
  tibble::tibble(
    reference_id = info[["reference_id"]],
    batch = info[["batch"]],
    floor_type = floor_type,
    unit = unit,
    marker = markers,
    split = as.integer(split),
    value = unname(value),
    status = if (anyNA(value)) "not_assessed" else "assessed",
    reason = reason
  )
}
