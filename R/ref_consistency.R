# R/ref_consistency.R
# Is a repeated reference the same material across batches? Per reference id,
# marker and batch: pair_emd, rest_emd, rest_emd_z, median_shift_z, order_rho and
# an event-count ratio. Built on the per-group rows of the batch scores. Own code.

#' Check that a reference looks the same in every batch
#'
#' @description
#' For each reference id (the same material measured in several batches), give
#' per marker and batch the numbers that show whether the material moved. The
#' result is a long table of numbers with a `status`. No reference, batch or
#' marker is marked consistent or inconsistent. Rank within a check and a
#' marker, never over markers.
#'
#' * `pair_emd`: exact 1D earth mover's distance ([sw_emd_1d()]) between the
#'   pooled events of the reference in two batches (`batch` is `"b1 ~ b2"`).
#' * `rest_emd`: the reference in one batch against the same reference in all
#'   other batches pooled. `rest_emd_z` is its robust z over the batches of
#'   that reference and marker.
#' * `median_shift_z`: robust z over batches of the batch median (median over
#'   files of each file's marker median). It shows a plain shift of the centre.
#' * `order_rho`: Spearman correlation of file marker median with `order`, in
#'   a batch with at least 4 files of the reference. Only when `order` is given.
#' * `count_ratio`: median events per file in the batch over the median events
#'   per file of the whole reference (`marker` is `NA`). A batch far from 1 has a
#'   different amount of material.
#'
#' @details
#' **Cells are not replicates.** The unit that varies between batches is the
#' reference sample (aliquot), not the cell. With few batches the z-scores are
#' weak; `min_batches_per_reference` sets how many batches a reference needs
#' before `rest_emd_z` and `median_shift_z` are filled in. `pair_emd` and
#' `rest_emd` need only 2 batches.
#'
#' **Thin rows.** A file with fewer than `min_events` events is left out of
#' every score. A reference with fewer than 2 batches left, or a value that
#' cannot be computed (a robust z needs 3 batches and a non-zero spread), gets
#' `status = "not_assessed"`, `value = NA` and a `reason`. Nothing is dropped
#' silently.
#'
#' **Scale is yours (DR37).** `data` must already be transformed. This function
#' never transforms. `transform` and `cofactor` are recorded in attribute
#' `settings` so a reader can see which scale produced the numbers; they are not
#' checked against the data.
#'
#' **No defaults.** Every argument changes the numbers and must be passed.
#' Events are drawn per file without replacement with `seed` (a file with fewer
#' than `n_events` events but at least `min_events` uses all its events). The
#' caller's RNG state is restored.
#'
#' @param data A named list of event matrices or data frames, one per file.
#' @param manifest A manifest with `sample_id`, `batch`, `reference_role` and
#'   `reference_id` (see [sw_sample_manifest()]).
#' @param markers Character vector of marker (column) names to score.
#' @param n_events Whole number, events drawn per file (uncalibrated).
#' @param seed Whole number for the draws (uncalibrated).
#' @param order `NULL`, or a named numeric vector with the acquisition order of
#'   each reference file (uncalibrated; explicit skip when `NULL`).
#' @param min_batches_per_reference Whole number, batches a reference needs for
#'   the z-based rows (uncalibrated; the plan proposes 3, which is reasoning,
#'   not a sourced number).
#' @param min_events Whole number, events a file needs to be scored
#'   (uncalibrated).
#' @param roles Character vector of `reference_role` values that count as
#'   reference material (uncalibrated).
#' @param transform One string naming the transform the caller applied, for
#'   example `"asinh"`. Recorded only.
#' @param cofactor One positive number (the cofactor of `transform`), or `NULL`
#'   when the transform has none. Recorded only.
#'
#' @return A tibble of class `sw_ref_consistency` with columns `reference_id`,
#'   `check`, `marker`, `batch`, `value`, `n_files`, `status`
#'   (`"assessed"` or `"not_assessed"`) and `reason`. Attribute `settings`
#'   echoes every argument except `data` and `manifest`.
#'
#' @seealso [sw_ref_inventory()], [sw_ref_noise_floor()], [sw_batch_scores()],
#'   [sw_flag_iqr()]
#'
#' @examples
#' set.seed(1)
#' ids <- paste0("r", 1:3)
#' d <- lapply(ids, function(i) cbind(A = stats::rnorm(200), B = stats::rnorm(200)))
#' names(d) <- ids
#' man <- data.frame(
#'   sample_id = ids,
#'   batch = c("b1", "b2", "b3"),
#'   reference_role = "technical_replicate",
#'   reference_id = "ref1"
#' )
#' sw_ref_consistency(
#'   d,
#'   man,
#'   markers = c("A", "B"),
#'   n_events = 100, # illustrative, uncalibrated
#'   seed = 1,
#'   order = NULL,
#'   min_batches_per_reference = 3, # illustrative
#'   min_events = 50, # illustrative
#'   roles = "technical_replicate",
#'   transform = "asinh",
#'   cofactor = 6000
#' )
#'
#' @export
sw_ref_consistency <- function(
  data,
  manifest,
  markers,
  n_events,
  seed,
  order,
  min_batches_per_reference,
  min_events,
  roles,
  transform,
  cofactor
) {
  .sw_stop_missing(c(
    data = missing(data),
    manifest = missing(manifest),
    markers = missing(markers),
    n_events = missing(n_events),
    seed = missing(seed),
    order = missing(order),
    min_batches_per_reference = missing(min_batches_per_reference),
    min_events = missing(min_events),
    roles = missing(roles),
    transform = missing(transform),
    cofactor = missing(cofactor)
  ))
  .sw_check_count(n_events, "n_events")
  .sw_check_count(min_events, "min_events")
  .sw_check_count(min_batches_per_reference, "min_batches_per_reference")
  .sw_ref_check_seed(seed)
  if (!is.character(transform) || length(transform) != 1 || is.na(transform)) {
    stop(
      "`transform` must be one string, for example \"asinh\".",
      call. = FALSE
    )
  }
  if (
    !is.null(cofactor) &&
      (!is.numeric(cofactor) ||
        length(cofactor) != 1 ||
        is.na(cofactor) ||
        cofactor <= 0)
  ) {
    stop("`cofactor` must be NULL or one positive number.", call. = FALSE)
  }
  picked <- .sw_ref_select(data, manifest, roles)
  sel <- picked[["sel"]]
  mats <- .sw_ref_matrices(picked[["data"]], markers)
  if (
    !is.null(order) &&
      (!is.numeric(order) ||
        is.null(names(order)) ||
        !all(sel[["sample_id"]] %in% names(order)))
  ) {
    stop(
      "`order` must be NULL or a named numeric vector covering every ",
      "reference file.",
      call. = FALSE
    )
  }
  sel[["n_all"]] <- vapply(mats, nrow, integer(1))[sel[["sample_id"]]]
  used <- sel[sel[["n_all"]] >= min_events, ]

  # Scores per reference id, only from files that meet min_events.
  scores <- tibble::tibble()
  if (nrow(used) > 0) {
    files <- used[["sample_id"]]
    drawn <- .sw_with_seed(
      seed,
      lapply(mats[files], \(m) {
        m[sample.int(nrow(m), min(n_events, nrow(m))), , drop = FALSE]
      })
    )
    if (anyNA(do.call(rbind, drawn))) {
      stop("Events must not contain NA in the chosen markers.", call. = FALSE)
    }
    scores <- .sw_batch_scores_groups(
      drawn,
      files,
      stats::setNames(used[["batch"]], files)[files],
      markers,
      order,
      stats::setNames(used[["reference_id"]], files)[files]
    )
    if (nrow(scores) > 0) {
      scores <- scores |>
        dplyr::filter(
          check %in% c("pair_emd", "rest_emd", "rest_emd_z", "order_rho")
        ) |>
        dplyr::rename(reference_id = group)
    }
    shift <- .sw_ref_median_shift(drawn, used, markers)
    scores <- dplyr::bind_rows(scores, shift)
  }
  scores <- dplyr::bind_rows(
    tibble::tibble(
      reference_id = character(),
      check = character(),
      marker = character(),
      batch = character(),
      value = numeric(),
      n_files = integer()
    ),
    scores
  )

  # A z needs enough batches; a reference with fewer than min_batches_per_reference
  # usable batches gets NA there.
  n_used <- used |>
    dplyr::summarise(n_used = dplyr::n_distinct(batch), .by = reference_id)
  scores <- scores |>
    dplyr::left_join(
      n_used,
      by = "reference_id",
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      is_z = check %in% c("rest_emd_z", "median_shift_z"),
      reason = dplyr::case_when(
        is_z & n_used < min_batches_per_reference ~
          "fewer usable batches than min_batches_per_reference",
        is.na(value) & is_z ~ "robust z needs 3 batches and a non-zero spread",
        is.na(value) &
          check == "order_rho" ~ "correlation undefined (constant values)",
        is.na(value) ~ "not computable",
        .default = NA_character_
      ),
      value = dplyr::if_else(
        is_z & n_used < min_batches_per_reference,
        NA_real_,
        value
      )
    ) |>
    dplyr::select(reference_id, check, marker, batch, value, n_files, reason)

  # Per-batch checks with no row (a batch with no usable file, or a reference with
  # fewer than 2 usable batches) become explicit not_assessed rows.
  have <- scores |>
    dplyr::filter(check %in% c("rest_emd", "rest_emd_z", "median_shift_z")) |>
    dplyr::select(reference_id, check, marker, batch)
  batch_has_file <- used[c("reference_id", "batch")] |>
    dplyr::distinct() |>
    dplyr::mutate(has_file = TRUE)
  absent_rows <- .sw_ref_grid(sel, markers) |>
    dplyr::anti_join(
      have,
      by = c("reference_id", "check", "marker", "batch")
    ) |>
    dplyr::left_join(
      batch_has_file,
      by = c("reference_id", "batch"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      value = NA_real_,
      n_files = 0L,
      reason = dplyr::if_else(
        is.na(has_file),
        "no file of this reference in this batch has min_events events",
        "reference has fewer than 2 batches with min_events events"
      )
    ) |>
    dplyr::select(-has_file)

  out <- dplyr::bind_rows(scores, absent_rows, .sw_ref_count_ratio(sel)) |>
    dplyr::mutate(
      status = dplyr::if_else(is.na(value), "not_assessed", "assessed")
    ) |>
    dplyr::select(
      reference_id,
      check,
      marker,
      batch,
      value,
      n_files,
      status,
      reason
    ) |>
    dplyr::arrange(reference_id, check, marker, batch)
  structure(
    out,
    settings = list(
      markers = markers,
      n_events = n_events,
      seed = seed,
      order = !is.null(order),
      min_batches_per_reference = min_batches_per_reference,
      min_events = min_events,
      roles = roles,
      transform = transform,
      cofactor = cofactor
    ),
    class = c("sw_ref_consistency", class(out))
  )
}

# Full grid of per-batch checks: every reference x batch (from the manifest) x
# marker x check, so a batch that could not be scored shows up as a row.
.sw_ref_grid <- function(sel, markers) {
  rb <- unique(sel[c("reference_id", "batch")])
  chk <- c("rest_emd", "rest_emd_z", "median_shift_z")
  rb |>
    dplyr::cross_join(tibble::tibble(marker = markers)) |>
    dplyr::cross_join(tibble::tibble(check = chk))
}

# Robust z over batches of the batch median (median over files of file medians),
# within each reference id and marker.
.sw_ref_median_shift <- function(drawn, used, markers) {
  fmed <- t(vapply(
    drawn,
    \(m) apply(m, 2, stats::median),
    numeric(length(markers))
  ))
  if (length(markers) == 1) {
    fmed <- matrix(fmed, ncol = 1)
  }
  colnames(fmed) <- markers
  rows <- lapply(unique(used[["reference_id"]]), \(r) {
    u <- used[used[["reference_id"]] == r, ]
    if (dplyr::n_distinct(u[["batch"]]) < 2) {
      return(NULL)
    }
    lapply(markers, \(mk) {
      bmed <- tapply(fmed[u[["sample_id"]], mk], u[["batch"]], stats::median)
      nf <- as.integer(table(u[["batch"]])[names(bmed)])
      tibble::tibble(
        reference_id = r,
        check = "median_shift_z",
        marker = mk,
        batch = names(bmed),
        value = .sw_robust_z(as.numeric(bmed)),
        n_files = nf
      )
    }) |>
      dplyr::bind_rows()
  })
  dplyr::bind_rows(rows)
}

# Events per file in a batch over events per file in the whole reference, from the
# raw counts of every reference file (thin files included). NA with one batch.
.sw_ref_count_ratio <- function(sel) {
  sel |>
    dplyr::mutate(ref_med = stats::median(n_all), .by = reference_id) |>
    dplyr::summarise(
      value = stats::median(n_all) / dplyr::first(ref_med),
      n_files = dplyr::n(),
      .by = c(reference_id, batch)
    ) |>
    dplyr::mutate(
      check = "count_ratio",
      marker = NA_character_,
      reason = if (dplyr::n() < 2) {
        rep("reference has fewer than 2 batches", dplyr::n())
      } else {
        rep(NA_character_, dplyr::n())
      },
      value = if (dplyr::n() < 2) rep(NA_real_, dplyr::n()) else value,
      .by = reference_id
    )
}
