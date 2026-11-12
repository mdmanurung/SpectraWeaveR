# R/ref_vs_study.R
# Does a repeated reference move between batches by less than donors differ
# inside one batch? Own design from the plan DR17 and DR18; no literature source.
# Distances are exact 1D earth mover's distances (sw_emd_1d) per marker.

#' Reference movement between batches next to donor spread inside a batch
#'
#' @description
#' For each reference id measured in two or more batches, compute the 1D EMD of
#' each marker between every pair of batches (`ref_emd`). Put it next to
#' `donor_emd_median`: the median pairwise EMD between donors of the same
#' batch. The ratio `ref_emd / donor_emd_median` says how large the movement
#' of the reference is against ordinary donor differences. A ratio well below 1
#' means the reference moves less than donors differ.
#'
#' @details
#' **Scale is yours.** `data` must already be transformed (DR37). This function
#' never transforms. `transform` and `cofactor` are only recorded in the
#' result so a reader can see the scale.
#'
#' **Donor spread uses same-batch donors only** (donor is nested in batch).
#' All files of one donor in one batch are pooled into one distribution. A
#' batch with fewer than `min_donors_per_batch` donors gives `NA` for
#' `donor_emd_median` and the row is `not_assessed`. The donor spread of a
#' batch pair is the median over the donor pairs of both batches.
#'
#' **Units.** Donors and reference aliquots are the replicates. Events are the
#' measurements: the EMDs describe distributions, and no p-value is made from
#' the number of cells.
#'
#' **No roll-up.** Rows are per marker, reference id and batch pair. There is
#' no pass or fail and no summary over markers (Q14).
#'
#' **No defaults.** Events are drawn per file without replacement with `seed`.
#' A file with fewer events than `n_events` is an error. The caller's RNG state
#' is restored.
#'
#' @param data A named list of numeric matrices or data frames, one per file,
#'   one column per marker. Names are sample ids and must be in
#'   `manifest$sample_id`.
#' @param manifest A [sw_sample_manifest()] result, or a data frame with the
#'   columns `sample_id`, `batch`, `donor`, `reference_role` and
#'   `reference_id`.
#' @param markers Character vector of columns to score.
#' @param n_events Whole number, events drawn per file (uncalibrated: no
#'   sourced value).
#' @param seed Whole number for the draws.
#' @param reference_roles Character vector of `reference_role` values that mark
#'   reference samples, for example `c("technical_replicate", "anchor")`.
#'   Samples with any other role and a donor are study samples.
#' @param min_donors_per_batch Whole number >= 2. Fewest donors a batch needs
#'   for a donor spread (uncalibrated: the plan cites 4 from a dataset rule that
#'   was not verified).
#' @param noise_floor `NULL` (skip, `floor_emd` is `NA`), or a data frame with
#'   columns `reference_id`, `marker` and `floor_emd`, one row per key, for
#'   example the sampling-only floor of a reference. The unit is the EMD unit
#'   of the transformed scale (uncalibrated).
#' @param transform,cofactor Text note and number that describe the scale of
#'   `data`, for example `"asinh"` and `6000`. Recorded only; use `NA_real_`
#'   when the transform has no cofactor.
#'
#' @return A tibble of class `sw_ref_vs_study` with columns `marker`,
#'   `reference_id`, `batch_a`, `batch_b`, `ref_emd`, `donor_emd_median`,
#'   `floor_emd`, `ratio`, `n_donors_a`, `n_donors_b`, `n_donor_pairs`,
#'   `status` (`"ok"` or `"not_assessed"`) and `reason`. Attributes `n_events`,
#'   `seed`, `min_donors_per_batch`, `transform` and `cofactor` record the call.
#'
#' @seealso [sw_batch_scores()], [sw_emd_1d()], [sw_ref_composition()],
#'   [sw_ref_coverage()]
#'
#' @examples
#' set.seed(1)
#' ids <- c(paste0("r", 1:2), paste0("d", 1:8))
#' man <- data.frame(
#'   id = ids,
#'   b = c("b1", "b2", rep(c("b1", "b2"), each = 4)),
#'   donor = c(NA, NA, rep(paste0("p", 1:4), 2)),
#'   role = c("anchor", "anchor", rep("none", 8)),
#'   ref = c("R", "R", rep(NA, 8))
#' )
#' man <- sw_sample_manifest(
#'   man, sample_id = "id", batch = "b", donor = "donor",
#'   reference_role = "role", reference_id = "ref"
#' )
#' files <- lapply(ids, function(i) cbind(A = stats::rnorm(200)))
#' names(files) <- ids
#' sw_ref_vs_study(
#'   files, man, markers = "A",
#'   n_events = 100, # illustrative
#'   seed = 1,
#'   reference_roles = "anchor",
#'   min_donors_per_batch = 4, # illustrative
#'   noise_floor = NULL,
#'   transform = "none", cofactor = NA_real_
#' )
#'
#' @export
sw_ref_vs_study <- function(
  data,
  manifest,
  markers,
  n_events,
  seed,
  reference_roles,
  min_donors_per_batch,
  noise_floor,
  transform,
  cofactor
) {
  .sw_stop_missing(c(
    data = missing(data),
    manifest = missing(manifest),
    markers = missing(markers),
    n_events = missing(n_events),
    seed = missing(seed),
    reference_roles = missing(reference_roles),
    min_donors_per_batch = missing(min_donors_per_batch),
    noise_floor = missing(noise_floor),
    transform = missing(transform),
    cofactor = missing(cofactor)
  ))
  .sw_check_count(n_events, "n_events")
  .sw_check_count(min_donors_per_batch, "min_donors_per_batch")
  if (min_donors_per_batch < 2) {
    stop("`min_donors_per_batch` must be at least 2.", call. = FALSE)
  }
  .sw_check_seed(seed)
  .sw_check_scale_note(transform, cofactor)
  if (!is.character(markers) || length(markers) == 0) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(reference_roles) || length(reference_roles) == 0) {
    stop(
      "`reference_roles` must be a non-empty character vector.",
      call. = FALSE
    )
  }
  man <- .sw_ref_manifest(
    manifest,
    names(data),
    c("sample_id", "batch", "donor", "reference_role", "reference_id")
  )
  is_ref <- man[["reference_role"]] %in% reference_roles
  if (!any(is_ref)) {
    stop("No sample has a role in `reference_roles`.", call. = FALSE)
  }
  if (anyNA(man[["reference_id"]][is_ref])) {
    stop("Reference samples need a `reference_id`.", call. = FALSE)
  }
  if (!is.null(noise_floor)) {
    ok <- is.data.frame(noise_floor) &&
      all(c("reference_id", "marker", "floor_emd") %in% names(noise_floor)) &&
      !anyDuplicated(noise_floor[c("reference_id", "marker")])
    if (!ok) {
      stop(
        "`noise_floor` must be NULL or a data frame with `reference_id`, ",
        "`marker` and `floor_emd`, one row per reference and marker.",
        call. = FALSE
      )
    }
  }

  # Draw events per file, in sorted file order so the seed fixes the draw.
  man <- man[is_ref | !is.na(man[["donor"]]), ]
  man <- man[order(man[["sample_id"]]), ]
  is_ref <- man[["reference_role"]] %in% reference_roles
  drawn <- .sw_ref_draw(data, man[["sample_id"]], markers, n_events, seed)

  # Pool files into one distribution per reference id and batch, and per donor
  # and batch.
  ref_rows <- man[is_ref, ]
  ref_pool <- split(
    ref_rows[["sample_id"]],
    paste(ref_rows[["reference_id"]], ref_rows[["batch"]], sep = "\r")
  ) |>
    lapply(\(ids) do.call(rbind, drawn[ids]))
  study_rows <- man[!is_ref, ]
  donor_spread <- split(study_rows, study_rows[["batch"]]) |>
    lapply(\(rows) {
      pools <- split(rows[["sample_id"]], rows[["donor"]]) |>
        lapply(\(ids) do.call(rbind, drawn[ids]))
      .sw_donor_pair_emd(pools, markers)
    })

  # One row per marker, reference id and batch pair.
  ref_meta <- unique(ref_rows[c("reference_id", "batch")])
  rows <- list()
  for (rid in unique(ref_meta[["reference_id"]])) {
    bs <- sort(ref_meta[["batch"]][ref_meta[["reference_id"]] == rid])
    if (length(bs) < 2) {
      next
    }
    for (p in utils::combn(bs, 2, simplify = FALSE)) {
      xa <- ref_pool[[paste(rid, p[[1]], sep = "\r")]]
      xb <- ref_pool[[paste(rid, p[[2]], sep = "\r")]]
      sa <- donor_spread[[p[[1]]]]
      sb <- donor_spread[[p[[2]]]]
      nda <- if (is.null(sa)) 0L else sa[["n_donors"]]
      ndb <- if (is.null(sb)) 0L else sb[["n_donors"]]
      enough <- nda >= min_donors_per_batch && ndb >= min_donors_per_batch
      for (mk in markers) {
        rows[[length(rows) + 1]] <- tibble::tibble(
          marker = mk,
          reference_id = rid,
          batch_a = p[[1]],
          batch_b = p[[2]],
          ref_emd = sw_emd_1d(xa[, mk], xb[, mk]),
          donor_emd_median = if (enough) {
            stats::median(c(sa[["emd"]][, mk], sb[["emd"]][, mk]))
          } else {
            NA_real_
          },
          n_donors_a = nda,
          n_donors_b = ndb,
          n_donor_pairs = if (enough) {
            nrow(sa[["emd"]]) + nrow(sb[["emd"]])
          } else {
            0L
          }
        )
      }
    }
  }
  if (length(rows) == 0) {
    stop("No reference id is measured in two or more batches.", call. = FALSE)
  }

  floor_tbl <- if (is.null(noise_floor)) {
    tibble::tibble(
      reference_id = character(),
      marker = character(),
      floor_emd = double()
    )
  } else {
    tibble::as_tibble(noise_floor[c("reference_id", "marker", "floor_emd")])
  }
  out <- dplyr::bind_rows(rows) |>
    dplyr::left_join(
      floor_tbl,
      by = dplyr::join_by(reference_id, marker),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      ratio = dplyr::if_else(
        !is.na(donor_emd_median) & donor_emd_median > 0,
        ref_emd / donor_emd_median,
        NA_real_
      ),
      status = dplyr::if_else(is.na(ratio), "not_assessed", "ok"),
      reason = dplyr::case_when(
        n_donors_a < min_donors_per_batch | n_donors_b < min_donors_per_batch ~
          "fewer than min_donors_per_batch donors in a batch",
        is.na(ratio) ~ "donor spread is zero",
        .default = NA_character_
      )
    ) |>
    dplyr::select(
      marker,
      reference_id,
      batch_a,
      batch_b,
      ref_emd,
      donor_emd_median,
      floor_emd,
      ratio,
      n_donors_a,
      n_donors_b,
      n_donor_pairs,
      status,
      reason
    )
  structure(
    out,
    class = c("sw_ref_vs_study", class(out)),
    n_events = n_events,
    seed = seed,
    min_donors_per_batch = min_donors_per_batch,
    transform = transform,
    cofactor = cofactor
  )
}

# Pairwise EMD between pooled donor distributions of one batch: the donor count
# and a matrix with one row per donor pair and one column per marker.
.sw_donor_pair_emd <- function(pools, markers) {
  nm <- names(pools)
  emd <- matrix(numeric(), 0, length(markers), dimnames = list(NULL, markers))
  if (length(nm) >= 2) {
    pr <- utils::combn(nm, 2, simplify = FALSE)
    emd <- do.call(
      rbind,
      lapply(pr, \(p) {
        vapply(
          markers,
          \(mk) sw_emd_1d(pools[[p[[1]]]][, mk], pools[[p[[2]]]][, mk]),
          numeric(1)
        )
      })
    )
    colnames(emd) <- markers
  }
  list(n_donors = length(nm), emd = emd)
}

# Draw `n_events` rows per file without replacement; the seed fixes the draw and
# the caller's RNG state is restored.
.sw_ref_draw <- function(data, files, markers, n_events, seed) {
  mats <- lapply(files, \(f) {
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
    lapply(mats, \(m) m[sample.int(nrow(m), n_events), , drop = FALSE])
  )
  if (anyNA(do.call(rbind, drawn))) {
    stop("Events must not contain NA in the chosen markers.", call. = FALSE)
  }
  drawn
}

# Manifest rows (as character columns `cols`) for the files in `files`. Stops on
# a missing column or a file with no manifest row.
.sw_ref_manifest <- function(manifest, files, cols) {
  if (!is.data.frame(manifest)) {
    stop("`manifest` must be a data frame.", call. = FALSE)
  }
  if (is.null(files) || anyDuplicated(files)) {
    stop(
      "`data` must be a named list of files with unique names.",
      call. = FALSE
    )
  }
  absent <- setdiff(cols, names(manifest))
  if (length(absent) > 0) {
    stop(
      "`manifest` lacks column(s): ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  man <- tibble::as_tibble(manifest)[cols]
  man[] <- lapply(man, as.character)
  if (anyDuplicated(man[["sample_id"]])) {
    stop("`manifest$sample_id` must be unique.", call. = FALSE)
  }
  if (!all(files %in% man[["sample_id"]])) {
    stop(
      "Files not in `manifest$sample_id`: ",
      paste(setdiff(files, man[["sample_id"]]), collapse = ", "),
      call. = FALSE
    )
  }
  man[man[["sample_id"]] %in% files, ]
}

.sw_check_seed <- function(seed) {
  if (
    !is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed != round(seed)
  ) {
    stop("`seed` must be one whole number.", call. = FALSE)
  }
}

.sw_check_scale_note <- function(transform, cofactor) {
  if (!is.character(transform) || length(transform) != 1 || is.na(transform)) {
    stop(
      "`transform` must be one string, for example \"asinh\".",
      call. = FALSE
    )
  }
  if (!is.numeric(cofactor) || length(cofactor) != 1) {
    stop("`cofactor` must be one number (NA_real_ if none).", call. = FALSE)
  }
}
