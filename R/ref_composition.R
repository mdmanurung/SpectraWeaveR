# R/ref_composition.R
# Population shares of one reference across batches, with an analytic
# multinomial floor from event counts (plan DR19). Own design; no bootstrap and
# no seed. The Hellinger idea of CytofRUV (Trussart 2020) is not used here.

#' Population shares of a reference across batches
#'
#' @description
#' For each reference id, count the events per population in each batch, give
#' the share, and put next to it the floor that counting alone produces. If the
#' same material is measured in every batch, shares differ only by the floor
#' plus real technical or biological differences.
#'
#' @details
#' **Floor.** With `N` events in a reference and batch, and a population share
#' `p`, the multinomial standard error is `sqrt(p * (1 - p) / N)`. This is
#' computed from the plug-in `p`, so it is rough for shares near 0 or 1. It is
#' a sampling-only floor: events from one tube are not independent aliquots
#' (cells are not replicates), so real aliquot spread is larger. There is no
#' bootstrap and no seed.
#'
#' **Labels are yours.** `labels` holds one population label per event. Gate
#' labels are stable across batches. Cluster labels depend on the batch, so a
#' result built on them is for review only.
#'
#' **Not assessed.** A population with fewer than `min_population_events`
#' events summed over all batches of the reference gets
#' `status = "not_assessed"` and `NA` for `cv_across_batches`. Its share and
#' floor are still reported.
#'
#' **No roll-up.** Rows are per reference id, population and batch. There is
#' no pass or fail and no summary over populations (Q14).
#'
#' @param labels A named list, one element per file. Each element is a vector
#'   with one population label per event. Names are sample ids and must be in
#'   `manifest$sample_id`. Missing labels are an error.
#' @param manifest A [sw_sample_manifest()] result, or a data frame with the
#'   columns `sample_id`, `batch`, `reference_role` and `reference_id`.
#' @param reference_roles Character vector of `reference_role` values that
#'   mark reference samples. Other files are ignored.
#' @param min_population_events Whole number, fewest events a population needs
#'   over all batches of one reference (uncalibrated: no sourced value).
#' @param label_source Text note on where the labels come from, for example
#'   `"gate"` or `"kmeans"`. Recorded only.
#'
#' @return A tibble of class `sw_ref_composition` with columns `reference_id`,
#'   `population`, `batch`, `share`, `floor`, `n_population_events`,
#'   `n_events`, `cv_across_batches` (standard deviation of the shares over
#'   batches divided by their mean, the same on every batch row of a
#'   population), `status` and `reason`. The attribute `label_source` records
#'   the note.
#'
#' @seealso [sw_ref_vs_study()], [sw_ref_coverage()], [sw_cluster_scores()]
#'
#' @examples
#' set.seed(1)
#' man <- data.frame(
#'   id = c("r1", "r2", "r3"),
#'   b = c("b1", "b2", "b3"),
#'   role = "anchor",
#'   ref = "R"
#' )
#' man <- sw_sample_manifest(
#'   man, sample_id = "id", batch = "b",
#'   reference_role = "role", reference_id = "ref"
#' )
#' labs <- lapply(1:3, function(i) sample(c("T", "B"), 500, replace = TRUE, prob = c(0.7, 0.3)))
#' names(labs) <- man$sample_id
#' sw_ref_composition(
#'   labs, man,
#'   reference_roles = "anchor",
#'   min_population_events = 50, # illustrative
#'   label_source = "gate"
#' )
#'
#' @export
sw_ref_composition <- function(
  labels,
  manifest,
  reference_roles,
  min_population_events,
  label_source
) {
  .sw_stop_missing(c(
    labels = missing(labels),
    manifest = missing(manifest),
    reference_roles = missing(reference_roles),
    min_population_events = missing(min_population_events),
    label_source = missing(label_source)
  ))
  .sw_check_count(min_population_events, "min_population_events")
  if (
    !is.character(label_source) ||
      length(label_source) != 1 ||
      is.na(label_source)
  ) {
    stop("`label_source` must be one string.", call. = FALSE)
  }
  if (!is.list(labels) || is.null(names(labels))) {
    stop("`labels` must be a named list of label vectors.", call. = FALSE)
  }
  if (!is.character(reference_roles) || length(reference_roles) == 0) {
    stop(
      "`reference_roles` must be a non-empty character vector.",
      call. = FALSE
    )
  }
  man <- .sw_ref_manifest(
    manifest,
    names(labels),
    c("sample_id", "batch", "reference_role", "reference_id")
  )
  man <- man[man[["reference_role"]] %in% reference_roles, ]
  if (nrow(man) == 0) {
    stop("No labelled sample has a role in `reference_roles`.", call. = FALSE)
  }
  if (anyNA(man[["reference_id"]])) {
    stop("Reference samples need a `reference_id`.", call. = FALSE)
  }
  if (anyNA(unlist(labels[man[["sample_id"]]]))) {
    stop("`labels` must not contain NA.", call. = FALSE)
  }

  # One count per reference id, batch and population (files of a batch pooled).
  # A file with no events adds no rows here; its batch still comes from `man`.
  counts <- dplyr::bind_rows(lapply(man[["sample_id"]], \(id) {
    tab <- table(as.character(labels[[id]]))
    tibble::tibble(
      sample_id = rep(id, length(tab)),
      population = names(tab),
      n_population_events = as.integer(tab)
    )
  })) |>
    dplyr::left_join(
      man[c("sample_id", "batch", "reference_id")],
      by = "sample_id",
      relationship = "many-to-one"
    ) |>
    dplyr::summarise(
      n_population_events = sum(n_population_events),
      .by = c(reference_id, batch, population)
    )

  # Full grid: every population seen for a reference, in every batch of it.
  grid <- dplyr::inner_join(
    unique(counts[c("reference_id", "population")]),
    unique(man[c("reference_id", "batch")]),
    by = "reference_id",
    relationship = "many-to-many"
  )
  out <- grid |>
    dplyr::left_join(
      counts,
      by = c("reference_id", "batch", "population"),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      n_population_events = dplyr::coalesce(n_population_events, 0L)
    ) |>
    dplyr::mutate(
      n_events = sum(n_population_events),
      .by = c(reference_id, batch)
    ) |>
    dplyr::mutate(
      share = dplyr::if_else(
        n_events > 0,
        n_population_events / n_events,
        NA_real_
      ),
      floor = sqrt(share * (1 - share) / n_events)
    ) |>
    dplyr::mutate(
      total_pop = sum(n_population_events),
      n_batches = sum(!is.na(share)),
      cv_across_batches = dplyr::if_else(
        n_batches >= 2 & mean(share, na.rm = TRUE) > 0,
        stats::sd(share, na.rm = TRUE) / mean(share, na.rm = TRUE),
        NA_real_
      ),
      .by = c(reference_id, population)
    ) |>
    dplyr::mutate(
      status = dplyr::if_else(
        total_pop < min_population_events | n_events == 0,
        "not_assessed",
        "ok"
      ),
      cv_across_batches = dplyr::if_else(
        status == "ok",
        cv_across_batches,
        NA_real_
      ),
      reason = dplyr::case_when(
        n_events == 0 ~ "no events in this batch",
        status ==
          "not_assessed" ~ "fewer than min_population_events events over batches",
        n_batches < 2 ~ "reference measured in one batch",
        .default = NA_character_
      )
    ) |>
    dplyr::arrange(reference_id, population, batch) |>
    dplyr::select(
      reference_id,
      population,
      batch,
      share,
      floor,
      n_population_events,
      n_events,
      cv_across_batches,
      status,
      reason
    )
  structure(
    out,
    class = c("sw_ref_composition", class(out)),
    label_source = label_source
  )
}
