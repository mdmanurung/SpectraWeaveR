# R/batch_scores_population.R
# Population-level rows for sw_batch_scores() and the batch x population
# interaction check. Written from scratch from standard definitions.

# Check the `population` argument of sw_batch_scores().
.sw_check_population <- function(population, files, data) {
  if (
    !is.list(population) ||
      !all(c("labels", "label_source") %in% names(population))
  ) {
    stop(
      "`population` must be NULL or list(labels = <named list>, ",
      "label_source = \"gate\" or \"cluster\").",
      call. = FALSE
    )
  }
  src <- population[["label_source"]]
  if (
    !is.character(src) || length(src) != 1 || !src %in% c("gate", "cluster")
  ) {
    stop(
      "`population$label_source` must be \"gate\" or \"cluster\".",
      call. = FALSE
    )
  }
  lab <- population[["labels"]]
  if (!is.list(lab) || is.null(names(lab)) || !all(files %in% names(lab))) {
    stop(
      "`population$labels` must be a named list with one vector per file.",
      call. = FALSE
    )
  }
  for (f in files) {
    if (length(lab[[f]]) != nrow(as.matrix(data[[f]]))) {
      stop(
        "Population labels of file '",
        f,
        "' must have one label per event.",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

# Add population rows to the all-event rows `base`. `mats` holds each file's
# events (marker columns only), in the order of `files`.
.sw_batch_scores_population <- function(
  base,
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
) {
  labels <- lapply(population[["labels"]][files], as.character)
  src <- population[["label_source"]]
  pops <- sort(unique(unlist(labels)))
  names(batch) <- files
  summ <- if (group_summary == "mean") base::mean else stats::median

  # Events per file and population, for the share rows and the thin-cell rows.
  counts <- lapply(labels, function(l) {
    table(factor(l[!is.na(l)], levels = pops))
  })
  # One draw per file and population; NULL where the cell is too thin.
  pop_draw <- .sw_with_seed(
    seed,
    lapply(files, function(f) {
      cells <- lapply(pops, function(p) {
        idx <- which(labels[[f]] == p)
        if (length(idx) < min_events) {
          return(NULL)
        }
        take <- idx[sample.int(length(idx), min(n_events, length(idx)))]
        mats[[f]][take, , drop = FALSE]
      })
      names(cells) <- pops
      cells
    })
  )
  names(pop_draw) <- files
  if (anyNA(unlist(pop_draw))) {
    stop("Events must not contain NA in the chosen markers.", call. = FALSE)
  }

  rows_for <- function(p) {
    drawn_p <- lapply(pop_draw, `[[`, p)
    assessed <- files[!vapply(drawn_p, is.null, logical(1))]
    share <- vapply(
      files,
      function(f) {
        tot <- sum(counts[[f]])
        if (tot > 0) counts[[f]][[p]] / tot else NA_real_
      },
      numeric(1)
    )
    thin <- setdiff(files, assessed)
    thin_rows <- tibble::tibble(
      check = rep("cell_events", length(thin)),
      marker = NA_character_,
      batch = unname(batch[thin]),
      value = NA_real_,
      n_files = 0L,
      file = thin,
      n_events = as.integer(vapply(
        thin,
        function(f) as.numeric(counts[[f]][[p]]),
        numeric(1)
      )),
      status = "not_assessed"
    )

    # Scores for one set of files: all files, or the files of one group.
    one_set <- function(all_keep) {
      keep <- intersect(assessed, all_keep)
      parts <- list()
      if (length(unique(batch[keep])) >= 2) {
        parts$core <- .sw_batch_scores_core(
          drawn_p[keep],
          keep,
          unname(batch[keep]),
          markers,
          order,
          extra = !is.null(order)
        )
      }
      bs <- sort(unique(batch[all_keep]))
      if (length(bs) >= 2) {
        bshare <- vapply(
          bs,
          function(b) stats::median(share[all_keep][batch[all_keep] == b]),
          numeric(1)
        )
        nf <- vapply(bs, function(b) sum(batch[all_keep] == b), numeric(1))
        parts$share <- tibble::tibble(
          check = rep(c("pop_share", "share_rest_z"), each = length(bs)),
          marker = NA_character_,
          batch = rep(bs, 2),
          value = c(unname(bshare), .sw_robust_z(unname(bshare))),
          n_files = as.integer(rep(nf, 2))
        )
      }
      list(
        rows = dplyr::bind_rows(parts),
        scored = !is.null(parts[["core"]]),
        n_events = sum(vapply(drawn_p[keep], nrow, integer(1)))
      )
    }

    sets <- if (is.null(grp)) {
      list(one_set(files))
    } else {
      lapply(unique(grp), function(g) one_set(files[grp == g]))
    }
    rows <- dplyr::bind_rows(lapply(sets, `[[`, "rows"))
    n_ev <- as.integer(sum(vapply(sets, `[[`, numeric(1), "n_events")))
    if (!is.null(grp) && nrow(rows) > 0) {
      # Combine the groups per row.
      rows <- rows |>
        dplyr::summarise(
          value = summ(value),
          n_files = sum(n_files),
          .by = dplyr::all_of(c("check", "marker", "batch"))
        )
    }
    rows[["status"]] <- rep("assessed", nrow(rows))
    # No set had two batches with enough events: say so in one row.
    if (!any(vapply(sets, `[[`, logical(1), "scored"))) {
      rows <- dplyr::bind_rows(
        rows,
        tibble::tibble(
          check = "population_scores",
          marker = NA_character_,
          batch = NA_character_,
          value = NA_real_,
          n_files = length(assessed),
          status = "not_assessed"
        )
      )
    }
    rows[["file"]] <- NA_character_
    rows[["n_events"]] <- n_ev
    out <- dplyr::bind_rows(rows, thin_rows)
    out[["population"]] <- p
    out
  }

  pop_rows <- dplyr::bind_rows(lapply(pops, rows_for))
  pop_rows[["label_source"]] <- src
  pop_rows[["eligibility_cap"]] <- if (src == "cluster") {
    "eligible_for_review"
  } else {
    NA_character_
  }
  base[["population"]] <- NA_character_
  base[["file"]] <- NA_character_
  base[["n_events"]] <- NA_integer_
  base[["status"]] <- "assessed"
  base[["label_source"]] <- NA_character_
  base[["eligibility_cap"]] <- NA_character_
  out <- dplyr::bind_rows(base, pop_rows) |>
    dplyr::relocate(
      population,
      file,
      n_events,
      status,
      label_source,
      eligibility_cap,
      .after = n_files
    )
  attr(out, "label_source") <- src
  attr(out, "group_summary") <- group_summary
  out
}

#' Does a marker's batch shift differ between populations
#'
#' @description
#' Compare, per marker and batch, the score on all events with its spread over
#' populations. A ratio near 0 means the populations move together, so an
#' additive per-marker correction fits. A large ratio means the shift depends on
#' the population, so a correction fitted on all events can over- or
#' under-correct some of them. The output is numbers only; nothing is flagged.
#'
#' @details
#' `shift_pooled` is the all-event value of `which_check`.
#' `shift_sd_across_populations` is the standard deviation, over assessed
#' populations, of the value of `which_check` inside each population. `ratio` is
#' the second divided by the first, `NA` when the pooled value is 0 or missing.
#' EMD checks are unsigned, so the ratio does not say in which direction
#' populations differ. Populations are the units here, not cells.
#'
#' @param scores A result of [sw_batch_scores()] called with `population`.
#' @param which_check One check name with a marker and a batch per row, for
#'   example `"rest_emd"` (uncalibrated: the choice moves the numbers).
#' @param min_populations Whole number, at least 2: the fewest assessed
#'   populations a marker and batch need (uncalibrated: no data-derived value).
#'   Fewer gives `status = "not_assessed"` and `NA` values.
#'
#' @return A tibble with `marker`, `batch`, `check`, `shift_pooled`,
#'   `shift_sd_across_populations`, `ratio`, `n_populations` and `status`.
#'
#' @seealso [sw_batch_scores()]
#'
#' @examples
#' set.seed(1)
#' ids <- paste0("f", 1:6)
#' files <- lapply(ids, function(i) {
#'   matrix(stats::rnorm(800), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' })
#' names(files) <- ids
#' labs <- lapply(files, function(m) rep(c("p1", "p2"), each = nrow(m) / 2))
#' sc <- sw_batch_scores(
#'   files,
#'   batch = stats::setNames(rep(c("b1", "b2", "b3"), 2), ids),
#'   markers = c("A", "B"),
#'   n_events = 100, # illustrative
#'   seed = 1,
#'   order = NULL,
#'   reference_group = NULL,
#'   population = list(labels = labs, label_source = "gate"),
#'   min_events = 50, # illustrative
#'   group_summary = "median"
#' )
#' sw_batch_interaction(sc, which_check = "rest_emd", min_populations = 2)
#'
#' @export
sw_batch_interaction <- function(scores, which_check, min_populations) {
  .sw_stop_missing(c(
    scores = missing(scores),
    which_check = missing(which_check),
    min_populations = missing(min_populations)
  ))
  if (
    !inherits(scores, "sw_batch_scores") || !"population" %in% names(scores)
  ) {
    stop(
      "`scores` must come from sw_batch_scores() called with `population`.",
      call. = FALSE
    )
  }
  if (!is.character(which_check) || length(which_check) != 1) {
    stop("`which_check` must be one check name.", call. = FALSE)
  }
  .sw_check_count(min_populations, "min_populations")
  if (min_populations < 2) {
    stop("`min_populations` must be at least 2.", call. = FALSE)
  }
  rows <- scores[scores[["check"]] == which_check, ]
  if (
    nrow(rows) == 0 ||
      all(is.na(rows[["batch"]])) ||
      all(is.na(rows[["marker"]]))
  ) {
    stop(
      "`which_check` '",
      which_check,
      "' has no marker x batch rows in `scores`.",
      call. = FALSE
    )
  }
  pooled <- rows |>
    dplyr::filter(is.na(population)) |>
    dplyr::select(marker, batch, shift_pooled = value)
  by_pop <- rows |>
    dplyr::filter(!is.na(population), status == "assessed", !is.na(value)) |>
    dplyr::summarise(
      shift_sd_across_populations = if (dplyr::n() >= 2) {
        stats::sd(value)
      } else {
        NA_real_
      },
      n_populations = dplyr::n(),
      .by = dplyr::all_of(c("marker", "batch"))
    )
  pooled |>
    dplyr::left_join(
      by_pop,
      by = dplyr::join_by(marker, batch),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      n_populations = dplyr::coalesce(n_populations, 0L),
      enough = n_populations >= min_populations,
      shift_sd_across_populations = dplyr::if_else(
        enough,
        shift_sd_across_populations,
        NA_real_
      ),
      ratio = dplyr::if_else(
        enough & !is.na(shift_pooled) & shift_pooled != 0,
        shift_sd_across_populations / shift_pooled,
        NA_real_
      ),
      status = dplyr::if_else(enough, "assessed", "not_assessed"),
      check = which_check
    ) |>
    dplyr::select(
      marker,
      batch,
      check,
      shift_pooled,
      shift_sd_across_populations,
      ratio,
      n_populations,
      status
    )
}
