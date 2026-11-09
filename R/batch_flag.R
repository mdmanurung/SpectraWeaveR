# R/batch_flag.R
# Control-sample batch flags: an IQR (Tukey fence) assessment across batches.
# Clean-room implementation from a written method description; no upstream code.

#' Flag batches from a control sample with an IQR fence
#'
#' @description
#' Use a control sample (the same reference material measured in every batch)
#' to flag batches whose control looks different from the rest.
#'
#' For each control, marker and batch, three values are computed:
#'
#' * `neg_median`: median intensity of the negative cells (`<= threshold`).
#' * `pos_median`: median intensity of the positive cells (`> threshold`).
#' * `pct_positive`: percent of events above the threshold.
#'
#' For each control, marker and metric, the per-batch values form one panel.
#' Quartiles are taken across batches. A batch is flagged when its value lies
#' below `Q1 - k * IQR` or above `Q3 + k * IQR`, with `k = iqr_multiplier`.
#'
#' @details
#' **Scale.** `cells` must already be on the analysis scale (for example
#' arcsinh-transformed). The scale is the caller's choice. `thresholds` must be
#' on that same scale.
#'
#' **No defaults for the method knobs.** `iqr_multiplier`, `min_cells_per_side`
#' and `quantile_type` have no default because each one changes which batches
#' are flagged. The source methods use `iqr_multiplier = 1.5` and
#' `min_cells_per_side = 100`. Those are conventions, not validated cutoffs.
#' `quantile_type` matters most with 4 to 6 batches, where R's quantile types
#' give visibly different fences.
#'
#' **No automatic thresholds.** Positivity thresholds are analyst-supplied.
#' Different analysts pick different thresholds for the same marker, and the
#' threshold moves all three metrics.
#'
#' **Eligibility.** The fence is applied per panel and only when every batch
#' contributes a value. A panel is not assessed when any of these hold:
#'
#' * fewer than 4 batches exist in `cells` (`ineligible_n_batches`); with 3
#'   or fewer values the fence can never be exceeded;
#' * the control is absent from at least one batch
#'   (`ineligible_missing_control`);
#' * at least one batch has fewer than `min_cells_per_side` events on the
#'   side the metric uses (`ineligible_min_cells`). The side is the negative
#'   cells for `neg_median`, the positive cells for `pos_median` and all
#'   events for `pct_positive`;
#' * a value is `NA` after those checks, for example an empty side with
#'   `min_cells_per_side = 0` (`undefined`).
#'
#' When a panel is not assessed, every row of it carries that status, the
#' offending batches are named in `reason`, and `q1`, `q3`, `lower`, `upper`
#' and `flag` are `NA`. The fence is never applied to the remaining batches.
#' A large shift can push a batch's control entirely across the threshold, so
#' the short batch is often the shifted one. An `NA` flag means "not assessed",
#' never "no batch effect".
#'
#' **Repeated controls.** Event-level data cannot show two acquisitions of the
#' same control in one batch. Pass `sample_col` (one id per acquired sample)
#' to make that an error. Without it, repeated acquisitions pool silently.
#'
#' @param cells A data frame with one row per event. It holds the batch
#'   column, the control column and one numeric column per marker. Marker
#'   values must not be `NA`.
#' @param markers Character vector of marker column names to assess.
#' @param thresholds Named numeric vector of positivity thresholds, one per
#'   marker, on the same scale as `cells`. Required; there is no automatic
#'   thresholding.
#' @param iqr_multiplier Positive number `k` for the fence `Q1 - k * IQR`,
#'   `Q3 + k * IQR`. Required. The source methods use 1.5.
#' @param min_cells_per_side Non-negative whole number. A batch needs at least
#'   this many events on the side a metric uses. Required. The source methods
#'   use 100.
#' @param batch_col Name of the batch column.
#' @param control_col Name of the control-sample id column.
#' @param quantile_type Integer 1 to 9, passed as `type` to
#'   [stats::quantile()]. Required.
#' @param sample_col Optional name of a column holding one id per acquired
#'   sample. When given, a control acquired more than once in one batch is an
#'   error.
#'
#' @return A tibble of class `sw_flag_iqr`, one row per control, marker,
#'   metric and batch, with columns `control`, `marker`, `metric`, `batch`,
#'   `value`, `n_cells` (events on the side the metric uses), `q1`, `q3`,
#'   `lower`, `upper`, `flag`, `status` and `reason`. `control` and `batch`
#'   are character. Attribute `params` echoes every argument plus the batches
#'   found. Use [sw_flag_iqr_summary()] to count flags per batch or marker.
#'
#' @examples
#' set.seed(1)
#' cells <- data.frame(
#'   batch = rep(paste0("b", 1:5), each = 400),
#'   control = "ctrl1",
#'   CD3 = c(stats::rnorm(1000, 0, 0.3), stats::rnorm(1000, 3, 0.3))[
#'     sample(2000)
#'   ]
#' )
#' # Shift the positive cells of batch b5 upward.
#' shift <- cells$batch == "b5" & cells$CD3 > 1.5
#' cells$CD3[shift] <- cells$CD3[shift] + 2
#' res <- sw_flag_iqr(
#'   cells,
#'   markers = "CD3",
#'   thresholds = c(CD3 = 1.5),
#'   iqr_multiplier = 1.5,
#'   min_cells_per_side = 50,
#'   quantile_type = 7
#' )
#' res[res$flag %in% TRUE, c("batch", "marker", "metric", "value")]
#' sw_flag_iqr_summary(res, by = "batch")
#'
#' @export
sw_flag_iqr <- function(
  cells,
  markers,
  thresholds,
  iqr_multiplier,
  min_cells_per_side,
  batch_col = "batch",
  control_col = "control",
  quantile_type,
  sample_col = NULL
) {
  # Required arguments: each one changes the result, so none has a default.
  if (missing(cells)) {
    stop("`cells` is required and has no default.", call. = FALSE)
  }
  if (missing(markers)) {
    stop("`markers` is required and has no default.", call. = FALSE)
  }
  if (missing(thresholds)) {
    stop("`thresholds` is required and has no default.", call. = FALSE)
  }
  if (missing(iqr_multiplier)) {
    stop("`iqr_multiplier` is required and has no default.", call. = FALSE)
  }
  if (missing(min_cells_per_side)) {
    stop("`min_cells_per_side` is required and has no default.", call. = FALSE)
  }
  if (missing(quantile_type)) {
    stop("`quantile_type` is required and has no default.", call. = FALSE)
  }

  .sw_flag_check_args(
    cells,
    markers,
    thresholds,
    iqr_multiplier,
    min_cells_per_side,
    batch_col,
    control_col,
    quantile_type,
    sample_col
  )

  batch <- as.character(cells[[batch_col]])
  control <- as.character(cells[[control_col]])
  batches <- sort(unique(batch))
  controls <- sort(unique(control))

  # A control acquired twice in one batch would pool two samples into one value.
  if (!is.null(sample_col)) {
    .sw_flag_check_repeats(batch, control, as.character(cells[[sample_col]]))
  }

  # Per control x batch x marker: side counts and the three metrics.
  metrics <- c("neg_median", "pos_median", "pct_positive")
  stats_long <- lapply(markers, function(m) {
    one <- .sw_flag_metrics(batch, control, cells[[m]], thresholds[[m]])
    one$marker <- rep(m, nrow(one))
    one
  }) |>
    dplyr::bind_rows()

  # Full grid, so a control absent from a batch still gets a row.
  grid <- expand.grid(
    control = controls,
    marker = markers,
    metric = metrics,
    batch = batches,
    stringsAsFactors = FALSE
  )
  grid <- tibble::as_tibble(grid[
    order(
      grid$control,
      match(grid$marker, markers),
      match(grid$metric, metrics),
      grid$batch
    ),
  ])
  long <- dplyr::left_join(
    grid,
    stats_long,
    by = c("control", "marker", "metric", "batch"),
    relationship = "one-to-one"
  )
  long$present <- !is.na(long$n_cells)
  long$n_cells[!long$present] <- 0L

  # Fence and status per control x marker x metric panel.
  panel_key <- paste(long$control, long$marker, long$metric, sep = "\r")
  panels <- split(long, factor(panel_key, levels = unique(panel_key)))
  out <- lapply(panels, function(panel) {
    .sw_flag_panel(
      panel,
      n_batches = length(batches),
      iqr_multiplier = iqr_multiplier,
      min_cells_per_side = min_cells_per_side,
      quantile_type = quantile_type
    )
  }) |>
    dplyr::bind_rows()
  out <- out[, c(
    "control",
    "marker",
    "metric",
    "batch",
    "value",
    "n_cells",
    "q1",
    "q3",
    "lower",
    "upper",
    "flag",
    "status",
    "reason"
  )]

  attr(out, "params") <- list(
    markers = markers,
    thresholds = thresholds[markers],
    iqr_multiplier = iqr_multiplier,
    min_cells_per_side = min_cells_per_side,
    batch_col = batch_col,
    control_col = control_col,
    quantile_type = quantile_type,
    sample_col = sample_col,
    batches = batches,
    controls = controls
  )
  class(out) <- c("sw_flag_iqr", class(out))
  out
}

#' Count IQR flags per batch or per marker
#'
#' Count flagged, assessed and not-assessed rows of an [sw_flag_iqr()] result.
#' A row is assessed when its status is `"ok"`. Not-assessed rows are counted
#' apart, so a batch with no flags and no assessed rows is visible as such.
#'
#' @param x An `sw_flag_iqr` result.
#' @param by `"batch"` or `"marker"`.
#'
#' @return A tibble with the `by` column, `n_flagged`, `n_assessed` and
#'   `n_not_assessed`.
#'
#' @examples
#' cells <- data.frame(
#'   batch = rep(paste0("b", 1:4), each = 4),
#'   control = "c1",
#'   CD4 = rep(c(0, 0, 2, 2), 4)
#' )
#' res <- sw_flag_iqr(
#'   cells, "CD4", c(CD4 = 1),
#'   iqr_multiplier = 1.5, min_cells_per_side = 1, quantile_type = 7
#' )
#' sw_flag_iqr_summary(res, by = "marker")
#'
#' @export
sw_flag_iqr_summary <- function(x, by) {
  if (!inherits(x, "sw_flag_iqr")) {
    stop("`x` must be an sw_flag_iqr result.", call. = FALSE)
  }
  if (missing(by) || length(by) != 1L || !by %in% c("batch", "marker")) {
    stop("`by` must be \"batch\" or \"marker\".", call. = FALSE)
  }
  # One factor over the keys, so every count comes back in the same order.
  by_key <- factor(x[[by]], levels = sort(unique(x[[by]])))
  assessed <- x[["status"]] == "ok"
  flagged <- x[["flag"]] %in% TRUE
  res <- tibble::tibble(
    key = levels(by_key),
    n_flagged = as.integer(tapply(flagged, by_key, sum)),
    n_assessed = as.integer(tapply(assessed, by_key, sum)),
    n_not_assessed = as.integer(tapply(!assessed, by_key, sum))
  )
  names(res)[1] <- by
  res
}

#' @rdname sw_flag_iqr_summary
#' @param object An `sw_flag_iqr` result.
#' @param ... Unused.
#' @return `summary()` returns a list with `by_batch`, `by_marker` and
#'   `status` (row counts per status).
#' @export
summary.sw_flag_iqr <- function(object, ...) {
  list(
    by_batch = sw_flag_iqr_summary(object, by = "batch"),
    by_marker = sw_flag_iqr_summary(object, by = "marker"),
    status = tibble::tibble(
      status = sort(unique(object[["status"]])),
      n = as.integer(table(object[["status"]])[sort(unique(object[[
        "status"
      ]]))])
    )
  )
}

# Internal helpers -------------------------------------------------------------

.sw_flag_check_args <- function(
  cells,
  markers,
  thresholds,
  iqr_multiplier,
  min_cells_per_side,
  batch_col,
  control_col,
  quantile_type,
  sample_col
) {
  if (!is.data.frame(cells)) {
    stop("`cells` must be a data frame.", call. = FALSE)
  }
  if (
    !is.character(markers) || length(markers) == 0L || anyDuplicated(markers)
  ) {
    stop(
      "`markers` must be a non-empty character vector without repeats.",
      call. = FALSE
    )
  }
  needed <- c(batch_col, control_col, markers, sample_col)
  absent <- setdiff(needed, names(cells))
  if (length(absent) > 0L) {
    stop(
      sprintf(
        "Columns missing from `cells`: %s.",
        paste(absent, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!is.numeric(thresholds) || is.null(names(thresholds))) {
    stop("`thresholds` must be a named numeric vector.", call. = FALSE)
  }
  no_thr <- setdiff(markers, names(thresholds))
  if (length(no_thr) > 0L) {
    stop(
      sprintf(
        "No threshold for marker(s): %s.",
        paste(no_thr, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (anyNA(thresholds[markers])) {
    stop("`thresholds` must not be NA for an assessed marker.", call. = FALSE)
  }
  if (
    !is.numeric(iqr_multiplier) ||
      length(iqr_multiplier) != 1L ||
      !is.finite(iqr_multiplier) ||
      iqr_multiplier <= 0
  ) {
    stop("`iqr_multiplier` must be one positive finite number.", call. = FALSE)
  }
  if (
    !is.numeric(min_cells_per_side) ||
      length(min_cells_per_side) != 1L ||
      !is.finite(min_cells_per_side) ||
      min_cells_per_side < 0 ||
      min_cells_per_side != round(min_cells_per_side)
  ) {
    stop(
      "`min_cells_per_side` must be one non-negative whole number.",
      call. = FALSE
    )
  }
  if (length(quantile_type) != 1L || !quantile_type %in% 1:9) {
    stop("`quantile_type` must be one integer from 1 to 9.", call. = FALSE)
  }
  if (anyNA(cells[[batch_col]]) || anyNA(cells[[control_col]])) {
    stop("Batch and control columns must not contain NA.", call. = FALSE)
  }
  # NA intensities would miscount sides and make medians NA, so they fail here.
  for (m in markers) {
    if (!is.numeric(cells[[m]])) {
      stop(sprintf("Marker column `%s` is not numeric.", m), call. = FALSE)
    }
    n_na <- sum(is.na(cells[[m]]))
    if (n_na > 0L) {
      stop(
        sprintf("Marker `%s` has %d NA values; resolve them first.", m, n_na),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

.sw_flag_check_repeats <- function(batch, control, sample) {
  ids <- unique(data.frame(batch = batch, control = control, sample = sample))
  n_per_pair <- table(paste(ids$control, "in", ids$batch))
  pairs <- n_per_pair[n_per_pair > 1L]
  if (length(pairs) > 0L) {
    stop(
      sprintf(
        "Control acquired more than once in one batch: %s.",
        paste0(
          names(pairs),
          " (",
          as.integer(pairs),
          " samples)",
          collapse = "; "
        )
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Side counts and the three metrics for one marker, long by metric.
.sw_flag_metrics <- function(batch, control, x, threshold) {
  pos <- x > threshold
  # Row indices per control x batch, unnamed so results carry no names.
  idx <- unname(split(
    seq_along(x),
    list(control, batch),
    drop = TRUE,
    sep = "\r"
  ))
  per_batch <- tibble::tibble(
    control = vapply(idx, function(i) control[i[1]], character(1)),
    batch = vapply(idx, function(i) batch[i[1]], character(1)),
    n_total = lengths(idx),
    n_pos = vapply(idx, function(i) sum(pos[i]), integer(1)),
    n_neg = vapply(idx, function(i) sum(!pos[i]), integer(1)),
    # median() of an empty vector is NA, never 0.
    neg_median = vapply(
      idx,
      function(i) stats::median(x[i][!pos[i]]),
      numeric(1)
    ),
    pos_median = vapply(
      idx,
      function(i) stats::median(x[i][pos[i]]),
      numeric(1)
    ),
    pct_positive = vapply(idx, function(i) 100 * mean(pos[i]), numeric(1))
  )
  dplyr::bind_rows(
    tibble::tibble(
      control = per_batch$control,
      batch = per_batch$batch,
      metric = "neg_median",
      value = per_batch$neg_median,
      n_cells = per_batch$n_neg
    ),
    tibble::tibble(
      control = per_batch$control,
      batch = per_batch$batch,
      metric = "pos_median",
      value = per_batch$pos_median,
      n_cells = per_batch$n_pos
    ),
    tibble::tibble(
      control = per_batch$control,
      batch = per_batch$batch,
      metric = "pct_positive",
      value = per_batch$pct_positive,
      n_cells = per_batch$n_total
    )
  )
}

# Status, reason and fence for one control x marker x metric panel.
.sw_flag_panel <- function(
  panel,
  n_batches,
  iqr_multiplier,
  min_cells_per_side,
  quantile_type
) {
  absent <- panel$batch[!panel$present]
  short <- panel$batch[panel$present & panel$n_cells < min_cells_per_side]
  undefined <- panel$batch[panel$present & is.na(panel$value)]

  # Checks run in priority order; the first that fails names the panel status.
  if (n_batches < 4L) {
    status <- "ineligible_n_batches"
    reason <- sprintf(
      "%d batches; the IQR fence needs at least 4 values to be exceedable.",
      n_batches
    )
  } else if (length(absent) > 0L) {
    status <- "ineligible_missing_control"
    reason <- sprintf(
      "Control absent from batch(es): %s.",
      paste(absent, collapse = ", ")
    )
  } else if (length(short) > 0L) {
    status <- "ineligible_min_cells"
    reason <- sprintf(
      "Fewer than %d events on the side this metric uses in batch(es): %s.",
      as.integer(min_cells_per_side),
      paste(short, collapse = ", ")
    )
  } else if (length(undefined) > 0L) {
    status <- "undefined"
    reason <- sprintf(
      "Metric undefined in batch(es): %s.",
      paste(undefined, collapse = ", ")
    )
  } else {
    status <- "ok"
    reason <- NA_character_
  }

  panel$status <- status
  panel$reason <- reason
  panel$q1 <- NA_real_
  panel$q3 <- NA_real_
  panel$lower <- NA_real_
  panel$upper <- NA_real_
  panel$flag <- NA

  if (status == "ok") {
    q <- stats::quantile(
      panel$value,
      c(0.25, 0.75),
      type = quantile_type,
      names = FALSE
    )
    iqr <- q[2] - q[1]
    panel$q1 <- q[1]
    panel$q3 <- q[2]
    panel$lower <- q[1] - iqr_multiplier * iqr
    panel$upper <- q[2] + iqr_multiplier * iqr
    panel$flag <- panel$value < panel$lower | panel$value > panel$upper
  }
  panel
}
