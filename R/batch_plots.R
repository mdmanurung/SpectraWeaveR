# R/batch_plots.R
# Dot-matrix figure of one batch-by-marker check. No cutoff is drawn.

#' Plot a batch-by-marker check as a dot matrix
#'
#' @description
#' Draw one check of [sw_batch_matrix()]: markers on the rows, batches on the
#' columns, dot colour for the value and dot size for the rank inside the
#' marker (rank 1 is the largest dot). Markers are ordered by their largest
#' ranking key, so the marker with the most unusual batch is on top. No
#' cutoff or reference line is drawn.
#'
#' @details
#' The caption text and the plot frame are built first; the plot call only
#' maps them. The ordering key follows the `direction` stored in the matrix
#' (`"diag_key"`, `"high"`, `"low"` or `"abs"`). Needs ggplot2.
#'
#' @param matrix An `sw_batch_matrix` from [sw_batch_matrix()].
#' @param check One check in `matrix` to draw.
#' @param n_markers `NULL` to draw every marker, or a whole number: keep the
#'   first `n_markers` in the order above (uncalibrated: it decides which
#'   markers appear). Required.
#' @param population `NULL` when `matrix` has no `population` column, else
#'   the one population to draw (`NA` draws the all-event rows of
#'   [sw_batch_scores()] or [sw_batch_mixing()]). Required.
#'
#' @return A list with `plot` (a ggplot) and `source_data` (the plotted
#'   frame, one row per marker and batch).
#'
#' @seealso [sw_batch_matrix()]
#'
#' @examples
#' \dontrun{
#' m <- sw_batch_matrix(sc, check = "rest_emd", direction = c(rest_emd = "high"))
#' fig <- sw_plot_batch_matrix(m, check = "rest_emd", n_markers = NULL, population = NULL)
#' fig$plot
#' }
#'
#' @export
sw_plot_batch_matrix <- function(matrix, check, n_markers, population) {
  .sw_stop_missing(c(
    matrix = missing(matrix),
    check = missing(check),
    n_markers = missing(n_markers),
    population = missing(population)
  ))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required. Install it with ",
      "install.packages('ggplot2').",
      call. = FALSE
    )
  }
  if (!inherits(matrix, "sw_batch_matrix")) {
    stop("`matrix` must come from sw_batch_matrix().", call. = FALSE)
  }
  if (
    !is.character(check) || length(check) != 1 || !check %in% matrix[["check"]]
  ) {
    stop("`check` must be one check present in `matrix`.", call. = FALSE)
  }
  if (!is.null(n_markers)) {
    .sw_check_count(n_markers, "n_markers")
  }
  has_pop <- "population" %in% names(matrix)
  if (
    has_pop && (is.null(population) || !population %in% matrix[["population"]])
  ) {
    stop(
      "`population` must be one population present in `matrix`.",
      call. = FALSE
    )
  }
  if (!has_pop && !is.null(population)) {
    stop(
      "`matrix` has no population column; pass `population = NULL`.",
      call. = FALSE
    )
  }

  # Plot frame and caption first.
  direction <- attr(matrix, "params")[["direction"]]
  frame <- tibble::as_tibble(as.data.frame(matrix))
  frame <- frame[frame[["check"]] == check, , drop = FALSE]
  if (has_pop) {
    in_pop <- if (is.na(population)) {
      is.na(frame[["population"]])
    } else {
      !is.na(frame[["population"]]) & frame[["population"]] == population
    }
    frame <- frame[in_pop, , drop = FALSE]
  }
  frame[["key"]] <- .sw_matrix_key(
    frame[["check"]],
    frame[["value"]],
    direction
  )
  frame[["rank_size"]] <- stats::ave(
    frame[["rank"]],
    paste(frame[["marker"]], frame[["reference_id"]]),
    FUN = function(r) {
      if (all(is.na(r))) rep(1, length(r)) else max(r, na.rm = TRUE) + 1 - r
    }
  )
  top <- tapply(frame[["key"]], frame[["marker"]], function(k) {
    if (all(is.na(k))) NA_real_ else max(k, na.rm = TRUE)
  })
  marker_order <- names(sort(top, decreasing = TRUE, na.last = TRUE))
  if (!is.null(n_markers)) {
    marker_order <- utils::head(marker_order, n_markers)
  }
  frame <- frame[frame[["marker"]] %in% marker_order, , drop = FALSE]
  frame[["marker"]] <- factor(frame[["marker"]], levels = rev(marker_order))
  caption <- paste0(
    "Colour: ",
    check,
    " value. Dot size: rank within marker (largest dot = rank 1, direction '",
    direction[[check]],
    "'). Markers ordered by their largest ranking key. No cutoff drawn."
  )
  title <- if (has_pop) paste0(check, ", population ", population) else check

  plot <- ggplot2::ggplot(
    frame,
    ggplot2::aes(x = batch, y = marker, colour = value, size = rank_size)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_size_area(max_size = 4, guide = "none") +
    ggplot2::scale_colour_viridis_c(name = "Value") +
    ggplot2::labs(x = "Batch", y = "Marker", title = title, caption = caption) +
    .sw_theme_nature()

  if ("reference_id" %in% names(frame)) {
    plot <- plot + ggplot2::facet_wrap(ggplot2::vars(reference_id))
  }

  list(plot = plot, source_data = frame)
}
