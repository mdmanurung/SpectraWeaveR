# R/gate_plot.R
# Plot one gate over the parent-gated events of several samples.

#' Plot a gate on parent-gated events, one panel per sample
#'
#' @description
#' Draw the events inside the gate's parent chain for each sample, with the
#' gate overlaid. A polygon gate (two channels) gives a scatter with the
#' polygon outline. A range gate (one channel) gives a density with vertical
#' bounds. When `reference_gates` is given, the reference gate is also drawn
#' (dashed) on the same events, so you can see how far an adapted gate moved.
#'
#' @details
#' Events shown are those inside the parent chain of the gate in that
#' sample's own `gates` (adapted parents when you pass adapted gates). Data
#' are plotted as they are; this function never transforms.
#'
#' **No defaults.** `n_events` and `seed` change which events show, so both
#' must be passed. Events are drawn per sample without replacement; a sample
#' with fewer events shows all. The caller's RNG state is restored.
#'
#' Needs ggplot2 (Suggests).
#'
#' @param data A named list of numeric matrices, one per sample.
#' @param gates A list of `sw_gate` shared by all samples, or a named list
#'   with one list per sample (for example `$gates` from [sw_gate_transfer()]).
#' @param samples Character, the sample names to draw (order kept).
#' @param gate Character, the id of the gate to draw.
#' @param n_events Whole number, the most events drawn per sample.
#' @param seed Whole number for the draws.
#' @param reference_gates `NULL`, or a list of `sw_gate` (the reference gates)
#'   to draw dashed beside `gates`.
#'
#' @return A ggplot object.
#'
#' @seealso [sw_gate_transfer()], [sw_gate_qc()]
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- list(s1 = cbind(A = stats::rnorm(500), B = stats::rnorm(500)))
#'   g <- sw_gate_polygon(
#'     "pop", c("A", "B"),
#'     rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1)),
#'     parent = NA
#'   )
#'   sw_gate_plot(d, list(g), samples = "s1", gate = "pop",
#'     n_events = 300, seed = 1, reference_gates = NULL) # illustrative
#' }
#'
#' @export
sw_gate_plot <- function(
  data,
  gates,
  samples,
  gate,
  n_events,
  seed,
  reference_gates
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for sw_gate_plot().", call. = FALSE)
  }
  data <- .sw_gq_check_data(data)
  .sw_check_count(n_events, "n_events")
  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed != round(seed)) {
    stop("`seed` must be one whole number.", call. = FALSE)
  }
  if (!is.character(samples) || length(samples) == 0 || !all(samples %in% names(data))) {
    stop("`samples` must name samples present in `data`.", call. = FALSE)
  }
  if (!is.character(gate) || length(gate) != 1) stop("`gate` must be one gate id.", call. = FALSE)
  per_sample <- .sw_gq_per_sample(gates, samples)
  if (!all(vapply(per_sample, function(gl) gate %in% names(gl), logical(1)))) {
    stop("Gate '", gate, "' is not in `gates` for every sample.", call. = FALSE)
  }
  ref_gate <- NULL
  if (!is.null(reference_gates)) {
    rg <- .sw_gq_order(reference_gates)
    if (!gate %in% names(rg)) stop("Gate '", gate, "' is not in `reference_gates`.", call. = FALSE)
    ref_gate <- rg[[gate]]
  }
  g0 <- per_sample[[1]][[gate]]
  chans <- g0$channels
  for (gl in per_sample) {
    if (!identical(gl[[gate]]$channels, chans)) {
      stop("Gate '", gate, "' must have the same channels in all samples.", call. = FALSE)
    }
  }

  events <- .sw_with_seed(seed, {
    lapply(samples, function(s) {
      m <- data[[s]]
      pm <- .sw_gq_masks(per_sample[[s]], m)[[gate]]$parent
      e <- m[pm, chans, drop = FALSE]
      if (nrow(e) > n_events) e <- e[sample.int(nrow(e), n_events), , drop = FALSE]
      d <- as.data.frame(e)
      colnames(d) <- c("x", "y")[seq_along(chans)]
      d$sample <- rep(s, nrow(d))
      d
    })
  })
  ev <- do.call(rbind, events)
  ev$sample <- factor(ev$sample, levels = samples)

  outline <- function(g, label) {
    if (g$type == "polygon") {
      d <- as.data.frame(g$vertices)
      colnames(d) <- c("x", "y")
      d <- d[c(seq_len(nrow(d)), 1L), , drop = FALSE]
      d
    } else {
      data.frame(lo = g$vertices["min", 1], hi = g$vertices["max", 1])
    }
  }
  adapted_df <- do.call(rbind, lapply(samples, function(s) {
    d <- outline(per_sample[[s]][[gate]])
    d$sample <- rep(s, nrow(d))
    d
  }))
  adapted_df$sample <- factor(adapted_df$sample, levels = samples)
  ref_df <- NULL
  if (!is.null(ref_gate)) {
    ref_df <- do.call(rbind, lapply(samples, function(s) {
      d <- outline(ref_gate)
      d$sample <- rep(s, nrow(d))
      d
    }))
    ref_df$sample <- factor(ref_df$sample, levels = samples)
  }

  if (g0$type == "polygon") {
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = ev,
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        size = 0.2, alpha = 0.4, stroke = 0
      ) +
      ggplot2::geom_path(
        data = adapted_df,
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        colour = "#b2182b", linewidth = 0.6
      )
    if (!is.null(ref_df)) {
      p <- p + ggplot2::geom_path(
        data = ref_df,
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        colour = "black", linewidth = 0.5, linetype = "dashed"
      )
    }
  } else {
    bounds <- function(d) {
      d <- data.frame(sample = rep(d$sample, 2), x = c(d$lo, d$hi))
      d[is.finite(d$x), , drop = FALSE]
    }
    p <- ggplot2::ggplot() +
      ggplot2::geom_density(
        data = ev,
        mapping = ggplot2::aes(x = .data$x),
        fill = "grey80", colour = "grey30"
      ) +
      ggplot2::geom_vline(
        data = bounds(adapted_df),
        mapping = ggplot2::aes(xintercept = .data$x), colour = "#b2182b", linewidth = 0.6
      )
    if (!is.null(ref_df)) {
      p <- p + ggplot2::geom_vline(
        data = bounds(ref_df),
        mapping = ggplot2::aes(xintercept = .data$x), colour = "black",
        linewidth = 0.5, linetype = "dashed"
      )
    }
  }
  p +
    ggplot2::facet_wrap(~sample) +
    ggplot2::labs(title = paste0("Gate ", gate), x = chans[1], y = if (length(chans) == 2) chans[2] else "density") +
    ggplot2::theme_bw()
}
