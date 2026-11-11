# Figures F6 to F8: hexbin grid per setup and sample, and the spike-in summaries
# (plan v2 steps 14 and 15; user grill R12, R3.4, Q20, 2026-11-11). Each builder
# returns list(panel, plot, source, stem) for sw_unmix_save_figures().

# ---- Internal helpers --------------------------------------------------------

# Donor shapes: filled where possible, distinct in greyscale.
.sw_fig_donor_shapes <- c(16, 17, 15, 18, 8, 4, 3, 25)

# Setup linetypes and colours for F8; a setup without one is an error, never a
# silently undrawn curve (the 88b/88e benchmark bug).
.sw_fig_setup_linetypes <- c(
  "solid",
  "dashed",
  "dotted",
  "dotdash",
  "longdash",
  "twodash"
)
.sw_fig_setup_colours <- c(
  "#0072B2",
  "#D55E00",
  "#009E73",
  "#CC79A7",
  "#E69F00",
  "#56B4E9"
)

# Thin grid on the value axis only (user rule r-0014).
.sw_fig_value_grid_x <- function() {
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(
      colour = "grey85",
      linewidth = 0.3
    ),
    panel.grid.minor.x = ggplot2::element_line(
      colour = "grey92",
      linewidth = 0.2
    )
  )
}

# Black or white text for a strip fill, by luminance.
.sw_fig_text_on <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  ifelse(lum > 0.5, "black", "white")
}

# Donor shape scale; more donors than shapes is an error.
.sw_fig_donor_scale <- function(donors) {
  if (length(donors) > length(.sw_fig_donor_shapes)) {
    stop(
      sprintf(
        "%d donors; the figure distinguishes at most %d by shape.",
        length(donors),
        length(.sw_fig_donor_shapes)
      ),
      call. = FALSE
    )
  }
  ggplot2::scale_shape_manual(
    values = stats::setNames(.sw_fig_donor_shapes[seq_along(donors)], donors),
    name = "Donor"
  )
}

# Spike score check shared by F7 and F8; returns the requested level's summary rows.
.sw_fig_check_score <- function(score, readout_level) {
  if (!inherits(score, "sw_spike_score")) {
    stop("'score' must be an sw_spike_score.", call. = FALSE)
  }
  if (!is.numeric(readout_level) || length(readout_level) != 1L) {
    stop("'readout_level' must be one spike level.", call. = FALSE)
  }
  s <- score$summary[score$summary$level == readout_level, ]
  if (nrow(s) == 0L) {
    stop(
      sprintf(
        "No spiked cells at level %s; levels present: %s.",
        format(readout_level),
        paste(sort(unique(score$summary$level)), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  s
}

# ---- F6: hexbin grid ---------------------------------------------------------

#' F6: hexbin grid per setup and sample
#'
#' One square hexbin panel per setup (columns) and sample (rows), on the
#' arcsinh scale with equal axes. The saved height follows the grid.
#'
#' * `mode = "scatter"`: channel `y` against channel `x` from `events`, a red
#'   line at `cut` and the percentage of events above it per panel.
#' * `mode = "spike"`: one marker before and after the spike, from an
#'   `sw_spike_score`. x is the marker in the unspiked (base) unmix
#'   (`base_value`); y is the same marker in the spiked unmix
#'   (`base_value + recovered`). The dashed curve is full recovery
#'   (`asinh(sinh(x) + sinh(level))`), the solid diagonal is a lost spike, and
#'   each panel prints the percentage of cells with recovery below the score's
#'   `loss_cut`.
#'
#' Strips are filled with `setup_fill` when ggh4x is installed (the colour is
#' a key of the setups); plain otherwise.
#'
#' @param events Scatter mode: a data frame with `setup`, `sample` and the
#'   channels named by `x` and `y` (linear scale). Spike mode: an
#'   `sw_spike_score` from [sw_unmix_spike_score()].
#' @param x Scatter mode: x channel. Spike mode: the spiked marker.
#' @param y Scatter mode: y channel. Not used in spike mode (an error if
#'   given).
#' @param mode `"scatter"` or `"spike"`.
#' @param cut Cut on the arcsinh scale. Required in scatter mode.
#' @param cofactor arcsinh cofactor. In spike mode it must match the score's
#'   `k = cofactor * sinh(level)`.
#' @param setup_fill Named colours, one per setup; a setup without a colour is
#'   an error.
#' @param stem File stem for [sw_unmix_save_figures()], for example
#'   `"scatter-grid_TNF_CD8T"`.
#' @param level Spike level on the arcsinh scale. Required in spike mode; it
#'   selects the spiked cells and places the full-recovery curve.
#' @param pop Population to draw in spike mode. Required when the marker was
#'   spiked into more than one population.
#' @param bins Hexagons per axis. Display only.
#' @param limits Optional `c(min, max)` on the arcsinh scale for both axes.
#'   Display only.
#' @param font_family Font family. Display only.
#'
#' @return A list with `panel` (`"F6"`), `plot`, `source` (one row per setup
#'   and sample: `n`, `pct`), `stem` and `size` (`c(width, height)` in mm, the
#'   height set by the number of sample rows).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' ev <- data.frame(
#'   setup = rep(c("a", "b"), each = 400),
#'   sample = "d1",
#'   GLUT1 = rnorm(800, 8000, 6000),
#'   TNF = rnorm(800, 3000, 6000)
#' )
#' f6 <- sw_unmix_figure_scatter_grid(
#'   ev, "GLUT1", "TNF", mode = "scatter", cut = 1, cofactor = 6000,
#'   setup_fill = c(a = "#0072B2", b = "#D55E00"), stem = "scatter-grid_TNF"
#' )
#' }
#'
#' @export
sw_unmix_figure_scatter_grid <- function(
  events,
  x,
  y,
  mode,
  cut,
  cofactor,
  setup_fill,
  stem,
  level = NULL,
  pop = NULL,
  bins = 50,
  limits = NULL,
  font_family = "sans"
) {
  .sw_fig_require(missing(events), "events")
  .sw_fig_require(missing(x), "x")
  .sw_fig_require(missing(mode), "mode")
  .sw_fig_require(missing(cofactor), "cofactor")
  .sw_fig_require(missing(setup_fill), "setup_fill")
  .sw_fig_require(missing(stem), "stem")
  mode <- match.arg(mode, c("scatter", "spike"))
  .sw_fig_need_ggplot2()
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    stop(
      "Package 'hexbin' is required for hexbin panels. Install it with install.packages('hexbin').",
      call. = FALSE
    )
  }

  # Linear x and y per event, plus the per-panel share, for either mode.
  if (mode == "scatter") {
    .sw_fig_require(missing(y), "y")
    .sw_fig_require(missing(cut), "cut")
    absent <- setdiff(c("setup", "sample", x, y), names(events))
    if (length(absent) > 0L) {
      stop(
        sprintf(
          "`events` lacks column(s): %s.",
          paste(absent, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    ev <- tibble::tibble(
      setup = as.character(events$setup),
      sample = as.character(events$sample),
      x_lin = events[[x]],
      y_lin = events[[y]]
    )
    ev$hit <- asinh(ev$y_lin / cofactor) > cut
    x_lab <- sprintf("%s (arcsinh, cofactor %s)", x, format(cofactor))
    y_lab <- sprintf("%s (arcsinh, cofactor %s)", y, format(cofactor))
    label_fmt <- "above cut %.1f%%"
  } else {
    if (!missing(y)) {
      stop(
        "`y` is not used in spike mode: x and y are the marker before and after the spike.",
        call. = FALSE
      )
    }
    if (is.null(level)) {
      stop(
        "`level` is required in spike mode and has no default.",
        call. = FALSE
      )
    }
    ev <- .sw_fig_spike_cells(events, x, level, pop, cofactor)
    x_lab <- sprintf(
      "%s, unspiked unmix (arcsinh, cofactor %s)",
      x,
      format(cofactor)
    )
    y_lab <- sprintf("%s, +%s arcsinh spiked unmix", x, format(level))
    label_fmt <- "lost %.1f%%"
  }
  setups <- unique(ev$setup)
  no_fill <- setdiff(setups, names(setup_fill))
  if (length(no_fill) > 0L) {
    stop(
      sprintf(
        "Setup(s) without a colour in `setup_fill`: %s.",
        paste(no_fill, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Plot frame on the arcsinh scale; setups and samples in the order given.
  samples <- unique(ev$sample)
  df <- tibble::tibble(
    setup = factor(ev$setup, levels = setups),
    sample = factor(ev$sample, levels = samples),
    x = asinh(ev$x_lin / cofactor),
    y = asinh(ev$y_lin / cofactor),
    hit = ev$hit
  )
  src <- df |>
    dplyr::summarise(
      n = dplyr::n(),
      pct = 100 * mean(.data$hit),
      .by = c("setup", "sample")
    )
  lim <- if (is.null(limits)) range(c(df$x, df$y), finite = TRUE) else limits
  size <- .sw_fig_grid_size(n_cols = length(setups), n_rows = length(samples))

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_hex(bins = bins) +
    ggplot2::scale_fill_viridis_c(trans = "sqrt", name = "Events")
  if (mode == "scatter") {
    p <- p +
      ggplot2::geom_hline(yintercept = cut, colour = "red", linewidth = 0.4)
  } else {
    grid_x <- seq(lim[1], lim[2], length.out = 200)
    curves <- rbind(
      data.frame(
        x = grid_x,
        y = asinh(sinh(grid_x) + sinh(level)),
        curve = "full recovery"
      ),
      data.frame(x = grid_x, y = grid_x, curve = "spike lost")
    )
    p <- p +
      ggplot2::geom_line(
        data = curves,
        ggplot2::aes(x, y, linetype = .data$curve),
        inherit.aes = FALSE,
        colour = "grey20",
        linewidth = 0.4
      ) +
      ggplot2::scale_linetype_manual(
        values = c("full recovery" = "dashed", "spike lost" = "solid"),
        name = NULL
      )
  }
  p <- p +
    ggplot2::geom_text(
      data = src,
      ggplot2::aes(x = -Inf, y = Inf, label = sprintf(label_fmt, .data$pct)),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.3,
      size = 6,
      size.unit = "pt",
      family = font_family
    ) +
    ggplot2::coord_equal(xlim = lim, ylim = lim) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    .sw_fig_theme()

  # Strips filled with the setup colours when ggh4x is present (a data key, rule r-0010).
  if (requireNamespace("ggh4x", quietly = TRUE)) {
    fills <- unname(setup_fill[setups])
    p <- p +
      ggh4x::facet_grid2(
        sample ~ setup,
        axes = "all",
        remove_labels = "all",
        strip = ggh4x::strip_themed(
          background_x = ggh4x::elem_list_rect(fill = fills, colour = NA),
          text_x = ggh4x::elem_list_text(colour = .sw_fig_text_on(fills))
        )
      )
  } else {
    p <- p +
      ggplot2::facet_grid(sample ~ setup, axes = "all", axis.labels = "margins")
  }

  list(panel = "F6", plot = p, source = src, stem = stem, size = size)
}

# Spiked cells of one marker as linear before (base unmix) and after (spiked
# unmix) values, with `hit` = spike lost (recovery below the score's loss cut).
.sw_fig_spike_cells <- function(score, marker, level, pop, cofactor) {
  if (!inherits(score, "sw_spike_score")) {
    stop(
      "In spike mode `events` must be an sw_spike_score; x and y are the marker before and after the spike.",
      call. = FALSE
    )
  }
  cells <- score$cells[
    score$cells$marker == marker & score$cells$level == level,
  ]
  if (nrow(cells) == 0L) {
    stop(
      sprintf(
        "No spiked cells for marker '%s' at level %s.",
        marker,
        format(level)
      ),
      call. = FALSE
    )
  }
  pops <- unique(cells$pop)
  if (is.null(pop) && length(pops) > 1L) {
    stop(
      sprintf(
        "Marker '%s' was spiked into several populations (%s); choose one with `pop`.",
        marker,
        paste(pops, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!is.null(pop)) {
    cells <- cells[cells$pop == pop, ]
    if (nrow(cells) == 0L) {
      stop(sprintf("No spiked cells in population '%s'.", pop), call. = FALSE)
    }
  }
  # The score's spike size must match this cofactor, or the curve is misplaced.
  if (!isTRUE(all.equal(unique(cells$k), cofactor * sinh(level)))) {
    stop(
      sprintf(
        "The score's spike size k does not equal cofactor * sinh(level) for cofactor %s; use the cofactor the spike was built with.",
        format(cofactor)
      ),
      call. = FALSE
    )
  }
  tibble::tibble(
    setup = as.character(cells$setup),
    sample = as.character(cells$sample),
    x_lin = cells$base_value,
    y_lin = cells$base_value + cells$recovered,
    hit = cells$recovery < score$loss_cut
  )
}

# Saved size of a grid of square panels: Nature double width, height from the
# rows (panel side = free width / columns, plus strips, axes and legend room).
.sw_fig_grid_size <- function(n_cols, n_rows) {
  width <- 183
  side <- (width - 12 - 30 - 5 - 2 * (n_cols - 1)) / n_cols
  height <- n_rows * side + 2 * (n_rows - 1) + 5 + 14
  c(width, max(height, 60))
}

# ---- F7: spike summary -------------------------------------------------------

#' F7: spike-in loss and recovery per setup, donors as points
#'
#' Two readouts at one spike level, setups down the y axis (first setup on
#' top) and one row of panels per marker and population: the percentage of
#' cells losing the spike (recovery below the score's `loss_cut`) and each
#' donor's median recovery with its interquartile range. Donors are shapes; a
#' vertical bar marks the median across donors.
#'
#' @param score An `sw_spike_score` from [sw_unmix_spike_score()].
#' @param readout_level Spike level to show, as in the design's `levels`.
#' @param stem File stem for [sw_unmix_save_figures()].
#' @param font_family Font family. Display only.
#'
#' @return A list with `panel` (`"F7"`), `plot`, `source` (the plotted frame)
#'   and `stem`.
#'
#' @export
sw_unmix_figure_spike_summary <- function(
  score,
  readout_level,
  stem,
  font_family = "sans"
) {
  .sw_fig_require(missing(score), "score")
  .sw_fig_require(missing(readout_level), "readout_level")
  .sw_fig_require(missing(stem), "stem")
  .sw_fig_need_ggplot2()
  s <- .sw_fig_check_score(score, readout_level)
  setups <- unique(as.character(score$summary$setup))
  donors <- sort(unique(as.character(s$sample)))
  readouts <- c(
    lost = sprintf(
      "Cells losing the spike (%%, recovery < %s)",
      format(score$loss_cut)
    ),
    recovery = "Recovery per donor (median, IQR)"
  )

  # Long frame: one row per setup, donor, marker, population and readout;
  # setups as numbers so donors can be offset within a setup's row.
  offsets <- stats::setNames(
    seq(-0.25, 0.25, length.out = max(length(donors), 2))[seq_along(donors)],
    donors
  )
  base <- s |>
    dplyr::mutate(
      row = paste(.data$marker, "in", .data$pop),
      y = match(.data$setup, rev(setups)) + offsets[.data$sample]
    )
  df <- dplyr::bind_rows(
    base |>
      dplyr::transmute(
        .data$setup,
        .data$sample,
        .data$row,
        .data$y,
        readout = readouts[["lost"]],
        value = 100 * .data$share_lost,
        lo = NA_real_,
        hi = NA_real_
      ),
    base |>
      dplyr::transmute(
        .data$setup,
        .data$sample,
        .data$row,
        .data$y,
        readout = readouts[["recovery"]],
        value = .data$q50,
        lo = .data$q25,
        hi = .data$q75
      )
  ) |>
    dplyr::mutate(readout = factor(.data$readout, levels = unname(readouts)))
  medians <- df |>
    dplyr::summarise(
      value = stats::median(.data$value),
      .by = c("setup", "row", "readout")
    ) |>
    dplyr::mutate(y = match(.data$setup, rev(setups)))
  refs <- data.frame(
    readout = factor(readouts[["recovery"]], levels = unname(readouts)),
    xintercept = c(score$loss_cut, 1)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$value, .data$y)) +
    ggplot2::geom_vline(
      data = refs,
      ggplot2::aes(xintercept = .data$xintercept),
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.3
    ) +
    ggplot2::geom_segment(
      data = medians,
      ggplot2::aes(
        x = .data$value,
        xend = .data$value,
        y = .data$y - 0.35,
        yend = .data$y + 0.35
      ),
      linewidth = 0.6,
      colour = "grey30"
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
      na.rm = TRUE,
      linewidth = 0.3
    ) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$sample), size = 1.2) +
    .sw_fig_donor_scale(donors) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(setups),
      labels = rev(setups),
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::facet_grid(
      row ~ readout,
      scales = "free_x",
      axes = "all",
      axis.labels = "margins"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      caption = sprintf(
        "Spike level %s; bar = median across donors; dashed = loss cut and full recovery.",
        format(readout_level)
      )
    ) +
    .sw_fig_theme() +
    .sw_fig_value_grid_x()

  list(panel = "F7", plot = p, source = df, stem = stem)
}

# ---- F8: recovery ECDF -------------------------------------------------------

#' F8: cumulative distribution of per-cell spike recovery per setup
#'
#' One curve per setup (colour and linetype), one panel per marker and
#' population, at one spike level. Every setup in the score is drawn; a setup
#' without a linetype is an error. The curves are evaluated on a grid over
#' `xlim`; cells outside it still count in the cumulative share.
#'
#' @inheritParams sw_unmix_figure_spike_summary
#' @param xlim Recovery range shown, `c(min, max)`. Display only.
#'
#' @return A list with `panel` (`"F8"`), `plot`, `source` (the evaluated
#'   curves) and `stem`.
#'
#' @export
sw_unmix_figure_spike_ecdf <- function(
  score,
  readout_level,
  stem,
  xlim = c(-0.5, 1.5),
  font_family = "sans"
) {
  .sw_fig_require(missing(score), "score")
  .sw_fig_require(missing(readout_level), "readout_level")
  .sw_fig_require(missing(stem), "stem")
  .sw_fig_need_ggplot2()
  .sw_fig_check_score(score, readout_level)
  setups <- unique(as.character(score$summary$setup))
  linetypes <- stats::setNames(
    .sw_fig_setup_linetypes[seq_len(min(
      length(setups),
      length(.sw_fig_setup_linetypes)
    ))],
    setups[seq_len(min(length(setups), length(.sw_fig_setup_linetypes)))]
  )
  colours <- stats::setNames(
    .sw_fig_setup_colours[seq_along(linetypes)],
    names(linetypes)
  )
  no_line <- setdiff(setups, names(linetypes))
  if (length(no_line) > 0L) {
    stop(
      sprintf(
        "Setup(s) without a linetype: %s. F8 draws at most %d setups.",
        paste(no_line, collapse = ", "),
        length(.sw_fig_setup_linetypes)
      ),
      call. = FALSE
    )
  }

  # Curves on a grid: one ECDF per setup, marker and population.
  cells <- score$cells[score$cells$level == readout_level, ]
  grid <- seq(xlim[1], xlim[2], length.out = 201)
  df <- cells |>
    # cum_share first: a `recovery` output column would mask the cells' recovery.
    dplyr::reframe(
      cum_share = stats::ecdf(.data$recovery)(grid),
      n_cells = dplyr::n(),
      grid_x = grid,
      .by = c("setup", "marker", "pop")
    ) |>
    dplyr::rename(recovery = "grid_x") |>
    dplyr::mutate(
      row = paste(.data$marker, "in", .data$pop),
      setup = factor(.data$setup, levels = setups)
    )
  drawn <- unique(as.character(df$setup))
  if (!setequal(drawn, setups)) {
    stop(
      sprintf(
        "Setup(s) with no cells at level %s: %s.",
        format(readout_level),
        paste(setdiff(setups, drawn), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      .data$recovery,
      .data$cum_share,
      colour = .data$setup,
      linetype = .data$setup
    )
  ) +
    ggplot2::geom_vline(
      xintercept = c(score$loss_cut, 1),
      colour = "grey60",
      linewidth = 0.3,
      linetype = "dashed"
    ) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = colours, name = "Setup") +
    ggplot2::scale_linetype_manual(values = linetypes, name = "Setup") +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$row),
      axes = "all",
      axis.labels = "margins"
    ) +
    ggplot2::labs(
      x = "Recovered fraction of the spike",
      y = "Cumulative share of spiked cells",
      caption = sprintf(
        "Spike level %s; dashed = loss cut (%s) and full recovery.",
        format(readout_level),
        format(score$loss_cut)
      )
    ) +
    .sw_fig_theme() +
    .sw_fig_value_grid_x()

  list(
    panel = "F8",
    plot = p,
    source = dplyr::mutate(df, setup = as.character(.data$setup)),
    stem = stem
  )
}
