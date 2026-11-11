#' @title Unmixing Report Figures (F1 to F4)
#'
#' @description
#' Build the four figures approved in
#' `plans/beta_unmix_report_figure_plan.md`. Each builder returns the ggplot
#' object and the source table it was drawn from. Nothing is drawn or written
#' until you call [sw_unmix_save_figures()].
#'
#' * [sw_unmix_sample_events()] samples events per file for F3 and F4.
#' * [sw_unmix_figure_spectra()] is F1, the peak-normalised spectra overlay.
#' * [sw_unmix_figure_similarity()] is F2, the cosine-similarity matrix.
#' * [sw_unmix_figure_distributions()] is F3, per-channel densities.
#' * [sw_unmix_figure_biplot()] is F4, a two-channel biplot per method.
#' * [sw_unmix_save_figures()] writes PNG, PDF and source CSV per panel.
#'
#' @details
#' F5 (condition number and AF library size) stays a table in
#' [sw_unmix_report()] and has no builder here.
#'
#' Methods are coloured from one fixed palette keyed by method label
#' (Okabe-Ito colours). The accepted labels are `ours`, `spectracle`,
#' `spectroflo`, `ols`, `autospectral_fast` and `autospectral_slow`. Any other
#' label is an error in every builder.
#'
#' The figures use the package's Nature-spec theme: a prism-style
#' black-and-white theme at 7 pt, plain text with bold titles, strips and axis
#' titles, line width 0.4 and the legend on the right, so every drawn text sits
#' in 5 to 7 pt (user decision 2026-11-11, grill Q20). Faceted panels draw axis
#' lines and ticks on every panel, with tick labels on the outer edges only.
#' Saved widths are Nature column widths, 89 or 183 mm. The theme is added to
#' each plot; the session theme is not changed.
#'
#' @name unmix_figures
#' @keywords internal
NULL

# ---- Internal helpers: panels, palette, theme --------------------------------

# One named vector drives plot keys, file stems and source-data names.
.sw_fig_panels <- c(
  F1 = "spectra-overlay",
  F2 = "similarity-matrix",
  F3 = "unmixed-distributions",
  F4 = "biplot"
)

# Default saved size per panel, c(width, height) in mm. Widths are Nature column widths
# (89 single, 183 double); faceted panels take the double column.
.sw_fig_sizes <- list(
  F1 = c(183, 140),
  F2 = c(89, 85),
  F3 = c(183, 110),
  F4 = c(183, 90)
)

# Nature final artwork: column widths in mm and the drawn text range in pt.
.sw_fig_nature_widths <- c(single = 89, double = 183)
.sw_fig_nature_text_pt <- c(min = 5, max = 7)

# Okabe-Ito colours keyed by method label.
.sw_fig_palette_values <- c(
  ours = "#0072B2",
  spectracle = "#E69F00",
  spectroflo = "#009E73",
  ols = "#CC79A7",
  autospectral_fast = "#D55E00",
  autospectral_slow = "#56B4E9",
  # correction benchmark: the corrected version is a darker shade of its method's colour
  original = "#0072B2",
  corrected = "#D55E00",
  ols_corrected = "#7A3B5A",
  autospectral_fast_corrected = "#7A2E00"
)

# Palette for the given labels; errors on an unknown label.
.sw_fig_palette <- function(methods) {
  methods <- unique(as.character(methods))
  unknown <- setdiff(methods, names(.sw_fig_palette_values))
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown method label(s): %s. Accepted labels: %s.",
        paste(unknown, collapse = ", "),
        paste(names(.sw_fig_palette_values), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  .sw_fig_palette_values[methods]
}

.sw_fig_need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for unmixing figures. ",
      "Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
}

# House theme, added to each plot rather than set globally: the Nature-spec theme
# (R/theme_prism.R), legend on the right, caption left-aligned to the plot.
.sw_fig_theme <- function() {
  .sw_theme_nature() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0),
      plot.caption.position = "plot",
      # Extra right margin so the last tick label is not cut off.
      plot.margin = ggplot2::margin(5.5, 12, 5.5, 5.5)
    )
}

# Facet axes on every panel, tick labels on the outer edges only (user rule r-0008).
.sw_fig_facet_wrap <- function(var) {
  ggplot2::facet_wrap(
    ggplot2::vars({{ var }}),
    axes = "all",
    axis.labels = "margins"
  )
}

# Wrap a caption so it is not clipped at the saved width.
.sw_fig_wrap <- function(text, width = 80) {
  paste(strwrap(text, width = width), collapse = "\n")
}

# Legend on the right: one column of method keys.
.sw_fig_legend <- function(n_keys) {
  ggplot2::guide_legend(ncol = 1)
}

# Error naming the argument when a required argument is missing.
.sw_fig_require <- function(is_missing, arg) {
  if (is_missing) {
    stop(sprintf("`%s` is required and has no default.", arg), call. = FALSE)
  }
}

# Symmetric-log transform; sigma sets the width of the linear zone.
.sw_fig_pseudo_log <- function(sigma) {
  if (
    !is.numeric(sigma) ||
      length(sigma) != 1L ||
      !is.finite(sigma) ||
      sigma <= 0
  ) {
    stop("`pseudo_log_sigma` must be one positive number.", call. = FALSE)
  }
  scales::transform_pseudo_log(sigma = sigma)
}

# Minimum gap between labelled breaks, as a fraction of the axis width.
.sw_fig_min_break_gap <- 0.08

# Breaks at 0 and signed powers of ten inside the data range. Positions are
# computed through the axis transform; going outward from 0, a break closer
# than 8% of the axis width to a kept break is dropped. With sigma = 1000
# this drops +-10^3, which sit almost on 0.
.sw_fig_symlog_breaks <- function(trans) {
  function(limits) {
    limits <- limits[is.finite(limits)]
    if (length(limits) < 2L) {
      return(numeric(0))
    }
    limits <- range(limits)
    cand <- c(-10^(8:1), 0, 10^(1:8))
    cand <- cand[cand >= limits[1] & cand <= limits[2]]
    pos <- trans$transform(cand)
    min_gap <- .sw_fig_min_break_gap * diff(trans$transform(limits))
    kept <- logical(length(cand))
    # Visit 0 first, then outward by magnitude.
    for (i in order(abs(cand))) {
      if (all(abs(pos[i] - pos[kept]) >= min_gap)) {
        kept[i] <- TRUE
      }
    }
    cand[kept]
  }
}

# Labels 0, 10^5, -10^5 as plotmath; empty for a break ggplot2 dropped.
.sw_fig_symlog_labels <- function(b) {
  txt <- rep("''", length(b))
  ok <- is.finite(b)
  txt[ok] <- ifelse(
    b[ok] == 0,
    "0",
    sprintf(
      "%s10^%d",
      ifelse(b[ok] < 0, "-", ""),
      as.integer(round(log10(pmax(abs(b[ok]), 1))))
    )
  )
  parse(text = txt)
}

# Pseudo-log axis scale with the spaced breaks.
.sw_fig_symlog_scale <- function(axis, trans) {
  scale_fn <- if (axis == "x") {
    ggplot2::scale_x_continuous
  } else {
    ggplot2::scale_y_continuous
  }
  scale_fn(
    transform = trans,
    breaks = .sw_fig_symlog_breaks(trans),
    labels = .sw_fig_symlog_labels
  )
}

# Events frame with method and file columns, one file per method.
.sw_fig_check_events <- function(events, channels) {
  if (!is.data.frame(events)) {
    stop("`events` must be a data frame.", call. = FALSE)
  }
  missing_cols <- setdiff(c("method", "file", channels), names(events))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "`events` lacks column(s): %s.",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  # Pooling two files under one method would be silent, so it is an error.
  files_per_method <- tapply(
    events[["file"]],
    events[["method"]],
    \(f) length(unique(f))
  )
  if (any(files_per_method > 1L)) {
    stop(
      sprintf(
        "Each method must come from one file; more than one file for: %s.",
        paste(names(files_per_method)[files_per_method > 1L], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  .sw_fig_palette(events[["method"]])
}


# ---- Event sampling ----------------------------------------------------------

#' Sample events per file for unmixing figures F3 and F4
#'
#' Draw up to `n_events` events from each file, without replacement.
#'
#' @details
#' One `set.seed(seed)` call precedes the draws. Files are sampled in the
#' order of `x`, so the same `x`, `n_events` and `seed` give the same rows.
#' The caller's random-number state is restored afterwards, also on error.
#' A file with fewer than `n_events` events contributes all of them; the
#' `sample_log` attribute records `n_available` and `n_sampled` per file.
#'
#' The sample size and seed are impactful choices and have no default. The
#' Mendeley benchmark values are recorded in
#' `experiments/benchmarks/01_mendeley-7color/README.md`.
#'
#' Every file must have the same channel names. Channels are aligned by name
#' to the order of the first file.
#'
#' @param x A named list of numeric matrices or data frames (events x
#'   channels, with column names), one per file. Or a character vector of
#'   FCS paths, read with `flowCore::read.FCS(truncate_max_range = FALSE)`.
#' @param n_events Positive integer. Events to draw per file. Required.
#' @param seed Integer seed for the draws. Required.
#' @param method Character vector of method labels, one per file. When
#'   `NULL`, `names(x)` is used.
#'
#' @return A tibble with columns `file`, `method`, `event_id` (row index in
#'   the source file) and one column per channel. Attribute `sample_log` is a
#'   tibble with `file`, `method`, `n_available` and `n_sampled`.
#'
#' @examples
#' set.seed(10)
#' x <- list(
#'   ours = matrix(rnorm(300), ncol = 3, dimnames = list(NULL, c("A", "B", "C"))),
#'   spectroflo = matrix(rnorm(150), ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
#' )
#' ev <- sw_unmix_sample_events(x, n_events = 60, seed = 1)
#' attr(ev, "sample_log")
#'
#' @export
sw_unmix_sample_events <- function(x, n_events, seed, method = NULL) {
  .sw_fig_require(missing(x), "x")
  .sw_fig_require(missing(n_events), "n_events")
  .sw_fig_require(missing(seed), "seed")
  if (
    !is.numeric(n_events) ||
      length(n_events) != 1L ||
      !is.finite(n_events) ||
      n_events < 1 ||
      n_events != round(n_events)
  ) {
    stop("`n_events` must be one positive whole number.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("`seed` must be one finite number.", call. = FALSE)
  }

  # Read FCS paths into matrices; file ids are the path basenames.
  if (is.character(x)) {
    if (!requireNamespace("flowCore", quietly = TRUE)) {
      stop(
        "Package 'flowCore' is required to read FCS paths. ",
        "Pass a named list of matrices instead.",
        call. = FALSE
      )
    }
    file_ids <- if (is.null(names(x))) basename(x) else names(x)
    mats <- lapply(x, \(p) {
      flowCore::exprs(flowCore::read.FCS(
        p,
        transformation = FALSE,
        truncate_max_range = FALSE
      ))
    })
  } else if (is.list(x)) {
    file_ids <- names(x)
    mats <- lapply(x, as.matrix)
  } else {
    stop(
      "`x` must be a named list of matrices or a character vector of FCS paths.",
      call. = FALSE
    )
  }
  if (is.null(file_ids) || any(!nzchar(file_ids)) || anyDuplicated(file_ids)) {
    stop("`x` must have unique, non-empty names (one per file).", call. = FALSE)
  }
  if (is.null(method)) {
    method <- file_ids
  }
  if (length(method) != length(mats)) {
    stop("`method` must have one label per file in `x`.", call. = FALSE)
  }
  .sw_fig_palette(method)

  # Channels must match by name across files; order follows the first file.
  channels <- colnames(mats[[1]])
  if (is.null(channels) || any(!nzchar(channels))) {
    stop("Every file in `x` needs column (channel) names.", call. = FALSE)
  }
  for (i in seq_along(mats)) {
    cn <- colnames(mats[[i]])
    if (!setequal(cn, channels) || length(cn) != length(channels)) {
      stop(
        sprintf(
          "Channels of '%s' differ from '%s': %s.",
          file_ids[i],
          file_ids[1],
          paste(
            union(setdiff(cn, channels), setdiff(channels, cn)),
            collapse = ", "
          )
        ),
        call. = FALSE
      )
    }
    mats[[i]] <- mats[[i]][, channels, drop = FALSE]
  }

  # Save the caller's RNG state and restore it on exit, also on error.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    },
    add = TRUE
  )
  set.seed(seed)

  # Draw sorted row indices per file, in the order of `x`.
  idx <- lapply(mats, \(m) {
    n <- nrow(m)
    sort(sample.int(n, min(n, n_events)))
  })

  parts <- lapply(seq_along(mats), \(i) {
    tibble::tibble(
      file = file_ids[i],
      method = method[i],
      event_id = idx[[i]],
      tibble::as_tibble(mats[[i]][idx[[i]], , drop = FALSE])
    )
  })
  out <- dplyr::bind_rows(parts)
  attr(out, "sample_log") <- tibble::tibble(
    file = file_ids,
    method = method,
    n_available = vapply(mats, nrow, integer(1), USE.NAMES = FALSE),
    n_sampled = lengths(idx, use.names = FALSE)
  )
  out
}


# ---- F1: spectra overlay -----------------------------------------------------

#' F1: overlay reference spectra across methods
#'
#' One panel per fluorophore, detectors on x, peak-normalised intensity on y,
#' one line per method.
#'
#' @details
#' Matrices are aligned by fluorophore and detector names, never by position.
#' Every matrix must have the same fluorophore and detector sets as the
#' first. Detector order and panel order follow the first matrix.
#'
#' Each row is peak-normalised as `row / max(abs(row))`. This is a positive
#' scalar per row, so the angle is the same on raw and normalised spectra.
#' The angle between each method and the first method comes from
#' [sw_unmix_compare_spectra()] and is shown in the panel strip.
#'
#' @param spectra_list Named list of numeric matrices (fluorophores x
#'   detectors), keyed by method label. The first entry is the reference.
#'
#' @return A list with `panel` (`"F1"`), `plot` (a ggplot) and `source` (a
#'   tibble with `method`, `fluorophore`, `detector`, `value_raw`,
#'   `value_norm`, `angle_deg` and `panel_label`).
#'
#' @examples
#' det <- paste0("D", 1:6)
#' ours <- rbind(
#'   FITC = c(0.1, 1, 0.5, 0.2, 0.1, 0),
#'   PE = c(0, 0.1, 0.3, 1, 0.6, 0.2)
#' )
#' colnames(ours) <- det
#' spectroflo <- ours[c("PE", "FITC"), rev(det)] * 1.02
#' f1 <- sw_unmix_figure_spectra(list(ours = ours, spectroflo = spectroflo))
#' f1$source
#'
#' @export
sw_unmix_figure_spectra <- function(spectra_list) {
  .sw_fig_require(missing(spectra_list), "spectra_list")
  .sw_fig_need_ggplot2()
  methods <- names(spectra_list)
  if (
    !is.list(spectra_list) ||
      length(spectra_list) < 1L ||
      is.null(methods) ||
      any(!nzchar(methods)) ||
      anyDuplicated(methods)
  ) {
    stop(
      "`spectra_list` must be a list of matrices with unique method names.",
      call. = FALSE
    )
  }
  palette <- .sw_fig_palette(methods)
  for (m in methods) {
    .sw_correct_check_spectra(
      spectra_list[[m]],
      sprintf("spectra_list[[\"%s\"]]", m)
    )
  }

  # Align every matrix to the first by fluorophore and detector names.
  ref <- spectra_list[[1]]
  fluors <- rownames(ref)
  detectors <- colnames(ref)
  for (m in methods[-1]) {
    s <- spectra_list[[m]]
    if (!setequal(rownames(s), fluors) || !setequal(colnames(s), detectors)) {
      stop(
        sprintf(
          "Fluorophore or detector names of '%s' differ from '%s': %s.",
          m,
          methods[1],
          paste(
            c(
              union(setdiff(rownames(s), fluors), setdiff(fluors, rownames(s))),
              union(
                setdiff(colnames(s), detectors),
                setdiff(detectors, colnames(s))
              )
            ),
            collapse = ", "
          )
        ),
        call. = FALSE
      )
    }
    spectra_list[[m]] <- s[fluors, detectors, drop = FALSE]
  }

  # Angle of each method to the first, per fluorophore.
  angles <- lapply(methods, \(m) {
    a <- sw_unmix_compare_spectra(ref, spectra_list[[m]])
    stats::setNames(a[["angle_deg"]], a[["fluorophore"]])
  })
  names(angles) <- methods

  # Strip label: fluorophore plus the angle of every non-reference method.
  panel_label <- vapply(
    fluors,
    \(f) {
      if (length(methods) < 2L) {
        return(f)
      }
      parts <- vapply(
        methods[-1],
        \(m) sprintf("%s %.1f deg", m, angles[[m]][[f]]),
        character(1)
      )
      # One line per method keeps the strip narrow.
      paste(c(f, parts), collapse = "\n")
    },
    character(1)
  )

  # Long plot frame, one row per method x fluorophore x detector.
  source <- dplyr::bind_rows(lapply(methods, \(m) {
    s <- spectra_list[[m]]
    peak <- apply(abs(s), 1, max)
    if (any(peak == 0)) {
      stop(
        sprintf(
          "Spectra of '%s' have an all-zero row: %s.",
          m,
          paste(names(peak)[peak == 0], collapse = ", ")
        ),
        call. = FALSE
      )
    }
    norm <- s / peak
    tibble::tibble(
      method = m,
      fluorophore = rep(fluors, times = length(detectors)),
      detector = rep(detectors, each = length(fluors)),
      value_raw = as.vector(s),
      value_norm = as.vector(norm),
      angle_deg = rep(unname(angles[[m]][fluors]), times = length(detectors))
    )
  }))
  source[["panel_label"]] <- unname(panel_label[source[["fluorophore"]]])

  plot_df <- source
  plot_df[["method"]] <- factor(plot_df[["method"]], levels = methods)
  plot_df[["detector"]] <- factor(plot_df[["detector"]], levels = detectors)
  plot_df[["panel_label"]] <- factor(
    plot_df[["panel_label"]],
    levels = unname(panel_label)
  )

  # Label about 8 detectors per axis so names stay legible.
  detector_breaks <- detectors[seq(
    1L,
    length(detectors),
    by = ceiling(length(detectors) / 8)
  )]

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = detector,
      y = value_norm,
      colour = method,
      group = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::geom_point(size = 0.6) +
    .sw_fig_facet_wrap(panel_label) +
    ggplot2::scale_x_discrete(breaks = detector_breaks) +
    ggplot2::scale_colour_manual(
      values = palette,
      name = "Method",
      guide = .sw_fig_legend(length(palette))
    ) +
    ggplot2::labs(
      x = "Detector",
      y = "Peak-normalised intensity (row / max |row|)",
      caption = .sw_fig_wrap(sprintf(
        "%d fluorophores x %d detectors; angle (deg) to %s in each strip.",
        length(fluors),
        length(detectors),
        methods[1]
      ))
    ) +
    .sw_fig_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  list(panel = "F1", plot = p, source = source)
}


# ---- F2: similarity matrix ---------------------------------------------------

#' F2: cosine-similarity matrix of reference spectra
#'
#' A tile plot of the cosine similarity between every pair of fluorophores,
#' diagonal included.
#'
#' @details
#' Similarity is computed with base R as `u %*% t(u)`, where `u` holds the
#' rows of `spectra` scaled to unit length. The fill uses a sequential
#' viridis scale with limits 0 to 1. A negative cosine is an error, because
#' it would fall outside the approved scale and show as blank.
#' Fluorophores keep the row order of `spectra`.
#'
#' @param spectra Numeric matrix, fluorophores x detectors.
#'
#' @return A list with `panel` (`"F2"`), `plot`, `source` (a long tibble of
#'   the full matrix: `fluorophore_x`, `fluorophore_y`, `cosine`) and
#'   `matrix` (the square similarity matrix).
#'
#' @examples
#' s <- rbind(
#'   FITC = c(0.1, 1, 0.5, 0.2),
#'   BB515 = c(0.1, 0.9, 0.6, 0.2),
#'   PE = c(0, 0.1, 0.3, 1)
#' )
#' colnames(s) <- paste0("D", 1:4)
#' f2 <- sw_unmix_figure_similarity(s)
#' round(f2$matrix, 2)
#'
#' @export
sw_unmix_figure_similarity <- function(spectra) {
  .sw_fig_require(missing(spectra), "spectra")
  .sw_fig_need_ggplot2()
  .sw_correct_check_spectra(spectra, "spectra")
  sim <- .sw_row_cosine(spectra)
  if (anyNA(sim)) {
    stop("`spectra` has an all-zero row; cosine is undefined.", call. = FALSE)
  }
  # Rounding can push a self cosine to 1 + 2e-10, outside the 0..1 scale.
  sim <- pmin(pmax(sim, -1), 1)
  # Noisy real spectra give small negative cosines; the approved scale is 0 to 1,
  # so they are drawn as 0. The raw values stay in the source table and the
  # caption says how many were clamped.
  sim_raw <- sim
  n_clamped <- sum(sim_raw < -1e-9)
  min_raw <- min(sim_raw)
  sim <- pmax(sim, 0)
  fluors <- rownames(spectra)

  # Long form of the full matrix, diagonal included.
  source <- tibble::tibble(
    fluorophore_x = rep(fluors, times = length(fluors)),
    fluorophore_y = rep(fluors, each = length(fluors)),
    cosine = as.vector(sim),
    cosine_raw = as.vector(sim_raw)
  )

  # Caption: closest off-diagonal pair, computed from the matrix.
  off <- sim
  diag(off) <- NA_real_
  caption <- if (length(fluors) < 2L) {
    "1 fluorophore; no pairs."
  } else {
    top <- which(off == max(off, na.rm = TRUE), arr.ind = TRUE)[1, ]
    sprintf(
      "%d fluorophores; highest off-diagonal cosine %.3f (%s vs %s).",
      length(fluors),
      off[top[1], top[2]],
      fluors[top[1]],
      fluors[top[2]]
    )
  }
  if (n_clamped > 0L) {
    caption <- sprintf(
      "%s %d value(s) below 0 (min %.3f) are drawn as 0.",
      caption,
      n_clamped,
      min_raw
    )
  }

  plot_df <- source
  plot_df[["fluorophore_x"]] <- factor(plot_df[["fluorophore_x"]], fluors)
  plot_df[["fluorophore_y"]] <- factor(plot_df[["fluorophore_y"]], rev(fluors))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = fluorophore_x, y = fluorophore_y, fill = cosine)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(
      limits = c(0, 1),
      name = "Cosine\nsimilarity"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = "Fluorophore",
      y = "Fluorophore",
      caption = .sw_fig_wrap(caption, width = 70)
    ) +
    .sw_fig_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "right"
    )

  list(panel = "F2", plot = p, source = source, matrix = sim)
}


# ---- F3: per-channel distributions -------------------------------------------

#' F3: per-channel distribution of unmixed values by method
#'
#' One density per method, overlaid, one panel per channel, on a
#' symmetric-log x axis.
#'
#' @details
#' Each method must come from one file. The x axis uses
#' `scales::transform_pseudo_log(sigma = pseudo_log_sigma)`; the density is
#' estimated on that transformed scale. `pseudo_log_sigma` is a display
#' choice that changes how the negative population looks, so it has no
#' default. Panels share both axes. Non-finite values are dropped and the
#' caption gives the finite events per method per channel.
#'
#' @param events Tibble from [sw_unmix_sample_events()] (columns `file`,
#'   `method` and the channels).
#' @param channels Character vector of channel columns to show.
#' @param pseudo_log_sigma Positive number, the linear-zone width of the
#'   pseudo-log axis. Required.
#'
#' @return A list with `panel` (`"F3"`), `plot`, `source` (long tibble:
#'   `method`, `file`, `event_id`, `channel`, `value`; finite values only)
#'   and `counts` (finite events per `method` and `channel`).
#'
#' @examples
#' set.seed(2)
#' x <- list(
#'   ours = matrix(rnorm(400, sd = 50), ncol = 2, dimnames = list(NULL, c("A", "B"))),
#'   spectroflo = matrix(rnorm(400, sd = 60), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' )
#' ev <- sw_unmix_sample_events(x, n_events = 100, seed = 1)
#' f3 <- sw_unmix_figure_distributions(ev, c("A", "B"), pseudo_log_sigma = 10)
#' f3$counts
#'
#' @export
sw_unmix_figure_distributions <- function(events, channels, pseudo_log_sigma) {
  .sw_fig_require(missing(events), "events")
  .sw_fig_require(missing(channels), "channels")
  .sw_fig_require(missing(pseudo_log_sigma), "pseudo_log_sigma")
  .sw_fig_need_ggplot2()
  trans <- .sw_fig_pseudo_log(pseudo_log_sigma)
  palette <- .sw_fig_check_events(events, channels)
  methods <- names(palette)
  ids <- if ("event_id" %in% names(events)) {
    events[["event_id"]]
  } else {
    seq_len(nrow(events))
  }

  # Long frame, one row per event x channel; non-finite values dropped.
  long <- dplyr::bind_rows(lapply(channels, \(ch) {
    tibble::tibble(
      method = as.character(events[["method"]]),
      file = as.character(events[["file"]]),
      event_id = ids,
      channel = ch,
      value = as.numeric(events[[ch]])
    )
  }))
  source <- long[is.finite(long[["value"]]), ]

  # Finite events per method x channel, for the caption.
  counts <- dplyr::count(source, method, channel, name = "n_events")
  n_by_method <- vapply(
    methods,
    \(m) {
      n <- counts[["n_events"]][counts[["method"]] == m]
      if (min(n) == max(n)) {
        format(min(n), big.mark = ",")
      } else {
        sprintf(
          "%s to %s",
          format(min(n), big.mark = ","),
          format(max(n), big.mark = ",")
        )
      }
    },
    character(1)
  )
  caption <- sprintf(
    "Events per channel: %s.",
    paste(sprintf("%s %s", methods, n_by_method), collapse = "; ")
  )

  plot_df <- source
  plot_df[["method"]] <- factor(plot_df[["method"]], levels = methods)
  plot_df[["channel"]] <- factor(plot_df[["channel"]], levels = channels)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = value, colour = method)
  ) +
    ggplot2::geom_density(linewidth = 0.4, key_glyph = "path") +
    .sw_fig_facet_wrap(channel) +
    .sw_fig_symlog_scale("x", trans) +
    ggplot2::scale_colour_manual(
      values = palette,
      name = "Method",
      guide = .sw_fig_legend(length(palette))
    ) +
    ggplot2::labs(
      x = sprintf(
        "Unmixed intensity (a.u., pseudo-log, sigma = %s)",
        format(pseudo_log_sigma)
      ),
      y = "Density (on the pseudo-log scale)",
      caption = .sw_fig_wrap(caption)
    ) +
    .sw_fig_theme() +
    ggplot2::theme(panel.spacing.x = ggplot2::unit(6, "mm"))

  list(panel = "F3", plot = p, source = source, counts = counts)
}


# ---- F4: biplot --------------------------------------------------------------

#' F4: two-channel biplot, methods side by side
#'
#' A 2-D histogram of `y_channel` against `x_channel`, one panel per method,
#' on shared symmetric-log axes.
#'
#' @details
#' Both axes use `scales::transform_pseudo_log(sigma = pseudo_log_sigma)`.
#' Bins are counted with `ggplot2::geom_bin2d(bins = bins)` on the
#' transformed scale; the fill shows the count on a log10 colour scale.
#' Panels are square because both axes share units. Events with a
#' non-finite value in either channel are dropped; the caption gives the
#' events plotted per method.
#'
#' @param events Tibble from [sw_unmix_sample_events()].
#' @param x_channel,y_channel Channel column names for x and y.
#' @param pseudo_log_sigma Positive number, the linear-zone width of the
#'   pseudo-log axes. Required.
#' @param bins Positive whole number of bins per axis. Required.
#'
#' @return A list with `panel` (`"F4"`), `plot`, `source` (tibble: `method`,
#'   `file`, `event_id`, `x`, `y`; finite pairs only) and `counts` (events per
#'   method).
#'
#' @examples
#' set.seed(3)
#' x <- list(
#'   ours = matrix(rnorm(400, sd = 50), ncol = 2, dimnames = list(NULL, c("A", "B"))),
#'   spectroflo = matrix(rnorm(400, sd = 60), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' )
#' ev <- sw_unmix_sample_events(x, n_events = 100, seed = 1)
#' f4 <- sw_unmix_figure_biplot(ev, "A", "B", pseudo_log_sigma = 10, bins = 30)
#' f4$counts
#'
#' @export
sw_unmix_figure_biplot <- function(
  events,
  x_channel,
  y_channel,
  pseudo_log_sigma,
  bins
) {
  .sw_fig_require(missing(events), "events")
  .sw_fig_require(missing(x_channel), "x_channel")
  .sw_fig_require(missing(y_channel), "y_channel")
  .sw_fig_require(missing(pseudo_log_sigma), "pseudo_log_sigma")
  .sw_fig_require(missing(bins), "bins")
  .sw_fig_need_ggplot2()
  trans <- .sw_fig_pseudo_log(pseudo_log_sigma)
  if (
    !is.numeric(bins) || length(bins) != 1L || bins < 1 || bins != round(bins)
  ) {
    stop("`bins` must be one positive whole number.", call. = FALSE)
  }
  palette <- .sw_fig_check_events(events, c(x_channel, y_channel))
  methods <- names(palette)
  ids <- if ("event_id" %in% names(events)) {
    events[["event_id"]]
  } else {
    seq_len(nrow(events))
  }

  # Plot frame; events with a non-finite value in either channel dropped.
  all_pairs <- tibble::tibble(
    method = as.character(events[["method"]]),
    file = as.character(events[["file"]]),
    event_id = ids,
    x = as.numeric(events[[x_channel]]),
    y = as.numeric(events[[y_channel]])
  )
  source <- all_pairs[
    is.finite(all_pairs[["x"]]) & is.finite(all_pairs[["y"]]),
  ]

  # Events per method, for the caption.
  counts <- dplyr::count(source, method, name = "n_events")
  counts <- counts[match(methods, counts[["method"]]), ]
  caption <- sprintf(
    "Events plotted: %s. %d bins per axis; pseudo-log sigma = %s.",
    paste(
      sprintf(
        "%s %s",
        counts[["method"]],
        format(counts[["n_events"]], big.mark = ",")
      ),
      collapse = "; "
    ),
    as.integer(bins),
    format(pseudo_log_sigma)
  )

  plot_df <- source
  plot_df[["method"]] <- factor(plot_df[["method"]], levels = methods)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bin2d(bins = bins) +
    .sw_fig_facet_wrap(method) +
    .sw_fig_symlog_scale("x", trans) +
    .sw_fig_symlog_scale("y", trans) +
    ggplot2::scale_fill_viridis_c(transform = "log10", name = "Events") +
    ggplot2::labs(
      x = sprintf("%s (a.u., pseudo-log)", x_channel),
      y = sprintf("%s (a.u., pseudo-log)", y_channel),
      caption = .sw_fig_wrap(caption, width = 120) # double column holds ~120 characters at 6 pt
    ) +
    .sw_fig_theme() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "right",
      panel.spacing = ggplot2::unit(4, "mm")
    )

  list(panel = "F4", plot = p, source = source, counts = counts)
}


# ---- Saving ------------------------------------------------------------------

#' Save unmixing figures and their source data
#'
#' Write each figure as PNG (and PDF when cairo is available) and its source
#' table as CSV.
#'
#' @details
#' Files are `<dir>/<stem>.png`, `<dir>/<stem>.pdf` and
#' `<dir>/source_data/<stem>.csv`, where the stem is the panel name
#' (`spectra-overlay`, `similarity-matrix`, `unmixed-distributions`,
#' `biplot`). Default sizes (width x height, mm) are F1 183 x 140, F2
#' 89 x 85, F3 183 x 110 and F4 183 x 110; the widths are Nature column
#' widths. The PNG is 300 dpi and uses the cairo `png()` type when available;
#' the PDF uses [grDevices::cairo_pdf()]. If any target file already exists,
#' nothing is written.
#'
#' With `check = "nature"`, each figure must meet the Nature final artwork
#' spec before anything is written: width 89 or 183 mm, and every text the
#' plot draws between 5 and 7 pt (font sizes read from the built plot on a
#' cairo device). After the PDFs are written, `pdffonts` (poppler) checks
#' that every font is embedded; when `pdffonts` or cairo is absent, that check
#' is skipped and `fonts_embedded` is `NA`. A breach stops with the panel
#' name; files already written stay.
#'
#' @param figures Named list of figure results keyed by panel (`F1` to
#'   `F4`), as returned by the `sw_unmix_figure_*()` builders.
#' @param dir Output folder. Required. Created if absent.
#' @param sizes Optional named list of `c(width, height)` in mm, keyed by
#'   panel, replacing the default size of those panels.
#' @param check `"none"` (default) or `"nature"`, the journal spec check
#'   described in Details.
#'
#' @return A tibble with `panel`, `png`, `pdf` (`NA` without cairo),
#'   `source_csv`, `n_rows` (rows in the source table), `width_mm`,
#'   `height_mm`, and with `check = "nature"` the drawn text range
#'   `text_pt_min`, `text_pt_max` and `fonts_embedded` (`NA` when not checked
#'   or not checkable).
#'
#' @examples
#' \dontrun{
#' s <- rbind(FITC = c(0.1, 1, 0.5), PE = c(0, 0.2, 1))
#' colnames(s) <- c("D1", "D2", "D3")
#' figs <- list(F2 = sw_unmix_figure_similarity(s))
#' sw_unmix_save_figures(figs, dir = tempfile("figs"))
#' }
#'
#' @export
sw_unmix_save_figures <- function(
  figures,
  dir,
  sizes = NULL,
  check = c("none", "nature")
) {
  .sw_fig_require(missing(figures), "figures")
  .sw_fig_require(missing(dir), "dir")
  .sw_fig_need_ggplot2()
  check <- match.arg(check)
  keys <- names(figures)
  if (!is.list(figures) || is.null(keys) || anyDuplicated(keys)) {
    stop("`figures` must be a list with unique panel names.", call. = FALSE)
  }
  unknown <- setdiff(keys, names(.sw_fig_panels))
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown panel(s): %s. Accepted: %s.",
        paste(unknown, collapse = ", "),
        paste(names(.sw_fig_panels), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  # Each result must carry the panel it is filed under.
  carried <- vapply(figures, \(f) as.character(f[["panel"]]), character(1))
  if (!identical(unname(carried), keys)) {
    stop(
      sprintf(
        "Figure names do not match their `panel` fields: %s.",
        paste(sprintf("%s=%s", keys, carried), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Per-panel sizes: defaults, replaced by any panel named in `sizes`.
  bad_size <- setdiff(names(sizes), names(.sw_fig_panels))
  if (
    length(bad_size) > 0L ||
      (length(sizes) > 0L && is.null(names(sizes))) ||
      !all(vapply(
        sizes,
        \(s) is.numeric(s) && length(s) == 2L && all(s > 0),
        logical(1)
      ))
  ) {
    stop(
      "`sizes` must be a list named by panel (F1 to F4) of c(width, height) in mm.",
      call. = FALSE
    )
  }
  sizes <- utils::modifyList(.sw_fig_sizes, as.list(sizes))

  # Journal spec: widths and drawn text sizes, checked before anything is written.
  text_pt <- matrix(NA_real_, nrow = length(keys), ncol = 2)
  if (check == "nature") {
    for (i in seq_along(keys)) {
      text_pt[i, ] <- .sw_fig_check_nature(
        figures[[i]][["plot"]],
        keys[i],
        sizes[[keys[i]]][1]
      )
    }
  }

  use_cairo <- isTRUE(capabilities("cairo"))
  stems <- unname(.sw_fig_panels[keys])
  png_paths <- file.path(dir, paste0(stems, ".png"))
  pdf_paths <- if (use_cairo) {
    file.path(dir, paste0(stems, ".pdf"))
  } else {
    rep(NA_character_, length(stems))
  }
  csv_paths <- file.path(dir, "source_data", paste0(stems, ".csv"))

  # Refuse before writing anything if a target exists.
  targets <- c(png_paths, csv_paths, pdf_paths[!is.na(pdf_paths)])
  exists_now <- targets[file.exists(targets)]
  if (length(exists_now) > 0L) {
    stop(
      sprintf(
        "Refusing to overwrite existing file(s): %s.",
        paste(exists_now, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  dir.create(
    file.path(dir, "source_data"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (i in seq_along(keys)) {
    p <- figures[[i]][["plot"]]
    size <- sizes[[keys[i]]]
    if (use_cairo) {
      ggplot2::ggsave(
        png_paths[i],
        plot = p,
        width = size[1],
        height = size[2],
        units = "mm",
        dpi = 300,
        device = "png",
        type = "cairo",
        bg = "white"
      )
      ggplot2::ggsave(
        pdf_paths[i],
        plot = p,
        width = size[1],
        height = size[2],
        units = "mm",
        device = grDevices::cairo_pdf
      )
    } else {
      ggplot2::ggsave(
        png_paths[i],
        plot = p,
        width = size[1],
        height = size[2],
        units = "mm",
        dpi = 300,
        device = "png",
        bg = "white"
      )
    }
    utils::write.csv(figures[[i]][["source"]], csv_paths[i], row.names = FALSE)
  }

  # Font embedding, read from the written PDFs (needs cairo and poppler's pdffonts).
  fonts_embedded <- rep(NA, length(keys))
  if (check == "nature") {
    for (i in which(!is.na(pdf_paths))) {
      fonts_embedded[i] <- .sw_fig_fonts_embedded(pdf_paths[i])
      if (isFALSE(fonts_embedded[i])) {
        stop(
          sprintf(
            "Panel %s: a font is not embedded in %s.",
            keys[i],
            pdf_paths[i]
          ),
          call. = FALSE
        )
      }
    }
  }

  tibble::tibble(
    panel = keys,
    png = png_paths,
    pdf = pdf_paths,
    source_csv = csv_paths,
    n_rows = vapply(
      figures,
      \(f) nrow(f[["source"]]),
      integer(1),
      USE.NAMES = FALSE
    ),
    width_mm = vapply(keys, \(k) sizes[[k]][1], numeric(1), USE.NAMES = FALSE),
    height_mm = vapply(keys, \(k) sizes[[k]][2], numeric(1), USE.NAMES = FALSE),
    text_pt_min = text_pt[, 1],
    text_pt_max = text_pt[, 2],
    fonts_embedded = fonts_embedded
  )
}

# Font sizes (pt) of every non-empty text grob the built plot draws, measured on a cairo
# device (as the saved PDF) so the font family resolves through fontconfig.
.sw_fig_drawn_font_sizes <- function(p) {
  path <- tempfile(fileext = ".pdf")
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(path)
  } else {
    grDevices::pdf(path)
  }
  on.exit({
    grDevices::dev.off()
    unlink(path)
  })
  g <- ggplot2::ggplotGrob(p)
  sizes <- numeric(0)
  walk <- function(x) {
    label <- x[["label"]]
    fontsize <- x[["gp"]][["fontsize"]]
    if (
      !is.null(label) && any(nzchar(as.character(label))) && !is.null(fontsize)
    ) {
      cex <- x[["gp"]][["cex"]]
      sizes <<- c(sizes, fontsize * if (is.null(cex)) 1 else cex)
    }
    if (inherits(x, "gtable")) {
      lapply(x[["grobs"]], walk)
    } else if (inherits(x, "gTree")) {
      lapply(x[["children"]], walk)
    }
    invisible(NULL)
  }
  walk(g)
  sizes
}

# Stop unless the figure meets the Nature width and text-size spec; returns c(min, max) pt.
.sw_fig_check_nature <- function(p, panel, width_mm) {
  if (!width_mm %in% .sw_fig_nature_widths) {
    stop(
      sprintf(
        "Panel %s: width %s mm; Nature needs %s mm.",
        panel,
        format(width_mm),
        paste(.sw_fig_nature_widths, collapse = " or ")
      ),
      call. = FALSE
    )
  }
  pt <- round(.sw_fig_drawn_font_sizes(p), 2)
  if (
    length(pt) > 0L &&
      (min(pt) < .sw_fig_nature_text_pt[["min"]] ||
        max(pt) > .sw_fig_nature_text_pt[["max"]])
  ) {
    stop(
      sprintf(
        "Panel %s: drawn text %.1f to %.1f pt; Nature needs %s to %s pt.",
        panel,
        min(pt),
        max(pt),
        .sw_fig_nature_text_pt[["min"]],
        .sw_fig_nature_text_pt[["max"]]
      ),
      call. = FALSE
    )
  }
  if (length(pt) == 0L) c(NA_real_, NA_real_) else range(pt)
}

# Run poppler's pdffonts; NULL when it is not installed (wrapped so tests can mock it).
.sw_pdffonts <- function(path) {
  if (!nzchar(Sys.which("pdffonts"))) {
    return(NULL)
  }
  system2("pdffonts", shQuote(path), stdout = TRUE, stderr = FALSE)
}

# TRUE when every font in the PDF is embedded, NA when pdffonts is unavailable.
.sw_fig_fonts_embedded <- function(path) {
  out <- .sw_pdffonts(path)
  if (is.null(out)) {
    return(NA)
  }
  rows <- out[-(1:2)]
  if (length(rows) == 0L) {
    return(TRUE)
  }
  # pdffonts rows end in "emb sub uni object ID": emb is the fifth field from the end
  emb <- vapply(
    strsplit(trimws(rows), "\\s+"),
    \(x) x[length(x) - 4L],
    character(1)
  )
  all(emb == "yes")
}

# Bare column names used inside ggplot2::aes() and dplyr::count().
utils::globalVariables(c(
  "detector",
  "value_norm",
  "method",
  "panel_label",
  "fluorophore_x",
  "fluorophore_y",
  "cosine",
  "value",
  "channel",
  "x",
  "y"
))
