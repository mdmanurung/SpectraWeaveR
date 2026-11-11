# Figure theme and panel export, shared by every figure script.
# Nature spec (user, 2026-11-11; lucid-coder r-0016, r-0017): width 89 mm (single column) or
# 183 mm (double), all drawn text 5-7 pt, Helvetica (Nimbus Sans, metric clone) embedded.
library(ggplot2)

# Journal spec -------------------------------------------------------------------------------
JOURNAL_WIDTHS_MM <- c(single = 89, double = 183) # Nature final artwork
JOURNAL_TEXT_PT <- c(min = 5, max = 7) # Nature: 5-7 pt
FIG_FONT <- "Nimbus Sans" # Helvetica clone (URW), in /usr/share/fonts
FIG_LINE <- 0.4 # axis line width at 7 pt (user-approved demo, 2026-11-11)

# Theme ----------------------------------------------------------------------------------------
# theme_prism with plain text and bold titles, strips and axis titles (r-0017); legend right (r-0015).
theme_set(
  ggprism::theme_prism(
    base_size = 7,
    base_family = FIG_FONT,
    base_fontface = "plain",
    base_line_size = FIG_LINE
  ) +
    theme(
      plot.title = element_text(face = "bold", size = 7),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      strip.clip = "off"
    )
)

# Thin grid along the value axis only, for plots that compare values along that axis (r-0014).
value_grid <- function(axis = "x") {
  major <- element_line(colour = "grey85", linewidth = FIG_LINE / 2)
  minor <- element_line(colour = "grey92", linewidth = FIG_LINE / 3)
  if (axis == "x") {
    theme(panel.grid.major.x = major, panel.grid.minor.x = minor)
  } else {
    theme(panel.grid.major.y = major, panel.grid.minor.y = minor)
  }
}

# Font sizes (pt) of every text grob the built plot draws.
drawn_font_sizes <- function(p) {
  # a cairo device, as for the saved PDF, so Nimbus Sans resolves through fontconfig
  grDevices::cairo_pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off())
  g <- if (inherits(p, "ggplot")) ggplotGrob(p) else p
  sizes <- c()
  walk <- function(x) {
    if (!is.null(x$label) && any(nzchar(x$label)) && !is.null(x$gp$fontsize)) {
      sizes <<- c(sizes, x$gp$fontsize * (x$gp$cex %||% 1))
    }
    if (inherits(x, "gtable")) {
      for (ch in x$grobs) walk(ch)
    } else if (inherits(x, "gTree")) {
      for (ch in x$children) walk(ch)
    }
  }
  suppressWarnings(walk(g))
  sizes
}

# Write one panel as <run_dir>/<name>.pdf (vector, for the figure), .png (for review) and .csv
# (the source data, the same frame that was passed to ggplot()). Stops if the panel breaks the
# journal spec: width not a column width, drawn text outside 5-7 pt, or a font not embedded.
save_panel <- function(p, data, run_dir, name, width_mm, height_mm) {
  if (!dir.exists(run_dir)) {
    stop("save_panel: run_dir ", run_dir, " does not exist; create it first.", call. = FALSE)
  }
  if (!width_mm %in% JOURNAL_WIDTHS_MM) {
    stop(
      sprintf("save_panel %s: width %s mm; use 89 (single) or 183 (double).", name, width_mm),
      call. = FALSE
    )
  }
  sizes <- round(drawn_font_sizes(p), 2)
  if (min(sizes) < JOURNAL_TEXT_PT[["min"]] || max(sizes) > JOURNAL_TEXT_PT[["max"]]) {
    stop(
      sprintf("save_panel %s: drawn text %.1f-%.1f pt; Nature needs 5-7 pt.", name, min(sizes), max(sizes)),
      call. = FALSE
    )
  }
  stem <- file.path(run_dir, name)
  # No dpi on the PDF: vector output ignores it.
  ggsave(
    paste0(stem, ".pdf"),
    p,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = cairo_pdf
  )
  fonts <- system2("pdffonts", paste0(stem, ".pdf"), stdout = TRUE)[-(1:2)]
  # pdffonts rows end in "emb sub uni object ID": emb is the fifth field from the end
  embedded <- vapply(strsplit(trimws(fonts), "\\s+"), \(x) x[length(x) - 4], character(1))
  if (length(fonts) == 0 || any(embedded != "yes")) {
    stop(sprintf("save_panel %s: a font is not embedded in the PDF.", name), call. = FALSE)
  }
  ggsave(
    paste0(stem, ".png"),
    p,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 300,
    bg = "white"
  )
  write.csv(data, paste0(stem, ".csv"), row.names = FALSE)
  invisible(stem)
}
