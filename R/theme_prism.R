# Theme code supplied by the package author (Mikhael D. Manurung), who states they hold its
# copyright, included under MIT (user decision 2026-11-11). The code matches `theme_prism()` in
# ggprism (github.com/csdaw/ggprism).
#
# Changes from the supplied code, so it runs without ggprism (plan v2 step 2):
# * the `palette` argument is gone; the `black_and_white` colours are hard-coded;
# * the custom element `prism.ticks.length` is dropped (only ggprism registers it);
# * `legend.text.align` and `legend.title.align` are dropped (removed from ggplot2 >= 3.5);
# * `deframe()` and `is_bool()` are replaced by base R;
# * the parent is `theme_grey(base_size, base_family)` instead of ggprism's `all_null` list.

# ggprism `black_and_white` palette colours (read from ggprism 1.0.7 `ggprism_data`).
.sw_prism_colours <- c(
  pageBackgroundColor = "#FFFFFF",
  plottingAreaColor = "#FFFFFF",
  axisColor = "#000000",
  axisLabelColor = "#000000",
  axisTitleColor = "#000000",
  graphTitleColor = "#000000"
)

#' Prism-style ggplot2 theme (internal)
#'
#' A black-and-white theme with bold axis lines, outward ticks, no grid and
#' no strip boxes.
#'
#' @param base_size Base font size, in pt.
#' @param base_family Base font family.
#' @param base_fontface Base font face.
#' @param base_line_size Base linewidth for line elements.
#' @param base_rect_size Base linewidth for rect elements.
#' @param axis_text_angle Angle of x axis text: one of 0, 45, 90, 270.
#' @param border Draw a border around the panel instead of axis lines.
#'
#' @return A complete ggplot2 theme.
#' @keywords internal
#' @noRd
.sw_theme_prism <- function(
  base_size,
  base_family = "sans",
  base_fontface = "bold",
  base_line_size = base_size / 14,
  base_rect_size = base_size / 14,
  axis_text_angle = 0,
  border = FALSE
) {
  .sw_fig_require(missing(base_size), "base_size")
  .sw_fig_need_ggplot2()

  # Ensure x axis text is at a sensible angle
  angle <- axis_text_angle[1]
  if (!angle %in% c(0, 45, 90, 270)) {
    stop(
      sprintf(
        "'axis_text_angle' must be one of [%s]",
        paste(c(0, 45, 90, 270), collapse = ", ")
      ),
      ".\nFor other angles, use the guide_axis() function in ggplot2 instead",
      call. = FALSE
    )
  }
  colours <- .sw_prism_colours

  # Draw border or not
  if (!(is.logical(border) && length(border) == 1L && !is.na(border))) {
    stop("border must be either: TRUE or FALSE", call. = FALSE)
  }
  if (border) {
    panel.border <- ggplot2::element_rect(fill = NA)
    axis.line <- ggplot2::element_blank()
  } else {
    panel.border <- ggplot2::element_blank()
    axis.line <- ggplot2::element_line()
  }

  margin <- ggplot2::margin
  unit <- ggplot2::unit
  rel <- ggplot2::rel
  element_text <- ggplot2::element_text
  element_blank <- ggplot2::element_blank
  element_rect <- ggplot2::element_rect

  t <- ggplot2::theme(
    # Base elements (to be inherited by other elements)
    line = ggplot2::element_line(
      colour = colours[["axisColor"]],
      linewidth = base_line_size,
      linetype = 1,
      lineend = "square"
    ),
    rect = element_rect(
      fill = "white",
      colour = colours[["axisColor"]],
      linewidth = base_rect_size,
      linetype = 1
    ),
    text = element_text(
      family = base_family,
      face = base_fontface,
      colour = colours[["graphTitleColor"]],
      size = base_size,
      lineheight = 0.9,
      hjust = 0.5,
      vjust = 0.5,
      angle = 0,
      margin = margin(),
      debug = FALSE
    ),

    # Normal ggplot2 theme elements
    axis.line = axis.line,
    axis.line.x = NULL,
    axis.line.y = NULL,
    axis.text = element_text(
      size = rel(0.95),
      colour = colours[["axisLabelColor"]]
    ),
    axis.text.x = element_text(
      margin = margin(t = 0.8 * base_size / 4),
      angle = axis_text_angle,
      hjust = ifelse(axis_text_angle %in% c(45, 90, 270), 1, 0.5),
      vjust = ifelse(axis_text_angle %in% c(0, 90, 270), 0.5, 1)
    ),
    axis.text.x.top = element_text(
      margin = margin(b = 0.8 * base_size / 4),
      vjust = 0
    ),
    axis.text.y = element_text(
      margin = margin(r = 0.5 * base_size / 4),
      hjust = 1
    ),
    axis.text.y.right = element_text(
      margin = margin(l = 0.5 * base_size / 4),
      hjust = 0
    ),
    axis.ticks = ggplot2::element_line(),
    axis.ticks.length = unit(base_size / 2.5, "pt"),
    axis.ticks.length.x = NULL,
    axis.ticks.length.x.top = NULL,
    axis.ticks.length.x.bottom = NULL,
    axis.ticks.length.y = NULL,
    axis.ticks.length.y.left = NULL,
    axis.ticks.length.y.right = NULL,
    axis.title = element_text(colour = colours[["axisTitleColor"]]),
    axis.title.x = element_text(
      margin = margin(t = base_size * 0.6),
      vjust = 1
    ),
    axis.title.x.top = element_text(
      margin = margin(b = base_size * 0.6),
      vjust = 0
    ),
    axis.title.y = element_text(
      angle = 90,
      margin = margin(r = base_size * 0.6),
      vjust = 1
    ),
    axis.title.y.right = element_text(
      angle = -90,
      margin = margin(l = base_size * 0.6),
      vjust = 0
    ),

    legend.background = element_blank(),
    legend.spacing = unit(base_size, "pt"),
    legend.spacing.x = NULL,
    legend.spacing.y = NULL,
    legend.margin = margin(
      base_size / 2,
      base_size / 2,
      base_size / 2,
      base_size / 2
    ),
    legend.key = element_blank(),
    legend.key.size = unit(1.2, "lines"),
    legend.key.height = NULL,
    legend.key.width = unit(base_size * 1.8, "pt"),
    legend.text = element_text(size = rel(0.8), face = "plain"),
    legend.title = element_blank(),
    legend.position = "right",
    legend.direction = NULL,
    legend.justification = "center",
    legend.box = NULL,
    legend.box.margin = margin(0, 0, 0, 0, "cm"),
    legend.box.background = element_blank(),
    legend.box.spacing = unit(base_size, "pt"),

    panel.background = element_rect(fill = NA, colour = NA),
    panel.border = panel.border,
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(base_size / 2, "pt"),
    panel.spacing.x = NULL,
    panel.spacing.y = NULL,
    panel.ontop = FALSE,

    strip.background = element_blank(),
    strip.text = element_text(
      colour = colours[["axisTitleColor"]],
      size = rel(0.8),
      margin = margin(
        base_size / 2.5,
        base_size / 2.5,
        base_size / 2.5,
        base_size / 2.5
      )
    ),
    strip.text.x = element_text(margin = margin(b = base_size / 3)),
    strip.text.y = element_text(
      angle = -90,
      margin = margin(l = base_size / 3)
    ),
    strip.text.y.left = element_text(angle = 90),
    strip.placement = "inside",
    strip.placement.x = NULL,
    strip.placement.y = NULL,
    strip.switch.pad.grid = unit(base_size / 4, "pt"),
    strip.switch.pad.wrap = unit(base_size / 4, "pt"),

    plot.background = element_rect(
      fill = colours[["pageBackgroundColor"]],
      colour = NA
    ),
    plot.title = element_text(
      size = rel(1.2),
      hjust = 0.5,
      vjust = 1,
      margin = margin(b = base_size)
    ),
    plot.title.position = "panel",
    plot.subtitle = element_text(
      hjust = 0.5,
      vjust = 1,
      margin = margin(b = base_size / 2)
    ),
    plot.caption = element_text(
      size = rel(0.8),
      hjust = 1,
      vjust = 1,
      margin = margin(t = base_size / 2)
    ),
    plot.caption.position = "panel",
    plot.tag = element_text(size = rel(1.2), hjust = 0.5, vjust = 0.5),
    plot.tag.position = "topleft",
    plot.margin = margin(
      base_size / 2,
      base_size / 2,
      base_size / 2,
      base_size / 2
    ),

    complete = TRUE
  )

  # Elements the code above leaves unset come from theme_grey (ggprism uses its all_null list)
  parent <- ggplot2::theme_grey(
    base_size = base_size,
    base_family = base_family
  )
  ggplot2::`%+replace%`(parent, t)
}

#' Nature-spec figure theme (internal)
#'
#' The prism-style theme at 7 pt with plain text, bold plot title, strips,
#' axis titles and legend title, line width 0.4 and the legend on the right,
#' so every drawn text sits in 5 to 7 pt (Nature final artwork guide; user
#' rule r-0017, 2026-11-11).
#'
#' @param font_family Font family. Display only.
#'
#' @return A complete ggplot2 theme.
#' @keywords internal
#' @noRd
.sw_theme_nature <- function(font_family = "sans") {
  .sw_theme_prism(
    base_size = 7,
    base_family = font_family,
    base_fontface = "plain",
    base_line_size = 0.4
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 7),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      strip.clip = "off"
    )
}
