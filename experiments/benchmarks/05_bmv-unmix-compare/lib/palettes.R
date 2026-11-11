# Project palettes and figure choices the house style leaves open.
# Fill each field once per project; every figure script reads this file.
#   dotplot defaults:    [fill in: size range, colour scale, dot.min]
#   volcano thresholds:  [fill in; one seen pair is Q < 0.1 and |log2FC| > 0.5]
#   heatmap anchors:     [fill in; one seen example is -4 / 0 / 4]
#   font family:         [fill in: e.g. Helvetica]
#
# Palettes are named vectors keyed by biological label, never by position.
# The two below are placeholders; the project replaces them with its own labels.
pal_condition <- c(healthy = "#4477AA", disease = "#EE6677")
pal_cell_type <- c(t_cell = "#228833", b_cell = "#CCBB44", monocyte = "#AA3377")

# Colour-blind-safe categorical sets (fig.colourblind). Hex values from the
# sources: Paul Tol (tol-colors 2.2.0 colors.json), Okabe-Ito
# (grDevices::palette.colors), Tableau "Color Blind" (ggthemes); standard
# Tableau 10 is not safe (red against green) and is not listed.
TOL_BRIGHT <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
OKABE_ITO <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7",
  "#000000"
)
TOL_MUTED <- c(
  "#CC6677",
  "#332288",
  "#DDCC77",
  "#117733",
  "#88CCEE",
  "#882255",
  "#44AA99",
  "#999933",
  "#AA4499"
)
TABLEAU_CB10 <- c(
  "#1170AA",
  "#FC7D0B",
  "#A3ACB9",
  "#57606C",
  "#5FA2CE",
  "#C85200",
  "#7B848F",
  "#A3CCE9",
  "#FFBC79",
  "#C8D0D9"
)
# Tol discrete rainbow: his 29-colour ramp and, for each n (11-23), the
# indexes of the colours to use.
TOL_RAINBOW_HEX <- c(
  "#E8ECFB",
  "#D9CCE3",
  "#D1BBD7",
  "#CAACCB",
  "#BA8DB4",
  "#AE76A3",
  "#AA6F9E",
  "#994F88",
  "#882E72",
  "#1965B0",
  "#437DBF",
  "#5289C7",
  "#6195CF",
  "#7BAFDE",
  "#4EB265",
  "#90C987",
  "#CAE0AB",
  "#F7F056",
  "#F7CB45",
  "#F6C141",
  "#F4A736",
  "#F1932D",
  "#EE8026",
  "#E8601C",
  "#E65518",
  "#DC050C",
  "#A5170E",
  "#72190E",
  "#42150A"
)
TOL_RAINBOW_INDEX <- list(
  `11` = c(9, 10, 12, 14, 15, 17, 18, 21, 24, 26, 28),
  `12` = c(3, 6, 9, 10, 12, 14, 15, 17, 18, 21, 24, 26),
  `13` = c(3, 6, 9, 10, 12, 14, 15, 16, 17, 18, 21, 24, 26),
  `14` = c(3, 6, 9, 10, 12, 14, 15, 16, 17, 18, 20, 22, 24, 26),
  `15` = c(3, 6, 9, 10, 12, 14, 15, 16, 17, 18, 20, 22, 24, 26, 28),
  `16` = c(3, 5, 7, 9, 10, 12, 14, 15, 16, 17, 18, 20, 22, 24, 26, 28),
  `17` = c(3, 5, 7, 8, 9, 10, 12, 14, 15, 16, 17, 18, 20, 22, 24, 26, 28),
  `18` = c(3, 5, 7, 8, 9, 10, 12, 14, 15, 16, 17, 18, 20, 22, 24, 26, 27, 28),
  `19` = c(2, 4, 5, 7, 8, 9, 10, 12, 14, 15, 16, 17, 18, 20, 22, 24, 26, 27, 28),
  `20` = c(2, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 20, 22, 24, 26, 27, 28),
  `21` = c(2, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 21, 23, 25, 26, 27, 28),
  `22` = c(2, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 21, 23, 25, 26, 27, 28, 29),
  `23` = c(1, 2, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 21, 23, 25, 26, 27, 28, 29)
)

# A colour-blind-safe palette keyed by label, from the smallest set that fits:
# up to 7 Tol bright, 8 Okabe-Ito, 9 Tol muted, 10 Tableau Color Blind, 11-23
# Tol discrete rainbow. Past 23 no set stays distinguishable, so it stops.
categorical_palette <- function(labels) {
  labels <- unique(as.character(labels))
  n <- length(labels)
  for (colours in list(TOL_BRIGHT, OKABE_ITO, TOL_MUTED, TABLEAU_CB10)) {
    if (n <= length(colours)) {
      return(setNames(colours[seq_len(n)], labels))
    }
  }
  if (n > 23) {
    stop(
      sprintf(
        "categorical_palette: %d categories; no colour set stays colour-blind safe past 23. ",
        n
      ),
      "Facet, label directly or lump rare categories.",
      call. = FALSE
    )
  }
  setNames(TOL_RAINBOW_HEX[TOL_RAINBOW_INDEX[[as.character(n)]]], labels)
}

# Black or white text for a hex fill, whichever contrasts more (WCAG luminance).
text_colour <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  lin <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  lum <- colSums(lin * c(0.2126, 0.7152, 0.0722))
  ifelse((lum + 0.05) / 0.05 >= 1.05 / (lum + 0.05), "#000000", "#FFFFFF")
}

# Colours for a vector of labels. An unknown label stops the script, so a
# mislabelled or renamed group can never fall through as an unplanned colour.
palette_for <- function(labels, pal) {
  # Character first: a factor would index pal by level number, not by name.
  labels <- as.character(labels)
  missing <- setdiff(unique(labels), names(pal))
  if (length(missing) > 0) {
    stop(
      "palette_for: labels not in the palette: ",
      toString(missing),
      "; add them to palettes.R or fix the labels.",
      call. = FALSE
    )
  }
  pal[labels]
}

if (sys.nframe() == 0) {
  stopifnot(identical(unname(palette_for("disease", pal_condition)), "#EE6677"))
  stopifnot(inherits(try(palette_for("unknown", pal_condition), silent = TRUE), "try-error"))
  stopifnot(identical(unname(categorical_palette(c("a", "b"))), TOL_BRIGHT[1:2]))
  stopifnot(length(unique(categorical_palette(1:23))) == 23)
  stopifnot(identical(text_colour(c("#000000", "#F0E442")), c("#FFFFFF", "#000000")))
}
