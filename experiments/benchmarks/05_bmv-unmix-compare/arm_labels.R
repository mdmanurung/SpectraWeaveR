# Purpose: display labels and the version-row page layout for the unmixing arms, shared by the
#          figure scripts (10, 19, 31-34, 38, 40). Arm IDs stay the folder and table keys; only plot
#          titles and tick labels use these labels.
# Run:     source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") from the repo root. It does not
#          read SW_DATASET or SW_RUN_DATE: arm IDs, labels and version rows are the same for 3D, 24H and EXV
#          (each zip ships the same four arms; the self-unmixed arms use the same IDs).
# Inputs:  none. Settings per arm come from the zip scripts and unmix_cells/<arm>/args.txt.
# Outputs: ARM_LABELS, a named character vector (arm ID -> label)
#          ZIP_ARMS, ZIP_FILL, the arms shipped in the zip and their panel tint
#          VERSION_ROWS, arms per code-version row
#          PANEL_THEME, text sizes for one panel in a version grid
#          version_grid(), lays out one ggcyto panel per arm in labelled version rows (cowplot)

# Labels ----------------------------------------------------------------------------------
# One short line per arm, used as panel titles and boxplot tick labels. "zip" = your original run.
ARM_LABELS <- c(
  spectroflo = "SpectroFlo (zip)",
  as102_basic = "AS 1.0.2, no variants (zip)",
  as102_var = "AS 1.0.2, variants (zip)",
  best102_neg = "AS 1.0.2 best, neg. spectra",
  best102_af = "AS 1.0.2 best, AF-clean spectra",
  as188_joint = "AS 1.8.8 joint (zip, ref.)",
  as188_fast_csv = "AS 1.8.8 joint, rerun",
  best188_p1 = "AS 1.8.8 best, 1 pass",
  best188_p5 = "AS 1.8.8 best, 5 passes",
  s188_c102 = "1.8.8 spectra, 1.0.2 code",
  s102_c188 = "1.0.2 spectra, 1.8.8 code"
)

# Zip arms ------------------------------------------------------------------------------
# The four arms the user shipped in each dataset zip (3D, 24H, EXV); their panels get a tinted background
# (user, 2026-11-10) so they stand out from the reruns.
ZIP_ARMS <- c("spectroflo", "as102_basic", "as102_var", "as188_joint")
ZIP_FILL <- "#FFE9B3" # light amber; display only

# Version rows --------------------------------------------------------------------------
# Figure rows: one row per code version (user, 2026-11-10); crossed arms are named by the code they ran
VERSION_ROWS <- list(
  `SpectroFlo` = "spectroflo",
  `AutoSpectral 1.0.2` = c(
    "as102_basic",
    "as102_var",
    "best102_neg",
    "best102_af"
  ),
  `AutoSpectral 1.8.8` = c(
    "as188_joint",
    "as188_fast_csv",
    "best188_p1",
    "best188_p5"
  ),
  `Crossed: spectra x code` = c("s188_c102", "s102_c188")
)

# Panel theme -----------------------------------------------------------------------------
# Text sizes for one ggcyto panel inside version_grid(), legible at the 14 x 10 inch page.
PANEL_THEME <- ggplot2::theme(
  strip.text = ggplot2::element_blank(),
  legend.position = "none",
  plot.title = ggplot2::element_text(size = 10, face = "bold"),
  plot.subtitle = ggplot2::element_text(size = 9),
  axis.title = ggplot2::element_text(size = 9),
  axis.text = ggplot2::element_text(size = 8)
)

# Version grid: one labelled row per code version, panels named by arm -------------------------
version_grid <- function(panels) {
  n_col <- max(lengths(VERSION_ROWS))
  VERSION_ROWS |>
    imap(\(row_arms, label) {
      cells <- c(
        map(row_arms, \(a) {
          g <- as.ggplot(panels[[a]])
          if (a %in% ZIP_ARMS) {
            g <- g +
              theme(
                plot.background = element_rect(fill = ZIP_FILL, colour = NA)
              )
          }
          g
        }),
        rep(list(NULL), n_col - length(row_arms))
      )
      c(
        list(
          cowplot::ggdraw() +
            cowplot::draw_label(
              stringr::str_wrap(label, width = 12),
              angle = 90,
              fontface = "bold",
              size = 12
            )
        ),
        cells
      )
    }) |>
    list_flatten() |>
    cowplot::plot_grid(
      plotlist = _,
      ncol = n_col + 1,
      rel_widths = c(0.3, rep(1, n_col))
    )
}
