# Purpose: figures for 92_umap_knn.R. (1) per-donor median kNN overlap of each arm's UMAP with the
#          reference UMAP (1.8.8 zip or SpectroFlo), next to the 1.8.8 zip seed noise floor; (2) UMAP per
#          arm coloured by TNF and Th2, one figure per marker set and scaling. Lucid-coder style (lib/).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/92b_plot_umap_knn.R \
#            experiments/runs/2026-11-10_bmv-umap-knn   (repo root; writes into that run dir)
# Inputs:  <run_dir>/knn_overlap.csv, embeddings.rds  [92_umap_knn.R]
# Outputs: <run_dir>/fig_knn_overlap.{pdf,png,csv}            overlap per donor, arms and noise floor
#          <run_dir>/fig_umap_<marker set>_<scaling>.{pdf,png,csv}  UMAPs coloured by TNF and Th2
#          <run_dir>/sessionInfo_figures.txt

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
})
LIB <- Sys.getenv("SW_LIB", "experiments/benchmarks/05_bmv-unmix-compare/lib")
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))
source(file.path(LIB, "plots.R"))

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
UMAP_SEED <- 42 # as 92
ARM_ORDER <- c("SpectroFlo", "1.8.8 zip", "opt188", "refined AF", "1.0.2 zip")
NOISE <- paste0("1.8.8 zip, seed ", 1:3)
COLOUR_LIMITS <- c(-0.5, 3) # arcsinh display window for TNF and Th2, shared by every panel; display only
MARKER_SET_LABEL <- c(
  lineage = "14 lineage markers",
  lineage_functional = "lineage + 17 functional"
)
SCALING_LABEL <- c(zscore = "z-scored", none = "unscaled")

stopifnot(
  "92 outputs present" = file.exists(file.path(RUN_DIR, "knn_overlap.csv"))
)

# Figure 1: neighbour overlap per donor -----------------------------------------------------------------
overlap_df <- read_csv(
  file.path(RUN_DIR, "knn_overlap.csv"),
  show_col_types = FALSE
) |>
  mutate(
    group = if_else(comparison %in% NOISE, "1.8.8 zip, other seed", comparison),
    group = factor(group, levels = c(ARM_ORDER, "1.8.8 zip, other seed")),
    reference = factor(
      paste("vs", reference),
      levels = c("vs 1.8.8 zip", "vs SpectroFlo")
    ),
    panel = factor(
      paste0(MARKER_SET_LABEL[marker_set], ", ", SCALING_LABEL[scaling]),
      levels = as.vector(outer(
        MARKER_SET_LABEL,
        SCALING_LABEL,
        paste,
        sep = ", "
      ))
    )
  ) |>
  filter(!is.na(group))
pal_group <- categorical_palette(levels(overlap_df$group))
p_overlap <- ggplot(
  overlap_df,
  aes(group, median_overlap, colour = group, shape = donor)
) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.6,
    colour = "black",
    linewidth = 0.3,
    aes(group = group),
    show.legend = FALSE
  ) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0, seed = 1),
    size = 1.4
  ) +
  scale_colour_manual(values = pal_group, guide = "none") +
  facet_grid(reference ~ panel) +
  labs(
    x = NULL,
    y = "Median share of 30 UMAP neighbours kept",
    shape = "Donor",
    title = "CD3+ T cells: UMAP neighbourhoods vs reference"
  ) +
  rotate_x()
save_panel(
  p_overlap,
  overlap_df,
  RUN_DIR,
  "fig_knn_overlap",
  width_mm = 183,
  height_mm = 120
)

# Figure 2: UMAP per arm coloured by TNF and Th2 ----------------------------------------------------------
embeddings <- readRDS(file.path(RUN_DIR, "embeddings.rds")) |>
  filter(seed == UMAP_SEED)
umap_panel <- \(df, marker, title) {
  df <- df[order(df[[marker]]), ]
  ggplot(df, aes(dim1, dim2, colour = .data[[marker]])) +
    embedding_points() +
    scale_colour_viridis_c(limits = COLOUR_LIMITS, oob = scales::squish) +
    coord_equal() +
    corner_arrows(df$dim1, df$dim2, "UMAP") +
    labs(colour = sprintf("%s\n(arcsinh)", marker), title = title) +
    no_axes()
}
umap_figure <- \(ms, sc) {
  df <- embeddings |>
    filter(marker_set == ms, scaling == sc) |>
    select(arm, umap) |>
    tidyr::unnest(umap) |>
    mutate(arm = factor(arm, levels = ARM_ORDER))
  # one row per marker: five arms without legends, then that marker's shared colour bar
  marker_row <- \(marker) {
    panels <- map(ARM_ORDER, \(a) {
      umap_panel(filter(df, arm == a), marker, a) +
        theme(legend.position = "none")
    })
    legend <- cowplot::get_legend(
      umap_panel(
        filter(df, arm == ARM_ORDER[1]),
        marker,
        ""
      ) +
        theme(legend.position = "right", legend.direction = "vertical")
    )
    cowplot::plot_grid(
      plotlist = c(panels, list(legend)),
      nrow = 1,
      rel_widths = c(rep(1, 5), 0.6)
    )
  }
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      sprintf(
        "CD3+ T cells, UMAP on %s, %s",
        MARKER_SET_LABEL[[ms]],
        SCALING_LABEL[[sc]]
      ),
      size = 9,
      x = 0,
      hjust = 0
    )
  p <- cowplot::plot_grid(
    title,
    marker_row("TNF"),
    marker_row("Th2"),
    ncol = 1,
    rel_heights = c(0.1, 1, 1)
  )
  # cowplot draws on a transparent canvas; a white one keeps the PNG readable
  p <- cowplot::ggdraw(p) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  save_panel(
    p,
    df,
    RUN_DIR,
    sprintf("fig_umap_%s_%s", ms, sc),
    width_mm = 183,
    height_mm = 90
  )
}
umap_figure("lineage", "zscore")
umap_figure("lineage", "none")
umap_figure("lineage_functional", "zscore")
umap_figure("lineage_functional", "none")

writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_figures.txt")
)
