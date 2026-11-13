# Purpose: copy of 92b_plot_umap_knn.R for the run of 92_umap_knn_common.R (common gates, plans/common_gates_plan.md).
#          Figures for it. (1) per-donor median kNN overlap of each arm's UMAP with the
#          reference UMAP (1.8.8 zip or SpectroFlo), next to the 1.8.8 zip seed noise floor; (2) UMAP per
#          arm coloured by the config colour markers (3D: TNF and Th2), one figure per marker set and scaling.
#          With conditions, the overlap figure has one row block per condition and shapes mark the donor id.
#          Lucid-coder style (lib/).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> \
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/92b_plot_umap_knn_common.R [<run dir>]
#            (repo root; writes into that run dir; without the argument it is run_in("common-umap-knn"))
# Inputs:  datasets/<SW_DATASET>.json (cfg("umap_markers")), <SW_DATASET>_runs.csv
#          <run_dir>/knn_overlap.csv, embeddings.rds  [92_umap_knn_common.R]; run_dir is registry key common-umap-knn
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
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(RUN_DIR)) {
  RUN_DIR <- run_in("common-umap-knn")
}
UMAP_SEED <- 42 # as 92
ARM_ORDER <- c("SpectroFlo", "1.8.8 zip", "opt188", "refined AF", "1.0.2 zip")
NOISE <- paste0("1.8.8 zip, seed ", 1:3)
COLOUR_LIMITS <- c(-0.5, 3) # arcsinh display window for the colour markers, shared by every panel; display only
UMAP_MARKERS <- cfg("umap_markers") # as 92
COLOUR_MARKERS <- UMAP_MARKERS$colour
MARKER_SET_LABEL <- c(
  lineage = sprintf("%d lineage markers", length(UMAP_MARKERS$lineage)),
  lineage_functional = sprintf(
    "lineage + %d functional",
    length(UMAP_MARKERS$functional)
  )
)
HEIGHT_PER_ROW_MM <- 60 # overlap figure: one row per reference (and condition); 2 rows give the 120 mm of 3D
# shapes for up to 10 donor ids when there are conditions (3D keeps the ggplot2 default shapes)
DONOR_SHAPES <- c(16, 17, 15, 18, 8, 3, 4, 7, 9, 10)
# replicate column of knn_overlap.csv: the sample in 3D and EXV, the donor id when conditions split the samples
SHAPE_COL <- if (is.null(DS$conditions)) "donor" else "donor_id"
FACET_ROWS <- c("reference", if (!is.null(DS$conditions)) "condition")
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
n_shapes <- n_distinct(overlap_df[[SHAPE_COL]])
stopifnot("a shape for every donor id" = n_shapes <= length(DONOR_SHAPES))
p_overlap <- ggplot(
  overlap_df,
  aes(group, median_overlap, colour = group, shape = .data[[SHAPE_COL]])
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
  facet_grid(rows = vars(!!!rlang::syms(FACET_ROWS)), cols = vars(panel)) +
  labs(
    x = NULL,
    y = "Median share of 30 UMAP neighbours kept",
    shape = "Donor",
    title = "CD3+ T cells: UMAP neighbourhoods vs reference"
  ) +
  rotate_x()
if (!is.null(DS$conditions)) {
  p_overlap <- p_overlap +
    scale_shape_manual(values = DONOR_SHAPES[seq_len(n_shapes)])
}
save_panel(
  p_overlap,
  overlap_df,
  RUN_DIR,
  "fig_knn_overlap",
  width_mm = 183,
  height_mm = HEIGHT_PER_ROW_MM * n_distinct(overlap_df[FACET_ROWS])
)

# Figure 2: UMAP per arm coloured by the colour markers ----------------------------------------------------------
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
      size = 7,
      x = 0,
      hjust = 0
    )
  p <- cowplot::plot_grid(
    plotlist = c(list(title), map(COLOUR_MARKERS, marker_row)),
    ncol = 1,
    rel_heights = c(0.1, rep(1, length(COLOUR_MARKERS)))
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
    height_mm = 45 * length(COLOUR_MARKERS) # 2 colour markers give the 90 mm of 3D
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
