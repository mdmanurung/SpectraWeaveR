# Purpose: redraw the report's section 19 figures (as188_joint vs opt188 downstream; UMAP neighbourhoods)
#          in the Nature figure style (user, 2026-11-11; lucid-coder r-0008, r-0012, r-0014 to r-0017) from
#          the tables and embeddings their runs saved. Nothing is recomputed; values are the originals.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_downstream.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-downstream-opt188-v2/fig_pop_pct.csv    [89b_pop_pct_fixed_cut.R]
#          experiments/runs/2026-11-10_bmv-downstream-opt188/fig_cluster_abundance.csv   [89_downstream_as188_vs_opt188.R]
#          experiments/runs/2026-11-10_bmv-downstream-opt188/fig_marker_flags.csv        [89]
#          experiments/runs/2026-11-10_bmv-downstream-opt188/fig_ari.csv                 [89]
#          experiments/runs/2026-11-10_bmv-umap-knn/knn_overlap.csv     neighbour overlap per donor [92_umap_knn.R]
#          experiments/runs/2026-11-10_bmv-umap-knn/embeddings.rds      UMAP per arm, seed 42 [92_umap_knn.R]
# Outputs: experiments/runs/2026-11-11_bmv-restyle/19_downstream_umap/fig_pop_pct.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/19_downstream_umap/fig_cluster_abundance.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/19_downstream_umap/fig_marker_flags.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/19_downstream_umap/fig_ari.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/19_downstream_umap/fig_knn_overlap.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/19_downstream_umap/fig_umap_lineage_zscore.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/19_downstream_umap/sessionInfo_95_downstream.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
})
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))
source(file.path(LIB, "plots.R"))

# Parameters ---------------------------------------------------------------------------------------
RUNS <- "experiments/runs"
OUT <- "experiments/runs/2026-11-11_bmv-restyle/19_downstream_umap"
ARMS <- c("as188_joint", "opt188") # 89, grill Q1
FLAG_LIMIT <- 0.25 # 89: symmetric limit of the diverging squeezed-share scale (display only)
DONORS <- c("A3", "A6", "B11", "B2", "D1")
DONOR_SHAPES <- setNames(c(16, 17, 15, 18, 8), DONORS) # as the approved demo (2026-11-11)
UMAP_SEED <- 42 # 92
ARM_ORDER <- c("SpectroFlo", "1.8.8 zip", "opt188", "refined AF", "1.0.2 zip") # 92b
NOISE <- paste0("1.8.8 zip, seed ", 1:3) # 92b: noise-floor comparisons
COLOUR_LIMITS <- c(-0.5, 3) # 92b: shared arcsinh display window for TNF and Th2 (display only)
MARKER_SET_LABEL <- c(
  lineage = "14 lineage markers",
  lineage_functional = "lineage + 17 functional"
)
SCALING_LABEL <- c(zscore = "z-scored", none = "unscaled")

stopifnot("F4 section 19 dir is new" = !dir.exists(OUT))
dir.create(OUT, recursive = TRUE)
pal_arm <- categorical_palette(c(ARMS, "as188_joint, other seed")) # 89
shape_arm <- c(as188_joint = 16, opt188 = 17, "as188_joint, other seed" = 1) # 89: encode twice
# "as188\n(n = 5)": arm short name and donor count, as 89
arm_n_labels <- \(donor) {
  n <- n_distinct(donor)
  \(arm) sprintf("%s\n(n = %d)", sub("_joint", "", arm), n)
}

# Paired arm comparison per facet: donors joined by lines, median crossbar ---------------------------
# Used by the population and metacluster figures (the two figures differ only in data and labels).
paired_arms <- function(d, value, facet, y_label, title) {
  ggplot(d, aes(arm, .data[[value]])) +
    geom_line(aes(group = donor), colour = "grey60", linewidth = FIG_LINE / 2) +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.5,
      linewidth = FIG_LINE / 2,
      colour = "black"
    ) +
    geom_point(aes(colour = arm, shape = arm), size = 1) +
    scale_colour_manual(values = pal_arm, guide = "none") +
    scale_shape_manual(values = shape_arm, guide = "none") +
    scale_x_discrete(labels = arm_n_labels(d$donor)) +
    facet_wrap(
      vars(.data[[facet]]),
      scales = "free_y",
      ncol = 4,
      axes = "all",
      axis.labels = "margins"
    ) +
    labs(x = NULL, y = y_label, title = title) +
    value_grid("y")
}

# Figure: gated populations, % of parent -----------------------------------------------------------------
pop_df <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-downstream-opt188-v2/fig_pop_pct.csv"),
  show_col_types = FALSE
) |>
  mutate(
    arm = factor(arm, levels = ARMS),
    population = factor(population, levels = unique(population)) # gate-tree order, as 89b
  )
p_pop <- paired_arms(
  pop_df,
  "pct_of_parent",
  "population",
  "% of parent",
  "Gated populations, % of parent"
)
save_panel(p_pop, pop_df, OUT, "fig_pop_pct", width_mm = 183, height_mm = 120)

# Figure: metacluster abundance, matched clusters -------------------------------------------------------
abund_df <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-downstream-opt188/fig_cluster_abundance.csv"),
  show_col_types = FALSE
) |>
  mutate(
    arm = factor(arm, levels = ARMS),
    cluster = factor(cluster, levels = sort(unique(cluster)))
  )
p_abund <- paired_arms(
  abund_df,
  "pct",
  "cluster",
  "% of live events",
  "Metacluster abundance, matched clusters"
)
save_panel(
  p_abund,
  abund_df,
  OUT,
  "fig_cluster_abundance",
  width_mm = 183,
  height_mm = 200
)

# Figure: marker flags, donor-median squeezed-share change (opt188 - as188_joint) ------------------------
# Rows keep 89's marker order (first seen); NA (CD3 has no negatives in T cells) drawn grey.
flag_df <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-downstream-opt188/fig_marker_flags.csv"),
  col_types = cols(text = col_character()) # keep "0.00", not 0
) |>
  mutate(
    label = factor(label, levels = unique(label)),
    text = if_else(is.na(median_d_squeezed), "NA", text)
  )
p_flag <- plot_heatmap(
  flag_df,
  row = "label",
  col = "pop",
  value = "median_d_squeezed",
  limit = FLAG_LIMIT,
  title = "Squeezed share change",
  x_label = NULL,
  y_label = NULL,
  fill_label = "opt188 minus\nas188, donor\nmedian"
) +
  geom_text(aes(label = text), size = 5.5, size.unit = "pt") +
  guides(
    fill = guide_colourbar(barheight = unit(25, "mm"), barwidth = unit(3, "mm"))
  )
save_panel(
  p_flag,
  flag_df,
  OUT,
  "fig_marker_flags",
  width_mm = 89,
  height_mm = 130
)

# Figure: ARI per donor, opt188 and the as188_joint seed noise floor ---------------------------------------
ari_df <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-downstream-opt188/fig_ari.csv"),
  show_col_types = FALSE
) |>
  mutate(
    comparison = factor(
      comparison,
      levels = c("opt188", "as188_joint, other seed")
    )
  )
p_ari <- ggplot(ari_df, aes(comparison, ari)) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.5,
    linewidth = FIG_LINE / 2,
    colour = "black"
  ) +
  geom_point(
    aes(colour = comparison, shape = comparison),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    size = 1.2
  ) +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_shape_manual(values = shape_arm, guide = "none") +
  scale_x_discrete(labels = n_labels(ari_df$comparison)) +
  labs(
    x = NULL,
    y = "ARI vs as188_joint (seed 42)",
    title = "Cluster agreement per donor"
  ) +
  value_grid("y")
save_panel(p_ari, ari_df, OUT, "fig_ari", width_mm = 89, height_mm = 60)

# Figure: UMAP neighbour overlap per donor, flipped so arm names read horizontally (r-0012) -------------
overlap_df <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-umap-knn/knn_overlap.csv"),
  show_col_types = FALSE
) |>
  mutate(
    group = if_else(comparison %in% NOISE, "1.8.8 zip, other seed", comparison),
    group = factor(group, levels = rev(c(ARM_ORDER, "1.8.8 zip, other seed"))),
    reference = factor(
      paste("vs", reference),
      levels = c("vs 1.8.8 zip", "vs SpectroFlo")
    ),
    panel = factor(
      paste0(MARKER_SET_LABEL[marker_set], ",\n", SCALING_LABEL[scaling]),
      levels = as.vector(outer(
        MARKER_SET_LABEL,
        SCALING_LABEL,
        paste,
        sep = ",\n"
      ))
    )
  )
pal_group <- categorical_palette(rev(levels(overlap_df$group)))
p_overlap <- ggplot(overlap_df, aes(median_overlap, group)) +
  stat_summary(
    aes(group = group),
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.6,
    colour = "black",
    linewidth = FIG_LINE / 2,
    orientation = "y"
  ) +
  geom_point(
    aes(colour = group, shape = donor),
    position = position_jitter(width = 0, height = 0.15, seed = 1),
    size = 1
  ) +
  scale_colour_manual(values = pal_group, guide = "none") +
  scale_shape_manual(values = DONOR_SHAPES, name = "Donor") +
  # free_y + space: each reference row lists only the arms scored against it
  facet_grid(
    reference ~ panel,
    scales = "free_y",
    space = "free_y",
    axes = "all",
    axis.labels = "margins"
  ) +
  labs(
    x = "Median share of 30 UMAP neighbours kept",
    y = NULL,
    title = "CD3+ T cells: UMAP neighbourhoods vs reference"
  ) +
  value_grid("x")
save_panel(
  p_overlap,
  overlap_df,
  OUT,
  "fig_knn_overlap",
  width_mm = 183,
  height_mm = 85
)

# Figure: UMAP per arm coloured by TNF and Th2, lineage markers z-scored -----------------------------------
# Each arm has its own embedding, so each panel is its own plot (equal aspect, no frame, corner arrows);
# TNF and Th2 share one viridis window, so one colour bar serves both rows.
umap_df <- readRDS(file.path(RUNS, "2026-11-10_bmv-umap-knn/embeddings.rds")) |>
  filter(seed == UMAP_SEED, marker_set == "lineage", scaling == "zscore") |>
  select(arm, umap) |>
  tidyr::unnest(umap) |>
  mutate(arm = factor(arm, levels = ARM_ORDER))
umap_long <- umap_df |>
  tidyr::pivot_longer(c(TNF, Th2), names_to = "marker", values_to = "arcsinh")
umap_panel <- function(d, title, show_arrows) {
  d <- d[order(d$arcsinh), ] # high values on top
  arrow_len <- 0.22 * c(diff(range(d$dim1)), diff(range(d$dim2)))
  origin <- c(min(d$dim1), min(d$dim2))
  p <- ggplot(d, aes(dim1, dim2, colour = arcsinh)) +
    scattermore::geom_scattermore(pointsize = 1.2) +
    scale_colour_viridis_c(
      limits = COLOUR_LIMITS,
      oob = scales::squish,
      name = "arcsinh\n(cofactor\n6000)"
    ) +
    coord_equal(clip = "off") +
    labs(title = title) +
    theme_void(base_size = 7, base_family = FIG_FONT) +
    theme(
      plot.title = element_text(face = "bold", size = 7, hjust = 0.5),
      legend.title = element_text(face = "bold", size = 6),
      legend.text = element_text(size = 6),
      legend.position = "none"
    )
  if (!show_arrows) {
    return(p)
  }
  head <- arrow(length = unit(0.8, "mm"), type = "closed")
  p +
    annotate(
      "segment",
      x = origin[1],
      y = origin[2],
      xend = origin[1] + arrow_len[1],
      yend = origin[2],
      arrow = head,
      linewidth = FIG_LINE / 2
    ) +
    annotate(
      "segment",
      x = origin[1],
      y = origin[2],
      xend = origin[1],
      yend = origin[2] + arrow_len[2],
      arrow = head,
      linewidth = FIG_LINE / 2
    ) +
    annotate(
      "text",
      x = origin[1],
      y = origin[2] - 0.06 * diff(range(d$dim2)),
      label = "UMAP 1",
      hjust = 0,
      size = 5,
      size.unit = "pt",
      family = FIG_FONT
    ) +
    annotate(
      "text",
      x = origin[1] - 0.06 * diff(range(d$dim1)),
      y = origin[2],
      label = "UMAP 2",
      angle = 90,
      hjust = 0,
      size = 5,
      size.unit = "pt",
      family = FIG_FONT
    )
}
# one row per marker: five arms, arrows on the first panel only
marker_row <- function(which_marker) {
  panels <- imap(ARM_ORDER, \(a, i) {
    umap_panel(
      filter(umap_long, .data$arm == a, .data$marker == which_marker),
      "",
      show_arrows = i == 1
    )
  })
  label <- cowplot::ggdraw() +
    cowplot::draw_label(
      which_marker,
      fontface = "bold",
      size = 7,
      fontfamily = FIG_FONT,
      angle = 90
    )
  cowplot::plot_grid(
    plotlist = c(list(label), panels),
    nrow = 1,
    rel_widths = c(0.15, rep(1, 5))
  )
}
legend <- cowplot::get_legend(
  umap_panel(filter(umap_long, .data$arm == ARM_ORDER[1]), "", FALSE) +
    theme(legend.position = "right", legend.key.height = unit(6, "mm"))
)
title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "CD3+ T cells, UMAP on 14 lineage markers, z-scored",
    fontface = "bold",
    size = 7,
    fontfamily = FIG_FONT
  )
# arm names in their own row, so they align whatever each embedding's aspect
arm_header <- cowplot::plot_grid(
  plotlist = c(
    list(cowplot::ggdraw()),
    map(ARM_ORDER, \(a) {
      cowplot::ggdraw() +
        cowplot::draw_label(
          a,
          fontface = "bold",
          size = 7,
          fontfamily = FIG_FONT
        )
    })
  ),
  nrow = 1,
  rel_widths = c(0.15, rep(1, 5))
)
grid_rows <- cowplot::plot_grid(
  arm_header,
  marker_row("TNF"),
  marker_row("Th2"),
  ncol = 1,
  rel_heights = c(0.1, 1, 1)
)
p_umap <- cowplot::plot_grid(
  title,
  cowplot::plot_grid(grid_rows, legend, nrow = 1, rel_widths = c(1, 0.12)),
  ncol = 1,
  rel_heights = c(0.07, 1)
)
# cowplot draws on a transparent canvas; a white one keeps the PNG readable
p_umap <- cowplot::ggdraw(p_umap) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
save_panel(
  p_umap,
  umap_long,
  OUT,
  "fig_umap_lineage_zscore",
  width_mm = 183,
  height_mm = 85
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT, "sessionInfo_95_downstream.txt")
)
