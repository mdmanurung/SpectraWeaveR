# Purpose: redraw the section 19 figures (opt188 vs as188_joint downstream; UMAP neighbourhoods) from the
#          common-gate runs (95_restyle_downstream.R repointed; plans/common_gates_plan.md G11, G13). Donors
#          are encoded by colour and shape as in the report; no line joins a donor across arms, and the
#          median is a bar drawn behind the points. Nothing is recomputed.
# Run:     SW_DATASET=3d SW_RUN_DATE=YYYY-MM-DD pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_downstream.R
#          (repo root)
# Inputs:  run_in("common-downstream-opt188-v2")/fig_pop_pct.csv           [89b common]
#          run_in("common-downstream-opt188")/fig_{cluster_abundance,marker_flags,ari}.csv [89 common]
#          run_in("common-umap-knn")/knn_overlap.csv, embeddings.rds        [92 common]
# Outputs: run_out("restyle-common") (shared by the 95c scripts)/19_downstream_umap/
#            fig_{pop_pct,cluster_abundance,marker_flags,ari,knn_overlap,umap_lineage_zscore}.{pdf,png,csv},
#            sessionInfo_95c_downstream.txt
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")
SCRIPT <- "95c_restyle_downstream.R"
TNF_LAB <- cfg("targets")$tnf_label
TH2_LAB <- need_target(SCRIPT, "th2_label") # the UMAP figure colours by TNF and Th2

# Parameters ---------------------------------------------------------------------------------------
ARMS <- c("as188_joint", "opt188") # 89, grill Q1
FLAG_LIMIT <- 0.25 # 89: symmetric limit of the diverging squeezed-share scale (display only)
UMAP_SEED <- 42 # 92
ARM_ORDER <- c("SpectroFlo", "1.8.8 zip", "opt188", "refined AF", "1.0.2 zip") # 92b
NOISE <- paste0("1.8.8 zip, seed ", 1:3) # 92b: noise-floor comparisons
COLOUR_LIMITS <- c(-0.5, 3) # 92b: shared arcsinh display window for TNF and Th2 (display only)
MARKER_SET_LABEL <- c(
  lineage = "14 lineage markers",
  lineage_functional = "lineage + 17 functional"
)
SCALING_LABEL <- c(zscore = "z-scored", none = "unscaled")
RUN_89 <- run_in("common-downstream-opt188")
RUN_92 <- run_in("common-umap-knn")

OUT <- fresh_dir("19_downstream_umap")
OLD <- "2026-11-11_bmv-restyle/19_downstream_umap"
log_fig <- function(name, source) {
  add_manifest(
    SCRIPT,
    file.path(OLD, paste0(name, ".png")),
    file.path(OUT, paste0(name, ".png")),
    source
  )
}
# "as188\n(n = 5)": arm short name and donor count, as 89
arm_n_labels <- \(donor) {
  n <- n_distinct(donor)
  \(arm) sprintf("%s\n(n = %d)", sub("_joint", "", arm), n)
}

# Arms along x, donors as colour and shape over a median bar, no line across arms --------------------
paired_arms <- function(d, value, facet, y_label, title) {
  need_cond(d, SCRIPT)
  d <- mutate(d, donor = factor_by_pal(.data$donor, DONOR_PAL))
  ggplot(d, aes(arm, .data[[value]])) +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.6,
      linewidth = FIG_LINE / 2,
      colour = "grey30"
    ) +
    geom_point(aes(colour = donor, shape = donor), size = 1) +
    donor_scales() +
    scale_x_discrete(labels = arm_n_labels(d$donor)) +
    facet_wrap(
      vars(.data[[facet]], !!!rlang::syms(COND)),
      scales = "free_y",
      ncol = 4,
      axes = "all",
      axis.labels = "margins"
    ) +
    labs(x = NULL, y = y_label, title = title) +
    value_grid("y")
}

# Figure: gated populations, % of parent -----------------------------------------------------------------
src_pop <- file.path(run_in("common-downstream-opt188-v2"), "fig_pop_pct.csv")
pop_df <- read_csv(src_pop, show_col_types = FALSE) |>
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
log_fig("fig_pop_pct", src_pop)

# Figure: metacluster abundance, matched clusters -------------------------------------------------------
src_abund <- file.path(RUN_89, "fig_cluster_abundance.csv")
abund_df <- read_csv(src_abund, show_col_types = FALSE) |>
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
log_fig("fig_cluster_abundance", src_abund)

# Figure: marker flags, donor-median squeezed-share change (opt188 - as188_joint) ------------------------
src_flag <- file.path(RUN_89, "fig_marker_flags.csv")
flag_df <- read_csv(src_flag, col_types = cols(text = col_character())) |>
  mutate(
    # PE has no marker in the control file (marker NA in fcs_control_file.csv)
    label = if_else(label == "NA (PE)", "PE (no marker)", label),
    label = factor(label, levels = unique(label)),
    pop = factor(pop, levels = unique(pop))
  ) |>
  complete_grid("label", "pop") |>
  mutate(text = if_else(is.na(median_d_squeezed), "", text))
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
  heat_fill(FLAG_LIMIT, "opt188 minus\nas188, donor\nmedian") +
  labs(caption = NOT_TESTED) +
  theme(plot.caption = element_text(hjust = 0)) +
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
log_fig("fig_marker_flags", src_flag)

# Figure: ARI per donor, opt188 and the as188_joint seed noise floor ---------------------------------------
src_ari <- file.path(RUN_89, "fig_ari.csv")
ari_df <- read_csv(src_ari, show_col_types = FALSE) |>
  need_cond(SCRIPT) |>
  mutate(
    comparison = factor(
      comparison,
      levels = c("opt188", "as188_joint, other seed")
    ),
    donor = factor_by_pal(donor, DONOR_PAL)
  )
p_ari <- ggplot(ari_df, aes(comparison, ari)) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.6,
    linewidth = FIG_LINE / 2,
    colour = "grey30"
  ) +
  geom_point(
    aes(colour = donor, shape = donor),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    size = 1.2
  ) +
  donor_scales() +
  scale_x_discrete(labels = scales::label_wrap(12)) +
  cond_facet() +
  labs(
    x = NULL,
    y = "ARI vs as188_joint (seed 42)",
    title = "Cluster agreement per donor"
  ) +
  value_grid("y")
save_panel(p_ari, ari_df, OUT, "fig_ari", width_mm = 89, height_mm = 60)
log_fig("fig_ari", src_ari)

# Figure: UMAP neighbour overlap per donor, flipped so arm names read horizontally (r-0012) -------------
src_knn <- file.path(RUN_92, "knn_overlap.csv")
overlap_df <- read_csv(src_knn, show_col_types = FALSE) |>
  need_cond(SCRIPT) |>
  mutate(
    group = if_else(comparison %in% NOISE, "1.8.8 zip, other seed", comparison),
    group = factor(group, levels = rev(c(ARM_ORDER, "1.8.8 zip, other seed"))),
    reference = factor(
      paste("vs", reference),
      levels = c("vs 1.8.8 zip", "vs SpectroFlo")
    ),
    donor = factor_by_pal(donor, DONOR_PAL),
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
p_overlap <- ggplot(overlap_df, aes(median_overlap, group)) +
  stat_summary(
    aes(group = group),
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.6,
    colour = "grey30",
    linewidth = FIG_LINE / 2,
    orientation = "y"
  ) +
  geom_point(
    aes(colour = donor, shape = donor),
    position = position_jitter(width = 0, height = 0.15, seed = 1),
    size = 1
  ) +
  donor_scales() +
  # free_y + space: each reference row lists only the arms scored against it
  facet_grid(
    rows = vars(reference, !!!rlang::syms(COND)),
    cols = vars(panel),
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
log_fig("fig_knn_overlap", src_knn)

# Figure: UMAP per arm coloured by TNF and Th2, lineage markers z-scored -----------------------------------
# Each arm has its own embedding, so each panel is its own plot (equal aspect, no frame, corner arrows);
# TNF and Th2 share one viridis window, so one colour bar serves both rows.
umap_df <- readRDS(file.path(RUN_92, "embeddings.rds")) |>
  filter(seed == UMAP_SEED, marker_set == "lineage", scaling == "zscore") |>
  select(arm, umap) |>
  tidyr::unnest(umap) |>
  mutate(arm = factor(arm, levels = ARM_ORDER))
umap_long <- umap_df |>
  tidyr::pivot_longer(
    all_of(c(TNF_LAB, TH2_LAB)),
    names_to = "marker",
    values_to = "arcsinh"
  )
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
  marker_row(TNF_LAB),
  marker_row(TH2_LAB),
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
log_fig("fig_umap_lineage_zscore", file.path(RUN_92, "embeddings.rds"))

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT, "sessionInfo_95c_downstream.txt")
)
