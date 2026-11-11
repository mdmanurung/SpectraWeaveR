# Purpose: redraw the ranking, cytokine-benchmark and G1/G2 figures from their saved source tables,
#          with no caption text on the figure and legible fonts (user, 2026-11-10). Plots only: no
#          value is recomputed; arm keys and methods live in the summary INDEX.md, not on the figure.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/41_replot_group_b.R
#          from the repo root.
# Inputs:  experiments/runs/2026-11-10_bmv-ranking3/fig_{spread,ssm,ssm_heatmap,tilt}_source.csv [30_rank_metrics2.R]
#          experiments/runs/2026-11-10_bmv-ranking-sym/fig_<metric>_sym_source.csv [34_rank_symmetric.R]
#          experiments/runs/2026-11-10_bmv-cytokine/fig_<metric>_source.csv [31_cytokine_benchmark.R]
#          experiments/runs/2026-11-10_bmv-g1g2-figures/fig_g1_heldout_source.csv,
#            fig_g2_simulation_source.csv, fig_g2_simulation_seeds.csv [35_plot_g1_g2.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  VERSION_ROWS (arm order and groups)
# Outputs: experiments/runs/2026-11-10_bmv-figures-v2/<source run>/<same file name>.png, one per input
#          figure: ranking3 (3), ranking-sym (6), cytokine (6), g1g2-figures (2)
#          experiments/runs/2026-11-10_bmv-figures-v2/sessionInfo_41.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # VERSION_ROWS: arm order, groups

# Parameters ------------------------------------------------------------------------------
BASE_SIZE <- 11 # font size in pt at the saved size (user, 2026-11-10: legible labels)
DPI <- 200
JITTER_SEED <- 42 # jitter only moves points sideways, as in 35_plot_g1_g2.R
DONOR_COLOURS <- c(
  A3 = "#0072B2",
  A6 = "#D55E00",
  B11 = "#009E73",
  B2 = "#CC79A7",
  D1 = "#E69F00"
) # Okabe-Ito, one colour per donor (31_cytokine_benchmark.R)

RUNS <- "experiments/runs"
OUT <- file.path(RUNS, "2026-11-10_bmv-figures-v2")
SOURCE_RUNS <- c(
  "2026-11-10_bmv-ranking3",
  "2026-11-10_bmv-ranking-sym",
  "2026-11-10_bmv-cytokine",
  "2026-11-10_bmv-g1g2-figures"
)
out_dirs <- set_names(file.path(OUT, SOURCE_RUNS), SOURCE_RUNS)
stopifnot("41 outputs must not exist" = !any(dir.exists(out_dirs)))
walk(out_dirs, \(d) dir.create(d, recursive = TRUE))

theme_set(
  theme_bw(base_size = BASE_SIZE) +
    theme(
      axis.text = element_text(size = BASE_SIZE - 2),
      strip.text = element_text(size = BASE_SIZE - 1),
      legend.text = element_text(size = BASE_SIZE - 2),
      panel.grid.minor = element_blank()
    )
)
arm_order <- unlist(VERSION_ROWS, use.names = FALSE)
group_end <- cumsum(lengths(VERSION_ROWS))
group_breaks <- head(group_end, -1) + 0.5 # dashed lines between code-version groups
group_mid <- group_end - lengths(VERSION_ROWS) / 2 + 0.5
group_names <- c("SpectroFlo", "AS 1.0.2", "AS 1.8.8", "Crossed")
x_arm <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Ranking3: donors as points over a median bar, arms ordered by the median -------------------
ranking3 <- function(file, ylab) {
  d <- read_csv(
    file.path(RUNS, SOURCE_RUNS[[1]], file),
    show_col_types = FALSE
  ) |>
    mutate(arm = reorder(arm, value, FUN = median))
  med <- summarise(d, value = median(value), .by = arm)
  ggplot(d, aes(arm, value)) +
    geom_col(data = med, fill = "grey80", width = 0.7) +
    geom_point(
      aes(colour = sample),
      position = position_jitter(width = 0.15, height = 0, seed = 1),
      size = 2
    ) +
    scale_colour_manual(values = DONOR_COLOURS) +
    labs(x = NULL, y = ylab, colour = "Donor") +
    x_arm
}
p_spread <- ranking3("fig_spread_source.csv", "Median Nguyen rSD of negatives")
ggsave(
  file.path(out_dirs[[1]], "fig_spread.png"),
  p_spread,
  width = 180,
  height = 110,
  units = "mm",
  dpi = DPI
)
p_tilt <- ranking3("fig_tilt_source.csv", "Median |Spearman rho|")
ggsave(
  file.path(out_dirs[[1]], "fig_tilt.png"),
  p_tilt,
  width = 180,
  height = 110,
  units = "mm",
  dpi = DPI
)
# SSM: arm summary (a) and the reference arm's spillover-spread heatmap (b)
p_ssm_a <- ranking3("fig_ssm_source.csv", "Median spillover spread")
p_ssm_b <- read_csv(
  file.path(RUNS, SOURCE_RUNS[[1]], "fig_ssm_heatmap_source.csv"),
  show_col_types = FALSE
) |>
  ggplot(aes(target, source, fill = ss)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "sqrt", na.value = "white", name = "SS") +
  labs(x = "Target channel", y = "Source channel", title = "as188_joint") +
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7)
  )
p_ssm <- cowplot::plot_grid(
  p_ssm_a,
  p_ssm_b,
  ncol = 2,
  rel_widths = c(1, 1.25),
  labels = c("a", "b")
)
ggsave(
  file.path(out_dirs[[1]], "fig_ssm.png"),
  p_ssm,
  width = 330,
  height = 160,
  units = "mm",
  dpi = DPI,
  bg = "white"
)

# Ranking-sym: boxplot per arm, donors joined by lines, arms in version-group order -------------
sym_ylabs <- c(
  spread = "Median Nguyen rSD of negatives",
  stain_index = "Median stain index",
  hyperneg = "Median hypernegatives (%)",
  streak = "TNF vs CD8 Spearman rho",
  tilt = "Median |Spearman rho|",
  ssm = "Median spillover spread"
)
iwalk(sym_ylabs, \(ylab, metric) {
  d <- read_csv(
    file.path(RUNS, SOURCE_RUNS[[2]], sprintf("fig_%s_sym_source.csv", metric)),
    show_col_types = FALSE
  ) |>
    mutate(arm = factor(arm, levels = arm_order))
  p <- ggplot(d, aes(arm, value)) +
    geom_vline(
      xintercept = group_breaks,
      linetype = "dashed",
      colour = "grey60"
    ) +
    geom_boxplot(
      outlier.shape = NA,
      fill = "grey92",
      colour = "grey40",
      width = 0.55
    ) +
    geom_line(
      aes(group = sample, colour = sample),
      linewidth = 0.4,
      alpha = 0.7
    ) +
    geom_point(aes(colour = sample), size = 2) +
    annotate(
      "text",
      x = group_mid,
      y = Inf,
      label = group_names,
      vjust = 1.4,
      size = 3,
      fontface = "bold"
    ) +
    scale_colour_manual(values = DONOR_COLOURS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = ylab, colour = "Donor") +
    x_arm
  ggsave(
    file.path(out_dirs[[2]], sprintf("fig_%s_sym.png", metric)),
    p,
    width = 215,
    height = 115,
    units = "mm",
    dpi = DPI,
    bg = "white"
  )
})

# Cytokine benchmark: median crossbar, donors as points, open = flagged -----------------------
cyto_ylabs <- c(
  pct_pos_own = "% positive (own gate)",
  stain_index = "Stain index",
  neg_rsd = "Nguyen rSD of negatives",
  pileup = "Share of negatives near 0",
  emd = "EMD vs as188_joint (arcsinh)",
  biology = "% or ratio"
)
iwalk(cyto_ylabs, \(ylab, metric) {
  d <- read_csv(
    file.path(RUNS, SOURCE_RUNS[[3]], sprintf("fig_%s_source.csv", metric)),
    show_col_types = FALSE
  ) |>
    mutate(arm = factor(arm, levels = arm_order))
  facets <- if (n_distinct(d$row) == 1) {
    facet_wrap(vars(col), scales = "free_y")
  } else {
    facet_grid(row ~ col, scales = "free_y")
  }
  p <- ggplot(d, aes(arm, value)) +
    geom_vline(
      xintercept = group_breaks,
      colour = "grey70",
      linetype = "dashed"
    ) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.7,
      linewidth = 0.3,
      colour = "grey20"
    ) +
    geom_point(
      aes(colour = donor, shape = open),
      position = position_jitter(width = 0.15, height = 0, seed = 1),
      size = 2,
      stroke = 0.7,
      na.rm = TRUE
    ) +
    facets +
    scale_colour_manual(values = DONOR_COLOURS, name = "Donor") +
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 1),
      labels = c(`FALSE` = "no", `TRUE` = "yes"),
      name = "Flagged"
    ) +
    labs(x = NULL, y = ylab) +
    x_arm
  ggsave(
    file.path(out_dirs[[3]], sprintf("fig_%s.png", metric)),
    p,
    width = 250,
    height = if (n_distinct(d$row) == 1) 120 else 200,
    units = "mm",
    dpi = DPI,
    bg = "white"
  )
})

# G1 held-out controls: log2(arm / as188_joint) per control ----------------------------------
g1_arms <- c(
  as102_basic = "as102_basic",
  as102_var = "as102_var",
  s188_c102 = "s188_c102",
  s102_c188 = "s102_c188"
)
g1_groups <- c(
  `cell single-stain` = "Cell single stains",
  `bead single-stain` = "Bead single stains",
  `Unstained (Cells)` = "Unstained cells",
  `Negative_1 (Beads)` = "Unstained beads",
  `H12 Unstained_BMV_Aurora_3dCD328_Batch1` = "H12 unstained\n(not held out)"
) # facet labels as in 35_plot_g1_g2.R, shortened
g1 <- read_csv(
  file.path(RUNS, SOURCE_RUNS[[4]], "fig_g1_heldout_source.csv"),
  show_col_types = FALSE
) |>
  mutate(
    arm = factor(unname(g1_arms[arm]), levels = rev(g1_arms)),
    group = factor(unname(g1_groups[group]), levels = g1_groups),
    metric = factor(
      metric,
      levels = c("Off-target spread", "Off-target |bias|")
    )
  )
stopifnot(
  "every G1 arm and group has a label" = !anyNA(g1$arm) && !anyNA(g1$group)
)
p_g1 <- ggplot(g1, aes(log2_ratio, arm)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_boxplot(
    data = filter(g1, n_controls >= 5),
    outlier.shape = NA,
    width = 0.5,
    fill = NA
  ) +
  geom_point(
    position = position_jitter(height = 0.15, width = 0, seed = JITTER_SEED),
    alpha = 0.7,
    size = 1.8
  ) +
  facet_grid(metric ~ group, scales = "free_x") +
  labs(x = "log2(arm / as188_joint) per control", y = NULL)
ggsave(
  file.path(out_dirs[[4]], "fig_g1_heldout.png"),
  p_g1,
  width = 330,
  height = 150,
  units = "mm",
  dpi = DPI,
  bg = "white"
)

# G2 simulation: per-dye error ratio to the 1.8.8 joint recipe --------------------------------
g2_arms <- c(
  ols = "OLS",
  as102_basic = "1.0.2 code, no variants",
  as102_var = "1.0.2 code, variants"
) # all G2 arms use the 1.8.8 spectra and AF (35_plot_g1_g2.R)
g2_levels <- c(
  "RMSE (all events)",
  "|bias| (all events)",
  "|bias| (dim layer)",
  "Negative spread"
)
g2_read <- function(file) {
  read_csv(file.path(RUNS, SOURCE_RUNS[[4]], file), show_col_types = FALSE) |>
    mutate(
      arm = factor(unname(g2_arms[arm]), levels = rev(g2_arms)),
      metric = factor(metric, levels = g2_levels)
    )
}
g2_dyes <- g2_read("fig_g2_simulation_source.csv")
g2_seeds <- g2_read("fig_g2_simulation_seeds.csv")
p_g2 <- ggplot(g2_dyes, aes(ratio, arm)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(
    data = g2_seeds,
    colour = "grey60",
    alpha = 0.3,
    size = 0.9,
    position = position_jitter(height = 0.2, width = 0, seed = JITTER_SEED)
  ) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA) +
  geom_point(
    position = position_jitter(height = 0.2, width = 0, seed = JITTER_SEED),
    size = 1.8,
    alpha = 0.8
  ) +
  scale_x_log10() +
  facet_wrap(vars(metric), nrow = 1, scales = "free_x") +
  labs(x = "Error ratio to 1.8.8 joint (log scale)", y = NULL)
ggsave(
  file.path(out_dirs[[4]], "fig_g2_simulation.png"),
  p_g2,
  width = 330,
  height = 110,
  units = "mm",
  dpi = DPI,
  bg = "white"
)

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_41.txt"))
