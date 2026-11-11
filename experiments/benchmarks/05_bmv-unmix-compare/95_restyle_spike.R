# Purpose: redraw the spike-in summary figures of report sections 18, 20, 22 and 23 at Nature spec
#          (user, 2026-11-11: theme_prism 7 pt Nimbus Sans, 89 or 183 mm, facet axes on every panel,
#          value-axis grid, long labels flipped, legend right). Plot only: reads the tables the spike
#          runs saved; nothing is unmixed or rescored. Fixes the 88b/88e ECDF bug, where the linetype
#          scale held 2 values keyed to the first 2 of 5 setup labels, so "af405_afrefine" was not
#          drawn.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_spike.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-spike-in-v2/spike_summary.csv, spike_cells.csv
#            TNF, Th2, IL-10 in CD4T and CD8T; 5 setups [88_spike_in_markers.R]
#          experiments/runs/2026-11-10_bmv-spike-in-panel2/spike_summary.csv, spike_cells.csv
#            9 cytokines and enzymes in CD4T and CD8T; 2 setups [88b_spike_in_panel2.R]
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg/spike_summary.csv, spike_cells.csv
#            12 markers in CD3- cells; 2 setups [88e_spike_cd3neg.R]
#          experiments/runs/2026-11-11_bmv-cospike-tnf-glut1/cospike_summary.csv
#            TNF +1 with GLUT1 +0/+1/+2 [88d_cospike_tnf_glut1.R]
#          experiments/runs/2026-11-11_bmv-tnf-spike-by-glut1/summary_by_glut1.csv
#            TNF +1 loss per natural GLUT1 quartile [88f_tnf_spike_by_glut1.R]
# Outputs: experiments/runs/2026-11-11_bmv-restyle/18_spike_in/fig_spike_<marker>.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/18_spike_in/fig_spike_recovery_ecdf.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/20_spike_in_panel2/fig_spike_<marker>.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/20_spike_in_panel2/fig_spike_recovery_ecdf.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/22_spike_cd3neg/fig_spike_{GLUT1,TNF}.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/22_spike_cd3neg/fig_spike_recovery_ecdf.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/23_cospike/fig_cospike_tnf_loss.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/23_cospike/fig_cospike_glut1_recovery.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/23_cospike/by_glut1_fig_loss_by_glut1.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/manifest_F3.csv   summary path, new path, size
#          experiments/runs/2026-11-11_bmv-restyle/sessionInfo_F3.txt
#          The ECDF .csv holds the ECDF evaluated on ECDF_GRID per setup, marker and population (the
#          plotted frame), not the per-cell table, which stays in the spike run folders.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))
source(file.path(LIB, "plots.R"))

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- "experiments/runs/2026-11-11_bmv-restyle"
RUNS <- "experiments/runs"
FIG_LEVEL <- 1 # the level the figures show (user, 2026-11-10: +1 arcsinh), as 88
LOSS_FRAC <- 0.5 # "losing most of the spike" = recovery below 0.5, as 88
SEED <- 42 # figure jitter, as 88
ECDF_GRID <- seq(-0.5, 1.5, by = 0.0025) # recovery values where the ECDF is evaluated; = 88's x window
# marker = dye, from the control file (as 88, 88b, 88e)
MARKERS <- c(
  TNF = "BV510",
  Th2 = "eFluor 450",
  IL10 = "BB700",
  IFNg = "BV570",
  IL2 = "BV605",
  GLUT1 = "Alexa Fluor 405",
  ATP5A = "Alexa Fluor 488",
  CS = "Alexa Fluor 594",
  CPT1A = "PE-Cy5.5",
  ATGL = "PE-Cy7",
  G6PD = "APC-Cy7",
  CD71 = "Spark YG 593"
)
# setup labels and their order, as 88/88b/88e; one palette over all five keeps colours fixed
SETUP_LABELS <- c(
  c188 = "as188_joint",
  s102_c102 = "as102_var",
  s188_c102 = "s188_c102",
  opt188 = "opt188",
  afrefine = "af405_afrefine"
)
SETUP_LINES <- c(
  as188_joint = "solid",
  as102_var = "dashed",
  s188_c102 = "dotted",
  opt188 = "dotdash",
  af405_afrefine = "longdash"
) # second encoding of setup in the ECDFs; every label has its own
# 88d and 88f label the setups by name
COSPIKE_LABELS <- c(
  c188 = "1.8.8 zip",
  afrefine = "refined AF",
  s102_c102 = "1.0.2 zip"
)
BY_GLUT1_LABELS <- c(
  c188 = "1.8.8 zip",
  afrefine = "refined AF",
  opt188 = "opt188",
  s188_c102 = "1.0.2 code,\n1.8.8 spectra",
  s102_c102 = "1.0.2 zip"
)
DONOR_SHAPES <- c(A3 = 16, A6 = 17, B11 = 15, B2 = 18, D1 = 8)
# the figures report section 18, 20, 22 and 23 cites (90_report.qmd show_kept), plus 22's ECDF
SPIKE_FIGURES <- tribble(
  ~section             , ~run                             , ~marker ,
  "18_spike_in"        , "2026-11-10_bmv-spike-in-v2"     , "TNF"   ,
  "18_spike_in"        , "2026-11-10_bmv-spike-in-v2"     , "Th2"   ,
  "18_spike_in"        , "2026-11-10_bmv-spike-in-v2"     , "IL10"  ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" , "GLUT1" ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" , "CS"    ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" , "IFNg"  ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" , "IL2"   ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" , "ATP5A" ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" , "CPT1A" ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" , "ATGL"  ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" , "G6PD"  ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" , "CD71"  ,
  "22_spike_cd3neg"    , "2026-11-11_bmv-spike-in-cd3neg" , "GLUT1" ,
  "22_spike_cd3neg"    , "2026-11-11_bmv-spike-in-cd3neg" , "TNF"
)
ECDF_FIGURES <- tribble(
  ~section             , ~run                             ,
  "18_spike_in"        , "2026-11-10_bmv-spike-in-v2"     ,
  "20_spike_in_panel2" , "2026-11-10_bmv-spike-in-panel2" ,
  "22_spike_cd3neg"    , "2026-11-11_bmv-spike-in-cd3neg"
)

# Output folders: other restyle forks write next to these, so check files, not folders
planned <- c(
  file.path(
    OUT_DIR,
    SPIKE_FIGURES$section,
    paste0("fig_spike_", SPIKE_FIGURES$marker, ".png")
  ),
  file.path(OUT_DIR, ECDF_FIGURES$section, "fig_spike_recovery_ecdf.png"),
  file.path(
    OUT_DIR,
    "23_cospike",
    c(
      "fig_cospike_tnf_loss.png",
      "fig_cospike_glut1_recovery.png",
      "by_glut1_fig_loss_by_glut1.png"
    )
  )
)
stopifnot("restyled figures must not exist yet" = !any(file.exists(planned)))
walk(unique(dirname(planned)), \(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
})

# Per-marker figures: level-1 spike, donors as points over the median, setups down the y axis -----------
# Two readouts side by side; flipped (r-0012) so the setup names read horizontally.
draw_spike <- function(section, run, which_marker) {
  summary_tbl <- read_csv(
    file.path(RUNS, run, "spike_summary.csv"),
    show_col_types = FALSE
  )
  setup_pal <- categorical_palette(unname(SETUP_LABELS))
  # keep the shown level and marker, one row per readout
  panel_data <- summary_tbl |>
    filter(level == FIG_LEVEL, marker == which_marker) |>
    mutate(
      setup = factor_by_pal(SETUP_LABELS[setup], setup_pal),
      `Cells losing >50% of spike (%)` = 100 * share_lost,
      `Median recovery (fraction of spike)` = median_recovery
    ) |>
    pivot_longer(
      c(
        `Cells losing >50% of spike (%)`,
        `Median recovery (fraction of spike)`
      ),
      names_to = "readout",
      values_to = "value"
    ) |>
    mutate(setup = droplevels(setup))
  n_per_setup <- filter(
    panel_data,
    readout == first(readout),
    pop == first(pop)
  )$setup
  n_setups <- nlevels(panel_data$setup)
  n_pops <- n_distinct(panel_data$pop)
  # 89 mm fits two setups; five setups with their n labels need the double column
  width <- if (n_setups <= 2) 89 else 183
  n_breaks <- if (width == 89) 3 else 5 # x ticks per panel
  # display floor per readout: % lost starts at 0, recovery spans at least 0.99-1.01, so panels whose
  # values differ by 1e-7 still get readable ticks; no data point moves
  span_floor <- tibble(
    readout = c(
      "Cells losing >50% of spike (%)",
      "Median recovery (fraction of spike)",
      "Median recovery (fraction of spike)"
    ),
    value = c(0, 0.99, 1.01)
  )
  p <- ggplot(panel_data, aes(value, setup)) +
    geom_blank(data = span_floor, aes(x = value), inherit.aes = FALSE) +
    stat_summary(
      aes(colour = setup),
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.5,
      linewidth = FIG_LINE / 2,
      orientation = "y"
    ) +
    geom_point(
      aes(colour = setup, shape = sample),
      position = position_jitter(width = 0, height = 0.15, seed = SEED),
      size = 1.1
    ) +
    facet_grid(
      pop ~ readout,
      scales = "free_x",
      axes = "all",
      axis.labels = "margins",
      labeller = labeller(readout = label_wrap_gen(22))
    ) +
    scale_colour_manual(values = setup_pal, guide = "none") +
    scale_shape_manual(values = DONOR_SHAPES, name = "Donor") +
    # first setup on top
    scale_y_discrete(limits = rev, labels = n_labels(n_per_setup)) +
    labs(
      x = NULL,
      y = NULL,
      title = sprintf(
        "%s (%s) spike +%d arcsinh",
        which_marker,
        MARKERS[[which_marker]],
        FIG_LEVEL
      )
    ) +
    # n_breaks ticks per panel; each tick label gets the fewest significant digits that keep the ticks
    # of that panel distinct (recovery panels can span 1e-7), and panels sit 4 mm apart
    scale_x_continuous(
      breaks = scales::breaks_extended(n = n_breaks),
      labels = \(x) {
        shown <- x[!is.na(x)]
        digits <- 2
        while (digits < 9 && anyDuplicated(signif(shown, digits))) {
          digits <- digits + 1
        }
        as.character(signif(x, digits))
      },
      guide = guide_axis(check.overlap = TRUE)
    ) +
    theme(panel.spacing.x = unit(4, "mm")) +
    value_grid("x")
  # rows set the height; at least 55 mm so the five-donor legend fits
  height <- max(55, 18 + n_pops * (6 + 7 * n_setups))
  save_panel(
    p,
    panel_data,
    file.path(OUT_DIR, section),
    paste0("fig_spike_", which_marker),
    width,
    height
  )
  tibble(
    summary_path = file.path(
      section,
      paste0("fig_spike_", which_marker, ".png")
    ),
    width_mm = width,
    height_mm = height
  )
}
spike_manifest <- pmap(SPIKE_FIGURES, \(section, run, marker) {
  draw_spike(section, run, marker)
}) |>
  list_rbind()

# Recovery ECDFs: every level-1 cell, pooled over donors, one curve per setup ------------------------------
# The ECDF is evaluated on ECDF_GRID so the source table holds the plotted curve, not millions of cells.
draw_ecdf <- function(section, run) {
  spiked <- read_csv(
    file.path(RUNS, run, "spike_cells.csv"),
    show_col_types = FALSE,
    col_select = c(setup, marker, pop, level, recovery)
  ) |>
    filter(level == FIG_LEVEL)
  setup_pal <- categorical_palette(unname(SETUP_LABELS))
  stopifnot(
    "every setup has a label and a linetype" = all(
      SETUP_LABELS[unique(spiked$setup)] %in%
        names(SETUP_LINES)
    )
  )
  # ECDF of recovery per setup, marker and population, read off on the grid
  ecdf_frame <- spiked |>
    # the grid gets its own name so ecdf() reads the cells' recovery, not the grid
    reframe(
      n_cells = n(),
      cumulative = ecdf(recovery)(ECDF_GRID),
      recovery_grid = ECDF_GRID,
      .by = c(setup, marker, pop)
    ) |>
    mutate(
      setup = droplevels(factor_by_pal(SETUP_LABELS[setup], setup_pal)),
      marker = factor(
        marker,
        levels = intersect(names(MARKERS), unique(marker))
      )
    )
  n_markers <- nlevels(ecdf_frame$marker)
  n_pops <- n_distinct(ecdf_frame$pop)
  # one population: wrap the markers four across; two populations: markers down, populations across
  facets <- if (n_pops == 1) {
    facet_wrap(~marker, ncol = 4, axes = "all", axis.labels = "margins")
  } else {
    facet_grid(marker ~ pop, axes = "all", axis.labels = "margins")
  }
  n_rows <- if (n_pops == 1) ceiling(n_markers / 4) else n_markers
  p <- ggplot(
    ecdf_frame,
    aes(recovery_grid, cumulative, colour = setup, linetype = setup)
  ) +
    geom_vline(
      xintercept = c(LOSS_FRAC, 1),
      linetype = "dashed",
      colour = "grey60",
      linewidth = FIG_LINE / 2
    ) +
    geom_step(linewidth = FIG_LINE) +
    facets +
    scale_colour_manual(values = setup_pal, name = "Setup") +
    scale_linetype_manual(values = SETUP_LINES, name = "Setup") +
    scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5)) +
    labs(
      x = "Recovered fraction of the +1 spike",
      y = "Cumulative fraction of spiked cells",
      title = sprintf(
        "Spike recovery per cell (%s)",
        paste(unique(ecdf_frame$pop), collapse = ", ")
      )
    ) +
    value_grid("x")
  height <- 20 + n_rows * 22
  save_panel(
    p,
    ecdf_frame,
    file.path(OUT_DIR, section),
    "fig_spike_recovery_ecdf",
    183,
    height
  )
  tibble(
    summary_path = file.path(section, "fig_spike_recovery_ecdf.png"),
    width_mm = 183,
    height_mm = height
  )
}
ecdf_manifest <- pmap(ECDF_FIGURES, draw_ecdf) |>
  list_rbind()

# Co-spike: TNF +1 with GLUT1 +0/+1/+2, donors as points, median line per setup --------------------------
cospike <- read_csv(
  file.path(RUNS, "2026-11-11_bmv-cospike-tnf-glut1", "cospike_summary.csv"),
  show_col_types = FALSE
) |>
  mutate(
    setup_label = factor(COSPIKE_LABELS[setup], levels = COSPIKE_LABELS),
    glut1_label = factor(paste0("+", glut1_level), levels = paste0("+", 0:2)),
    donor = sample
  )
pal_cospike <- categorical_palette(levels(cospike$setup_label))
p_tnf <- ggplot(cospike, aes(glut1_label, pct_tnf_lost, colour = setup_label)) +
  stat_summary(
    aes(group = setup_label),
    fun = median,
    geom = "line",
    linewidth = FIG_LINE,
    show.legend = FALSE
  ) +
  geom_point(
    aes(shape = donor),
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 1.1
  ) +
  facet_grid(pop ~ setup_label, axes = "all", axis.labels = "margins") +
  scale_colour_manual(values = pal_cospike, guide = "none") +
  scale_shape_manual(values = DONOR_SHAPES) +
  labs(
    x = "GLUT1 (Alexa Fluor 405) added with TNF +1 arcsinh (arcsinh units)",
    y = "Cells losing >50% of the TNF spike (%)",
    shape = "Donor",
    title = "TNF + GLUT1 co-spike: TNF loss"
  ) +
  value_grid("y")
save_panel(
  p_tnf,
  cospike,
  file.path(OUT_DIR, "23_cospike"),
  "fig_cospike_tnf_loss",
  183,
  80
)

glut1_df <- filter(cospike, glut1_level > 0)
p_glut1 <- ggplot(
  glut1_df,
  aes(glut1_label, median_glut1_recovery, colour = setup_label)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "grey50",
    linewidth = FIG_LINE / 2
  ) +
  stat_summary(
    aes(group = setup_label),
    fun = median,
    geom = "line",
    linewidth = FIG_LINE,
    show.legend = FALSE
  ) +
  geom_point(
    aes(shape = donor),
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 1.1
  ) +
  facet_grid(pop ~ setup_label, axes = "all", axis.labels = "margins") +
  scale_colour_manual(values = pal_cospike, guide = "none") +
  scale_shape_manual(values = DONOR_SHAPES) +
  labs(
    x = "GLUT1 (Alexa Fluor 405) added with TNF +1 arcsinh (arcsinh units)",
    y = "Median GLUT1 recovery (fraction of spike)",
    shape = "Donor",
    title = "TNF + GLUT1 co-spike: GLUT1 recovery"
  ) +
  value_grid("y")
save_panel(
  p_glut1,
  glut1_df,
  file.path(OUT_DIR, "23_cospike"),
  "fig_cospike_glut1_recovery",
  183,
  80
)

# TNF +1 loss by natural GLUT1 quartile, donors as points, median line per setup -------------------------
# Tick labels shortened to Q1-Q4 (the axis title names the quartiles) so they need no rotation.
by_glut1 <- read_csv(
  file.path(RUNS, "2026-11-11_bmv-tnf-spike-by-glut1", "summary_by_glut1.csv"),
  show_col_types = FALSE
) |>
  mutate(
    setup_label = factor(BY_GLUT1_LABELS[setup], levels = BY_GLUT1_LABELS),
    glut1_bin = factor(sub("GLUT1 ", "", glut1_bin), levels = paste0("Q", 1:4)),
    donor = sample
  )
pal_by_glut1 <- categorical_palette(levels(by_glut1$setup_label))
p_loss <- ggplot(by_glut1, aes(glut1_bin, pct_lost, colour = setup_label)) +
  stat_summary(
    aes(group = setup_label),
    fun = median,
    geom = "line",
    linewidth = FIG_LINE,
    show.legend = FALSE
  ) +
  geom_point(
    aes(shape = donor),
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 1.1
  ) +
  facet_grid(pop ~ setup_label, axes = "all", axis.labels = "margins") +
  scale_colour_manual(values = pal_by_glut1, guide = "none") +
  scale_shape_manual(values = DONOR_SHAPES) +
  labs(
    x = "Natural GLUT1 (Alexa Fluor 405), quartile per donor",
    y = "Cells losing >50% of the TNF spike (%)",
    shape = "Donor",
    title = "TNF +1 arcsinh spike by natural GLUT1"
  ) +
  value_grid("y")
save_panel(
  p_loss,
  by_glut1,
  file.path(OUT_DIR, "23_cospike"),
  "by_glut1_fig_loss_by_glut1",
  183,
  85
)

# Manifest: summary copy each figure replaces, and its new path ---------------------------------------
manifest <- bind_rows(
  spike_manifest,
  ecdf_manifest,
  tibble(
    summary_path = file.path(
      "23_cospike",
      c(
        "fig_cospike_tnf_loss.png",
        "fig_cospike_glut1_recovery.png",
        "by_glut1_fig_loss_by_glut1.png"
      )
    ),
    width_mm = 183,
    height_mm = c(80, 80, 85)
  )
) |>
  mutate(
    new_path = file.path(OUT_DIR, summary_path),
    script = "experiments/benchmarks/05_bmv-unmix-compare/95_restyle_spike.R",
    .after = summary_path
  )
write_csv(manifest, file.path(OUT_DIR, "manifest_F3.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo_F3.txt")
)
