# Purpose: redraw the co-spike raw hexbin the report cites (CD8T) from the common-gate co-spike run (plans/common_gates_plan.md
#          G5, G10b, G11). Copy of 95_restyle_events_cospike_raw.R with the source repointed to bmv-cospike-tnf-glut1-cg,
#          setups named through arm_label(), a density colourbar on the hexagons, and axis limits taken from the data.
#          TNF before against after the TNF +1 spike with GLUT1 co-spiked at +0, +1, +2; donors pooled (per condition
#          when the dataset has conditions); the % is the donor median of cells losing more than half of the TNF spike.
#          TNF is the dataset's tnf_label (the BV510 slot).
# Run:     SW_DATASET=3d SW_RUN_DATE=YYYY-MM-DD REDRAW_ONLY=spike pixi run Rscript
#            experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_events_cospike_raw.R from the repo root.
# Inputs:  run_in("cospike-tnf-glut1-cg")/cospike_cells.csv     per co-spiked cell
#          run_in("cospike-tnf-glut1-cg")/cospike_summary.csv   % TNF lost per donor
#          datasets/<SW_DATASET>.json   lib_dirs (ggplot2 3.5 entry), targets, conditions
# Outputs: run_out("restyle-common") (shared by the 95c scripts)/23_cospike/fig_cospike_raw_CD8T_v2.{pdf,png,csv}   the .csv holds the donor-median %
#          run_out("restyle-common")/manifest_parts/95c_restyle_events_cospike_raw.R.csv

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
# ggplot2 3.5.2 for ggcyto: only the ggplot2_35 entry of lib_dirs, as before
.libPaths(c(
  grep("lib_ggplot2_35$", cfg("lib_dirs"), value = TRUE),
  .libPaths()
))
suppressPackageStartupMessages({
  library(flowCore)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/95c_spike_common.R")
SCRIPT <- "95c_restyle_events_cospike_raw.R"
STEM <- "fig_cospike_raw_CD8T_v2"
TNF_LAB <- cfg("targets")$tnf_label
if (!want(STEM)) {
  quit(save = "no")
}

# Parameters ---------------------------------------------------------------------------------------
POP <- "CD8T" # the population the report cites
GLUT1_LEVELS <- c(0, 1, 2) # GLUT1 co-spike, arcsinh units (88d)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
AXIS_STEP <- 1 # axis limits are the data range rounded outward to this step (arcsinh); no cell is cut
BINS <- 60 # hexagons per panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
WIDTH_MM <- 183 # Nature double column, three setups
cells_file <- SRC("cospike-tnf-glut1-cg", "cospike_cells.csv")
summary_file <- SRC("cospike-tnf-glut1-cg", "cospike_summary.csv")
out_dir <- fresh_dir("23_cospike")

# Cells and the donor-median % lost per setup and GLUT1 level (as 88d) -------------------------------------
glut1_levels <- paste0("GLUT1 +", GLUT1_LEVELS)
setup_levels <- SETUP_ORDER[
  SETUP_ORDER %in% setup_label(c("c188", "afrefine", "s102_c102"))
]
cospike_summary <- read_csv(summary_file, show_col_types = FALSE) |>
  need_cond(SCRIPT) |>
  dplyr::filter(.data$pop == POP)
N_DONORS <- n_distinct(sample_keys(cospike_summary)$donor_key) # donors, not samples (24H)
labels <- cospike_summary |>
  mutate(
    setup_label = factor(setup_label(.data$setup), levels = setup_levels),
    glut1_label = factor(
      paste0("GLUT1 +", .data$glut1_level),
      levels = glut1_levels
    ),
    across(any_of(COND), \(x) forcats::fct(x, levels = DS$conditions))
  ) |>
  summarise(
    pct = median(.data$pct_tnf_lost),
    .by = c(setup_label, glut1_label, any_of(COND))
  )
cells <- data.table::fread(
  cells_file,
  select = c("pop", "tnf_before", "tnf_recovered", "setup", "glut1_level", COND)
) |>
  as_tibble() |>
  dplyr::filter(.data$pop == POP) |>
  mutate(
    across(any_of(COND), \(x) forcats::fct(x, levels = DS$conditions)),
    before = asinh(tnf_before / COFACTOR),
    after = asinh((tnf_before + tnf_recovered) / COFACTOR),
    setup_label = factor(setup_label(.data$setup), levels = setup_levels),
    glut1_label = paste0("GLUT1 +", .data$glut1_level)
  )
frames <- cells |>
  mutate(
    name = paste(
      .data$setup_label,
      .data$glut1_label,
      !!!rlang::syms(COND),
      sep = "|"
    )
  ) |>
  split(~name) |>
  map(\(x) flowFrame(cbind(before = x$before, after = x$after)))
fs <- flowSet(frames)
pd <- tibble(name = sampleNames(fs)) |>
  tidyr::separate_wider_delim(
    name,
    "|",
    names = c("setup_label", "glut1_label", COND),
    cols_remove = FALSE
  ) |>
  mutate(
    setup_label = factor(.data$setup_label, levels = setup_levels),
    glut1_label = factor(.data$glut1_label, levels = glut1_levels),
    across(any_of(COND), \(x) forcats::fct(x, levels = DS$conditions))
  )
pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)

# Reference curves: full recovery adds sinh(1) in linear units; loss stays on the diagonal ----------------
axis_lim <- c(
  floor(min(cells$before, cells$after) / AXIS_STEP) * AXIS_STEP,
  ceiling(max(cells$before, cells$after) / AXIS_STEP) * AXIS_STEP
)
grid_x <- seq(axis_lim[1], axis_lim[2], length.out = 200)
curves <- bind_rows(
  tibble(
    x = grid_x,
    y = asinh(sinh(grid_x) + sinh(FIG_LEVEL)),
    curve = "full recovery"
  ),
  tibble(x = grid_x, y = grid_x, curve = "spike lost")
)

# Figure ---------------------------------------------------------------------------------------------------
p <- ggcyto(fs, aes(x = before, y = after)) +
  geom_hex(bins = BINS) +
  geom_line(
    data = curves,
    aes(x, y, linetype = curve),
    inherit.aes = FALSE,
    colour = "grey30",
    linewidth = FIG_LINE
  ) +
  geom_label(
    data = labels,
    aes(
      x = -Inf,
      y = Inf,
      label = sprintf("%s lost >50%%: %.1f%%", TNF_LAB, pct)
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = LABEL_PT,
    size.unit = "pt",
    family = FIG_FONT,
    fill = "white",
    alpha = 0.85,
    label.size = 0,
    label.padding = unit(0.1, "lines")
  ) +
  facet_grid(
    rows = vars(glut1_label, !!!rlang::syms(COND)),
    cols = vars(setup_label),
    switch = "y",
    axes = "all",
    axis.labels = "margins",
    labeller = labeller(setup_label = label_wrap_gen(16))
  ) +
  ggcyto_par_set(limits = list(x = axis_lim, y = axis_lim)) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    trans = "sqrt",
    name = "Cells per\nhexagon",
    guide = guide_colourbar()
  ) +
  scale_linetype_manual(
    values = c("full recovery" = "dashed", "spike lost" = "solid"),
    name = NULL
  ) +
  labs(
    x = paste0(TNF_LAB, ", unspiked unmix (arcsinh)"),
    y = paste0(TNF_LAB, ", co-spiked unmix (arcsinh)"),
    title = sprintf(
      "%s: %s +1 with GLUT1 co-spike, per cell\n%d donors pooled; %% = donor median",
      POP,
      TNF_LAB,
      N_DONORS
    )
  ) +
  theme(
    strip.text.y.left = element_text(face = "bold", angle = 0),
    strip.placement = "outside",
    panel.spacing = unit(1.2, "mm"),
    # ggcyto resets the legend position, so it is set here
    legend.position = "right",
    legend.key.width = unit(3, "mm"),
    legend.key.height = unit(14, "mm"),
    # both axes are TNF in arcsinh, so the panels are square
    aspect.ratio = 1
  )
# build on a cairo device, as save_panel() draws, so Nimbus Sans resolves through fontconfig
cairo_pdf(tempfile(fileext = ".pdf"))
g <- ggplotGrob(as.ggplot(p))
dev.off()
height_mm <- round(
  (WIDTH_MM - 60) /
    length(setup_levels) *
    length(GLUT1_LEVELS) *
    max(1, length(DS$conditions)) +
    40
)
save_panel(g, labels, out_dir, STEM, WIDTH_MM, height_mm)
add_manifest(
  SCRIPT,
  file.path("2026-11-11_bmv-restyle/23_cospike", paste0(STEM, ".png")),
  file.path(out_dir, paste0(STEM, ".png")),
  cells_file,
  note = sprintf(
    "axis window %g to %g arcsinh from the data",
    axis_lim[1],
    axis_lim[2]
  )
)
