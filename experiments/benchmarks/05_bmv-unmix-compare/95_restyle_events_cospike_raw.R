# Purpose: Nature-spec redraw (user, 2026-11-11; lucid-coder r-0006, r-0008, r-0015, r-0016, r-0017) of the
#          co-spike raw figure from 88d_cospike_tnf_glut1.R (mode figures_raw, CD8T): TNF before vs after the
#          TNF +1 spike with GLUT1 co-spiked at +0, +1, +2, five donors pooled, % = donor median of cells
#          losing >50% of the TNF spike. Same cells and numbers as 88d (its saved tables); only the layout
#          changes: theme_prism 7 pt Nimbus Sans, axes on every panel, labels on the outer edges, square
#          panels, legend right, 183 mm wide.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_events_cospike_raw.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-11_bmv-cospike-tnf-glut1/cospike_cells.csv     per co-spiked cell [88d]
#          experiments/runs/2026-11-11_bmv-cospike-tnf-glut1/cospike_summary.csv   % lost per donor [88d]
# Outputs: experiments/runs/2026-11-11_bmv-restyle/23_cospike/fig_cospike_raw_CD8T_v2.pdf, .png, .csv
#            the .csv holds the donor-median % lost per panel

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto
suppressPackageStartupMessages({
  library(flowCore)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters ---------------------------------------------------------------------------------------
RUN_88D <- "experiments/runs/2026-11-11_bmv-cospike-tnf-glut1"
OUT_DIR <- "experiments/runs/2026-11-11_bmv-restyle/23_cospike"
STEM <- "fig_cospike_raw_CD8T_v2"
POP <- "CD8T" # the population the report cites
TNF_LEVEL <- 1 # TNF spike, arcsinh units (88d)
GLUT1_LEVELS <- c(0, 1, 2) # GLUT1 co-spike, arcsinh units (88d)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
SETUP_LABELS <- c(
  c188 = "1.8.8 zip",
  afrefine = "refined AF",
  s102_c102 = "1.0.2 zip"
)
AXIS <- c(-1.5, 4) # display window only
BINS <- 50 # hexagons per panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
WIDTH_MM <- 183 # Nature double column, three setups

stopifnot(
  "no output file exists yet" = !file.exists(file.path(
    OUT_DIR,
    paste0(STEM, ".png")
  ))
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Cells and the donor-median % lost per setup and GLUT1 level (as 88d) -------------------------------------
glut1_levels <- paste0("GLUT1 +", GLUT1_LEVELS)
labels <- read_csv(
  file.path(RUN_88D, "cospike_summary.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(pop == POP) |>
  mutate(
    setup_label = factor(SETUP_LABELS[setup], levels = SETUP_LABELS),
    glut1_label = factor(paste0("GLUT1 +", glut1_level), levels = glut1_levels)
  ) |>
  summarise(pct = median(pct_tnf_lost), .by = c(setup_label, glut1_label))
cells <- read_csv(
  file.path(RUN_88D, "cospike_cells.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(pop == POP) |>
  mutate(
    before = asinh(tnf_before / COFACTOR),
    after = asinh((tnf_before + tnf_recovered) / COFACTOR),
    setup_label = factor(SETUP_LABELS[setup], levels = SETUP_LABELS),
    glut1_label = paste0("GLUT1 +", glut1_level)
  )
frames <- cells |>
  mutate(name = paste(setup_label, glut1_label, sep = "|")) |>
  split(~name) |>
  map(\(x) flowFrame(cbind(before = x$before, after = x$after)))
fs <- flowSet(frames)
pd <- tibble(name = sampleNames(fs)) |>
  tidyr::separate_wider_delim(
    name,
    "|",
    names = c("setup_label", "glut1_label"),
    cols_remove = FALSE
  ) |>
  mutate(
    setup_label = factor(setup_label, levels = SETUP_LABELS),
    glut1_label = factor(glut1_label, levels = glut1_levels)
  )
pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)

# Reference curves: full recovery adds sinh(1) in linear units; loss stays on the diagonal ----------------
grid_x <- seq(AXIS[1], AXIS[2], length.out = 200)
curves <- bind_rows(
  tibble(
    x = grid_x,
    y = asinh(sinh(grid_x) + sinh(TNF_LEVEL)),
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
    aes(x = -Inf, y = Inf, label = sprintf("TNF lost >50%%: %.1f%%", pct)),
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
    glut1_label ~ setup_label,
    switch = "y",
    axes = "all",
    axis.labels = "margins"
  ) +
  ggcyto_par_set(limits = list(x = AXIS, y = AXIS)) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    trans = "sqrt",
    guide = "none"
  ) +
  scale_linetype_manual(
    values = c("full recovery" = "dashed", "spike lost" = "solid"),
    name = NULL
  ) +
  labs(
    x = "TNF, unspiked unmix (arcsinh)",
    y = "TNF, co-spiked unmix (arcsinh)",
    title = sprintf(
      "%s: TNF +1 with GLUT1 co-spike, per cell\n5 donors pooled; %% = donor median",
      POP
    )
  ) +
  theme(
    strip.text.y.left = element_text(face = "bold", angle = 0),
    strip.placement = "outside",
    panel.spacing = unit(1.2, "mm"),
    # ggcyto resets the legend position, so it is set here (r-0015)
    legend.position = "right",
    legend.key.width = unit(6, "mm"),
    # both axes are TNF in arcsinh, so the panels are square (r-0006)
    aspect.ratio = 1
  )
# build on a cairo device, as save_panel() draws, so Nimbus Sans resolves through fontconfig
cairo_pdf(tempfile(fileext = ".pdf"))
g <- ggplotGrob(as.ggplot(p))
dev.off()
height_mm <- round(
  (WIDTH_MM - 55) / length(SETUP_LABELS) * length(GLUT1_LEVELS) + 32
)
save_panel(g, labels, OUT_DIR, STEM, WIDTH_MM, height_mm)
