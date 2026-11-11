# Purpose: Nature-spec redraw (user, 2026-11-11; lucid-coder r-0006, r-0008, r-0015, r-0016, r-0017) of the
#          per-cell spike-in views from 88c_spike_raw_scatter.R (CD4T/CD8T) and 88g_spike_raw_scatter_cd3neg.R
#          (live CD3-), for the figures the report cites. Same cells (+1 arcsinh level), same % losing >50% of
#          the spike; only the layout changes: theme_prism 7 pt Nimbus Sans, axes on every panel, labels on
#          the outer edges, square panels (both axes are the same marker in arcsinh), legend right.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_events_spike_raw.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-10_bmv-spike-in-v2/spike_cells.csv         TNF, Th2, IL-10; 5 setups [88]
#          experiments/runs/2026-11-10_bmv-spike-in-panel2/spike_cells.csv     9 markers; 2 setups [88b]
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg/spike_cells.csv     12 markers; 2 setups [88e]
# Outputs: experiments/runs/2026-11-11_bmv-restyle/21_raw_companions/spike_raw/fig_raw_<marker>_<pop>.pdf,
#            .png, .csv   CD4T/CD8T figures; the .csv holds the % losing >50% per panel
#          experiments/runs/2026-11-11_bmv-restyle/22_spike_cd3neg/raw/fig_raw_<marker>_CD3neg.pdf, .png, .csv

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
OUT_ROOT <- "experiments/runs/2026-11-11_bmv-restyle"
# the figures 90_report.qmd cites: marker, population, source table, output folder
FIGURES <- tribble(
  ~marker , ~pop     , ~source                                                           , ~folder                       ,
  "TNF"   , "CD8T"   , "experiments/runs/2026-11-10_bmv-spike-in-v2/spike_cells.csv"     , "21_raw_companions/spike_raw" ,
  "Th2"   , "CD8T"   , "experiments/runs/2026-11-10_bmv-spike-in-v2/spike_cells.csv"     , "21_raw_companions/spike_raw" ,
  "CS"    , "CD8T"   , "experiments/runs/2026-11-10_bmv-spike-in-panel2/spike_cells.csv" , "21_raw_companions/spike_raw" ,
  "TNF"   , "CD3neg" , "experiments/runs/2026-11-11_bmv-spike-in-cd3neg/spike_cells.csv" , "22_spike_cd3neg/raw"         ,
  "GLUT1" , "CD3neg" , "experiments/runs/2026-11-11_bmv-spike-in-cd3neg/spike_cells.csv" , "22_spike_cd3neg/raw"
)
SETUP_LABELS <- c(
  c188 = "1.8.8 zip",
  afrefine = "refined AF",
  opt188 = "opt188",
  s188_c102 = "1.0.2 code,\n1.8.8 spectra",
  s102_c102 = "1.0.2 zip"
)
POP_LABELS <- c(CD4T = "CD4T", CD8T = "CD8T", CD3neg = "Live CD3-")
LEVEL <- 1 # the +1 arcsinh cells, as the summary figures (88c)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
LOSS_FRAC <- 0.5 # "losing most of the spike", as 88 and 88c
AXIS <- c(-1.5, 4) # display window only
BINS <- 50 # hexagons per panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt

stems <- file.path(
  OUT_ROOT,
  FIGURES$folder,
  sprintf("fig_raw_%s_%s.png", FIGURES$marker, FIGURES$pop)
)
stopifnot("no output file exists yet" = !any(file.exists(stems)))

# Reference curves: full recovery adds sinh(1) in linear units; loss stays on the diagonal ----------------
grid_x <- seq(AXIS[1], AXIS[2], length.out = 200)
curves <- bind_rows(
  tibble(
    x = grid_x,
    y = asinh(sinh(grid_x) + sinh(LEVEL)),
    curve = "full recovery"
  ),
  tibble(x = grid_x, y = grid_x, curve = "spike lost")
)

# One ggcyto figure per cited marker and population -------------------------------------------------------
pwalk(FIGURES, \(marker, pop, source, folder) {
  # +1 arcsinh cells of this marker and population, on arcsinh (as 88c)
  cells <- read_csv(source, show_col_types = FALSE) |>
    dplyr::filter(
      level == LEVEL,
      .data$marker == .env$marker,
      .data$pop == .env$pop
    ) |>
    mutate(
      before = asinh(base_value / COFACTOR),
      after = asinh((base_value + recovered) / COFACTOR),
      lost = recovery < LOSS_FRAC,
      setup_label = factor(
        SETUP_LABELS[setup],
        levels = intersect(SETUP_LABELS, SETUP_LABELS[unique(setup)])
      )
    )
  pct_lost <- cells |>
    summarise(
      pct_lost = 100 * mean(lost),
      n = n(),
      .by = c(marker, pop, sample, setup, setup_label)
    )
  # one flowFrame per setup and donor, annotated for the facets
  frames <- cells |>
    mutate(name = paste(setup_label, sample, sep = "|")) |>
    split(~name) |>
    map(\(x) flowFrame(cbind(before = x$before, after = x$after)))
  fs <- flowSet(frames)
  pd <- tibble(name = sampleNames(fs)) |>
    tidyr::separate_wider_delim(
      name,
      "|",
      names = c("setup_label", "sample"),
      cols_remove = FALSE
    ) |>
    mutate(
      setup_label = factor(setup_label, levels = levels(cells$setup_label))
    )
  pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)

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
      data = pct_lost,
      aes(x = -Inf, y = Inf, label = sprintf("lost >50%%: %.1f%%", pct_lost)),
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
      sample ~ setup_label,
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
      x = sprintf("%s, unspiked unmix (arcsinh, cofactor 6000)", marker),
      y = sprintf("%s, +1 arcsinh spiked unmix", marker),
      title = sprintf("%s: %s spike, per cell", POP_LABELS[[pop]], marker)
    ) +
    theme(
      strip.text.y.left = element_text(face = "bold", angle = 0),
      strip.placement = "outside",
      panel.spacing = unit(1.2, "mm"),
      legend.position = "right",
      legend.key.width = unit(6, "mm"),
      # both axes are the same marker in arcsinh, so the panels are square (r-0006)
      aspect.ratio = 1
    )
  # Nature single column for up to 2 setups, double for more; ~40 mm for the y axis and the legend
  n_setups <- nlevels(pd$setup_label)
  width_mm <- if (n_setups <= 2) 89 else 183
  height_mm <- round((width_mm - 40) / n_setups * n_distinct(pd$sample) + 28)
  out_dir <- file.path(OUT_ROOT, folder)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # build on a cairo device, as save_panel() draws, so Nimbus Sans resolves through fontconfig
  cairo_pdf(tempfile(fileext = ".pdf"))
  g <- ggplotGrob(as.ggplot(p))
  dev.off()
  save_panel(
    g,
    select(pct_lost, -setup_label),
    out_dir,
    sprintf("fig_raw_%s_%s", marker, pop),
    width_mm,
    height_mm
  )
})
