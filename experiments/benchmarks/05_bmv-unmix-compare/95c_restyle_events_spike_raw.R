# Purpose: redraw the per-cell spike hexbins the report cites from the common-gate spike runs (plans/common_gates_plan.md
#          G5, G10b, G11, G13). Copy of 95_restyle_events_spike_raw.R with sources repointed to the *-cg runs, setups named
#          through arm_label(), a density colourbar on the hexagons, and axis limits taken from the data so no cell is
#          cut (the old window -1.5 to 4 dropped about 0.2% of TNF cells). Same cells (+1 arcsinh level), same % losing
#          more than half of the spike. Also covers the n.passes figure of section 25 (fig_raw_Th2_CD8T).
# Run:     REDRAW_ONLY=spike pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_events_spike_raw.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-11_bmv-spike-in-cg/spike_cells.csv          TNF, Th2 in CD8T; 5 setups
#          experiments/runs/2026-11-11_bmv-spike-in-panel2-cg/spike_cells.csv    CS in CD8T; 2 setups
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg-cg/spike_cells.csv    TNF, GLUT1 in live CD3-; 2 setups
#          experiments/runs/2026-11-11_bmv-spike-passes-cg/spike_cells.csv       Th2 in CD8T; 4 pass settings
# Outputs: experiments/runs/2026-11-11_bmv-restyle-common/21_raw_companions/spike_raw/fig_raw_<marker>_CD8T.{pdf,png,csv}
#          .../22_spike_cd3neg/raw/fig_raw_<marker>_CD3neg.{pdf,png,csv}
#          .../25_spike_passes/raw/fig_raw_Th2_CD8T.{pdf,png,csv}   each .csv holds the % losing >50% per panel
#          .../manifest_parts/95c_restyle_events_spike_raw.R.csv

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto
suppressPackageStartupMessages({
  library(flowCore)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/95c_spike_common.R")
SCRIPT <- "95c_restyle_events_spike_raw.R"

# Parameters ---------------------------------------------------------------------------------------
# the figures 90_report.qmd cites: marker, population, run, column that names the unmix, folder, old figure
FIGURES <- tribble(
  ~marker , ~pop     , ~run                 , ~by       , ~folder                       , ~old                                                                      ,
  "TNF"   , "CD8T"   , "spike-in-cg"        , "setup"   , "21_raw_companions/spike_raw" , "2026-11-11_bmv-restyle/21_raw_companions/spike_raw/fig_raw_TNF_CD8T.png" ,
  "Th2"   , "CD8T"   , "spike-in-cg"        , "setup"   , "21_raw_companions/spike_raw" , "2026-11-11_bmv-restyle/21_raw_companions/spike_raw/fig_raw_Th2_CD8T.png" ,
  "CS"    , "CD8T"   , "spike-in-panel2-cg" , "setup"   , "21_raw_companions/spike_raw" , "2026-11-11_bmv-restyle/21_raw_companions/spike_raw/fig_raw_CS_CD8T.png"  ,
  "TNF"   , "CD3neg" , "spike-in-cd3neg-cg" , "setup"   , "22_spike_cd3neg/raw"         , "2026-11-11_bmv-restyle/22_spike_cd3neg/raw/fig_raw_TNF_CD3neg.png"       ,
  "GLUT1" , "CD3neg" , "spike-in-cd3neg-cg" , "setup"   , "22_spike_cd3neg/raw"         , "2026-11-11_bmv-restyle/22_spike_cd3neg/raw/fig_raw_GLUT1_CD3neg.png"     ,
  "Th2"   , "CD8T"   , "spike-passes-cg"    , "setting" , "25_spike_passes/raw"         , "2026-11-11_bmv-spike-passes/raw/fig_raw_Th2_CD8T.png"
)
POP_LABELS <- c(CD4T = "CD4T", CD8T = "CD8T", CD3neg = "Live CD3-")
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
AXIS_STEP <- 1 # axis limits are the data range rounded outward to this step (arcsinh); no cell is cut
BINS <- 60 # hexagons per panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt

cache <- new.env()
read_cells <- \(run, by) {
  key <- paste(run, by)
  if (is.null(cache[[key]])) {
    cache[[key]] <- data.table::fread(
      SRC(run, "spike_cells.csv"),
      select = c(
        by,
        "marker",
        "pop",
        "sample",
        "level",
        "base_value",
        "recovered",
        "recovery"
      )
    ) |>
      as_tibble() |>
      dplyr::filter(.data$level == FIG_LEVEL)
  }
  cache[[key]]
}

# One ggcyto figure per cited marker, population and run ---------------------------------------------------
pwalk(FIGURES, \(marker, pop, run, by, folder, old) {
  stem <- sprintf("fig_raw_%s_%s", marker, pop)
  if (!want(c(stem, paste(folder, stem)))) {
    return(invisible())
  }
  out_dir <- fresh_dir(folder)
  levels_shown <- if (by == "setup") SETUP_ORDER else unname(SETTING_ORDER)
  cells <- read_cells(run, by) |>
    dplyr::filter(.data$marker == .env$marker, .data$pop == .env$pop) |>
    mutate(
      before = asinh(base_value / COFACTOR),
      after = asinh((base_value + recovered) / COFACTOR),
      group = if (by == "setup") {
        setup_label(.data$setup)
      } else {
        unname(SETTING_ORDER[.data$setting])
      },
      group = factor(
        .data$group,
        levels = intersect(levels_shown, unique(.data$group))
      )
    )
  # square display window that holds every cell of the figure
  axis_lim <- c(
    floor(min(cells$before, cells$after) / AXIS_STEP) * AXIS_STEP,
    ceiling(max(cells$before, cells$after) / AXIS_STEP) * AXIS_STEP
  )
  # reference curves: full recovery adds sinh(1) in linear units; loss stays on the diagonal
  grid_x <- seq(axis_lim[1], axis_lim[2], length.out = 200)
  curves <- bind_rows(
    tibble(
      x = grid_x,
      y = asinh(sinh(grid_x) + sinh(FIG_LEVEL)),
      curve = "full recovery"
    ),
    tibble(x = grid_x, y = grid_x, curve = "spike lost")
  )
  pct_lost <- cells |>
    summarise(
      pct_lost = 100 * mean(.data$recovery < LOSS_FRAC),
      n = n(),
      .by = c(marker, pop, sample, group)
    )
  # one flowFrame per unmix and donor, annotated for the facets
  frames <- cells |>
    mutate(name = paste(.data$group, .data$sample, sep = "|")) |>
    split(~name) |>
    map(\(x) flowFrame(cbind(before = x$before, after = x$after)))
  fs <- flowSet(frames)
  pd <- tibble(name = sampleNames(fs)) |>
    tidyr::separate_wider_delim(
      name,
      "|",
      names = c("group", "sample"),
      cols_remove = FALSE
    ) |>
    mutate(group = factor(.data$group, levels = levels(cells$group)))
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
      sample ~ group,
      switch = "y",
      axes = "all",
      axis.labels = "margins",
      labeller = labeller(group = label_wrap_gen(14))
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
      x = sprintf("%s, unspiked unmix (arcsinh, cofactor 6000)", marker),
      y = sprintf("%s, +1 arcsinh spiked unmix", marker),
      title = sprintf("%s: %s spike, per cell", POP_LABELS[[pop]], marker)
    ) +
    theme(
      strip.text.y.left = element_text(face = "bold", angle = 0),
      strip.placement = "outside",
      panel.spacing = unit(1.2, "mm"),
      legend.position = "right",
      legend.key.width = unit(3, "mm"),
      legend.key.height = unit(14, "mm"),
      # both axes are the same marker in arcsinh, so the panels are square
      aspect.ratio = 1
    )
  # Nature single column for up to 2 unmixes, double for more
  n_groups <- nlevels(pd$group)
  width_mm <- if (n_groups <= 2) 89 else 183
  height_mm <- round((width_mm - 48) / n_groups * n_distinct(pd$sample) + 28)
  # build on a cairo device, as save_panel() draws, so Nimbus Sans resolves through fontconfig
  cairo_pdf(tempfile(fileext = ".pdf"))
  g <- ggplotGrob(as.ggplot(p))
  dev.off()
  save_panel(g, pct_lost, out_dir, stem, width_mm, height_mm)
  add_manifest(
    SCRIPT,
    old,
    file.path(out_dir, paste0(stem, ".png")),
    SRC(run, "spike_cells.csv"),
    note = sprintf(
      "axis window %g to %g arcsinh from the data",
      axis_lim[1],
      axis_lim[2]
    )
  )
})
