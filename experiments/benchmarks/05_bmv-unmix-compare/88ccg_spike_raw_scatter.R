# Copy of 88c_spike_raw_scatter.R for the common SpectroFlo labels (plans/common_gates_plan.md G5, G10b). Only the inputs (the score tables of 88cg and 88bcg)
#          and the output dirs differ; method, levels, seed 42, loss cut and setups are inherited unchanged.
# Purpose: raw per-cell view behind the spike-in summary figures (user, 2026-11-10: "for each summary
#          figure show raw data like ggcyto geom hex"). For each spiked marker and T cell population, the
#          +1 arcsinh cells: marker value in the unspiked unmix (x) against the spiked unmix (y), ggcyto
#          hexbin with ggcyto's default pseudocolor. Dashed curve = full recovery of the spike; grey
#          diagonal = spike fully lost. Rows donors, columns setups; % of cells losing >50% per panel.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_spike_cg.sbatch 88ccg_spike_raw_scatter.R
#          (repo root). SW_DATASET and SW_RUN_DATE are required; the dataset sets the run folders and lib dirs.
# Inputs:  run_in("spike-in-cg")/spike_cells.csv         TNF, Th2, IL-10; 5 setups [88cg]
#          run_in("spike-in-panel2-cg")/spike_cells.csv  9 markers; 2 setups [88bcg]
#          (registry keys of datasets/<ds>_runs.csv; the folders are written by 88cg and 88bcg)
# Outputs: run_out("spike-raw-scatter-cg") = experiments/runs/<SW_RUN_DATE>_bmv<run_tag>-spike-raw-scatter-cg/
#          <out_dir>/fig_raw_<marker>_<pop>.png, .pdf   one figure per marker and population
#          <out_dir>/pct_lost.csv                       % of level-1 cells losing >50%, per panel
#          <out_dir>/sessionInfo.txt

# the dataset config comes first: the lib dir it names must be on .libPaths before ggcyto loads
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
filter <- dplyr::filter # dataset.R calls filter() from the search path; flowCore, attached later, would mask it
.libPaths(c(grep("lib_ggplot2_35", cfg("lib_dirs"), value = TRUE), .libPaths())) # ggplot2 3.5.2 for ggcyto
suppressPackageStartupMessages({
  library(flowCore)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
})

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- run_out("spike-raw-scatter-cg")
SOURCES <- file.path(
  c(run_in("spike-in-cg"), run_in("spike-in-panel2-cg")),
  "spike_cells.csv"
)
SETUP_LABELS <- c(
  c188 = "1.8.8 zip",
  afrefine = "refined AF",
  opt188 = "opt188",
  s188_c102 = "1.0.2 code,\n1.8.8 spectra",
  s102_c102 = "1.0.2 zip"
)
LEVEL <- 1 # the +1 arcsinh cells, as the summary figures (user, 2026-11-10)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
LOSS_FRAC <- 0.5 # "losing most of the spike", as 88 and 64
AXIS <- c(-1.5, 4) # display window only
BINS <- 60 # hexagons per panel; display only

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Level-1 cells, both spike runs, on arcsinh ------------------------------------------------------------
cells <- map(SOURCES, \(f) read_csv(f, show_col_types = FALSE)) |>
  list_rbind() |>
  dplyr::filter(level == LEVEL) |>
  mutate(
    before = asinh(base_value / COFACTOR),
    after = asinh((base_value + recovered) / COFACTOR),
    lost = recovery < LOSS_FRAC,
    setup_label = factor(SETUP_LABELS[setup], levels = SETUP_LABELS)
  )
pct_lost <- cells |>
  summarise(
    pct_lost = 100 * mean(lost),
    n = n(),
    .by = c(marker, pop, sample, setup, setup_label)
  )
write_csv(select(pct_lost, -setup_label), file.path(OUT_DIR, "pct_lost.csv"))

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

# One ggcyto figure per marker and population ------------------------------------------------------------
raw_figure <- \(which_marker, which_pop) {
  d <- dplyr::filter(cells, marker == which_marker, pop == which_pop)
  frames <- d |>
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
      setup_label = factor(
        setup_label,
        levels = intersect(SETUP_LABELS, unique(d$setup_label))
      )
    )
  pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)
  labels <- dplyr::filter(pct_lost, marker == which_marker, pop == which_pop)
  p <- ggcyto(fs, aes(x = before, y = after)) +
    geom_hex(bins = BINS) +
    geom_line(
      data = curves,
      aes(x, y, linetype = curve),
      inherit.aes = FALSE,
      colour = "grey30",
      linewidth = 0.4
    ) +
    geom_label(
      data = labels,
      aes(x = -Inf, y = Inf, label = sprintf("lost >50%%: %.1f%%", pct_lost)),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 7,
      size.unit = "pt",
      fontface = "bold",
      fill = "white",
      alpha = 0.85,
      label.size = 0
    ) +
    facet_grid(sample ~ setup_label, switch = "y") +
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
      x = sprintf("%s, unspiked unmix (arcsinh, cofactor 6000)", which_marker),
      y = sprintf("%s, +1 arcsinh spiked unmix", which_marker),
      title = sprintf("%s: %s spike, per cell", which_pop, which_marker)
    ) +
    ggpubr::theme_pubr(base_size = 8, border = TRUE) +
    theme(
      strip.text.x = element_text(size = 7),
      strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(1, "mm")
    )
  n_setups <- nlevels(pd$setup_label)
  stem <- file.path(OUT_DIR, sprintf("fig_raw_%s_%s", which_marker, which_pop))
  ggsave(
    paste0(stem, ".pdf"),
    as.ggplot(p),
    width = 30 * n_setups + 25,
    height = 175,
    units = "mm",
    device = cairo_pdf
  )
  ggsave(
    paste0(stem, ".png"),
    as.ggplot(p),
    width = 30 * n_setups + 25,
    height = 175,
    units = "mm",
    dpi = 200,
    bg = "white"
  )
}
distinct(cells, marker, pop) |>
  pwalk(\(marker, pop) raw_figure(marker, pop))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
