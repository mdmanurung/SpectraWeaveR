# Purpose: render F1 to F4 from the unit-test fixtures for the before/after approval of the Nature
#          theme migration (plan v2 step 4, plot-first rows F1-F4 migrated; user, 2026-11-11).
# Run:     pixi run Rscript experiments/integration/fig_migration_render.R <pkg_dir> <label> <out_dir>
#            pkg_dir  package checkout to load (main for "before", the feature branch for "after")
#            label    "before" or "after"; appended to each file name
#            out_dir  run folder; files are never overwritten
# Inputs:  the package at <pkg_dir>; fixtures as in tests/testthat/test-unmix_figures.R
# Outputs: <out_dir>/F1_<label>.png .. F4_<label>.png   one PNG per panel at its saved size
#          <out_dir>/sizes_<label>.csv                    panel, width_mm, height_mm

suppressPackageStartupMessages({
  pkgload::load_all(
    commandArgs(trailingOnly = TRUE)[1],
    quiet = TRUE,
    export_all = FALSE
  )
})

# Parameters ---------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
LABEL <- args[2]
OUT_DIR <- args[3]

# Fixtures, as the unit tests build them --------------------------------------------------------
spectra <- rbind(
  FITC = c(0.1, 1, 0.5, 0.2, 0.1, 0),
  BB515 = c(0.1, 0.9, 0.6, 0.2, 0.1, 0),
  PE = c(0, 0.1, 0.3, 1, 0.6, 0.2)
)
colnames(spectra) <- paste0("D", 1:6)
set.seed(42)
events_list <- list(
  ours = matrix(
    rnorm(3000, sd = 50),
    ncol = 3,
    dimnames = list(NULL, c("A", "B", "C"))
  ),
  spectroflo = matrix(
    rnorm(1500, sd = 60),
    ncol = 3,
    dimnames = list(NULL, c("A", "B", "C"))
  )
)
events <- sw_unmix_sample_events(events_list, n_events = 1000, seed = 1)

# Build and save the four panels -----------------------------------------------------------------
figures <- list(
  F1 = sw_unmix_figure_spectra(list(
    ours = spectra,
    spectroflo = spectra * 0.9
  )),
  F2 = sw_unmix_figure_similarity(spectra),
  F3 = sw_unmix_figure_distributions(
    events,
    c("A", "B", "C"),
    pseudo_log_sigma = 10
  ),
  F4 = sw_unmix_figure_biplot(
    events,
    "A",
    "B",
    pseudo_log_sigma = 10,
    bins = 30
  )
)
tmp <- tempfile("figs")
saved <- sw_unmix_save_figures(figures, tmp)

# Copy each PNG under its panel and label ----------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
targets <- file.path(OUT_DIR, sprintf("%s_%s.png", saved$panel, LABEL))
stopifnot("figures must not exist" = !any(file.exists(targets)))
file.copy(saved$png, targets)
utils::write.csv(
  saved[c("panel", "width_mm", "height_mm")],
  file.path(OUT_DIR, sprintf("sizes_%s.csv", LABEL)),
  row.names = FALSE
)
