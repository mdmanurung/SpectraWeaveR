#!/usr/bin/env Rscript
# Purpose: run the AutoSpectral workflow through SpectraWeaveR on one dataset, then
#          write diagnostics, a report and figures F1 to F4. Shared by all benchmarks.
# Run:     not called directly; each benchmark's run.R sets CONFIG_PATH and sources this file.
#          sbatch experiments/slurm/run_script.sbatch experiments/benchmarks/<NN_name>/run.R experiments/runs/<date>_<slug>
# Inputs:  the dataset config (config.R of the benchmark) and the data it names  [experiments/datasets/download.sh]
# Outputs (in <run_dir>):
#   unmix/                   AutoSpectral outputs: spectra, AF, variants, control tables
#   unmixed/<method>/        unmixed FCS: ols, autospectral_fast, autospectral_slow
#   report/                  report.html, checks.csv, spectra.csv, provenance.csv, files.csv
#   report/dashboard.html    interactive dashboard over the report tables and figure source data
#   figures/                 F1 to F4 as png/pdf, with source_data/<panel>.csv
#   spectra_compare.csv      angle between our spectra and spectracle's, per fluorophore
#   sessionInfo.txt          R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE) # exported functions only

# Parameters (from the benchmark's config.R) -----------------------------------
cfg <- source(CONFIG_PATH, local = TRUE)$value
SEED <- cfg$seed
N_EVENTS <- cfg$n_events
PSEUDO_LOG_SIGMA <- cfg$pseudo_log_sigma
BIPLOT_BINS <- cfg$biplot_bins
BIPLOT_PAIR <- cfg$biplot_pair
CYTOMETER <- cfg$cytometer
stopifnot("config names its diagnostic limits" = is.list(cfg$limits))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run dir must not exist yet" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)
set.seed(SEED)

data_dir <- Sys.getenv("SW_DATA_DIR", "/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata")
root <- file.path(data_dir, cfg$dataset)
data_path <- \(rel) if (is.null(rel)) NULL else file.path(root, rel)
control_dir <- data_path(cfg$control_dir)
control_file <- if (!is.null(cfg$control_file_repo)) cfg$control_file_repo else data_path(cfg$control_file) # NULL: AutoSpectral builds and checks one
unstained <- data_path(cfg$unstained)
stained <- data_path(cfg$stained)
spectroflo <- data_path(cfg$spectroflo)
stopifnot(
  "inputs exist" = all(file.exists(c(control_file, unstained, stained, spectroflo))),
  "control dir exists" = dir.exists(control_dir)
)

# Optional spectra set, for one control file that holds cell and bead controls (AutoSpectral article 02b).
# Bead rows are named "<dye> beads". The set file keeps the cell rows ("cells"), or the bead row where one
# exists and the cell row otherwise ("beads"), and names every row by its dye so that both sets give the
# same channel names. Each row is extracted on its own, so a set gives the spectra a subset of one run would.
if (!is.null(cfg$spectra_set)) {
  stopifnot("spectra_set is cells or beads" = cfg$spectra_set %in% c("cells", "beads"))
  combined <- read.csv(control_file, check.names = FALSE, stringsAsFactors = FALSE)
  is_bead <- grepl(" beads$", combined$fluorophore)
  dye <- sub(" beads$", "", combined$fluorophore)
  keep <- if (cfg$spectra_set == "cells") {
    !is_bead & combined$control.type == "cells"
  } else {
    is_bead | (combined$control.type == "cells" & !dye %in% dye[is_bead])
  }
  if (cfg$spectra_set == "beads") keep <- keep | combined$fluorophore == "Negative"
  set_file <- combined[keep, , drop = FALSE]
  set_file$fluorophore <- sub(" beads$", "", set_file$fluorophore)
  stopifnot("one row per fluorophore in the set file" = !anyDuplicated(set_file$fluorophore))
  control_file <- file.path(run_dir, "control_file_set.csv")
  write.csv(set_file, control_file, row.names = FALSE)
}

# Setup and reference spectra --------------------------------------------------
setup <- sw_unmix_setup(
  control_dir = control_dir,
  cytometer = CYTOMETER,
  control_file = control_file,
  output_dir = file.path(run_dir, "unmix")
)
controls <- sw_unmix_prepare(setup)

# Autofluorescence and variants (AutoSpectral defaults) -------------------------
af_spectra <- sw_unmix_extract_af(unstained, setup, controls$spectra)
variants <- sw_unmix_extract_variants(setup, controls$spectra)

# Unmix the stained lung file with each method ---------------------------------
unmix_one <- function(method, speed = NULL) {
  label <- if (is.null(speed)) method else paste0(method, "_", speed)
  out <- file.path(run_dir, "unmixed", tolower(label))
  dir.create(out, recursive = TRUE)
  fn_args <- list(
    input = stained,
    spectra = controls$spectra,
    setup = setup,
    flow_control = controls$flow_control,
    method = method,
    output_dir = out,
    file_suffix = label
  )
  if (method == "AutoSpectral") {
    fn_args <- c(fn_args, list(af_spectra = af_spectra, spectra_variants = variants, speed = speed))
  }
  do.call(sw_unmix_run, fn_args)
}
run_ols <- unmix_one("OLS")
run_fast <- unmix_one("AutoSpectral", "fast")
run_slow <- unmix_one("AutoSpectral", "slow")
unmixed_files <- c(
  ols = run_ols$output_files,
  autospectral_fast = run_fast$output_files,
  autospectral_slow = run_slow$output_files
)

# Diagnostics and report (tables; verdicts only where a limit is given) ---------
result <- sw_unmix_result(
  spectra = controls$spectra,
  af_spectra = af_spectra,
  spectra_variants = variants,
  unmixed_files = unname(unmixed_files),
  setup = setup,
  flow_control = controls$flow_control,
  method = "OLS, AutoSpectral fast, AutoSpectral slow",
  args = list(workflow = "automated", dataset = cfg$dataset, methods = cfg$methods),
  seed = SEED,
  inputs = c(control_file, unstained, stained, spectroflo)
)
diagnostics <- sw_unmix_diagnose(result, limits = cfg$limits, unmixed_files = unname(unmixed_files))
sw_unmix_report(
  result,
  diagnostics,
  path = file.path(run_dir, "report"),
  title = cfg$title,
  notes = cfg$notes
)

# Comparison figures and tables (reads the outputs written above) ---------------
source("experiments/benchmarks/_common/render_figures.R")

# Interactive dashboard over the report tables and the figure source data -------
# Same pseudo-log sigma and bin count as F3 and F4, so the binned panels match the figures.
sw_unmix_dashboard(
  stats::setNames(run_dir, basename(run_dir)),
  path = file.path(run_dir, "report", "dashboard.html"),
  pseudo_log_sigma = PSEUDO_LOG_SIGMA,
  bins = BIPLOT_BINS,
  title = cfg$title,
  figures_dir = out_name,
  axis_labels = BIPLOT_PAIR
)
