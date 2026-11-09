#!/usr/bin/env Rscript
# Purpose: test AutoSpectral's unmixing-error correction through SpectraWeaveR on the Mendeley lung file.
# Run:     sbatch experiments/slurm/run_script.sbatch experiments/benchmarks/04_mendeley-correction/run.R experiments/runs/<date>_mendeley-correction
# Inputs:  settings from experiments/benchmarks/01_mendeley-7color/config.R and the data it names
# Outputs (in <run_dir>):
#   unmix/                  AutoSpectral outputs
#   unmixed/<run>/          ols_original, ols_corrected, autospectral_fast_original, autospectral_fast_corrected
#   correction_changes.csv  per fluorophore: angle between original and corrected spectrum
#   report/                 report for the corrected spectra
#   figures/                F1 (original vs corrected), F3, F4 with source_data/
#   sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters (benchmark 01's config; nothing new is chosen here) ----------------
cfg <- source("experiments/benchmarks/01_mendeley-7color/config.R", local = TRUE)$value
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run dir must not exist yet" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)
set.seed(cfg$seed)

root <- file.path(Sys.getenv("SW_DATA_DIR", "/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata"), cfg$dataset)
control_dir <- file.path(root, cfg$control_dir)
control_file <- file.path(root, cfg$control_file)
unstained <- file.path(root, cfg$unstained)
stained <- file.path(root, cfg$stained)
spectroflo <- file.path(root, cfg$spectroflo)
stopifnot(all(file.exists(c(control_file, unstained, stained, spectroflo))))

# Baseline: spectra, AF and variants exactly as in benchmark 01 -----------------
setup <- sw_unmix_setup(control_dir, cytometer = cfg$cytometer, control_file = control_file, output_dir = file.path(run_dir, "unmix"))
controls <- sw_unmix_prepare(setup)
af_spectra <- sw_unmix_extract_af(unstained, setup, controls$spectra)
variants <- sw_unmix_extract_variants(setup, controls$spectra)

# Correction: signatures, then spillover (AutoSpectral defaults for all tuning) --
correction <- sw_unmix_correct(
  stained = stained,
  spectra = controls$spectra,
  asp = setup$asp,
  steps = c("signatures", "spillover"),
  signatures_args = list(
    unstained.sample = unstained,
    flow.control = controls$flow_control,
    variants = variants
  ),
  spillover_args = list(
    unstained.sample = unstained,
    flow.control = controls$flow_control,
    variants = variants,
    af.spectra = af_spectra
  )
)
write.csv(correction$changes, file.path(run_dir, "correction_changes.csv"), row.names = FALSE)
variants_corrected <- sw_unmix_extract_variants(setup, correction$spectra)

# Unmix the stained file with original and corrected spectra --------------------
unmix_one <- function(label, spectra, method, variants = NULL) {
  out <- file.path(run_dir, "unmixed", label)
  dir.create(out, recursive = TRUE)
  args <- list(input = stained, spectra = spectra, setup = setup, flow_control = controls$flow_control,
               method = method, output_dir = out, file_suffix = label)
  if (method == "AutoSpectral") {
    args <- c(args, list(af_spectra = af_spectra, spectra_variants = variants, speed = "fast"))
  }
  do.call(sw_unmix_run, args)$output_files
}
unmixed <- c(
  ols = unmix_one("ols", controls$spectra, "OLS"),
  ols_corrected = unmix_one("ols_corrected", correction$spectra, "OLS"),
  autospectral_fast = unmix_one("autospectral_fast", controls$spectra, "AutoSpectral", variants),
  autospectral_fast_corrected = unmix_one("autospectral_fast_corrected", correction$spectra, "AutoSpectral", variants_corrected)
)

# Report for the corrected spectra (tables; no limits supplied) -----------------
result <- sw_unmix_result(
  spectra = correction$spectra, af_spectra = af_spectra, spectra_variants = variants_corrected,
  unmixed_files = unname(unmixed), setup = setup, flow_control = controls$flow_control,
  method = "OLS and AutoSpectral fast, original vs corrected spectra",
  args = list(steps = c("signatures", "spillover")), seed = cfg$seed,
  inputs = c(control_file, unstained, stained, spectroflo)
)
diagnostics <- sw_unmix_diagnose(result, limits = list(), unmixed_files = unname(unmixed))
sw_unmix_report(
  result, diagnostics, path = file.path(run_dir, "report"),
  title = "Mendeley 7-colour lung: unmixing-error correction",
  notes = "Candidate correction judged against unmodified AutoSpectral; no improvement is claimed. No pass/fail limits were set."
)

# Figures F1, F3, F4 -----------------------------------------------------------
read_exprs <- function(path) flowCore::exprs(flowCore::read.FCS(path, truncate_max_range = FALSE))
event_list <- c(lapply(unmixed, read_exprs), list(spectroflo = read_exprs(spectroflo)))
names(event_list) <- c(names(unmixed), "spectroflo")
common <- Reduce(intersect, lapply(event_list, colnames))
common <- common[grepl("-A$", common) & !grepl("^(FSC|SSC|Time)", common)]
events <- sw_unmix_sample_events(lapply(event_list, \(m) m[, common, drop = FALSE]), n_events = cfg$n_events, seed = cfg$seed)
panels <- list(
  F1 = sw_unmix_figure_spectra(list(original = controls$spectra, corrected = correction$spectra)),
  F3 = sw_unmix_figure_distributions(events, channels = common, pseudo_log_sigma = cfg$pseudo_log_sigma),
  F4 = sw_unmix_figure_biplot(events, x_channel = cfg$biplot_pair[[1]], y_channel = cfg$biplot_pair[[2]],
                              pseudo_log_sigma = cfg$pseudo_log_sigma, bins = cfg$biplot_bins)
)
sw_unmix_save_figures(panels, file.path(run_dir, "figures"))

writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
