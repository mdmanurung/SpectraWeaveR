# Purpose: AutoSpectral 1.8.8 in its documented best form: AF spectra rebuilt with the package default
#          refine = TRUE, joint pipeline, with a chosen n.passes.
# Run:     pixi run Rscript 07_best_c188.R <n_passes> <out_dir> [af_rds]   (from the repo root)
#            af_rds: reuse AF spectra built by an earlier run, so arms differ only in n.passes
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (Reference Group Clean controls, run-143803 spectra,
#          Spectral_variants.rds built with defaults, unstained_AF, Raw) [user zip]
# Outputs: <out_dir>/*.fcs, <out_dir>/af_spectra.rds, <out_dir>/args.txt, <out_dir>/sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
AF_REFINE <- TRUE # package default in 1.8.8 (get_af_spectra.R:242); the original run used FALSE
PIPELINE <- "joint"

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
n_passes <- as.integer(args[[1]])
out_dir <- args[[2]]
af_rds <- if (length(args) >= 3) args[[3]] else NA
stopifnot(
  "out dir must not exist" = !dir.exists(out_dir),
  "n_passes is a positive integer" = !is.na(n_passes) && n_passes >= 1
)
dir.create(out_dir, recursive = TRUE)
out_dir <- normalizePath(out_dir)
suppressPackageStartupMessages(library(AutoSpectral))
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)

# Inputs ----------------------------------------------------------------------------------
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
control_dir <- file.path(d188, "Reference Group Clean")
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(
  control_dir,
  file.path(control_dir, "fcs_control_file.csv"),
  asp
)
spectra <- read.spectra(
  "Automated_autospectral_spectra_20261105_143803.csv",
  spectra.dir = file.path(d188, "table_spectra")
)
spectra.variants <- readRDS(file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
))
unstained <- list.files(
  file.path(d188, "Reference Group", "unstained_AF"),
  pattern = "\\.fcs$",
  full.names = TRUE
)
stopifnot("one unstained AF file" = length(unstained) == 1)

# AF spectra: rebuilt with refine = TRUE, or reused ----------------------------------------
af.spectra <- if (is.na(af_rds)) {
  get.af.spectra(
    unstained,
    asp,
    spectra,
    figures = FALSE,
    plot.dir = out_dir,
    table.dir = out_dir,
    refine = AF_REFINE
  )
} else {
  readRDS(af_rds)
}
saveRDS(af.spectra, file.path(out_dir, "af_spectra.rds"))

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("af refine:", AF_REFINE),
    paste("af reused from:", af_rds),
    paste("pipeline:", PIPELINE),
    paste("n.passes:", n_passes)
  ),
  file.path(out_dir, "args.txt")
)
unmix.folder(
  fcs.dir = file.path(d188, "Raw"),
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af.spectra,
  spectra.variants = spectra.variants,
  pipeline = PIPELINE,
  n.passes = n_passes,
  speed = "slow", # always slow (user, 2026-11-10); no effect under pipeline "joint"
  output.dir = out_dir,
  parallel = TRUE,
  verbose = TRUE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
