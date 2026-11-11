# Purpose: unmix the 5 raw BMV files with AutoSpectral 1.8.8 code and the 1.0.2 spectra, AF and variants.
# Run:     pixi run Rscript 06_unmix_c188.R <cell> <speed> <out_dir>   (from the repo root; the pixi env has AutoSpectral 1.8.8)
#            cell = s102 (1.0.2 spectra, AF, variants) or s188 (1.8.8 run 143803 spectra, AF, variants)
#            speed = slow only (user, 2026-11-10); 10 variants, the 1.0.2 default
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2 or autoSpectral_v1.8.8 (spectra, AF, variants, controls, raw) [user zip]
# Outputs: <out_dir>/*.fcs unmixed files, <out_dir>/args.txt (what was run), <out_dir>/sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
PIPELINE <- "joint" # as in AutoSpectral_unmix_v1.8.8.R (inherited: user zip); other 1.8.8 arguments are its defaults, as there

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
cell <- args[[1]]
speed <- args[[2]]
out_dir <- args[[3]]
stopifnot(
  "speed is slow (user, 2026-11-10: always slow)" = speed == "slow"
)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
suppressPackageStartupMessages(library(AutoSpectral))
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)

# Inputs by cell --------------------------------------------------------------------------
d102 <- file.path(DATA, "autoSpectral_v1.0.2")
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
if (cell == "s102") {
  control_dir <- file.path(d102, "Reference Group")
  spectra_dir <- file.path(control_dir, "table_spectra")
  spectra_file <- "Cells Matching Negative_autospectral_spectra.csv"
  af_file <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv"
  variants_rds <- file.path(
    control_dir,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  )
  raw_dir <- file.path(d102, "Raw")
} else {
  stopifnot(cell == "s188")
  control_dir <- file.path(d188, "Reference Group Clean")
  spectra_dir <- file.path(d188, "table_spectra")
  spectra_file <- "Automated_autospectral_spectra_20261105_143803.csv" # run 143803: variants file at 14:40, unmix at 14:41
  af_file <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
  variants_rds <- file.path(
    d188,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  )
  raw_dir <- file.path(d188, "Raw")
}
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(
  control_dir,
  file.path(control_dir, "fcs_control_file.csv"),
  asp
)
spectra <- read.spectra(spectra_file, spectra.dir = spectra_dir)
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in ReRunUnmix.R; no AF row in the 1.8.8 file
af.spectra <- read.spectra(af_file, spectra.dir = spectra_dir)
spectra.variants <- readRDS(variants_rds)

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("cell:", cell),
    paste("speed:", speed),
    paste("pipeline:", PIPELINE),
    paste("spectra:", file.path(spectra_dir, spectra_file)),
    paste("af:", file.path(spectra_dir, af_file)),
    paste("variants:", variants_rds)
  ),
  file.path(out_dir, "args.txt")
)
unmix.folder(
  fcs.dir = raw_dir,
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af.spectra,
  spectra.variants = spectra.variants,
  pipeline = PIPELINE,
  speed = speed,
  output.dir = out_dir,
  parallel = TRUE,
  verbose = TRUE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
