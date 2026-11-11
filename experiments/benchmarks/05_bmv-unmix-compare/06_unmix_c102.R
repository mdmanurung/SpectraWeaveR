# Purpose: unmix the 5 raw BMV files with AutoSpectral 1.0.2 code and a chosen spectra source.
# Run:     pixi run Rscript 06_unmix_c102.R <cell> <out_dir>
#            cell = s102 (reproduction check: 1.0.2 spectra, AF, variants) or s188_<stamp> (1.8.8 spectra)
#          from the repo root. Uses the 1.0.2 library installed in data/lib_as102.
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2 and autoSpectral_v1.8.8 [user zip]
# Outputs: <out_dir>/*.fcs unmixed files, <out_dir>/args.txt (what was run), <out_dir>/sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
LIB_102 <- file.path(DATA, "lib_as102") # AutoSpectral and AutoSpectralRcpp v1.0.2 (separate library)
WEIGHTED <- FALSE # as in ReRunUnmix.R (inherited: user zip)
# 1.8.8 AF file stamps follow the spectra stamps by order (04_spectra.R)
AF_STAMP_188 <- c("20261105_142441" = "20261105_142610", "20261105_143455" = "20261105_143604", "20261105_143803" = "20261105_143933")

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
cell <- args[[1]]
out_dir <- args[[2]]
stopifnot("out dir must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
.libPaths(c(LIB_102, .libPaths()))
suppressPackageStartupMessages({
  library(AutoSpectral)
  library(AutoSpectralRcpp)
})
stopifnot("1.0.2 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.0.2")

# Cell definition -------------------------------------------------------------------------
d102 <- file.path(DATA, "autoSpectral_v1.0.2")
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
if (cell == "s102") {
  control_dir <- file.path(d102, "Reference Group")
  control_file <- file.path(control_dir, "fcs_control_file.csv")
  spectra_file <- "Cells Matching Negative_autospectral_spectra.csv"
  spectra_dir <- file.path(control_dir, "table_spectra")
  af_file <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv"
  variants_rds <- file.path(control_dir, "figure_spectral_variants", "Spectral_variants.rds")
} else {
  stamp <- sub("^s188_", "", cell)
  control_dir <- file.path(d188, "Reference Group Clean")
  control_file <- file.path(control_dir, "fcs_control_file.csv")
  spectra_dir <- file.path(d188, "table_spectra")
  spectra_file <- paste0("Automated_autospectral_spectra_", stamp, ".csv")
  af_file <- paste0("H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra ", AF_STAMP_188[[stamp]], ".csv")
  variants_rds <- file.path(d188, "figure_spectral_variants", "Spectral_variants.rds")
}

asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(control_dir, control_file, asp)
spectra <- read.spectra(spectra_file, spectra.dir = spectra_dir)
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in ReRunUnmix.R
af.spectra <- read.spectra(af_file, spectra.dir = spectra_dir)
spectra.variants <- readRDS(variants_rds)

# Unmix -----------------------------------------------------------------------------------
writeLines(c(paste("cell:", cell), paste("spectra:", file.path(spectra_dir, spectra_file)), paste("af:", file.path(spectra_dir, af_file)), paste("variants:", variants_rds), paste("controls:", control_dir), paste("weighted:", WEIGHTED)), file.path(out_dir, "args.txt"))
unmix.folder(
  fcs.dir = file.path(d102, "Raw"),
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af.spectra,
  weighted = WEIGHTED,
  spectra.variants = spectra.variants,
  speed = "slow", # always slow (user, 2026-11-10); the 1.0.2 default
  output.dir = out_dir,
  parallel = TRUE,
  threads = NULL
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
