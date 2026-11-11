# Purpose: AutoSpectral 1.0.2 in its documented best form: AF spectra and spectral variants rebuilt with
#          refine = TRUE, then unmixing at speed = "slow" (10 variants, the 1.0.2 default).
# Run:     pixi run Rscript 07_best_c102.R <spectra: neg | af> <out_dir>   (from the repo root)
#            neg = "Cells Matching Negative" spectra (used in the original run), af = "Cells AF clean up" spectra
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2 (controls, spectra, unstained_AF, Raw) [user zip]
# Outputs: <out_dir>/*.fcs, <out_dir>/af_spectra.rds, <out_dir>/spectral_variants.rds, <out_dir>/args.txt, sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
LIB_102 <- file.path(DATA, "lib_as102")
SPECTRA_FILES <- c(
  neg = "Cells Matching Negative_autospectral_spectra.csv",
  af = "Cells AF clean up_autospectral_spectra.csv"
)
AF_REFINE <- TRUE # grill Q16: documented quality lever
VARIANTS_REFINE <- TRUE # grill Q16: "may help reduce spillover spread", no unmixing-time cost

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
spectra_choice <- args[[1]]
out_dir <- args[[2]]
stopifnot(
  "out dir must not exist" = !dir.exists(out_dir),
  "spectra is neg or af" = spectra_choice %in% names(SPECTRA_FILES)
)
dir.create(out_dir, recursive = TRUE)
out_dir <- normalizePath(out_dir)
.libPaths(c(LIB_102, .libPaths()))
suppressPackageStartupMessages({
  library(AutoSpectral)
  library(AutoSpectralRcpp)
})
stopifnot(
  "1.0.2 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.0.2"
)

# Inputs ----------------------------------------------------------------------------------
d102 <- file.path(DATA, "autoSpectral_v1.0.2")
control_dir <- file.path(d102, "Reference Group")
setwd(control_dir) # 1.0.2 helpers read control files relative to the control folder, as in Autospectral_Setup.R
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(control_dir, "fcs_control_file.csv", asp)
spectra <- read.spectra(
  SPECTRA_FILES[[spectra_choice]],
  spectra.dir = file.path(control_dir, "table_spectra")
)
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE]
unstained <- list.files(
  file.path(control_dir, "unstained_AF"),
  pattern = "\\.fcs$",
  full.names = TRUE
)
stopifnot("one unstained AF file" = length(unstained) == 1)

# Rebuild AF spectra and variants with refine = TRUE ---------------------------------------
af.spectra <- get.af.spectra(
  unstained,
  asp,
  spectra,
  figures = FALSE,
  plot.dir = out_dir,
  table.dir = out_dir,
  refine = AF_REFINE
)
saveRDS(af.spectra, file.path(out_dir, "af_spectra.rds"))
spectra.variants <- get.spectral.variants(
  control_dir,
  "fcs_control_file.csv",
  asp,
  spectra,
  af.spectra,
  figures = FALSE,
  output.dir = out_dir,
  parallel = TRUE,
  refine = VARIANTS_REFINE
)
saveRDS(spectra.variants, file.path(out_dir, "spectral_variants.rds"))

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("spectra:", SPECTRA_FILES[[spectra_choice]]),
    paste("af refine:", AF_REFINE),
    paste("variants refine:", VARIANTS_REFINE),
    "speed: slow (1.0.2 default, 10 variants)",
    "weighted: FALSE"
  ),
  file.path(out_dir, "args.txt")
)
unmix.folder(
  fcs.dir = file.path(d102, "Raw"),
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af.spectra,
  weighted = FALSE,
  spectra.variants = spectra.variants,
  output.dir = out_dir,
  parallel = TRUE,
  threads = NULL
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
