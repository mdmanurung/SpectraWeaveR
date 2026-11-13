# Purpose: AutoSpectral 1.0.2 in its documented best form: AF spectra and spectral variants rebuilt with
#          refine = TRUE, then unmixing at speed = "slow" (10 variants, the 1.0.2 default).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript 07_best_c102.R <spectra: neg | af> <out_dir> [condition]
#            (from the repo root)
#            neg = "Cells Matching Negative" spectra (used in the original run), af = "Cells AF clean up" spectra
#            condition = last argument, required when the dataset has conditions (24H: one task per condition,
#              out_dir = <data_dir>/unmix_cells/<arm>/<condition>), refused otherwise
#            SW_RUN_DATE is only read by scripts that write a run folder; this one writes to <out_dir>.
# Inputs:  <data_dir>/autoSpectral_v1.0.2 (controls, spectra, unstained_AF, Raw[/<condition>]) [user zip]
#          spectra file name for neg from config stamps$as102_spectra; lib dir from config lib_dirs
# Outputs: <out_dir>/*.fcs, <out_dir>/af_spectra.rds, <out_dir>/spectral_variants.rds, <out_dir>/args.txt, sessionInfo.txt

# Dataset config --------------------------------------------------------------------------
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath(data_path())
LIB_102 <- normalizePath(cfg("lib_dirs")[
  basename(cfg("lib_dirs")) == "lib_as102"
])
SPECTRA_FILES <- c(
  neg = cfg("stamps")$as102_spectra,
  af = "Cells AF clean up_autospectral_spectra.csv" # a 1.0.2 setup file name, not a stamp
)
AF_REFINE <- TRUE # grill Q16: documented quality lever
VARIANTS_REFINE <- TRUE # grill Q16: "may help reduce spillover spread", no unmixing-time cost

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
spectra_choice <- args[[1]]
out_dir <- args[[2]]
# 24H: one task per condition subfolder (argument 3); 3D and EXV have no condition level
cond <- if (length(args) >= 3) args[[3]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "07_best_c102.R: dataset ",
    DS$dataset,
    " has no conditions; drop argument 3 ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "07_best_c102.R: needs argument 3, a condition of dataset ",
    DS$dataset,
    ": ",
    paste(DS$conditions, collapse = ", "),
    call. = FALSE
  )
}
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
  fcs.dir = paste(c(file.path(d102, "Raw"), cond), collapse = "/"),
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
