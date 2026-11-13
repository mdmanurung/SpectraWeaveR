# Purpose: copy of 06_unmix_c102.R that loads the AF csv RAW instead of through read.spectra().
#          Why: read.spectra() rescales each row to max 1, but AF row 1 of the 1.8.8 AF csv peaks
#          at 0.9828, so 06_unmix_c102.R used it 1.7% brighter than the user's original run.
#          Raw loading reproduced the user's 1.8.8 A3 output to float noise (verified 2026-11-10).
#          This script measures whether that offset changes the s188_c102 arm. Only the AF
#          load line and the out dir differ from 06_unmix_c102.R.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/06b_unmix_c102_rawaf.R <cell> <out_dir> [condition]
#            cell = s102 (reproduction check: 1.0.2 spectra, AF, variants) or s188 / s188_<stamp> (the config
#              1.8.8 spectra; the stamp must be the one in config stamps$as188_spectra)
#            condition = last argument, required when the dataset has conditions (24H: one task per condition,
#              out_dir = <data_dir>/unmix_cells/<arm>/<condition>), refused otherwise
#          from the repo root. Uses the 1.0.2 library from config lib_dirs (lib_as102).
#          SW_RUN_DATE is only read by scripts that write a run folder; this one writes to <out_dir>.
# Inputs:  <data_dir>/autoSpectral_v1.0.2 and autoSpectral_v1.8.8 (Raw[/<condition>], spectra, AF, variants, controls) [user zip]
#          stamp file names (spectra, AF) from config stamps; lib dir from config lib_dirs
# Outputs: <out_dir>/*.fcs unmixed files, <out_dir>/args.txt (what was run), <out_dir>/sessionInfo.txt

# Dataset config --------------------------------------------------------------------------
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath(data_path())
LIB_102 <- normalizePath(cfg("lib_dirs")[
  basename(cfg("lib_dirs")) == "lib_as102"
]) # AutoSpectral and AutoSpectralRcpp v1.0.2 (separate library)
WEIGHTED <- FALSE # as in ReRunUnmix.R (inherited: user zip)
STAMP_188 <- sub(
  "^Automated_autospectral_spectra_(.*)\\.csv$",
  "\\1",
  cfg("stamps")$as188_spectra
) # the only 1.8.8 run this script can load

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
cell <- args[[1]]
out_dir <- args[[2]]
# 24H: one task per condition subfolder (argument 3); 3D and EXV have no condition level
cond <- if (length(args) >= 3) args[[3]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "06b_unmix_c102_rawaf.R: dataset ",
    DS$dataset,
    " has no conditions; drop argument 3 ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "06b_unmix_c102_rawaf.R: needs argument 3, a condition of dataset ",
    DS$dataset,
    ": ",
    paste(DS$conditions, collapse = ", "),
    call. = FALSE
  )
}
stopifnot("out dir must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
.libPaths(c(LIB_102, .libPaths()))
suppressPackageStartupMessages({
  library(AutoSpectral)
  library(AutoSpectralRcpp)
})
stopifnot(
  "1.0.2 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.0.2"
)

# Cell definition -------------------------------------------------------------------------
d102 <- file.path(DATA, "autoSpectral_v1.0.2")
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
if (cell == "s102") {
  control_dir <- file.path(d102, "Reference Group")
  control_file <- file.path(control_dir, "fcs_control_file.csv")
  spectra_file <- cfg("stamps")$as102_spectra
  spectra_dir <- file.path(control_dir, "table_spectra")
  af_file <- cfg("stamps")$as102_af
  variants_rds <- file.path(
    control_dir,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  )
} else {
  stopifnot(
    "cell is s102, s188 or s188_<stamp> of config stamps$as188_spectra" = cell %in%
      c("s188", paste0("s188_", STAMP_188))
  )
  control_dir <- file.path(d188, "Reference Group Clean")
  control_file <- file.path(control_dir, "fcs_control_file.csv")
  spectra_dir <- file.path(d188, "table_spectra")
  spectra_file <- cfg("stamps")$as188_spectra
  af_file <- cfg("stamps")$as188_af
  variants_rds <- file.path(
    d188,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  )
}

asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(control_dir, control_file, asp)
spectra <- read.spectra(spectra_file, spectra.dir = spectra_dir)
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in ReRunUnmix.R
af.spectra <- as.matrix(read.csv(
  file.path(spectra_dir, af_file),
  row.names = 1,
  check.names = FALSE
)) # raw, no row rescale
spectra.variants <- readRDS(variants_rds)

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("cell:", cell),
    paste("spectra:", file.path(spectra_dir, spectra_file)),
    paste("af:", file.path(spectra_dir, af_file)),
    paste("variants:", variants_rds),
    paste("controls:", control_dir),
    paste("weighted:", WEIGHTED)
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
  weighted = WEIGHTED,
  spectra.variants = spectra.variants,
  speed = "slow", # always slow (user, 2026-11-10); the 1.0.2 default
  output.dir = out_dir,
  parallel = TRUE,
  threads = NULL
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
