# Purpose: unmix the raw BMV files of the dataset with AutoSpectral 1.8.8 code and the 1.0.2 spectra, AF and variants.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript 06_unmix_c188.R <cell> <speed> <out_dir> [condition]
#            (from the repo root; the pixi env has AutoSpectral 1.8.8)
#            cell = s102 (1.0.2 spectra, AF, variants) or s188 (1.8.8 config-stamp spectra, AF, variants)
#            speed = slow only (user, 2026-11-10); 10 variants, the 1.0.2 default
#            condition = last argument, required when the dataset has conditions (24H: one task per condition,
#              out_dir = <data_dir>/unmix_cells/<arm>/<condition>), refused otherwise
#            SW_RUN_DATE is only read by scripts that write a run folder; this one writes to <out_dir>.
# Inputs:  <data_dir>/autoSpectral_v1.0.2 or autoSpectral_v1.8.8 (spectra, AF, variants, controls, Raw[/<condition>]) [user zip]
#          stamp file names (spectra, AF) from config stamps; samples from datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <out_dir>/*.fcs unmixed files, <out_dir>/args.txt (what was run), <out_dir>/sessionInfo.txt

# Dataset config --------------------------------------------------------------------------
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath(data_path())
PIPELINE <- "joint" # as in AutoSpectral_unmix_v1.8.8.R (inherited: user zip); other 1.8.8 arguments are its defaults, as there

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
cell <- args[[1]]
speed <- args[[2]]
out_dir <- args[[3]]
stopifnot(
  "speed is slow (user, 2026-11-10: always slow)" = speed == "slow"
)
# 24H: one task per condition subfolder (4th argument); 3D and EXV have no condition level
cond <- if (length(args) >= 4) args[[4]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "06_unmix_c188.R: dataset ",
    DS$dataset,
    " has no conditions; drop the 4th argument ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "06_unmix_c188.R: needs a 4th argument, a condition of dataset ",
    DS$dataset,
    ": ",
    paste(DS$conditions, collapse = ", "),
    call. = FALSE
  )
}
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
  spectra_file <- cfg("stamps")$as102_spectra
  af_file <- cfg("stamps")$as102_af
  variants_rds <- file.path(
    control_dir,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  )
  raw_dir <- paste(c(file.path(d102, "Raw"), cond), collapse = "/")
} else {
  stopifnot(cell == "s188")
  control_dir <- file.path(d188, "Reference Group Clean")
  spectra_dir <- file.path(d188, "table_spectra")
  spectra_file <- cfg("stamps")$as188_spectra # 3D: run 143803, variants file at 14:40, unmix at 14:41
  af_file <- cfg("stamps")$as188_af
  variants_rds <- file.path(
    d188,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  )
  raw_dir <- paste(c(file.path(d188, "Raw"), cond), collapse = "/")
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
