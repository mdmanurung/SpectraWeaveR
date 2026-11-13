# Purpose: AutoSpectral 1.8.8 in its documented best form: AF spectra rebuilt with the package default
#          refine = TRUE, joint pipeline, with a chosen n.passes.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript 07_best_c188.R <n_passes> <out_dir> [af_rds|NA] [condition]
#            (from the repo root)
#            af_rds: reuse AF spectra built by an earlier run, so arms differ only in n.passes; NA = rebuild
#            condition = last argument, required when the dataset has conditions (24H: one task per condition,
#              out_dir = <data_dir>/unmix_cells/<arm>/<condition>, so af_rds needs "NA" or the earlier condition's rds),
#              refused otherwise
#            SW_RUN_DATE is only read by scripts that write a run folder; this one writes to <out_dir>.
# Inputs:  <data_dir>/autoSpectral_v1.8.8 (Reference Group Clean controls, config-stamp spectra,
#          Spectral_variants.rds built with defaults, unstained_AF, Raw[/<condition>]) [user zip]
#          stamp file name (spectra) from config stamps
# Outputs: <out_dir>/*.fcs, <out_dir>/af_spectra.rds, <out_dir>/args.txt, <out_dir>/sessionInfo.txt

# Dataset config --------------------------------------------------------------------------
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath(data_path())
AF_REFINE <- TRUE # package default in 1.8.8 (get_af_spectra.R:242); the original run used FALSE
PIPELINE <- "joint"

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
n_passes <- as.integer(args[[1]])
out_dir <- args[[2]]
af_rds <- if (length(args) >= 3 && args[[3]] != "NA") args[[3]] else NA
# 24H: one task per condition subfolder (argument 4); 3D and EXV have no condition level
cond <- if (length(args) >= 4) args[[4]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "07_best_c188.R: dataset ",
    DS$dataset,
    " has no conditions; drop argument 4 ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "07_best_c188.R: needs argument 4, a condition of dataset ",
    DS$dataset,
    ": ",
    paste(DS$conditions, collapse = ", "),
    call. = FALSE
  )
}
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
  cfg("stamps")$as188_spectra,
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
  fcs.dir = paste(c(file.path(d188, "Raw"), cond), collapse = "/"),
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
