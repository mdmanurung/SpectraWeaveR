# Purpose: find which part of the AutoSpectral 1.0.2 spectra set makes BV510 (TNF) collapse to ~0 in many
#          T cells. Unmix the 5 raw BMV files with 1.8.8 code and the 1.8.8 dye spectra, swapping one
#          component of the 1.8.8 set for the 1.0.2 one: the AF spectra (swap_af102) or the variants
#          (swap_var102). Scoring is in 37b_score_swap.R.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/37_swap_af_variants.R <arm> <out_dir>
#            from the repo root. arm = swap_af102 or swap_var102. <out_dir> must not exist.
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (run 143803 dye spectra, AF, variants,
#            Reference Group Clean controls, Raw, AutoSpectral_unmixed) [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2 (AF csv, Spectral_variants.rds) [user zip]
# Outputs: <out_dir>/*.fcs      unmixed files
#          <out_dir>/args.txt   what was run
#          <out_dir>/verify.csv event count and channel names vs as188_joint, per file
#          <out_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
PIPELINE <- "joint" # as in 06_unmix_c188.R and AutoSpectral_unmix_v1.8.8.R (inherited: user zip)
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint" (repro188 = s188_c188_slow); was "fast" to match as188_joint

# Setup -----------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
arm <- cli_args[[1]]
out_dir <- cli_args[[2]]
stopifnot(
  "arm is swap_af102 or swap_var102" = arm %in% c("swap_af102", "swap_var102")
)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)

# Inputs: 1.8.8 set (the s188 cell of 06) and the 1.0.2 AF and variants (the s102 cell of 06) ----
d102 <- file.path(DATA, "autoSpectral_v1.0.2")
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
control_dir <- file.path(d188, "Reference Group Clean")
spectra_file <- file.path(
  d188,
  "table_spectra",
  "Automated_autospectral_spectra_20261105_143803.csv"
)
af188_file <- file.path(
  d188,
  "table_spectra",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
)
af102_file <- file.path(
  d102,
  "Reference Group",
  "table_spectra",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv"
)
var188_rds <- file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)
var102_rds <- file.path(
  d102,
  "Reference Group",
  "figure_spectral_variants",
  "Spectral_variants.rds"
)

# Raw csv read: read.spectra() rescales rows and caused an offset (user, 2026-11-10). The dye rows
# already have max 1 (checked), so the dye spectra equal 06's.
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra <- read_raw(spectra_file)
af188 <- read_raw(af188_file)
af102 <- read_raw(af102_file)
variants188 <- readRDS(var188_rds)
variants102 <- readRDS(var102_rds)

# Detector set and fluor set must match; stop, do not coerce
stopifnot(
  "AF 1.8.8 and 1.0.2 detectors equal dye spectra detectors" = identical(
    colnames(af188),
    colnames(spectra)
  ) &&
    identical(colnames(af102), colnames(spectra)),
  "variants 1.0.2 fluors equal dye spectra fluors" = identical(
    names(variants102$variants),
    rownames(spectra)
  ),
  "variants 1.8.8 fluors equal dye spectra fluors" = identical(
    names(variants188$variants),
    rownames(spectra)
  ),
  "variants 1.0.2 detectors equal dye spectra detectors" = all(map_lgl(
    variants102$variants,
    \(m) identical(colnames(m), colnames(spectra))
  ))
)
af_spectra <- if (arm == "swap_af102") af102 else af188
spectra_variants <- if (arm == "swap_var102") variants102 else variants188

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste("speed:", SPEED),
    paste("pipeline:", PIPELINE),
    paste("spectra:", spectra_file),
    paste("af:", if (arm == "swap_af102") af102_file else af188_file),
    paste("variants:", if (arm == "swap_var102") var102_rds else var188_rds),
    paste("controls:", control_dir),
    "weighted: default (as 06)",
    "csv read: raw as.matrix(read.csv(row.names = 1, check.names = FALSE))"
  ),
  file.path(out_dir, "args.txt")
)
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(
  control_dir,
  file.path(control_dir, "fcs_control_file.csv"),
  asp
)
t0 <- Sys.time()
unmix.folder(
  fcs.dir = file.path(d188, "Raw"),
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af_spectra,
  spectra.variants = spectra_variants,
  pipeline = PIPELINE,
  speed = SPEED,
  output.dir = out_dir,
  parallel = TRUE,
  verbose = TRUE
)
message(
  "unmix minutes: ",
  round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
)

# Verify: same files, event counts and channel names as as188_joint -----------------------------
ref_dir <- file.path(d188, "AutoSpectral_unmixed")
read_hdr <- \(path) {
  ff <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  list(n = nrow(ff), channels = colnames(ff))
}
verify <- tibble(file = list.files(ref_dir, pattern = "\\.fcs$")) |>
  mutate(
    ref = map(file.path(ref_dir, file), read_hdr),
    new = map(file.path(out_dir, file), read_hdr),
    n_ref = map_int(ref, "n"),
    n_new = map_int(new, "n"),
    channels_equal = map2_lgl(ref, new, \(r, n) {
      identical(r$channels, n$channels)
    })
  ) |>
  select(file, n_ref, n_new, channels_equal)
write_csv(verify, file.path(out_dir, "verify.csv"))
stopifnot(
  "5 files unmixed" = nrow(verify) == 5,
  "event counts equal as188_joint" = all(verify$n_ref == verify$n_new),
  "channel names equal as188_joint" = all(verify$channels_equal)
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
