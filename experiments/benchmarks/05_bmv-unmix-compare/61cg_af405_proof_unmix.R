# Copy of 61_af405_proof_unmix.R for the common-gates spike (plans/common_gates_plan.md G5): the unmix is unchanged;
#          verify.csv covers the files of the given raw_dir (one donor in a smoke run), not the five of as188_joint.
# Purpose: unmix one proof-chain arm (user, 2026-11-10): the 1.8.8 zip setup (1.8.8 code, joint,
#          speed "slow" (user rule), run-143803 dye spectra, raw-csv 1.8.8 AF, 1.8.8 thresholds) with a
#          given variants object, on a given folder of raw FCS files (the zip Raw folder, or spiked copies).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/61_af405_proof_unmix.R \
#            <arm> <variants_rds|188> [raw_dir]
#          from the repo root; on SLURM via experiments/slurm/run_spike_cg.sbatch. "188" = the zip 1.8.8
#          variants object unchanged. raw_dir defaults to the zip Raw folder. Output goes to
#          experiments/bmv_test/compare_3D/data/unmix_cells/af405proof_<arm>, which must not exist.
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (spectra, AF, variants, controls, Raw)
#          <variants_rds> [60_af405_proof_objects.R]
# Outputs: <out_dir>/*.fcs   unmixed files
#          <out_dir>/args.txt, verify.csv (event count and channels vs as188_joint), sessionInfo.txt

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Setup -----------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
cli_args <- commandArgs(trailingOnly = TRUE)
arm <- cli_args[[1]]
variants_arg <- cli_args[[2]]
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
raw_dir <- if (length(cli_args) >= 3) {
  normalizePath(cli_args[[3]])
} else {
  file.path(d188, "Raw")
}
out_dir <- file.path(DATA, "unmix_cells", paste0("af405proof_", arm))
stopifnot("out dir must not exist" = !dir.exists(out_dir))
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)
dir.create(out_dir, recursive = TRUE)

# Inputs: the zip 1.8.8 objects; variants from the arm ---------------------------------------
control_dir <- file.path(d188, "Reference Group Clean")
spectra_file <- file.path(
  d188,
  "table_spectra",
  "Automated_autospectral_spectra_20261105_143803.csv"
)
af_file <- file.path(
  d188,
  "table_spectra",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
)
var188_rds <- file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
} # no row rescale
spectra <- read_raw(spectra_file)
af <- read_raw(af_file)
variants_file <- if (variants_arg == "188") {
  var188_rds
} else {
  normalizePath(variants_arg)
}
spectra_variants <- readRDS(variants_file)

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste("variants:", variants_file),
    paste("raw dir:", raw_dir),
    "pipeline: joint",
    "speed: slow (user rule, 2026-11-10; no effect under joint)",
    paste("spectra:", spectra_file),
    paste("af:", af_file, "(raw csv, no row rescale)"),
    paste("controls:", control_dir),
    "thresholds: 1.8.8 (from the variants object)"
  ),
  file.path(out_dir, "args.txt")
)
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(
  control_dir,
  file.path(control_dir, "fcs_control_file.csv"),
  asp
)
unmix.folder(
  fcs.dir = raw_dir,
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af,
  spectra.variants = spectra_variants,
  pipeline = "joint",
  speed = "slow",
  output.dir = out_dir,
  parallel = TRUE,
  verbose = TRUE
)

# Verify: same files, event counts and channel names as as188_joint -----------------------------
ref_dir <- file.path(d188, "AutoSpectral_unmixed")
read_hdr <- \(path) {
  ff <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  list(n = nrow(ff), channels = colnames(ff))
}
verify <- tibble(
  file = sub(
    "\\.fcs$",
    " AutoSpectral.fcs",
    list.files(raw_dir, pattern = "\\.fcs$")
  )
) |>
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
  "every raw file unmixed" = nrow(verify) ==
    length(list.files(raw_dir, pattern = "\\.fcs$")),
  "event counts equal as188_joint" = all(verify$n_ref == verify$n_new),
  "channel names equal as188_joint" = all(verify$channels_equal)
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
