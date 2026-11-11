# Purpose: test 8 (README decision table, user 2026-11-10): split SBUV510's TNF effect into its reference
#          spectrum and its variants. Arm sbuv510_refonly = the zip 1.8.8 setup with only SBUV510's spectrum
#          row replaced by the zip 1.0.2 'Cells Matching Negative' row; all variants stay 1.8.8. The
#          variants-only arm is var102_only_SBUV510 (45_variant_per_dye.R), reused, not rerun. Everything
#          else is the 45 setup (1.8.8 code, joint, speed slow, 1.8.8 AF). No subsampling, so no seed.
#          Score with 51.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/84_sbuv510_decompose.R  (repo root)
#            On SLURM: experiments/slurm/run_script.sbatch with this script.
#            Output goes to experiments/bmv_test/compare_3D/data/unmix_cells/sbuv510_refonly (must not exist).
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (run 143803 dye spectra, AF, variants,
#            Reference Group Clean controls, Raw, AutoSpectral_unmixed) [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/table_spectra/
#            Cells Matching Negative_autospectral_spectra.csv  1.0.2 dye spectra [user zip]
# Outputs: <out_dir>/*.fcs               unmixed files
#          <out_dir>/spectra_used.csv    dye spectra passed to unmix.folder()
#          <out_dir>/variants_used.rds   the variants object passed to unmix.folder() (1.8.8, unchanged)
#          <out_dir>/args.txt            what was run
#          <out_dir>/verify.csv          event count and channel names vs as188_joint, per file
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
PIPELINE <- "joint" # as in 39_variant_count.R (inherited: user zip)
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint" (repro188 = s188_c188_slow); was "fast" to match as188_joint
DYE <- "SBUV510" # the CD3 dye; second-largest TNF effect in the per-dye test (0.19)

# Setup -----------------------------------------------------------------------------------
dye <- DYE
arm <- "sbuv510_refonly"
out_dir <- file.path(DATA, "unmix_cells", arm)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)

# Inputs: 1.8.8 dye spectra, AF and variants; 1.0.2 dye spectra ------------------------------------
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
var188_rds <- file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)
spectra102_file <- file.path(
  DATA,
  "autoSpectral_v1.0.2",
  "Reference Group",
  "table_spectra",
  "Cells Matching Negative_autospectral_spectra.csv"
)
# raw csv read: read.spectra() rescales rows (user, 2026-11-10)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra <- read_raw(spectra_file)
af188 <- read_raw(af188_file)
variants188 <- readRDS(var188_rds)
spectra102 <- read_raw(spectra102_file)
stopifnot(
  "1.0.2 spectra hold the dye and the 1.8.8 detectors" = dye %in%
    rownames(spectra102) &&
    all(colnames(spectra) %in% colnames(spectra102)),
  "variants 1.8.8 fluors equal dye spectra fluors" = identical(
    names(variants188$variants),
    rownames(spectra)
  )
)
dir.create(out_dir, recursive = TRUE)

# Spectra and variants: 1.8.8, with SBUV510's spectrum row from 1.0.2 (selected by detector name) --------
spectra_variants <- variants188
spectra[dye, ] <- spectra102[dye, colnames(spectra)]
write.csv(spectra, file.path(out_dir, "spectra_used.csv"))
saveRDS(spectra_variants, file.path(out_dir, "variants_used.rds"))

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste("dye with the 1.0.2 spectrum row:", dye),
    paste("speed:", SPEED),
    paste("pipeline:", PIPELINE),
    paste("spectra:", spectra_file, "; SBUV510 row from", spectra102_file),
    paste("af:", af188_file),
    paste("variants: 1.8.8", var188_rds, "(unchanged)"),
    paste(
      "variants per dye:",
      paste(map_int(spectra_variants$variants, nrow), collapse = ",")
    ),
    paste("controls:", control_dir),
    "weighted: default (as 06)",
    "csv read: raw as.matrix(read.csv(row.names = 1, check.names = FALSE))",
    "seed: none (no subsampling)"
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
  af.spectra = af188,
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
