# Purpose: test 4 (README decision table, user 2026-11-10): is the widened raw-file bead gate the source of
#          SBUV510's BV510 part? Arm sbuv510_beadgate = the zip 1.8.8 setup with only SBUV510's spectrum row
#          and its per-dye variant lists (variants, delta.list, delta.norms, as swapped in 45) taken from
#          chain5_bead (49_ablation_chain.R: bead controls re-gated from the raw files, widened rectangle).
#          Thresholds and the other variant slots stay 1.8.8, as in the per-dye arms. Everything else is the
#          45 setup (1.8.8 code, joint, speed slow, 1.8.8 AF). No subsampling, so no seed. Score with 51.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/78_sbuv510_beadgate_arm.R  (repo root)
#            On SLURM: experiments/slurm/run_script.sbatch with this script.
#            Output goes to experiments/bmv_test/compare_3D/data/unmix_cells/sbuv510_beadgate (must not exist).
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (run 143803 dye spectra, AF, variants,
#            Reference Group Clean controls, Raw, AutoSpectral_unmixed) [user zip]
#          experiments/bmv_test/compare_3D/data/unmix_cells/chain5_bead/objects/spectra.rds,
#            spectral_variants.rds [49_ablation_chain.R]
# Outputs: <out_dir>/*.fcs               unmixed files
#          <out_dir>/spectra_used.csv    dye spectra passed to unmix.folder()
#          <out_dir>/variants_used.rds   the variants object passed to unmix.folder()
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
SOURCE_ARM <- "chain5_bead" # first chain step with the widened raw-file bead gate
VARIANT_SLOTS <- c("variants", "delta.list", "delta.norms") # per-dye lists swapped together, as in 39

# Setup -----------------------------------------------------------------------------------
dye <- DYE
arm <- "sbuv510_beadgate"
out_dir <- file.path(DATA, "unmix_cells", arm)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)

# Inputs: 1.8.8 dye spectra, AF and variants; chain5_bead spectra and variants ----------------------
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
src_dir <- file.path(DATA, "unmix_cells", SOURCE_ARM, "objects")
# raw csv read: read.spectra() rescales rows (user, 2026-11-10)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra <- read_raw(spectra_file)
af188 <- read_raw(af188_file)
variants188 <- readRDS(var188_rds)
spectra_src <- readRDS(file.path(src_dir, "spectra.rds"))
variants_src <- readRDS(file.path(src_dir, "spectral_variants.rds"))
stopifnot(
  "source spectra fluors and detectors equal 1.8.8" = identical(
    dimnames(spectra_src),
    dimnames(spectra)
  ),
  "source variants fluors equal dye spectra fluors" = identical(
    names(variants_src$variants),
    rownames(spectra)
  ),
  "variants 1.8.8 fluors equal dye spectra fluors" = identical(
    names(variants188$variants),
    rownames(spectra)
  )
)
dir.create(out_dir, recursive = TRUE)

# Spectra and variants: 1.8.8, with SBUV510's row and per-dye lists from chain5_bead ----------------
spectra_variants <- modifyList(
  variants188,
  map(set_names(VARIANT_SLOTS), \(slot) {
    out <- variants188[[slot]]
    out[[dye]] <- variants_src[[slot]][[dye]]
    out
  })
)
spectra[dye, ] <- spectra_src[dye, ]
write.csv(spectra, file.path(out_dir, "spectra_used.csv"))
saveRDS(spectra_variants, file.path(out_dir, "variants_used.rds"))

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste("dye with chain5_bead spectrum row and variants:", dye),
    paste("speed:", SPEED),
    paste("pipeline:", PIPELINE),
    paste("spectra:", spectra_file, "; SBUV510 row from", src_dir),
    paste("af:", af188_file),
    paste("variants: 1.8.8", var188_rds, "; SBUV510 lists from", src_dir),
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
