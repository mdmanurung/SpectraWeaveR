# Purpose: causal tests of why Alexa Fluor 405's AutoSpectral 1.0.2 variants squeeze BV510 (TNF). Each arm is
#          the 1.8.8 variants object (as in 45_variant_per_dye.R) with one change to the Alexa Fluor 405 entry:
#            af405_102_<subset>  a subset of the 65 1.0.2 AF405 variants (rows from 54: median splits on the
#                                delta's similarity to AF or to BV510, or 26 random rows)
#            af405_102_thr       the 1.8.8 variants, with only the AF405 positivity threshold set to 1.0.2's
#            af405_188_noattr    the 1.8.8 variants, with the attributes on the AF405 delta matrix removed
#                                (the swapped-in 1.0.2 delta has none; tests whether their loss matters)
#          Everything else as 45: 1.8.8 code, joint pipeline, 1.8.8 dye spectra and raw AF, speed "slow".
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/55_af405_subset_unmix.R <arm>
#            from the repo root, after 54. Output goes to
#            experiments/bmv_test/compare_3D/data/unmix_cells/<arm>, which must not exist.
#            On SLURM: experiments/slurm/run_af405_arm.sbatch (array over the arm list).
# Inputs:  experiments/runs/2026-11-10_bmv-af405-mechanism/variant_subsets.rds [54_af405_variant_geometry.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 and autoSpectral_v1.0.2 [user zip]
# Outputs: <out_dir>/*.fcs               unmixed files
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
SUBSETS_RDS <- "experiments/runs/2026-11-10_bmv-af405-mechanism/variant_subsets.rds"
DYE <- "Alexa Fluor 405"
PIPELINE <- "joint" # as 45 (inherited: user zip)
SPEED <- "slow" # every AutoSpectral unmix call (user, 2026-11-10); no effect under joint
VARIANT_SLOTS <- c("variants", "delta.list", "delta.norms") # per-dye lists swapped together, as in 45

# Setup -----------------------------------------------------------------------------------
arm <- commandArgs(trailingOnly = TRUE)[[1]]
out_dir <- file.path(DATA, "unmix_cells", arm)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)
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
# raw csv read: read.spectra() rescales rows (user, 2026-11-10)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra <- read_raw(spectra_file)
af188 <- read_raw(af188_file)
variants188 <- readRDS(var188_rds)
variants102 <- readRDS(var102_rds)
subsets <- readRDS(SUBSETS_RDS)
stopifnot(
  "arm is a defined AF405 arm" = arm %in%
    c(names(subsets), "af405_102_thr", "af405_188_noattr")
)

# Variants object for this arm -------------------------------------------------------------
spectra_variants <- variants188
change <- if (arm %in% names(subsets)) {
  # subset rows of the 1.0.2 AF405 entry, the same rows in all three per-dye slots
  rows <- subsets[[arm]]
  for (slot in VARIANT_SLOTS) {
    x <- variants102[[slot]][[DYE]]
    spectra_variants[[slot]][[DYE]] <- if (is.matrix(x)) {
      x[rows, , drop = FALSE]
    } else {
      x[rows]
    }
  }
  paste0(
    "1.0.2 AF405 variant rows ",
    paste(rows, collapse = ","),
    " (",
    length(rows),
    " of 65)"
  )
} else if (arm == "af405_102_thr") {
  spectra_variants$thresholds[[DYE]] <- variants102$thresholds[[DYE]]
  sprintf(
    "AF405 threshold %.1f (1.0.2) instead of %.1f (1.8.8)",
    variants102$thresholds[[DYE]],
    variants188$thresholds[[DYE]]
  )
} else {
  d <- spectra_variants$delta.list[[DYE]]
  spectra_variants$delta.list[[DYE]] <- matrix(
    d,
    nrow(d),
    ncol(d),
    dimnames = dimnames(d)
  )
  paste(
    "AF405 delta.list attributes removed:",
    paste(setdiff(names(attributes(d)), c("dim", "dimnames")), collapse = ",")
  )
}
dir.create(out_dir, recursive = TRUE)
saveRDS(spectra_variants, file.path(out_dir, "variants_used.rds"))

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste("change:", change),
    paste("speed:", SPEED),
    paste("pipeline:", PIPELINE),
    paste("spectra:", spectra_file),
    paste("af:", af188_file),
    paste("variants: 1.8.8", var188_rds, "; AF405 change from", var102_rds),
    paste(
      "variants per dye:",
      paste(map_int(spectra_variants$variants, nrow), collapse = ",")
    ),
    paste("controls:", control_dir),
    "csv read: raw as.matrix(read.csv(row.names = 1, check.names = FALSE))",
    "seed: random-26 rows drawn in 54 with set.seed(42); no other randomness here"
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
