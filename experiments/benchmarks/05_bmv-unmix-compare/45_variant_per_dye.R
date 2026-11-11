# Purpose: find which dye's AutoSpectral 1.0.2 spectral variants make BV510 (TNF) collapse toward 0.
#          One arm per non-BV510 dye: the 1.8.8 variants object (as in var102_bv510only of
#          39_variant_count.R) with that one dye's variants replaced by its full 1.0.2 set. Several dyes
#          joined by "+" swap together in one arm (SBUV510 + AF405, user 2026-11-10). Everything else
#          is the 37/39 setup (1.8.8 code, joint, speed slow, 1.8.8 dye spectra, raw 1.8.8 AF). No subsampling,
#          so no seed. Scoring is in 45b_score_variant_per_dye.R.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/45_variant_per_dye.R "<dye>"
#            from the repo root; dye is a row name of the 1.8.8 spectra other than BV510, e.g. "Alexa Fluor 405",
#            or several joined by "+", e.g. "SBUV510+Alexa Fluor 405" (arm var102_only_SBUV510__Alexa_Fluor_405).
#            Output goes to experiments/bmv_test/compare_3D/data/unmix_cells/var102_only_<dye_id>, which must
#            not exist (dye_id = dye with spaces and dashes as "_").
#            On SLURM: experiments/slurm/run_variant_per_dye.sbatch (array over the dye list).
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (run 143803 dye spectra, AF, variants,
#            Reference Group Clean controls, Raw, AutoSpectral_unmixed) [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2 (Spectral_variants.rds) [user zip]
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
PIPELINE <- "joint" # as in 39_variant_count.R (inherited: user zip)
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint" (repro188 = s188_c188_slow); was "fast" to match as188_joint
TARGET_DYE <- "BV510" # the TNF dye; its own variants stay 1.8.8 here (var102_bv510only showed no effect)
VARIANT_SLOTS <- c("variants", "delta.list", "delta.norms") # per-dye lists swapped together, as in 39

# Setup -----------------------------------------------------------------------------------
dyes <- strsplit(commandArgs(trailingOnly = TRUE)[[1]], "+", fixed = TRUE)[[1]]
dye_id <- paste(gsub("[ -]", "_", dyes), collapse = "__")
arm <- paste0("var102_only_", dye_id)
out_dir <- file.path(DATA, "unmix_cells", arm)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)

# Inputs: 1.8.8 dye spectra, AF and variants; 1.0.2 variants ----------------------------------
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
stopifnot(
  "dyes are non-BV510 fluors of the spectra" = all(
    dyes %in% setdiff(rownames(spectra), TARGET_DYE)
  ),
  "variants 1.0.2 fluors equal dye spectra fluors" = identical(
    names(variants102$variants),
    rownames(spectra)
  ),
  "variants 1.8.8 fluors equal dye spectra fluors" = identical(
    names(variants188$variants),
    rownames(spectra)
  )
)
dir.create(out_dir, recursive = TRUE)

# Variants object: 1.8.8 object, the chosen dyes' three per-dye lists from 1.0.2 --------------------
spectra_variants <- modifyList(
  variants188,
  map(set_names(VARIANT_SLOTS), \(slot) {
    out <- variants188[[slot]]
    out[dyes] <- variants102[[slot]][dyes]
    out
  })
)
saveRDS(spectra_variants, file.path(out_dir, "variants_used.rds"))

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste("dyes with 1.0.2 variants:", paste(dyes, collapse = ", ")),
    paste("speed:", SPEED),
    paste("pipeline:", PIPELINE),
    paste("spectra:", spectra_file),
    paste("af:", af188_file),
    paste("variants: 1.8.8", var188_rds, "; dye from 1.0.2", var102_rds),
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
