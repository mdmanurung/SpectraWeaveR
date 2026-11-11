# Purpose: unmix the 5 BMV files with the AutoSpectral 1.8.8 noise model and GLS (estimate.noise.model()
#          then unmix.gls()), the Honeychrome-style noise weighting that unmix.folder() never calls.
#          Same 1.8.8 inputs as as188_joint: run 143803 dye spectra, raw-csv AF, Spectral_variants.rds,
#          Reference Group Clean controls. Per-cell AF comes from assign.af.joint.cov(); per-cell spectral
#          variants are committed inside unmix.gls() (spectra.variants), as in the joint pipeline.
# Run:     from the repo root, in this order:
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/42_unmix_gls_noise.R noise
#              fits the noise model once; both arms load it
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/42_unmix_gls_noise.R <arm> <mode>
#              arm = gls188_noise (include.spillover FALSE) or gls188_noise_spill (TRUE)
#              mode = smoke: donor A3 only, written to <RUN_DIR>/smoke_<arm>/, then verified
#                     A3, A6, B11, B2 or D1: that donor into the arm folder (one SLURM job each;
#                       unmix.gls() takes about 2 core-minutes per 1000 events)
#                     verify: after all 5 donors, checks the arm folder against as188_joint
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (Raw, Reference Group Clean,
#            table_spectra, figure_spectral_variants, AutoSpectral_unmixed) [user zip]
# Outputs: <RUN_DIR>/noise_model/noise_model.rds        fitted noise model (noise step)
#          <RUN_DIR>/noise_model/detectors.csv          sigma0^2, kappa and flags per detector
#          <RUN_DIR>/noise_model/calibration_sources.csv held-out calibration, one row per source
#          <RUN_DIR>/noise_model/summary.txt            in-band shares, clamped flag, notes
#          <out_dir>/*.fcs           unmixed files, same channels as as188_joint
#          <out_dir>/args.txt        what was run
#          <out_dir>/noise_model.rds copy of the noise model used
#          <out_dir>/verify.csv      event count and channel names vs as188_joint, per file
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
# All noise-model and GLS settings are the AutoSpectral 1.8.8 package defaults (user, 2026-11-10).
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
RUN_DIR <- "experiments/runs/2026-11-10_bmv-gls-noise"
ARMS <- c(gls188_noise = FALSE, gls188_noise_spill = TRUE) # include.spillover per arm (user, 2026-11-10)
SMOKE_DONOR <- "A3" # whole donor for the smoke test (all events)
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "6")) # parallel chunks; unmix.gls() is per cell, so chunking does not change results
DONORS <- c("A3", "A6", "B11", "B2", "D1")
# estimate.noise.model() (user, 2026-11-10: package defaults)
NOISE_SELECTION <- "auto" # scatter core and FSC-A/FSC-H singlets of every control file
NOISE_UNSTAINED_TRAIN_FRAC <- 0.5 # training share of the unstained cells
NOISE_CONTROL_TRAIN_FRAC <- 0.7 # training share of each control; the rest checks the fit
NOISE_N_AF_ROWS <- 4L # AF design rows of a cell negative
NOISE_FIT_QUADRATIC <- FALSE # no per-detector c * mu^2 term
NOISE_FLOOR_FRAC <- 0.01 # lower clamp on sigma0^2, as a fraction of its median
NOISE_BEAD_FACTOR <- 10 # minimum positive / internal bead negative ratio
NOISE_CLAMP_MODE <- "auto" # detect data clamped at zero; then clamped.weights are recommended
NOISE_SEED <- NULL # NULL = asp$bird.seed, the package default (user, 2026-11-10)
# build.variant.basis() (package defaults; built with the spectra, as unmix.gls() docs advise)
BASIS_RANK <- 6L # maximum components per fluorophore
BASIS_VAR_EXPLAINED <- 0.99 # components kept until this share of delta variance
# unmix.gls() (user, 2026-11-10: package defaults; include.spillover from ARMS)
GLS_SPREAD_KAPPA <- 2 # spread SDs added to the positivity threshold for variant commitment
GLS_N_SELECT_SWEEPS <- 2L # variant selection sweeps per cell
GLS_N_ITER <- 2L # GLS refinement iterations
GLS_GAIN_CV <- 0 # multiplicative gain term off
GLS_ACTIVE_THRESHOLD <- 1e-3 # fluorophores below this share of a cell's largest abundance leave Sigma

# Setup -----------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
step <- cli_args[[1]]
mode <- if (length(cli_args) > 1) cli_args[[2]] else ""
smoke <- mode == "smoke"
stopifnot(
  "step is noise or an arm" = step %in% c("noise", names(ARMS)),
  "arm mode is smoke, a donor or verify" = step == "noise" ||
    mode %in% c("smoke", DONORS, "verify"),
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
control_dir <- file.path(d188, "Reference Group Clean")
control_file <- file.path(control_dir, "fcs_control_file.csv")
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
variants_rds <- file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)
noise_dir <- file.path(RUN_DIR, "noise_model")

# Raw csv read: read.spectra() rescales rows and caused an offset (user, 2026-11-10)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra <- read_raw(spectra_file)
af_spectra <- read_raw(af_file)
variants <- readRDS(variants_rds)
stopifnot(
  "AF detectors equal dye spectra detectors" = identical(
    colnames(af_spectra),
    colnames(spectra)
  )
)
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)

# Noise step: fit the noise model once, shared by both arms ------------------------------------
if (step == "noise") {
  stopifnot("noise model dir must not exist" = !dir.exists(noise_dir))
  dir.create(noise_dir, recursive = TRUE)
  t0 <- Sys.time()
  noise_model <- estimate.noise.model(
    control.dir = control_dir,
    control.def.file = control_file,
    asp = asp,
    spectra = spectra,
    spectral.variants = variants,
    output.dir = noise_dir,
    save = FALSE,
    selection = NOISE_SELECTION,
    unstained.train.frac = NOISE_UNSTAINED_TRAIN_FRAC,
    control.train.frac = NOISE_CONTROL_TRAIN_FRAC,
    n.af.rows = NOISE_N_AF_ROWS,
    fit.quadratic = NOISE_FIT_QUADRATIC,
    floor.frac = NOISE_FLOOR_FRAC,
    bead.factor = NOISE_BEAD_FACTOR,
    clamp.mode = NOISE_CLAMP_MODE,
    seed = if (is.null(NOISE_SEED)) asp$bird.seed else NOISE_SEED
  )
  saveRDS(noise_model, file.path(noise_dir, "noise_model.rds"))
  tibble(
    detector = names(noise_model$read.var),
    sigma0_sq = unname(noise_model$read.var),
    kappa = unname(noise_model$counts.per.unit[names(noise_model$read.var)]),
    flags = unname(noise_model$flags[names(noise_model$read.var)])
  ) |>
    write_csv(file.path(noise_dir, "detectors.csv"))
  calibration_sources <- noise_model$calibration |>
    keep(is.data.frame) |>
    pluck(length(keep(noise_model$calibration, is.data.frame)))
  write_csv(
    as_tibble(calibration_sources),
    file.path(noise_dir, "calibration_sources.csv")
  )
  writeLines(
    c(
      paste(
        "fit minutes:",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
      ),
      paste("seed:", if (is.null(NOISE_SEED)) asp$bird.seed else NOISE_SEED),
      paste(
        "in.band:",
        paste(format(noise_model$calibration$in.band), collapse = " ")
      ),
      paste(
        "in.band.controls:",
        paste(format(noise_model$calibration$in.band.controls), collapse = " ")
      ),
      paste("clamped:", paste(format(noise_model$clamped), collapse = " ")),
      paste("flag.text:", paste(noise_model$flag.text, collapse = " | ")),
      paste("notes:", paste(noise_model$notes, collapse = " | "))
    ),
    file.path(noise_dir, "summary.txt")
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(noise_dir, "sessionInfo.txt")
  )
  quit(save = "no")
}

# Arm step: per-cell AF, then GLS with the shared noise model --------------------------------
arm <- step
out_dir <- if (smoke) {
  file.path(RUN_DIR, paste0("smoke_", arm))
} else {
  file.path(DATA, "unmix_cells", arm)
}
# smoke: a new folder; donor jobs share the arm folder but never overwrite a file (checked per file)
if (smoke) {
  stopifnot("smoke dir must not exist" = !dir.exists(out_dir))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
noise_model <- readRDS(file.path(noise_dir, "noise_model.rds"))
file.copy(file.path(noise_dir, "noise_model.rds"), out_dir, overwrite = FALSE)
variant_basis <- build.variant.basis(
  variants,
  spectra = spectra,
  rank = BASIS_RANK,
  var.explained = BASIS_VAR_EXPLAINED
)
writeLines(
  c(
    paste("arm:", arm, if (smoke) "(smoke, donor A3)" else ""),
    "method: estimate.noise.model() + assign.af.joint.cov() + unmix.gls()",
    paste(
      "include.spillover:",
      ARMS[[arm]],
      " spillover.kappa: default (kappa.pooled)"
    ),
    paste("spectra:", spectra_file),
    paste("af:", af_file, "(raw csv read)"),
    paste("variants:", variants_rds),
    paste("controls:", control_dir),
    paste("noise model:", file.path(noise_dir, "noise_model.rds")),
    paste(
      "gls: spread.kappa",
      GLS_SPREAD_KAPPA,
      "n.select.sweeps",
      GLS_N_SELECT_SWEEPS,
      "n.iter",
      GLS_N_ITER,
      "gain.cv",
      GLS_GAIN_CV,
      "active.threshold",
      GLS_ACTIVE_THRESHOLD
    ),
    paste(
      "variant basis: rank",
      BASIS_RANK,
      "var.explained",
      BASIS_VAR_EXPLAINED
    ),
    "af basis: none (raw-csv AF has no af.model)"
  ),
  file.path(out_dir, "args.txt")
)

raw_files <- list.files(
  file.path(d188, "Raw"),
  pattern = "\\.fcs$",
  full.names = TRUE
)
if (mode != "verify") {
  donor <- if (smoke) SMOKE_DONOR else mode
  raw_files <- raw_files[startsWith(basename(raw_files), paste0(donor, " "))]
}
ref_dir <- file.path(d188, "AutoSpectral_unmixed")
for (raw_file in if (mode == "verify") character() else raw_files) {
  t0 <- Sys.time()
  ff <- read.FCS(raw_file, transformation = FALSE, truncate_max_range = FALSE)
  raw_data <- exprs(ff)[, colnames(spectra)]
  # per-cell AF spectrum, then GLS on chunks of cells (per-cell fits, so chunking is exact)
  af_index <- assign.af.joint.cov(raw_data, spectra, af_spectra)
  chunks <- split(
    seq_len(nrow(raw_data)),
    cut(seq_len(nrow(raw_data)), N_CORES, labels = FALSE)
  )
  unmixed <- parallel::mclapply(
    chunks,
    \(rows) {
      unmix.gls(
        raw.data = raw_data[rows, , drop = FALSE],
        spectra = spectra,
        noise.model = noise_model,
        variant.basis = variant_basis,
        spectra.variants = variants,
        spread.kappa = GLS_SPREAD_KAPPA,
        n.select.sweeps = GLS_N_SELECT_SWEEPS,
        af.spectra = af_spectra,
        af.index = af_index[rows],
        n.iter = GLS_N_ITER,
        include.spillover = ARMS[[arm]],
        gain.cv = GLS_GAIN_CV,
        active.threshold = GLS_ACTIVE_THRESHOLD,
        verbose = FALSE
      )
    },
    mc.cores = N_CORES
  )
  failed <- map_lgl(unmixed, \(x) inherits(x, "try-error") || !is.matrix(x))
  if (any(failed)) {
    stop(sprintf(
      "unmix.gls failed in %d chunk(s): %s",
      sum(failed),
      unmixed[failed][[1]]
    ))
  }
  unmixed <- do.call(rbind, unmixed)

  # Write with the as188_joint channel order: scatter and Time from the raw file, dyes, AF, AF Index
  ref_file <- file.path(
    ref_dir,
    sub("\\.fcs$", " AutoSpectral.fcs", basename(raw_file))
  )
  ref_channels <- colnames(read.FCS(
    ref_file,
    which.lines = 1:2,
    truncate_max_range = FALSE
  ))
  af_col <- setdiff(colnames(unmixed), rownames(spectra))
  out <- cbind(
    exprs(ff)[, setdiff(colnames(exprs(ff)), colnames(spectra))],
    unmixed[, rownames(spectra), drop = FALSE],
    `AF-A` = if (length(af_col) == 1) unmixed[, af_col] else NA_real_,
    `AF Index` = af_index
  )
  colnames(out)[colnames(out) %in% rownames(spectra)] <- paste0(
    colnames(out)[colnames(out) %in% rownames(spectra)],
    "-A"
  )
  stopifnot(
    "output channels equal as188_joint" = setequal(colnames(out), ref_channels)
  )
  out_file <- file.path(out_dir, basename(ref_file))
  stopifnot("output file must not exist" = !file.exists(out_file))
  write.FCS(
    flowFrame(out[, ref_channels]),
    out_file
  )
  message(
    basename(raw_file),
    ": ",
    nrow(out),
    " events, AF column '",
    paste(af_col, collapse = ","),
    "', minutes ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  )
}

# a donor job stops here; the verify step checks the arm once all 5 donors are written
if (mode %in% DONORS) {
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, paste0("sessionInfo_", mode, ".txt"))
  )
  quit(save = "no")
}

# Verify: same files, event counts and channel names as as188_joint -----------------------------
read_hdr <- \(path) {
  ff <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  list(n = nrow(ff), channels = colnames(ff))
}
verify <- tibble(file = list.files(out_dir, pattern = "\\.fcs$")) |>
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
  "all files unmixed" = nrow(verify) == length(raw_files),
  "event counts equal as188_joint" = all(verify$n_ref == verify$n_new),
  "channel names equal as188_joint" = all(verify$channels_equal)
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
