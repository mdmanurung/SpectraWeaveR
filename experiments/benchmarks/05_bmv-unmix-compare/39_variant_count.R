# Purpose: decide whether the NUMBER or the CONTENT of the AutoSpectral 1.0.2 spectral variants makes BV510
#          (TNF) collapse toward 0. Same setup as 37_swap_af_variants.R (1.8.8 code, joint, speed slow, 1.8.8
#          dye spectra, raw 1.8.8 AF); only the variants object changes. Scoring is in 39b_score_variant_count.R.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/39_variant_count.R <arm>
#            from the repo root. arm = var102_n26, var188_n100, var102_bv510only or var102_except_bv510.
#            Output goes to experiments/bmv_test/compare_3D/data/unmix_cells/<arm>, which must not exist.
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (run 143803 dye spectra, AF, variants,
#            Reference Group Clean controls, Raw, AutoSpectral_unmixed) [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2 (Spectral_variants.rds) [user zip]
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
SEED <- 42 # seed for the subsample and the resample; a choice the user must confirm
N_CAP <- 26 # 1.8.8 variants per dye; the cap for var102_n26 and the threshold for var188_n100
SWAP_DYE <- "BV510" # dye whose variants are swapped in the two bv510 arms
ARMS <- c(
  "var102_n26",
  "var188_n100",
  "var102_bv510only",
  "var102_except_bv510"
)

# Setup -----------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
arm <- cli_args[[1]]
out_dir <- file.path(DATA, "unmix_cells", arm)
stopifnot("arm is one of ARMS" = arm %in% ARMS)
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
variants188 <- readRDS(var188_rds)
variants102 <- readRDS(var102_rds)

# Detector set and fluor set must match; stop, do not coerce
stopifnot(
  "AF 1.8.8 detectors equal dye spectra detectors" = identical(
    colnames(af188),
    colnames(spectra)
  ),
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

# Variants object per arm. Slots the 1.0.2 object lacks (neg.thresholds etc.) stay missing, as in
# swap_var102 of 37 (that arm ran with the unfilled 1.0.2 object); 1.8.8-based arms keep their slots.
# Subset the three per-dye lists together; keep each matrix's attributes (noise.events etc.)
pick_rows <- \(v, idx, new_names = FALSE) {
  v$variants <- imap(v$variants, \(m, dye) {
    out <- m[idx[[dye]], , drop = FALSE]
    keep_attr <- attributes(m)[setdiff(
      names(attributes(m)),
      c("dim", "dimnames")
    )]
    attributes(out) <- c(attributes(out), keep_attr)
    if (new_names) {
      rownames(out) <- paste0(dye, "_", seq_len(nrow(out)))
    }
    out
  })
  v$delta.list <- imap(v$delta.list, \(m, dye) {
    out <- m[idx[[dye]], , drop = FALSE]
    keep_attr <- attributes(m)[setdiff(
      names(attributes(m)),
      c("dim", "dimnames")
    )]
    attributes(out) <- c(attributes(out), keep_attr)
    if (new_names) {
      rownames(out) <- paste0(dye, "_", seq_len(nrow(out)))
    }
    out
  })
  v$delta.norms <- imap(v$delta.norms, \(x, dye) {
    out <- x[idx[[dye]]]
    if (new_names) {
      names(out) <- paste0(dye, "_", seq_along(out))
    }
    out
  })
  v
}
set.seed(SEED)
n102 <- map_int(variants102$variants, nrow)
n188 <- map_int(variants188$variants, nrow)
spectra_variants <- switch(
  arm,
  # 1.0.2 variants, at most N_CAP per dye (random rows, original order)
  var102_n26 = pick_rows(
    variants102,
    map(n102, \(n) if (n > N_CAP) sort(sample.int(n, N_CAP)) else seq_len(n))
  ),
  # 1.8.8 variants resampled with replacement up to the 1.0.2 count where that exceeds N_CAP
  var188_n100 = pick_rows(
    variants188,
    map2(n188, n102, \(a, b) {
      if (b > N_CAP) sample.int(a, b, replace = TRUE) else seq_len(a)
    }),
    new_names = TRUE
  ),
  # 1.8.8 object, SWAP_DYE from 1.0.2
  var102_bv510only = modifyList(
    variants188,
    map(
      c(
        variants = "variants",
        delta.list = "delta.list",
        delta.norms = "delta.norms"
      ),
      \(slot) {
        out <- variants188[[slot]]
        out[[SWAP_DYE]] <- variants102[[slot]][[SWAP_DYE]]
        out
      }
    )
  ),
  # 1.0.2 object, SWAP_DYE from 1.8.8
  var102_except_bv510 = modifyList(
    variants102,
    map(
      c(
        variants = "variants",
        delta.list = "delta.list",
        delta.norms = "delta.norms"
      ),
      \(slot) {
        out <- variants102[[slot]]
        out[[SWAP_DYE]] <- variants188[[slot]][[SWAP_DYE]]
        out
      }
    )
  )
)
af_spectra <- af188
saveRDS(spectra_variants, file.path(out_dir, "variants_used.rds"))

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste("speed:", SPEED),
    paste("pipeline:", PIPELINE),
    paste("spectra:", spectra_file),
    paste("af:", af188_file),
    paste("variants:", var102_rds, "and", var188_rds, "combined per arm"),
    paste("seed:", SEED, " n_cap:", N_CAP, " swap_dye:", SWAP_DYE),
    paste(
      "variants per dye:",
      paste(map_int(spectra_variants$variants, nrow), collapse = ",")
    ),
    paste("controls:", control_dir),
    "weighted: default (as 06)",
    "csv read: raw as.matrix(read.csv(row.names = 1, check.names = FALSE))",
    "extra slots: 1.0.2-based arms lack neg.thresholds etc. (as swap_var102 of 37)"
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
