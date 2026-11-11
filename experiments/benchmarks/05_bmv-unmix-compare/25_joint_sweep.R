# Purpose: one-factor sweep around the AutoSpectral 1.8.8 joint pipeline (as188_joint): change exactly one
#          unmix argument per arm, run on the 5 raw BMV files. Arm repro188 changes nothing (baseline check).
# Run:     pixi run Rscript 25_joint_sweep.R <arm>   (from the repo root; run each arm as its own process)
#            arm = repro188, jpr_false, coll03, coll07, alpha025 or alpha075
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (Reference Group Clean, run-143803 spectra, its AF
#          spectra, Spectral_variants.rds, Raw) [user zip]; reference files in unmix_cells/repro_s188_c188_fast (as188_joint)
# Outputs: unmix_cells/<arm>/ with 5 FCS, args.txt, verify.csv, sessionInfo.txt (repro188 also repro_check.csv)

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
PIPELINE <- "joint" # as in as188_joint (inherited: unmix_cells/repro_s188_c188_fast/args.txt)
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint" (repro188 = s188_c188_slow); was "fast" to match as188_joint
N_PASSES <- 1L # as188_joint: 1 (1.8.8 default); fixed in every arm
THREADS <- 1L # 6 arms share 8 cores, so 1 thread each; this does not change per-cell results (repro188 checks it)
REF_DIR <- file.path(DATA, "unmix_cells", "repro_s188_c188_fast") # as188_joint files
# one changed argument per arm; the 1.8.8 defaults are joint.pair.resolution TRUE, collinear.threshold 0.5, alpha 0.5
ARMS <- list(
  repro188 = list(),
  jpr_false = list(joint.pair.resolution = FALSE),
  coll03 = list(collinear.threshold = 0.3),
  coll07 = list(collinear.threshold = 0.7),
  alpha025 = list(alpha = 0.25),
  alpha075 = list(alpha = 0.75)
)

# Setup -----------------------------------------------------------------------------------
arm <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("arm is known" = arm %in% names(ARMS))
out_dir <- file.path(DATA, "unmix_cells", arm)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
suppressPackageStartupMessages(library(AutoSpectral))
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)
defaults <- formals(unmix.fcs)[c(
  "joint.pair.resolution",
  "collinear.threshold",
  "alpha",
  "n.passes"
)]
stopifnot(
  "defaults as expected" = identical(defaults$joint.pair.resolution, TRUE),
  identical(defaults$collinear.threshold, 0.5),
  identical(defaults$alpha, 0.5)
)

# Inputs ----------------------------------------------------------------------------------
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
control_dir <- file.path(d188, "Reference Group Clean")
spectra_dir <- file.path(d188, "table_spectra")
spectra_file <- "Automated_autospectral_spectra_20261105_143803.csv"
af_file <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
variants_rds <- file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(
  control_dir,
  file.path(control_dir, "fcs_control_file.csv"),
  asp
)
spectra <- read.spectra(spectra_file, spectra.dir = spectra_dir)
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in 06_unmix_c188.R
af.spectra <- read.spectra(af_file, spectra.dir = spectra_dir)
spectra.variants <- readRDS(variants_rds)

# Unmix -----------------------------------------------------------------------------------
changed <- ARMS[[arm]]
writeLines(
  c(
    paste("arm:", arm),
    paste(
      "changed:",
      if (length(changed)) {
        paste(names(changed), unlist(changed), sep = " = ", collapse = "; ")
      } else {
        "none"
      }
    ),
    paste(
      "defaults 1.8.8: joint.pair.resolution",
      defaults$joint.pair.resolution,
      "collinear.threshold",
      defaults$collinear.threshold,
      "alpha",
      defaults$alpha,
      "n.passes",
      defaults$n.passes
    ),
    paste("n.passes:", N_PASSES),
    paste("speed:", SPEED),
    paste("pipeline:", PIPELINE),
    paste("threads:", THREADS),
    paste("spectra:", file.path(spectra_dir, spectra_file)),
    paste("af:", file.path(spectra_dir, af_file)),
    paste("variants:", variants_rds)
  ),
  file.path(out_dir, "args.txt")
)
t0 <- Sys.time()
do.call(
  unmix.folder,
  c(
    list(
      fcs.dir = file.path(d188, "Raw"),
      spectra = spectra,
      asp = asp,
      flow.control = flow.control,
      method = "AutoSpectral",
      af.spectra = af.spectra,
      spectra.variants = spectra.variants,
      pipeline = PIPELINE,
      speed = SPEED,
      n.passes = N_PASSES,
      output.dir = out_dir,
      parallel = TRUE,
      threads = THREADS,
      verbose = TRUE
    ),
    changed
  )
)
runtime <- difftime(Sys.time(), t0, units = "mins")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# Verify against as188_joint files --------------------------------------------------------
read_ff <- function(f) {
  flowCore::read.FCS(f, truncate_max_range = FALSE, emptyValue = FALSE)
}
ref_files <- sort(list.files(REF_DIR, pattern = "\\.fcs$"))
new_files <- sort(list.files(out_dir, pattern = "\\.fcs$"))
stopifnot(
  "same 5 file names" = identical(ref_files, new_files),
  length(ref_files) == 5
)
rows <- list()
chan <- list()
for (f in ref_files) {
  r <- read_ff(file.path(REF_DIR, f))
  n <- read_ff(file.path(out_dir, f))
  same_chan <- identical(colnames(r), colnames(n))
  rows[[f]] <- data.frame(
    arm = arm,
    file = f,
    n_events_ref = nrow(r),
    n_events_arm = nrow(n),
    channels_equal = same_chan,
    ok = same_chan && nrow(r) == nrow(n)
  )
  if (same_chan && nrow(r) == nrow(n)) {
    d <- abs(flowCore::exprs(r) - flowCore::exprs(n))
    chan[[f]] <- data.frame(
      file = f,
      channel = colnames(d),
      max_abs_diff = apply(d, 2, max)
    )
  }
}
verify <- do.call(rbind, rows)
verify$runtime_min <- round(as.numeric(runtime), 1)
write.csv(verify, file.path(out_dir, "verify.csv"), row.names = FALSE)
chan <- do.call(rbind, chan)
write.csv(
  chan,
  file.path(out_dir, "channel_diff_vs_ref.csv"),
  row.names = FALSE
)
if (arm == "repro188") {
  # baseline check: max abs difference per channel, over all 5 files
  per_channel <- aggregate(max_abs_diff ~ channel, chan, max)
  write.csv(
    per_channel,
    file.path(out_dir, "repro_check.csv"),
    row.names = FALSE
  )
  message(
    "repro188 max abs diff over all channels: ",
    max(per_channel$max_abs_diff)
  )
}
stopifnot("all files verified" = all(verify$ok))
message("arm ", arm, " OK in ", round(as.numeric(runtime), 1), " min")
