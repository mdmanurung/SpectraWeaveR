# Purpose: unmix the RAW single-stain controls (35 files) and the unstained AF file with the exact spectra,
#          AF spectra, variants, AutoSpectral version and settings of one of 10 benchmark arms, so the
#          controls can be compared across arms.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/16_unmix_controls.R <arm> [smoke_out_dir]
#          from the repo root. <arm> is one of the names in ARMS. With smoke_out_dir, only SMOKE_FILES are
#          unmixed and the output goes to that folder instead of unmix_controls/<arm>.
# Settings come from (nothing invented; each arm lists its source in ARMS$source):
#   as102_basic     = 06_unmix_c102.R cell s102 (folder unmix_cells/repro_s102_c102), WITHOUT variants
#                     (README: "no variants"; ReRunUnmix.R first call). AMBIGUOUS: the repro folder's args.txt
#                     lists variants; see the report to the user.
#   as102_var       = the same, with the 1.0.2 per-cell variants (06_unmix_c102.R cell s102; speed = 1.0.2 default "slow")
#   as188_joint     = settings of unmix_cells/repro_s188_c188_fast (06_unmix_c188.R cell s188, speed fast, joint);
#                     identical to as188_fast_csv
#   as188_fast_csv  = unmix_cells/repro_s188_c188_fast
#   s102_c188       = 06_unmix_c188.R cell s102; args.txt has no speed line, so the 1.8.8 default "fast" is used
#   s188_c102       = 06_unmix_c102.R cell s188_20261105_143803 (unmix_cells/s188_143803_c102)
#   best102_neg/af  = 07_best_c102.R; the saved af_spectra.rds and spectral_variants.rds are reloaded, not rebuilt
#   best188_p1/p5   = 07_best_c188.R; the saved af_spectra.rds is reloaded, n.passes 1 / 5
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group/*.fcs (35 raw controls)
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group/unstained_AF/*.fcs (1 raw file)
#          spectra, AF, variants, control files and lib_as102 as in the 06/07 scripts [user zip]
#          experiments/bmv_test/compare_3D/data/unmix_cells/best*/ (af_spectra.rds, spectral_variants.rds)
# Outputs: experiments/bmv_test/compare_3D/data/unmix_controls/<arm>/*.fcs unmixed controls
#          <arm>/args.txt inputs and settings, <arm>/verify.csv event counts and channel names vs the arm's
#          unmixed samples, <arm>/sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
LIB_102 <- file.path(DATA, "lib_as102") # AutoSpectral and AutoSpectralRcpp v1.0.2 (separate library)
OUT_ROOT <- file.path(DATA, "unmix_controls") # git-ignored data dir, one folder per arm
SMOKE_FILES <- c("Zombie NIR (Cells).fcs", "PE (Beads).fcs") # one cell and one bead control
WEIGHTED_102 <- FALSE # as in 06_unmix_c102.R and 07_best_c102.R
PIPELINE_188 <- "joint" # as in 06_unmix_c188.R and 07_best_c188.R
SPECTRA_102 <- "Cells Matching Negative_autospectral_spectra.csv"
AF_102 <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv"
SPECTRA_188 <- "Automated_autospectral_spectra_20261105_143803.csv"
AF_188 <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"

# Arm table -------------------------------------------------------------------------------
d102 <- file.path(DATA, "autoSpectral_v1.0.2")
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
cells <- file.path(DATA, "unmix_cells")
ctrl_102 <- file.path(d102, "Reference Group")
ctrl_188 <- file.path(d188, "Reference Group Clean")
var_102 <- file.path(
  ctrl_102,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)
var_188 <- file.path(d188, "figure_spectral_variants", "Spectral_variants.rds")

# code = AutoSpectral version; pipeline / n_passes NA = not passed (package default); speed always slow
arm_row <- function(
  code,
  ctrl,
  spectra_dir,
  spectra,
  af_csv = NA,
  af_rds = NA,
  variants = NA,
  speed = "slow", # always slow (user, 2026-11-10); 1.0.2 default is slow, joint ignores speed
  n_passes = NA,
  ref_dir,
  setwd_ctrl = FALSE
) {
  list(
    code = code,
    ctrl = ctrl,
    spectra_dir = spectra_dir,
    spectra = spectra,
    af_csv = af_csv,
    af_rds = af_rds,
    variants = variants,
    speed = speed,
    n_passes = n_passes,
    ref_dir = ref_dir,
    setwd_ctrl = setwd_ctrl
  )
}
tab_102 <- file.path(ctrl_102, "table_spectra")
tab_188 <- file.path(d188, "table_spectra")
ARMS <- list(
  as102_basic = arm_row(
    "1.0.2",
    ctrl_102,
    tab_102,
    SPECTRA_102,
    af_csv = AF_102,
    ref_dir = file.path(d102, "AutospectralUnmix")
  ),
  as102_var = arm_row(
    "1.0.2",
    ctrl_102,
    tab_102,
    SPECTRA_102,
    af_csv = AF_102,
    variants = var_102,
    ref_dir = file.path(d102, "AutospectralUnmix_optimized")
  ),
  as188_joint = arm_row(
    "1.8.8",
    ctrl_188,
    tab_188,
    SPECTRA_188,
    af_csv = AF_188,
    variants = var_188,
    speed = "slow",
    ref_dir = file.path(d188, "AutoSpectral_unmixed")
  ),
  as188_fast_csv = arm_row(
    "1.8.8",
    ctrl_188,
    tab_188,
    SPECTRA_188,
    af_csv = AF_188,
    variants = var_188,
    speed = "slow",
    ref_dir = file.path(cells, "repro_s188_c188_fast")
  ),
  s102_c188 = arm_row(
    "1.8.8",
    ctrl_102,
    tab_102,
    SPECTRA_102,
    af_csv = AF_102,
    variants = var_102,
    speed = "slow",
    ref_dir = file.path(cells, "s102_c188")
  ),
  s188_c102 = arm_row(
    "1.0.2",
    ctrl_188,
    tab_188,
    SPECTRA_188,
    af_csv = AF_188,
    variants = var_188,
    ref_dir = file.path(cells, "s188_143803_c102")
  ),
  best102_neg = arm_row(
    "1.0.2",
    ctrl_102,
    tab_102,
    SPECTRA_102,
    af_rds = file.path(cells, "best102_neg", "af_spectra.rds"),
    variants = file.path(cells, "best102_neg", "spectral_variants.rds"),
    ref_dir = file.path(cells, "best102_neg"),
    setwd_ctrl = TRUE
  ),
  best102_af = arm_row(
    "1.0.2",
    ctrl_102,
    tab_102,
    "Cells AF clean up_autospectral_spectra.csv",
    af_rds = file.path(cells, "best102_af", "af_spectra.rds"),
    variants = file.path(cells, "best102_af", "spectral_variants.rds"),
    ref_dir = file.path(cells, "best102_af"),
    setwd_ctrl = TRUE
  ),
  best188_p1 = arm_row(
    "1.8.8",
    ctrl_188,
    tab_188,
    SPECTRA_188,
    af_rds = file.path(cells, "best188_p1", "af_spectra.rds"),
    variants = var_188,
    n_passes = 1L,
    ref_dir = file.path(cells, "best188_p1")
  ),
  best188_p5 = arm_row(
    "1.8.8",
    ctrl_188,
    tab_188,
    SPECTRA_188,
    af_rds = file.path(cells, "best188_p5", "af_spectra.rds"),
    variants = var_188,
    n_passes = 5L,
    ref_dir = file.path(cells, "best188_p5")
  )
)

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
arm <- args[[1]]
smoke_dir <- if (length(args) >= 2) args[[2]] else NA
stopifnot("arm is one of ARMS" = arm %in% names(ARMS))
a <- ARMS[[arm]]
out_dir <- if (is.na(smoke_dir)) file.path(OUT_ROOT, arm) else smoke_dir
stopifnot("out dir must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
out_dir <- normalizePath(out_dir)
if (a$code == "1.0.2") {
  .libPaths(c(LIB_102, .libPaths()))
}
suppressPackageStartupMessages({
  library(AutoSpectral)
  library(dplyr)
  library(purrr)
  library(tibble)
})
if (a$code == "1.0.2") {
  suppressPackageStartupMessages(library(AutoSpectralRcpp))
}
stopifnot(
  "requested AutoSpectral version loaded" = as.character(packageVersion(
    "AutoSpectral"
  )) ==
    a$code
)

# Inputs ----------------------------------------------------------------------------------
raw_dir <- file.path(d188, "Reference Group")
raw_files <- c(
  list.files(raw_dir, pattern = "\\.fcs$", full.names = TRUE),
  list.files(
    file.path(raw_dir, "unstained_AF"),
    pattern = "\\.fcs$",
    full.names = TRUE
  )
)
if (!is.na(smoke_dir)) {
  raw_files <- raw_files[basename(raw_files) %in% SMOKE_FILES]
}
# unmix.folder takes a folder, so link the raw files into one staging folder
stage_dir <- file.path(tempdir(), "controls_raw")
dir.create(stage_dir)
stopifnot("control file names are unique" = !anyDuplicated(basename(raw_files)))
file.symlink(raw_files, file.path(stage_dir, basename(raw_files)))

if (a$setwd_ctrl) {
  setwd(a$ctrl)
} # 07_best_c102.R does this
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(
  a$ctrl,
  file.path(a$ctrl, "fcs_control_file.csv"),
  asp
)
spectra <- read.spectra(a$spectra, spectra.dir = a$spectra_dir)
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in the 06/07 scripts
af.spectra <- if (is.na(a$af_rds)) {
  read.spectra(a$af_csv, spectra.dir = a$spectra_dir)
} else {
  readRDS(a$af_rds)
}
spectra.variants <- if (is.na(a$variants)) NULL else readRDS(a$variants)

# Unmix -----------------------------------------------------------------------------------
# settings the arm's script passes; speed / n.passes only when the arm sets them
unmix_args <- list(
  fcs.dir = stage_dir,
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af.spectra,
  spectra.variants = spectra.variants,
  output.dir = out_dir,
  parallel = TRUE
)
if (a$code == "1.0.2") {
  unmix_args <- c(unmix_args, list(weighted = WEIGHTED_102, threads = NULL))
} else {
  unmix_args <- c(unmix_args, list(pipeline = PIPELINE_188, verbose = TRUE))
}
if (!is.na(a$speed)) {
  unmix_args$speed <- a$speed
}
if (!is.na(a$n_passes)) {
  unmix_args$n.passes <- a$n_passes
}

writeLines(
  c(
    paste("arm:", arm),
    paste("AutoSpectral:", a$code),
    paste("controls unmixed (n):", length(raw_files)),
    paste("controls:", raw_files),
    paste("flow.control from:", a$ctrl),
    paste("spectra:", file.path(a$spectra_dir, a$spectra)),
    paste(
      "af:",
      if (is.na(a$af_rds)) file.path(a$spectra_dir, a$af_csv) else a$af_rds
    ),
    paste("variants:", if (is.na(a$variants)) "none" else a$variants),
    paste("speed:", if (is.na(a$speed)) "package default" else a$speed),
    paste(
      "n.passes:",
      if (is.na(a$n_passes)) "package default" else a$n_passes
    ),
    paste("weighted:", if (a$code == "1.0.2") WEIGHTED_102 else "n/a"),
    paste("pipeline:", if (a$code == "1.0.2") "n/a" else PIPELINE_188),
    paste("setwd to control folder:", a$setwd_ctrl),
    "parallel: TRUE, threads: NULL where the arm's script sets it"
  ),
  file.path(out_dir, "args.txt")
)
t0 <- Sys.time()
do.call(unmix.folder, unmix_args)
cat(
  "elapsed_min:",
  format(round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)),
  "\n",
  file = file.path(out_dir, "args.txt"),
  append = TRUE
)

# Verify against the raw inputs and the arm's unmixed samples ------------------------------
ref_file <- list.files(
  a$ref_dir,
  pattern = "AutoSpectral.*\\.fcs$",
  full.names = TRUE
)[[1]]
fcs_names <- function(f) {
  colnames(flowCore::read.FCS(f, truncate_max_range = FALSE, which.lines = 1))
}
fcs_n <- function(f) as.integer(flowCore::read.FCSheader(f)[[1]][["$TOT"]])
ref_names <- fcs_names(ref_file)
out_files <- list.files(out_dir, pattern = "\\.fcs$", full.names = TRUE)
verify <- tibble(raw = raw_files) |>
  mutate(
    # 1.0.2 names the output after the $FIL keyword, 1.8.8 after the file on disk: accept either
    out = map_chr(raw, \(x) {
      stems <- c(
        tools::file_path_sans_ext(basename(x)),
        tools::file_path_sans_ext(flowCore::read.FCSheader(x)[[1]][["$FIL"]])
      )
      paste(
        out_files[basename(out_files) %in% paste(stems, "AutoSpectral.fcs")],
        collapse = ";"
      )
    }),
    n_raw = map_int(raw, fcs_n),
    n_out = map_int(out, \(x) if (x == "") NA_integer_ else fcs_n(x)),
    same_names = map_lgl(out, \(x) {
      if (x == "") NA else identical(fcs_names(x), ref_names)
    }),
    ref_file = ref_file
  )
write.csv(verify, file.path(out_dir, "verify.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
