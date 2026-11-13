# Purpose: tune AutoSpectral 1.8.8 n.passes and n.af.passes on a 3 x 3 grid (1, 5, 10 each; user,
#          2026-11-10). Every arm uses the zip 1.8.8 objects exactly as as188_joint (run 143803 dye spectra,
#          raw-csv AF, zip variants) and differs only in the two pass counts. Scoring is 51_score_ablation.R;
#          the grid figure is 52b_plot_passes_grid.R.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/52_passes_grid.R <arm> [condition]
#            from the repo root; arm = p<n.passes>_af<n.af.passes>, e.g. p5_af10.
#            Output goes to <data_dir>/unmix_cells/<arm>[/<condition>], which must not exist.
#            condition = last argument, required when the dataset has conditions (24H: one task per condition),
#            refused otherwise. SW_RUN_DATE is only read by scripts that write a run folder.
#            SLURM: sbatch --array=1-9 experiments/slurm/run_passes_grid.sbatch
# Inputs:  <data_dir>/autoSpectral_v1.8.8 (config-stamp dye spectra and AF, variants,
#            Reference Group Clean controls, Raw[/<condition>], AutoSpectral_unmixed[/<condition>]) [user zip]
#          stamp file names (spectra, AF) from config stamps; samples from datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <out_dir>/*.fcs      unmixed files
#          <out_dir>/args.txt   what was run, incl. unmix minutes
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

# Dataset config --------------------------------------------------------------------------
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath(data_path())
PIPELINE <- "joint" # as as188_joint (inherited: user zip AutoSpectral_unmix_v1.8.8.R)
SPEED <- "slow" # user rule, 2026-11-10: always slow; no effect under joint (repro188 = s188_c188_slow)
PASS_GRID <- c(1, 5, 10) # n.passes and n.af.passes values (user, 2026-11-10)
ARMS <- tidyr::expand_grid(n_passes = PASS_GRID, n_af_passes = PASS_GRID) |>
  mutate(arm = paste0("p", n_passes, "_af", n_af_passes))

# Setup -----------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
arm_id <- cli_args[[1]]
settings <- dplyr::filter(ARMS, arm == arm_id)
stopifnot("arm is one of the grid arms" = nrow(settings) == 1)
# 24H: one task per condition subfolder (argument 2); 3D and EXV have no condition level
cond <- if (length(cli_args) >= 2) cli_args[[2]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "52_passes_grid.R: dataset ",
    DS$dataset,
    " has no conditions; drop argument 2 ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "52_passes_grid.R: needs argument 2, a condition of dataset ",
    DS$dataset,
    ": ",
    paste(DS$conditions, collapse = ", "),
    call. = FALSE
  )
}
samples <- samples_tbl()
if (!is.null(cond)) {
  samples <- filter(samples, .data$condition == .env$cond)
}
out_dir <- paste(
  c(file.path(DATA, "unmix_cells", arm_id), cond),
  collapse = "/"
)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)

# Inputs: the zip 1.8.8 objects, as as188_joint ------------------------------------------------
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
control_dir <- file.path(d188, "Reference Group Clean")
spectra_file <- file.path(d188, "table_spectra", cfg("stamps")$as188_spectra)
af_file <- file.path(d188, "table_spectra", cfg("stamps")$as188_af)
variants_rds <- file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)

# Raw csv read: read.spectra() rescales rows and caused an offset (user, 2026-11-10)
spectra <- as.matrix(read.csv(spectra_file, row.names = 1, check.names = FALSE))
af_spectra <- as.matrix(read.csv(af_file, row.names = 1, check.names = FALSE))
spectra_variants <- readRDS(variants_rds)
stopifnot(
  "AF detectors equal dye spectra detectors" = identical(
    colnames(af_spectra),
    colnames(spectra)
  )
)

# Unmix -----------------------------------------------------------------------------------
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(
  control_dir,
  file.path(control_dir, "fcs_control_file.csv"),
  asp
)
t0 <- Sys.time()
unmix.folder(
  fcs.dir = paste(c(file.path(d188, "Raw"), cond), collapse = "/"),
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af_spectra,
  spectra.variants = spectra_variants,
  pipeline = PIPELINE,
  speed = SPEED,
  n.passes = settings$n_passes,
  n.af.passes = settings$n_af_passes,
  output.dir = out_dir,
  parallel = TRUE,
  verbose = TRUE
)
unmix_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
message("unmix minutes: ", unmix_minutes)

writeLines(
  c(
    paste("arm:", arm_id),
    paste("n.passes:", settings$n_passes),
    paste("n.af.passes:", settings$n_af_passes),
    paste("speed:", SPEED),
    paste("pipeline:", PIPELINE),
    paste("spectra:", spectra_file),
    paste("af:", af_file),
    paste("variants:", variants_rds),
    paste("controls:", control_dir),
    "csv read: raw as.matrix(read.csv(row.names = 1, check.names = FALSE))",
    paste("unmix minutes:", unmix_minutes)
  ),
  file.path(out_dir, "args.txt")
)

# Verify: same files, event counts and channel names as as188_joint -----------------------------
ref_dir <- paste(
  c(file.path(d188, "AutoSpectral_unmixed"), cond),
  collapse = "/"
)
read_hdr <- \(path) {
  ff <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  list(n = nrow(ff), channels = colnames(ff))
}
# one unmixed file per sample: the stem is the delivered file name without its " AutoSpectral*.fcs" suffix (EXV: _v1-8-8);
# the join stops on any other file
verify <- tibble(
  ref_file = list.files(ref_dir, pattern = "\\.fcs$", recursive = TRUE)
) |>
  mutate(stem = sub(" AutoSpectral[^ ]*\\.fcs$", "", ref_file)) |>
  inner_join(
    select(samples, stem),
    by = "stem",
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  mutate(
    file = paste0(stem, " AutoSpectral.fcs"), # the arm's own file name
    ref = map(file.path(ref_dir, ref_file), read_hdr),
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
  "one unmixed file per sample" = nrow(verify) == nrow(samples),
  "event counts and channels equal as188_joint" = all(
    verify$n_ref == verify$n_new & verify$channels_equal
  )
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
