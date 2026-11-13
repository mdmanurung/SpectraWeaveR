# Copy of 61_af405_proof_unmix.R for the common-gates spike (plans/common_gates_plan.md G5): the unmix is unchanged;
#          verify.csv covers the files of the given raw_dir (one donor in a smoke run), not the five of as188_joint.
# Purpose: unmix one proof-chain arm (user, 2026-11-10): the 1.8.8 zip setup (1.8.8 code, joint,
#          speed "slow" (user rule), run-143803 dye spectra, raw-csv 1.8.8 AF, 1.8.8 thresholds) with a
#          given variants object, on a given folder of raw FCS files (the zip Raw folder, or spiked copies).
# Run:     SW_DATASET=<3d|24h|exv> pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/61cg_af405_proof_unmix.R \
#            <arm> <variants_rds|188> [raw_dir|-] [condition]
#          from the repo root; on SLURM via experiments/slurm/run_spike_cg.sbatch. SW_DATASET is required (3d, 24h, exv):
#          it sets the data dir and the 1.8.8 spectra and AF names (cfg("stamps")). "188" = the zip 1.8.8
#          variants object unchanged. raw_dir defaults to the zip Raw folder ("-" = default). Output goes to
#          data_path("unmix_cells")/af405proof_<arm>, which must not exist.
#          Datasets with conditions (24H): <condition> (iRBC, SEB or uRBC) is required, one task per condition;
#          the default raw_dir is Raw/<condition> and the output goes to af405proof_<arm>/<condition>.
# Inputs:  data_path("autoSpectral_v1.8.8") (spectra, AF, variants, controls, Raw)
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
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Setup -----------------------------------------------------------------------------------
DATA <- normalizePath(data_path())
cli_args <- commandArgs(trailingOnly = TRUE)
arm <- cli_args[[1]]
variants_arg <- cli_args[[2]]
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
# datasets with conditions keep one raw and one unmixed folder per condition; one task per condition
condition <- if (is.null(DS$conditions)) NULL else cli_args[4]
if (!is.null(DS$conditions) && !condition %in% DS$conditions) {
  stop(
    "61cg_af405_proof_unmix.R: dataset ",
    DS$dataset,
    " needs a 4th argument <condition>, one of ",
    paste(DS$conditions, collapse = ", "),
    call. = FALSE
  )
}
raw_dir <- if (length(cli_args) >= 3 && cli_args[[3]] != "-") {
  normalizePath(cli_args[[3]])
} else if (is.null(condition)) {
  file.path(d188, "Raw")
} else {
  file.path(d188, "Raw", condition)
}
if (!is.null(condition) && basename(raw_dir) != condition) {
  stop(
    "61cg_af405_proof_unmix.R: raw_dir must end in the condition ",
    condition,
    "; got ",
    raw_dir,
    call. = FALSE
  )
}
out_dir <- file.path(DATA, "unmix_cells", paste0("af405proof_", arm))
if (!is.null(condition)) {
  out_dir <- file.path(out_dir, condition)
}
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
  cfg("stamps")$as188_spectra
)
af_file <- file.path(
  d188,
  "table_spectra",
  cfg("stamps")$as188_af
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
    if (!is.null(condition)) paste("condition:", condition),
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
ref_dir <- file.path(d188, "AutoSpectral_unmixed") # 24H: one subfolder per condition, found recursively
read_hdr <- \(path) {
  ff <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  list(n = nrow(ff), channels = colnames(ff))
}
# raw file -> registered sample (join on the full file stem) -> its as188_joint file, one per sample
raw_names <- list.files(raw_dir, pattern = "\\.fcs$")
raw_stems <- sub("\\.fcs$", "", raw_names)
stopifnot(
  "61cg: every raw file is a registered sample (samples_tbl stem)" = all(
    raw_stems %in% samples_tbl()$stem
  )
)
raw_files <- tibble(raw_file = raw_names, stem = raw_stems) |>
  inner_join(
    samples_tbl() |>
      filter(stem %in% raw_stems) |>
      select(sample, stem),
    by = "stem",
    relationship = "one-to-one",
    unmatched = "error"
  )
ref_files <- tibble(
  ref_path = list.files(
    ref_dir,
    pattern = " AutoSpectral\\.fcs$",
    recursive = TRUE,
    full.names = TRUE
  )
) |>
  mutate(stem = sub(" AutoSpectral\\.fcs$", "", basename(ref_path))) |>
  filter(stem %in% raw_files$stem)
stopifnot(
  "61cg: one as188_joint file per raw file under ref_dir" = setequal(
    ref_files$stem,
    raw_files$stem
  )
)
verify <- raw_files |>
  inner_join(
    ref_files,
    by = "stem",
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  mutate(file = sub("\\.fcs$", " AutoSpectral.fcs", raw_file)) |>
  mutate(
    ref = map(ref_path, read_hdr),
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
