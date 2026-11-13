# Purpose: ablation Block A (user, 2026-11-10, README rows D1, D8). Object-swap arms that complete the
#          spectra x AF x variants cube between the zip 1.8.8 set and the 1.0.2 set, plus one arm that swaps
#          only the per-dye positivity thresholds of the variants object. Same setup as 39_variant_count.R:
#          1.8.8 code, joint pipeline (legacy for two arms), speed slow, Reference Group Clean
#          flow.control, raw csv loading.
#          Scoring is in 51_score_ablation.R.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/48_ablation_swap.R <arm> [condition]
#            from the repo root. arm = one of ARM_SOURCES. Output goes to
#            <data_dir>/unmix_cells/<arm>[/<condition>], which must not exist.
#            condition = last argument, required when the dataset has conditions (24H: one task per condition),
#            refused otherwise. SW_RUN_DATE is only read by scripts that write a run folder.
# Inputs:  <data_dir>/autoSpectral_v1.8.8 (config-stamp dye spectra and AF,
#            Spectral_variants.rds, Reference Group Clean, Raw[/<condition>], AutoSpectral_unmixed[/<condition>]) [user zip]
#          <data_dir>/autoSpectral_v1.0.2/Reference Group (config-stamp 1.0.2 spectra and AF csv,
#            Spectral_variants.rds) [user zip]
#          stamp file names (spectra, AF) from config stamps; samples from datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <out_dir>/*.fcs              unmixed files
#          <out_dir>/args.txt           what was run
#          <out_dir>/variants_used.rds  the variants object passed to unmix.folder()
#          <out_dir>/verify.csv         event count and channel names vs as188_joint, per file
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
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint" (repro188 vs s188_c188_slow)
# Which set each object comes from, per arm ("188" = zip 1.8.8, "102" = zip 1.0.2), and the unmixing
# pipeline: "joint" as in AutoSpectral_unmix_v1.8.8.R (user zip); "legacy" arms asked by the user
# (2026-11-10). The thresholds column swaps only the variants object's `thresholds` slot, the one
# extra slot the joint unmixer reads.
ARM_SOURCES <- tribble(
  ~arm                , ~spectra , ~af   , ~variants , ~thresholds , ~pipeline ,
  "swap_spec102"      , "102"    , "188" , "188"     , "188"       , "joint"   ,
  "swap_af102_var102" , "188"    , "102" , "102"     , "102"       , "joint"   ,
  "rev_var188"        , "102"    , "102" , "188"     , "188"       , "joint"   ,
  "rev_af188"         , "102"    , "188" , "102"     , "102"       , "joint"   ,
  "swap_thr102"       , "188"    , "188" , "188"     , "102"       , "joint"   ,
  "legacy188"         , "188"    , "188" , "188"     , "188"       , "legacy"  ,
  "legacy102in"       , "102"    , "102" , "102"     , "102"       , "legacy"
)

# Setup -----------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
arm_name <- cli_args[[1]]
src <- ARM_SOURCES |>
  dplyr::filter(arm == arm_name)
stopifnot("arm is one of ARM_SOURCES" = nrow(src) == 1)
# 24H: one task per condition subfolder (argument 2); 3D and EXV have no condition level
cond <- if (length(cli_args) >= 2) cli_args[[2]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "48_ablation_swap.R: dataset ",
    DS$dataset,
    " has no conditions; drop argument 2 ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "48_ablation_swap.R: needs argument 2, a condition of dataset ",
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
  c(file.path(DATA, "unmix_cells", arm_name), cond),
  collapse = "/"
)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)

# Inputs: both zip sets ------------------------------------------------------------------------
d102 <- file.path(DATA, "autoSpectral_v1.0.2", "Reference Group")
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
control_dir <- file.path(d188, "Reference Group Clean")
files <- list(
  spectra = c(
    "188" = file.path(d188, "table_spectra", cfg("stamps")$as188_spectra),
    "102" = file.path(d102, "table_spectra", cfg("stamps")$as102_spectra)
  ),
  af = c(
    "188" = file.path(d188, "table_spectra", cfg("stamps")$as188_af),
    "102" = file.path(d102, "table_spectra", cfg("stamps")$as102_af)
  ),
  variants = c(
    "188" = file.path(
      d188,
      "figure_spectral_variants",
      "Spectral_variants.rds"
    ),
    "102" = file.path(d102, "figure_spectral_variants", "Spectral_variants.rds")
  )
)

# Raw csv read: read.spectra() rescales rows and caused an offset (user, 2026-11-10)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra <- read_raw(files$spectra[[src$spectra]])
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # 1.0.2 file has an AF row, as ReRunUnmix.R
af_spectra <- read_raw(files$af[[src$af]])
spectra_variants <- readRDS(files$variants[[src$variants]])
thresholds_src <- readRDS(files$variants[[src$thresholds]])$thresholds
spectra_variants$thresholds <- thresholds_src[names(
  spectra_variants$thresholds
)]

# Fluor and detector sets must match; stop, do not coerce
stopifnot(
  "AF detectors equal dye spectra detectors" = identical(
    colnames(af_spectra),
    colnames(spectra)
  ),
  "variants fluors equal dye spectra fluors" = identical(
    names(spectra_variants$variants),
    rownames(spectra)
  ),
  "thresholds named for every fluor" = !anyNA(spectra_variants$thresholds)
)
saveRDS(spectra_variants, file.path(out_dir, "variants_used.rds"))

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm_name),
    paste("speed:", SPEED),
    paste("pipeline:", src$pipeline),
    paste(
      "spectra:",
      files$spectra[[src$spectra]],
      "(AF row dropped if present)"
    ),
    paste("af:", files$af[[src$af]]),
    paste("variants:", files$variants[[src$variants]]),
    paste("thresholds slot from:", files$variants[[src$thresholds]]),
    paste(
      paste(cfg("targets")$tnf_dye, "threshold:"),
      round(spectra_variants$thresholds[[cfg("targets")$tnf_dye]], 1)
    ),
    paste("AF rows:", nrow(af_spectra)),
    paste(
      "variants per dye:",
      paste(map_int(spectra_variants$variants, nrow), collapse = ",")
    ),
    paste("controls (flow.control):", control_dir),
    "csv read: raw as.matrix(read.csv(row.names = 1, check.names = FALSE))"
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
  fcs.dir = paste(c(file.path(d188, "Raw"), cond), collapse = "/"),
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af_spectra,
  spectra.variants = spectra_variants,
  pipeline = src$pipeline,
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
  "event counts equal as188_joint" = all(verify$n_ref == verify$n_new),
  "channel names equal as188_joint" = all(verify$channels_equal)
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
