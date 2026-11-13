# Purpose: unmix one opt188 arm from the objects one route's extraction already built (user, 2026-11-10,
#          grill Q8/Q9): spectra, AF spectra and variants are extracted once per route by 53_opt188.R
#          (stage extract), and every arm of that route unmixes from them, so arms differ only in
#          pipeline and passes. Legacy takes no passes (AutoSpectralRcpp sends n.passes and n.af.passes
#          only to the joint call), so legacy arms carry no passes in their name.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_opt188_unmix.sbatch <route> <pipeline> [n_passes n_af_passes] [condition]
#          route = pregate | legacyspec; pipeline = joint (needs n_passes, n_af_passes) | legacy
#          condition = last argument, required when the dataset has conditions (24H: one task per condition),
#          refused otherwise. SW_RUN_DATE is only read by scripts that write a run folder.
# Inputs:  run key opt188-smoke: opt188_<route>_neg_same_fsc4e5_p1/
#            spectra.csv, af_spectra.csv, variants.rds, controls/ (gated control copies and control file)
#            [53_opt188.R stage extract; register the run folder in datasets/<ds>_runs.csv first]
#          <data_dir>/autoSpectral_v1.8.8/Raw[/<condition>]/*.fcs, AutoSpectral_unmixed[/<condition>]/*.fcs   samples [user zip]
#          samples from datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <data_dir>/unmix_cells/<arm>[/<condition>]/*.fcs   unmixed samples
#          <arm>[/<condition>]/args.txt (settings and the object folder), verify.csv, sessionInfo.txt

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
OBJECT_ROOT <- normalizePath(run_in("opt188-smoke")) # objects built by 53_opt188.R stage extract
NEG <- "neg_same_fsc4e5" # the only cell negative kept (user, 2026-11-10)
SPEED <- "slow" # always slow (user, 2026-11-10); under legacy it sets the variants tested per cell

# Setup -----------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
route <- cli_args[[1]]
pipeline <- cli_args[[2]]
stopifnot(
  "route is pregate or legacyspec" = route %in% c("pregate", "legacyspec"),
  "pipeline is joint or legacy" = pipeline %in% c("joint", "legacy")
)
# 24H: one task per condition subfolder (the argument after the pass counts); 3D and EXV have no condition level
n_fixed <- if (pipeline == "joint") 4 else 2
cond <- if (length(cli_args) > n_fixed) cli_args[[n_fixed + 1]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "53b_opt188_unmix.R: dataset ",
    DS$dataset,
    " has no conditions; drop the argument after the pass counts ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "53b_opt188_unmix.R: needs a last argument, a condition of dataset ",
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
stopifnot(
  "joint needs n_passes and n_af_passes; legacy takes none" = length(
    cli_args
  ) ==
    n_fixed + length(cond),
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)
n_passes <- if (pipeline == "joint") as.integer(cli_args[[3]]) else 1L
n_af_passes <- if (pipeline == "joint") as.integer(cli_args[[4]]) else 1L
arm <- paste0(
  "opt188_",
  route,
  "_",
  pipeline,
  if (pipeline == "joint") paste0("_p", n_passes, "af", n_af_passes)
)
obj_dir <- file.path(OBJECT_ROOT, paste0("opt188_", route, "_", NEG, "_p1"))
out_dir <- paste(c(file.path(DATA, "unmix_cells", arm), cond), collapse = "/")
stopifnot("out dir must not exist" = !dir.exists(out_dir))
work_dir <- file.path(out_dir, "work")
dir.create(work_dir, recursive = TRUE)
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
dir_fields <- names(asp)[grepl("dir", names(asp)) & map_lgl(asp, is.character)]
asp[dir_fields] <- map(asp[dir_fields], \(d) file.path(work_dir, d))

# Objects from the route's extraction --------------------------------------------------------------
spectra <- as.matrix(read.csv(
  file.path(obj_dir, "spectra.csv"),
  row.names = 1,
  check.names = FALSE
))
af_spectra <- as.matrix(read.csv(
  file.path(obj_dir, "af_spectra.csv"),
  row.names = 1,
  check.names = FALSE
))
spectra_variants <- readRDS(file.path(obj_dir, "variants.rds"))
ctrl_dir <- file.path(obj_dir, "controls")
flow.control <- reload.flow.control(
  ctrl_dir,
  file.path(ctrl_dir, "fcs_control_file.csv"),
  asp
)

writeLines(
  c(
    paste("arm:", arm),
    paste("objects from:", obj_dir),
    paste("route:", route, " cell negative:", NEG),
    paste(
      "pipeline:",
      pipeline,
      " speed:",
      SPEED,
      " n.passes:",
      if (pipeline == "joint") n_passes else "not used by legacy",
      " n.af.passes:",
      if (pipeline == "joint") n_af_passes else "not used by legacy"
    )
  ),
  file.path(out_dir, "args.txt")
)

# Unmix ------------------------------------------------------------------------------------------
t1 <- Sys.time()
unmix.folder(
  fcs.dir = paste(c(file.path(d188, "Raw"), cond), collapse = "/"),
  spectra = spectra,
  asp = asp,
  flow.control = flow.control,
  method = "AutoSpectral",
  af.spectra = af_spectra,
  spectra.variants = spectra_variants,
  pipeline = pipeline,
  speed = SPEED,
  n.passes = n_passes,
  n.af.passes = n_af_passes,
  output.dir = out_dir,
  parallel = TRUE,
  verbose = TRUE
)
message(
  "unmix minutes: ",
  round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1)
)

# Verify against as188_joint: event counts and channel names -------------------------------------
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
