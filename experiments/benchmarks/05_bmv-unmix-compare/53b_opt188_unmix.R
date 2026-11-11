# Purpose: unmix one opt188 arm from the objects one route's extraction already built (user, 2026-11-10,
#          grill Q8/Q9): spectra, AF spectra and variants are extracted once per route by 53_opt188.R
#          (stage extract), and every arm of that route unmixes from them, so arms differ only in
#          pipeline and passes. Legacy takes no passes (AutoSpectralRcpp sends n.passes and n.af.passes
#          only to the joint call), so legacy arms carry no passes in their name.
# Run:     sbatch experiments/slurm/run_opt188_unmix.sbatch <route> <pipeline> [n_passes n_af_passes]
#          route = pregate | legacyspec; pipeline = joint (needs n_passes, n_af_passes) | legacy
# Inputs:  experiments/runs/2026-11-10_bmv-opt188-smoke/opt188_<route>_neg_same_fsc4e5_p1/
#            spectra.csv, af_spectra.csv, variants.rds, controls/ (gated control copies and control file)
#            [53_opt188.R stage extract, jobs 25992962 (pregate) and 25991442 (legacyspec)]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Raw/*.fcs   5 raw samples [user zip]
# Outputs: experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/*.fcs   unmixed samples
#          <arm>/args.txt (settings and the object folder), <arm>/verify.csv, <arm>/sessionInfo.txt

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
OBJECT_ROOT <- normalizePath("experiments/runs/2026-11-10_bmv-opt188-smoke")
NEG <- "neg_same_fsc4e5" # the only cell negative kept (user, 2026-11-10)
SPEED <- "slow" # always slow (user, 2026-11-10); under legacy it sets the variants tested per cell

# Setup -----------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
route <- cli_args[[1]]
pipeline <- cli_args[[2]]
stopifnot(
  "route is pregate or legacyspec" = route %in% c("pregate", "legacyspec"),
  "pipeline is joint or legacy" = pipeline %in% c("joint", "legacy"),
  "joint needs n_passes and n_af_passes; legacy takes none" = length(
    cli_args
  ) ==
    if (pipeline == "joint") 4 else 2,
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
out_dir <- file.path(DATA, "unmix_cells", arm)
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
  fcs.dir = file.path(d188, "Raw"),
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
