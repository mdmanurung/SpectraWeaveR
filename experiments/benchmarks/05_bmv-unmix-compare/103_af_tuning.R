# Purpose: stage 1 of the AS 1.8.8 negative-cell tuning on 24H. One arm = one setting of AF refine, n.passes (= n.af.passes) and
#          the variants necessity threshold, extracted and unmixed from the zip's own controls, gates and control file
#          (Reference Group Clean, no 3D gates), speed slow, joint pipeline, on the smoke donor only (3 files).
#          Grid and guardrail approved 2026-11-14: refine {FALSE, TRUE} x passes {1, 5} x threshold {0.01, 0.001, 0.1}.
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> SW_TUNE_SCOPE=<smoke|all> sbatch experiments/slurm/run_af_tuning.sbatch <refine> <n_passes> <threshold>
# Inputs:  <data_dir>/autoSpectral_v1.8.8/Reference Group Clean/ + fcs_control_file.csv   cleaned controls and control file [zip]
#          <data_dir>/autoSpectral_v1.8.8/Reference Group/unstained_AF/                    the one unstained file for AF [zip]
#          <data_dir>/autoSpectral_v1.8.8/Raw/<condition>/<stem>.fcs                       raw files of the smoke donor [zip]
# Outputs: <root> = experiments/runs/<SW_RUN_DATE>_bmv24h-af-tuning (shared by the arms; the arm folder must not exist)
#          <root>/<arm>/spectra.csv, af_spectra.csv, variants.rds, args.txt, sessionInfo.txt
#          <root>/<arm>/*.fcs   unmixed smoke-donor files, "<stem> AutoSpectral.fcs"
#          arm = a188_r<refine>_p<passes>_t<threshold>, e.g. a188_rFALSE_p1_t0.01 (the zip's own settings, up to the seed)

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # cfg(), data_path(), samples_tbl()

# Parameters ------------------------------------------------------------------------------
EXTRACT_SEED <- 42 # set.seed() before extraction, as 53_opt188.R (user, 2026-11-10)
PIPELINE <- "joint" # zip script
SPEED <- "slow" # always slow (user)
cli <- commandArgs(trailingOnly = TRUE)
stopifnot("args: refine n_passes threshold" = length(cli) == 3)
REFINE <- as.logical(cli[[1]])
N_PASSES <- as.integer(cli[[2]])
THRESHOLD <- as.numeric(cli[[3]])
stopifnot("refine is TRUE or FALSE" = !is.na(REFINE), "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8")
arm <- sprintf("a188_r%s_p%d_t%s", REFINE, N_PASSES, format(THRESHOLD))
SCOPE <- Sys.getenv("SW_TUNE_SCOPE") # "smoke" = smoke donor (stage 1); "all" = all 30 files (confirmation, user 2026-11-14)
stopifnot("SW_TUNE_SCOPE is smoke or all" = SCOPE %in% c("smoke", "all"))
ROOT <- file.path("experiments/runs", paste0(Sys.getenv("SW_RUN_DATE"), "_bmv", cfg("run_tag"), c(smoke = "-af-tuning", all = "-af-confirm")[[SCOPE]]))
out_dir <- file.path(normalizePath(ROOT, mustWork = FALSE), arm)
stopifnot("arm folder must not exist" = !dir.exists(out_dir))
dir.create(file.path(out_dir, "work", "raw"), recursive = TRUE)

# Inputs ----------------------------------------------------------------------------------
d188 <- file.path(normalizePath(data_path()), "autoSpectral_v1.8.8")
clean_dir <- file.path(d188, "Reference Group Clean")
control_file <- file.path(clean_dir, "fcs_control_file.csv")
unstained <- list.files(file.path(d188, "Reference Group", "unstained_AF"), pattern = "\\.fcs$", full.names = TRUE)
stopifnot("one unstained AF file" = length(unstained) == 1)
smoke <- samples_tbl() |> dplyr::filter(SCOPE == "all" | .data$donor == cfg("smoke_donor"))
stopifnot("smoke donor has 3 files" = SCOPE == "all" || nrow(smoke) == 3)
raw_files <- file.path(d188, "Raw", smoke$condition, paste0(smoke$stem, ".fcs"))
stopifnot("smoke raw files exist" = all(file.exists(raw_files)))
file.symlink(raw_files, file.path(out_dir, "work", "raw", basename(raw_files)))
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
dir_fields <- names(asp)[grepl("dir", names(asp)) & map_lgl(asp, is.character)]
asp[dir_fields] <- map(asp[dir_fields], \(d) file.path(out_dir, "work", d))

# Extraction and unmixing --------------------------------------------------------------------
set.seed(EXTRACT_SEED)
spectra <- get.spectra.automated(
  control.dir = clean_dir, control.def.file = control_file, asp = asp,
  n.candidates = 10000, n.spectral = 500, cosine.threshold = 0.95
)
af_spectra <- get.af.spectra(unstained.sample = unstained, asp = asp, spectra = spectra, refine = REFINE)
spectra_variants <- get.spectral.variants(
  control.dir = clean_dir, control.def.file = control_file, asp = asp, spectra = spectra,
  output.dir = file.path(out_dir, "work"), optimize.necessity.threshold = THRESHOLD
)
write.csv(spectra, file.path(out_dir, "spectra.csv"))
write.csv(af_spectra, file.path(out_dir, "af_spectra.csv"))
saveRDS(spectra_variants, file.path(out_dir, "variants.rds"))
flow.control <- reload.flow.control(clean_dir, control_file, asp)
unmix.folder(
  fcs.dir = file.path(out_dir, "work", "raw"), spectra = spectra, asp = asp, flow.control = flow.control,
  method = "AutoSpectral", af.spectra = af_spectra, spectra.variants = spectra_variants, pipeline = PIPELINE,
  speed = SPEED, n.passes = N_PASSES, n.af.passes = N_PASSES, output.dir = out_dir, parallel = TRUE, verbose = TRUE
)
out_files <- list.files(out_dir, pattern = "\\.fcs$")
stopifnot("one unmixed file per sample" = length(out_files) == nrow(smoke))
writeLines(c(sprintf("arm: %s", arm), sprintf("refine=%s n.passes=n.af.passes=%d necessity.threshold=%s seed=%d speed=%s pipeline=%s",
  REFINE, N_PASSES, format(THRESHOLD), EXTRACT_SEED, SPEED, PIPELINE)), file.path(out_dir, "args.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
