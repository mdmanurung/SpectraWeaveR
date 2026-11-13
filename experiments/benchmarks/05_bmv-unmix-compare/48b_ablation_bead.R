# Purpose: regenerate the zip AutoSpectral 1.8.8 spectra set from the controls and unmix, for two arms
#          (the unmix runs on the dataset's raw samples).
#          bead_wide: bead pre-gate widened to FSC 2e5-9e5, SSC 1.3e6-2.25e6 (user, 2026-11-10, README O1);
#            the zip gate (FSC 3e5-8e5, SSC 1.3e6-1.9e6) clipped the brightest stained beads.
#          regen188: zip bead gate, everything else identical; the nuisance control for bead_wide, because
#            bead_wide is regenerated and as188_joint is not.
#          Everything else follows AutoSpectral_unmix_v1.8.8.R: cell rectangle FSC 6e5-4e6, SSC 1e5-1.6e6
#          (cell files copied from Reference Group Clean, including the 84,528-event Zombie NIR file, user
#          D13), get.spectra.automated(10000, 500, 0.95), get.af.spectra(refine = FALSE) on H12,
#          get.spectral.variants() defaults, unmix.folder(joint). Scoring is in 51_score_ablation.R.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/48b_ablation_bead.R <arm> [condition]
#            from the repo root; arm = bead_wide or regen188. Output goes to
#            <data_dir>/unmix_cells/<arm>[/<condition>], which must not exist.
#            condition = last argument, required when the dataset has conditions (24H: one task per condition;
#            the spectra are regenerated in each condition folder), refused otherwise.
#            SW_RUN_DATE is only read by scripts that write a run folder.
# Inputs:  <data_dir>/autoSpectral_v1.8.8/Reference Group (raw bead controls,
#            unstained_AF), Reference Group Clean (cell controls, fcs_control_file.csv), Raw[/<condition>],
#            AutoSpectral_unmixed[/<condition>], table_spectra (config-stamp spectra),
#            figure_spectral_variants/Spectral_variants.rds [user zip]
#          bead rectangles from config control_gates (bead_zip, bead_wide); stamp file name (spectra) from
#            config stamps; samples from datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <out_dir>/*.fcs                 unmixed files
#          <out_dir>/controls/             control files used (beads re-gated, cells copied)
#          <out_dir>/work/                 AutoSpectral tables and figures from regeneration
#          <out_dir>/spectra.csv, af_spectra.csv, variants.rds  regenerated objects
#          <out_dir>/bead_events.csv       bead events kept per control, this arm vs Reference Group Clean
#          <out_dir>/object_check.csv      regenerated objects vs the zip objects
#          <out_dir>/args.txt, verify.csv, sessionInfo.txt

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
PIPELINE <- "joint" # as in AutoSpectral_unmix_v1.8.8.R (inherited: user zip)
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint"
EXTRACT_SEED <- 42 # set.seed() around get.spectra.automated (user, 2026-11-10, D5); inside, asp$bird.seed
# bead rectangles are dataset config (3D: zip gate FSC 3e5-8e5, SSC 1.3e6-1.9e6, AutoSpectral_unmix_v1.8.8.R:62-65;
# bead_wide FSC 2e5-9e5, SSC 1.3e6-2.25e6, user 2026-11-10 O1)
BEAD_GATES <- list(
  regen188 = cfg("control_gates")$bead_zip,
  bead_wide = cfg("control_gates")$bead_wide
)

# Setup -----------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
arm <- cli_args[[1]]
stopifnot("arm is bead_wide or regen188" = arm %in% names(BEAD_GATES))
# 24H: one task per condition subfolder (argument 2); 3D and EXV have no condition level
cond <- if (length(cli_args) >= 2) cli_args[[2]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "48b_ablation_bead.R: dataset ",
    DS$dataset,
    " has no conditions; drop argument 2 ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "48b_ablation_bead.R: needs argument 2, a condition of dataset ",
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
out_dir <- paste(c(file.path(DATA, "unmix_cells", arm), cond), collapse = "/")
stopifnot("out dir must not exist" = !dir.exists(out_dir))
ctrl_dir <- file.path(out_dir, "controls")
work_dir <- file.path(out_dir, "work")
dir.create(ctrl_dir, recursive = TRUE)
dir.create(work_dir)
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
clean_dir <- file.path(d188, "Reference Group Clean")
raw_ctrl_dir <- file.path(d188, "Reference Group")

# Controls: cell files and control file copied, bead files re-gated from the raw beads -------------
cell_files <- list.files(clean_dir, pattern = "\\(Cells\\)\\.fcs$")
file.copy(file.path(clean_dir, cell_files), ctrl_dir)
file.copy(file.path(clean_dir, "fcs_control_file.csv"), ctrl_dir)
gate <- BEAD_GATES[[arm]]
bead_gate <- rectangleGate("FSC-A" = gate$fsc, "SSC-A" = gate$ssc)
bead_files <- list.files(raw_ctrl_dir, pattern = "\\(Beads\\)\\.fcs$")
bead_events <- bead_files |>
  map(\(f) {
    ff <- read.FCS(
      file.path(raw_ctrl_dir, f),
      transformation = FALSE,
      truncate_max_range = FALSE
    )
    kept <- Subset(ff, bead_gate)
    write.FCS(kept, file.path(ctrl_dir, f))
    zip_n <- nrow(read.FCS(
      file.path(clean_dir, f),
      transformation = FALSE,
      truncate_max_range = FALSE
    ))
    tibble(
      file = f,
      n_raw = nrow(ff),
      n_kept = nrow(kept),
      n_zip_clean = zip_n
    )
  }) |>
  list_rbind()
write_csv(bead_events, file.path(out_dir, "bead_events.csv"))

# Regenerate spectra, AF and variants as the zip script -----------------------------------
# AutoSpectral writes tables and figures to the asp *.dir fields; point them into work/
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
dir_fields <- names(asp)[grepl("dir", names(asp)) & map_lgl(asp, is.character)]
asp[dir_fields] <- map(asp[dir_fields], \(d) file.path(work_dir, d))
control_file <- file.path(ctrl_dir, "fcs_control_file.csv")
set.seed(EXTRACT_SEED)
spectra <- get.spectra.automated(
  control.dir = ctrl_dir,
  control.def.file = control_file,
  asp = asp,
  n.candidates = 10000,
  n.spectral = 500,
  cosine.threshold = 0.95
)
unstained <- list.files(
  file.path(raw_ctrl_dir, "unstained_AF"),
  pattern = "\\.fcs$",
  full.names = TRUE
)
stopifnot("one unstained AF file" = length(unstained) == 1)
af_spectra <- get.af.spectra(
  unstained.sample = unstained,
  asp = asp,
  spectra = spectra,
  refine = FALSE
)
spectra_variants <- get.spectral.variants(
  control.dir = ctrl_dir,
  control.def.file = control_file,
  asp = asp,
  spectra = spectra,
  output.dir = work_dir
)
write.csv(spectra, file.path(out_dir, "spectra.csv"))
write.csv(af_spectra, file.path(out_dir, "af_spectra.csv"))
saveRDS(spectra_variants, file.path(out_dir, "variants.rds"))

# Object check: regenerated objects vs the zip objects ---------------------------------------
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra_zip <- read_raw(file.path(
  d188,
  "table_spectra",
  cfg("stamps")$as188_spectra
))
variants_zip <- readRDS(file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
))
angle_deg <- \(a, b) {
  acos(min(1, sum(a * b) / sqrt(sum(a^2) * sum(b^2)))) * 180 / pi
}
stopifnot(
  "regenerated fluors equal zip fluors" = identical(
    rownames(spectra),
    rownames(spectra_zip)
  )
)
object_check <- tibble(dye = rownames(spectra)) |>
  mutate(
    spectra_angle_vs_zip = map_dbl(dye, \(d) {
      angle_deg(spectra[d, ], spectra_zip[d, ])
    }),
    n_variants = map_int(dye, \(d) nrow(spectra_variants$variants[[d]])),
    n_variants_zip = map_int(dye, \(d) nrow(variants_zip$variants[[d]])),
    threshold = unname(spectra_variants$thresholds[dye]),
    threshold_zip = unname(variants_zip$thresholds[dye]),
    variant_angle_median = map_dbl(dye, \(d) {
      median(apply(
        spectra_variants$variants[[d]],
        1,
        angle_deg,
        b = spectra_zip[d, ]
      ))
    }),
    variant_angle_median_zip = map_dbl(dye, \(d) {
      median(apply(
        variants_zip$variants[[d]],
        1,
        angle_deg,
        b = spectra_zip[d, ]
      ))
    }),
    af_rows = nrow(af_spectra)
  )
write_csv(object_check, file.path(out_dir, "object_check.csv"))

# Unmix -----------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste(
      "bead gate FSC-A:",
      paste(gate$fsc, collapse = "-"),
      " SSC-A:",
      paste(gate$ssc, collapse = "-")
    ),
    "cell gate: zip FSC-A 6e5-4e6, SSC-A 1e5-1.6e6 (cell files copied from Reference Group Clean)",
    "Zombie NIR: 84,528-event Reference Group Clean file (user, 2026-11-10, D13)",
    paste(
      "extraction seed (set.seed):",
      EXTRACT_SEED,
      " asp$bird.seed:",
      asp$bird.seed
    ),
    "get.spectra.automated(n.candidates = 10000, n.spectral = 500, cosine.threshold = 0.95)",
    "get.af.spectra(refine = FALSE); get.spectral.variants() defaults",
    paste("pipeline:", PIPELINE, " speed:", SPEED),
    paste("AF rows:", nrow(af_spectra))
  ),
  file.path(out_dir, "args.txt")
)
flow.control <- reload.flow.control(ctrl_dir, control_file, asp)
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
