# Purpose: test 3 (README decision table, user 2026-11-10): does the AutoSpectral 1.0.2 variant code, run
#          on the CLEAN controls (zip rectangle, 'Reference Group Clean'), give Alexa Fluor 405 (GLUT1)
#          variants with the BV510 spread of the zip 1.0.2 variants? Separates "1.0.2 variant code" from
#          "debris-heavy raw controls" as the source of the spread.
#          Step variants (1.0.2 library): get.spectral.variants exactly as the zip 1.0.2 setup called it
#            (Autospectral_Setup.R: spectra = 'Cells Matching Negative', AF = H12 unstained AF set, defaults
#            n.cells 2000, som.dim 10, parallel FALSE), but control.dir = the clean controls; seed 42.
#          Step unmix (1.8.8 library): measure the AF405 spread with the 74 metric against zip_102 and zip_188;
#            build arm var102code_clean_af405 = the 1.8.8 setup with only AF405's variants replaced (as in
#            45_variant_per_dye.R) and unmix the 5 raw files (joint, speed slow).
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/77_af405_variants_102code_clean.R <run_dir> variants
#          then the same with <step> = unmix (run dir must exist and hold the variants step's output).
#          run_dir = experiments/runs/2026-11-10_bmv-af405-102code-clean
#          Scoring: experiments/slurm/run_score_ablation.sbatch
#            experiments/runs/2026-11-10_bmv-af405-102code-clean-score var102code_clean_af405
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/  clean controls +
#            fcs_control_file.csv [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/table_spectra/
#            'Cells Matching Negative' spectra and H12 AF csv; figure_spectral_variants/Spectral_variants.rds
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/ spectra 143803, AF 143933, variants, Raw
#          experiments/bmv_test/compare_3D/data/lib_as102  AutoSpectral 1.0.2 library (as in 06_unmix_c102.R)
# Outputs: <run_dir>/variants_102code_clean.rds       1.0.2-code variants object (variants step)
#          <run_dir>/figure_spectral_variants/         1.0.2 variant figures (variants step)
#          <run_dir>/args_variants.txt                 what the variants step ran
#          <run_dir>/af405_spread.csv                  per set and variant: BV510 coefficient (unmix step)
#          <run_dir>/af405_spread_summary.csv          per set: n, median, 90th percentile
#          <run_dir>/fig_af405_spread.png (+ _source.csv)
#          experiments/bmv_test/compare_3D/data/unmix_cells/var102code_clean_af405/  unmixed FCS, variants_used.rds,
#            args.txt, verify.csv, sessionInfo.txt
#          <run_dir>/sessionInfo_<step>.txt

args <- commandArgs(trailingOnly = TRUE)
RUN_DIR <- args[[1]]
STEP <- args[[2]]
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
if (STEP == "variants") {
  .libPaths(c(file.path(DATA, "lib_as102"), .libPaths()))
}

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ---------------------------------------------------------------------------------------
DYE <- "Alexa Fluor 405"
SEED <- 42 # README seed rule (bird.seed + 42)
PIPELINE <- "joint" # as 45_variant_per_dye.R (inherited: user zip)
SPEED <- "slow" # always slow (user, 2026-11-10)
VARIANT_SLOTS <- c("variants", "delta.list", "delta.norms") # per-dye lists swapped together, as in 45
ARM <- "var102code_clean_af405"

# Inputs -------------------------------------------------------------------------------------------
d102 <- file.path(DATA, "autoSpectral_v1.0.2", "Reference Group")
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
clean_dir <- file.path(d188, "Reference Group Clean")
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
} # no rescaling
s102_file <- file.path(
  d102,
  "table_spectra",
  "Cells Matching Negative_autospectral_spectra.csv"
)
af102_file <- file.path(
  d102,
  "table_spectra",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv"
)

# Step variants: 1.0.2 get.spectral.variants on the clean controls -----------------------------------
if (STEP == "variants") {
  stopifnot(
    "1.0.2 code loaded" = as.character(packageVersion("AutoSpectral")) ==
      "1.0.2"
  )
  stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
  dir.create(RUN_DIR, recursive = TRUE)
  RUN_DIR <- normalizePath(RUN_DIR)
  s102 <- read_raw(s102_file)
  af102 <- read_raw(af102_file)
  writeLines(
    c(
      "code: AutoSpectral 1.0.2 (lib_as102)",
      paste("control.dir:", clean_dir),
      "control.def.file: fcs_control_file.csv in control.dir",
      paste("spectra:", s102_file),
      paste("af.spectra:", af102_file, "rows", nrow(af102)),
      "n.cells 2000, som.dim 10, parallel FALSE (zip defaults)",
      paste("seed:", SEED)
    ),
    file.path(RUN_DIR, "args_variants.txt")
  )
  asp <- get.autospectral.param(cytometer = "aurora", figures = TRUE)
  set.seed(SEED)
  # the zip call used relative paths from inside control.dir; the function reads control.def.file as given
  old_wd <- setwd(clean_dir)
  variants <- get.spectral.variants(
    clean_dir,
    "fcs_control_file.csv",
    asp,
    s102,
    af102,
    output.dir = file.path(RUN_DIR, "figure_spectral_variants"),
    parallel = FALSE,
    verbose = TRUE
  )
  setwd(old_wd)
  saveRDS(variants, file.path(RUN_DIR, "variants_102code_clean.rds"))
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_variants.txt")
  )
}

# Step unmix: measure the AF405 spread, then build and unmix the swap arm ------------------------------
if (STEP == "unmix") {
  stopifnot(
    "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) ==
      "1.8.8"
  )
  s188 <- read_raw(file.path(
    d188,
    "table_spectra",
    "Automated_autospectral_spectra_20261105_143803.csv"
  ))
  s102 <- read_raw(s102_file)
  var188 <- readRDS(file.path(
    d188,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  ))
  var102 <- readRDS(file.path(
    d102,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  ))
  var_new <- readRDS(file.path(RUN_DIR, "variants_102code_clean.rds"))
  b510 <- s188["BV510", ]

  # BV510 coefficient of each AF405 variant minus its own set's reference (script 74 metric)
  proj_coef <- \(d, b) as.vector(d %*% b) / sum(b^2)
  sets <- list(
    zip_188 = list(v = var188$variants[[DYE]], ref = s188[DYE, ]),
    zip_102 = list(v = var102$variants[[DYE]], ref = s102[DYE, ]),
    code102_clean = list(v = var_new$variants[[DYE]], ref = s102[DYE, ])
  )
  spread <- imap(sets, \(x, set) {
    v <- x$v[, names(b510)]
    tibble(
      variant = seq_len(nrow(v)),
      spread_coef = proj_coef(sweep(v, 2, x$ref[names(b510)]), b510)
    )
  }) |>
    list_rbind(names_to = "set") |>
    mutate(set = factor(set, levels = names(sets)))
  write_csv(spread, file.path(RUN_DIR, "af405_spread.csv"))
  write_csv(
    summarise(
      spread,
      n_variants = n(),
      median = median(spread_coef),
      p90 = quantile(spread_coef, 0.9),
      .by = set
    ),
    file.path(RUN_DIR, "af405_spread_summary.csv")
  )
  p <- ggplot(spread, aes(set, spread_coef)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_point(
      position = position_jitter(width = 0.2, height = 0, seed = 1),
      size = 0.8,
      alpha = 0.6
    ) +
    geom_boxplot(outlier.shape = NA, fill = NA, width = 0.5) +
    labs(
      x = NULL,
      y = "BV510 coefficient of (variant - own reference)",
      title = "Alexa Fluor 405 variants: zip 1.8.8, zip 1.0.2, 1.0.2 code on clean controls"
    ) +
    theme_bw(base_size = 8)
  ggsave(
    file.path(RUN_DIR, "fig_af405_spread.png"),
    p,
    width = 4.5,
    height = 3.5,
    dpi = 150,
    bg = "white"
  )
  write_csv(spread, file.path(RUN_DIR, "fig_af405_spread_source.csv"))

  # swap arm: 1.8.8 variants object with AF405's three per-dye lists from the new 1.0.2-code object
  out_dir <- file.path(DATA, "unmix_cells", ARM)
  stopifnot("arm out dir must not exist" = !dir.exists(out_dir))
  dir.create(out_dir)
  spectra_variants <- modifyList(
    var188,
    map(set_names(VARIANT_SLOTS), \(slot) {
      out <- var188[[slot]]
      out[[DYE]] <- var_new[[slot]][[DYE]]
      out
    })
  )
  saveRDS(spectra_variants, file.path(out_dir, "variants_used.rds"))
  af188_file <- file.path(
    d188,
    "table_spectra",
    "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
  )
  writeLines(
    c(
      paste("arm:", ARM),
      paste("dye with 1.0.2-code clean-control variants:", DYE),
      paste(
        "variants source:",
        file.path(RUN_DIR, "variants_102code_clean.rds")
      ),
      paste("speed:", SPEED),
      paste("pipeline:", PIPELINE),
      "spectra: 1.8.8 143803; af: 1.8.8 143933; other variants: 1.8.8 zip",
      paste("controls:", clean_dir)
    ),
    file.path(out_dir, "args.txt")
  )
  asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
  flow.control <- reload.flow.control(
    clean_dir,
    file.path(clean_dir, "fcs_control_file.csv"),
    asp
  )
  unmix.folder(
    fcs.dir = file.path(d188, "Raw"),
    spectra = s188,
    asp = asp,
    flow.control = flow.control,
    method = "AutoSpectral",
    af.spectra = read_raw(af188_file),
    spectra.variants = spectra_variants,
    pipeline = PIPELINE,
    speed = SPEED,
    output.dir = out_dir,
    parallel = TRUE,
    verbose = TRUE
  )

  # verify: same files, event counts and channel names as as188_joint
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
    "5 files, counts and channels equal as188_joint" = nrow(verify) == 5 &&
      all(verify$n_ref == verify$n_new) &&
      all(verify$channels_equal)
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_unmix.txt")
  )
}
