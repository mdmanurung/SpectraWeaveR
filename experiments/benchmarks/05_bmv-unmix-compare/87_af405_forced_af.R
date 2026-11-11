# Purpose: grill Q4/Q9 (README rows "Forced AF removal on AF405", user 2026-11-10): force the AF handling the
#          zip 1.8.8 run skipped and see whether it changes the AF405 (GLUT1) objects and the TNF / Th2 squeeze.
#            af405_afrefine  the zip 1.8.8 setup with get.af.spectra(refine = TRUE) instead of the zip's
#                            refine = FALSE; everything else as the zip script (AutoSpectral_unmix_v1.8.8.R):
#                            zip dye spectra (143803), zip variants, joint, Reference Group Clean controls.
#                            The zip built the AF405 spectrum (get.spectra.automated) and variants
#                            (get.spectral.variants takes no AF input in 1.8.8) without the AF set, so only the
#                            AF set used in the unmix changes; the AF405 objects stay the zip's.
#            af405_afremove  NOT BUILT: clean.controls() needs a define.flow.control() gating object, and the
#                            zip's automated route never made one; how the cleaned AF405 events reach
#                            get.spectra.automated / get.fluorophore.spectra and get.spectral.variants (which
#                            reads FCS files, not a flow.control) is not set by the zip. User choice pending.
# Run:     step build (SLURM: experiments/slurm/run_script.sbatch):
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/87_af405_forced_af.R build
#              -> experiments/bmv_test/compare_3D/data/unmix_cells/af405_afrefine (must not exist)
#          then score with 51: sbatch experiments/slurm/run_score_ablation.sbatch \
#              experiments/runs/2026-11-10_bmv-af405-forced-af-score af405_afrefine
#          step report (after the scorer): ... 87_af405_forced_af.R report
#              -> experiments/runs/2026-11-10_bmv-af405-forced-af (must not exist)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/  zip spectra 143803, AF 143933, variants,
#            Reference Group Clean, Reference Group/unstained_AF/H12 *.fcs, Raw, AutoSpectral_unmixed [user zip]
#          experiments/runs/2026-11-10_bmv-af405-forced-af-score/gap_closed.csv   [51_score_ablation.R]
#          experiments/runs/2026-11-10_bmv-ablation-gated/gs/as188_joint, 2026-11-10_bmv-gated-pooled/manifest.csv
#            [19a, 12_cleanup.R]; helper75_marker-scores.R; lib/theme.R, lib/palettes.R, lib/plots.R
# Outputs: build: <out_dir>/*.fcs unmixed files; af_refined.csv the AF set; args.txt; verify.csv;
#            figure_autofluorescence/ (get.af.spectra plots); sessionInfo.txt
#          report: <report_dir>/af_bv510_coef.csv        BV510 coefficient of every AF spectrum, zip vs refined
#                  <report_dir>/af_bv510_coef.pdf/.png/.csv   panel: AF spectra BV510 coefficient per AF set
#                  <report_dir>/marker_scores.csv        TNF and Th2 readouts per donor and population
#                  <report_dir>/q2_flags.csv             grill Q2 flag per marker and population
#                  <report_dir>/gap_closed.pdf/.png/.csv panel: CD8T TNF gap closed per donor
#                  <report_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ---------------------------------------------------------------------------------------
STEP <- commandArgs(trailingOnly = TRUE)[1]
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
ARM <- "af405_afrefine"
PIPELINE <- "joint" # zip script
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under joint
REPORT_DIR <- "experiments/runs/2026-11-10_bmv-af405-forced-af"
SCORE_DIR <- "experiments/runs/2026-11-10_bmv-af405-forced-af-score"
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
CUT <- 1 # negatives below / positives above arcsinh 1 (user, 2026-11-10)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| < this fraction of as188_joint's rSD (39b; user, 2026-11-10)
Q2_DELTA <- 0.10 # grill Q2: donor-median rise in squeezed share that flags a marker (user, 2026-11-10)
Q2_RSD <- 0.8 # grill Q2: and donor-median negative rSD ratio at or below this
MARKERS <- c(TNF = "BV510-A", Th2 = "eFluor 450-A")
POPS <- c("CD4T", "CD8T")
JITTER_WIDTH <- 0.15 # display only
JITTER_SEED <- 1 # display only

stopifnot("step is build or report" = STEP %in% c("build", "report"))
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
out_dir <- file.path(DATA, "unmix_cells", ARM)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra <- read_raw(file.path(
  d188,
  "table_spectra",
  "Automated_autospectral_spectra_20261105_143803.csv"
))
af188 <- read_raw(file.path(
  d188,
  "table_spectra",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
))

# Build: refined AF set, then unmix with everything else from the zip ---------------------------------
if (STEP == "build") {
  stopifnot(
    "out dir must not exist" = !dir.exists(out_dir),
    "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) ==
      "1.8.8"
  )
  dir.create(out_dir, recursive = TRUE)
  control_dir <- file.path(d188, "Reference Group Clean")
  unstained <- list.files(
    file.path(d188, "Reference Group", "unstained_AF"),
    "\\.fcs$",
    full.names = TRUE
  )
  stopifnot(
    "one unstained AF file, as the zip script requires" = length(unstained) == 1
  )
  variants188 <- readRDS(file.path(
    d188,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  ))
  asp <- get.autospectral.param(cytometer = "aurora", figures = TRUE)

  # the zip call with refine = TRUE; output folders point into the arm dir, not the repo root
  af_refined <- get.af.spectra(
    unstained.sample = unstained,
    asp = asp,
    spectra = spectra,
    refine = TRUE,
    plot.dir = file.path(out_dir, "figure_autofluorescence"),
    table.dir = out_dir
  )
  write.csv(af_refined, file.path(out_dir, "af_refined.csv"))
  writeLines(
    c(
      paste("arm:", ARM),
      "af: get.af.spectra(unstained.sample = H12, asp, spectra = zip 143803, refine = TRUE); other args default",
      paste(
        "af spectra:",
        nrow(af_refined),
        "(zip refine = FALSE:",
        nrow(af188),
        ")"
      ),
      "spectra: zip 143803 (unchanged); variants: zip Spectral_variants.rds (unchanged)",
      paste("pipeline:", PIPELINE, "; speed:", SPEED),
      paste("controls:", control_dir),
      "seed: none set (as the zip script)"
    ),
    file.path(out_dir, "args.txt")
  )
  flow.control <- reload.flow.control(
    control_dir,
    file.path(control_dir, "fcs_control_file.csv"),
    asp
  )
  unmix.folder(
    fcs.dir = file.path(d188, "Raw"),
    spectra = spectra,
    asp = asp,
    flow.control = flow.control,
    method = "AutoSpectral",
    af.spectra = af_refined,
    spectra.variants = variants188,
    pipeline = PIPELINE,
    speed = SPEED,
    output.dir = out_dir,
    parallel = TRUE,
    verbose = TRUE
  )

  # Verify: same files, event counts and channel names as as188_joint ----------------------------------
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
    "5 files, same events and channels as as188_joint" = nrow(verify) == 5 &&
      all(verify$n_ref == verify$n_new) &&
      all(verify$channels_equal)
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )
}

# Report: AF set BV510 part, TNF and Th2 readouts with the Q2 flag, gap-closed panel -------------------
if (STEP == "report") {
  stopifnot("report dir must not exist" = !dir.exists(REPORT_DIR))
  dir.create(REPORT_DIR, recursive = TRUE)
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/plots.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/helper75_marker-scores.R")

  # BV510 coefficient of each AF spectrum (script 74 metric), zip AF set vs refined AF set
  b510 <- spectra["BV510", ]
  af_refined <- read_raw(file.path(out_dir, "af_refined.csv"))
  af_coef <- list(
    `zip AF` = af188,
    `refined AF` = af_refined
  ) |>
    map(\(m) {
      tibble(
        af_spectrum = rownames(m),
        bv510_coef = as.vector(m[, names(b510)] %*% b510) / sum(b510^2)
      )
    }) |>
    list_rbind(names_to = "af_set")
  write_csv(af_coef, file.path(REPORT_DIR, "af_bv510_coef.csv"))
  af_pal <- categorical_palette(unique(af_coef$af_set))
  p_af <- plot_box(
    af_coef,
    group = "af_set",
    value = "bv510_coef",
    pal = af_pal,
    jitter_width = JITTER_WIDTH,
    seed = JITTER_SEED,
    value_label = "BV510 coefficient of AF spectrum",
    title = "AF spectra: BV510 part"
  )
  save_panel(
    p_af,
    af_coef,
    REPORT_DIR,
    "af_bv510_coef",
    width_mm = 89,
    height_mm = 70
  )

  # TNF and Th2 readouts per donor on as188_joint's CD4T and CD8T rows (helper75, 39b metrics)
  pooled <- read_csv(POOLED_MANIFEST, show_col_types = FALSE)
  gs <- load_gs("experiments/runs/2026-11-10_bmv-ablation-gated/gs/as188_joint")
  donors <- sort(pooled$sample[pooled$arm == "as188_joint"])
  marker_scores <- donors |>
    set_names() |>
    map(\(s) {
      ref <- arm_matrix("as188_joint", s, unname(MARKERS), pooled)
      x <- arm_matrix(ARM, s, unname(MARKERS), pooled)
      map(set_names(POPS), \(p) {
        score_channels(x, ref, which(gh_pop_get_indices(gs[[s]], p)))
      }) |>
        list_rbind(names_to = "pop")
    }) |>
    list_rbind(names_to = "donor")
  ref_scores <- donors |>
    set_names() |>
    map(\(s) {
      ref <- arm_matrix("as188_joint", s, unname(MARKERS), pooled)
      map(set_names(POPS), \(p) {
        score_channels(ref, ref, which(gh_pop_get_indices(gs[[s]], p)))
      }) |>
        list_rbind(names_to = "pop")
    }) |>
    list_rbind(names_to = "donor") |>
    select(donor, pop, channel, ref_squeezed = squeezed_share)
  marker_scores <- marker_scores |>
    left_join(
      ref_scores,
      by = join_by(donor, pop, channel),
      relationship = "one-to-one",
      unmatched = "error"
    ) |>
    mutate(
      marker = names(MARKERS)[match(channel, MARKERS)],
      d_squeezed = squeezed_share - ref_squeezed
    )
  write_csv(marker_scores, file.path(REPORT_DIR, "marker_scores.csv"))

  # grill Q2 flag: median rise >= Q2_DELTA, median rSD ratio <= Q2_RSD, same sign of the rise in all donors
  q2_flags <- marker_scores |>
    summarise(
      median_d_squeezed = median(d_squeezed),
      median_rsd_ratio = median(rsd_ratio),
      same_sign = all(sign(d_squeezed) == sign(median(d_squeezed))),
      .by = c(marker, pop)
    ) |>
    mutate(
      flagged = median_d_squeezed >= Q2_DELTA &
        median_rsd_ratio <= Q2_RSD &
        same_sign
    )
  write_csv(q2_flags, file.path(REPORT_DIR, "q2_flags.csv"))

  # CD8T TNF gap closed per donor (51 scorer output), donors as points over the median
  gap <- read_csv(
    file.path(SCORE_DIR, "gap_closed.csv"),
    show_col_types = FALSE
  ) |>
    filter(arm == ARM) |>
    mutate(
      readout = if_else(
        readout == "squeezed_share",
        "TNF squeezed near 0",
        "TNF above arcsinh 1"
      )
    )
  gap_pal <- categorical_palette(c(
    "TNF squeezed near 0",
    "TNF above arcsinh 1"
  ))
  p_gap <- plot_box(
    gap,
    group = "readout",
    value = "frac_closed",
    pal = gap_pal,
    jitter_width = JITTER_WIDTH,
    seed = JITTER_SEED,
    value_label = "Fraction of 1.8.8 to 1.0.2 gap closed",
    title = "Refined AF: CD8T TNF gap closed"
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")
  save_panel(
    p_gap,
    gap,
    REPORT_DIR,
    "gap_closed",
    width_mm = 89,
    height_mm = 70
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(REPORT_DIR, "sessionInfo.txt")
  )
}
