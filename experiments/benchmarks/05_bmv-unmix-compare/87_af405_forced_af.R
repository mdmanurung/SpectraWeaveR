# Purpose: grill Q4/Q9 (README rows "Forced AF removal on AF405", user 2026-11-10): force the AF handling the
#          zip 1.8.8 run skipped and see whether it changes the AF405 (GLUT1) objects and the TNF / Th2 squeeze.
#            af405_afrefine  the zip 1.8.8 setup with get.af.spectra(refine = TRUE) instead of the zip's
#                            refine = FALSE; everything else as the zip script (AutoSpectral_unmix_v1.8.8.R):
#                            zip dye spectra (config stamp), zip variants, joint, Reference Group Clean controls.
#                            The zip built the AF405 spectrum (get.spectra.automated) and variants
#                            (get.spectral.variants takes no AF input in 1.8.8) without the AF set, so only the
#                            AF set used in the unmix changes; the AF405 objects stay the zip's.
#            af405_afremove  NOT BUILT: clean.controls() needs a define.flow.control() gating object, and the
#                            zip's automated route never made one; how the cleaned AF405 events reach
#                            get.spectra.automated / get.fluorophore.spectra and get.spectral.variants (which
#                            reads FCS files, not a flow.control) is not set by the zip. User choice pending.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD>, then
#          step build (SLURM: experiments/slurm/run_script.sbatch):
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/87_af405_forced_af.R build [condition]
#              -> <data_dir>/unmix_cells/af405_afrefine[/<condition>] (must not exist)
#          then score with 51: sbatch experiments/slurm/run_score_ablation.sbatch \
#              <run folder af405-forced-af-score> af405_afrefine
#          step report (after the scorer; register the score run folder first): ... 87_af405_forced_af.R report [condition]
#              -> run_out("af405-forced-af") = experiments/runs/<SW_RUN_DATE>_bmv<run_tag>-af405-forced-af (must not exist)
#          condition = last argument, required when the dataset has conditions (24H: one task per condition; report
#          reads that condition's AF set), refused otherwise.
# Inputs:  <data_dir>/autoSpectral_v1.8.8/  config-stamp zip spectra and AF, variants,
#            Reference Group Clean, Reference Group/unstained_AF/*.fcs, Raw[/<condition>], AutoSpectral_unmixed[/<condition>] [user zip]
#          run key af405-forced-af-score: gap_closed.csv   [51_score_ablation.R]
#          run key ablation-gated: gs/as188_joint; run key gated-pooled: manifest.csv (via helper75)
#            [19a, 12_cleanup.R]; helper75_marker-scores.R; lib/theme.R, lib/palettes.R, lib/plots.R
#          stamp file names (spectra, AF) from config stamps; samples from datasets/<ds>_samples.csv [00_samples.R];
#            config channels (TNF_slot, Th2) and targets (tnf_dye, tnf_label); the report stops when the config has no Th2 channel
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

# Dataset config -----------------------------------------------------------------------------------
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ---------------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
STEP <- cli_args[1]
DATA <- normalizePath(data_path())
ARM <- "af405_afrefine"
PIPELINE <- "joint" # zip script
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under joint
TNF <- cfg("targets")$tnf_label # marker label of the TNF slot (3D: TNF)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
CUT <- 1 # negatives below / positives above arcsinh 1 (user, 2026-11-10; 3D only)
if (DS$dataset != "3d") {
  stop(
    "87_af405_forced_af.R: CUT = 1 is the 3D value; the cut rule (cut_of()) is the user's decision for dataset ",
    DS$dataset,
    call. = FALSE
  )
}
SQUEEZE_FRAC <- 0.5 # squeezed = |x| < this fraction of as188_joint's rSD (39b; user, 2026-11-10)
Q2_DELTA <- 0.10 # grill Q2: donor-median rise in squeezed share that flags a marker (user, 2026-11-10)
Q2_RSD <- 0.8 # grill Q2: and donor-median negative rSD ratio at or below this
TNF_CHANNEL <- cfg("channels")$TNF_slot # 3D: BV510-A
TH2_CHANNEL <- cfg("channels")$Th2 # 3D: eFluor 450-A; NULL when the dataset has no Th2 channel
POPS <- c("CD4T", "CD8T")
JITTER_WIDTH <- 0.15 # display only
JITTER_SEED <- 1 # display only

stopifnot("step is build or report" = STEP %in% c("build", "report"))
# 24H: one task per condition subfolder (argument 2); 3D and EXV have no condition level
cond <- if (length(cli_args) >= 2) cli_args[[2]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "87_af405_forced_af.R: dataset ",
    DS$dataset,
    " has no conditions; drop argument 2 ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "87_af405_forced_af.R: needs argument 2, a condition of dataset ",
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
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
out_dir <- paste(c(file.path(DATA, "unmix_cells", ARM), cond), collapse = "/")
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra <- read_raw(file.path(
  d188,
  "table_spectra",
  cfg("stamps")$as188_spectra
))
af188 <- read_raw(file.path(d188, "table_spectra", cfg("stamps")$as188_af))

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
      "af: get.af.spectra(unstained.sample = the unstained AF file, asp, spectra = zip config stamp, refine = TRUE); other args default",
      paste(
        "af spectra:",
        nrow(af_refined),
        "(zip refine = FALSE:",
        nrow(af188),
        ")"
      ),
      "spectra: zip config stamp (unchanged); variants: zip Spectral_variants.rds (unchanged)",
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
    fcs.dir = paste(c(file.path(d188, "Raw"), cond), collapse = "/"),
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
    "one file per sample, same events and channels as as188_joint" = nrow(
      verify
    ) ==
      nrow(samples) &&
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
  if (is.null(TH2_CHANNEL)) {
    stop(
      "87_af405_forced_af.R: needs the Th2 channel (config channels$Th2), missing in dataset ",
      DS$dataset,
      call. = FALSE
    )
  }
  MARKERS <- set_names(c(TNF_CHANNEL, TH2_CHANNEL), c(TNF, "Th2")) # marker label = channel
  REPORT_DIR <- run_out("af405-forced-af")
  SCORE_DIR <- run_in("af405-forced-af-score")
  stopifnot("report dir must not exist" = !dir.exists(REPORT_DIR))
  dir.create(REPORT_DIR, recursive = TRUE)
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/plots.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/helper75_marker-scores.R")

  # BV510 coefficient of each AF spectrum (script 74 metric), zip AF set vs refined AF set
  b510 <- spectra[cfg("targets")$tnf_dye, ] # TNF dye signature (3D: BV510)
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
    value_label = paste(cfg("targets")$tnf_dye, "coefficient of AF spectrum"),
    title = paste("AF spectra:", cfg("targets")$tnf_dye, "part")
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
  gs <- load_gs(file.path(run_in("ablation-gated"), "gs/as188_joint"))
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
  # 24H: one summary per condition (marker_scores gains a condition column); 3D and EXV keep the old columns
  q2_by <- c("marker", "pop")
  if (!is.null(DS$conditions)) {
    marker_scores <- marker_scores |>
      left_join(
        select(samples_tbl(), donor = sample, condition),
        by = join_by(donor),
        relationship = "many-to-one",
        unmatched = "error"
      )
    q2_by <- c(q2_by, "condition")
  }
  write_csv(marker_scores, file.path(REPORT_DIR, "marker_scores.csv"))

  # grill Q2 flag: median rise >= Q2_DELTA, median rSD ratio <= Q2_RSD, same sign of the rise in all donors
  q2_flags <- marker_scores |>
    summarise(
      median_d_squeezed = median(d_squeezed),
      median_rsd_ratio = median(rsd_ratio),
      same_sign = all(sign(d_squeezed) == sign(median(d_squeezed))),
      .by = all_of(q2_by)
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
        paste(TNF, "squeezed near 0"),
        paste(TNF, "above arcsinh", CUT)
      )
    )
  gap_pal <- categorical_palette(c(
    paste(TNF, "squeezed near 0"),
    paste(TNF, "above arcsinh", CUT)
  ))
  p_gap <- plot_box(
    gap,
    group = "readout",
    value = "frac_closed",
    pal = gap_pal,
    jitter_width = JITTER_WIDTH,
    seed = JITTER_SEED,
    value_label = "Fraction of 1.8.8 to 1.0.2 gap closed",
    title = paste("Refined AF: CD8T", TNF, "gap closed")
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
