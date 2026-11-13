# Purpose: grill Q4/Q9, arm af405_afremove (README row "af405_afremove arm (script 87b)", user 2026-11-10).
#          Run AutoSpectral's cell AF removal on the AF405 (GLUT1) control, rebuild only AF405's spectrum row and
#          variants from the cleaned positives, swap them into the zip 1.8.8 setup, unmix and score.
#            1. landmark gate: define.gate.landmarks (config control_gates$landmark: 3D n.cells 2000, percentile 70)
#               on Reference Group Clean with the landmark control file, then define.flow.control(gating.system = "landmarks")
#            2. clean.controls(af.remove = TRUE, downsample = FALSE, positive.n = 1e6); other args default
#            3. cleaned AF405 positives (rows not copied from the universal negative) are taken from the
#               original AF405 FCS by event key and written as the AF405 control in a copy of Reference Group Clean
#            4. get.spectra.automated with the zip arguments and 1.8.8 get.spectral.variants on that copy
#            5. zip setup with AF405's spectrum row and its variant slots (variants, delta.list, delta.norms, as
#               45_variant_per_dye.R) replaced; zip AF set; joint, slow
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD>, then
#          build:  sbatch experiments/slurm/run_script_lib.sbatch \
#                    experiments/benchmarks/05_bmv-unmix-compare/87b_af405_afremove.R build [condition]
#                  -> <data_dir>/unmix_cells/af405_afremove[/<condition>] (must not exist)
#          score:  sbatch experiments/slurm/run_score_ablation.sbatch \
#                    <run folder af405-afremove-score> af405_afremove
#          report: sbatch experiments/slurm/run_script_lib.sbatch <this script> report [condition]
#                  -> run_out("af405-afremove") = experiments/runs/<SW_RUN_DATE>_bmv<run_tag>-af405-afremove (must not exist);
#                  register the score run folder first
#          condition = last argument, required when the dataset has conditions (24H: one task per condition; the AF405
#          control is rebuilt in each condition folder; report reads that condition's objects), refused otherwise.
# Inputs:  <data_dir>/autoSpectral_v1.8.8/  Reference Group Clean, config-stamp zip spectra and AF,
#            Spectral_variants.rds, Raw[/<condition>], AutoSpectral_unmixed[/<condition>] [user zip]
#          run key cell-gating: fcs_control_file_landmark.csv  landmark control file [50_cell_gate_plots.R, user-chosen gates]
#          run key af405-afremove-score: gap_closed.csv  [51_score_ablation.R]
#          run key ablation-gated: gs/as188_joint; run key gated-pooled: manifest.csv (via helper75)
#            [19a, 12_cleanup.R]; helper75_marker-scores.R; lib/theme.R, lib/palettes.R, lib/plots.R
#          stamp file names (spectra, AF) from config stamps; samples from datasets/<ds>_samples.csv [00_samples.R];
#            config control_gates$landmark, channels (TNF_slot, Th2) and targets (glut1_dye, tnf_dye, tnf_label);
#            the report stops when the config has no Th2 channel
# Outputs: build:  <out_dir>/controls/            Reference Group Clean copy with the cleaned AF405 control
#                  <out_dir>/objects/             landmark_gate.rds, spectra_auto.csv, variants_auto.rds
#                  <out_dir>/clean_log.csv        AF405 events: gated, cleaned rows, positives kept
#                  <out_dir>/*.fcs, spectra_used.csv, variants_used.rds, args.txt, verify.csv, sessionInfo.txt
#          report: <report_dir>/offset_spread.csv            AF405 BV510 offset and variant spread per set
#                  <report_dir>/variant_spread.pdf/.png/.csv  panel: BV510 coefficient per AF405 variant
#                  <report_dir>/marker_scores.csv, q2_flags.csv   TNF and Th2 readouts; revised Q2 flag
#                  <report_dir>/gap_closed.pdf/.png/.csv      panel: CD8T TNF gap closed per donor
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
DATA <- normalizePath(data_path()) # absolute: the build step changes its working directory
ARM <- "af405_afremove"
DYE <- cfg("targets")$glut1_dye # the GLUT1 dye (3D: Alexa Fluor 405)
TNF <- cfg("targets")$tnf_label # marker label of the TNF slot (3D: TNF)
LANDMARK_FILE <- if (STEP == "build") {
  normalizePath(
    file.path(run_in("cell-gating"), "fcs_control_file_landmark.csv")
  ) # written by 50_cell_gate_plots.R
}
LANDMARK_N_CELLS <- if (STEP == "build") cfg("control_gates")$landmark$n_cells # define.gate.landmarks, as chain7 (3D: 2000, user 2026-11-10)
LANDMARK_PERCENTILE <- if (STEP == "build") {
  cfg("control_gates")$landmark$percentile
} # same (3D: 70)
POSITIVE_N <- 1e6 # above the AF405 positive count, so every cleaned positive is kept (user, 2026-11-10)
SEED <- 42 # repo seed (README: bird.seed + 42)
VARIANT_SLOTS <- c("variants", "delta.list", "delta.norms") # per-dye slots swapped together, as in 45
PIPELINE <- "joint" # zip script
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under joint
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
CUT <- 1 # negatives below / positives above arcsinh 1 (user, 2026-11-10; 3D only)
if (DS$dataset != "3d") {
  stop(
    "87b_af405_afremove.R: CUT = 1 is the 3D value; the cut rule (cut_of()) is the user's decision for dataset ",
    DS$dataset,
    call. = FALSE
  )
}
SQUEEZE_FRAC <- 0.5 # squeezed = |x| < this fraction of as188_joint's rSD (39b; user, 2026-11-10)
Q2_DELTA <- 0.10 # revised Q2 flag: donor-median rise in squeezed share, same sign in all donors (user)
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
    "87b_af405_afremove.R: dataset ",
    DS$dataset,
    " has no conditions; drop argument 2 ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "87b_af405_afremove.R: needs argument 2, a condition of dataset ",
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
clean_dir <- file.path(d188, "Reference Group Clean")
out_dir <- paste(c(file.path(DATA, "unmix_cells", ARM), cond), collapse = "/")
obj_dir <- file.path(out_dir, "objects")
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra_zip <- read_raw(file.path(
  d188,
  "table_spectra",
  cfg("stamps")$as188_spectra
))
af_zip <- read_raw(file.path(d188, "table_spectra", cfg("stamps")$as188_af))
variants_zip <- readRDS(file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
))

# Build ----------------------------------------------------------------------------------------------
if (STEP == "build") {
  stopifnot(
    "out dir must not exist" = !dir.exists(out_dir),
    "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) ==
      "1.8.8"
  )
  dir.create(obj_dir, recursive = TRUE)
  setwd(obj_dir) # AutoSpectral writes tables and figures relative to the working directory
  set.seed(SEED)
  asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)

  # 1. landmark gate on Reference Group Clean, then landmark gating for the flow control
  landmark_gate <- define.gate.landmarks(
    control.file = LANDMARK_FILE,
    control.dir = clean_dir,
    asp = asp,
    gate.name = "lymphocyte",
    n.cells = LANDMARK_N_CELLS,
    percentile = LANDMARK_PERCENTILE,
    output.dir = file.path(obj_dir, "figure_gate")
  )
  saveRDS(landmark_gate, file.path(obj_dir, "landmark_gate.rds"))
  flow_control <- define.flow.control(
    clean_dir,
    LANDMARK_FILE,
    asp,
    gating.system = "landmarks",
    gate.list = list(lymphocyte = landmark_gate)
  )

  # 2. AF removal and universal-negative matching, every cleaned positive kept
  cleaned <- clean.controls(
    flow_control,
    asp,
    af.remove = TRUE,
    downsample = FALSE,
    positive.n = POSITIVE_N,
    main.figures = FALSE
  )

  # 3. cleaned AF405 positives: drop rows copied from the universal negative, then take the original
  #    events by key so the written control keeps every FCS column
  dye_sample <- cleaned$sample[cleaned$fluorophore == DYE]
  neg_sample <- cleaned$sample[cleaned$fluorophore == "AF"]
  dye_rows <- cleaned$clean.expr[
    cleaned$clean.event.sample == dye_sample,
    ,
    drop = FALSE
  ]
  key_cols <- intersect(
    c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "V2-A", "V7-A"),
    colnames(dye_rows)
  )
  key_of <- \(m) {
    do.call(paste, c(as.data.frame(m[, key_cols, drop = FALSE]), sep = "|"))
  }
  neg_raw <- exprs(read.FCS(
    file.path(clean_dir, "Unstained (Cells).fcs"),
    transformation = FALSE,
    truncate_max_range = FALSE
  ))
  pos_rows <- dye_rows[!key_of(dye_rows) %in% key_of(neg_raw), , drop = FALSE]
  dye_file <- paste0(DYE, " (Cells).fcs") # control file name of the GLUT1 dye
  dye_ff <- read.FCS(
    file.path(clean_dir, dye_file),
    transformation = FALSE,
    truncate_max_range = FALSE
  )
  keep <- key_of(exprs(dye_ff)) %in% key_of(pos_rows)
  stopifnot(
    "every cleaned positive found in the original AF405 file" = sum(keep) >=
      nrow(pos_rows)
  )
  ctrl_dir <- file.path(out_dir, "controls")
  dir.create(ctrl_dir)
  file.copy(list.files(clean_dir, full.names = TRUE), ctrl_dir)
  write.FCS(dye_ff[keep, ], file.path(ctrl_dir, dye_file))
  write_csv(
    tibble(
      dye_sample,
      neg_sample,
      events_in_file = nrow(dye_ff),
      cleaned_rows = nrow(dye_rows),
      negative_rows_dropped = nrow(dye_rows) - nrow(pos_rows),
      positives_kept = sum(keep)
    ),
    file.path(out_dir, "clean_log.csv")
  )

  # 4. automated spectra (zip arguments) and 1.8.8 variants on the copy
  ctrl_csv <- file.path(ctrl_dir, "fcs_control_file.csv")
  spectra_auto <- get.spectra.automated(
    control.dir = ctrl_dir,
    control.def.file = ctrl_csv,
    asp = asp,
    n.candidates = 10000,
    n.spectral = 500,
    cosine.threshold = 0.95
  )
  write.csv(spectra_auto, file.path(obj_dir, "spectra_auto.csv"))
  variants_auto <- get.spectral.variants(
    control.dir = ctrl_dir,
    control.def.file = ctrl_csv,
    asp = asp,
    spectra = spectra_auto
  )
  saveRDS(variants_auto, file.path(obj_dir, "variants_auto.rds"))

  # 5. zip setup with only AF405's spectrum row and variant slots replaced
  spectra <- spectra_zip
  spectra[DYE, ] <- spectra_auto[DYE, colnames(spectra)]
  spectra_variants <- variants_zip
  for (slot in VARIANT_SLOTS) {
    spectra_variants[[slot]][[DYE]] <- variants_auto[[slot]][[DYE]]
  }
  write.csv(spectra, file.path(out_dir, "spectra_used.csv"))
  saveRDS(spectra_variants, file.path(out_dir, "variants_used.rds"))
  writeLines(
    c(
      paste("arm:", ARM),
      paste(
        "landmark gate: define.gate.landmarks(n.cells =",
        LANDMARK_N_CELLS,
        ", percentile =",
        LANDMARK_PERCENTILE,
        ") on Reference Group Clean"
      ),
      paste(
        "clean.controls(af.remove = TRUE, downsample = FALSE, positive.n =",
        POSITIVE_N,
        "); other args default"
      ),
      "spectra: zip config stamp with the AF405 row from get.spectra.automated(n.candidates 10000, n.spectral 500, cosine 0.95) on the cleaned copy",
      paste(
        "variants: zip, with AF405",
        paste(VARIANT_SLOTS, collapse = "/"),
        "from get.spectral.variants (defaults) on the cleaned copy"
      ),
      paste("AF405 variants:", nrow(variants_auto$variants[[DYE]])),
      "af: zip config stamp",
      paste("pipeline:", PIPELINE, "; speed:", SPEED, "; seed:", SEED)
    ),
    file.path(out_dir, "args.txt")
  )
  setwd(Sys.getenv("SW_REPO"))
  flow_control_unmix <- reload.flow.control(
    clean_dir,
    file.path(clean_dir, "fcs_control_file.csv"),
    asp
  )
  unmix.folder(
    fcs.dir = paste(c(file.path(d188, "Raw"), cond), collapse = "/"),
    spectra = spectra,
    asp = asp,
    flow.control = flow_control_unmix,
    method = "AutoSpectral",
    af.spectra = af_zip,
    spectra.variants = spectra_variants,
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

# Report: offset and spread, TNF and Th2 with the revised Q2 flag, gap closed ---------------------------
if (STEP == "report") {
  if (is.null(TH2_CHANNEL)) {
    stop(
      "87b_af405_afremove.R: needs the Th2 channel (config channels$Th2), missing in dataset ",
      DS$dataset,
      call. = FALSE
    )
  }
  MARKERS <- set_names(c(TNF_CHANNEL, TH2_CHANNEL), c(TNF, "Th2")) # marker label = channel
  REPORT_DIR <- run_out("af405-afremove")
  SCORE_DIR <- run_in("af405-afremove-score")
  stopifnot("report dir must not exist" = !dir.exists(REPORT_DIR))
  dir.create(REPORT_DIR, recursive = TRUE)
  lib <- Sys.getenv("SW_LIB", "experiments/benchmarks/05_bmv-unmix-compare/lib")
  source(file.path(lib, "theme.R"))
  source(file.path(lib, "palettes.R"))
  source(file.path(lib, "plots.R"))
  source("experiments/benchmarks/05_bmv-unmix-compare/helper75_marker-scores.R")

  # BV510 offset (reference - zip 1.8.8 reference) and variant spread (variant - own reference), script 74 metric
  b510 <- spectra_zip[cfg("targets")$tnf_dye, ] # TNF dye signature (3D: BV510)
  proj_coef <- \(d) as.vector(d %*% b510) / sum(b510^2)
  spectra_used <- read_raw(file.path(out_dir, "spectra_used.csv"))
  variants_used <- readRDS(file.path(out_dir, "variants_used.rds"))
  spread <- list(
    `zip 1.8.8` = list(
      ref = spectra_zip[DYE, ],
      v = variants_zip$variants[[DYE]]
    ),
    `AF removed` = list(
      ref = spectra_used[DYE, ],
      v = variants_used$variants[[DYE]]
    )
  ) |>
    map(\(s) {
      tibble(
        spread_coef = proj_coef(sweep(
          s$v[, names(b510)],
          2,
          s$ref[names(b510)]
        )),
        offset_coef = sum(
          (s$ref[names(b510)] - spectra_zip[DYE, names(b510)]) * b510
        ) /
          sum(b510^2)
      )
    }) |>
    list_rbind(names_to = "set")
  write_csv(
    summarise(
      spread,
      n_variants = n(),
      median_spread = median(spread_coef),
      p90_spread = quantile(spread_coef, 0.9),
      offset = first(offset_coef),
      .by = set
    ),
    file.path(REPORT_DIR, "offset_spread.csv")
  )
  set_pal <- categorical_palette(unique(spread$set))
  p_spread <- plot_box(
    spread,
    group = "set",
    value = "spread_coef",
    pal = set_pal,
    jitter_width = JITTER_WIDTH,
    seed = JITTER_SEED,
    value_label = paste(
      cfg("targets")$tnf_dye,
      "coefficient of variant spread"
    ),
    title = paste("AF405 variants:", cfg("targets")$tnf_dye, "part")
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")
  save_panel(
    p_spread,
    spread,
    REPORT_DIR,
    "variant_spread",
    width_mm = 89,
    height_mm = 70
  )

  # TNF and Th2 readouts per donor on as188_joint's CD4T and CD8T rows (helper75, 39b metrics)
  pooled <- read_csv(POOLED_MANIFEST, show_col_types = FALSE)
  gs <- load_gs(file.path(run_in("ablation-gated"), "gs/as188_joint"))
  donors <- sort(pooled$sample[pooled$arm == "as188_joint"])
  score_arm <- \(arm_name) {
    donors |>
      set_names() |>
      map(\(s) {
        ref <- arm_matrix("as188_joint", s, unname(MARKERS), pooled)
        x <- arm_matrix(arm_name, s, unname(MARKERS), pooled)
        map(set_names(POPS), \(p) {
          score_channels(x, ref, which(gh_pop_get_indices(gs[[s]], p)))
        }) |>
          list_rbind(names_to = "pop")
      }) |>
      list_rbind(names_to = "donor")
  }
  ref_scores <- score_arm("as188_joint") |>
    select(donor, pop, channel, ref_squeezed = squeezed_share)
  marker_scores <- score_arm(ARM) |>
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

  # revised Q2 flag: donor-median rise >= Q2_DELTA and the same sign in all donors (no rSD criterion)
  q2_flags <- marker_scores |>
    summarise(
      median_d_squeezed = median(d_squeezed),
      median_rsd_ratio = median(rsd_ratio),
      same_sign = all(sign(d_squeezed) == sign(median(d_squeezed))),
      .by = all_of(q2_by)
    ) |>
    mutate(flagged = median_d_squeezed >= Q2_DELTA & same_sign)
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
    title = paste("AF removed: CD8T", TNF, "gap closed")
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
