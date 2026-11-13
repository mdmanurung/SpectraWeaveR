# Copy of 88d_cospike_tnf_glut1.R for the common SpectroFlo labels (plans/common_gates_plan.md G5, G10b). Only the inputs (the spike keys and scores of the 88cg family)
#          and the output dirs differ; method, levels, seed 42, loss cut and setups are inherited unchanged.
# Purpose: TNF + GLUT1 co-spike (grill 2026-11-11 Q6, Q7; README row "Co-spike (grill Q6, Q7; script
#          88d)"). The TNF +1 arcsinh cells of the first spike-in (88 spike_key_<tnf_label>.csv: CD4T and CD8T
#          cells, seed 42) get, in the same cells, TNF (BV510) +1 arcsinh plus GLUT1 (Alexa Fluor 405)
#          at +0, +1 or +2 arcsinh, each dye as its 1.8.8 spectrum added to the raw detector values. One
#          spiked set per GLUT1 level; the unspiked quarter (level 0 in the key) stays unspiked; the other
#          cells keep their raw values. Setups c188 (as188_joint), afrefine (af405_afrefine) and s102_c102
#          (as102_var, positive control) unmix each set; recovered = spiked unmix - unspiked unmix of the
#          same cell, with the unspiked unmixes spk_<setup>_base of 88 reused. Readouts per cell: TNF
#          recovery = recovered TNF / TNF added; GLUT1 recovery = recovered GLUT1 / GLUT1 added (GLUT1 > 0).
#          Unmix mode is the 88b unmix block with the set names changed.
# Run:     from the repo root, on SLURM (experiments/slurm/run_spike_cg.sbatch 88dcg_cospike_tnf_glut1.R <mode>), with
#          SW_DATASET (3d | 24h | exv) and SW_RUN_DATE (YYYY-MM-DD) set; the dataset comes from datasets/<SW_DATASET>.json:
#            build                 spike the raw files, one set per GLUT1 level
#            unmix <setup> <set> [condition]   setup = c188 | afrefine | s102_c102; set = cospike_g0 | cospike_g1 | cospike_g2;
#                                  condition = 4th argument, required when the dataset has conditions (24H), refused otherwise
#            score                 recovery tables
#            figures               summary and ggcyto raw figures from the saved tables
#          SW_SMOKE_DONOR: any value marks a smoke run; the cells come from the smoke build of 88cg (run key
#          spike-in-cg-build-smoke), so the donor choice is made there.
# Inputs:  <data_dir>/autoSpectral_v1.8.8/{Raw, table_spectra, Reference Group Clean, figure_spectral_variants}  [user zip]
#          <data_dir>/autoSpectral_v1.0.2/Reference Group/  [user zip]
#          <data_dir>/lib_as102/, <data_dir>/lib_ggplot2_35/  code libraries (config lib_dirs) [06_unmix_c102.R setup]
#          stamp file names (spectra, AF) from config stamps; dyes and marker labels from config targets;
#            samples from datasets/<ds>_samples.csv [00_samples.R]
#          run key spike-in-cg-build: spike_key_<targets$tnf_label>.csv  which cells, and their level [88cg_spike_in_markers.R build]
#          run key opt188-smoke: opt188 objects, setup opt188 only [53b_opt188_unmix.R]
#          run key cospike-tnf-glut1-cg: build run, read by score and figures [this script, mode build]
#          <data_dir>/unmix_cells/spk_<setup>_base/[<condition>/]  unspiked unmixes [88bcg_spike_in_panel2.R unmix <setup> base]
#          <data_dir>/unmix_cells/af405_afrefine/af_refined.csv  refined AF set [87]
# Outputs: <data_dir>/spiked/<SPIKED_PREFIX><set>/[<condition>/]*.fcs   spiked raw files
#          <data_dir>/unmix_cells/<UNMIX_PREFIX><setup>_<set>/[<condition>/]*.fcs  unmixed files
#          <RUN_DIR> = run_out("cospike-tnf-glut1-cg") (build), the same registered folder in score and figures:
#          <RUN_DIR>/cospike_key.csv          sample, raw_row, pop, raw_file, level, TNF and GLUT1 amounts
#          <RUN_DIR>/cospike_cells.csv        per co-spiked cell, setup and set: TNF and GLUT1 before, recovered,
#                                             recovery (24H adds condition)
#          <RUN_DIR>/cospike_summary.csv      per setup, set, pop, sample: % lost >50% TNF, median recoveries
#          <RUN_DIR>/baseline_check.csv       unspiked cells: max |spiked - base| per channel
#          <RUN_DIR>/fig_cospike_tnf_loss.{pdf,png,csv}       % TNF lost >50% by GLUT1 level, donors as points
#          <RUN_DIR>/fig_cospike_glut1_recovery.{pdf,png,csv} median GLUT1 recovery by level, donors as points
#          <RUN_DIR>/fig_cospike_raw_<pop>[_<condition>].{pdf,png}      ggcyto hexbin TNF before vs after, levels x setups
#          <RUN_DIR>/fig_cospike_raw_<pop>[_<condition>]_v2.{pdf,png}   the same, two-line title (mode figures_raw)
#          <RUN_DIR>/sessionInfo_<mode>.txt (build, score, figures)
#          Register the build run folder in datasets/<ds>_runs.csv before score reads it.

# Dataset config -----------------------------------------------------------------------------------
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ---------------------------------------------------------------------------------------
DATA <- data_path()
SMOKE_DONOR <- Sys.getenv("SW_SMOKE_DONOR", "") # any value = smoke run (cells come from the 88cg smoke build); empty = full run
SMOKE_TAG <- if (nzchar(SMOKE_DONOR)) "smoke_" else "" # smoke runs write to their own dirs
SMOKE_SUFFIX <- if (nzchar(SMOKE_DONOR)) "-smoke" else ""
TNF_LABEL <- cfg("targets")$tnf_label # marker label of the TNF-slot dye in titles, axes and key file names
TNF_KEY_RUN <- paste0("spike-in-cg-build", SMOKE_SUFFIX) # run key of the 88cg build: the cells, as 88cg chose them
RUN_KEY <- paste0("cospike-tnf-glut1-cg", SMOKE_SUFFIX) # run key of this script's run folder
SPIKED_PREFIX <- paste0(SMOKE_TAG, "2026-11-11_cg_") # spiked raw sets, apart from the old 2026-11-11_cospike_g*
UNMIX_PREFIX <- paste0(SMOKE_TAG, "spkcg_") # unmixed spiked sets; the base unmixes keep spk_<setup>_base
BASE_PREFIX <- "spk_" # unspiked (base) unmixes read by score; written by 88bcg_spike_in_panel2.R unmix <setup> base, no smoke tag
# opt188 objects and the refined AF set, as 88b (the unmix block reads both; opt188 is not a setup here)
OPT188_ARM <- "opt188_pregate_neg_same_fsc4e5_p1" # folder of the opt188 objects inside run key opt188-smoke (resolved in the unmix block)
AFREFINE_AF <- data_path("unmix_cells/af405_afrefine/af_refined.csv")
TNF_DYE <- cfg("targets")$tnf_dye
GLUT1_DYE <- cfg("targets")$glut1_dye
TNF_LEVEL <- 1 # TNF spike, arcsinh units (user, 2026-11-11, grill Q6)
SETS <- c(cospike_g0 = 0, cospike_g1 = 1, cospike_g2 = 2) # GLUT1 co-spike, arcsinh units (grill Q6)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
LOSS_FRAC <- 0.5 # "losing most of the spike" = recovery below 0.5 (64 readout)
DONOR_SHAPES <- c(16, 17, 15, 18, 8, 3, 4, 7, 9, 10) # 24H figure point shapes in sorted donor order (3D keeps the ggplot default)
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across unmixes of one raw file
SETUPS <- c("c188", "afrefine", "s102_c102") # grill Q7: s102_c102 is the positive control
SETUPS_188 <- c("c188", "opt188", "afrefine") # run on the AutoSpectral 1.8.8 library
SETUP_LABELS <- c(
  c188 = "1.8.8 zip",
  afrefine = "refined AF",
  s102_c102 = "1.0.2 zip"
)

if (nzchar(SMOKE_DONOR)) {
  SETS <- SETS["cospike_g1"] # smoke: one GLUT1 level only
}

cli_args <- commandArgs(trailingOnly = TRUE)
mode <- cli_args[[1]]
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
d102 <- file.path(DATA, "autoSpectral_v1.0.2")
spectra_188_file <- file.path(
  d188,
  "table_spectra",
  cfg("stamps")$as188_spectra
)
read_matrix <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}

# Mode build: one spiked raw set per GLUT1 level ------------------------------------------------------
if (mode == "build") {
  TNF_KEY_DIR <- run_in(TNF_KEY_RUN) # resolved before library(flowCore): its filter() masks dplyr's inside run_in()
  suppressPackageStartupMessages({
    library(flowCore)
    library(dplyr)
    library(purrr)
    library(readr)
  })
  RUN_DIR <- run_out(RUN_KEY) # stops when the folder exists
  dir.create(RUN_DIR, recursive = TRUE)
  s188 <- read_matrix(spectra_188_file)
  # level-1 cells of the TNF key get the co-spike; level-0 cells stay unspiked as the check
  key <- read_csv(
    file.path(TNF_KEY_DIR, paste0("spike_key_", TNF_LABEL, ".csv")),
    show_col_types = FALSE
  ) |>
    dplyr::filter(level %in% c(0, TNF_LEVEL)) |>
    mutate(
      spiked = level == TNF_LEVEL,
      k_tnf = if_else(spiked, COFACTOR * sinh(TNF_LEVEL), 0)
    )
  write_csv(key, file.path(RUN_DIR, "cospike_key.csv"))
  iwalk(SETS, \(glut1_level, set) {
    out_fcs <- file.path(DATA, "spiked", paste0(SPIKED_PREFIX, set))
    stopifnot("spiked fcs dir must not exist" = !dir.exists(out_fcs))
    dir.create(out_fcs, recursive = TRUE)
    k_glut1 <- COFACTOR * sinh(glut1_level)
    split(key, key$raw_file) |>
      walk(\(key_df) {
        ff <- read.FCS(
          file.path(d188, "Raw", key_df$raw_file[1]),
          transformation = FALSE,
          truncate_max_range = FALSE
        )
        rows <- key_df$raw_row[key_df$spiked]
        x <- exprs(ff)
        x[rows, colnames(s188)] <- x[rows, colnames(s188)] +
          outer(
            rep(COFACTOR * sinh(TNF_LEVEL), length(rows)),
            s188[TNF_DYE, ]
          ) +
          outer(rep(k_glut1, length(rows)), s188[GLUT1_DYE, ])
        exprs(ff) <- x
        out_file <- file.path(out_fcs, key_df$raw_file[1])
        dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE) # condition subfolder in 24H
        write.FCS(ff, out_file)
      })
  })
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_build.txt")
  )
}

# Mode unmix: one setup on one spiked set (88b unmix block) -------------------------------------------
if (mode == "unmix") {
  setup <- cli_args[[2]]
  set <- cli_args[[3]]
  stopifnot(
    "known setup and set" = setup %in%
      SETUPS &&
      set %in% names(SETS)
  )
  # 24H: one task per condition subfolder (4th job argument); 3D and EXV have no condition level
  cond <- if (length(cli_args) >= 4) cli_args[[4]] else NULL
  if (is.null(DS$conditions) && !is.null(cond)) {
    stop(
      "88dcg_cospike_tnf_glut1.R: dataset ",
      DS$dataset,
      " has no conditions; drop the 4th argument ",
      cond,
      call. = FALSE
    )
  }
  if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
    stop(
      "88dcg_cospike_tnf_glut1.R: unmix needs a 4th argument, a condition of dataset ",
      DS$dataset,
      ": ",
      paste(DS$conditions, collapse = ", "),
      call. = FALSE
    )
  }
  raw_dir <- file.path(DATA, "spiked", paste0(SPIKED_PREFIX, set))
  raw_dir <- paste(c(raw_dir, cond), collapse = "/")
  out_dir <- normalizePath(file.path(DATA, "unmix_cells"), mustWork = TRUE) |>
    file.path(paste0(UNMIX_PREFIX, setup, "_", set))
  out_dir <- paste(c(out_dir, cond), collapse = "/")
  stopifnot("out dir must not exist" = !dir.exists(out_dir))
  if (!setup %in% SETUPS_188) {
    .libPaths(c(
      cfg("lib_dirs")[basename(cfg("lib_dirs")) == "lib_as102"],
      .libPaths()
    ))
  }
  OPT188_RUN <- if (setup == "opt188") run_in("opt188-smoke") # before library(flowCore), as in build
  suppressPackageStartupMessages({
    library(AutoSpectral)
    library(flowCore)
  })
  stopifnot(
    "AutoSpectral version matches setup" = as.character(packageVersion(
      "AutoSpectral"
    )) ==
      if (setup %in% SETUPS_188) "1.8.8" else "1.0.2"
  )
  dir.create(out_dir, recursive = TRUE)
  asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
  if (setup == "opt188") {
    # opt188 exactly as 53b_opt188_unmix.R (route pregate, joint, 1 pass, 1 AF pass), on the chosen raw set
    work_dir <- file.path(out_dir, "work")
    dir.create(work_dir)
    dir_fields <- names(asp)[
      grepl("dir", names(asp)) & vapply(asp, is.character, logical(1))
    ]
    asp[dir_fields] <- lapply(asp[dir_fields], \(d) file.path(work_dir, d))
    OPT188_OBJECTS <- file.path(OPT188_RUN, OPT188_ARM)
    ctrl_dir <- file.path(OPT188_OBJECTS, "controls")
    unmix.folder(
      fcs.dir = raw_dir,
      spectra = read_matrix(file.path(OPT188_OBJECTS, "spectra.csv")),
      asp = asp,
      flow.control = reload.flow.control(
        ctrl_dir,
        file.path(ctrl_dir, "fcs_control_file.csv"),
        asp
      ),
      method = "AutoSpectral",
      af.spectra = read_matrix(file.path(OPT188_OBJECTS, "af_spectra.csv")),
      spectra.variants = readRDS(file.path(OPT188_OBJECTS, "variants.rds")),
      pipeline = "joint",
      speed = "slow",
      n.passes = 1,
      n.af.passes = 1,
      output.dir = out_dir,
      parallel = TRUE,
      verbose = TRUE
    )
  } else if (setup == "afrefine") {
    # af405_afrefine exactly as 87_af405_forced_af.R: the c188 setup with its refined AF set
    control_dir <- file.path(d188, "Reference Group Clean")
    unmix.folder(
      fcs.dir = raw_dir,
      spectra = read_matrix(spectra_188_file),
      asp = asp,
      flow.control = reload.flow.control(
        control_dir,
        file.path(control_dir, "fcs_control_file.csv"),
        asp
      ),
      method = "AutoSpectral",
      af.spectra = read_matrix(AFREFINE_AF),
      spectra.variants = readRDS(file.path(
        d188,
        "figure_spectral_variants",
        "Spectral_variants.rds"
      )),
      pipeline = "joint",
      speed = "slow",
      output.dir = out_dir,
      parallel = TRUE,
      verbose = TRUE
    )
  } else if (setup == "c188") {
    # 1.8.8 setup exactly as 61_af405_proof_unmix.R with the zip variants
    control_dir <- file.path(d188, "Reference Group Clean")
    flow.control <- reload.flow.control(
      control_dir,
      file.path(control_dir, "fcs_control_file.csv"),
      asp
    )
    af <- read_matrix(file.path(
      d188,
      "table_spectra",
      cfg("stamps")$as188_af
    ))
    unmix.folder(
      fcs.dir = raw_dir,
      spectra = read_matrix(spectra_188_file),
      asp = asp,
      flow.control = flow.control,
      method = "AutoSpectral",
      af.spectra = af,
      spectra.variants = readRDS(file.path(
        d188,
        "figure_spectral_variants",
        "Spectral_variants.rds"
      )),
      pipeline = "joint",
      speed = "slow",
      output.dir = out_dir,
      parallel = TRUE,
      verbose = TRUE
    )
  } else {
    # 1.0.2 setups exactly as 06_unmix_c102.R (cells s102 and s188_143803), on the chosen raw set
    if (setup == "s102_c102") {
      control_dir <- file.path(d102, "Reference Group")
      spectra_dir <- file.path(control_dir, "table_spectra")
      spectra_file <- cfg("stamps")$as102_spectra
      af_file <- cfg("stamps")$as102_af
      variants_rds <- file.path(
        control_dir,
        "figure_spectral_variants",
        "Spectral_variants.rds"
      )
    } else {
      control_dir <- file.path(d188, "Reference Group Clean")
      spectra_dir <- file.path(d188, "table_spectra")
      spectra_file <- cfg("stamps")$as188_spectra
      af_file <- cfg("stamps")$as188_af
      variants_rds <- file.path(
        d188,
        "figure_spectral_variants",
        "Spectral_variants.rds"
      )
    }
    flow.control <- reload.flow.control(
      control_dir,
      file.path(control_dir, "fcs_control_file.csv"),
      asp
    )
    spectra <- read.spectra(spectra_file, spectra.dir = spectra_dir)
    spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in ReRunUnmix.R and 06
    unmix.folder(
      fcs.dir = raw_dir,
      spectra = spectra,
      asp = asp,
      flow.control = flow.control,
      method = "AutoSpectral",
      af.spectra = read.spectra(af_file, spectra.dir = spectra_dir),
      weighted = FALSE, # as in 06 (inherited: user zip)
      spectra.variants = readRDS(variants_rds),
      speed = "slow",
      output.dir = out_dir,
      parallel = TRUE,
      threads = NULL
    )
  }
  stopifnot(
    "every raw file unmixed" = setequal(
      list.files(out_dir, pattern = "\\.fcs$"),
      sub(
        "\\.fcs$",
        " AutoSpectral.fcs",
        list.files(raw_dir, pattern = "\\.fcs$")
      )
    )
  )
  writeLines(
    c(paste("setup:", setup), paste("set:", set), paste("raw dir:", raw_dir)),
    file.path(out_dir, "args.txt")
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )
}

# Mode score: TNF and GLUT1 recovery per co-spiked cell -----------------------------------------------
if (mode == "score") {
  RUN_DIR <- run_in(RUN_KEY) # build run: holds cospike_key.csv, score adds its tables here; resolved before library(flowCore): its filter() masks dplyr's inside run_in()
  suppressPackageStartupMessages({
    library(flowCore)
    library(dplyr)
    library(purrr)
    library(readr)
    library(tidyr)
  })
  stopifnot(
    "score outputs must not exist" = !file.exists(file.path(
      RUN_DIR,
      "cospike_cells.csv"
    ))
  )
  key <- read_csv(file.path(RUN_DIR, "cospike_key.csv"), show_col_types = FALSE)
  read_x <- \(path) {
    exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
  }
  unmixed_name <- \(raw_file) sub("\\.fcs$", " AutoSpectral.fcs", raw_file)
  tnf_ch <- paste0(TNF_DYE, "-A")
  glut1_ch <- paste0(GLUT1_DYE, "-A")
  cells <- expand_grid(which_setup = SETUPS, which_set = names(SETS)) |>
    pmap(\(which_setup, which_set) {
      split(key, key$raw_file) |>
        map(\(key_df) {
          f <- unmixed_name(key_df$raw_file[1])
          x_spk <- read_x(file.path(
            DATA,
            "unmix_cells",
            paste0(UNMIX_PREFIX, which_setup, "_", which_set),
            f
          ))
          x_base <- read_x(file.path(
            DATA,
            "unmix_cells",
            paste0(BASE_PREFIX, which_setup, "_base"),
            f
          ))
          stopifnot(
            "spiked and base unmix rows match" = identical(
              x_spk[, KEY_COLS],
              x_base[, KEY_COLS]
            ),
            "unmixed rows pair with the labels (FSC-A)" = max(abs(
              x_base[key_df$raw_row, "FSC-A"] - key_df$FSC_A
            )) <
              0.5
          )
          r <- key_df$raw_row
          key_df |>
            mutate(
              tnf_before = x_base[r, tnf_ch],
              tnf_recovered = x_spk[r, tnf_ch] - x_base[r, tnf_ch],
              glut1_before = x_base[r, glut1_ch],
              glut1_recovered = x_spk[r, glut1_ch] - x_base[r, glut1_ch]
            )
        }) |>
        list_rbind() |>
        mutate(
          setup = which_setup,
          set = which_set,
          glut1_level = SETS[[which_set]]
        )
    }) |>
    list_rbind()
  if (!is.null(DS$conditions)) {
    cells <- inner_join(
      cells,
      select(samples_tbl(), sample, condition),
      by = "sample",
      relationship = "many-to-one",
      unmatched = c("error", "drop")
    )
  }

  # unspiked cells must not move (unmixing is per cell)
  baseline_check <- cells |>
    dplyr::filter(!spiked) |>
    summarise(
      max_abs_diff_tnf = max(abs(tnf_recovered)),
      max_abs_diff_glut1 = max(abs(glut1_recovered)),
      n = n(),
      .by = c(setup, set)
    )
  write_csv(baseline_check, file.path(RUN_DIR, "baseline_check.csv"))

  spiked_cells <- cells |>
    dplyr::filter(spiked) |>
    mutate(
      k_glut1 = COFACTOR * sinh(glut1_level),
      tnf_recovery = tnf_recovered / k_tnf,
      glut1_recovery = if_else(k_glut1 > 0, glut1_recovered / k_glut1, NA_real_)
    )
  write_csv(spiked_cells, file.path(RUN_DIR, "cospike_cells.csv"))
  summary_tbl <- spiked_cells |>
    summarise(
      pct_tnf_lost = 100 * mean(tnf_recovery < LOSS_FRAC),
      median_tnf_recovery = median(tnf_recovery),
      median_glut1_recovery = median(glut1_recovery),
      n = n(),
      .by = c(setup, set, glut1_level, pop, sample, any_of("condition"))
    )
  write_csv(summary_tbl, file.path(RUN_DIR, "cospike_summary.csv"))
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_score.txt")
  )
}

# Mode figures: summary per GLUT1 level and ggcyto raw scatter ---------------------------------------------
# figures_raw: only the per-cell ggcyto figures, as *_v2 files with a two-line title (the first title was clipped)
if (mode %in% c("figures", "figures_raw")) {
  LIB_GGPLOT2 <- cfg("lib_dirs")[basename(cfg("lib_dirs")) == "lib_ggplot2_35"]
  stopifnot(
    "exactly one lib_ggplot2_35 in config lib_dirs" = length(LIB_GGPLOT2) == 1
  )
  RUN_DIR <- run_in(RUN_KEY) # build and score run; before library(flowCore): its filter() masks dplyr's inside run_in()
  .libPaths(c(LIB_GGPLOT2, .libPaths())) # ggplot2 3.5.2 for ggcyto
  suppressPackageStartupMessages({
    library(flowCore)
    library(ggcyto)
    library(ggplot2)
    library(dplyr, exclude = "filter")
    library(purrr)
    library(readr)
  })
  source("experiments/benchmarks/05_bmv-unmix-compare/style_theme.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/style_palettes.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/style_plots.R")
  summary_tbl <- read_csv(
    file.path(RUN_DIR, "cospike_summary.csv"),
    show_col_types = FALSE
  ) |>
    mutate(
      setup_label = factor(SETUP_LABELS[setup], levels = SETUP_LABELS),
      glut1_label = factor(
        paste0("GLUT1 +", glut1_level),
        levels = paste0("GLUT1 +", SETS)
      ),
      donor = sample
    )
  # 24H: donor is the donor id (3 samples each), condition facets the figures; 3D and EXV have no condition level
  tnf_facets <- pop ~ setup_label
  shape_scale <- NULL
  if (!is.null(DS$conditions)) {
    summary_tbl <- summary_tbl |>
      select(-donor) |>
      inner_join(
        select(samples_tbl(), sample, donor),
        by = "sample",
        relationship = "many-to-one",
        unmatched = c("error", "drop")
      )
    donor_ids <- sort(unique(samples_tbl()$donor))
    stopifnot(
      "enough point shapes for the donors" = length(donor_ids) <=
        length(DONOR_SHAPES)
    )
    shape_scale <- scale_shape_manual(
      values = setNames(DONOR_SHAPES[seq_along(donor_ids)], donor_ids)
    )
    tnf_facets <- pop + condition ~ setup_label
  }
  pal_setup <- categorical_palette(levels(summary_tbl$setup_label))
  if (mode == "figures") {
    # Figure 1: % of cells losing >50% of the TNF spike by GLUT1 co-spike level
    p_tnf <- ggplot(
      summary_tbl,
      aes(glut1_label, pct_tnf_lost, colour = setup_label)
    ) +
      stat_summary(
        aes(group = setup_label),
        fun = median,
        geom = "line",
        show.legend = FALSE
      ) +
      geom_point(
        aes(shape = donor),
        position = position_jitter(width = 0.12, height = 0, seed = 1)
      ) +
      scale_colour_manual(values = pal_setup, guide = "none") +
      facet_grid(tnf_facets) +
      shape_scale +
      labs(
        x = paste0(
          "GLUT1 (",
          GLUT1_DYE,
          ") added with ",
          TNF_LABEL,
          " +1 arcsinh"
        ),
        y = paste0("Cells losing >50% of the ", TNF_LABEL, " spike (%)"),
        shape = "Donor",
        title = paste0(TNF_LABEL, " + GLUT1 co-spike: ", TNF_LABEL, " loss")
      )
    save_panel(p_tnf, summary_tbl, RUN_DIR, "fig_cospike_tnf_loss", 183, 100)

    # Figure 2: GLUT1 recovery at the levels where GLUT1 was added
    glut1_df <- dplyr::filter(summary_tbl, glut1_level > 0)
    p_glut1 <- ggplot(
      glut1_df,
      aes(glut1_label, median_glut1_recovery, colour = setup_label)
    ) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
      stat_summary(
        aes(group = setup_label),
        fun = median,
        geom = "line",
        show.legend = FALSE
      ) +
      geom_point(
        aes(shape = donor),
        position = position_jitter(width = 0.12, height = 0, seed = 1)
      ) +
      scale_colour_manual(values = pal_setup, guide = "none") +
      facet_grid(tnf_facets) +
      shape_scale +
      labs(
        x = paste0(
          "GLUT1 (",
          GLUT1_DYE,
          ") added with ",
          TNF_LABEL,
          " +1 arcsinh"
        ),
        y = "Median GLUT1 recovery (fraction of spike)",
        shape = "Donor",
        title = paste0(TNF_LABEL, " + GLUT1 co-spike: GLUT1 recovery")
      )
    save_panel(
      p_glut1,
      glut1_df,
      RUN_DIR,
      "fig_cospike_glut1_recovery",
      183,
      100
    )
  }

  # Figure 3 per population (24H: and condition): ggcyto hexbin TNF before vs after, GLUT1 levels x setups
  cells <- read_csv(
    file.path(RUN_DIR, "cospike_cells.csv"),
    show_col_types = FALSE
  ) |>
    mutate(
      before = asinh(tnf_before / COFACTOR),
      after = asinh((tnf_before + tnf_recovered) / COFACTOR),
      setup_label = factor(SETUP_LABELS[setup], levels = SETUP_LABELS),
      glut1_label = paste0("GLUT1 +", glut1_level)
    )
  axis_lim <- c(-1.5, 4) # display window only
  grid_x <- seq(axis_lim[1], axis_lim[2], length.out = 200)
  curves <- bind_rows(
    tibble(
      x = grid_x,
      y = asinh(sinh(grid_x) + sinh(TNF_LEVEL)),
      curve = "full recovery"
    ),
    tibble(x = grid_x, y = grid_x, curve = "spike lost")
  )
  fig3_panels <- tidyr::expand_grid(
    which_pop = sort(unique(cells$pop)),
    which_cond = if (is.null(DS$conditions)) NA_character_ else DS$conditions
  )
  for (i_panel in seq_len(nrow(fig3_panels))) {
    which_pop <- fig3_panels$which_pop[i_panel]
    which_cond <- fig3_panels$which_cond[i_panel]
    d <- dplyr::filter(cells, pop == which_pop)
    labels_src <- dplyr::filter(summary_tbl, pop == which_pop)
    if (!is.na(which_cond)) {
      d <- dplyr::filter(d, condition == which_cond)
      labels_src <- dplyr::filter(labels_src, condition == which_cond)
    }
    frames <- d |>
      mutate(name = paste(setup_label, glut1_label, sep = "|")) |>
      split(~name) |>
      map(\(x) flowFrame(cbind(before = x$before, after = x$after)))
    fs <- flowSet(frames)
    pd <- tibble(name = sampleNames(fs)) |>
      tidyr::separate_wider_delim(
        name,
        "|",
        names = c("setup_label", "glut1_label"),
        cols_remove = FALSE
      ) |>
      mutate(setup_label = factor(setup_label, levels = SETUP_LABELS))
    pData(fs) <- as.data.frame(pd) |>
      `rownames<-`(pd$name)
    labels <- labels_src |>
      summarise(pct = median(pct_tnf_lost), .by = c(setup_label, glut1_label))
    p <- ggcyto(fs, aes(x = before, y = after)) +
      geom_hex(bins = 50) +
      geom_line(
        data = curves,
        aes(x, y, linetype = curve),
        inherit.aes = FALSE,
        colour = "grey30",
        linewidth = 0.4
      ) +
      geom_label(
        data = labels,
        aes(
          x = -Inf,
          y = Inf,
          label = sprintf("%s lost >50%%: %.1f%%", TNF_LABEL, pct)
        ),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = 7,
        size.unit = "pt",
        fontface = "bold",
        fill = "white",
        alpha = 0.85,
        label.size = 0
      ) +
      facet_grid(glut1_label ~ setup_label, switch = "y") +
      ggcyto_par_set(limits = list(x = axis_lim, y = axis_lim)) +
      scale_fill_gradientn(
        colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
        trans = "sqrt",
        guide = "none"
      ) +
      scale_linetype_manual(
        values = c("full recovery" = "dashed", "spike lost" = "solid"),
        name = NULL
      ) +
      labs(
        x = paste0(TNF_LABEL, ", unspiked unmix (arcsinh)"),
        y = paste0(TNF_LABEL, ", co-spiked unmix (arcsinh)"),
        title = sprintf(
          "%s: %s +1 with GLUT1 co-spike, per cell\n%d donors pooled; %% = donor median",
          if (is.na(which_cond)) which_pop else paste(which_pop, which_cond),
          TNF_LABEL,
          n_distinct(d$sample)
        )
      ) +
      ggpubr::theme_pubr(base_size = 8, border = TRUE) +
      theme(
        strip.text.x = element_text(size = 7),
        strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
        strip.background.y = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(1, "mm")
      )
    stem <- file.path(
      RUN_DIR,
      paste0(
        "fig_cospike_raw_",
        which_pop,
        if (is.na(which_cond)) "" else paste0("_", which_cond),
        if (mode == "figures_raw") "_v2" else ""
      )
    )
    ggsave(
      paste0(stem, ".pdf"),
      as.ggplot(p),
      width = 140,
      height = 130,
      units = "mm",
      device = cairo_pdf
    )
    ggsave(
      paste0(stem, ".png"),
      as.ggplot(p),
      width = 140,
      height = 130,
      units = "mm",
      dpi = 200,
      bg = "white"
    )
  }
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, paste0("sessionInfo_", mode, ".txt"))
  )
}
