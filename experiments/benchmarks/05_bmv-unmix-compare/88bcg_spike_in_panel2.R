# Copy of 88b_spike_in_panel2.R for the common SpectroFlo labels (plans/common_gates_plan.md G5, G10b): only the
#          population source, run dirs and the spiked and unmixed dir prefixes differ. Spiked cells come from
#          run key common-gates-clean3 (labels.csv.gz) (columns CD4T, CD8T, CD3neg_live), keyed by raw_row;
#          the old 88b_spike_in_panel2.R picked them from the per-arm as188_joint gates. Inherited unchanged: +1 arcsinh grid,
#          seed 42, loss cut recovery < 0.5, setups, unmix settings. The spk_<setup>_base unmixes are reused.
# Purpose: spike-in panel 2 (user, 2026-11-10: "extend the spike-in to IFNg, IL-2 and the 7 metabolic
#          markers" on 1.8.8 zip and refined AF). A copy of 88_spike_in_markers.R with a new marker panel and
#          two setups; method, level +1 arcsinh, CD4T and CD8T cells and seed 42 unchanged (README row
#          "Spike-in (grill Q6, script 88)"). The base unmixes spk_c188_base and spk_afrefine_base from 88 are
#          reused (same raw files, same setups). Copied rather than parameterised so 88 stays as it ran.
#          Setups: c188 = as188_joint (1.8.8 zip); afrefine = af405_afrefine (87).
# Run:     from the repo root, on SLURM (experiments/slurm/run_spike_cg.sbatch 88bcg_spike_in_panel2.R <mode>), with SW_DATASET
#          (3d | 24h | exv) and SW_RUN_DATE (YYYY-MM-DD) set; the dataset comes from datasets/<SW_DATASET>.json:
#            build                 spike the raw files, one set per marker; writes the build run dir
#            unmix <setup> <set> [condition]   setup = c188 | afrefine (SETUPS_188 also lists opt188); set = base | a marker name of the config panel;
#                                  condition = 4th argument, required when the dataset has conditions (24H), refused otherwise
#            score                 recovery tables
#            figures               figures from the saved tables
#          SW_SMOKE_DONOR: a donor or sample id (its samples only; 24H: all its conditions); empty = full run.
# Inputs:  <data_dir>/autoSpectral_v1.8.8/{Raw, table_spectra, Reference Group Clean, figure_spectral_variants}  [user zip]
#          <data_dir>/autoSpectral_v1.0.2/Reference Group/  [user zip]
#          <data_dir>/lib_as102/  AutoSpectral 1.0.2 library (config lib_dirs) [06_unmix_c102.R setup]
#          stamp file names (spectra, AF) from config stamps; samples from datasets/<ds>_samples.csv [00_samples.R]
#          run key common-gates-clean3: labels.csv.gz [97d_spectroflo_labels_clean3.R]
#          run key opt188-smoke: opt188 objects, setup opt188 only [53b_opt188_unmix.R]
#          run key spike-in-panel2-cg-build: build run, read by score [this script, mode build]
#          run key spike-in-panel2-cg: score run, read by figures [this script, mode score]
#          config spike_panels: the marker = dye panel [dataset config]
# Outputs: <data_dir>/spiked/<SPIKED_PREFIX><marker>/[<condition>/]*.fcs   spiked raw files
#          <data_dir>/unmix_cells/<UNMIX_PREFIX><setup>_<set>/[<condition>/]*.fcs   unmixed files
#          <data_dir>/unmix_cells/spk_<setup>_base/[<condition>/]*.fcs   unspiked unmixes (set base, no smoke tag); the only script that writes them, 88cg, 88ecg and 88dcg score read them
#          run_out("spike-in-panel2-cg-build")/spike_key_<marker>.csv  sample, raw_row, pop, level, k (raw_file is the path under Raw/)
#          run_out("spike-in-panel2-cg")/spike_cells.csv     per spiked cell: setup, marker, recovered, recovery (24H adds condition)
#          run_out("spike-in-panel2-cg")/spike_summary.csv   per setup, marker, pop, sample, level: median recovery, share losing >50%
#          run_out("spike-in-panel2-cg")/baseline_check.csv  unspiked cells: max |spiked - base|
#          run_out("spike-in-panel2-cg")/fig_spike_<marker>.{pdf,png,csv}  level 1, donors as points
#          run_out("spike-in-panel2-cg")/fig_spike_recovery_ecdf.{pdf,png,csv}  level 1 cells, ECDF
#          run_out("spike-in-panel2-cg")/sessionInfo_<mode>.txt (build, score, figures)
#          Run folders are experiments/runs/<SW_RUN_DATE>_bmv<run_tag>-<key>; register them in datasets/<ds>_runs.csv
#          before the next mode reads them.

# Dataset config -----------------------------------------------------------------------------------
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ---------------------------------------------------------------------------------------
DATA <- data_path()
SMOKE_DONOR <- Sys.getenv("SW_SMOKE_DONOR", "") # donor or sample id for a smoke run: all its samples; empty = full run
SMOKE_TAG <- if (nzchar(SMOKE_DONOR)) "smoke_" else "" # smoke runs write to their own dirs
SMOKE_SUFFIX <- if (nzchar(SMOKE_DONOR)) "-smoke" else ""
KEY_KEY <- paste0("spike-in-panel2-cg-build", SMOKE_SUFFIX) # build outputs: spike keys
RUN_KEY <- paste0("spike-in-panel2-cg", SMOKE_SUFFIX) # score and figures
SPIKED_PREFIX <- paste0(SMOKE_TAG, "2026-11-11_cg_spike_") # spiked raw sets, apart from the old 2026-11-10_spike_*
UNMIX_PREFIX <- paste0(SMOKE_TAG, "spkcg_") # unmixed spiked sets; the base unmixes keep spk_<setup>_base
BASE_PREFIX <- "spk_" # unspiked (base) unmixes read by score; written by 88bcg_spike_in_panel2.R unmix <setup> base, no smoke tag
# opt188 objects (53b_opt188_unmix.R arm opt188_pregate_joint_p1af1) and the refined AF set (87, af405_afrefine)
OPT188_ARM <- "opt188_pregate_neg_same_fsc4e5_p1" # folder of the opt188 objects inside run key opt188-smoke (resolved in the unmix block)
AFREFINE_AF <- data_path("unmix_cells/af405_afrefine/af_refined.csv")
LABELS_KEY <- cfg("label_key") # 3D: SpectroFlo labels of clean3 (97d); 24H and ex vivo: pooled AS 1.8.8 (99a)
# marker = dye (control file); IFNg, IL-2 and the metabolic markers (user, 2026-11-10). The panel is dataset config
# (spike_panels$panel2: 3D = IFNg BV570, IL2 BV605, GLUT1 Alexa Fluor 405, ATP5A Alexa Fluor 488, CS Alexa Fluor 594,
# CPT1A PE-Cy5.5, ATGL PE-Cy7, G6PD APC-Cy7, CD71 Spark YG 593)
MARKERS <- unlist(cfg("spike_panels")$panel2)
if (is.null(MARKERS)) {
  stop(
    "88bcg_spike_in_panel2.R: config spike_panels$panel2 is missing for dataset ",
    DS$dataset,
    call. = FALSE
  )
}
LEVELS <- c(0, 1, 2, 3) # arcsinh units added; 0 = unspiked quarter (62, user-approved grid)
FIG_LEVEL <- 1 # the level the figures show (user, 2026-11-10: +1 arcsinh)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
SEED <- 42 # which cells get which level (62, user-approved); also the figure jitter
LOSS_FRAC <- 0.5 # "losing most of the spike" = recovery below 0.5 (64 readout)
DONOR_SHAPES <- c(16, 17, 15, 18, 8, 3, 4, 7, 9, 10) # figure point shapes in sorted id order; the first five are the 3D shapes, 10 for the 24H donors
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across setups; maps raw rows to cleaned rows
POPS <- c(CD4T = "CD4T", CD8T = "CD8T") # pop name = label column
SETUPS <- c("c188", "afrefine") # user, 2026-11-10
SETUPS_188 <- c("c188", "opt188", "afrefine") # run on the AutoSpectral 1.8.8 library
SETUP_LABELS <- c(
  c188 = "as188_joint",
  s102_c102 = "as102_var",
  s188_c102 = "s188_c102",
  opt188 = "opt188",
  afrefine = "af405_afrefine"
)

if (nzchar(SMOKE_DONOR)) {
  MARKERS <- MARKERS[1] # smoke: first marker only
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

# Mode build: spike the raw files, one set per marker ------------------------------------------------
if (mode == "build") {
  LABELS_DIR <- run_in(LABELS_KEY) # resolved before library(flowCore): its filter() masks dplyr's inside run_in()
  suppressPackageStartupMessages({
    library(flowCore)
    library(flowWorkspace)
    library(dplyr)
    library(purrr)
    library(readr)
    library(tibble)
  })
  KEY_DIR <- run_out(KEY_KEY) # stops when the folder exists
  dir.create(KEY_DIR, recursive = TRUE)
  s188 <- read_matrix(spectra_188_file)

  labels <- read_labels(c(
    "sample",
    "raw_row",
    "FSC_A",
    "CD4T",
    "CD8T",
    "CD3neg_live"
  ))
  # Raw files of the 1.8.8 arm, joined to the sample table on the full file stem (one file per sample)
  raw_files <- tibble(
    raw_file = list.files(
      file.path(d188, "Raw"),
      pattern = "\\.fcs$",
      recursive = TRUE
    )
  ) |>
    mutate(stem = sub("\\.fcs$", "", basename(raw_file))) |>
    inner_join(
      samples_tbl(),
      by = "stem",
      relationship = "one-to-one",
      unmatched = "error"
    )
  if (nzchar(SMOKE_DONOR)) {
    raw_files <- dplyr::filter(
      raw_files,
      donor == SMOKE_DONOR | sample == SMOKE_DONOR
    )
    if (nrow(raw_files) == 0) {
      stop(
        "88bcg_spike_in_panel2.R: SW_SMOKE_DONOR ",
        SMOKE_DONOR,
        " matches no donor or sample of dataset ",
        DS$dataset,
        call. = FALSE
      )
    }
    labels <- dplyr::filter(labels, sample %in% raw_files$sample)
  }

  # CD4T and CD8T raw rows per donor from the SpectroFlo labels, the same for every marker
  pop_rows <- map(set_names(unique(labels$sample)), \(s) {
    raw_file <- raw_files$raw_file[raw_files$sample == s]
    stopifnot("one raw file per sample" = length(raw_file) == 1)
    donor <- labels[labels$sample == s, ]
    imap(POPS, \(label_col, p) {
      tibble(
        raw_row = donor$raw_row[donor[[label_col]]],
        FSC_A = donor$FSC_A[donor[[label_col]]],
        pop = p
      )
    }) |>
      list_rbind() |>
      distinct(raw_row, .keep_all = TRUE) |>
      mutate(raw_file = raw_file)
  })

  # One spiked set per marker; set.seed per marker so every marker spikes the same cells
  iwalk(MARKERS, \(dye, marker) {
    out_fcs <- file.path(DATA, "spiked", paste0(SPIKED_PREFIX, marker))
    stopifnot("spiked fcs dir must not exist" = !dir.exists(out_fcs))
    dir.create(out_fcs, recursive = TRUE)
    set.seed(SEED)
    spike_key <- imap(pop_rows, \(rows, s) {
      ff <- read.FCS(
        file.path(d188, "Raw", rows$raw_file[1]),
        transformation = FALSE,
        truncate_max_range = FALSE
      )
      stopifnot(
        "label raw_row pairs with the raw file (FSC-A)" = max(abs(
          exprs(ff)[rows$raw_row, "FSC-A"] - rows$FSC_A
        )) <
          0.5 # float32 spacing at 1e6 is 0.06 to 0.12
      )
      key <- rows |>
        select(raw_row, FSC_A, pop, raw_file) |>
        mutate(
          level = sample(rep_len(LEVELS, n())),
          k = COFACTOR * sinh(level),
          sample = s
        )
      x <- exprs(ff)
      x[key$raw_row, colnames(s188)] <- x[key$raw_row, colnames(s188)] +
        outer(key$k, s188[dye, ])
      exprs(ff) <- x
      out_file <- file.path(out_fcs, rows$raw_file[1])
      dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE) # condition subfolder in 24H
      write.FCS(ff, out_file)
      key
    }) |>
      list_rbind() |>
      mutate(marker = marker, dye = dye)
    write_csv(
      spike_key,
      file.path(KEY_DIR, paste0("spike_key_", marker, ".csv"))
    )
  })
  writeLines(
    capture.output(sessionInfo()),
    file.path(KEY_DIR, "sessionInfo_build.txt")
  )
}

# Mode unmix: one setup on one raw set ----------------------------------------------------------------
if (mode == "unmix") {
  setup <- cli_args[[2]]
  set <- cli_args[[3]]
  stopifnot(
    "known setup and set" = setup %in%
      SETUPS &&
      set %in% c("base", names(MARKERS))
  )
  # 24H: one task per condition subfolder (4th job argument); 3D and EXV have no condition level
  cond <- if (length(cli_args) >= 4) cli_args[[4]] else NULL
  if (is.null(DS$conditions) && !is.null(cond)) {
    stop(
      "88bcg_spike_in_panel2.R: dataset ",
      DS$dataset,
      " has no conditions; drop the 4th argument ",
      cond,
      call. = FALSE
    )
  }
  if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
    stop(
      "88bcg_spike_in_panel2.R: unmix needs a 4th argument, a condition of dataset ",
      DS$dataset,
      ": ",
      paste(DS$conditions, collapse = ", "),
      call. = FALSE
    )
  }
  raw_dir <- if (set == "base") {
    file.path(d188, "Raw")
  } else {
    file.path(DATA, "spiked", paste0(SPIKED_PREFIX, set))
  }
  raw_dir <- paste(c(raw_dir, cond), collapse = "/")
  out_dir <- normalizePath(file.path(DATA, "unmix_cells"), mustWork = TRUE) |>
    file.path(
      if (set == "base") {
        paste0(BASE_PREFIX, setup, "_base")
      } else {
        paste0(UNMIX_PREFIX, setup, "_", set)
      }
    )
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

# Mode score: recovery per spiked cell, tables and figures ----------------------------------------------
if (mode == "score") {
  KEY_DIR <- run_in(KEY_KEY) # build run: spike keys; resolved before library(flowCore): its filter() masks dplyr's inside run_in()
  suppressPackageStartupMessages({
    library(flowCore)
    library(dplyr)
    library(ggplot2)
    library(purrr)
    library(readr)
    library(tidyr)
    library(tibble)
  })
  source("experiments/benchmarks/05_bmv-unmix-compare/style_theme.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/style_palettes.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/style_plots.R")
  RUN_DIR <- run_out(RUN_KEY) # stops when the folder exists
  dir.create(RUN_DIR, recursive = TRUE)
  read_x <- \(path) {
    exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
  }
  unmixed_name <- \(raw_file) sub("\\.fcs$", " AutoSpectral.fcs", raw_file)

  # Spiked minus unspiked marker value, per cell, setup and marker
  cells <- expand_grid(setup = SETUPS, marker = names(MARKERS)) |>
    pmap(\(setup, marker) {
      which_setup <- setup
      which_marker <- marker
      channel <- paste0(MARKERS[[which_marker]], "-A")
      key <- read_csv(
        file.path(KEY_DIR, paste0("spike_key_", which_marker, ".csv")),
        show_col_types = FALSE
      )
      split(key, key$raw_file) |>
        map(\(key_df) {
          f <- unmixed_name(key_df$raw_file[1])
          x_spk <- read_x(file.path(
            DATA,
            "unmix_cells",
            paste0(UNMIX_PREFIX, which_setup, "_", which_marker),
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
          key_df |>
            mutate(
              recovered = x_spk[raw_row, channel] - x_base[raw_row, channel],
              base_value = x_base[raw_row, channel]
            )
        }) |>
        list_rbind() |>
        mutate(setup = which_setup)
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

  # Unspiked quarter: the spike must not move them (unmixing is per cell)
  baseline_check <- cells |>
    dplyr::filter(level == 0) |>
    summarise(
      max_abs_diff = max(abs(recovered)),
      n = n(),
      .by = c(setup, marker)
    )
  write_csv(baseline_check, file.path(RUN_DIR, "baseline_check.csv"))

  spiked <- cells |>
    dplyr::filter(level > 0) |>
    mutate(recovery = recovered / k)
  write_csv(spiked, file.path(RUN_DIR, "spike_cells.csv"))
  summary_tbl <- spiked |>
    summarise(
      n_cells = n(),
      median_recovery = median(recovery),
      share_lost = mean(recovery < LOSS_FRAC),
      .by = c(setup, marker, pop, sample, level, any_of("condition"))
    )
  write_csv(summary_tbl, file.path(RUN_DIR, "spike_summary.csv"))

  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_score.txt")
  )
}

# Mode figures: draw from the saved tables, so a redraw does not rescore ---------------------------------
if (mode == "figures") {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(purrr)
    library(readr)
    library(tidyr)
  })
  source("experiments/benchmarks/05_bmv-unmix-compare/style_theme.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/style_palettes.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/style_plots.R")
  RUN_DIR <- run_in(RUN_KEY) # score run
  summary_tbl <- read_csv(
    file.path(RUN_DIR, "spike_summary.csv"),
    show_col_types = FALSE
  )
  spiked <- read_csv(
    file.path(RUN_DIR, "spike_cells.csv"),
    show_col_types = FALSE
  )

  # Figure per marker: level 1, donors as points, setup on x, two readouts by population
  setup_pal <- categorical_palette(unname(SETUP_LABELS))
  setup_lines <- setNames(
    c("solid", "dashed", "dotted", "dotdash", "longdash"),
    names(setup_pal)
  ) # second encoding of setup
  if (!is.null(DS$conditions)) {
    summary_tbl <- inner_join(
      summary_tbl,
      select(samples_tbl(), sample, donor),
      by = "sample",
      relationship = "many-to-one",
      unmatched = c("error", "drop")
    )
  }
  shape_col <- if (is.null(DS$conditions)) "sample" else "donor" # point shape per sample (3D, EXV) or per donor (24H)
  shape_ids <- sort(unique(samples_tbl()[[shape_col]]))
  donor_shapes <- setNames(DONOR_SHAPES[seq_along(shape_ids)], shape_ids)
  fig_facets <- if (is.null(DS$conditions)) {
    readout ~ pop
  } else {
    readout ~ pop + condition
  }
  ecdf_facets <- if (is.null(DS$conditions)) {
    marker ~ pop
  } else {
    marker ~ pop + condition
  }
  fig_frame <- summary_tbl |>
    filter(level == FIG_LEVEL) |>
    mutate(
      setup = factor_by_pal(SETUP_LABELS[setup], setup_pal),
      `Median recovery (fraction of spike)` = median_recovery,
      `Cells losing >50% of spike (%)` = 100 * share_lost
    ) |>
    pivot_longer(
      c(
        `Median recovery (fraction of spike)`,
        `Cells losing >50% of spike (%)`
      ),
      names_to = "readout",
      values_to = "value"
    )
  stopifnot(
    "every donor has a shape" = length(shape_ids) <= length(DONOR_SHAPES) &&
      all(fig_frame[[shape_col]] %in% names(donor_shapes))
  )
  n_per_setup <- fig_frame |>
    filter(readout == first(readout), pop == first(pop))
  iwalk(MARKERS, \(dye, which_marker) {
    panel_data <- filter(fig_frame, marker == which_marker)
    p <- ggplot(panel_data, aes(setup, value)) +
      stat_summary(
        aes(colour = setup),
        fun = median,
        fun.min = median,
        fun.max = median,
        geom = "crossbar",
        width = 0.5
      ) +
      geom_point(
        aes(colour = setup, shape = .data[[shape_col]]),
        position = position_jitter(width = 0.15, height = 0, seed = SEED),
        size = 1.5
      ) +
      facet_grid(fig_facets, scales = "free_y", switch = "y") +
      scale_colour_manual(values = setup_pal, guide = "none") +
      scale_shape_manual(values = donor_shapes, name = "Donor") +
      scale_x_discrete(
        labels = n_labels(filter(n_per_setup, marker == which_marker)$setup)
      ) +
      labs(
        x = NULL,
        y = NULL,
        title = sprintf(
          "%s (%s) spike +%d arcsinh",
          which_marker,
          dye,
          FIG_LEVEL
        )
      ) +
      theme(strip.placement = "outside", strip.background = element_blank())
    save_panel(
      p,
      panel_data,
      RUN_DIR,
      paste0("fig_spike_", which_marker),
      183,
      110
    )
  })

  # Distribution figure: cumulative recovery of every level-1 cell, pooled over donors, by setup.
  # An ECDF, because many 1.8.8 cells recover exactly 1.0 and a density spike hides the other setups.
  dist_frame <- spiked |>
    filter(level == FIG_LEVEL) |>
    mutate(
      setup = factor_by_pal(SETUP_LABELS[setup], setup_pal),
      marker = factor(marker, levels = names(MARKERS))
    )
  p_dist <- ggplot(
    dist_frame,
    aes(recovery, colour = setup, linetype = setup)
  ) +
    geom_vline(
      xintercept = c(LOSS_FRAC, 1),
      linetype = "dashed",
      colour = "grey60"
    ) +
    stat_ecdf(geom = "step", linewidth = 0.4) +
    facet_grid(ecdf_facets, switch = "y") +
    scale_colour_manual(values = setup_pal, name = NULL) +
    scale_linetype_manual(values = setup_lines, name = NULL) +
    coord_cartesian(xlim = c(-0.5, 1.5)) +
    labs(
      x = "Recovered fraction of the +1 spike",
      y = "Cumulative fraction of spiked cells",
      title = "Spike recovery per cell"
    ) +
    theme(strip.placement = "outside", strip.background = element_blank())
  save_panel(p_dist, dist_frame, RUN_DIR, "fig_spike_recovery_ecdf", 183, 150)
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_figures.txt")
  )
}
