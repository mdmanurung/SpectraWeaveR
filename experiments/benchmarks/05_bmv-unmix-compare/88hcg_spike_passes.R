# Copy of 88h_spike_passes.R for the common SpectroFlo labels (plans/common_gates_plan.md G5, G10b). Only the inputs (the spike keys, spiked sets and 1/1 unmixes of 88cg and 88ecg)
#          and the output dirs differ; method, levels, seed 42, loss cut and setups are inherited unchanged.
# Purpose: does raising AutoSpectral n.passes / n.af.passes change spike recovery for Th2, IL-10 and TNF?
#          Grill 2026-11-11 (README row "Spike passes (grill 2026-11-11; script 88h)"): pass settings 5/1,
#          1/5 and 5/5 (n.passes/n.af.passes) on the as188_joint (1.8.8 zip) setup, joint pipeline, speed
#          "slow"; the 1/1 reference is the existing spk_c188_* / spkneg_c188_* unmixes of 88 and 88e (Q4).
#          Spiked raw files, levels 0-3, seed 42 and the loss cut (recovery < 0.5) are inherited from 88
#          (CD4T, CD8T) and 88e (live CD3-) (Q6). Markers Th2, IL-10 and TNF as control (Q2). Descriptive:
#          the final setting stays 1/1 unless the user decides after seeing the result (Q5).
# Run:     from the repo root, on SLURM, with SW_DATASET=<3d|24h|exv> and SW_RUN_DATE=<YYYY-MM-DD> set (both required;
#          the dataset sets the data dir, stamps, run folders and spike markers):
#          experiments/slurm/run_spike_cg.sbatch 88hcg_spike_passes.R <mode>
#            unmix <setting> <set> [<condition>]
#                                    setting = p5_af1 | p1_af5 | p5_af5;
#                                    set = base | <marker> | neg_<marker> for the markers of cfg("spike_panels")$panel1 (the spike keys and spiked sets exist only for panel1);
#                                    condition = iRBC | SEB | uRBC, required for datasets with conditions (24H):
#                                    one task per condition folder, raw and unmixed folders get /<condition>
#            score                   recovery tables, all settings incl. the 1/1 reference
#            figures                 summary figures from the saved tables (Nature spec, lib/theme.R)
#            raw                     per-cell hexbin companions (ggcyto, private ggplot2 3.5.2 lib)
#            drift                   unspiked-cell shift table (does a cell unmix independently of the file?)
# Inputs:  data_path("autoSpectral_v1.8.8")/{Raw, table_spectra, Reference Group Clean,
#            figure_spectral_variants}  [user zip]; spectra and AF file names from cfg("stamps")
#          data_path("spiked")/<tag>2026-11-11_cg_spike_<marker>/, ..._cg_spike_cd3neg_<marker>/  spiked raw [88cg, 88ecg]
#          data_path("unmix_cells")/spk_c188_{base,<marker>}/, spkcg_c188_<marker>/, spkcgneg_c188_<marker>/
#            1/1 reference unmixes [88, 88cg, 88ecg]; found recursively and matched to samples_tbl() on the
#            file stem (24H keeps one folder per condition)
#          run_in("spike-in-cg-build")/spike_key_<marker>.csv        T-cell keys [88cg build]
#          run_in("spike-in-cd3neg-cg-build")/spike_key_<marker>.csv CD3- keys [88ecg build]
#          (registry keys of datasets/<ds>_runs.csv; "-smoke" is added to both keys and to the run key in a smoke run)
# Outputs: data_path("unmix_cells")/spkcgpass_<setting>_<set>[/<condition>]/*.fcs  unmixed files
#          run_out("spike-passes-cg") = experiments/runs/<SW_RUN_DATE>_bmv<run_tag>-spike-passes-cg/, created by
#          mode score; drift, figures and raw add to it and stop when it is missing
#          <run dir>/baseline_check.csv  unspiked cells: max |spiked - base|
#          <run dir>/drift_summary.csv   unspiked cells: max shift and % moved
#                                        >0.05 / >0.5 arcsinh, own and any channel
#          <run dir>/spike_cells.csv     per spiked cell: setting, marker,
#                                        pop, level, k, recovered, recovery
#          <run dir>/spike_summary.csv   per setting, marker, pop, donor,
#                                        level: median recovery, share lost
#          <run dir>/fig_passes_<marker>.{pdf,png,csv}   level 1, donors
#          <run dir>/raw/fig_raw_<marker>_<pop>.{pdf,png,csv}  per cell
#          <run dir>/sessionInfo_<mode>.txt
#          Datasets with conditions (24H) add donor and condition columns to the figure source csv files and facet by
#          condition.

# the dataset config; sourced before the mode libraries, so dplyr is attached first and every filter below is
# written dplyr::filter() (flowCore, loaded later, would otherwise mask it)
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
filter <- dplyr::filter # dataset.R calls filter() from the search path; flowCore, attached in the modes, would mask it

# Parameters ---------------------------------------------------------------------------------------
DATA <- data_path()
cli_args <- commandArgs(trailingOnly = TRUE)
mode <- cli_args[[1]]
SMOKE_DONOR <- Sys.getenv("SW_SMOKE_DONOR", "") # one whole donor (all its rows) for a smoke run; empty = full run
SMOKE_TAG <- if (nzchar(SMOKE_DONOR)) "smoke_" else "" # smoke runs write to their own dirs
SMOKE_SUFFIX <- if (nzchar(SMOKE_DONOR)) "-smoke" else ""
RUN_KEY <- paste0("spike-passes-cg", SMOKE_SUFFIX)
# score creates the run dir (run_out stops when it exists); drift, figures and raw add to the same dir
RUN_DIR <- if (mode == "score") {
  run_out(RUN_KEY)
} else {
  file.path(
    "experiments/runs",
    paste0(Sys.getenv("SW_RUN_DATE"), "_bmv", cfg("run_tag"), "-", RUN_KEY)
  )
}
if (mode %in% c("drift", "figures", "raw") && !dir.exists(RUN_DIR)) {
  stop(
    "88hcg_spike_passes.R: run dir ",
    RUN_DIR,
    " is missing; run mode score first with the same SW_RUN_DATE",
    call. = FALSE
  )
}
KEY_DIRS <- if (mode %in% c("score", "drift")) {
  c(
    t = run_in(paste0("spike-in-cg-build", SMOKE_SUFFIX)), # CD4T, CD8T keys [88cg]
    neg = run_in(paste0("spike-in-cd3neg-cg-build", SMOKE_SUFFIX)) # CD3- keys [88ecg]
  )
}
# grill Q1: n.passes / n.af.passes; the 1/1 reference reuses the c188 unmixes (Q4)
SETTINGS <- list(
  p1_af1 = c(n_passes = 1, n_af_passes = 1),
  p5_af1 = c(n_passes = 5, n_af_passes = 1),
  p1_af5 = c(n_passes = 1, n_af_passes = 5),
  p5_af5 = c(n_passes = 5, n_af_passes = 5)
)
SETTING_LABELS <- c(
  p1_af1 = "1 / 1 (as188_joint)",
  p5_af1 = "5 / 1",
  p1_af5 = "1 / 5",
  p5_af5 = "5 / 5"
)
MARKERS <- unlist(cfg("spike_panels")$panel1) # grill Q2; same panel as 88cg (3D: TNF BV510, Th2 eFluor 450, IL10 BB700); marker label = dye
if (is.null(MARKERS)) {
  stop(
    "88hcg_spike_passes.R: needs spike_panels$panel1 (marker label = dye), missing in dataset ",
    DS$dataset,
    call. = FALSE
  )
}
DONOR_SHAPES <- c(16, 17, 15, 18, 8, 0, 1, 2, 5, 6) # first five as before; 24H has 10 donors
if (nzchar(SMOKE_DONOR)) {
  MARKERS <- MARKERS[1] # smoke: first marker only
}
FIG_LEVEL <- 1 # inherited (Q6): figures show the +1 arcsinh cells
LOSS_FRAC <- 0.5 # inherited (Q6): "lost" = recovery below 0.5
SEED <- 42 # inherited (Q6): figure jitter (the spiked files already fix the cell assignment)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across unmixes of the same raw rows

d188 <- file.path(DATA, "autoSpectral_v1.8.8")
read_matrix <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
# unmixed folder of one setting and set. Base unmixes (whole Raw files) are the old ones, reused: spk_c188_base
# for 1/1 and spkpass_<setting>_base. The 1/1 spiked unmixes are 88cg's and 88ecg's c188 output.
unmix_dir <- \(setting, set) {
  if (set == "base") {
    sub_dir <- if (setting == "p1_af1") {
      "spk_c188_base"
    } else {
      paste0("spkpass_", setting, "_base")
    }
  } else if (setting == "p1_af1") {
    sub_dir <- if (startsWith(set, "neg_")) {
      paste0(SMOKE_TAG, "spkcgneg_c188_", sub("^neg_", "", set))
    } else {
      paste0(SMOKE_TAG, "spkcg_c188_", set)
    }
  } else {
    sub_dir <- paste0(SMOKE_TAG, "spkcgpass_", setting, "_", set)
  }
  file.path(DATA, "unmix_cells", sub_dir)
}
# unmixed file of every raw file in a spike key, matched on the full file stem and to samples_tbl().
# Recursive, because 24H keeps one folder per condition; one file per sample or the join stops.
unmixed_paths <- \(setting, set, key) {
  want <- distinct(key, sample, raw_file) |>
    mutate(stem = sub("\\.fcs$", "", basename(raw_file)))
  stopifnot(
    "88hcg: every key raw file is a registered sample (samples_tbl stem)" = all(
      want$stem %in% samples_tbl()$stem
    )
  )
  want <- want |>
    inner_join(
      samples_tbl() |>
        dplyr::filter(sample %in% key$sample) |>
        select(sample, stem),
      by = c("sample", "stem"),
      relationship = "one-to-one",
      unmatched = "error"
    )
  found <- tibble(
    path = list.files(
      unmix_dir(setting, set),
      pattern = " AutoSpectral\\.fcs$",
      recursive = TRUE,
      full.names = TRUE
    )
  ) |>
    mutate(stem = sub(" AutoSpectral\\.fcs$", "", basename(path))) |>
    dplyr::filter(stem %in% want$stem)
  stopifnot(
    "88hcg: one unmixed file per raw file in the unmix dir" = setequal(
      found$stem,
      want$stem
    )
  )
  inner_join(
    found,
    want,
    by = "stem",
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
    select(raw_file, path)
}

# Mode unmix: one pass setting on one raw set -------------------------------------------------------
if (mode == "unmix") {
  setting <- cli_args[[2]]
  set <- cli_args[[3]]
  stopifnot(
    "known setting and set" = setting %in%
      setdiff(names(SETTINGS), "p1_af1") &&
      set %in% c("base", names(MARKERS), paste0("neg_", names(MARKERS)))
  )
  # datasets with conditions keep one raw and one unmixed folder per condition; one task per condition
  condition <- if (is.null(DS$conditions)) NULL else cli_args[4]
  if (!is.null(DS$conditions) && !condition %in% DS$conditions) {
    stop(
      "88hcg_spike_passes.R: dataset ",
      DS$dataset,
      " needs a 4th argument <condition>, one of ",
      paste(DS$conditions, collapse = ", "),
      call. = FALSE
    )
  }
  raw_dir <- if (set == "base") {
    file.path(d188, "Raw")
  } else if (startsWith(set, "neg_")) {
    file.path(
      DATA,
      "spiked",
      paste0(SMOKE_TAG, "2026-11-11_cg_spike_cd3neg_", sub("^neg_", "", set))
    )
  } else {
    file.path(DATA, "spiked", paste0(SMOKE_TAG, "2026-11-11_cg_spike_", set))
  }
  out_dir <- unmix_dir(setting, set)
  if (!is.null(condition)) {
    raw_dir <- file.path(raw_dir, condition)
    out_dir <- file.path(out_dir, condition)
  }
  stopifnot("out dir must not exist" = !dir.exists(out_dir))
  suppressPackageStartupMessages({
    library(AutoSpectral)
    library(flowCore)
  })
  stopifnot(
    "AutoSpectral 1.8.8" = as.character(packageVersion("AutoSpectral")) ==
      "1.8.8"
  )
  dir.create(out_dir, recursive = TRUE)
  asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
  # the c188 setup exactly as 88 (61_af405_proof_unmix.R with the zip variants), passes changed only
  control_dir <- file.path(d188, "Reference Group Clean")
  unmix.folder(
    fcs.dir = raw_dir,
    spectra = read_matrix(file.path(
      d188,
      "table_spectra",
      cfg("stamps")$as188_spectra
    )),
    asp = asp,
    flow.control = reload.flow.control(
      control_dir,
      file.path(control_dir, "fcs_control_file.csv"),
      asp
    ),
    method = "AutoSpectral",
    af.spectra = read_matrix(file.path(
      d188,
      "table_spectra",
      cfg("stamps")$as188_af
    )),
    spectra.variants = readRDS(file.path(
      d188,
      "figure_spectral_variants",
      "Spectral_variants.rds"
    )),
    pipeline = "joint",
    speed = "slow",
    n.passes = SETTINGS[[setting]][["n_passes"]],
    n.af.passes = SETTINGS[[setting]][["n_af_passes"]],
    output.dir = out_dir,
    parallel = TRUE,
    verbose = TRUE
  )
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
    c(
      paste("setting:", setting),
      paste("set:", set),
      if (!is.null(condition)) paste("condition:", condition),
      paste("raw dir:", raw_dir)
    ),
    file.path(out_dir, "args.txt")
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )
}

# Mode score: recovery per spiked cell, every setting ----------------------------------------------------
if (mode == "score") {
  suppressPackageStartupMessages({
    library(flowCore)
    library(dplyr)
    library(purrr)
    library(readr)
    library(tidyr)
  })
  stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
  dir.create(RUN_DIR, recursive = TRUE)
  read_x <- \(path) {
    exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
  }

  # spiked minus unspiked marker value of the same cell, per setting, marker and population group
  cells <- expand_grid(
    setting = names(SETTINGS),
    marker = names(MARKERS),
    group = names(KEY_DIRS)
  ) |>
    pmap(\(setting, marker, group) {
      which_marker <- marker
      channel <- paste0(MARKERS[[which_marker]], "-A")
      set <- if (group == "neg") paste0("neg_", which_marker) else which_marker
      key <- read_csv(
        file.path(
          KEY_DIRS[[group]],
          paste0("spike_key_", which_marker, ".csv")
        ),
        show_col_types = FALSE
      )
      files_spk <- unmixed_paths(setting, set, key)
      files_base <- unmixed_paths(setting, "base", key)
      split(key, key$raw_file) |>
        map(\(key_df) {
          f <- key_df$raw_file[1]
          x_spk <- read_x(files_spk$path[files_spk$raw_file == f])
          x_base <- read_x(files_base$path[files_base$raw_file == f])
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
        mutate(setting = setting)
    }) |>
    list_rbind()

  # unspiked quarter: the spike must not move them (unmixing is per cell)
  baseline_check <- cells |>
    dplyr::filter(level == 0) |>
    summarise(
      max_abs_diff = max(abs(recovered)),
      n = n(),
      .by = c(setting, marker, pop)
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
      .by = c(setting, marker, pop, sample, level)
    )
  write_csv(summary_tbl, file.path(RUN_DIR, "spike_summary.csv"))
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_score.txt")
  )
}

# Mode drift: do unspiked cells move when other cells in the file are spiked? ----------------------------
# The spike method assumes each cell unmixes independently of the others in its file. For every unspiked
# (level 0) cell: the largest shift over all fluorescence channels and the shift in the spiked dye's own
# channel, in arcsinh units, spiked-set unmix vs base unmix. Thresholds below are reporting bins only.
if (mode == "drift") {
  suppressPackageStartupMessages({
    library(flowCore)
    library(dplyr)
    library(purrr)
    library(readr)
    library(tidyr)
  })
  DRIFT_BINS <- c(0.05, 0.5) # arcsinh shift bins for the share of cells moved; reporting only
  out_file <- file.path(RUN_DIR, "drift_summary.csv")
  stopifnot("drift table must not exist" = !file.exists(out_file))
  read_x <- \(path) {
    exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
  }
  drift <- expand_grid(
    setting = names(SETTINGS),
    marker = names(MARKERS),
    group = names(KEY_DIRS)
  ) |>
    pmap(\(setting, marker, group) {
      which_marker <- marker
      set <- if (group == "neg") paste0("neg_", which_marker) else which_marker
      key <- read_csv(
        file.path(
          KEY_DIRS[[group]],
          paste0("spike_key_", which_marker, ".csv")
        ),
        show_col_types = FALSE
      ) |>
        dplyr::filter(level == 0)
      files_spk <- unmixed_paths(setting, set, key)
      files_base <- unmixed_paths(setting, "base", key)
      split(key, key$raw_file) |>
        map(\(key_df) {
          f <- key_df$raw_file[1]
          x_spk <- read_x(files_spk$path[files_spk$raw_file == f])
          x_base <- read_x(files_base$path[files_base$raw_file == f])
          channels <- setdiff(
            grep("-A$", colnames(x_spk), value = TRUE),
            c("FSC-A", "SSC-A")
          )
          shift <- abs(
            asinh(x_spk[key_df$raw_row, channels, drop = FALSE] / COFACTOR) -
              asinh(x_base[key_df$raw_row, channels, drop = FALSE] / COFACTOR)
          )
          tibble(
            sample = key_df$sample,
            pop = key_df$pop,
            own = shift[, paste0(MARKERS[[which_marker]], "-A")],
            any = apply(shift, 1, max)
          )
        }) |>
        list_rbind() |>
        mutate(setting = setting, marker = which_marker)
    }) |>
    list_rbind()
  drift_summary <- drift |>
    summarise(
      n_unspiked = n(),
      max_any = max(any),
      pct_own_gt_bin1 = 100 * mean(own > DRIFT_BINS[1]),
      pct_own_gt_bin2 = 100 * mean(own > DRIFT_BINS[2]),
      pct_any_gt_bin1 = 100 * mean(any > DRIFT_BINS[1]),
      pct_any_gt_bin2 = 100 * mean(any > DRIFT_BINS[2]),
      .by = c(setting, marker, pop)
    )
  write_csv(drift_summary, out_file)
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_drift.txt")
  )
}

# Mode figures: % lost and median recovery at +1, one figure per marker ---------------------------------
if (mode == "figures") {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tidyr)
  })
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")
  summary_tbl <- read_csv(
    file.path(RUN_DIR, "spike_summary.csv"),
    show_col_types = FALSE
  )
  # datasets with conditions: add donor and condition; shapes then follow the donor (24H has 3 samples per donor)
  if (!is.null(DS$conditions)) {
    summary_tbl <- summary_tbl |>
      inner_join(
        samples_tbl() |>
          dplyr::filter(sample %in% summary_tbl$sample) |>
          select(sample, donor, condition),
        by = "sample",
        relationship = "many-to-one",
        unmatched = "error"
      )
  }
  shape_by <- if (is.null(DS$conditions)) "sample" else "donor"
  fig_rows <- if (is.null(DS$conditions)) vars(pop) else vars(pop, condition)
  donors <- sort(unique(summary_tbl[[shape_by]]))
  stopifnot(
    "enough point shapes for the donors" = length(donors) <=
      length(DONOR_SHAPES)
  )
  # long format, one row per setting, population, donor and readout; reference on top
  plot_df <- summary_tbl |>
    dplyr::filter(level == FIG_LEVEL) |>
    mutate(
      `Cells losing >50% of spike (%)` = 100 * share_lost,
      `Median recovery (fraction of spike)` = median_recovery
    ) |>
    pivot_longer(
      c(
        `Cells losing >50% of spike (%)`,
        `Median recovery (fraction of spike)`
      ),
      names_to = "readout"
    ) |>
    mutate(
      setting = factor(SETTING_LABELS[setting], levels = rev(SETTING_LABELS)),
      pop = factor(pop, levels = c("CD4T", "CD8T", "CD3neg")),
      across(all_of(shape_by), \(x) factor(x, levels = donors))
    )
  setting_pal <- categorical_palette(unname(SETTING_LABELS))
  donor_shapes <- setNames(DONOR_SHAPES[seq_along(donors)], donors)
  for (which_marker in names(MARKERS)) {
    d <- dplyr::filter(plot_df, marker == which_marker)
    p <- ggplot(d, aes(value, setting)) +
      stat_summary(
        aes(colour = setting),
        fun = median,
        fun.min = median,
        fun.max = median,
        geom = "crossbar",
        width = 0.5,
        linewidth = FIG_LINE / 2,
        orientation = "y"
      ) +
      geom_point(
        aes(colour = setting, shape = .data[[shape_by]]),
        position = position_jitter(width = 0, height = 0.15, seed = SEED),
        size = 1.3
      ) +
      facet_grid(
        rows = fig_rows,
        cols = vars(readout),
        scales = "free_x",
        axes = "all",
        axis.labels = "margins",
        labeller = labeller(readout = label_wrap_gen(28))
      ) +
      scale_colour_manual(values = setting_pal, guide = "none") +
      scale_shape_manual(values = donor_shapes, name = "Donor") +
      labs(
        x = NULL,
        y = "n.passes / n.af.passes",
        title = sprintf(
          "%s (%s) spike +%d arcsinh",
          which_marker,
          MARKERS[[which_marker]],
          FIG_LEVEL
        )
      ) +
      value_grid("x")
    save_panel(
      p,
      d,
      RUN_DIR,
      paste0("fig_passes_", which_marker),
      width_mm = 183,
      height_mm = 110 * max(1, length(DS$conditions)) # one block of rows per condition
    )
  }
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_figures.txt")
  )
}

# Mode raw: per-cell before/after hexbins behind each summary figure ------------------------------------
if (mode == "raw") {
  .libPaths(c(
    grep("lib_ggplot2_35", cfg("lib_dirs"), value = TRUE),
    .libPaths()
  )) # ggcyto needs 3.5.2
  suppressPackageStartupMessages({
    library(flowCore)
    library(ggcyto)
    library(dplyr)
    library(purrr)
    library(readr)
  })
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
  raw_out <- file.path(RUN_DIR, "raw")
  stopifnot("raw dir must not exist" = !dir.exists(raw_out))
  dir.create(raw_out)
  AXIS <- c(-1.5, 4) # display window only (as 88c)
  BINS <- 50 # hexagons per panel; display only
  cells <- read_csv(
    file.path(RUN_DIR, "spike_cells.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(level == FIG_LEVEL) |>
    mutate(
      before = asinh(base_value / COFACTOR),
      after = asinh((base_value + recovered) / COFACTOR),
      setting_label = factor(SETTING_LABELS[setting], levels = SETTING_LABELS)
    )
  pct_lost <- cells |>
    summarise(
      pct_lost = 100 * mean(recovery < LOSS_FRAC),
      .by = c(marker, pop, sample, setting_label)
    )
  # full recovery adds sinh(1) in linear units; a lost spike stays on the diagonal
  grid_x <- seq(AXIS[1], AXIS[2], length.out = 200)
  curves <- bind_rows(
    tibble(
      x = grid_x,
      y = asinh(sinh(grid_x) + sinh(FIG_LEVEL)),
      curve = "full recovery"
    ),
    tibble(x = grid_x, y = grid_x, curve = "spike lost")
  )
  distinct(cells, marker, pop) |>
    pwalk(\(marker, pop) {
      which_marker <- marker
      which_pop <- pop
      d <- dplyr::filter(
        cells,
        .data$marker == which_marker,
        .data$pop == which_pop
      )
      frames <- d |>
        mutate(name = paste(setting_label, sample, sep = "|")) |>
        split(~name) |>
        map(\(x) flowFrame(cbind(before = x$before, after = x$after)))
      fs <- flowSet(frames)
      pd <- tibble(name = sampleNames(fs)) |>
        tidyr::separate_wider_delim(
          name,
          "|",
          names = c("setting_label", "sample"),
          cols_remove = FALSE
        ) |>
        mutate(setting_label = factor(setting_label, levels = SETTING_LABELS))
      pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)
      labels <- dplyr::filter(
        pct_lost,
        .data$marker == which_marker,
        .data$pop == which_pop
      )
      p <- ggcyto(fs, aes(x = before, y = after)) +
        geom_hex(bins = BINS) +
        geom_line(
          data = curves,
          aes(x, y, linetype = curve),
          inherit.aes = FALSE,
          colour = "grey20",
          linewidth = FIG_LINE
        ) +
        geom_text(
          data = labels,
          aes(x = -Inf, y = Inf, label = sprintf("lost %.1f%%", pct_lost)),
          inherit.aes = FALSE,
          hjust = -0.05,
          vjust = 1.3,
          size = 6,
          size.unit = "pt"
        ) +
        facet_grid(
          sample ~ setting_label,
          axes = "all",
          axis.labels = "margins"
        ) +
        ggcyto_par_set(limits = list(x = AXIS, y = AXIS)) +
        scale_fill_gradientn(
          colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
          trans = "sqrt",
          guide = "none"
        ) +
        scale_linetype_manual(
          values = c("full recovery" = "dashed", "spike lost" = "solid"),
          name = NULL
        ) +
        coord_fixed() +
        labs(
          x = sprintf(
            "%s, unspiked unmix (arcsinh, cofactor 6000)",
            which_marker
          ),
          y = sprintf("%s, +%d arcsinh spiked unmix", which_marker, FIG_LEVEL),
          title = sprintf(
            "%s: %s spike per cell, by n.passes / n.af.passes",
            which_pop,
            which_marker
          )
        ) +
        theme(legend.position = "right")
      save_panel(
        as.ggplot(p),
        labels,
        raw_out,
        sprintf("fig_raw_%s_%s", which_marker, which_pop),
        width_mm = 183,
        height_mm = 200
      )
    })
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_raw.txt")
  )
}
