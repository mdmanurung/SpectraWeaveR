# Copy of 88d_cospike_tnf_glut1.R for the common SpectroFlo labels (plans/common_gates_plan.md G5, G10b). Only the inputs (the spike keys and scores of the 88cg family)
#          and the output dirs differ; method, levels, seed 42, loss cut and setups are inherited unchanged.
# Purpose: TNF + GLUT1 co-spike (grill 2026-11-11 Q6, Q7; README row "Co-spike (grill Q6, Q7; script
#          88d)"). The TNF +1 arcsinh cells of the first spike-in (88 spike_key_TNF.csv: CD4T and CD8T
#          cells, seed 42) get, in the same cells, TNF (BV510) +1 arcsinh plus GLUT1 (Alexa Fluor 405)
#          at +0, +1 or +2 arcsinh, each dye as its 1.8.8 spectrum added to the raw detector values. One
#          spiked set per GLUT1 level; the unspiked quarter (level 0 in the key) stays unspiked; the other
#          cells keep their raw values. Setups c188 (as188_joint), afrefine (af405_afrefine) and s102_c102
#          (as102_var, positive control) unmix each set; recovered = spiked unmix - unspiked unmix of the
#          same cell, with the unspiked unmixes spk_<setup>_base of 88 reused. Readouts per cell: TNF
#          recovery = recovered TNF / TNF added; GLUT1 recovery = recovered GLUT1 / GLUT1 added (GLUT1 > 0).
#          Unmix mode is the 88b unmix block with the set names changed.
# Run:     from the repo root, on SLURM (experiments/slurm/run_spike_cg.sbatch 88dcg_cospike_tnf_glut1.R <mode>):
#            build                 spike the raw files, one set per GLUT1 level
#            unmix <setup> <set>   setup = c188 | afrefine | s102_c102; set = cospike_g0 | cospike_g1 | cospike_g2
#            score                 recovery tables
#            figures               summary and ggcyto raw figures from the saved tables
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/{Raw, table_spectra, Reference Group Clean,
#            figure_spectral_variants}  [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/  [user zip]
#          experiments/bmv_test/compare_3D/data/lib_as102/  AutoSpectral 1.0.2 library [06_unmix_c102.R setup]
#          experiments/runs/2026-11-10_bmv-spike-in/spike_key_TNF.csv  which cells, and their level [88 build]
#          experiments/bmv_test/compare_3D/data/unmix_cells/spk_<setup>_base/  unspiked unmixes [88 unmix]
#          experiments/bmv_test/compare_3D/data/unmix_cells/af405_afrefine/af_refined.csv  refined AF set [87]
# Outputs: experiments/bmv_test/compare_3D/data/spiked/2026-11-11_cospike_g<level>/*.fcs   spiked raw files
#          experiments/bmv_test/compare_3D/data/unmix_cells/spk_<setup>_cospike_g<level>/*.fcs  unmixed files
#          <RUN_DIR>/cospike_key.csv          sample, raw_row, pop, raw_file, level, TNF and GLUT1 amounts
#          <RUN_DIR>/cospike_cells.csv        per co-spiked cell, setup and set: TNF and GLUT1 before, recovered,
#                                             recovery
#          <RUN_DIR>/cospike_summary.csv      per setup, set, pop, donor: % lost >50% TNF, median recoveries
#          <RUN_DIR>/baseline_check.csv       unspiked cells: max |spiked - base| per channel
#          <RUN_DIR>/fig_cospike_tnf_loss.{pdf,png,csv}       % TNF lost >50% by GLUT1 level, donors as points
#          <RUN_DIR>/fig_cospike_glut1_recovery.{pdf,png,csv} median GLUT1 recovery by level, donors as points
#          <RUN_DIR>/fig_cospike_raw_<pop>.{pdf,png}          ggcyto hexbin TNF before vs after, levels x setups
#          <RUN_DIR>/fig_cospike_raw_<pop>_v2.{pdf,png}       the same, two-line title (mode figures_raw)
#          <RUN_DIR>/sessionInfo_<mode>.txt (build, score, figures)

# Parameters ---------------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
SMOKE_DONOR <- Sys.getenv("SW_SMOKE_DONOR", "") # one whole donor (all its rows) for a smoke run; empty = full run
SMOKE_TAG <- if (nzchar(SMOKE_DONOR)) "smoke_" else "" # smoke runs write to their own dirs
SMOKE_SUFFIX <- if (nzchar(SMOKE_DONOR)) "-smoke" else ""
TNF_KEY <- paste0(
  "experiments/runs/2026-11-11_bmv-spike-in-cg-build",
  SMOKE_SUFFIX,
  "/spike_key_TNF.csv"
) # the cells, as 88cg chose them
RUN_DIR <- paste0(
  "experiments/runs/2026-11-11_bmv-cospike-tnf-glut1-cg",
  SMOKE_SUFFIX
)
SPIKED_PREFIX <- paste0(SMOKE_TAG, "2026-11-11_cg_") # spiked raw sets, apart from the old 2026-11-11_cospike_g*
UNMIX_PREFIX <- paste0(SMOKE_TAG, "spkcg_") # unmixed spiked sets; the base unmixes keep spk_<setup>_base
# opt188 objects and the refined AF set, as 88b (the unmix block reads both; opt188 is not a setup here)
OPT188_OBJECTS <- "experiments/runs/2026-11-10_bmv-opt188-smoke/opt188_pregate_neg_same_fsc4e5_p1"
AFREFINE_AF <- "experiments/bmv_test/compare_3D/data/unmix_cells/af405_afrefine/af_refined.csv"
TNF_DYE <- "BV510"
GLUT1_DYE <- "Alexa Fluor 405"
TNF_LEVEL <- 1 # TNF spike, arcsinh units (user, 2026-11-11, grill Q6)
SETS <- c(cospike_g0 = 0, cospike_g1 = 1, cospike_g2 = 2) # GLUT1 co-spike, arcsinh units (grill Q6)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
LOSS_FRAC <- 0.5 # "losing most of the spike" = recovery below 0.5 (64 readout)
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
  "Automated_autospectral_spectra_20261105_143803.csv"
)
read_matrix <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}

# Mode build: one spiked raw set per GLUT1 level ------------------------------------------------------
if (mode == "build") {
  suppressPackageStartupMessages({
    library(flowCore)
    library(dplyr)
    library(purrr)
    library(readr)
  })
  stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
  dir.create(RUN_DIR, recursive = TRUE)
  s188 <- read_matrix(spectra_188_file)
  # level-1 cells of the TNF key get the co-spike; level-0 cells stay unspiked as the check
  key <- read_csv(TNF_KEY, show_col_types = FALSE) |>
    filter(level %in% c(0, TNF_LEVEL)) |>
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
        write.FCS(ff, file.path(out_fcs, key_df$raw_file[1]))
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
  raw_dir <- file.path(DATA, "spiked", paste0(SPIKED_PREFIX, set))
  out_dir <- normalizePath(file.path(DATA, "unmix_cells"), mustWork = TRUE) |>
    file.path(paste0(UNMIX_PREFIX, setup, "_", set))
  stopifnot("out dir must not exist" = !dir.exists(out_dir))
  if (!setup %in% SETUPS_188) {
    .libPaths(c(file.path(DATA, "lib_as102"), .libPaths()))
  }
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
  dir.create(out_dir)
  asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
  if (setup == "opt188") {
    # opt188 exactly as 53b_opt188_unmix.R (route pregate, joint, 1 pass, 1 AF pass), on the chosen raw set
    work_dir <- file.path(out_dir, "work")
    dir.create(work_dir)
    dir_fields <- names(asp)[
      grepl("dir", names(asp)) & vapply(asp, is.character, logical(1))
    ]
    asp[dir_fields] <- lapply(asp[dir_fields], \(d) file.path(work_dir, d))
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
      "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
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
      spectra_file <- "Cells Matching Negative_autospectral_spectra.csv"
      af_file <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv"
      variants_rds <- file.path(
        control_dir,
        "figure_spectral_variants",
        "Spectral_variants.rds"
      )
    } else {
      control_dir <- file.path(d188, "Reference Group Clean")
      spectra_dir <- file.path(d188, "table_spectra")
      spectra_file <- "Automated_autospectral_spectra_20261105_143803.csv"
      af_file <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
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
            paste0("spk_", which_setup, "_base"),
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

  # unspiked cells must not move (unmixing is per cell)
  baseline_check <- cells |>
    filter(!spiked) |>
    summarise(
      max_abs_diff_tnf = max(abs(tnf_recovered)),
      max_abs_diff_glut1 = max(abs(glut1_recovered)),
      n = n(),
      .by = c(setup, set)
    )
  write_csv(baseline_check, file.path(RUN_DIR, "baseline_check.csv"))

  spiked_cells <- cells |>
    filter(spiked) |>
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
      .by = c(setup, set, glut1_level, pop, sample)
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
  .libPaths(c(file.path(DATA, "lib_ggplot2_35"), .libPaths())) # ggplot2 3.5.2 for ggcyto
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
      facet_grid(pop ~ setup_label) +
      labs(
        x = "GLUT1 (Alexa Fluor 405) added with TNF +1 arcsinh",
        y = "Cells losing >50% of the TNF spike (%)",
        shape = "Donor",
        title = "TNF + GLUT1 co-spike: TNF loss"
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
      facet_grid(pop ~ setup_label) +
      labs(
        x = "GLUT1 (Alexa Fluor 405) added with TNF +1 arcsinh",
        y = "Median GLUT1 recovery (fraction of spike)",
        shape = "Donor",
        title = "TNF + GLUT1 co-spike: GLUT1 recovery"
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

  # Figure 3 per population: ggcyto hexbin TNF before vs after, GLUT1 levels x setups
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
  for (which_pop in c("CD4T", "CD8T")) {
    d <- dplyr::filter(cells, pop == which_pop)
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
    labels <- summary_tbl |>
      dplyr::filter(pop == which_pop) |>
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
        aes(x = -Inf, y = Inf, label = sprintf("TNF lost >50%%: %.1f%%", pct)),
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
        x = "TNF, unspiked unmix (arcsinh)",
        y = "TNF, co-spiked unmix (arcsinh)",
        title = sprintf(
          "%s: TNF +1 with GLUT1 co-spike, per cell\n5 donors pooled; %% = donor median",
          which_pop
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
