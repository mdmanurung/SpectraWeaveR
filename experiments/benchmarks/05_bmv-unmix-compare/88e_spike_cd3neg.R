# Purpose: CD3- spike (grill 2026-11-10 Q8, Q9; README row "CD3- spike (grill Q8, Q9; script 88e)"). The
#          12 spiked markers (88 and 88b) at +1 arcsinh in live CD3- cells (as188_joint gates singlets/live
#          minus live/CD3pos), one marker per spiked set, 1.8.8 dye spectra, seed 42. A copy of
#          88b_spike_in_panel2.R with a new population and marker list; method unchanged. Setups: c188 =
#          as188_joint (1.8.8 zip); afrefine = af405_afrefine (87). The base unmixes spk_c188_base and
#          spk_afrefine_base from 88 are reused: they unmix the same Raw files, which spiking does not alter.
# Run:     from the repo root, on SLURM (experiments/slurm/run_af405_proof.sbatch 88e_spike_cd3neg.R <mode>):
#            build                 spike the raw files, one set per marker; writes the run dir
#            unmix <setup> <set>   setup = c188 | afrefine; set = a marker of MARKERS (base: reused from 88)
#            score                 recovery tables
#            figures               figures from the saved tables
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/{Raw, AutoSpectral_unmixed, table_spectra,
#            Reference Group Clean, figure_spectral_variants}  [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/  [user zip]
#          experiments/bmv_test/compare_3D/data/lib_as102/  AutoSpectral 1.0.2 library [06_unmix_c102.R setup]
#          experiments/runs/2026-11-10_bmv-gated/manifest.csv, gs/as188_joint  [12_cleanup.R, 19a_gate_arm.R]
# Outputs: experiments/bmv_test/compare_3D/data/spiked/2026-11-11_spike_cd3neg_<marker>/*.fcs   spiked raw files
#          experiments/bmv_test/compare_3D/data/unmix_cells/spkneg_<setup>_<marker>/*.fcs      unmixed files
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg/spike_key_<marker>.csv  sample, raw_row, pop, level, k
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg/spike_cells.csv     per spiked cell: setup, marker, added,
#                                                                       recovered, recovery
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg/spike_summary.csv   per setup, marker, pop, donor, level:
#                                                                       median recovery, share losing >50%
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg/baseline_check.csv  unspiked cells: max |spiked - base|
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg/fig_spike_<marker>.{pdf,png,csv}  level 1, donors as points
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg/fig_spike_recovery_ecdf.{pdf,png,csv}  level 1 cells, ECDF
#          experiments/runs/2026-11-11_bmv-spike-in-cd3neg/sessionInfo_<mode>.txt (build, score, figures)

# Parameters ---------------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
KEY_DIR <- "experiments/runs/2026-11-11_bmv-spike-in-cd3neg-build" # build outputs: spike keys
RUN_DIR <- "experiments/runs/2026-11-11_bmv-spike-in-cd3neg" # score and figures
SPIKED_PREFIX <- "2026-11-11_spike_cd3neg_" # spiked raw sets, apart from the T-cell sets of 88 and 88b
UNMIX_PREFIX <- "spkneg_" # unmixed spiked sets; the base unmixes keep 88's spk_<setup>_base
# opt188 objects (53b_opt188_unmix.R arm opt188_pregate_joint_p1af1) and the refined AF set (87, af405_afrefine)
OPT188_OBJECTS <- "experiments/runs/2026-11-10_bmv-opt188-smoke/opt188_pregate_neg_same_fsc4e5_p1"
AFREFINE_AF <- "experiments/bmv_test/compare_3D/data/unmix_cells/af405_afrefine/af_refined.csv"
GATED <- "experiments/runs/2026-11-10_bmv-gated" # the gates 62 used
# marker = dye (control file); the 12 markers of 88 and 88b (user, 2026-11-11)
MARKERS <- c(
  TNF = "BV510",
  Th2 = "eFluor 450",
  IL10 = "BB700",
  IFNg = "BV570",
  IL2 = "BV605",
  GLUT1 = "Alexa Fluor 405",
  ATP5A = "Alexa Fluor 488",
  CS = "Alexa Fluor 594",
  CPT1A = "PE-Cy5.5",
  ATGL = "PE-Cy7",
  G6PD = "APC-Cy7",
  CD71 = "Spark YG 593"
)
LEVELS <- c(0, 1, 2, 3) # arcsinh units added; 0 = unspiked quarter (62, user-approved grid)
FIG_LEVEL <- 1 # the level the figures show (user, 2026-11-10: +1 arcsinh)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
SEED <- 42 # which cells get which level (62, user-approved); also the figure jitter
LOSS_FRAC <- 0.5 # "losing most of the spike" = recovery below 0.5 (64 readout)
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across setups; maps raw rows to cleaned rows
POPS <- "CD3neg" # live cells outside the CD3pos gate (grill Q8)
SETUPS <- c("c188", "afrefine") # user, 2026-11-10
SETUPS_188 <- c("c188", "opt188", "afrefine") # run on the AutoSpectral 1.8.8 library
SETUP_LABELS <- c(
  c188 = "as188_joint",
  s102_c102 = "as102_var",
  s188_c102 = "s188_c102",
  opt188 = "opt188",
  afrefine = "af405_afrefine"
)

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

# Mode build: spike the raw files, one set per marker ------------------------------------------------
if (mode == "build") {
  suppressPackageStartupMessages({
    library(flowCore)
    library(flowWorkspace)
    library(dplyr)
    library(purrr)
    library(readr)
    library(tibble)
  })
  stopifnot("run dir must not exist" = !dir.exists(KEY_DIR))
  dir.create(KEY_DIR, recursive = TRUE)
  s188 <- read_matrix(spectra_188_file)
  manifest <- read_csv(
    file.path(GATED, "manifest.csv"),
    show_col_types = FALSE
  ) |>
    filter(arm == "as188_joint")
  read_x <- \(path) {
    exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
  }
  key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))
  gs <- load_gs(file.path(GATED, "gs/as188_joint"))

  # live CD3- raw rows per donor, the same for every marker
  pop_rows <- map2(set_names(manifest$sample), manifest$file, \(s, clean_file) {
    unmixed_name <- basename(clean_file)
    unmixed_key <- key_of(read_x(file.path(
      d188,
      "AutoSpectral_unmixed",
      unmixed_name
    )))
    kept_rows <- which(unmixed_key %in% key_of(read_x(clean_file)))
    imap(
      list(
        CD3neg = kept_rows[
          gh_pop_get_indices(gs[[s]], "live") &
            !gh_pop_get_indices(gs[[s]], "CD3pos")
        ]
      ),
      \(r, p) {
        tibble(raw_row = r, pop = p)
      }
    ) |>
      list_rbind() |>
      distinct(raw_row, .keep_all = TRUE) |>
      mutate(
        raw_file = sub(" AutoSpectral\\.fcs$", ".fcs", unmixed_name),
        unmixed_key = list(unmixed_key)
      )
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
        "raw rows equal unmixed rows" = identical(
          key_of(exprs(ff)),
          rows$unmixed_key[[1]]
        )
      )
      key <- rows |>
        select(raw_row, pop, raw_file) |>
        mutate(
          level = sample(rep_len(LEVELS, n())),
          k = COFACTOR * sinh(level),
          sample = s
        )
      x <- exprs(ff)
      x[key$raw_row, colnames(s188)] <- x[key$raw_row, colnames(s188)] +
        outer(key$k, s188[dye, ])
      exprs(ff) <- x
      write.FCS(ff, file.path(out_fcs, rows$raw_file[1]))
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
      set %in% names(MARKERS)
  )
  raw_dir <- if (set == "base") {
    file.path(d188, "Raw")
  } else {
    file.path(DATA, "spiked", paste0(SPIKED_PREFIX, set))
  }
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
  stopifnot("score run dir must not exist" = !dir.exists(RUN_DIR))
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
            paste0("spk_", which_setup, "_base"),
            f
          ))
          stopifnot(
            "spiked and base unmix rows match" = identical(
              x_spk[, KEY_COLS],
              x_base[, KEY_COLS]
            )
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

  # Unspiked quarter: the spike must not move them (unmixing is per cell)
  baseline_check <- cells |>
    filter(level == 0) |>
    summarise(
      max_abs_diff = max(abs(recovered)),
      n = n(),
      .by = c(setup, marker)
    )
  write_csv(baseline_check, file.path(RUN_DIR, "baseline_check.csv"))

  spiked <- cells |>
    filter(level > 0) |>
    mutate(recovery = recovered / k)
  write_csv(spiked, file.path(RUN_DIR, "spike_cells.csv"))
  summary_tbl <- spiked |>
    summarise(
      n_cells = n(),
      median_recovery = median(recovery),
      share_lost = mean(recovery < LOSS_FRAC),
      .by = c(setup, marker, pop, sample, level)
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
    c("solid", "dashed")[seq_along(setup_pal)],
    names(setup_pal)
  ) # second encoding of setup
  donor_shapes <- c(A3 = 16, A6 = 17, B11 = 15, B2 = 18, D1 = 8)
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
    "every donor has a shape" = all(fig_frame$sample %in% names(donor_shapes))
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
        aes(colour = setup, shape = sample),
        position = position_jitter(width = 0.15, height = 0, seed = SEED),
        size = 1.5
      ) +
      facet_grid(readout ~ pop, scales = "free_y", switch = "y") +
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
    facet_grid(marker ~ pop, switch = "y") +
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
