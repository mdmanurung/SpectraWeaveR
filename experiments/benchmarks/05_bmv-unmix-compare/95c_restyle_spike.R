# Purpose: redraw the certified spike summary figures (report sections 10, 18, 20, 22, 23 and 25) from the
#          common-gate spike runs (plans/common_gates_plan.md G5, G10b, G11, G13). Copy of 95_restyle_spike.R
#          with the sources repointed to the *-cg runs, setups named through arm_label(), the donor colour
#          and shape of report_figures.R, and no lines across setups. Nothing is unmixed or rescored.
#          Spike level +1 arcsinh, loss cut recovery < 0.5 and seed 42 are inherited (G5).
#          The markers are the keys of spike_marker_dyes in the dataset config; a figure whose marker is not in
#          that list stops the script (EXV has no Th2, IL-2 or IL-10); REDRAW_ONLY names the figures to draw.
# Run:     SW_DATASET=3d SW_RUN_DATE=YYYY-MM-DD REDRAW_ONLY=spike pixi run Rscript
#            experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_spike.R
#            from the repo root. REDRAW_ONLY=spike (or a list of figure stems) adds only these figures to
#            MANIFEST_restyle_common.csv and reuses the existing report_figures and 10_af405_proof folders.
# Inputs:  run_in("spike-in-cg")/{spike_summary,spike_cells}.csv         TNF, Th2, IL10; 5 setups
#          run_in("spike-in-panel2-cg")/{spike_summary,spike_cells}.csv   9 markers; 2 setups
#          run_in("spike-in-cd3neg-cg")/{spike_summary,spike_cells}.csv   12 markers in CD3-
#          run_in("cospike-tnf-glut1-cg")/cospike_summary.csv             TNF +1, GLUT1 +0/+1/+2
#          run_in("tnf-spike-by-glut1-cg")/summary_by_glut1.csv           TNF +1 by GLUT1 quartile
#          run_in("spike-passes-cg")/spike_summary.csv                    n.passes / n.af.passes
#          run_in("af405-proof-cg")/{spike_summary,fig5b_..._source,fig6_..._source}.csv
#          datasets/<SW_DATASET>.json   targets, spike_marker_dyes, conditions
# Outputs: run_out("restyle-common") (shared by the 95c scripts)/18_spike_in/fig_spike_<marker>.{pdf,png,csv}, fig_spike_recovery_ecdf
#          .../20_spike_in_panel2/, .../22_spike_cd3neg/   same names
#          .../23_cospike/fig_cospike_tnf_loss, fig_cospike_glut1_recovery, by_glut1_fig_loss_by_glut1
#          .../25_spike_passes/fig_passes_<marker>
#          .../10_af405_proof/fig5_spike_recovery, fig5b_recovery_by_glut1, fig6_anatomy
#          .../report_figures/s10_spike, plus its row in captions.csv
#          .../manifest_parts/95c_restyle_spike.R.csv   the ECDF .csv holds the curve on ECDF_GRID, not the cells

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
source("experiments/benchmarks/05_bmv-unmix-compare/95c_spike_common.R")
SCRIPT <- "95c_restyle_spike.R"
TNF_LAB <- cfg("targets")$tnf_label
IFNG_LAB <- cfg("targets")$ifng_label
TNF_DYE <- cfg("targets")$tnf_dye
GLUT1_DYE <- cfg("targets")$glut1_dye

# Parameters ---------------------------------------------------------------------------------------
ECDF_GRID <- seq(-0.5, 1.5, by = 0.0025) # display window of the ECDF; cells outside it are counted in the .csv
RECOVERY_FLOOR <- c(0.99, 1.01) # recovery axes span at least this, so panels differing by 1e-7 keep readable ticks
LOST_FLOOR <- c(0, 5) # % lost axes span at least this; display only, no point moves
SETUP_LINES <- setNames(
  c("solid", "dashed", "dotted", "dotdash", "longdash"),
  setup_label(c("c188", "s102_c102", "s188_c102", "opt188", "afrefine"))
) # second encoding of setup in the ECDFs
PASS_MARKERS <- c(TNF_LAB, "Th2", "IL10") # markers of the n.passes figures (section 25)

# the figures report sections 18, 20, 22 and 23 cite (90_report.qmd show_kept), plus each section's ECDF
SPIKE_FIGURES <- tribble(
  ~section             , ~run                 , ~marker  ,
  "18_spike_in"        , "spike-in-cg"        , TNF_LAB  ,
  "18_spike_in"        , "spike-in-cg"        , "Th2"    ,
  "18_spike_in"        , "spike-in-cg"        , "IL10"   ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" , "GLUT1"  ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" , "CS"     ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" , IFNG_LAB ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" , "IL2"    ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" , "ATP5A"  ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" , "CPT1A"  ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" , "ATGL"   ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" , "G6PD"   ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" , "CD71"   ,
  "22_spike_cd3neg"    , "spike-in-cd3neg-cg" , "GLUT1"  ,
  "22_spike_cd3neg"    , "spike-in-cd3neg-cg" , TNF_LAB
)
# every figure to draw needs its marker in the dataset's spike panel (stop, never skip or reuse a 3D value)
drawn_markers <- c(
  SPIKE_FIGURES$marker[map_lgl(
    paste0("fig_spike_", SPIKE_FIGURES$marker),
    want
  )],
  PASS_MARKERS[map_lgl(paste0("fig_passes_", PASS_MARKERS), want)]
)
absent <- setdiff(drawn_markers, names(MARKERS))
if (length(absent) > 0) {
  stop(
    SCRIPT,
    ": needs ",
    paste(absent, collapse = ", "),
    ", missing in dataset ",
    DS$dataset,
    " (name other figures in REDRAW_ONLY)",
    call. = FALSE
  )
}
ECDF_FIGURES <- tribble(
  ~section             , ~run                 ,
  "18_spike_in"        , "spike-in-cg"        ,
  "20_spike_in_panel2" , "spike-in-panel2-cg" ,
  "22_spike_cd3neg"    , "spike-in-cd3neg-cg"
)
# old path of each figure, as audit_freshness.csv writes it (under experiments/runs/)
OLD_SPIKE_DIR <- c(
  `18_spike_in` = "2026-11-11_bmv-restyle/18_spike_in",
  `20_spike_in_panel2` = "2026-11-11_bmv-restyle/20_spike_in_panel2",
  `22_spike_cd3neg` = "2026-11-11_bmv-restyle/22_spike_cd3neg"
)

# Helpers used by every figure ---------------------------------------------------------------------
# Tick labels with the fewest significant digits that keep the ticks of one panel distinct.
tick_labels <- \(x) {
  shown <- x[!is.na(x)]
  digits <- 2
  while (digits < 9 && anyDuplicated(signif(shown, digits))) {
    digits <- digits + 1
  }
  as.character(signif(x, digits))
}
donor_factor <- \(x) factor_by_pal(x, DONOR_PAL)

# Setups down the y axis, donors as points over the median bar, two readouts across, populations down.
draw_summary <- function(
  d,
  group,
  group_levels,
  title,
  y_title,
  out_dir,
  stem
) {
  panel_data <- d |>
    sample_keys() |>
    need_cond(SCRIPT) |>
    mutate(
      group = factor(.data[[group]], levels = group_levels),
      donor = donor_factor(.data$donor_key),
      donor_key = NULL,
      pop = factor(.data$pop, levels = c("CD4T", "CD8T", "CD3neg")),
      `Cells losing >50% of spike (%)` = 100 * .data$share_lost,
      `Median recovery (fraction of spike)` = .data$median_recovery
    ) |>
    pivot_longer(
      c(
        `Cells losing >50% of spike (%)`,
        `Median recovery (fraction of spike)`
      ),
      names_to = "readout",
      values_to = "value"
    ) |>
    mutate(group = droplevels(.data$group), pop = droplevels(.data$pop))
  n_groups <- nlevels(panel_data$group)
  n_pops <- nlevels(panel_data$pop) * max(1, length(DS$conditions))
  floors <- tibble(
    readout = c(
      rep("Cells losing >50% of spike (%)", 2),
      rep("Median recovery (fraction of spike)", 2)
    ),
    value = c(LOST_FLOOR, RECOVERY_FLOOR)
  )
  p <- ggplot(panel_data, aes(value, group)) +
    geom_blank(data = floors, aes(x = value), inherit.aes = FALSE) +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      colour = "grey30",
      width = 0.5,
      linewidth = FIG_LINE / 2,
      orientation = "y"
    ) +
    geom_point(
      aes(colour = donor, shape = donor),
      position = position_jitter(width = 0, height = 0.15, seed = SEED),
      size = 1.2
    ) +
    facet_grid(
      rows = vars(pop, !!!rlang::syms(COND)),
      cols = vars(readout),
      scales = "free_x",
      axes = "all",
      axis.labels = "margins",
      labeller = labeller(readout = label_wrap_gen(28))
    ) +
    donor_scales() +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(
      breaks = scales::breaks_extended(n = 4),
      labels = tick_labels,
      guide = guide_axis(check.overlap = TRUE)
    ) +
    labs(x = NULL, y = y_title, title = title) +
    theme(panel.spacing.x = unit(4, "mm")) +
    value_grid("x")
  height <- max(55, 22 + n_pops * (6 + 7 * n_groups))
  save_panel(p, panel_data, out_dir, stem, 183, height)
}

# Ordered x (spike level, quartile, decile): donors as points, median bar per x and panel, no lines.
plot_donor_levels <- function(
  d,
  x,
  y,
  facets,
  x_title,
  y_title,
  title,
  y_floor = NULL,
  hline = NULL
) {
  p <- ggplot(d, aes(.data[[x]], .data[[y]]))
  if (!is.null(hline)) {
    p <- p +
      geom_hline(
        yintercept = hline,
        linetype = "dashed",
        colour = "grey50",
        linewidth = FIG_LINE / 2
      )
  }
  if (!is.null(y_floor)) {
    p <- p +
      geom_blank(data = tibble(y = y_floor), aes(y = y), inherit.aes = FALSE)
  }
  p +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      colour = "grey30",
      width = 0.6,
      linewidth = FIG_LINE / 2
    ) +
    geom_point(
      aes(colour = donor, shape = donor),
      position = position_jitter(width = 0.12, height = 0, seed = SEED),
      size = 1.2
    ) +
    facets +
    donor_scales() +
    scale_y_continuous(labels = tick_labels) +
    labs(x = x_title, y = y_title, title = title) +
    value_grid("y")
}
FACET_ALL <- \(...) facet_grid(..., axes = "all", axis.labels = "margins")

# Per-marker summary figures: sections 18, 20, 22 ------------------------------------------------------
pwalk(SPIKE_FIGURES, \(section, run, marker) {
  stem <- paste0("fig_spike_", marker)
  if (!want(stem)) {
    return(invisible())
  }
  out_dir <- fresh_dir(section)
  src <- SRC(run, "spike_summary.csv")
  d <- read_csv(src, show_col_types = FALSE) |>
    filter(.data$level == FIG_LEVEL, .data$marker == .env$marker) |>
    mutate(setup = setup_label(.data$setup))
  levels_shown <- intersect(SETUP_ORDER, unique(d$setup))
  draw_summary(
    d,
    "setup",
    levels_shown,
    sprintf(
      "%s (%s) spike +%d arcsinh",
      marker,
      MARKERS[[marker]],
      FIG_LEVEL
    ),
    NULL,
    out_dir,
    stem
  )
  add_manifest(
    SCRIPT,
    file.path(OLD_SPIKE_DIR[[section]], paste0(stem, ".png")),
    file.path(out_dir, paste0(stem, ".png")),
    src
  )
})

# Recovery ECDFs: every level-1 cell pooled over donors, one curve per setup -----------------------------
pwalk(ECDF_FIGURES, \(section, run) {
  stem <- "fig_spike_recovery_ecdf"
  key <- paste0(stem, "_", section)
  if (!want(stem) && !want(key)) {
    return(invisible())
  }
  out_dir <- fresh_dir(section)
  src <- SRC(run, "spike_cells.csv")
  spiked <- data.table::fread(
    src,
    select = c("setup", "marker", "pop", "level", "recovery", COND)
  ) |>
    as_tibble() |>
    filter(.data$level == FIG_LEVEL)
  setup_pal <- categorical_palette(SETUP_ORDER)
  # ECDF per setup, marker and population on the grid; cells outside the window are counted, not dropped
  ecdf_frame <- spiked |>
    reframe(
      n_cells = n(),
      n_below_window = sum(.data$recovery < min(ECDF_GRID)),
      n_above_window = sum(.data$recovery > max(ECDF_GRID)),
      cumulative = ecdf(.data$recovery)(ECDF_GRID),
      recovery_grid = ECDF_GRID,
      .by = c(setup, marker, pop, any_of(COND))
    ) |>
    mutate(
      setup = factor(setup_label(.data$setup), levels = SETUP_ORDER),
      marker = factor(
        .data$marker,
        levels = intersect(names(MARKERS), unique(.data$marker))
      )
    ) |>
    mutate(setup = droplevels(.data$setup))
  n_pops <- n_distinct(ecdf_frame$pop)
  n_rows <- if (n_pops == 1) {
    ceiling(nlevels(ecdf_frame$marker) / 4)
  } else {
    nlevels(ecdf_frame$marker)
  }
  facets <- if (n_pops == 1) {
    facet_wrap(
      vars(marker, !!!rlang::syms(COND)),
      ncol = 4,
      axes = "all",
      axis.labels = "margins"
    )
  } else {
    facet_grid(
      rows = vars(marker),
      cols = vars(pop, !!!rlang::syms(COND)),
      axes = "all",
      axis.labels = "margins"
    )
  }
  p <- ggplot(
    ecdf_frame,
    aes(recovery_grid, cumulative, colour = setup, linetype = setup)
  ) +
    geom_vline(
      xintercept = c(LOSS_FRAC, 1),
      linetype = "dashed",
      colour = "grey60",
      linewidth = FIG_LINE / 2
    ) +
    geom_step(linewidth = FIG_LINE) +
    facets +
    scale_colour_manual(values = setup_pal, name = "Setup") +
    scale_linetype_manual(values = SETUP_LINES, name = "Setup") +
    scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5)) +
    labs(
      x = "Recovered fraction of the +1 spike (window -0.5 to 1.5)",
      y = "Cumulative fraction of spiked cells",
      title = sprintf(
        "Spike recovery per cell (%s), donors pooled",
        paste(unique(ecdf_frame$pop), collapse = ", ")
      )
    ) +
    value_grid("x")
  save_panel(p, ecdf_frame, out_dir, stem, 183, 24 + n_rows * 22)
  add_manifest(
    SCRIPT,
    file.path(OLD_SPIKE_DIR[[section]], paste0(stem, ".png")),
    file.path(out_dir, paste0(stem, ".png")),
    src,
    note = "curve on ECDF_GRID; cells below or above the window are counted in the csv"
  )
})

# Co-spike and TNF loss by natural GLUT1: section 23 ---------------------------------------------------
cospike_src <- SRC("cospike-tnf-glut1-cg", "cospike_summary.csv")
by_glut1_src <- SRC("tnf-spike-by-glut1-cg", "summary_by_glut1.csv")
COSPIKE_X <- paste0(
  "GLUT1 (",
  GLUT1_DYE,
  ") added with ",
  TNF_LAB,
  " +1 arcsinh (arcsinh units)"
)
# rows of the pop-by-setup grids: population, then condition when the dataset has them
POP_BY_SETUP <- \(...) {
  FACET_ALL(
    rows = vars(pop, !!!rlang::syms(COND)),
    cols = vars(setup_label),
    ...
  )
}
if (any(want(c("fig_cospike_tnf_loss", "fig_cospike_glut1_recovery")))) {
  out_dir <- fresh_dir("23_cospike")
  cospike <- read_csv(cospike_src, show_col_types = FALSE) |>
    sample_keys() |>
    need_cond(SCRIPT) |>
    mutate(
      setup_label = factor(
        setup_label(.data$setup),
        levels = intersect(SETUP_ORDER, setup_label(unique(.data$setup)))
      ),
      glut1_label = factor(
        paste0("+", .data$glut1_level),
        levels = paste0("+", 0:2)
      ),
      donor = donor_factor(.data$donor_key),
      donor_key = NULL
    )
  facets <- POP_BY_SETUP(labeller = labeller(setup_label = label_wrap_gen(22)))
  p_tnf <- plot_donor_levels(
    cospike,
    "glut1_label",
    "pct_tnf_lost",
    facets,
    COSPIKE_X,
    paste0("Cells losing >50% of the ", TNF_LAB, " spike (%)"),
    paste0(TNF_LAB, " + GLUT1 co-spike: ", TNF_LAB, " loss"),
    y_floor = LOST_FLOOR
  )
  save_panel(p_tnf, cospike, out_dir, "fig_cospike_tnf_loss", 183, 85)
  glut1_df <- filter(cospike, .data$glut1_level > 0)
  p_glut1 <- plot_donor_levels(
    glut1_df,
    "glut1_label",
    "median_glut1_recovery",
    facets,
    COSPIKE_X,
    "Median GLUT1 recovery (fraction of spike)",
    paste0(TNF_LAB, " + GLUT1 co-spike: GLUT1 recovery"),
    y_floor = RECOVERY_FLOOR,
    hline = 1
  )
  save_panel(p_glut1, glut1_df, out_dir, "fig_cospike_glut1_recovery", 183, 85)
  walk(c("fig_cospike_tnf_loss", "fig_cospike_glut1_recovery"), \(stem) {
    add_manifest(
      SCRIPT,
      file.path("2026-11-11_bmv-restyle/23_cospike", paste0(stem, ".png")),
      file.path(out_dir, paste0(stem, ".png")),
      cospike_src
    )
  })
}
if (want("by_glut1_fig_loss_by_glut1")) {
  out_dir <- fresh_dir("23_cospike")
  by_glut1 <- read_csv(by_glut1_src, show_col_types = FALSE) |>
    sample_keys() |>
    need_cond(SCRIPT) |>
    mutate(
      setup_label = factor(
        setup_label(.data$setup),
        levels = intersect(SETUP_ORDER, setup_label(unique(.data$setup)))
      ),
      glut1_bin = factor(
        sub("GLUT1 ", "", .data$glut1_bin),
        levels = paste0("Q", 1:4)
      ),
      donor = donor_factor(.data$donor_key),
      donor_key = NULL
    )
  p_loss <- plot_donor_levels(
    by_glut1,
    "glut1_bin",
    "pct_lost",
    POP_BY_SETUP(labeller = labeller(setup_label = label_wrap_gen(22))),
    paste0("Natural GLUT1 (", GLUT1_DYE, "), quartile per donor"),
    paste0("Cells losing >50% of the ", TNF_LAB, " spike (%)"),
    paste0(TNF_LAB, " +1 arcsinh spike by natural GLUT1"),
    y_floor = LOST_FLOOR
  )
  save_panel(p_loss, by_glut1, out_dir, "by_glut1_fig_loss_by_glut1", 183, 90)
  add_manifest(
    SCRIPT,
    "2026-11-11_bmv-restyle/23_cospike/by_glut1_fig_loss_by_glut1.png",
    file.path(out_dir, "by_glut1_fig_loss_by_glut1.png"),
    by_glut1_src
  )
}

# Pass settings: section 25 ----------------------------------------------------------------------------
pass_src <- SRC("spike-passes-cg", "spike_summary.csv")
walk(PASS_MARKERS, \(marker) {
  stem <- paste0("fig_passes_", marker)
  if (!want(stem)) {
    return(invisible())
  }
  out_dir <- fresh_dir("25_spike_passes")
  d <- read_csv(pass_src, show_col_types = FALSE) |>
    filter(.data$level == FIG_LEVEL, .data$marker == .env$marker) |>
    mutate(setting = unname(SETTING_ORDER[.data$setting]))
  draw_summary(
    d,
    "setting",
    unname(SETTING_ORDER),
    sprintf("%s (%s) spike +%d arcsinh", marker, MARKERS[[marker]], FIG_LEVEL),
    "n.passes / n.af.passes",
    out_dir,
    stem
  )
  add_manifest(
    SCRIPT,
    file.path("2026-11-11_bmv-spike-passes", paste0(stem, ".png")),
    file.path(out_dir, paste0(stem, ".png")),
    pass_src
  )
})

# AF405 proof: spike recovery, GLUT1 deciles, anatomy -----------------------------------------------------
OLD_PROOF <- "2026-11-10_bmv-af405-proof"
proof_summary <- SRC("af405-proof-cg", "spike_summary.csv")
PROOF_SETUPS <- c("1.8.8", "1.8.8 + 1.0.2 AF405 variants")
proof_levels <- \(x) {
  factor(paste0("+", x, " arcsinh"), levels = paste0("+", 1:3, " arcsinh"))
}
proof <- read_csv(proof_summary, show_col_types = FALSE) |>
  sample_keys() |>
  need_cond(SCRIPT) |>
  mutate(
    setup = factor(.data$setup, levels = PROOF_SETUPS),
    level_label = proof_levels(.data$level),
    donor = donor_factor(.data$donor_key),
    donor_key = NULL
  )
# the AF405 proof spikes the BV510 slot (TNF in 3D)
PROOF_TITLE <- paste0(
  "AF405 proof: recovery of a known ",
  TNF_DYE,
  " (",
  TNF_LAB,
  ") spike"
)

# s10_spike: the report figure, level across, setup in the columns
if (want("s10_spike")) {
  out_dir <- fresh_dir("report_figures")
  p10 <- plot_donor_levels(
    proof,
    "level_label",
    "median_recovery",
    FACET_ALL(rows = vars(pop, !!!rlang::syms(COND)), cols = vars(setup)),
    "Spike level",
    paste0("Median spike recovery (", TNF_DYE, ")"),
    PROOF_TITLE,
    y_floor = RECOVERY_FLOOR,
    hline = 1
  )
  save_panel(p10, proof, out_dir, "s10_spike", 183, 85)
  caption <- tibble(
    key = "s10_spike",
    caption = paste0(
      "Median recovery of the added ",
      TNF_DYE,
      " (",
      TNF_LAB,
      ") spike per spike level, population (rows) and setup (columns). ",
      PAIRED
    )
  )
  captions_file <- file.path(out_dir, "captions.csv")
  captions <- read_csv(captions_file, show_col_types = FALSE) |>
    filter(.data$key != "s10_spike") |>
    bind_rows(caption)
  write_csv(captions, captions_file)
  add_manifest(
    SCRIPT,
    "2026-11-10_bmv-report/figures-v7/s10_spike.png",
    file.path(out_dir, "s10_spike.png"),
    proof_summary
  )
}

# fig5: setups down the y axis, spike level across
if (want("fig5_spike_recovery")) {
  out_dir <- fresh_dir("10_af405_proof")
  p5 <- ggplot(proof, aes(median_recovery, setup)) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "grey50",
      linewidth = FIG_LINE / 2
    ) +
    geom_blank(
      data = tibble(median_recovery = RECOVERY_FLOOR),
      inherit.aes = FALSE,
      aes(x = median_recovery)
    ) +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      colour = "grey30",
      width = 0.5,
      linewidth = FIG_LINE / 2,
      orientation = "y"
    ) +
    geom_point(
      aes(colour = donor, shape = donor),
      position = position_jitter(width = 0, height = 0.15, seed = SEED),
      size = 1.2
    ) +
    FACET_ALL(
      rows = vars(pop, !!!rlang::syms(COND)),
      cols = vars(level_label),
      scales = "free_x"
    ) +
    donor_scales() +
    scale_y_discrete(limits = rev, labels = scales::label_wrap(22)) +
    scale_x_continuous(
      breaks = scales::breaks_extended(n = 3),
      labels = tick_labels,
      guide = guide_axis(check.overlap = TRUE)
    ) +
    labs(
      x = paste0("Median share of added ", TNF_DYE, " recovered"),
      y = NULL,
      title = PROOF_TITLE
    ) +
    theme(panel.spacing.x = unit(4, "mm")) +
    value_grid("x")
  save_panel(p5, proof, out_dir, "fig5_spike_recovery", 183, 80)
  add_manifest(
    SCRIPT,
    file.path(OLD_PROOF, "fig5_spike_recovery.png"),
    file.path(out_dir, "fig5_spike_recovery.png"),
    proof_summary
  )
}

# fig5b: median recovery per GLUT1 decile within donor
fig5b_src <- SRC("af405-proof-cg", "fig5b_recovery_by_glut1_source.csv")
if (want("fig5b_recovery_by_glut1")) {
  out_dir <- fresh_dir("10_af405_proof")
  fig5b <- read_csv(fig5b_src, show_col_types = FALSE) |>
    sample_keys() |>
    need_cond(SCRIPT) |>
    mutate(
      setup = factor(.data$setup, levels = PROOF_SETUPS),
      level_label = proof_levels(.data$level),
      decile = factor(.data$glut1_decile, levels = 1:10),
      donor = donor_factor(.data$donor_key),
      donor_key = NULL
    )
  p5b <- plot_donor_levels(
    fig5b,
    "decile",
    "median_recovery",
    FACET_ALL(
      rows = vars(level_label, !!!rlang::syms(COND)),
      cols = vars(pop, setup),
      labeller = labeller(setup = label_wrap_gen(16))
    ),
    paste0("GLUT1 (", GLUT1_DYE, ") decile within donor"),
    "Median share recovered",
    paste0("AF405 proof: ", TNF_DYE, " recovery by GLUT1 decile"),
    y_floor = RECOVERY_FLOOR,
    hline = 1
  )
  save_panel(p5b, fig5b, out_dir, "fig5b_recovery_by_glut1", 183, 110)
  add_manifest(
    SCRIPT,
    file.path(OLD_PROOF, "fig5b_recovery_by_glut1.png"),
    file.path(out_dir, "fig5b_recovery_by_glut1.png"),
    fig5b_src
  )
}

# fig6: raw and fitted spectra of 10 squeezed CD8T cells (detectors are an ordered axis, so lines stay)
fig6_src <- SRC("af405-proof-cg", "fig6_anatomy_source.csv")
if (want("fig6_anatomy")) {
  out_dir <- fresh_dir("10_af405_proof")
  fig6 <- read_csv(fig6_src, show_col_types = FALSE)
  detectors <- unique(fig6$detector)
  SERIES <- c("raw", "fit 1.8.8", "fit 1.0.2 AF405 variants")
  fig6 <- fig6 |>
    mutate(
      detector = factor(.data$detector, levels = detectors),
      series = factor(.data$series, levels = SERIES)
    )
  p6 <- ggplot(fig6, aes(detector, value, colour = series, group = series)) +
    geom_line(linewidth = FIG_LINE) +
    facet_wrap(
      ~cell,
      ncol = 2,
      scales = "free_y",
      axes = "all",
      axis.labels = "margins"
    ) +
    scale_colour_manual(values = categorical_palette(SERIES), name = NULL) +
    scale_x_discrete(breaks = detectors[seq(1, length(detectors), by = 8)]) +
    labs(
      x = "Detector (every 8th labelled)",
      y = "Raw signal",
      title = "AF405 proof: 10 CD8T cells, raw spectrum against the two fits"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  save_panel(p6, fig6, out_dir, "fig6_anatomy", 183, 170)
  add_manifest(
    SCRIPT,
    file.path(OLD_PROOF, "fig6_anatomy.png"),
    file.path(out_dir, "fig6_anatomy.png"),
    fig6_src
  )
}
