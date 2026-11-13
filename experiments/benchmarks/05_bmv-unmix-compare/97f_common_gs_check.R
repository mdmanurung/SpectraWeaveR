# Purpose: sanity check of the common GatingSets (97e): CD4 vs CD8 of the PLOT_DONOR sample(s) in the P arms
#          that 97e gated, every cell coloured by its SpectroFlo label (G11). Gate lines appear on the SpectroFlo
#          panel only. Also writes the % of CD8T above the TNF cut per P arm and sample, with the median per
#          arm. The TNF table needs a dataset that has TNF; it stops for EXV (BV510 is CD1c there).
# Run:     SW_DATASET=3d SW_RUN_DATE=2026-11-11 pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/97f_common_gs_check.R
#          from the repo root, after 97e_common_gs.R. SW_DATASET is 3d, 24h or exv and SW_RUN_DATE is
#          YYYY-MM-DD (neither has a default).
# Inputs:  registry keys of datasets/<SW_DATASET>_runs.csv:
#          run_in("common-gated")/{gs,gate_counts}/<arm>, manifest.csv   [97e]
#          run_in("common-gates-clean3")/cuts.csv   SpectroFlo gate bounds [97d]
#          datasets/<SW_DATASET>_samples.csv, <SW_DATASET>.json   samples, channels, targets, lib_dirs
#          experiments/benchmarks/05_bmv-unmix-compare/{arm_labels.R,lib/theme.R}
# Outputs: <out_dir> = run_out("common-gated-check")
#          <out_dir>/cd4_cd8_<sample>_p_arms.{png,pdf,csv}   the figure and its source data (events drawn), one
#                                                           per sample of PLOT_DONOR (3D: A3)
#          <out_dir>/cd8t_tnf_pos_p_arms.csv           arm, donor (= sample), n, n_parent, pct_pos, condition (NA
#                                                      without conditions), median_pct per arm and condition
#          <out_dir>/sessionInfo_97f.txt

# dataset.R first: it only loads dplyr and readr, and cfg("lib_dirs") must be known before the other libraries
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths())) # ggplot2 3.5.2 for ggcyto
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS, VERSION_ROWS
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set, save_panel, FIG_LINE

# Parameters ---------------------------------------------------------------------------------------
GATED <- run_in("common-gated")
CUTS <- file.path(run_in(cfg("label_key")), "cuts.csv") # 3D: 97d cuts per sample; 24H and ex vivo: 99a pooled cuts
OUT_DIR <- run_out("common-gated-check") # new dated dir; must not exist
PLOT_DONOR <- "BMV_009C" # donor drawn: sample A3 in 3D and EXV; 24H has one sample per condition, each gets a figure
P_ARMS <- unlist(VERSION_ROWS, use.names = FALSE) # the 11 P arms, in version-group order
CD4 <- cfg("channels")$CD4 # channels as in 19a_gate_arm.R
CD8 <- cfg("channels")$CD8
TNF_LABEL <- cfg("targets")$tnf_label # gate name TNFpos; the CD8T table needs TNF itself
AXIS_T <- c(-3, 6) # display window, arcsinh; labels and counts use every CD3+ cell
N_DRAW <- 15000 # cells drawn per panel; the same rows in every arm, so panels show the same cells
SEED <- 42 # for the row subsample
QUAD_FILL <- c(
  CD4T = "#2166ac",
  CD8T = "#b2182b",
  DP = "#762a83",
  DN = "#7f7f7f"
) # as 97d
T_PARENT <- "/nonDebris/singlets/live/CD3pos"

stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE)

# Samples and arms ---------------------------------------------------------------------------------
smp <- samples_tbl()
plot_samples <- smp$sample[smp$donor == PLOT_DONOR]
if (length(plot_samples) == 0) {
  stop(
    "97f_common_gs_check.R: donor ",
    PLOT_DONOR,
    " is not in the sample table of dataset ",
    DS$dataset,
    call. = FALSE
  )
}
# the P arms that 97e gated (a dataset without the self-unmixed arms yet has the four zip arms only)
ARMS <- intersect(
  P_ARMS,
  unique(read_csv(file.path(GATED, "manifest.csv"), show_col_types = FALSE)$arm)
)
if (length(ARMS) < length(P_ARMS)) {
  message(
    "97f_common_gs_check.R: P arms not in the 97e manifest, not drawn: ",
    paste(setdiff(P_ARMS, ARMS), collapse = ", ")
  )
}
stopifnot("spectroflo is among the arms" = "spectroflo" %in% ARMS)
cuts_97d <- read_csv(CUTS, show_col_types = FALSE)
if (!"version" %in% names(cuts_97d)) {
  # 99a cuts.csv: one pooled cut per node; polygon bounds are the quadrant edges. Give every sample the same cut
  cuts_97d <- cuts_97d |>
    dplyr::filter(type == "polygon", node %in% names(QUAD_FILL)) |>
    transmute(
      version = "clean3_pooled",
      gate = node,
      dim,
      lo = min,
      hi = max
    ) |>
    tidyr::crossing(sample = samples_tbl()$sample)
}

# One figure per sample: CD3+ cells with their SpectroFlo quadrant label ----------------------------
for (s in plot_samples) {
  events <- map(set_names(ARMS), \(a) {
    gh <- load_gs(file.path(GATED, "gs", a))[[s]]
    root <- exprs(gh_pop_get_data(gh, "root"))
    in_cd3 <- gh_pop_get_indices(gh, T_PARENT)
    quad <- map(names(QUAD_FILL), \(q) {
      gh_pop_get_indices(gh, paste0(T_PARENT, "/", q))
    }) |>
      set_names(names(QUAD_FILL))
    quad_mat <- do.call(cbind, quad)
    stopifnot(
      "every CD3+ cell is in exactly one quadrant" = all(
        rowSums(quad_mat)[in_cd3] == 1
      )
    )
    label <- names(QUAD_FILL)[max.col(quad_mat, ties.method = "first")]
    tibble(
      row = which(in_cd3),
      cd4 = root[in_cd3, CD4],
      cd8 = root[in_cd3, CD8],
      label = label[in_cd3]
    )
  }) |>
    list_rbind(names_to = "arm")
  # every CD3+ cell sits in one quadrant, and the rows are the same cells in every arm
  stopifnot(
    "rows are the same in every arm" = n_distinct(map(
      split(events$row, events$arm),
      identity
    )) ==
      1
  )
  set.seed(SEED)
  rows_drawn <- sample(unique(events$row), N_DRAW)
  plot_data <- events |>
    dplyr::filter(row %in% rows_drawn) |>
    mutate(
      arm_label = factor(
        str_wrap(ARM_LABELS[arm], 24),
        levels = str_wrap(ARM_LABELS[ARMS], 24)
      ),
      label = factor(label, levels = names(QUAD_FILL))
    ) |>
    arrange(label %in% c("DP", "DN")) # rare populations drawn last so they stay visible

  # Gate lines: SpectroFlo panel only, from the pooled clean3 bounds of this sample's gate group
  gate_rects <- cuts_97d |>
    dplyr::filter(
      version == "clean3_pooled",
      sample == s,
      gate %in% names(QUAD_FILL),
      dim %in% c(CD4, CD8)
    ) |>
    select(gate, dim, lo, hi) |>
    tidyr::pivot_wider(names_from = dim, values_from = c(lo, hi)) |>
    transmute(
      quadrant = gate,
      xmin = pmax(.data[[paste0("lo_", CD4)]], AXIS_T[[1]]),
      xmax = pmin(.data[[paste0("hi_", CD4)]], AXIS_T[[2]]),
      ymin = pmax(.data[[paste0("lo_", CD8)]], AXIS_T[[1]]),
      ymax = pmin(.data[[paste0("hi_", CD8)]], AXIS_T[[2]]),
      arm_label = factor(
        str_wrap(ARM_LABELS[["spectroflo"]], 24),
        levels = levels(plot_data$arm_label)
      )
    )

  # Figure
  p <- ggplot(plot_data, aes(cd4, cd8, colour = label)) +
    geom_point(size = 0.12, alpha = 0.5, stroke = 0, shape = 16) +
    geom_rect(
      data = gate_rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = quadrant),
      inherit.aes = FALSE,
      fill = NA,
      colour = "black",
      linewidth = FIG_LINE
    ) +
    facet_wrap(vars(arm_label), ncol = 4) +
    scale_colour_manual(
      values = QUAD_FILL,
      name = "SpectroFlo label",
      guide = guide_legend(override.aes = list(size = 1.5, alpha = 1))
    ) +
    coord_cartesian(xlim = AXIS_T, ylim = AXIS_T) +
    labs(
      x = sprintf("CD4, %s (arcsinh)", CD4),
      y = sprintf("CD8, %s (arcsinh)", CD8),
      title = sprintf(
        "Donor %s, CD3+ cells coloured by SpectroFlo label (%d cells, same in every panel)",
        s,
        N_DRAW
      )
    ) +
    theme(aspect.ratio = 1, legend.position = "bottom")
  save_panel(p, plot_data, OUT_DIR, paste0("cd4_cd8_", s, "_p_arms"), 183, 160)
}

# CD8T above the TNF cut per P arm and sample ------------------------------------------------------
# The table is about TNF in CD8T; a dataset whose BV510 slot is another marker (EXV: CD1c) has no such table.
if (!identical(TNF_LABEL, "TNF")) {
  stop(
    "97f_common_gs_check.R: the CD8T TNF table needs target TNF, missing in dataset ",
    DS$dataset,
    " (the TNF slot is labelled ",
    TNF_LABEL,
    "); the CD4 vs CD8 figures are written",
    call. = FALSE
  )
}
tnf <- map(ARMS, \(a) {
  read_csv(
    file.path(GATED, "gate_counts", paste0(a, ".csv")),
    show_col_types = FALSE
  )
}) |>
  list_rbind() |>
  dplyr::filter(Population == paste0(T_PARENT, "/CD8T/", TNF_LABEL, "pos")) |>
  transmute(
    arm,
    donor = name,
    n = Count,
    n_parent = ParentCount,
    pct_pos = 100 * Count / ParentCount
  ) |>
  left_join(
    select(smp, sample, condition),
    by = join_by(donor == sample),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(median_pct = median(pct_pos), .by = c(arm, condition))
stopifnot(
  "one value per arm and sample" = nrow(tnf) == length(ARMS) * nrow(smp)
)
write_csv(tnf, file.path(OUT_DIR, "cd8t_tnf_pos_p_arms.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo_97f.txt")
)
