# Purpose: Nature-spec redraw on the common gates (plans/common_gates_plan.md G1a, G10b, G11) of the CD4 vs CD8
#          gate check (old 95_restyle_events_gate_cd4_cd8.R) for every donor: CD3+ cells of the 11 main arms,
#          each cell coloured by its SpectroFlo label (CD4T, CD8T, DP, DN). The same rows sit in every arm,
#          so panels differ only by unmixing. Gate lines appear on the SpectroFlo panel only, with % of CD3+
#          per label. Also draws the check figure of 97f_common_gs_check.R (donor A3, no % labels).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_events_gate_cd4_cd8.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-11_bmv-common-gated/gs/<arm>/   common-gate GatingSets [97e]
#          experiments/runs/2026-11-11_bmv-common-gates-clean3/cuts.csv   SpectroFlo pooled gate bounds [97d]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R   ARM_LABELS, VERSION_ROWS
# Outputs: experiments/runs/2026-11-11_bmv-restyle-common/02_gates/cd4_cd8_tmix_pooled/gate_cd4_cd8_<donor>.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle-common/02_gates/check/cd4_cd8_A3_p_arms.{pdf,png,csv}

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths()))
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
source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS, VERSION_ROWS

# Parameters ---------------------------------------------------------------------------------------
COMMON_GS <- "experiments/runs/2026-11-11_bmv-common-gated/gs"
CUTS <- "experiments/runs/2026-11-11_bmv-common-gates-clean3/cuts.csv"
CD4 <- "Spark Blue 574-A" # channels as in 19a_gate_arm.R
CD8 <- "SBUV445-A"
AXIS_T <- c(-3, 6) # display window, arcsinh; % use every CD3+ cell
N_DRAW <- 15000 # cells drawn per panel; the same rows in every arm
SEED <- 42 # row subsample
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
LABEL_INSET <- 0.06 # % labels sit this fraction of the axis range in from the panel edges, clear of the axes and cut lines
T_PARENT <- "/nonDebris/singlets/live/CD3pos"
# label colours as 97d / 97f; a label key, so no other use of these colours
QUAD_FILL <- c(
  CD4T = "#2166ac",
  CD8T = "#b2182b",
  DP = "#762a83",
  DN = "#7f7f7f"
)
ARMS <- unlist(VERSION_ROWS, use.names = FALSE) # the 11 main arms, in figure order

OUT_GATES <- fresh_dir("02_gates/cd4_cd8_tmix_pooled")
OUT_CHECK <- fresh_dir("02_gates/check")

# CD3+ events with their SpectroFlo label, per arm and donor ----------------------------------------
events <- map(set_names(ARMS), \(a) {
  gs <- load_gs(file.path(COMMON_GS, a))
  map(set_names(sampleNames(gs)), \(d) {
    gh <- gs[[d]]
    root <- exprs(gh_pop_get_data(gh, "root"))
    in_cd3 <- gh_pop_get_indices(gh, T_PARENT)
    quad <- map(names(QUAD_FILL), \(q) {
      gh_pop_get_indices(gh, paste0(T_PARENT, "/", q))
    }) |>
      set_names(names(QUAD_FILL)) |>
      do.call(what = cbind)
    stopifnot(
      "every CD3+ cell is in one label" = all(rowSums(quad)[in_cd3] == 1)
    )
    tibble(
      row = which(in_cd3),
      cd4 = root[in_cd3, CD4],
      cd8 = root[in_cd3, CD8],
      label = names(QUAD_FILL)[max.col(quad, ties.method = "first")][in_cd3]
    )
  }) |>
    list_rbind(names_to = "donor")
}) |>
  list_rbind(names_to = "arm")
# the rows are the same cells in every arm (G8)
stopifnot(
  "rows are the same in every arm" = n_distinct(map(
    split(events$row, paste(events$arm, events$donor)),
    identity
  )) ==
    n_distinct(events$donor)
)

# Gate lines from the SpectroFlo pooled cuts: one CD4T/CD8T boundary, same for every donor --------
cuts <- read_csv(CUTS, show_col_types = FALSE) |>
  dplyr::filter(
    version == "clean3_pooled",
    gate %in% c("CD4T", "CD8T"),
    dim %in% c(CD4, CD8)
  ) |>
  select(sample, gate, dim, lo, hi)
cut_at <- \(d, g, dm, side) {
  cuts[[side]][cuts$sample == d & cuts$gate == g & cuts$dim == dm]
}
lines_of <- \(d) {
  x_lo <- cut_at(d, "CD4T", CD4, "lo") # CD4T lower CD4 bound
  x_hi <- cut_at(d, "CD8T", CD4, "hi") # CD8T upper CD4 bound
  y_cut <- cut_at(d, "CD8T", CD8, "lo") # CD8 bound between CD4T and CD8T
  tibble(
    x = c(AXIS_T[[1]], AXIS_T[[1]] + 0, x_lo, x_hi),
    xend = c(AXIS_T[[2]], x_lo, x_lo, x_hi),
    y = c(y_cut, y_cut, AXIS_T[[1]], y_cut),
    yend = c(y_cut, y_cut, y_cut, AXIS_T[[2]])
  ) |>
    dplyr::slice(-2)
}

# Subtitle: the pooled cut values (arcsinh), the same for every donor
cut_text <- \(d) {
  sprintf(
    "Pooled gate cuts (arcsinh): CD4 %.3f / %.3f, CD8 %.3f",
    cut_at(d, "CD4T", CD4, "lo"),
    cut_at(d, "CD8T", CD4, "hi"),
    cut_at(d, "CD8T", CD8, "lo")
  )
}

arm_levels <- str_wrap(ARM_LABELS[ARMS], 24)
ref_level <- str_wrap(ARM_LABELS[["spectroflo"]], 24)

draw_donor <- \(d, pct, out_dir, stem, title) {
  set.seed(SEED)
  rows_d <- unique(events$row[events$donor == d])
  rows_drawn <- sample(rows_d, N_DRAW)
  plot_df <- events |>
    dplyr::filter(donor == d, row %in% rows_drawn) |>
    mutate(
      arm_label = factor(str_wrap(ARM_LABELS[arm], 24), levels = arm_levels),
      label = factor(label, levels = names(QUAD_FILL))
    ) |>
    arrange(label %in% c("DP", "DN")) # rare labels drawn last so they stay visible
  seg <- lines_of(d) |>
    mutate(arm_label = factor(ref_level, levels = arm_levels))
  p <- ggplot(plot_df, aes(cd4, cd8, colour = label)) +
    geom_point(size = 0.3, alpha = 0.6, stroke = 0, shape = 16) +
    geom_segment(
      data = seg,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = FIG_LINE
    ) +
    facet_wrap(
      vars(arm_label),
      ncol = 4,
      axes = "all",
      axis.labels = "margins"
    ) +
    scale_colour_manual(
      values = QUAD_FILL,
      name = "SpectroFlo label",
      guide = guide_legend(override.aes = list(size = 1.5, alpha = 1))
    ) +
    coord_cartesian(xlim = AXIS_T, ylim = AXIS_T) +
    labs(
      subtitle = cut_text(d),
      x = sprintf("CD4, %s (arcsinh)", CD4),
      y = sprintf("CD8, %s (arcsinh)", CD8),
      title = title
    ) +
    theme(aspect.ratio = 1, panel.spacing = unit(1.2, "mm"))
  if (!is.null(pct)) {
    p <- p +
      geom_text(
        data = pct |>
          mutate(arm_label = factor(ref_level, levels = arm_levels)),
        aes(x = x, y = y, label = txt, hjust = hjust, vjust = vjust),
        inherit.aes = FALSE,
        size = LABEL_PT,
        size.unit = "pt",
        family = FIG_FONT
      )
  }
  save_panel(p, plot_df, out_dir, stem, 183, 150)
  invisible(plot_df)
}

# One figure per donor, with % of CD3+ per label printed on the SpectroFlo panel --------------------
walk(sort(unique(events$donor)), \(d) {
  ref <- events |> dplyr::filter(arm == "spectroflo", donor == d)
  pct <- ref |>
    summarise(n = n(), .by = label) |>
    mutate(
      pct = 100 * n / sum(n),
      txt = sprintf("%s %.1f%%", label, pct),
      # CD8T top left, CD4T bottom right, DP top right, DN bottom left
      x = c(
        CD8T = AXIS_T[1] + LABEL_INSET * diff(AXIS_T),
        CD4T = AXIS_T[2] - LABEL_INSET * diff(AXIS_T),
        DP = AXIS_T[2] - LABEL_INSET * diff(AXIS_T),
        DN = AXIS_T[1] + LABEL_INSET * diff(AXIS_T)
      )[as.character(label)],
      y = c(
        CD8T = AXIS_T[2] - LABEL_INSET * diff(AXIS_T),
        CD4T = AXIS_T[1] + LABEL_INSET * diff(AXIS_T),
        DP = AXIS_T[2] - LABEL_INSET * diff(AXIS_T),
        DN = AXIS_T[1] + LABEL_INSET * diff(AXIS_T)
      )[as.character(label)],
      hjust = c(CD8T = 0, CD4T = 1, DP = 1, DN = 0)[as.character(label)],
      vjust = c(CD8T = 1, CD4T = 0, DP = 1, DN = 0)[as.character(label)]
    )
  stem <- paste0("gate_cd4_cd8_", d)
  draw_donor(
    d,
    pct,
    OUT_GATES,
    stem,
    sprintf(
      "CD3+ T cells, donor %s, coloured by SpectroFlo label (%d cells per panel)",
      d,
      N_DRAW
    )
  )
  add_manifest(
    "95c_restyle_events_gate_cd4_cd8",
    file.path(
      "2026-11-11_bmv-restyle/02_gates/cd4_cd8_tmix_pooled",
      paste0(stem, ".png")
    ),
    file.path(OUT_GATES, paste0(stem, ".png")),
    paste(COMMON_GS, CUTS, sep = "; ")
  )
})

# The 97f check figure: donor A3, labels only, x axes on the bottom panels of each column -----------
draw_donor(
  "A3",
  NULL,
  OUT_CHECK,
  "cd4_cd8_A3_p_arms",
  sprintf(
    "Check, donor A3: CD3+ cells coloured by SpectroFlo label (%d cells, same in every panel)",
    N_DRAW
  )
)
add_manifest(
  "95c_restyle_events_gate_cd4_cd8",
  NA_character_,
  file.path(OUT_CHECK, "cd4_cd8_A3_p_arms.png"),
  paste(COMMON_GS, CUTS, sep = "; "),
  note = "97f check figure; no old restyle counterpart"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_GATES, "sessionInfo.txt")
)
