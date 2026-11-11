# Purpose: sanity check of the common GatingSets (97e): CD4 vs CD8 of donor A3 in the 11 P arms, every cell
#          coloured by its SpectroFlo label (G11). Gate lines appear on the SpectroFlo panel only. Also writes
#          the % of CD8T above the fixed TNF cut (arcsinh 1.5) per P arm and donor, with the donor median.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/97f_common_gs_check.R
#          from the repo root, after 97e_common_gs.R
# Inputs:  experiments/runs/2026-11-11_bmv-common-gated/{gs,gate_counts}/<arm>   [97e]
#          experiments/runs/2026-11-11_bmv-common-gates-clean3/cuts.csv   SpectroFlo gate bounds [97d]
#          experiments/benchmarks/05_bmv-unmix-compare/{arm_labels.R,lib/theme.R}
# Outputs: <out_dir>/cd4_cd8_A3_p_arms.{png,pdf,csv}   the figure and its source data (events drawn)
#          <out_dir>/cd8t_tnf_pos_p_arms.csv           arm, donor, n, n_parent, pct_pos, median_pct per arm
#          <out_dir>/sessionInfo_97f.txt

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
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS, VERSION_ROWS
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set, save_panel, FIG_LINE

# Parameters ---------------------------------------------------------------------------------------
GATED <- "experiments/runs/2026-11-11_bmv-common-gated"
CUTS <- "experiments/runs/2026-11-11_bmv-common-gates-clean3/cuts.csv"
OUT_DIR <- "experiments/runs/2026-11-11_bmv-common-gated-check" # new dated dir; must not exist
DONOR <- "A3"
ARMS <- unlist(VERSION_ROWS, use.names = FALSE) # the 11 P arms, in version-group order
CD4 <- "Spark Blue 574-A" # channels as in 19a_gate_arm.R
CD8 <- "SBUV445-A"
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

# Events: CD3+ cells of A3 with their SpectroFlo quadrant label --------------------------------------
events <- map(set_names(ARMS), \(a) {
  gh <- load_gs(file.path(GATED, "gs", a))[[DONOR]]
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

# Gate lines: SpectroFlo panel only, from the pooled clean3 bounds ---------------------------------
gate_rects <- read_csv(CUTS, show_col_types = FALSE) |>
  dplyr::filter(
    version == "clean3_pooled",
    sample == DONOR,
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

# Figure -------------------------------------------------------------------------------------------
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
      DONOR,
      N_DRAW
    )
  ) +
  theme(aspect.ratio = 1, legend.position = "bottom")
save_panel(p, plot_data, OUT_DIR, "cd4_cd8_A3_p_arms", 183, 160)

# CD8T above the fixed TNF cut per P arm and donor -------------------------------------------------
tnf <- map(ARMS, \(a) {
  read_csv(
    file.path(GATED, "gate_counts", paste0(a, ".csv")),
    show_col_types = FALSE
  )
}) |>
  list_rbind() |>
  dplyr::filter(Population == paste0(T_PARENT, "/CD8T/TNFpos")) |>
  transmute(
    arm,
    donor = name,
    n = Count,
    n_parent = ParentCount,
    pct_pos = 100 * Count / ParentCount
  ) |>
  mutate(median_pct = median(pct_pos), .by = arm)
stopifnot("one value per arm and donor" = nrow(tnf) == length(ARMS) * 5)
write_csv(tnf, file.path(OUT_DIR, "cd8t_tnf_pos_p_arms.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo_97f.txt")
)
