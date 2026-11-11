# Purpose: redraw the dashboard figures of section 08 from the common-gate runs (95_restyle_dashboard_extras.R,
#          repointed; plans/common_gates_plan.md G7a, G11): BV510 density in the swap test, TNF+ among IFNg+
#          cells, rank correlation with activation markers, Th2 density, the zip 1.8.8 SSM, and two extras
#          (47c cytokine-pair rho). Cut lines and labels read "arcsinh 1.5". Rho x tick labels are turned 90
#          degrees so neighbours do not overlap. proof_subsets_gap_closed has no new source (see manifest).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_dashboard_extras.R (repo root)
# Inputs:  runs/2026-11-11_bmv-swap-test-common/bv510_events.csv.gz      BV510 per cell [37b common]
#          runs/2026-11-11_bmv-common-cut15-46c/{fig_a_tnf_among_ifng,fig_b_rho_activation}.csv  [46c]
#          runs/2026-11-11_bmv-common-cut15-47c/{fig_b_density_Th2,fig_d_rho_cytokine_pairs}.csv [47c]
#          runs/2026-11-10_bmv-ssm-188/ssm_long.csv                        zip SSM; unchanged source [81]
# Outputs: restyle-common/08_dashboard_extras/{fig_bv510_density,fig_coh_tnf_among_ifng,fig_coh_rho_activation,
#            fig_oc_density_th2,fig_ssm_188_zip,fig_d_rho_cytokine_pairs}.{pdf,png,csv}, sessionInfo.txt
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")
SCRIPT <- "95c_restyle_dashboard_extras.R"

# Parameters ---------------------------------------------------------------------------------------
DENSITY_WINDOW <- c(-1.5, 4) # display window, arcsinh; the density uses every event
CYTOKINE_PAL <- c(TNF = "#D55E00", IFNg = "#0072B2") # cytokine key for the rho figures
OUT_DIR <- fresh_dir("08_dashboard_extras")
OLD <- "2026-11-11_bmv-restyle/08_dashboard_extras"
arm_levels_all <- unlist(VERSION_ROWS, use.names = FALSE)
log_fig <- function(name, source, old = TRUE, ...) {
  add_manifest(
    SCRIPT,
    if (old) file.path(OLD, paste0(name, ".png")) else NA_character_,
    file.path(OUT_DIR, paste0(name, ".png")),
    source,
    ...
  )
}

# BV510 density in the swap test -------------------------------------------------------------------
src <- file.path(RUNS, "2026-11-11_bmv-swap-test-common/bv510_events.csv.gz")
swap_arms <- c(
  "as188_joint",
  "swap_af102",
  "swap_var102",
  "s188_c102",
  "s102_c188"
)
bv510_density <- read_csv(src, show_col_types = FALSE) |>
  filter(between(.data$arcsinh, DENSITY_WINDOW[1], DENSITY_WINDOW[2])) |>
  reframe(
    as_tibble(density(
      .data$arcsinh,
      from = DENSITY_WINDOW[1],
      to = DENSITY_WINDOW[2],
      n = 256
    )[c("x", "y")]),
    .by = c(arm, sample, pop)
  ) |>
  mutate(arm = factor(.data$arm, levels = swap_arms))
stopifnot("swap-test arms as expected" = !anyNA(bv510_density$arm))
p_bv510 <- ggplot(bv510_density, aes(x, y, colour = arm)) +
  geom_vline(
    xintercept = CUT,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.3
  ) +
  geom_line(linewidth = 0.35) +
  scale_colour_manual(
    values = categorical_palette(swap_arms),
    labels = arm_label,
    name = "Arm"
  ) +
  facet_grid(pop ~ sample, axes = "all", axis.labels = "margins") +
  labs(
    x = "TNF, BV510-A (arcsinh)",
    y = "Density",
    title = paste0(
      "BV510 (TNF) density, swap test; dashed line = cut at ",
      CUT_TXT
    )
  )
save_panel(p_bv510, bv510_density, OUT_DIR, "fig_bv510_density", 183, 80)
log_fig("fig_bv510_density", src)

# TNF+ among IFNg+ and IFNg- cells -----------------------------------------------------------------
src <- file.path(
  RUNS,
  "2026-11-11_bmv-common-cut15-46c/fig_a_tnf_among_ifng.csv"
)
coh <- read_csv(src, show_col_types = FALSE) |>
  mutate(among = if_else(.data$among == "IFNg+", "among IFNg+", "among IFNg-"))
coh_out <- panel_arm_donor(
  coh,
  "pct",
  paste0("% TNF+ (", CUT_TXT, " cut)"),
  "TNF+ among IFNg+ and IFNg- cells",
  intersect(arm_levels_all, unique(coh$arm)),
  rows = "among",
  cols = "parent"
)
save_panel(
  coh_out$plot,
  coh_out$data,
  OUT_DIR,
  "fig_coh_tnf_among_ifng",
  183,
  coh_out$height_mm
)
log_fig("fig_coh_tnf_among_ifng", src)

# TNF and IFNg rank correlation with activation markers --------------------------------------------
# donor points small, donor median a diamond; x labels turned 90 degrees (they overlapped in 46c)
src <- file.path(
  RUNS,
  "2026-11-11_bmv-common-cut15-46c/fig_b_rho_activation.csv"
)
rho <- read_csv(src, show_col_types = FALSE) |>
  mutate(
    arm = factor(
      .data$arm,
      levels = rev(intersect(arm_levels_all, unique(.data$arm)))
    )
  )
rho_median <- summarise(
  rho,
  rho = median(.data$rho),
  .by = c(arm, parent, cytokine, marker)
)
p_rho <- ggplot(rho, aes(rho, arm, colour = cytokine)) +
  value_grid("x") +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_point(size = 0.4, alpha = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(
    data = rho_median,
    shape = 18,
    size = 1.6,
    position = position_dodge(width = 0.5)
  ) +
  scale_colour_manual(values = CYTOKINE_PAL, name = "Cytokine") +
  scale_y_discrete(labels = arm_label) +
  scale_x_continuous(breaks = c(-0.4, 0, 0.4), guide = guide_axis(angle = 90)) +
  facet_grid(parent ~ marker, axes = "all", axis.labels = "margins") +
  labs(
    x = "Spearman rho (diamond = donor median)",
    y = NULL,
    title = paste0(
      "Cytokine rank correlation with activation markers (",
      CUT_TXT,
      ")"
    )
  )
save_panel(p_rho, rho, OUT_DIR, "fig_coh_rho_activation", 183, 150)
log_fig(
  "fig_coh_rho_activation",
  src,
  note = "46c fig_b_rho_activation; x tick labels rotated"
)

# Cytokine-pair rho (extra) ------------------------------------------------------------------------
src <- file.path(
  RUNS,
  "2026-11-11_bmv-common-cut15-47c/fig_d_rho_cytokine_pairs.csv"
)
pairs <- read_csv(src, show_col_types = FALSE)
# 10 pairs x 2 parents = 20 panels, 5 per row; fixed x scale with breaks at -0.2, 0, 0.2, 0.4
pair_arms <- intersect(arm_levels_all, unique(pairs$arm))
pairs <- mutate(
  pairs,
  donor = factor_by_pal(.data$donor, DONOR_PAL),
  arm = factor(.data$arm, levels = pair_arms)
)
PAIRS_HEIGHT_MM <- 235 # 4 rows of panels, 11 arms each; Nature page height limit
p_pairs <- ggplot(pairs, aes(rho, arm)) +
  value_grid("x") +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.2,
    orientation = "y"
  ) +
  geom_point(aes(colour = donor, shape = donor), size = 0.6) +
  donor_scales() +
  scale_y_discrete(limits = rev, labels = arm_label) +
  scale_x_continuous(
    breaks = c(-0.2, 0, 0.2, 0.4),
    guide = guide_axis(angle = 90)
  ) +
  facet_wrap(
    vars(parent, pair),
    ncol = 5,
    axes = "all",
    axis.labels = "margins"
  ) +
  labs(
    x = "Spearman rho",
    y = NULL,
    title = paste0("Rank correlation of cytokine pairs (", CUT_TXT, ")")
  ) +
  theme(panel.spacing.x = unit(3, "mm"))
save_panel(
  p_pairs,
  pairs,
  OUT_DIR,
  "fig_d_rho_cytokine_pairs",
  183,
  PAIRS_HEIGHT_MM
)
log_fig(
  "fig_d_rho_cytokine_pairs",
  src,
  old = FALSE,
  note = "extra from 47c fig_d; 5 panels per row, x breaks -0.2 to 0.4, smaller markers"
)

# Th2 density per arm ------------------------------------------------------------------------------
src <- file.path(RUNS, "2026-11-11_bmv-common-cut15-47c/fig_b_density_Th2.csv")
th2 <- read_csv(src, show_col_types = FALSE)
p_th2 <- ggplot(th2, aes(value, density, group = arm, colour = version)) +
  geom_vline(
    xintercept = CUT,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.3
  ) +
  geom_line(linewidth = 0.3) +
  scale_colour_manual(
    values = categorical_palette(unique(th2$version)),
    name = "Code version"
  ) +
  facet_grid(parent ~ donor, axes = "all", axis.labels = "margins") +
  labs(
    x = "Th2, eFluor 450-A (arcsinh)",
    y = "Density",
    title = paste0("Th2 density per arm; dashed line = cut at ", CUT_TXT)
  )
save_panel(p_th2, th2, OUT_DIR, "fig_oc_density_th2", 183, 80)
log_fig("fig_oc_density_th2", src)

# AF405 variant subsets: no new source -------------------------------------------------------------
add_manifest(
  SCRIPT,
  file.path(OLD, "proof_subsets_gap_closed.png"),
  NA_character_,
  "2026-11-10_bmv-af405-mechanism/scores/fig_gap_closed_source.csv",
  status = "no source",
  note = "the af405_102_* subset arms exist in common-gated/gs but no common-gate 51 score covers them"
)

# Zip 1.8.8 spillover spread matrix: no gate dependence, source unchanged ----------------------------
src <- file.path(RUNS, "2026-11-10_bmv-ssm-188/ssm_long.csv")
ssm <- read_csv(src, show_col_types = FALSE)
dye_levels <- unique(ssm$row_dye)
ssm <- mutate(
  ssm,
  row_dye = factor(.data$row_dye, levels = rev(dye_levels)),
  col_dye = factor(.data$col_dye, levels = dye_levels)
)
p_ssm <- ggplot(ssm, aes(col_dye, row_dye, fill = ssm)) +
  geom_tile(colour = "white", linewidth = 0.1) +
  scale_fill_viridis_c(name = "SSM", na.value = "grey90") +
  coord_equal() +
  labs(
    x = "Receiving dye",
    y = "Spreading dye",
    title = "AutoSpectral 1.8.8 spillover spread (zip)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.text.y = element_text(size = 5)
  )
save_panel(p_ssm, ssm, OUT_DIR, "fig_ssm_188_zip", 183, 175)
log_fig(
  "fig_ssm_188_zip",
  src,
  status = "unchanged source",
  note = "SSM is gate-independent"
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
