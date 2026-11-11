# Purpose: redraw the gap-closed and marker-flag figures of sections 16 and 17 from the common-gate
#          runs (95_restyle_gap_flags.R, repointed; plans/common_gates_plan.md G7a, G11). Readouts read
#          "arcsinh 1.5". Donors keep the report palette and shape; the donor median is a bar drawn
#          behind the points. Nothing is recomputed.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_gap_flags.R (repo root)
# Inputs:  runs/2026-11-11_bmv-common-ablation-sbuv510-af405/fig_gap_closed.csv   [51c]
#          runs/2026-11-11_bmv-common-af405-forced-af/gap_closed.csv              [87_af405_forced_af_common.R]
#          runs/2026-11-11_bmv-common-af405-afremove/gap_closed.csv               [87b_af405_afremove_common.R]
#          runs/2026-11-11_bmv-marker-flags-v2-common/fig_flags_screen_CD8T.csv   [86b_apply_flag_rule_v2_common.R]
#          runs/2026-11-11_bmv-marker-flags-v2-common/fig_flags_th2_proof_CD8T.csv [86b common]
# Outputs: restyle-common/16_scatter_spectroflo/sbuv510_af405_gap_closed.{pdf,png,csv}
#          restyle-common/17_flags_th2_af/{afrefine_gap_closed,afremove_gap_closed,flags_screen_CD8T,
#            flags_th2_proof_CD8T}.{pdf,png,csv}
#          restyle-common/17_flags_th2_af/sessionInfo_95c_gap_flags.txt
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")
SCRIPT <- "95c_restyle_gap_flags.R"

# Parameters ---------------------------------------------------------------------------------------
GAP_SOLID <- c(0, 1) # 51: no gap closed, full gap closed
GAP_DASHED <- c(0.1, 0.5, 0.8) # 51: CONTRIB_FRAC, RESEMBLE_MIN, RESEMBLE_MEDIAN
HEAT_LIMIT <- 0.5 # 86b: colour scale saturates at +/-0.5 squeezed share
JITTER_WIDTH <- 0.15 # display only
JITTER_SEED <- 1
READOUT_SQ <- "TNF squeezed near 0"
READOUT_POS <- paste("TNF above", CUT_TXT)

out_16 <- fresh_dir("16_scatter_spectroflo")
out_17 <- fresh_dir("17_flags_th2_af")

# Section 16: SBUV510 + AF405 arm ------------------------------------------------------------------
src_16 <- "2026-11-11_bmv-common-ablation-sbuv510-af405/fig_gap_closed.csv"
gap_16 <- read_csv(file.path(RUNS, src_16), show_col_types = FALSE) |>
  mutate(
    arm = factor(arm, levels = rev(unique(arm))),
    readout = factor(
      readout,
      levels = c("CD8T TNF squeezed near 0", paste("CD8T TNF above", CUT_TXT))
    ),
    donor = factor_by_pal(donor, DONOR_PAL)
  )
p_16 <- ggplot(gap_16, aes(frac_closed, arm)) +
  value_grid("x") +
  geom_vline(
    xintercept = GAP_SOLID,
    colour = "grey50",
    linewidth = FIG_LINE / 2
  ) +
  geom_vline(
    xintercept = GAP_DASHED,
    linetype = "dashed",
    colour = "grey70",
    linewidth = FIG_LINE / 2
  ) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.2,
    orientation = "y"
  ) +
  geom_point(aes(colour = donor, shape = donor), size = 1.3) +
  donor_scales() +
  scale_y_discrete(labels = \(x) str_wrap(arm_label(x), 30)) +
  facet_wrap(~readout, axes = "all", axis.labels = "margins") +
  labs(
    x = "gap closed (0 = 1.8.8 zip, 1 = 1.0.2 spectra + 1.8.8 code)",
    y = NULL,
    title = "SBUV510 + AF405 1.0.2 variants: CD8T TNF gap closed",
    caption = str_wrap(
      paste(
        "Dashed lines: 0.1 = smallest contribution counted, 0.5 = every donor must reach this share of the gap,",
        "0.8 = the donor median must reach this share (criteria of script 51)."
      ),
      130
    )
  ) +
  theme(plot.caption = element_text(hjust = 0))
save_panel(p_16, gap_16, out_16, "sbuv510_af405_gap_closed", 183, 55)
add_manifest(
  SCRIPT,
  "2026-11-11_bmv-restyle/16_scatter_spectroflo/sbuv510_af405_gap_closed.png",
  file.path(out_16, "sbuv510_af405_gap_closed.png"),
  file.path(RUNS, src_16)
)

# Section 17: refined AF and AF removed ------------------------------------------------------------
src_17 <- c(
  afrefine = "2026-11-11_bmv-common-af405-forced-af/gap_closed.csv",
  afremove = "2026-11-11_bmv-common-af405-afremove/gap_closed.csv"
)
gap_titles <- c(
  afrefine = "Refined AF: CD8T TNF gap closed",
  afremove = "AF removed: CD8T TNF gap closed"
)
for (which_gap in names(src_17)) {
  d <- read_csv(file.path(RUNS, src_17[[which_gap]]), show_col_types = FALSE) |>
    mutate(
      readout = factor(readout, levels = c(READOUT_SQ, READOUT_POS)),
      donor = factor_by_pal(donor, DONOR_PAL)
    )
  p <- ggplot(d, aes(readout, frac_closed)) +
    value_grid("y") +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey60",
      linewidth = FIG_LINE / 2
    ) +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.6,
      linewidth = 0.2
    ) +
    geom_point(
      aes(colour = donor, shape = donor),
      size = 1.3,
      position = position_jitter(
        width = JITTER_WIDTH,
        height = 0,
        seed = JITTER_SEED
      )
    ) +
    donor_scales() +
    scale_x_discrete(labels = scales::label_wrap(14)) +
    labs(
      x = NULL,
      y = "Fraction of 1.8.8 to 1.0.2 gap closed",
      title = gap_titles[[which_gap]]
    )
  save_panel(p, d, out_17, paste0(which_gap, "_gap_closed"), 89, 60)
  add_manifest(
    SCRIPT,
    paste0(
      "2026-11-11_bmv-restyle/17_flags_th2_af/",
      which_gap,
      "_gap_closed.png"
    ),
    file.path(out_17, paste0(which_gap, "_gap_closed.png")),
    file.path(RUNS, src_17[[which_gap]])
  )
}

# Section 17: marker-flag heatmaps, CD8T -----------------------------------------------------------
# Change in squeezed share vs as188_joint, diverging at 0; flagged cells outlined; NA cells grey.
heat_layers <- function(d) {
  list(
    geom_tile(
      data = filter(d, flag),
      colour = "black",
      fill = NA,
      linewidth = FIG_LINE * 1.2
    ),
    geom_text(aes(label = cell_text), size = 5.5, size.unit = "pt"),
    guides(
      fill = guide_colourbar(
        barheight = unit(25, "mm"),
        barwidth = unit(3, "mm")
      )
    )
  )
}
fill_label <- "Δ squeezed\nshare vs\nas188_joint"
read_flags <- function(f) {
  read_csv(file.path(RUNS, f), col_types = cols(cell_text = col_character())) |>
    mutate(
      cell_text = if_else(
        is.na(median_d_squeezed),
        "",
        coalesce(cell_text, "")
      )
    )
}

src_screen <- "2026-11-11_bmv-marker-flags-v2-common/fig_flags_screen_CD8T.csv"
# PE has no marker in the control file (marker NA in fcs_control_file.csv)
screen <- read_flags(src_screen) |>
  mutate(
    arm = arm_label(arm),
    arm = factor(arm, levels = unique(arm[order(.data$arm)])),
    label = if_else(label == "no marker (PE)", "PE (no marker)", label),
    label = factor(label, levels = sort(unique(label)))
  ) |>
  complete_grid("label", "arm") |>
  mutate(flag = coalesce(flag, FALSE), cell_text = coalesce(cell_text, ""))
p_screen <- plot_heatmap(
  screen,
  row = "label",
  col = "arm",
  value = "median_d_squeezed",
  limit = HEAT_LIMIT,
  title = "Marker squeeze, CD8 T cells",
  x_label = "Unmixing arm",
  y_label = NULL,
  fill_label = fill_label
) +
  heat_fill(HEAT_LIMIT, fill_label) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(caption = str_wrap(paste0(NOT_TESTED, "; ", BOX_RULE), 50)) +
  theme(plot.caption = element_text(hjust = 0)) +
  heat_layers(screen)
save_panel(p_screen, screen, out_17, "flags_screen_CD8T", 89, 175)
add_manifest(
  SCRIPT,
  "2026-11-11_bmv-restyle/17_flags_th2_af/flags_screen_CD8T.png",
  file.path(out_17, "flags_screen_CD8T.png"),
  file.path(RUNS, src_screen)
)

# Th2 proof: arms down y, three markers across x. A heatmap of donor medians; no donor lines exist.
src_proof <- "2026-11-11_bmv-marker-flags-v2-common/fig_flags_th2_proof_CD8T.csv"
proof <- read_flags(src_proof) |>
  mutate(
    arm = str_wrap(arm_label(arm), 30),
    arm = factor(arm, levels = unique(arm[order(.data$arm)])),
    label = gsub(" ", "\n", label, fixed = TRUE),
    label = factor(label, levels = sort(unique(label)))
  ) |>
  complete_grid("arm", "label") |>
  mutate(flag = coalesce(flag, FALSE), cell_text = coalesce(cell_text, ""))
p_proof <- plot_heatmap(
  proof,
  row = "arm",
  col = "label",
  value = "median_d_squeezed",
  limit = HEAT_LIMIT,
  title = "AF405 proof arms, CD8 T",
  x_label = NULL,
  y_label = NULL,
  fill_label = fill_label
) +
  heat_fill(HEAT_LIMIT, fill_label) +
  labs(caption = str_wrap(paste0(NOT_TESTED, "; ", BOX_RULE), 50)) +
  theme(plot.caption = element_text(hjust = 0)) +
  heat_layers(proof)
save_panel(p_proof, proof, out_17, "flags_th2_proof_CD8T", 89, 65)
add_manifest(
  SCRIPT,
  "2026-11-11_bmv-restyle/17_flags_th2_af/flags_th2_proof_CD8T.png",
  file.path(out_17, "flags_th2_proof_CD8T.png"),
  file.path(RUNS, src_proof)
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_17, "sessionInfo_95c_gap_flags.txt")
)
