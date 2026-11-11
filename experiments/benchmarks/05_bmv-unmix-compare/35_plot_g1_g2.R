# Purpose: make the G1 (held-out controls) and G2 (simulation) figures from the finished tables, so
#          the two benchmarks can be cited. Reads tables only; nothing is recomputed.
#          G1: per arm, log2(arm / as188_joint) of off-target spread and |bias|, one point per control
#              (the control's median over its off-target channels, as in heldout_summary.csv).
#          G2: per arm and dye, seed-mean error as a ratio to the 1.8.8 joint arm (33 dyes).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/35_plot_g1_g2.R <run_dir>
#          from the repo root. run_dir must not exist; use experiments/runs/2026-11-10_bmv-g1g2-figures.
# Inputs:  experiments/runs/2026-11-10_bmv-heldout/heldout_summary.csv [26_heldout_controls.R, 27_heldout_score.R]
#          experiments/runs/2026-11-10_bmv-simulation/sim_scores.csv [27_simulation.R]
# Outputs: <run_dir>/fig_g1_heldout.png          G1 figure
#          <run_dir>/fig_g1_heldout_source.csv   one row per arm, control and metric
#          <run_dir>/fig_g2_simulation.png       G2 figure
#          <run_dir>/fig_g2_simulation_source.csv  one row per arm, dye and metric: seed-mean ratio, seed range
#          <run_dir>/fig_g2_simulation_seeds.csv   same, one row per seed (the faint points)
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Parameters ------------------------------------------------------------------------------
REF_ARM <- "as188_joint" # reference arm (README, approved choices); its ratio is 1 (G2) or 0 on log2 (G1) by construction, so it is not drawn
MIN_BOX <- 5 # draw a boxplot only for groups with at least this many controls; smaller groups show points only
JITTER_SEED <- 42 # jitter only moves points sideways; same seed gives the same picture
G1_ARMS <- c(
  # arm ID -> label; real-data arms, G1
  as102_basic = "as102_basic\n1.0.2 code, own spectra, no variants",
  as102_var = "as102_var\n1.0.2 code, own spectra, variants",
  s188_c102 = "s188_c102\n1.8.8 spectra, 1.0.2 code",
  s102_c188 = "s102_c188\n1.0.2 spectra, 1.8.8 code"
)
# G2 arm folders are NOT the real-data arms of the same name (args.txt: all use 1.8.8 spectra and AF)
G2_ARMS <- c(
  ols = "OLS",
  as102_basic = "1.0.2 code + 1.8.8 spectra,\nno variants",
  as102_var = "1.0.2 code + 1.8.8 spectra\n+ variants (= s188_c102 recipe)"
)
G2_REF_LABEL <- "1.8.8 joint (as188_joint recipe)"
G1_GROUPS <- c(
  # group in heldout_summary.csv -> facet label
  `cell single-stain` = "Cell single stains",
  `bead single-stain` = "Bead single stains",
  `Unstained (Cells)` = "Unstained cells",
  `Negative_1 (Beads)` = "Unstained beads",
  `H12 Unstained_BMV_Aurora_3dCD328_Batch1` = "H12 unstained cells\n(NOT held out: AF circular)"
)

# Output folder ---------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
if (dir.exists(run_dir)) {
  stop("run_dir exists: ", run_dir)
}
dir.create(run_dir, recursive = TRUE)

# G1: one row per control, long over the two metrics --------------------------------------
heldout <- read_csv(
  "experiments/runs/2026-11-10_bmv-heldout/heldout_summary.csv",
  show_col_types = FALSE
)
g1 <- heldout |>
  filter(arm != REF_ARM) |>
  select(arm, control, group, log2_spread, log2_bias) |>
  pivot_longer(
    starts_with("log2_"),
    names_to = "metric",
    names_prefix = "log2_",
    values_to = "log2_ratio"
  ) |>
  mutate(
    metric = c(spread = "Off-target spread", bias = "Off-target |bias|")[
      metric
    ],
    n_controls = n(),
    .by = c(arm, group, metric)
  )
# Labels are looked up by name; an unknown arm or group gives NA and stops here
g1 <- g1 |>
  mutate(
    arm_label = unname(G1_ARMS[arm]),
    group_label = unname(G1_GROUPS[group])
  )
stopifnot(
  !anyNA(g1$arm_label),
  !anyNA(g1$group_label),
  all(is.finite(g1$log2_ratio))
)
g1 <- g1 |>
  mutate(
    arm_label = factor(arm_label, levels = rev(unname(G1_ARMS))),
    group_label = factor(group_label, levels = unname(G1_GROUPS)),
    metric = factor(
      metric,
      levels = c("Off-target spread", "Off-target |bias|")
    )
  )

g1_caption <- str_c(
  "Each point is one control: the median over its off-target channels of the held-out half, as log2(arm / as188_joint).\n",
  "Spectra were re-extracted from a training half of each control (seed 42); the other half was unmixed and scored. ",
  "Boxes (median, IQR) only for groups of ",
  MIN_BOX,
  " or more controls.\n",
  "Dashed line: no change from as188_joint; below 0 means lower than as188_joint. ",
  "H12 is NOT held out: its autofluorescence spectrum was extracted from the same file, so its score is circular.\n",
  "as188_joint is the reference (0 by construction) and is not drawn."
)
g1_plot <- ggplot(g1, aes(log2_ratio, arm_label)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_boxplot(
    data = filter(g1, n_controls >= MIN_BOX),
    outlier.shape = NA,
    width = 0.5,
    fill = NA
  ) +
  geom_point(
    position = position_jitter(height = 0.15, width = 0, seed = JITTER_SEED),
    alpha = 0.7,
    size = 1.6
  ) +
  facet_grid(metric ~ group_label, scales = "free_x") +
  labs(
    x = "log2(arm / as188_joint), per control",
    y = NULL,
    title = "G1 held-out controls: error relative to the 1.8.8 joint recipe",
    caption = g1_caption
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.caption = element_text(hjust = 0, size = 8),
    strip.text = element_text(size = 8)
  )
ggsave(
  file.path(run_dir, "fig_g1_heldout.png"),
  g1_plot,
  width = 15,
  height = 7,
  dpi = 150
)
g1 |>
  arrange(metric, group_label, arm_label, control) |>
  select(arm, control, group, n_controls, metric, log2_ratio) |>
  write_csv(file.path(run_dir, "fig_g1_heldout_source.csv"))

# G2: one row per seed, dye, arm and metric ------------------------------------------------
sim <- read_csv(
  "experiments/runs/2026-11-10_bmv-simulation/sim_scores.csv",
  show_col_types = FALSE
)
# Four metrics: RMSE (all events), |bias| (all events, dim layer), negative spread (negative cells)
g2_values <- bind_rows(
  sim |>
    filter(stratum == "all") |>
    transmute(seed, arm, fluor, metric = "RMSE (all events)", value = rmse),
  sim |>
    filter(stratum == "all") |>
    transmute(
      seed,
      arm,
      fluor,
      metric = "|bias| (all events)",
      value = abs(bias)
    ),
  sim |>
    filter(stratum == "dim") |>
    transmute(
      seed,
      arm,
      fluor,
      metric = "|bias| (dim layer)",
      value = abs(bias)
    ),
  sim |>
    filter(stratum == "neg") |>
    transmute(seed, arm, fluor, metric = "Negative spread", value = spread)
)
g2_ref <- g2_values |>
  filter(arm == REF_ARM) |>
  select(seed, fluor, metric, ref_value = value)
g2_seeds <- g2_values |>
  filter(arm != REF_ARM) |>
  left_join(
    g2_ref,
    by = c("seed", "fluor", "metric"),
    relationship = "many-to-one"
  ) |>
  mutate(
    arm_label = unname(G2_ARMS[arm]),
    ratio = value / ref_value
  )
stopifnot(!anyNA(g2_seeds$arm_label), all(is.finite(g2_seeds$ratio)))
# Dye point = ratio of seed means; faint points = per-seed ratios against the same seed of the reference
g2_dyes <- g2_seeds |>
  summarise(
    mean_value = mean(value),
    mean_ref = mean(ref_value),
    ratio_min_seed = min(ratio),
    ratio_max_seed = max(ratio),
    .by = c(arm, arm_label, fluor, metric)
  ) |>
  mutate(ratio = mean_value / mean_ref)
g2_levels <- c(
  "RMSE (all events)",
  "|bias| (all events)",
  "|bias| (dim layer)",
  "Negative spread"
)
g2_seeds <- g2_seeds |>
  mutate(
    arm_label = factor(arm_label, levels = rev(unname(G2_ARMS))),
    metric = factor(metric, levels = g2_levels)
  )
g2_dyes <- g2_dyes |>
  mutate(
    arm_label = factor(arm_label, levels = rev(unname(G2_ARMS))),
    metric = factor(metric, levels = g2_levels)
  )

g2_caption <- str_c(
  "Each bold point is one dye (",
  n_distinct(g2_dyes$fluor),
  "): ratio of the mean over 3 seeds (42-44) to the same mean for the ",
  G2_REF_LABEL,
  ". Faint points: per-seed ratios (seed range).\n",
  "Boxes: median and IQR over dyes. Dashed line at 1: same error as the reference; below 1 means lower error. ",
  "Log axis. 1e5 simulated cells per seed, known truth.\n",
  "Negative spread is scored on truly negative cells. Dim layer: dyes at 0.005 relative abundance. ",
  "The reference arm is 1 by construction and is not drawn.\n",
  "All arms use the 1.8.8 spectra and AF; the simulation folder names as102_basic and as102_var are not the real-data arms of those names."
)
g2_plot <- ggplot(g2_dyes, aes(ratio, arm_label)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(
    data = g2_seeds,
    aes(ratio),
    colour = "grey60",
    alpha = 0.3,
    size = 0.8,
    position = position_jitter(height = 0.2, width = 0, seed = JITTER_SEED)
  ) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA) +
  geom_point(
    position = position_jitter(height = 0.2, width = 0, seed = JITTER_SEED),
    size = 1.6,
    alpha = 0.8
  ) +
  scale_x_log10() +
  facet_wrap(vars(metric), nrow = 1, scales = "free_x") +
  labs(
    x = str_c("ratio to ", G2_REF_LABEL, " (log scale)"),
    y = NULL,
    title = "G2 simulation: error relative to the 1.8.8 joint recipe, per dye",
    caption = g2_caption
  ) +
  theme_bw(base_size = 10) +
  theme(plot.caption = element_text(hjust = 0, size = 8))
ggsave(
  file.path(run_dir, "fig_g2_simulation.png"),
  g2_plot,
  width = 15,
  height = 5,
  dpi = 150
)
g2_dyes |>
  arrange(metric, arm_label, fluor) |>
  select(
    arm,
    fluor,
    metric,
    ratio,
    ratio_min_seed,
    ratio_max_seed,
    mean_value,
    mean_ref
  ) |>
  write_csv(file.path(run_dir, "fig_g2_simulation_source.csv"))
g2_seeds |>
  arrange(metric, arm_label, fluor, seed) |>
  select(arm, fluor, metric, seed, value, ref_value, ratio) |>
  write_csv(file.path(run_dir, "fig_g2_simulation_seeds.csv"))

# Session ---------------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
