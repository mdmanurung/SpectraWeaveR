# Purpose: test whether "TNF+ among IFNg+" separates real TNF from artefact (user, 2026-11-10, grill Q11).
#          Per arm and donor, CD8T: TNF+ among IFNg+, TNF+ among IFNg- (derived from the scorer's
#          % TNF+, % IFNg+ and % TNF+ among IFNg+) and the TNF-IFNg odds ratio. Reads scorer CSVs only.
#          Common-gate copy (plans/common_gates_plan.md, batch 2): reads the 51c rescores at cut 1.5 (G7a).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/66_coherence_or_common.R   (repo root)
# Inputs:  experiments/runs/<run>/pct_pos.csv        % TNF+ and % IFNg+ per donor, arm, pop [51c_score_ablation_common.R]
#          experiments/runs/<run>/tnf_among_ifng.csv % TNF+ among IFNg+ per donor, arm, pop [51c_score_ablation_common.R]
#          for <run> in SOURCES below
# Outputs: experiments/runs/2026-11-11_bmv-common-coherence-or/coherence_or.csv   per donor and arm: the 3 readouts
#          experiments/runs/2026-11-11_bmv-common-coherence-or/fig_coherence_or.png  readouts per arm, donors as points
#          experiments/runs/2026-11-11_bmv-common-coherence-or/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))

# Parameters ---------------------------------------------------------------------------------------
DONORS <- c("A3", "A6", "B11", "B2", "D1")
DONOR_SHAPES <- setNames(c(16, 17, 15, 18, 8), DONORS) # as the approved demo (2026-11-11)
POP <- "CD8T" # the population where TNF differs most between versions (INDEX section 05)
# scorer run dirs and the experiment label each arm is shown under; base arms are taken from the first
SOURCES <- c(
  "Block A swaps" = "2026-11-11_bmv-common-ablation",
  "Rebuild chain" = "2026-11-11_bmv-common-chain-all",
  "Rebuild chain " = "2026-11-11_bmv-common-chain7a8",
  "AF405 proof" = "2026-11-11_bmv-common-af405-proof-score",
  "Passes grid" = "2026-11-11_bmv-common-passes-grid"
)
BASE_ARMS <- c("as188_joint", "s102_c188", "as102_var")
OUT_DIR <- "experiments/runs/2026-11-11_bmv-common-coherence-or"

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Read scorer tables ------------------------------------------------------------------------------
# one row per run, donor and arm; base arms kept once (from the first source)
read_run <- \(run) {
  inner_join(
    read_csv(
      file.path("experiments/runs", run, "pct_pos.csv"),
      show_col_types = FALSE
    ),
    read_csv(
      file.path("experiments/runs", run, "tnf_among_ifng.csv"),
      show_col_types = FALSE
    ),
    by = join_by(donor, arm, pop),
    relationship = "one-to-one",
    unmatched = "error"
  )
}
scored <- map(SOURCES, read_run) |>
  list_rbind(names_to = "experiment") |>
  mutate(
    experiment = if_else(arm %in% BASE_ARMS, "Reference", trimws(experiment))
  ) |>
  distinct(donor, arm, pop, .keep_all = TRUE) |>
  filter(pop == POP)

# Derive TNF+ among IFNg- and the odds ratio -----------------------------------------------------
# P(TNF+ | IFNg-) = (P(TNF+) - P(TNF+ | IFNg+) P(IFNg+)) / (1 - P(IFNg+))
coherence <- scored |>
  mutate(
    p_pos = pct_tnf_among_ifng / 100,
    p_neg = (pct_tnf / 100 - p_pos * pct_ifng / 100) / (1 - pct_ifng / 100),
    odds_ratio = (p_pos / (1 - p_pos)) / (p_neg / (1 - p_neg))
  ) |>
  select(
    experiment,
    arm,
    donor,
    pct_tnf,
    pct_ifng,
    n_ifng_pos,
    p_pos,
    p_neg,
    odds_ratio
  )
write_csv(coherence, file.path(OUT_DIR, "coherence_or.csv"))

# Plot frame -------------------------------------------------------------------------------------
# long format; arms ordered by the as-observed TNF+ among IFNg+ median within each experiment
readout_levels <- c(
  p_pos = "TNF+ among IFNg+ (%)",
  p_neg = "TNF+ among IFNg- (%)",
  log2_or = "TNF-IFNg log2 odds ratio"
)
# arm order: median TNF+ among IFNg+ (the readout under test), low to high
arm_order <- coherence |>
  summarise(m = median(p_pos), .by = arm) |>
  arrange(m) |>
  pull(arm)
plot_frame <- coherence |>
  mutate(
    p_pos = 100 * p_pos,
    p_neg = 100 * p_neg,
    log2_or = log2(odds_ratio)
  ) |>
  pivot_longer(
    c(p_pos, p_neg, log2_or),
    names_to = "readout",
    values_to = "value"
  ) |>
  mutate(
    readout = factor(readout_levels[readout], levels = unname(readout_levels)),
    experiment = factor(
      experiment,
      levels = c(
        "Reference",
        "Block A swaps",
        "Rebuild chain",
        "AF405 proof",
        "Passes grid"
      )
    ),
    arm = factor(arm, levels = arm_order)
  )
medians <- plot_frame |>
  summarise(value = median(value), .by = c(experiment, arm, readout))

# Figure ---------------------------------------------------------------------------------------------
p <- ggplot(plot_frame, aes(value, arm)) +
  geom_point(aes(colour = donor, shape = donor), size = 1.1, alpha = 0.8) +
  geom_point(data = medians, shape = 18, size = 2.2) +
  scale_colour_manual(values = categorical_palette(DONORS), name = "Donor") +
  scale_shape_manual(values = DONOR_SHAPES, name = "Donor") +
  facet_grid(experiment ~ readout, scales = "free", space = "free_y") +
  labs(
    x = NULL,
    y = NULL,
    title = "CD8T TNF coherence per arm, common gates (diamond = donor median)"
  ) +
  value_grid("x")
save_panel(
  p,
  plot_frame,
  OUT_DIR,
  "fig_coherence_or",
  width_mm = 183,
  height_mm = 230
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
