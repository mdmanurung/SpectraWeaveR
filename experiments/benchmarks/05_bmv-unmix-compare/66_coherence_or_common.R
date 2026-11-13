# Purpose: test whether "TNF+ among IFNg+" separates real TNF from artefact (user, 2026-11-10, grill Q11).
#          Per arm and donor, CD8T: TNF+ among IFNg+, TNF+ among IFNg- (derived from the scorer's
#          % TNF+, % IFNg+ and % TNF+ among IFNg+) and the TNF-IFNg odds ratio. Reads scorer CSVs only.
#          Common-gate copy (plans/common_gates_plan.md, batch 2): reads the 51c rescores at cut 1.5 (G7a).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/66_coherence_or_common.R   (repo root)
#          Needs the TNF and IFNg readouts: stops when the dataset's targets are not TNF and IFNg (EXV).
# Inputs:  run_in(<key>)/pct_pos.csv        % TNF+ and % IFNg+ per sample, arm, pop [51c_score_ablation_common.R]
#          run_in(<key>)/tnf_among_ifng.csv % TNF+ among IFNg+ per sample, arm, pop [51c_score_ablation_common.R]
#          for <key> in SOURCES below (registry keys of the 51c run folders)
#          datasets/<SW_DATASET>.json (targets), <ds>_samples.csv
# Outputs: run_out("common-coherence-or")/coherence_or.csv   per sample and arm: the 3 readouts (24H: + condition)
#          run_out("common-coherence-or")/fig_coherence_or.png  readouts per arm, donors as points (24H: by condition)
#          run_out("common-coherence-or")/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))

# Parameters ---------------------------------------------------------------------------------------
TNF_LABEL <- cfg("targets")$tnf_label
IFNG_LABEL <- cfg("targets")$ifng_label
if (!identical(c(TNF_LABEL, IFNG_LABEL), c("TNF", "IFNg"))) {
  stop(
    "66_coherence_or_common.R: needs TNF and IFNg, missing in dataset ",
    DS$dataset,
    " (the two dye slots are ",
    TNF_LABEL,
    " and ",
    IFNG_LABEL,
    ")",
    call. = FALSE
  )
}
# 51c writes the sample id in `donor`; 24H figures plot the donor (points) over conditions
SAMPLES <- samples_tbl()
COND_KEY <- if (is.null(DS$conditions)) character() else "condition"
PLOT_DONOR <- set_names(
  if (is.null(DS$conditions)) SAMPLES$sample else SAMPLES$donor,
  SAMPLES$sample
)
DONORS <- sort(unique(PLOT_DONOR))
DONOR_SHAPES <- setNames(rep_len(c(16, 17, 15, 18, 8), length(DONORS)), DONORS) # as the approved demo (2026-11-11)
add_condition <- \(d) {
  if (is.null(DS$conditions)) {
    d
  } else {
    left_join(
      d,
      select(SAMPLES, donor = sample, condition),
      by = "donor",
      relationship = "many-to-one",
      unmatched = c("error", "drop")
    )
  }
}
POP <- "CD8T" # the population where TNF differs most between versions (INDEX section 05)
# registry keys of the scorer run dirs and the experiment label each arm is shown under; base arms are taken from the first
SOURCES <- c(
  "Block A swaps" = "common-ablation",
  "Rebuild chain" = "common-chain-all",
  "Rebuild chain " = "common-chain7a8",
  "AF405 proof" = "common-af405-proof-score",
  "Passes grid" = "common-passes-grid"
)
BASE_ARMS <- c("as188_joint", "s102_c188", "as102_var")
OUT_DIR <- run_out("common-coherence-or")

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Read scorer tables ------------------------------------------------------------------------------
# one row per run, donor and arm; base arms kept once (from the first source)
read_run <- \(key) {
  inner_join(
    read_csv(
      file.path(run_in(key), "pct_pos.csv"),
      show_col_types = FALSE
    ),
    read_csv(
      file.path(run_in(key), "tnf_among_ifng.csv"),
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
  filter(pop == POP) |>
  add_condition()

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
    any_of(COND_KEY),
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
  p_pos = sprintf("%s+ among %s+ (%%)", TNF_LABEL, IFNG_LABEL),
  p_neg = sprintf("%s+ among %s- (%%)", TNF_LABEL, IFNG_LABEL),
  log2_or = sprintf("%s-%s log2 odds ratio", TNF_LABEL, IFNG_LABEL)
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
    arm = factor(arm, levels = arm_order),
    donor = unname(PLOT_DONOR[donor])
  )
medians <- plot_frame |>
  summarise(
    value = median(value),
    .by = c(experiment, arm, readout, any_of(COND_KEY))
  )

# Figure ---------------------------------------------------------------------------------------------
p <- ggplot(plot_frame, aes(value, arm)) +
  geom_point(aes(colour = donor, shape = donor), size = 1.1, alpha = 0.8) +
  geom_point(data = medians, shape = 18, size = 2.2) +
  scale_colour_manual(values = categorical_palette(DONORS), name = "Donor") +
  scale_shape_manual(values = DONOR_SHAPES, name = "Donor") +
  facet_grid(
    as.formula(paste(
      paste(c("experiment", COND_KEY), collapse = " + "),
      "~ readout"
    )),
    scales = "free",
    space = "free_y"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = sprintf(
      "%s %s coherence per arm, common gates (diamond = donor median)",
      POP,
      TNF_LABEL
    )
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
