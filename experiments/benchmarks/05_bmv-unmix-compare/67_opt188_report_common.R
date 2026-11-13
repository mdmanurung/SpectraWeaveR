# Purpose: describe the 6 opt188 arms (user, 2026-11-10, grill Q7/Q12). opt188 = pregate, joint, 1/1 by
#          design; the other 5 arms are sensitivity arms. No arm is picked by outcome. CD8T readouts, each
#          donor a point with the donor median: TNF squeezed share, TNF share above arcsinh 1, TNF+ among
#          IFNg+, TNF+ among IFNg-, TNF-IFNg log2 odds ratio, negative rSD ratio vs as188_joint for TNF,
#          IFNg, Th2 and IL-10.
#          Common-gate copy (plans/common_gates_plan.md, batch 2): reads the 51c opt188 rescore at cut 1.5.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/67_opt188_report_common.R   (repo root)
#          Needs the TNF and IFNg readouts: stops when the dataset's targets are not TNF and IFNg (EXV).
# Inputs:  run_in("common-opt188") is the 51c run folder of the opt188 arms (registry key common-opt188):
#          <run>/scores.csv          squeezed and pos share [51c_score_ablation_common.R]
#          <run>/pct_pos.csv         % TNF+, % IFNg+ [51c]
#          <run>/tnf_among_ifng.csv  % TNF+ among IFNg+ [51c]
#          <run>/neg_spread.csv      negative rSD ratios [51c]
#          datasets/<SW_DATASET>.json (targets, pos_rule), <ds>_samples.csv
# Outputs: <run>/report/opt188_readouts.csv  per sample, arm, readout (24H: + condition)
#          <run>/report/fig_opt188_readouts.png
#          <run>/report/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- run_in("common-opt188") # the 51c run folder; the report goes inside it
POP <- "CD8T" # grill Q12
DESIGN_ARM <- "opt188_pregate_joint_p1af1" # opt188 by design (grill Q7)
OUT_DIR <- file.path(RUN_DIR, "report")
TNF_LABEL <- cfg("targets")$tnf_label
IFNG_LABEL <- cfg("targets")$ifng_label
if (!identical(c(TNF_LABEL, IFNG_LABEL), c("TNF", "IFNg"))) {
  stop(
    "67_opt188_report_common.R: needs TNF and IFNg, missing in dataset ",
    DS$dataset,
    " (the two dye slots are ",
    TNF_LABEL,
    " and ",
    IFNG_LABEL,
    ")",
    call. = FALSE
  )
}
# label of the positive cut in the figure: the fixed value, or the per-donor rule
POS_LABEL <- if (identical(cfg("pos_rule")$type, "fixed")) {
  sprintf("arcsinh %s", cfg("pos_rule")$value)
} else {
  "the positive cut"
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

stopifnot("report dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)
read_run <- \(f) {
  read_csv(file.path(RUN_DIR, f), show_col_types = FALSE) |> filter(pop == POP)
}

# Readouts per donor and arm ---------------------------------------------------------------------
# squeezed and positive share of TNF from the scorer
tnf_shares <- read_run("scores.csv") |>
  filter(cytokine == TNF_LABEL) |>
  select(donor, arm, squeezed_share, pos_share)

# TNF+ among IFNg+ and IFNg- (derived), and the odds ratio
coherence <- inner_join(
  read_run("pct_pos.csv"),
  read_run("tnf_among_ifng.csv"),
  by = join_by(donor, arm, pop),
  relationship = "one-to-one",
  unmatched = "error"
) |>
  mutate(
    p_pos = pct_tnf_among_ifng / 100,
    p_neg = (pct_tnf / 100 - p_pos * pct_ifng / 100) / (1 - pct_ifng / 100),
    log2_or = log2((p_pos / (1 - p_pos)) / (p_neg / (1 - p_neg))),
    pct_among_ifng_pos = 100 * p_pos,
    pct_among_ifng_neg = 100 * p_neg
  ) |>
  select(donor, arm, pct_among_ifng_pos, pct_among_ifng_neg, log2_or)

# negative spread vs as188_joint, one column per cytokine
spread <- read_run("neg_spread.csv") |>
  select(donor, arm, cytokine, rsd_ratio) |>
  pivot_wider(
    names_from = cytokine,
    values_from = rsd_ratio,
    names_prefix = "rsd_ratio_"
  )

readouts <- tnf_shares |>
  inner_join(
    coherence,
    by = join_by(donor, arm),
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  inner_join(
    spread,
    by = join_by(donor, arm),
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  pivot_longer(-c(donor, arm), names_to = "readout", values_to = "value") |>
  add_condition()
write_csv(readouts, file.path(OUT_DIR, "opt188_readouts.csv"))

# Plot frame ---------------------------------------------------------------------------------------
RSD_MARKERS <- c(TNF_LABEL, IFNG_LABEL, "Th2", "IL-10")
readout_labels <- c(
  squeezed_share = sprintf("%s squeezed near 0 (share)", TNF_LABEL),
  pos_share = sprintf("%s above %s (share)", TNF_LABEL, POS_LABEL),
  pct_among_ifng_pos = sprintf("%s+ among %s+ (%%)", TNF_LABEL, IFNG_LABEL),
  pct_among_ifng_neg = sprintf("%s+ among %s- (%%)", TNF_LABEL, IFNG_LABEL),
  log2_or = sprintf("%s-%s log2 odds ratio", TNF_LABEL, IFNG_LABEL),
  set_names(
    sprintf("%s neg. rSD / as188", RSD_MARKERS),
    paste0("rsd_ratio_", RSD_MARKERS)
  )
)
arm_levels <- c(
  "as188_joint",
  "opt188_pregate_joint_p1af1",
  "opt188_pregate_joint_p5af5",
  "opt188_pregate_legacy",
  "opt188_legacyspec_joint_p1af1",
  "opt188_legacyspec_joint_p5af5",
  "opt188_legacyspec_legacy",
  "as102_var",
  "s102_c188"
)
plot_frame <- readouts |>
  filter(arm %in% arm_levels) |>
  mutate(
    readout = factor(readout_labels[readout], levels = unname(readout_labels)),
    arm = factor(arm, levels = rev(arm_levels)),
    design = arm == DESIGN_ARM,
    donor = unname(PLOT_DONOR[donor])
  )
stopifnot("every readout has a label" = !anyNA(plot_frame$readout))
medians <- plot_frame |>
  summarise(
    value = median(value),
    .by = c(arm, readout, design, any_of(COND_KEY))
  )

# Figure ---------------------------------------------------------------------------------------------
p <- ggplot(plot_frame, aes(value, arm)) +
  geom_point(aes(colour = donor, shape = donor), size = 1.1, alpha = 0.8) +
  geom_point(data = medians, aes(shape = design), size = 2.2) +
  scale_colour_manual(values = categorical_palette(DONORS), name = "Donor") +
  scale_shape_manual(
    values = c(DONOR_SHAPES, `FALSE` = 18, `TRUE` = 15),
    breaks = c(DONORS, "FALSE", "TRUE"),
    labels = c(DONORS, "median", "median, opt188 by design"),
    name = "Donor"
  ) +
  facet_wrap(c("readout", COND_KEY), scales = "free_x", ncol = 3) +
  labs(
    x = NULL,
    y = NULL,
    title = sprintf(
      "%s readouts of the opt188 arms, common gates (descriptive; no arm picked by outcome)",
      POP
    )
  ) +
  value_grid("x")
save_panel(
  p,
  plot_frame,
  OUT_DIR,
  "fig_opt188_readouts",
  width_mm = 183,
  height_mm = 150
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
