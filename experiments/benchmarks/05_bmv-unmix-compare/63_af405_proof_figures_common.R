# Purpose: proof chain steps 2 to 4 (user, 2026-11-10): plot the intervention arms' CD8T and CD4T TNF (BV510)
#          squeezed share and share above arcsinh 1 per donor, and test the pre-registered predictions:
#          ko_bv510 <= 0.15 (CD8T squeezed, donor median); ki_bv510 >= 0.25; dose 0 -> 0.25 -> 0.5 -> 1
#          monotonic in the median. Scores come from 51c_score_ablation_common.R (full gating, 39b metrics).
#          Common-gate copy (plans/common_gates_plan.md, batch 2): writes a NEW dir, never the old proof dir.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/63_af405_proof_figures_common.R <score_run_dir> <out_dir>
#            <out_dir> must not exist; it is run_out("common-af405-proof") of the caller
#            (3D: experiments/runs/2026-11-11_bmv-common-af405-proof)
# Inputs:  <score_run_dir>/scores.csv, gap_closed.csv   [51c_score_ablation_common.R on the proof arms;
#            registry key common-af405-proof-score]
#          datasets/<SW_DATASET>.json (targets, pos_rule), <ds>_samples.csv
# Outputs: <out_dir>/fig2_interventions.png  squeezed / above-1 per arm
#          <out_dir>/fig2_interventions_source.csv
#          <out_dir>/fig4_dose.png             dose-response, donors as lines
#          <out_dir>/fig4_dose_source.csv
#          <out_dir>/predictions_interventions.csv  pass / fail
#          <out_dir>/sessionInfo_63.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))

# Parameters ------------------------------------------------------------------------------
TNF_LABEL <- cfg("targets")$tnf_label
TNF_DYE <- cfg("targets")$tnf_dye
# label of the positive cut in the figure: the fixed value, or the per-donor rule
POS_LABEL <- if (identical(cfg("pos_rule")$type, "fixed")) {
  sprintf("arcsinh %s", cfg("pos_rule")$value)
} else {
  "the positive cut"
}
ARM_LABELS <- c(
  af405proof_base188 = "1.8.8 (all 1.8.8 variants)",
  af405proof_ki_bv510 = sprintf("knock-in: 1.8.8 AF405 + %s part", TNF_DYE),
  af405proof_dose_000 = "dose 0: AF405 = 1.0.2 reference only",
  af405proof_dose_025 = "dose 0.25",
  af405proof_dose_050 = "dose 0.5",
  var102_only_Alexa_Fluor_405 = "dose 1: 1.0.2 AF405 variants",
  af405proof_recentre = "1.0.2 spread around 1.8.8 reference",
  af405proof_ko_bv510 = sprintf(
    "knock-out: 1.0.2 AF405 minus %s part",
    TNF_DYE
  ),
  s102_c188 = "all 1.0.2 inputs (target)"
) # figure order
DOSE <- c(
  af405proof_dose_000 = 0,
  af405proof_dose_025 = 0.25,
  af405proof_dose_050 = 0.5,
  var102_only_Alexa_Fluor_405 = 1
)
KO_MAX <- 0.15 # pre-registered (user-approved plan, 2026-11-10); inherited, set on the old squeeze cut and cut 1
KI_MIN <- 0.25 # pre-registered; inherited like KO_MAX, not refitted on the common gates

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

score_dir <- commandArgs(trailingOnly = TRUE)[[1]]
RUN_DIR <- commandArgs(trailingOnly = TRUE)[[2]]
stopifnot("out dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
scores <- read_csv(
  file.path(score_dir, "scores.csv"),
  show_col_types = FALSE
) |>
  filter(cytokine == TNF_LABEL, arm %in% names(ARM_LABELS)) |>
  add_condition() |>
  mutate(arm_label = factor(ARM_LABELS[arm], levels = ARM_LABELS))
stopifnot(
  "every proof arm scored" = setequal(unique(scores$arm), names(ARM_LABELS))
)

# Figure 2: squeezed share and share above 1 per arm, donors as points ------------------------
plot_frame <- scores |>
  select(donor, arm, arm_label, pop, squeezed_share, pos_share) |>
  pivot_longer(
    c(squeezed_share, pos_share),
    names_to = "readout",
    values_to = "value"
  ) |>
  mutate(
    readout = factor(
      readout,
      levels = c("squeezed_share", "pos_share"),
      labels = c("Share squeezed near 0", paste("Share above", POS_LABEL))
    ),
    donor = unname(PLOT_DONOR[donor])
  )
p2 <- ggplot(plot_frame, aes(value, arm_label)) +
  geom_line(
    aes(group = donor, colour = donor),
    orientation = "y",
    alpha = 0.35,
    linewidth = FIG_LINE
  ) +
  geom_point(aes(colour = donor, shape = donor), size = 1.3) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 3,
    size = 2.5,
    colour = "black"
  ) +
  scale_colour_manual(values = categorical_palette(DONORS), name = "Donor") +
  scale_shape_manual(values = DONOR_SHAPES, name = "Donor") +
  scale_y_discrete(limits = rev) +
  facet_grid(
    as.formula(paste(paste(c("pop", COND_KEY), collapse = " + "), "~ readout")),
    scales = "free_x"
  ) +
  labs(x = sprintf("%s (%s)", TNF_LABEL, TNF_DYE), y = NULL) +
  value_grid("x")
save_panel(
  p2,
  plot_frame,
  RUN_DIR,
  "fig2_interventions",
  width_mm = 183,
  height_mm = 110
)

# Figure 4: dose-response, CD8T and CD4T squeezed share ----------------------------------------
dose_frame <- scores |>
  filter(arm %in% names(DOSE)) |>
  mutate(alpha = DOSE[arm], donor = unname(PLOT_DONOR[donor]))
base_line <- scores |>
  filter(arm == "af405proof_base188") |>
  summarise(base = median(squeezed_share), .by = c(pop, any_of(COND_KEY)))
p4 <- ggplot(dose_frame, aes(alpha, squeezed_share)) +
  geom_hline(
    data = base_line,
    aes(yintercept = base),
    linetype = "dashed",
    colour = "grey50",
    linewidth = FIG_LINE / 2
  ) +
  geom_line(aes(group = donor, colour = donor), linewidth = FIG_LINE) +
  geom_point(aes(colour = donor, shape = donor), size = 1.3) +
  scale_colour_manual(values = categorical_palette(DONORS), name = "Donor") +
  scale_shape_manual(values = DONOR_SHAPES, name = "Donor") +
  facet_wrap(c("pop", COND_KEY)) +
  labs(
    x = "Scale of the 1.0.2 AF405 spread (alpha)",
    y = sprintf("%s share squeezed near 0", TNF_LABEL)
  )
save_panel(p4, dose_frame, RUN_DIR, "fig4_dose", width_mm = 89, height_mm = 70)

# Pre-registered predictions ----------------------------------------------------------------
med <- scores |>
  filter(pop == "CD8T") |>
  summarise(
    median_squeezed = median(squeezed_share),
    .by = c(arm, any_of(COND_KEY))
  )
# the three predictions on the medians of one table (24H: one table per condition)
predict_from <- \(med) {
  med_of <- \(a) med$median_squeezed[med$arm == a]
  dose_medians <- map_dbl(names(DOSE), med_of)
  tibble(
    prediction = c(
      sprintf("knock-out CD8T squeezed median <= %s", KO_MAX),
      sprintf("knock-in CD8T squeezed median >= %s", KI_MIN),
      "dose 0, 0.25, 0.5, 1 medians increase monotonically"
    ),
    value = c(
      sprintf("%.3f", med_of("af405proof_ko_bv510")),
      sprintf("%.3f", med_of("af405proof_ki_bv510")),
      paste(sprintf("%.3f", dose_medians), collapse = " -> ")
    ),
    pass = c(
      med_of("af405proof_ko_bv510") <= KO_MAX,
      med_of("af405proof_ki_bv510") >= KI_MIN,
      all(diff(dose_medians) > 0)
    )
  )
}
predictions <- if (is.null(DS$conditions)) {
  predict_from(med)
} else {
  split(med, med$condition) |>
    map(predict_from) |>
    list_rbind(names_to = "condition")
}
write_csv(predictions, file.path(RUN_DIR, "predictions_interventions.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_63.txt")
)
