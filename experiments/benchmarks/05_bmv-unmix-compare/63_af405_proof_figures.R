# Purpose: proof chain steps 2 to 4 (user, 2026-11-10): plot the intervention arms' CD8T and CD4T TNF (BV510)
#          squeezed share and share above arcsinh 1 per donor, and test the pre-registered predictions:
#          ko_bv510 <= 0.15 (CD8T squeezed, donor median); ki_bv510 >= 0.25; dose 0 -> 0.25 -> 0.5 -> 1
#          monotonic in the median. Scores come from 51_score_ablation.R (full gating, 39b metrics).
# Run:     sbatch experiments/slurm/run_af405_proof.sbatch 63_af405_proof_figures.R <score_run_dir>
# Inputs:  <score_run_dir>/scores.csv, gap_closed.csv   [51_score_ablation.R on the proof arms]
# Outputs: experiments/runs/2026-11-10_bmv-af405-proof/fig2_interventions.png  squeezed / above-1 per arm
#          experiments/runs/2026-11-10_bmv-af405-proof/fig2_interventions_source.csv
#          experiments/runs/2026-11-10_bmv-af405-proof/fig4_dose.png             dose-response, donors as lines
#          experiments/runs/2026-11-10_bmv-af405-proof/fig4_dose_source.csv
#          experiments/runs/2026-11-10_bmv-af405-proof/predictions_interventions.csv  pass / fail
#          experiments/runs/2026-11-10_bmv-af405-proof/sessionInfo_63.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
RUN_DIR <- "experiments/runs/2026-11-10_bmv-af405-proof"
ARM_LABELS <- c(
  af405proof_base188 = "1.8.8 (all 1.8.8 variants)",
  af405proof_ki_bv510 = "knock-in: 1.8.8 AF405 + BV510 part",
  af405proof_dose_000 = "dose 0: AF405 = 1.0.2 reference only",
  af405proof_dose_025 = "dose 0.25",
  af405proof_dose_050 = "dose 0.5",
  var102_only_Alexa_Fluor_405 = "dose 1: 1.0.2 AF405 variants",
  af405proof_recentre = "1.0.2 spread around 1.8.8 reference",
  af405proof_ko_bv510 = "knock-out: 1.0.2 AF405 minus BV510 part",
  s102_c188 = "all 1.0.2 inputs (target)"
) # figure order
DOSE <- c(
  af405proof_dose_000 = 0,
  af405proof_dose_025 = 0.25,
  af405proof_dose_050 = 0.5,
  var102_only_Alexa_Fluor_405 = 1
)
KO_MAX <- 0.15 # pre-registered (user-approved plan, 2026-11-10)
KI_MIN <- 0.25 # pre-registered

score_dir <- commandArgs(trailingOnly = TRUE)[[1]]
scores <- read_csv(
  file.path(score_dir, "scores.csv"),
  show_col_types = FALSE
) |>
  filter(cytokine == "TNF", arm %in% names(ARM_LABELS)) |>
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
      labels = c("Share squeezed near 0", "Share above arcsinh 1")
    )
  )
write_csv(plot_frame, file.path(RUN_DIR, "fig2_interventions_source.csv"))
p2 <- ggplot(plot_frame, aes(value, arm_label)) +
  geom_line(
    aes(group = donor, colour = donor),
    orientation = "y",
    alpha = 0.35
  ) +
  geom_point(aes(colour = donor), size = 1.8) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 3,
    size = 4,
    colour = "black"
  ) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  scale_y_discrete(limits = rev) +
  facet_grid(pop ~ readout, scales = "free_x") +
  labs(x = "TNF (BV510)", y = NULL) +
  theme_bw(base_size = 11)
ggsave(
  file.path(RUN_DIR, "fig2_interventions.png"),
  p2,
  width = 11,
  height = 7,
  dpi = 130,
  bg = "white"
)

# Figure 4: dose-response, CD8T and CD4T squeezed share ----------------------------------------
dose_frame <- scores |>
  filter(arm %in% names(DOSE)) |>
  mutate(alpha = DOSE[arm])
write_csv(dose_frame, file.path(RUN_DIR, "fig4_dose_source.csv"))
base_line <- scores |>
  filter(arm == "af405proof_base188") |>
  summarise(base = median(squeezed_share), .by = pop)
p4 <- ggplot(dose_frame, aes(alpha, squeezed_share)) +
  geom_hline(
    data = base_line,
    aes(yintercept = base),
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line(aes(group = donor, colour = donor)) +
  geom_point(aes(colour = donor), size = 1.8) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  facet_wrap(~pop) +
  labs(
    x = "Scale of the 1.0.2 AF405 spread (alpha)",
    y = "TNF share squeezed near 0"
  ) +
  theme_bw(base_size = 11)
ggsave(
  file.path(RUN_DIR, "fig4_dose.png"),
  p4,
  width = 9,
  height = 4,
  dpi = 130,
  bg = "white"
)

# Pre-registered predictions ----------------------------------------------------------------
med <- scores |>
  filter(pop == "CD8T") |>
  summarise(median_squeezed = median(squeezed_share), .by = arm)
med_of <- \(a) med$median_squeezed[med$arm == a]
dose_medians <- map_dbl(names(DOSE), med_of)
predictions <- tibble(
  prediction = c(
    "knock-out CD8T squeezed median <= 0.15",
    "knock-in CD8T squeezed median >= 0.25",
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
write_csv(predictions, file.path(RUN_DIR, "predictions_interventions.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_63.txt")
)
