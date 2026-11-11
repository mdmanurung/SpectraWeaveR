# Purpose: describe the 6 opt188 arms (user, 2026-11-10, grill Q7/Q12). opt188 = pregate, joint, 1/1 by
#          design; the other 5 arms are sensitivity arms. No arm is picked by outcome. CD8T readouts, each
#          donor a point with the donor median: TNF squeezed share, TNF share above arcsinh 1, TNF+ among
#          IFNg+, TNF+ among IFNg-, TNF-IFNg log2 odds ratio, negative rSD ratio vs as188_joint for TNF,
#          IFNg, Th2 and IL-10.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/67_opt188_report.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-opt188/scores.csv          squeezed and pos share [51_score_ablation.R]
#          experiments/runs/2026-11-10_bmv-opt188/pct_pos.csv         % TNF+, % IFNg+ [51]
#          experiments/runs/2026-11-10_bmv-opt188/tnf_among_ifng.csv  % TNF+ among IFNg+ [51]
#          experiments/runs/2026-11-10_bmv-opt188/neg_spread.csv      negative rSD ratios [51]
# Outputs: experiments/runs/2026-11-10_bmv-opt188/report/opt188_readouts.csv  per donor, arm, readout
#          experiments/runs/2026-11-10_bmv-opt188/report/fig_opt188_readouts.png
#          experiments/runs/2026-11-10_bmv-opt188/report/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- "experiments/runs/2026-11-10_bmv-opt188"
POP <- "CD8T" # grill Q12
DESIGN_ARM <- "opt188_pregate_joint_p1af1" # opt188 by design (grill Q7)
OUT_DIR <- file.path(RUN_DIR, "report")

stopifnot("report dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)
read_run <- \(f) {
  read_csv(file.path(RUN_DIR, f), show_col_types = FALSE) |> filter(pop == POP)
}

# Readouts per donor and arm ---------------------------------------------------------------------
# squeezed and positive share of TNF from the scorer
tnf_shares <- read_run("scores.csv") |>
  filter(cytokine == "TNF") |>
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
  pivot_longer(-c(donor, arm), names_to = "readout", values_to = "value")
write_csv(readouts, file.path(OUT_DIR, "opt188_readouts.csv"))

# Plot frame ---------------------------------------------------------------------------------------
readout_labels <- c(
  squeezed_share = "TNF squeezed near 0 (share)",
  pos_share = "TNF above arcsinh 1 (share)",
  pct_among_ifng_pos = "TNF+ among IFNg+ (%)",
  pct_among_ifng_neg = "TNF+ among IFNg- (%)",
  log2_or = "TNF-IFNg log2 odds ratio",
  rsd_ratio_TNF = "TNF neg. rSD / as188",
  rsd_ratio_IFNg = "IFNg neg. rSD / as188",
  rsd_ratio_Th2 = "Th2 neg. rSD / as188",
  `rsd_ratio_IL-10` = "IL-10 neg. rSD / as188"
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
    design = arm == DESIGN_ARM
  )
stopifnot("every readout has a label" = !anyNA(plot_frame$readout))
medians <- plot_frame |>
  summarise(value = median(value), .by = c(arm, readout, design))

# Figure ---------------------------------------------------------------------------------------------
p <- ggplot(plot_frame, aes(value, arm)) +
  geom_point(aes(colour = donor), size = 1.3, alpha = 0.8) +
  geom_point(data = medians, aes(shape = design), size = 3) +
  scale_shape_manual(
    values = c(`FALSE` = 18, `TRUE` = 15),
    labels = c("median", "median, opt188 by design"),
    name = NULL
  ) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  facet_wrap(~readout, scales = "free_x", ncol = 3) +
  labs(
    x = NULL,
    y = NULL,
    title = "CD8T readouts of the opt188 arms (descriptive; no arm picked by outcome)"
  ) +
  theme_bw(base_size = 9)
ggsave(
  file.path(OUT_DIR, "fig_opt188_readouts.png"),
  p,
  width = 11,
  height = 8,
  dpi = 150,
  bg = "white"
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
