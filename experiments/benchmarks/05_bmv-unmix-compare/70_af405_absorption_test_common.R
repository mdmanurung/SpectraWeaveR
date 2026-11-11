# Purpose: Common-gates copy (plans/common_gates_plan.md, batch 3): gates are the common SpectroFlo labels and the
#          cytokine cut is arcsinh 1.5 (G7a). The predictions P1-P4 and P2_SHARE are unchanged.
#          Tests the proposed mechanism for the TNF zeroing (user, 2026-11-10): for a dim-TNF cell the per-cell
#          variant search picks an Alexa Fluor 405 (GLUT1) variant whose BV510-shaped part explains the cell's
#          BV510 signal, so TNF fits to ~0. In the AF405-only swap (var102_only_Alexa_Fluor_405) only the AF405
#          variants differ from 1.8.8 (as188_joint), so the BV510 a cell loses should equal its AF405
#          abundance x the chosen variant's BV510 coefficient:  r = (BV510_188 - BV510_swap) / AF405_swap.
# Predictions (written before the run, 2026-11-10):
#   P1  cells TNF+ in 1.8.8 and TNF- in the swap ("lost"): median r lies inside the range of the 65 1.0.2 AF405
#       variants' BV510 coefficients (overlap_coef.csv, "as used in 1.8.8"), in every donor and population.
#   P2  >= 80% of lost cells have r inside that range (pooled, per population).
#   P3  the same cells in the knock-out (af405proof_ko_bv510, BV510 part removed): median |r| < the smallest
#       coefficient, in every donor and population.
#   P4  cells TNF+ in both 1.8.8 and the swap ("kept"): median r below the smallest coefficient.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/70_af405_absorption_test_common.R   (repo root)
# Inputs:  experiments/runs/2026-11-11_bmv-common-gated/gs/<arm>/  common-gate GatingSets, rows aligned across arms [batch 0, 97d]
#          experiments/runs/2026-11-10_bmv-af405-proof/overlap_coef.csv  BV510 coefficient per AF405 variant [60]
# Outputs: experiments/runs/2026-11-11_bmv-af405-absorption-common/cells.csv.gz     per lost/kept cell: r in swap and KO
#          experiments/runs/2026-11-11_bmv-af405-absorption-common/predictions.csv  P1-P4 per donor/pop, pass/fail
#          experiments/runs/2026-11-11_bmv-af405-absorption-common/fig_r_distribution.{pdf,png,csv}  r per group, band = coef range
#          experiments/runs/2026-11-11_bmv-af405-absorption-common/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters ---------------------------------------------------------------------------------------
GATED <- "experiments/runs/2026-11-11_bmv-common-gated/gs" # batch 0: SpectroFlo labels as gates
REF_ARM <- "as188_joint"
SWAP_ARM <- "var102_only_Alexa_Fluor_405"
KO_ARM <- "af405proof_ko_bv510"
POPS <- c("CD4T", "CD8T")
CUT <- 1.5 # fixed cytokine cut, arcsinh (G7a, user 2026-11-11; was 1)
COFACTOR <- 6000 # BV510 is stored as arcsinh(x / 6000) in the GatingSets (19a)
P2_SHARE <- 0.8 # pre-registered above
OUT_DIR <- "experiments/runs/2026-11-11_bmv-af405-absorption-common"

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()

# Coefficient range of the 1.0.2 AF405 variants (BV510 part, as the 1.8.8 code uses them) -------------
coef_102 <- read_csv(
  "experiments/runs/2026-11-10_bmv-af405-proof/overlap_coef.csv",
  show_col_types = FALSE
) |>
  filter(version == "1.0.2 (as used in 1.8.8: v - s188)", target == "BV510") |>
  pull(coef)
coef_range <- range(coef_102)
stopifnot("65 AF405 variants" = length(coef_102) == 65)

# Per-cell BV510 and AF405 abundance, same rows in every arm ---------------------------------------
gating_sets <- c(REF_ARM, SWAP_ARM, KO_ARM) |>
  set_names() |>
  map(\(a) load_gs(file.path(GATED, a)))
cells <- sampleNames(gating_sets[[REF_ARM]]) |>
  set_names() |>
  map(\(s) {
    ev <- map(gating_sets, \(g) exprs(gh_pop_get_data(g[[s]], "root")))
    stopifnot(
      "rows align across arms" = all(map_int(ev, nrow) == nrow(ev[[REF_ARM]]))
    )
    POPS |>
      set_names() |>
      map(\(pop) {
        rows <- gh_pop_get_indices(gating_sets[[REF_ARM]][[s]], pop)
        tibble(
          tnf_188 = ev[[REF_ARM]][rows, "BV510-A"],
          tnf_swap = ev[[SWAP_ARM]][rows, "BV510-A"],
          tnf_ko = ev[[KO_ARM]][rows, "BV510-A"],
          af405_swap = ev[[SWAP_ARM]][rows, "Alexa Fluor 405-A"],
          af405_ko = ev[[KO_ARM]][rows, "Alexa Fluor 405-A"]
        )
      }) |>
      list_rbind(names_to = "pop")
  }) |>
  list_rbind(names_to = "donor")

# Groups and the absorption ratio r ------------------------------------------------------------------
# lost: TNF+ in 1.8.8, TNF- in the swap; kept: TNF+ in both
tested <- cells |>
  mutate(
    group = case_when(
      tnf_188 > CUT & tnf_swap < CUT ~ "lost",
      tnf_188 > CUT & tnf_swap > CUT ~ "kept",
      .default = NA_character_
    ),
    lin_188 = COFACTOR * sinh(tnf_188),
    r_swap = (lin_188 - COFACTOR * sinh(tnf_swap)) / af405_swap,
    r_ko = (lin_188 - COFACTOR * sinh(tnf_ko)) / af405_ko
  ) |>
  filter(!is.na(group), af405_swap > 0, af405_ko > 0)
write_csv(
  select(tested, donor, pop, group, tnf_188, tnf_swap, tnf_ko, r_swap, r_ko),
  file.path(OUT_DIR, "cells.csv.gz")
)

# Predictions -----------------------------------------------------------------------------------------
in_range <- \(r) r >= coef_range[1] & r <= coef_range[2]
per_donor <- tested |>
  summarise(
    n = n(),
    median_r_swap = median(r_swap),
    median_abs_r_ko = median(abs(r_ko)),
    share_in_range = mean(in_range(r_swap)),
    .by = c(donor, pop, group)
  )
predictions <- bind_rows(
  per_donor |>
    filter(group == "lost") |>
    transmute(
      prediction = "P1 lost: median r in coef range",
      donor,
      pop,
      value = median_r_swap,
      pass = in_range(value)
    ),
  tested |>
    filter(group == "lost") |>
    summarise(value = mean(in_range(r_swap)), .by = pop) |>
    transmute(
      prediction = "P2 lost: share r in coef range >= 0.8",
      donor = "pooled",
      pop,
      value,
      pass = value >= P2_SHARE
    ),
  per_donor |>
    filter(group == "lost") |>
    transmute(
      prediction = "P3 lost, KO: median |r| < min coef",
      donor,
      pop,
      value = median_abs_r_ko,
      pass = value < coef_range[1]
    ),
  per_donor |>
    filter(group == "kept") |>
    transmute(
      prediction = "P4 kept: median r < min coef",
      donor,
      pop,
      value = median_r_swap,
      pass = value < coef_range[1]
    )
) |>
  mutate(coef_min = coef_range[1], coef_max = coef_range[2])
write_csv(predictions, file.path(OUT_DIR, "predictions.csv"))

# Figure: r per group, swap and KO, band = coefficient range of the 1.0.2 AF405 variants ---------------
plot_frame <- tested |>
  pivot_longer(c(r_swap, r_ko), names_to = "arm", values_to = "r") |>
  mutate(
    arm = factor(
      arm,
      levels = c("r_swap", "r_ko"),
      labels = c("1.0.2 AF405 variants", "same, BV510 part removed")
    ),
    group = factor(
      group,
      levels = c("lost", "kept"),
      labels = c("lost (TNF+ in 1.8.8 only)", "kept (TNF+ in both)")
    ),
    r = pmin(pmax(r, -0.1), 0.35)
  )
p <- ggplot(plot_frame, aes(r, colour = donor)) +
  annotate(
    "rect",
    xmin = coef_range[1],
    xmax = coef_range[2],
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.15
  ) +
  geom_vline(
    xintercept = median(coef_102),
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_density(linewidth = 0.3) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  facet_grid(group + arm ~ pop, scales = "free_y") +
  labs(
    x = "r = BV510 lost vs 1.8.8 / AF405 abundance (clamped to -0.1..0.35)",
    y = "Density",
    title = "Grey band: BV510 coefficients of the 65 1.0.2 AF405 variants (dashed = median)"
  )
# the density frame is the per-cell table (r clamped as plotted); curves are drawn by geom_density
save_panel(p, plot_frame, OUT_DIR, "fig_r_distribution", 183, 160)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
