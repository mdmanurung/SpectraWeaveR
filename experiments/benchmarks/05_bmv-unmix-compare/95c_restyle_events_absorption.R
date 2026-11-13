# Purpose: Nature-spec redraw on the common gates (plans/common_gates_plan.md G11) of 70_af405_absorption_test.R's
#          figure, from its common-gate rerun: per-cell r = BV510 lost vs 1.8.8 / AF405 abundance, per donor, for
#          lost and kept TNF cells, in the AF405-only swap and its BV510 knock-out; grey band = BV510 coefficients
#          of the 65 1.0.2 AF405 variants, dashed = median. Lost and kept cells use TNF+ at arcsinh 1.5 (G7a) in
#          the cells table. Copy of 95_restyle_events_absorption.R with new inputs.
#          The cut wording reads the dataset's positive-cut rule (CUT_TXT in 95c_common.R).
# Run:     SW_DATASET=3d SW_RUN_DATE=YYYY-MM-DD pixi run Rscript
#            experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_events_absorption.R from the repo root.
# Inputs:  run_in("af405-absorption-common")/cells.csv.gz   r per lost/kept cell [70c]
#          run_in("af405-proof")/overlap_coef.csv   BV510 coefficient per variant [60];
#            a property of the spectra, not of the gates, so inherited unchanged
# Outputs: run_out("restyle-common") (shared by the 95c scripts)/10_af405_proof/absorption_fig_r_distribution.{pdf,png,csv}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")
SCRIPT <- "95c_restyle_events_absorption"
TNF_LAB <- cfg("targets")$tnf_label
TNF_DYE <- cfg("targets")$tnf_dye

# Parameters ---------------------------------------------------------------------------------------
STEM <- "absorption_fig_r_distribution"
R_CLAMP <- c(-0.1, 0.35) # display clamp, as 70

OUT_DIR <- fresh_dir("10_af405_proof")

# Coefficient range of the 1.0.2 AF405 variants (as 70) ---------------------------------------------------
coef_102 <- read_csv(
  file.path(run_in("af405-proof"), "overlap_coef.csv"),
  show_col_types = FALSE
) |>
  filter(version == "1.0.2 (as used in 1.8.8: v - s188)", target == TNF_DYE) |>
  pull(coef)
stopifnot("1.0.2 AF405 variants found" = length(coef_102) > 0)

# Plot frame: r per cell in the swap and the knock-out, clamped for display (as 70) ----------------------
plot_frame <- read_csv(
  file.path(run_in("af405-absorption-common"), "cells.csv.gz"),
  show_col_types = FALSE
) |>
  pivot_longer(c(r_swap, r_ko), names_to = "arm", values_to = "r") |>
  mutate(
    arm = factor(
      arm,
      levels = c("r_swap", "r_ko"),
      labels = c("1.0.2 AF405\nvariants", "same, BV510\npart removed")
    ),
    group = factor(
      group,
      levels = c("lost", "kept"),
      labels = c(
        paste0("lost (", TNF_LAB, "+ in\n1.8.8 only)"),
        paste0("kept (", TNF_LAB, "+\nin both)")
      )
    ),
    donor = forcats::fct(donor, levels = DONORS),
    r = pmin(pmax(r, R_CLAMP[1]), R_CLAMP[2])
  )

# Figure ---------------------------------------------------------------------------------------------------
p <- ggplot(plot_frame, aes(r, colour = donor)) +
  annotate(
    "rect",
    xmin = min(coef_102),
    xmax = max(coef_102),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.15
  ) +
  geom_vline(
    xintercept = median(coef_102),
    linetype = "dashed",
    linewidth = FIG_LINE / 2
  ) +
  geom_density(linewidth = FIG_LINE, key_glyph = "path") +
  scale_colour_manual(values = categorical_palette(DONORS), name = "Donor") +
  facet_grid(
    group + arm ~ pop,
    scales = "free_y",
    axes = "all",
    axis.labels = "margins"
  ) +
  labs(
    x = paste0(
      "r = ",
      TNF_DYE,
      " lost vs 1.8.8 / AF405 abundance (clamped to -0.1 to 0.35)"
    ),
    y = "Density",
    title = paste0(
      TNF_LAB,
      "+ = above ",
      CUT_TXT,
      ". Grey band: ",
      TNF_DYE,
      " coefficients of the ",
      length(coef_102),
      " 1.0.2 AF405 variants (dashed = median)"
    )
  ) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))
save_panel(p, plot_frame, OUT_DIR, STEM, 183, 150)
add_manifest(
  "95c_restyle_events_absorption",
  file.path("2026-11-11_bmv-restyle/10_af405_proof", paste0(STEM, ".png")),
  file.path(OUT_DIR, paste0(STEM, ".png")),
  sub(
    "^experiments/runs/",
    "",
    file.path(run_in("af405-absorption-common"), "cells.csv.gz")
  )
)
