# Purpose: Nature-spec redraw (user, 2026-11-11; lucid-coder r-0008, r-0015, r-0016, r-0017) of
#          70_af405_absorption_test.R's figure: per-cell r = BV510 lost vs 1.8.8 / AF405 abundance, per donor,
#          for lost and kept TNF cells, in the AF405-only swap and its BV510 knock-out; grey band = BV510
#          coefficients of the 65 1.0.2 AF405 variants, dashed = their median. Same cells and clamp as 70
#          (its saved per-cell table); only the layout changes: theme_prism 7 pt Nimbus Sans, axes on every
#          panel, labels on the outer edges, donor colours from the shared label-keyed palette, legend right.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_events_absorption.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-10_bmv-af405-absorption/cells.csv.gz   r per lost/kept cell [70]
#          experiments/runs/2026-11-10_bmv-af405-proof/overlap_coef.csv    BV510 coefficient per variant [60]
# Outputs: experiments/runs/2026-11-11_bmv-restyle/10_af405_proof/absorption_fig_r_distribution.pdf, .png,
#            .csv   the .csv is the plotted per-cell frame

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- "experiments/runs/2026-11-11_bmv-restyle/10_af405_proof"
STEM <- "absorption_fig_r_distribution"
R_CLAMP <- c(-0.1, 0.35) # display clamp, as 70
DONORS <- c("A3", "A6", "B11", "B2", "D1")

stopifnot(
  "no output file exists yet" = !file.exists(file.path(
    OUT_DIR,
    paste0(STEM, ".png")
  ))
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Coefficient range of the 1.0.2 AF405 variants (as 70) ---------------------------------------------------
coef_102 <- read_csv(
  "experiments/runs/2026-11-10_bmv-af405-proof/overlap_coef.csv",
  show_col_types = FALSE
) |>
  filter(version == "1.0.2 (as used in 1.8.8: v - s188)", target == "BV510") |>
  pull(coef)
stopifnot("65 AF405 variants" = length(coef_102) == 65)

# Plot frame: r per cell in the swap and the knock-out, clamped for display (as 70) ----------------------
plot_frame <- read_csv(
  "experiments/runs/2026-11-10_bmv-af405-absorption/cells.csv.gz",
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
      labels = c("lost (TNF+ in\n1.8.8 only)", "kept (TNF+\nin both)")
    ),
    donor = factor(donor, levels = DONORS),
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
    x = "r = BV510 lost vs 1.8.8 / AF405 abundance (clamped to -0.1 to 0.35)",
    y = "Density",
    title = "Grey band: BV510 coefficients of the 65 1.0.2 AF405 variants (dashed = median)"
  ) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))
save_panel(p, plot_frame, OUT_DIR, STEM, 183, 150)
