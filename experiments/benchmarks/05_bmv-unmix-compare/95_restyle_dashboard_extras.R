# Purpose: Nature-spec redraw (user, 2026-11-11; lucid-coder r-0008, r-0012, r-0014..r-0017) of the
#          dashboard figures that had no restyled version: BV510 density in the swap test, TNF+ among
#          IFNg+ cells, TNF/IFNg rank correlation with activation markers, Th2 density, the AF405
#          variant-subset gap closed, and the zip 1.8.8 spillover spread matrix. Reads saved tables only.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_dashboard_extras.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-10_bmv-swap-test/bv510_events.csv.gz   BV510 per cell [37]
#          experiments/runs/2026-11-10_bmv-tnf-coherence/fig_a_source.csv   TNF+ among IFNg+/- [46]
#          experiments/runs/2026-11-10_bmv-tnf-coherence/fig_b_source.csv   Spearman rho [46]
#          experiments/runs/2026-11-10_bmv-other-cytokines/fig_b_source.csv cytokine densities [47]
#          experiments/runs/2026-11-10_bmv-af405-mechanism/scores/fig_gap_closed_source.csv  [55 via 51]
#          experiments/runs/2026-11-10_bmv-ssm-188/ssm_long.csv             zip SSM, long [81]
# Outputs: experiments/runs/2026-11-11_bmv-restyle/08_dashboard_extras/<name>.pdf, .png, .csv
#            fig_bv510_density, fig_coh_tnf_among_ifng, fig_coh_rho_activation, fig_oc_density_th2,
#            proof_subsets_gap_closed, fig_ssm_188_zip

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/report_figures.R") # theme, save_panel, donor scales
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS, VERSION_ROWS

# Parameters ---------------------------------------------------------------------------------------
RUNS <- "experiments/runs"
OUT_DIR <- "experiments/runs/2026-11-11_bmv-restyle/08_dashboard_extras"
CUT <- 1 # fixed cytokine cut, arcsinh (user, 2026-11-10)
DENSITY_WINDOW <- c(-1.5, 4) # display window, arcsinh; the density uses every event

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE)
arm_levels_all <- unlist(VERSION_ROWS, use.names = FALSE)

# BV510 density in the swap test --------------------------------------------------------------------
# one density per arm, donor and population; arms keyed by the categorical palette
bv510 <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-swap-test/bv510_events.csv.gz"),
  show_col_types = FALSE
) |>
  filter(between(.data$arcsinh, DENSITY_WINDOW[1], DENSITY_WINDOW[2]))
swap_arms <- c(
  "as188_joint",
  "swap_af102",
  "swap_var102",
  "s188_c102",
  "s102_c188"
)
stopifnot("swap-test arms as expected" = setequal(unique(bv510$arm), swap_arms))
bv510_density <- bv510 |>
  reframe(
    as_tibble(density(
      .data$arcsinh,
      from = DENSITY_WINDOW[1],
      to = DENSITY_WINDOW[2],
      n = 256
    )[c("x", "y")]),
    .by = c(arm, sample, pop)
  ) |>
  mutate(arm = factor(.data$arm, levels = swap_arms))
p_bv510 <- ggplot(bv510_density, aes(x, y, colour = arm)) +
  geom_vline(
    xintercept = CUT,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.3
  ) +
  geom_line(linewidth = 0.35) +
  scale_colour_manual(values = categorical_palette(swap_arms), name = "Arm") +
  facet_grid(pop ~ sample, axes = "all", axis.labels = "margins") +
  labs(
    x = "TNF, BV510-A (arcsinh)",
    y = "Density",
    title = "BV510 (TNF) density, swap test"
  )
save_panel(p_bv510, bv510_density, OUT_DIR, "fig_bv510_density", 183, 80)

# TNF+ among IFNg+ and IFNg- cells ------------------------------------------------------------------
coh <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-tnf-coherence/fig_a_source.csv"),
  show_col_types = FALSE
) |>
  pivot_longer(
    c(pct_tnf_of_ifng, pct_tnf_of_ifng_neg),
    names_to = "among",
    values_to = "pct"
  ) |>
  mutate(
    among = if_else(
      .data$among == "pct_tnf_of_ifng",
      "among IFNg+",
      "among IFNg-"
    )
  )
coh_out <- panel_arm_donor(
  coh,
  "pct",
  "% TNF+ (arcsinh cut 1)",
  "TNF+ among IFNg+ and IFNg- cells",
  intersect(arm_levels_all, unique(coh$arm)),
  rows = "among",
  cols = "parent"
)
save_panel(
  coh_out$plot,
  coh_out$data,
  OUT_DIR,
  "fig_coh_tnf_among_ifng",
  183,
  coh_out$height_mm
)

# TNF and IFNg rank correlation with activation markers ----------------------------------------------
# donor points small, donor median as a larger point, coloured by cytokine
rho <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-tnf-coherence/fig_b_source.csv"),
  show_col_types = FALSE
) |>
  mutate(
    arm = factor(
      .data$arm,
      levels = rev(intersect(arm_levels_all, unique(.data$arm)))
    )
  )
rho_median <- summarise(
  rho,
  rho = median(.data$rho),
  .by = c(arm, parent, cytokine, marker)
)
cytokine_pal <- c(TNF = "#D55E00", IFNg = "#0072B2")
p_rho <- ggplot(rho, aes(rho, arm, colour = cytokine)) +
  value_grid("x") +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_point(size = 0.4, alpha = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(
    data = rho_median,
    shape = 18,
    size = 1.6,
    position = position_dodge(width = 0.5)
  ) +
  scale_colour_manual(values = cytokine_pal, name = "Cytokine") +
  scale_y_discrete(labels = arm_label) +
  # two labelled ticks per ~25 mm panel so neighbouring labels do not touch; 0 is the solid line
  scale_x_continuous(breaks = c(-0.2, 0.2)) +
  facet_grid(parent ~ marker, axes = "all", axis.labels = "margins") +
  labs(
    x = "Spearman rho (diamond = donor median)",
    y = NULL,
    title = "Cytokine rank correlation with activation markers"
  )
save_panel(p_rho, rho, OUT_DIR, "fig_coh_rho_activation", 183, 150)

# Th2 density per arm --------------------------------------------------------------------------------
th2 <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-other-cytokines/fig_b_source.csv"),
  show_col_types = FALSE
) |>
  filter(.data$cytokine == "Th2")
versions <- unique(th2$version)
p_th2 <- ggplot(th2, aes(value, density, group = arm, colour = version)) +
  geom_vline(
    xintercept = CUT,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.3
  ) +
  geom_line(linewidth = 0.3) +
  scale_colour_manual(
    values = categorical_palette(versions),
    name = "Code version"
  ) +
  facet_grid(parent ~ donor, axes = "all", axis.labels = "margins") +
  labs(
    x = "Th2, eFluor 450-A (arcsinh)",
    y = "Density",
    title = "Th2 density per arm"
  )
save_panel(p_th2, th2, OUT_DIR, "fig_oc_density_th2", 183, 80)

# AF405 variant subsets, gap closed ------------------------------------------------------------------
subsets <- redraw_panel(
  file.path(
    RUNS,
    "2026-11-10_bmv-af405-mechanism/scores/fig_gap_closed_source.csv"
  ),
  "gap",
  "AF405 1.0.2 variant subsets",
  OUT_DIR,
  "proof_subsets_gap_closed"
)

# Zip 1.8.8 spillover spread matrix ------------------------------------------------------------------
# rows = spreading dye (the zip SSM reads row = spreading dye; script 82); diagonal empty
ssm <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-ssm-188/ssm_long.csv"),
  show_col_types = FALSE
)
dye_levels <- unique(ssm$row_dye)
ssm <- mutate(
  ssm,
  row_dye = factor(.data$row_dye, levels = rev(dye_levels)),
  col_dye = factor(.data$col_dye, levels = dye_levels)
)
p_ssm <- ggplot(ssm, aes(col_dye, row_dye, fill = ssm)) +
  geom_tile(colour = "white", linewidth = 0.1) +
  scale_fill_viridis_c(name = "SSM", na.value = "grey90") +
  coord_equal() +
  labs(
    x = "Receiving dye",
    y = "Spreading dye",
    title = "AutoSpectral 1.8.8 spillover spread (zip)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.text.y = element_text(size = 5)
  )
save_panel(p_ssm, ssm, OUT_DIR, "fig_ssm_188_zip", 183, 175)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
