# Purpose: redraw the report's gap-closed and marker-flag figures (sections 16 and 17) in the Nature
#          figure style (user, 2026-11-11; lucid-coder r-0008, r-0012, r-0014 to r-0017) from the tables
#          their runs saved. Nothing is recomputed; values are the originals.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_gap_flags.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-ablation-sbuv510-af405/fig_gap_closed_source.csv   [51_score_ablation.R]
#          experiments/runs/2026-11-10_bmv-af405-forced-af/gap_closed.csv     refined AF arm [87_af405_forced_af.R]
#          experiments/runs/2026-11-10_bmv-af405-afremove/gap_closed.csv      AF-removed arm [87b_af405_afremove.R]
#          experiments/runs/2026-11-10_bmv-marker-flags-v2/fig_flags_screen_CD8T.csv      [86b_apply_flag_rule_v2.R]
#          experiments/runs/2026-11-10_bmv-marker-flags-v2/fig_flags_th2_proof_CD8T.csv   [86b_apply_flag_rule_v2.R]
# Outputs: experiments/runs/2026-11-11_bmv-restyle/16_scatter_spectroflo/sbuv510_af405_gap_closed.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/17_flags_th2_af/afrefine_gap_closed.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/17_flags_th2_af/afremove_gap_closed.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/17_flags_th2_af/flags_screen_CD8T.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/17_flags_th2_af/flags_th2_proof_CD8T.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle/17_flags_th2_af/sessionInfo_95_gap_flags.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))
source(file.path(LIB, "plots.R"))

# Parameters ---------------------------------------------------------------------------------------
RUNS <- "experiments/runs"
OUT <- "experiments/runs/2026-11-11_bmv-restyle"
DONORS <- c("A3", "A6", "B11", "B2", "D1")
DONOR_SHAPES <- setNames(c(16, 17, 15, 18, 8), DONORS) # as the approved demo (2026-11-11)
GAP_SOLID <- c(0, 1) # 51: no gap closed, full gap closed
GAP_DASHED <- c(0.1, 0.5, 0.8) # 51: CONTRIB_FRAC, RESEMBLE_MIN, RESEMBLE_MEDIAN
HEAT_LIMIT <- 0.5 # 86b: colour scale saturates at +/-0.5 squeezed share
JITTER_WIDTH <- 0.15 # 87, 87b: display only
JITTER_SEED <- 1

out_16 <- file.path(OUT, "16_scatter_spectroflo")
out_17 <- file.path(OUT, "17_flags_th2_af")
dir.create(out_16, recursive = TRUE, showWarnings = FALSE)
dir.create(out_17, recursive = TRUE, showWarnings = FALSE)
stopifnot(
  "no restyled F4 file exists yet" = !any(file.exists(c(
    file.path(out_16, "sbuv510_af405_gap_closed.png"),
    file.path(out_17, c("afrefine_gap_closed.png", "afremove_gap_closed.png")),
    file.path(out_17, c("flags_screen_CD8T.png", "flags_th2_proof_CD8T.png"))
  )))
)
donor_pal <- categorical_palette(DONORS)

# Section 16: gap closed by the SBUV510 + AF405 1.0.2 variant arm -----------------------------------
# Donors as points, donor median as a black diamond; thresholds as in script 51.
gap_16 <- read_csv(
  file.path(
    RUNS,
    "2026-11-10_bmv-ablation-sbuv510-af405/fig_gap_closed_source.csv"
  ),
  show_col_types = FALSE
) |>
  mutate(
    arm = factor(arm, levels = rev(unique(arm))),
    readout = factor(
      readout,
      levels = c("CD8T TNF squeezed near 0", "CD8T TNF above arcsinh 1")
    )
  )
gap_16_median <- gap_16 |>
  summarise(median_frac = median(frac_closed), .by = c(arm, readout))
p_16 <- ggplot(gap_16, aes(frac_closed, arm)) +
  geom_vline(
    xintercept = GAP_SOLID,
    colour = "grey50",
    linewidth = FIG_LINE / 2
  ) +
  geom_vline(
    xintercept = GAP_DASHED,
    linetype = "dashed",
    colour = "grey70",
    linewidth = FIG_LINE / 2
  ) +
  geom_point(aes(colour = donor, shape = donor), size = 1.3) +
  geom_point(
    data = gap_16_median,
    aes(median_frac, arm),
    shape = 18,
    size = 2.5
  ) +
  scale_colour_manual(values = donor_pal, name = "Donor") +
  scale_shape_manual(values = DONOR_SHAPES, name = "Donor") +
  facet_wrap(~readout, axes = "all", axis.labels = "margins") +
  labs(
    x = "Fraction of the as188_joint to s102_c188 gap closed",
    y = NULL,
    title = "SBUV510 + AF405 1.0.2 variants: CD8T TNF gap closed"
  ) +
  value_grid("x")
save_panel(
  p_16,
  gap_16,
  out_16,
  "sbuv510_af405_gap_closed",
  width_mm = 183,
  height_mm = 45
)

# Section 17: refined AF and AF removed, CD8T TNF gap closed per donor -------------------------------
# One readout per group, donors as points over the median (plot_box), 0 dashed; same layout for both arms.
gap_17 <- bind_rows(
  afrefine = read_csv(
    file.path(RUNS, "2026-11-10_bmv-af405-forced-af/gap_closed.csv"),
    show_col_types = FALSE
  ),
  afremove = read_csv(
    file.path(RUNS, "2026-11-10_bmv-af405-afremove/gap_closed.csv"),
    show_col_types = FALSE
  ),
  .id = "which_arm"
)
gap_pal <- categorical_palette(c("TNF squeezed near 0", "TNF above arcsinh 1"))
gap_titles <- c(
  afrefine = "Refined AF: CD8T TNF gap closed",
  afremove = "AF removed: CD8T TNF gap closed"
)
for (which_gap in names(gap_titles)) {
  d <- filter(gap_17, .data$which_arm == which_gap)
  p <- plot_box(
    d,
    group = "readout",
    value = "frac_closed",
    pal = gap_pal,
    jitter_width = JITTER_WIDTH,
    seed = JITTER_SEED,
    value_label = "Fraction of 1.8.8 to 1.0.2 gap closed",
    title = gap_titles[[which_gap]]
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey60",
      linewidth = FIG_LINE / 2
    ) +
    value_grid("y")
  save_panel(
    p,
    d,
    out_17,
    paste0(which_gap, "_gap_closed"),
    width_mm = 89,
    height_mm = 60
  )
}

# Section 17: marker-flag heatmaps, CD8T ---------------------------------------------------------------
# Change in squeezed share vs as188_joint, diverging at 0; flagged cells outlined; NA cells grey.
heat_layers <- function(d) {
  list(
    geom_tile(
      data = filter(d, flag),
      colour = "black",
      fill = NA,
      linewidth = FIG_LINE * 1.2
    ),
    geom_text(aes(label = cell_text), size = 5.5, size.unit = "pt"),
    guides(
      fill = guide_colourbar(
        barheight = unit(25, "mm"),
        barwidth = unit(3, "mm")
      )
    )
  )
}
fill_label <- "Δ squeezed\nshare vs\nas188_joint"

# Screen: markers down y; arm IDs wrapped at "_" so they read horizontally (r-0012)
screen <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-marker-flags-v2/fig_flags_screen_CD8T.csv"),
  col_types = cols(cell_text = col_character()) # keep "0.20", not 0.2
) |>
  mutate(
    cell_text = if_else(
      is.na(median_d_squeezed),
      "NA",
      coalesce(as.character(cell_text), "")
    ),
    arm = factor(
      sub("_", "\n", arm, fixed = TRUE),
      levels = sort(unique(sub("_", "\n", arm, fixed = TRUE)))
    ),
    label = factor(label, levels = sort(unique(label)))
  )
p_screen <- plot_heatmap(
  screen,
  row = "label",
  col = "arm",
  value = "median_d_squeezed",
  limit = HEAT_LIMIT,
  title = "Marker squeeze, CD8 T cells",
  x_label = "Unmixing arm",
  y_label = NULL,
  fill_label = fill_label
) +
  heat_layers(screen)
save_panel(
  p_screen,
  screen,
  out_17,
  "flags_screen_CD8T",
  width_mm = 89,
  height_mm = 150
)

# Th2 proof: long arm IDs run down y, the three markers across x (flipped, r-0012)
proof <- read_csv(
  file.path(
    RUNS,
    "2026-11-10_bmv-marker-flags-v2/fig_flags_th2_proof_CD8T.csv"
  ),
  col_types = cols(cell_text = col_character())
) |>
  mutate(
    cell_text = if_else(
      is.na(median_d_squeezed),
      "NA",
      coalesce(as.character(cell_text), "")
    ),
    arm = factor(arm, levels = sort(unique(arm))),
    # one word per line so the three marker labels read side by side at 89 mm
    label = gsub(" ", "\n", label, fixed = TRUE),
    label = factor(label, levels = sort(unique(label)))
  )
p_proof <- plot_heatmap(
  proof,
  row = "arm",
  col = "label",
  value = "median_d_squeezed",
  limit = HEAT_LIMIT,
  title = "AF405 proof arms, CD8 T",
  x_label = NULL,
  y_label = NULL, # the title names the arms
  fill_label = fill_label
) +
  heat_layers(proof)
save_panel(
  p_proof,
  proof,
  out_17,
  "flags_th2_proof_CD8T",
  width_mm = 89,
  height_mm = 65
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_17, "sessionInfo_95_gap_flags.txt")
)
