# Purpose: show the n.passes x n.af.passes grid (52_passes_grid.R) as one 3 x 3 tile panel per readout,
#          donor median with range, from the tables 51_score_ablation.R wrote. No winner is picked: the
#          selection criterion is the user's (not yet set, 2026-11-10).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/52b_plot_passes_grid.R <score_dir>
#            from the repo root, after 51_score_ablation.R <score_dir> p1_af1,...,p10_af10.
# Inputs:  <score_dir>/scores.csv, tnf_among_ifng.csv, neg_spread.csv  [51_score_ablation.R]
#          experiments/bmv_test/compare_3D/data/unmix_cells/p*_af*/args.txt  unmix minutes [52]
# Outputs: <score_dir>/grid_summary.csv          median, min, max over donors per arm and readout
#          <score_dir>/fig_passes_grid.png       3 x 3 tiles per readout (+ fig_passes_grid_source.csv)
#          <score_dir>/runtime.csv               unmix minutes per arm
#          <score_dir>/sessionInfo_52b.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Parameters ------------------------------------------------------------------------------
PASS_GRID <- c(1, 5, 10) # as 52_passes_grid.R
READOUTS <- c(
  squeezed = "TNF squeezed share",
  pos = "TNF share above 1",
  tnf_ifng = "% TNF+ among IFNg+",
  th2 = "Th2 neg. rSD / as188",
  il10 = "IL-10 neg. rSD / as188"
) # all CD8T; panel titles in this order
POP <- "CD8T"

score_dir <- commandArgs(trailingOnly = TRUE)[[1]]
arms_grid <- tidyr::expand_grid(
  n_passes = PASS_GRID,
  n_af_passes = PASS_GRID
) |>
  mutate(arm = paste0("p", n_passes, "_af", n_af_passes))

# Readouts per donor and arm, long --------------------------------------------------------------
scores <- read_csv(file.path(score_dir, "scores.csv"), show_col_types = FALSE)
tnf_ifng <- read_csv(
  file.path(score_dir, "tnf_among_ifng.csv"),
  show_col_types = FALSE
)
neg_spread <- read_csv(
  file.path(score_dir, "neg_spread.csv"),
  show_col_types = FALSE
)
long <- bind_rows(
  scores |>
    dplyr::filter(pop == POP, cytokine == "TNF") |>
    transmute(donor, arm, squeezed = squeezed_share, pos = pos_share) |>
    pivot_longer(c(squeezed, pos), names_to = "readout"),
  tnf_ifng |>
    dplyr::filter(pop == POP) |>
    transmute(donor, arm, readout = "tnf_ifng", value = pct_tnf_among_ifng),
  neg_spread |>
    dplyr::filter(pop == POP, cytokine %in% c("Th2", "IL-10")) |>
    transmute(
      donor,
      arm,
      readout = if_else(cytokine == "Th2", "th2", "il10"),
      value = rsd_ratio
    )
)

# Donor median and range per grid arm ---------------------------------------------------------
summary_tbl <- long |>
  inner_join(arms_grid, join_by(arm), relationship = "many-to-one") |>
  summarise(
    median = median(value),
    min = min(value),
    max = max(value),
    n_donors = n(),
    .by = c(arm, n_passes, n_af_passes, readout)
  ) |>
  mutate(
    readout = factor(READOUTS[readout], levels = READOUTS),
    label = sprintf("%.2f\n(%.2f-%.2f)", median, min, max)
  )
stopifnot(
  "every grid arm scored" = n_distinct(summary_tbl$arm) == nrow(arms_grid)
)
write_csv(select(summary_tbl, -label), file.path(score_dir, "grid_summary.csv"))

# Runtime from args.txt -----------------------------------------------------------------------
runtime <- arms_grid |>
  mutate(
    unmix_minutes = map_dbl(arm, \(a) {
      lines <- readLines(file.path(
        "experiments/bmv_test/compare_3D/data/unmix_cells",
        a,
        "args.txt"
      ))
      as.numeric(str_remove(
        str_subset(lines, "^unmix minutes:"),
        "^unmix minutes: "
      ))
    })
  )
write_csv(runtime, file.path(score_dir, "runtime.csv"))

# Figure: one 3 x 3 tile panel per readout, colour scaled within each panel ------------------------
fig <- summary_tbl |>
  mutate(
    fill_scaled = (median - min(median)) / (max(median) - min(median) + 1e-12),
    .by = readout
  )
write_csv(
  select(fig, -label),
  file.path(score_dir, "fig_passes_grid_source.csv")
)
p <- ggplot(
  fig,
  aes(factor(n_af_passes), factor(n_passes), fill = fill_scaled)
) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label), size = 3.2) +
  scale_fill_gradient(low = "grey95", high = "#4C9BD6", guide = "none") +
  facet_wrap(~readout, nrow = 1) +
  labs(x = "n.af.passes", y = "n.passes") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10))
ggsave(
  file.path(score_dir, "fig_passes_grid.png"),
  p,
  width = 15,
  height = 4,
  dpi = 130,
  bg = "white"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(score_dir, "sessionInfo_52b.txt")
)
