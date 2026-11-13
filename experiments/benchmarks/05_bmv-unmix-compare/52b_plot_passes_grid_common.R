# Purpose: show the n.passes x n.af.passes grid (52_passes_grid.R) as one 3 x 3 tile panel per readout,
#          donor median with range, from the tables 51c_score_ablation_common.R wrote. No winner is picked: the
#          selection criterion is the user's (not yet set, 2026-11-10).
#          Common-gate copy (plans/common_gates_plan.md, batch 2): <score_dir> is the 51c run dir; cut 1.5 (G7a).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/52b_plot_passes_grid_common.R <score_dir>
#            from the repo root, after 51c_score_ablation_common.R <score_dir> p1_af1,...,p10_af10.
#            <score_dir> is run_in("common-passes-grid") (the 51c run folder, written to here).
#          Needs the TNF and IFNg readouts: stops when the dataset's targets are not TNF and IFNg (EXV).
# Inputs:  <score_dir>/scores.csv, tnf_among_ifng.csv, neg_spread.csv  [51c_score_ablation_common.R]
#          data_path("unmix_cells", "p*_af*", [<condition>,] "args.txt")  unmix minutes [52]; 24H has one per condition
#          datasets/<SW_DATASET>.json (targets, pos_rule, conditions), <ds>_samples.csv
# Outputs: <score_dir>/grid_summary.csv          median, min, max over donors per arm and readout
#          <score_dir>/fig_passes_grid.png       3 x 3 tiles per readout (+ fig_passes_grid_source.csv)
#          <score_dir>/runtime.csv               unmix minutes per arm (24H: per arm and condition)
#          <score_dir>/sessionInfo_52b.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))

# Parameters ------------------------------------------------------------------------------
TNF_LABEL <- cfg("targets")$tnf_label
IFNG_LABEL <- cfg("targets")$ifng_label
if (!identical(c(TNF_LABEL, IFNG_LABEL), c("TNF", "IFNg"))) {
  stop(
    "52b_plot_passes_grid_common.R: needs TNF and IFNg, missing in dataset ",
    DS$dataset,
    " (the two dye slots are ",
    TNF_LABEL,
    " and ",
    IFNG_LABEL,
    ")",
    call. = FALSE
  )
}
# the positive cut in the panel title: the fixed value, or the per-donor rule
POS_VALUE <- if (identical(cfg("pos_rule")$type, "fixed")) {
  cfg("pos_rule")$value
} else {
  "the positive cut"
}
SAMPLES <- samples_tbl()
COND_KEY <- if (is.null(DS$conditions)) character() else "condition"
PASS_GRID <- c(1, 5, 10) # as 52_passes_grid.R
READOUTS <- c(
  squeezed = sprintf("%s squeezed share", TNF_LABEL),
  pos = sprintf("%s share above %s", TNF_LABEL, POS_VALUE),
  tnf_ifng = sprintf("%% %s+ among %s+", TNF_LABEL, IFNG_LABEL),
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
    dplyr::filter(pop == POP, cytokine == TNF_LABEL) |>
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
if (!is.null(DS$conditions)) {
  # the sample id in `donor` carries the condition in 24H
  long <- long |>
    left_join(
      select(SAMPLES, donor = sample, condition),
      by = "donor",
      relationship = "many-to-one",
      unmatched = c("error", "drop")
    )
}

# Donor median and range per grid arm ---------------------------------------------------------
summary_tbl <- long |>
  inner_join(arms_grid, join_by(arm), relationship = "many-to-one") |>
  summarise(
    median = median(value),
    min = min(value),
    max = max(value),
    n_donors = n(),
    .by = c(arm, n_passes, n_af_passes, readout, any_of(COND_KEY))
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
# condition is NULL for 3D and EXV; 24H unmixes one task per condition (unmix_cells/<arm>/<condition>/)
read_minutes <- \(a, condition) {
  lines <- readLines(do.call(
    file.path,
    as.list(c(data_path("unmix_cells"), a, condition, "args.txt"))
  ))
  as.numeric(str_remove(
    str_subset(lines, "^unmix minutes:"),
    "^unmix minutes: "
  ))
}
runtime <- if (is.null(DS$conditions)) {
  arms_grid |>
    mutate(unmix_minutes = map_dbl(arm, \(a) read_minutes(a, NULL)))
} else {
  expand_grid(arms_grid, condition = DS$conditions) |>
    mutate(unmix_minutes = map2_dbl(arm, condition, read_minutes))
}
write_csv(runtime, file.path(score_dir, "runtime.csv"))

# Figure: one 3 x 3 tile panel per readout, colour scaled within each panel ------------------------
fig <- summary_tbl |>
  mutate(
    fill_scaled = (median - min(median)) / (max(median) - min(median) + 1e-12),
    .by = c(readout, any_of(COND_KEY))
  )
p <- ggplot(
  fig,
  aes(factor(n_af_passes), factor(n_passes), fill = fill_scaled)
) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label), size = 5, size.unit = "pt", lineheight = 0.9) +
  scale_fill_gradient(low = "grey95", high = "#4C9BD6", guide = "none") +
  (if (is.null(DS$conditions)) {
    facet_wrap(~readout, nrow = 2)
  } else {
    facet_grid(condition ~ readout)
  }) +
  labs(x = "n.af.passes", y = "n.passes") +
  theme(panel.grid = element_blank())
save_panel(
  p,
  select(fig, -label),
  score_dir,
  "fig_passes_grid",
  width_mm = 183,
  height_mm = 120 * max(1, length(DS$conditions))
)
writeLines(
  capture.output(sessionInfo()),
  file.path(score_dir, "sessionInfo_52b.txt")
)
