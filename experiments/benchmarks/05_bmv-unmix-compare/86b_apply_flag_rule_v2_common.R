# Purpose: Common-gates copy (plans/common_gates_plan.md, batch 3): reads the common-gate runs of 75, 76 and 86; the
#          flag rule is unchanged. Original purpose: apply the REVISED affected-marker flag rule (README "Affected-marker flag, revised", user
#          2026-11-10; supersedes 86a) to the all-marker screen (75), the per-dye channel map (76) and the Th2
#          proof arms (86). A marker is flagged in a population if the donor-median squeezed share rises >= 0.10
#          vs as188_joint and rises in all 5 donors. The rSD criterion of 86a is dropped. The rSD ratio is kept
#          in flags.csv for reference only.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/86b_apply_flag_rule_v2_common.R   (repo root)
#          Stops on a dataset without the Th2 dye (config targets th2_dye): the Th2 proof is one of its sources.
# Inputs:  registry key all-marker-screen-common: scores.csv per donor, arm, pop, channel [75]; marker_parent.csv
#            channel -> marker, dye [75]
#          registry key per-dye-channel-map-common: scores.csv per donor, swapped dye [76]
#          registry key th2-proof-common: scores.csv per donor, proof arm [86]
#          datasets/<ds>.json (targets, marker_parent) and datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <run_out("marker-flags-v2-common")>/flags.csv   source, arm or swapped dye, pop, marker:
#                                                                  the three rule components and the flag
#          <run_out("marker-flags-v2-common")>/fig_flags_screen_<pop>.{pdf,png,csv}   screen heatmap per
#            population, flagged cells outlined (save_panel: PDF, PNG, plotted frame)
#          <run_out("marker-flags-v2-common")>/fig_flags_per_dye_<pop>.{pdf,png,csv}  per-dye heatmap
#          <run_out("marker-flags-v2-common")>/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), run_out(), samples_tbl()
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()
source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lib/plots.R")

# Parameters ---------------------------------------------------------------------------------------
if (is.null(cfg("targets")$th2_dye)) {
  stop(
    "86b_apply_flag_rule_v2_common.R: needs target th2_dye, missing in dataset ",
    DS$dataset,
    call. = FALSE
  )
}
D_SQUEEZE_MIN <- 0.10 # revised rule (user, 2026-11-10)
# every donor must rise (revised rule): samples per condition, one per donor (3D, EXV: no condition, all samples)
N_DONORS <- samples_tbl() |>
  count(condition) |>
  pull(n) |>
  unique()
stopifnot(
  "86b_apply_flag_rule_v2_common.R: unequal samples per condition" = length(
    N_DONORS
  ) ==
    1
)
REF_ARM <- "as188_joint"
HEAT_LIMIT <- 0.5 # colour scale saturates at +/-0.5 squeezed share (largest change 0.53); beyond is squished
# populations scored in 75 and 76, from config marker_parent (3D: CD4T, CD8T, live)
MARKER_PARENT <- cfg("marker_parent")
stopifnot(
  "marker_parent needs default and no_marker" = !is.null(
    MARKER_PARENT$default
  ) &&
    !is.null(MARKER_PARENT$no_marker)
)
POPS <- unique(c(
  MARKER_PARENT$default,
  unlist(MARKER_PARENT$by_marker),
  MARKER_PARENT$no_marker
))
COND_BY <- if (is.null(DS$conditions)) character() else "condition" # 24H only: rule per condition
TITLES <- list(
  screen = list(
    CD4T = "Marker squeeze, CD4 T cells",
    CD8T = "Marker squeeze, CD8 T cells",
    live = "Marker squeeze, live cells"
  ),
  per_dye = list(
    CD4T = "Per-dye 1.0.2 variants, CD4 T",
    CD8T = "Per-dye 1.0.2 variants, CD8 T",
    live = "Per-dye 1.0.2 variants, live"
  ),
  th2_proof = list(
    CD4T = "AF405 proof arms, CD4 T",
    CD8T = "AF405 proof arms, CD8 T"
  )
)
X_LABELS <- c(
  screen = "Unmixing arm",
  per_dye = "Dye given its 1.0.2 variants",
  th2_proof = "AF405 proof arm"
)
OUT_DIR <- run_out("marker-flags-v2-common")
stopifnot(
  "every scored population has a title in TITLES" = all(
    POPS %in% names(TITLES$screen)
  )
)
dir.create(OUT_DIR, recursive = TRUE)

# Per-donor change vs as188_joint, then the rule per arm, pop and channel ----------------------------
marker_of <- read_csv(
  file.path(run_in("all-marker-screen-common"), "marker_parent.csv"),
  show_col_types = FALSE
) |>
  select(channel, dye, marker)
inputs <- list(
  screen = read_csv(
    file.path(run_in("all-marker-screen-common"), "scores.csv"),
    show_col_types = FALSE
  ),
  per_dye = read_csv(
    file.path(run_in("per-dye-channel-map-common"), "scores.csv"),
    show_col_types = FALSE
  ) |>
    rename(arm = swapped_dye),
  th2_proof = read_csv(
    file.path(run_in("th2-proof-common"), "scores.csv"),
    show_col_types = FALSE
  ) |>
    select(-marker)
)
flags <- imap(inputs, \(scores, which_source) {
  ref <- scores |>
    filter(arm == REF_ARM) |>
    select(donor, pop, channel, ref_squeezed = squeezed_share)
  # each arm's donor rows joined to the same donor's as188_joint row
  scores |>
    filter(arm != REF_ARM) |>
    left_join(
      ref,
      by = join_by(donor, pop, channel),
      relationship = "many-to-one",
      unmatched = "error"
    ) |>
    mutate(d_squeezed = squeezed_share - ref_squeezed) |>
    # rule components across donors
    summarise(
      n_donors = n(),
      median_d_squeezed = median(d_squeezed),
      median_rsd_ratio = median(rsd_ratio),
      all_donors_up = all(d_squeezed > 0),
      .by = c(arm, pop, channel, all_of(COND_BY))
    ) |>
    mutate(
      source = which_source,
      flag = n_donors == N_DONORS &
        !is.na(median_d_squeezed) &
        median_d_squeezed >= D_SQUEEZE_MIN &
        all_donors_up
    )
}) |>
  list_rbind() |>
  left_join(
    marker_of,
    by = "channel",
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(
    flag = coalesce(flag, FALSE),
    # PE has no marker in the control file; say so instead of printing "NA"
    label = sprintf("%s (%s)", coalesce(marker, "no marker"), dye)
  )
write_csv(flags, file.path(OUT_DIR, "flags.csv"))

# Figures: one heatmap per source and population, house style ---------------------------------------
# Change in squeezed share on a diverging scale centred at 0; flagged cells outlined in black; cells with
# no score (CD3 has no negatives inside its own T-cell parents) drawn grey with "NA", never 0.
heat_panel <- function(which_source, which_pop, width_mm, height_mm) {
  d <- flags |>
    filter(source == which_source, pop == which_pop) |>
    mutate(
      arm = factor(arm, levels = sort(unique(arm))),
      label = factor(label, levels = sort(unique(label))),
      cell_text = case_when(
        is.na(median_d_squeezed) ~ "NA",
        abs(median_d_squeezed) >= D_SQUEEZE_MIN ~ sprintf(
          "%.2f",
          median_d_squeezed
        ),
        .default = ""
      )
    )
  p <- plot_heatmap(
    d,
    row = "label",
    col = "arm",
    value = "median_d_squeezed",
    limit = HEAT_LIMIT,
    title = TITLES[[which_source]][[which_pop]],
    x_label = X_LABELS[[which_source]],
    y_label = NULL,
    fill_label = "\u0394 squeezed share\nvs as188_joint"
  ) +
    geom_tile(
      data = filter(d, flag),
      colour = "black",
      fill = NA,
      linewidth = 0.6
    ) +
    geom_text(aes(label = cell_text), size = 6, size.unit = "pt") +
    # legend under the panel, title above the bar: keeps the columns wide at 89 mm
    theme(legend.position = "bottom") +
    guides(
      fill = guide_colourbar(title.position = "top", barwidth = unit(40, "mm"))
    ) +
    rotate_x()
  # 24H: one panel per condition
  if (!is.null(DS$conditions)) {
    p <- p + facet_wrap(vars(condition))
  }
  save_panel(
    p,
    d,
    OUT_DIR,
    sprintf("fig_flags_%s_%s", which_source, which_pop),
    width_mm,
    height_mm
  )
}
for (which_pop in POPS) {
  heat_panel("screen", which_pop, 89, 160)
}
for (which_pop in intersect(
  POPS,
  unique(flags$pop[flags$source == "per_dye"])
)) {
  heat_panel("per_dye", which_pop, 183, 183)
}

for (which_pop in MARKER_PARENT$default) {
  heat_panel("th2_proof", which_pop, 89, 89)
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
