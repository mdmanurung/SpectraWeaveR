# Purpose: apply the REVISED affected-marker flag rule (README "Affected-marker flag, revised", user
#          2026-11-10; supersedes 86a) to the all-marker screen (75), the per-dye channel map (76) and the Th2
#          proof arms (86). A marker is flagged in a population if the donor-median squeezed share rises >= 0.10
#          vs as188_joint and rises in all 5 donors. The rSD criterion of 86a is dropped. The rSD ratio is kept
#          in flags.csv for reference only.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/86a_apply_flag_rule.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-all-marker-screen/scores.csv        per donor, arm, pop, channel [75]
#          experiments/runs/2026-11-10_bmv-all-marker-screen/marker_parent.csv channel -> marker, dye [75]
#          experiments/runs/2026-11-10_bmv-per-dye-channel-map/scores.csv      per donor, swapped dye [76]
#          experiments/runs/2026-11-10_bmv-th2-proof/scores.csv                per donor, proof arm [86]
# Outputs: experiments/runs/2026-11-10_bmv-marker-flags-v2/flags.csv   source, arm or swapped dye, pop, marker:
#                                                                  the three rule components and the flag
#          experiments/runs/2026-11-10_bmv-marker-flags-v2/fig_flags_screen_<pop>.{pdf,png,csv}   screen heatmap per
#            population, flagged cells outlined (save_panel: PDF, PNG, plotted frame)
#          experiments/runs/2026-11-10_bmv-marker-flags-v2/fig_flags_per_dye_<pop>.{pdf,png,csv}  per-dye heatmap
#          experiments/runs/2026-11-10_bmv-marker-flags-v2/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
})
TEMPLATES <- "experiments/benchmarks/05_bmv-unmix-compare/lucid_templates" # lucid-coder 0.5.0 figure style
source(file.path(TEMPLATES, "theme.R"))
source(file.path(TEMPLATES, "palettes.R"))
source(file.path(TEMPLATES, "plots.R"))

# Parameters ---------------------------------------------------------------------------------------
D_SQUEEZE_MIN <- 0.10 # revised rule (user, 2026-11-10)
N_DONORS <- 5 # every donor must rise (revised rule)
REF_ARM <- "as188_joint"
HEAT_LIMIT <- 0.5 # colour scale saturates at +/-0.5 squeezed share (largest change 0.53); beyond is squished
POPS <- c("CD4T", "CD8T", "live")
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
OUT_DIR <- "experiments/runs/2026-11-10_bmv-marker-flags-v2"

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Per-donor change vs as188_joint, then the rule per arm, pop and channel ----------------------------
marker_of <- read_csv(
  "experiments/runs/2026-11-10_bmv-all-marker-screen/marker_parent.csv",
  show_col_types = FALSE
) |>
  select(channel, dye, marker)
inputs <- list(
  screen = read_csv(
    "experiments/runs/2026-11-10_bmv-all-marker-screen/scores.csv",
    show_col_types = FALSE
  ),
  per_dye = read_csv(
    "experiments/runs/2026-11-10_bmv-per-dye-channel-map/scores.csv",
    show_col_types = FALSE
  ) |>
    rename(arm = swapped_dye),
  th2_proof = read_csv(
    "experiments/runs/2026-11-10_bmv-th2-proof/scores.csv",
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
      .by = c(arm, pop, channel)
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

for (which_pop in c("CD4T", "CD8T")) {
  heat_panel("th2_proof", which_pop, 89, 89)
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
