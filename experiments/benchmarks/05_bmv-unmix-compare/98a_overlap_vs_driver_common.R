# Purpose: does the spectral-overlap metric (80) predict which swapped dye drives each channel's squeeze?
#          Common-gate rerun of the old claim "TNF is the only channel where the metric predicts the observed
#          driver". The metric (80) reads spectra only (no gate). The observed drivers (76 common) come from
#          the common SpectroFlo gates at cut 1.5. This script only joins the two; it unmixes and gates nothing.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/98a_overlap_vs_driver_common.R \
#            experiments/runs/2026-11-11_bmv-common-overlap-vs-driver      (repo root; run dir must not exist)
# Inputs:  experiments/runs/2026-11-10_bmv-variant-overlap-matrix/overlap.csv   part "swapped": coefficient of
#            channel X's 1.8.8 spectrum in dye D's 1.0.2 variant minus 1.8.8 reference [80, gate-free]
#          experiments/runs/2026-11-11_bmv-per-dye-channel-map-common/map.csv   donor-median rise of the squeezed
#            share per swapped dye and channel, vs as188_joint [76 common]
# Outputs: <run_dir>/pairs.csv       per swapped dye D and channel X (D not X): coefficient and observed rise
#          <run_dir>/by_channel.csv  per channel: predicted top dye, observed top dye, match, rank of the observed
#                                    driver in the predicted order
#          <run_dir>/fig_overlap_vs_driver.{pdf,png,csv}  coefficient against observed rise, one point per pair
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
OVERLAP <- "experiments/runs/2026-11-10_bmv-variant-overlap-matrix/overlap.csv"
OBSERVED <- "experiments/runs/2026-11-11_bmv-per-dye-channel-map-common/map.csv"
LABELLED_CHANNELS <- c("BV510", "eFluor 450") # TNF and Th2: the two channels the INDEX names as hit hardest

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)

# Join the metric to the observed rise -------------------------------------------------------------
metric <- read_csv(OVERLAP, show_col_types = FALSE) |>
  filter(part == "swapped") |>
  transmute(swapped_dye = variant_dye, channel_dye, coef)
observed <- read_csv(OBSERVED, show_col_types = FALSE) |>
  mutate(channel_dye = sub("-A$", "", channel)) |>
  filter(swapped_dye != channel_dye) |> # the dye's own channel is not a spill target
  select(swapped_dye, channel_dye, marker, pop, rise = median_d_squeezed)
pairs <- inner_join(
  observed,
  metric,
  by = c("swapped_dye", "channel_dye"),
  relationship = "one-to-one"
)
stopifnot("every observed pair has a metric" = nrow(pairs) == nrow(observed))
write_csv(pairs, file.path(RUN_DIR, "pairs.csv"))

# Per channel: does the top predicted dye equal the top observed dye? -------------------------------
by_channel <- pairs |>
  group_by(channel_dye, marker, pop) |>
  mutate(
    pred_rank = rank(-coef, ties.method = "first"),
    obs_rank = rank(-rise, ties.method = "first")
  ) |>
  summarise(
    predicted_dye = swapped_dye[pred_rank == 1],
    observed_dye = swapped_dye[obs_rank == 1],
    observed_rise = max(rise),
    rank_of_observed_in_prediction = pred_rank[obs_rank == 1],
    .groups = "drop"
  ) |>
  mutate(match = predicted_dye == observed_dye) |>
  arrange(desc(observed_rise))
write_csv(by_channel, file.path(RUN_DIR, "by_channel.csv"))
print(by_channel, n = Inf)

# Figure -------------------------------------------------------------------------------------------
plot_df <- pairs |>
  mutate(
    label = if_else(
      channel_dye %in%
        LABELLED_CHANNELS &
        rise == ave(rise, channel_dye, FUN = max),
      paste0(marker, ": ", swapped_dye),
      NA_character_
    )
  )
p <- ggplot(plot_df, aes(coef, rise)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = FIG_LINE) +
  geom_point(colour = "grey40", size = 0.8, alpha = 0.6) +
  geom_point(
    data = filter(plot_df, !is.na(label)),
    colour = "#B2182B",
    size = 1.4
  ) +
  geom_text_repel(
    aes(label = label),
    size = 5 / .pt,
    family = FIG_FONT,
    min.segment.length = 0,
    na.rm = TRUE
  ) +
  labs(
    x = "Overlap coefficient",
    y = "Squeezed-share rise (donor median)"
  )
save_panel(p, plot_df, RUN_DIR, "fig_overlap_vs_driver", 89, 70)

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
