# Purpose: test 2 (README decision table, user 2026-11-10): for each of the 32 per-dye arms (one dye's 1.8.8
#          variants replaced by its 1.0.2 set, 45_variant_per_dye.R), score every channel with the test 1
#          readouts on as188_joint's gates, to find which swapped dye squeezes which channel.
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/76_per_dye_channel_map.R \
#            experiments/runs/2026-11-10_bmv-per-dye-channel-map      (repo root; run dir must not exist)
# Inputs:  experiments/bmv_test/compare_3D/data/unmix_cells/var102_only_<dye>/   raw unmixed arms [45]
#          experiments/runs/2026-11-10_bmv-gated-pooled/manifest.csv, ablation-gated/gs/as188_joint  [12, 19a]
#          experiments/runs/2026-11-10_bmv-all-marker-screen/marker_parent.csv   channel parents [75]
# Outputs: <run_dir>/scores.csv        per swapped dye, donor, population, channel: readouts
#          <run_dir>/map.csv           per swapped dye and channel (designated parent): donor medians of the
#                                      difference vs as188_joint (squeezed, above 1) and of the rSD ratio
#          <run_dir>/drivers.csv       per channel: the swapped dye with the largest median squeezed increase
#          <run_dir>/fig_per_dye_channel_map.png (+ _source.csv)   swapped dye x channel heatmap
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
CUT <- 1 # negatives below / positives above arcsinh 1 (user, 2026-11-10)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
POPS <- c("CD4T", "CD8T", "live")
# one parent per channel for the map: CD8T for T-cell markers (where the TNF squeeze is largest), live otherwise
MAP_POP_T <- "CD8T"
MARKER_PARENT <- "experiments/runs/2026-11-10_bmv-all-marker-screen/marker_parent.csv"

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/helper75_marker-scores.R")

# Arms: the 32 per-dye folders, named back to their dye ------------------------------------------------
marker_parent <- read_csv(MARKER_PARENT, show_col_types = FALSE)
dye_arms <- marker_parent |>
  transmute(dye, arm = paste0("var102_only_", gsub("[ -]", "_", dye))) |>
  filter(dir.exists(file.path(UNMIX_DIR, arm)))
arms <- c(as188_joint = "as188_joint", set_names(dye_arms$arm, dye_arms$dye))

# Scores per donor: every arm on as188_joint's population rows ------------------------------------------
pooled <- read_csv(POOLED_MANIFEST, show_col_types = FALSE)
gs <- load_gs("experiments/runs/2026-11-10_bmv-ablation-gated/gs/as188_joint")
donors <- sort(pooled$sample[pooled$arm == "as188_joint"])
scores <- donors |>
  set_names() |>
  map(\(s) {
    mats <- map(arms, \(a) arm_matrix(a, s, marker_parent$channel, pooled))
    pop_rows <- map(set_names(POPS), \(p) which(gh_pop_get_indices(gs[[s]], p)))
    map(mats, \(x) {
      map(pop_rows, \(r) score_channels(x, mats[["as188_joint"]], r)) |>
        list_rbind(names_to = "pop")
    }) |>
      list_rbind(names_to = "swapped_dye")
  }) |>
  list_rbind(names_to = "donor")
write_csv(scores, file.path(RUN_DIR, "scores.csv"))

# Map: donor medians of the change vs as188_joint on each channel's map parent ---------------------------
map_pop <- marker_parent |>
  transmute(
    channel,
    marker,
    pop = if_else(parents == "live", "live", MAP_POP_T)
  )
map_tbl <- scores |>
  inner_join(map_pop, by = c("channel", "pop"), relationship = "many-to-one") |>
  mutate(
    d_squeezed = squeezed_share - squeezed_share[swapped_dye == "as188_joint"],
    d_pos = pos_share - pos_share[swapped_dye == "as188_joint"],
    .by = c(donor, channel)
  ) |>
  filter(swapped_dye != "as188_joint") |>
  summarise(
    median_d_squeezed = median(d_squeezed),
    min_d_squeezed = min(d_squeezed),
    median_rsd_ratio = median(rsd_ratio),
    median_d_pos = median(d_pos),
    .by = c(swapped_dye, channel, marker, pop)
  )
write_csv(map_tbl, file.path(RUN_DIR, "map.csv"))
drivers <- map_tbl |>
  slice_max(median_d_squeezed, n = 1, by = channel, with_ties = FALSE) |>
  arrange(desc(median_d_squeezed))
write_csv(drivers, file.path(RUN_DIR, "drivers.csv"))

# Figure: swapped dye x channel, fill = median squeezed-share increase ------------------------------------
plot_df <- map_tbl |>
  mutate(channel_label = paste0(marker, " (", sub("-A$", "", channel), ")"))
lim <- max(abs(plot_df$median_d_squeezed))
p <- ggplot(
  plot_df,
  aes(channel_label, swapped_dye, fill = median_d_squeezed)
) +
  geom_tile() +
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-lim, lim),
    name = "Δ squeezed\nshare"
  ) +
  labs(
    x = "Scored channel (CD8T for T-cell markers, live otherwise)",
    y = "Dye given its 1.0.2 variants",
    title = "Squeezed share near 0, per-dye arm minus as188_joint (donor median)"
  ) +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(
  file.path(RUN_DIR, "fig_per_dye_channel_map.png"),
  p,
  width = 10,
  height = 7.5,
  dpi = 150,
  bg = "white"
)
write_csv(plot_df, file.path(RUN_DIR, "fig_per_dye_channel_map_source.csv"))

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
