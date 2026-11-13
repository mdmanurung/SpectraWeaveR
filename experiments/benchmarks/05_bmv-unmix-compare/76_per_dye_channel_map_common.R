# Purpose: Common-gates copy (plans/common_gates_plan.md, batch 3): gates and rows are the common SpectroFlo labels
#          (clean3), the squeeze cut is refitted as 0.5 x the negative rSD of as188_joint on those cells, and the
#          above-cut readout uses arcsinh 1.5 (G7a). Negatives for the rSD are below arcsinh 1.5 (G13). Run after 75 (common).
#          Test 2 (README decision table, user 2026-11-10): for each of the 32 per-dye arms (one dye's 1.8.8
#          variants replaced by its 1.0.2 set, 45_variant_per_dye.R), score every channel with the test 1
#          readouts on as188_joint's gates, to find which swapped dye squeezes which channel. Cuts come from
#          cut_of() per arm, sample and channel (3D: arcsinh 1.5). 24H: scores carry a condition column and the map
#          is per condition.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/76_per_dye_channel_map_common.R [<run_dir>]
#            (repo root; run dir must not exist; without <run_dir> it is run_out("per-dye-channel-map-common"))
# Inputs:  <data_dir>/unmix_cells/var102_only_<dye>/   raw unmixed arms (24H: one subfolder per condition) [45]
#          registry key common-gated: manifest.csv, gs/as188_joint  [batch 0, 97c, 97d]
#          registry key all-marker-screen-common: marker_parent.csv   channel parents [75 common]
#          registry key common-cuts, only when the dataset's pos_rule or neg_rule is not fixed: cuts.csv [97e]
#          datasets/<ds>.json (marker_parent) and datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <run_dir>/scores.csv        per swapped dye, donor, population, channel: readouts
#          <run_dir>/map.csv           per swapped dye and channel (designated parent; 24H: and condition): donor
#                                      medians of the difference vs as188_joint (squeezed, above the cut) and of
#                                      the rSD ratio
#          <run_dir>/drivers.csv       per channel (24H: and condition): the swapped dye with the largest median
#                                      squeezed increase
#          <run_dir>/fig_per_dye_channel_map.{pdf,png,csv}   swapped dye x channel heatmap
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
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), data_path(), samples_tbl()
# shared scoring functions; sourcing it also sources dataset.R again, which is harmless
source("experiments/benchmarks/05_bmv-unmix-compare/helper98_marker_scores.R")

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(RUN_DIR)) {
  RUN_DIR <- run_out("per-dye-channel-map-common")
}
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
# negative cut (rSD, squeeze cut; G13) and positive cut (share above; G7a): cut_of(arm, sample, channel, kind)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
# config marker_parent: default = populations scored for every marker, by_marker, no_marker as in 75; map_pop =
# the one parent of a multi-parent channel in the map (3D: CD8T for T-cell markers, where the TNF squeeze is largest)
PARENTS <- marker_parent_cfg()
POPS <- unique(c(PARENTS$default, unlist(PARENTS$by_marker), PARENTS$no_marker))
MAP_POP_T <- PARENTS$map_pop
COND_BY <- if (is.null(DS$conditions)) character() else "condition" # 24H only: map per condition
MARKER_PARENT <- file.path(
  run_in("all-marker-screen-common"),
  "marker_parent.csv"
)

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()

# Arms: the 32 per-dye folders, named back to their dye ------------------------------------------------
marker_parent <- read_csv(MARKER_PARENT, show_col_types = FALSE)
dye_arms <- marker_parent |>
  transmute(dye, arm = paste0("var102_only_", gsub("[ -]", "_", dye))) |>
  filter(dir.exists(file.path(UNMIX_DIR, arm)))
arms <- c(as188_joint = "as188_joint", set_names(dye_arms$arm, dye_arms$dye))

# Scores per donor: every arm on as188_joint's population rows ------------------------------------------
pooled <- read_csv(POOLED_MANIFEST, show_col_types = FALSE)
gs <- load_gs(file.path(run_in("common-gated"), "gs/as188_joint"))
donors <- sort(pooled$sample[pooled$arm == "as188_joint"])
stopifnot(
  "scored samples are the sample table" = setequal(donors, samples_tbl()$sample)
)
scores <- donors |>
  set_names() |>
  map(\(s) {
    mats <- map(arms, \(a) arm_matrix(a, s, marker_parent$channel, pooled))
    neg_cuts <- map(arms, \(a) {
      cuts_by_channel(a, s, marker_parent$channel, "neg")
    })
    pos_cuts <- map(arms, \(a) {
      cuts_by_channel(a, s, marker_parent$channel, "pos")
    })
    pop_rows <- map(set_names(POPS), \(p) which(gh_pop_get_indices(gs[[s]], p)))
    imap(mats, \(x, dye) {
      map(pop_rows, \(r) {
        score_channels(
          x,
          mats[["as188_joint"]],
          r,
          neg_cut = neg_cuts[[dye]],
          pos_cut = pos_cuts[[dye]],
          ref_neg_cut = neg_cuts[["as188_joint"]]
        )
      }) |>
        list_rbind(names_to = "pop")
    }) |>
      list_rbind(names_to = "swapped_dye")
  }) |>
  list_rbind(names_to = "donor") # the sample id (3D: the well)
if (!is.null(DS$conditions)) {
  scores <- scores |>
    left_join(
      select(samples_tbl(), donor = sample, condition),
      by = "donor",
      relationship = "many-to-one",
      unmatched = "error"
    )
}
write_csv(scores, file.path(RUN_DIR, "scores.csv"))

# Map: donor medians of the change vs as188_joint on each channel's map parent ---------------------------
map_pop <- marker_parent |>
  transmute(
    channel,
    marker,
    pop = if_else(grepl(";", parents, fixed = TRUE), MAP_POP_T, parents)
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
    .by = c(swapped_dye, channel, marker, pop, all_of(COND_BY))
  )
write_csv(map_tbl, file.path(RUN_DIR, "map.csv"))
drivers <- map_tbl |>
  slice_max(
    median_d_squeezed,
    n = 1,
    by = c(channel, all_of(COND_BY)),
    with_ties = FALSE
  ) |>
  arrange(desc(median_d_squeezed))
write_csv(drivers, file.path(RUN_DIR, "drivers.csv"))

# Figure: swapped dye x channel, fill = median squeezed-share increase ------------------------------------
plot_df <- map_tbl |>
  mutate(channel_label = paste0(marker, " (", sub("-A$", "", channel), ")"))
lim <- max(abs(plot_df$median_d_squeezed))
# 24H: one panel per condition
facet_cond <- if (is.null(DS$conditions)) {
  NULL
} else {
  facet_wrap(vars(condition), nrow = 1)
}
p <- ggplot(
  plot_df,
  aes(channel_label, swapped_dye, fill = median_d_squeezed)
) +
  geom_tile() +
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-lim, lim),
    name = "\u0394 squeezed\nshare"
  ) +
  labs(
    x = "Scored channel (CD8T for T-cell markers, live otherwise)",
    y = "Dye given its 1.0.2 variants"
  ) +
  facet_cond +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    plot.margin = margin(5, 5, 16, 5) # room for the rotated labels
  )
save_panel(p, plot_df, RUN_DIR, "fig_per_dye_channel_map", 183, 160)

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
