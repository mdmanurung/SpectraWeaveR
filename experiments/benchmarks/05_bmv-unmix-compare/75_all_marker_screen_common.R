# Purpose: Common-gates copy (plans/common_gates_plan.md, batch 3): gates and rows are the common SpectroFlo labels
#          (clean3), the squeeze cut is refitted as 0.5 x the negative rSD of as188_joint on those cells, and the
#          above-cut readout uses the positive cut (G7a). Negatives for the rSD are below the negative cut (G13). Both
#          cuts come from cut_of() per arm, sample and channel (3D: arcsinh 1.5 everywhere).
#          Test 1 (README decision table, user 2026-11-10): check every fluorescence channel for the TNF-type
#          squeeze. For each arm, donor and channel: squeezed share near 0, negative rSD ratio vs as188_joint and
#          share above the positive cut (39b / 51 metrics, cofactor 6000 in memory), on as188_joint's gates (rows paired
#          across arms). Parents come from config marker_parent: CD4T and CD8T for T-cell markers, live cells
#          for the markers it lists (3D: Viability, CD56) and for channels without a marker (marker_parent.csv).
#          24H: scores carry a condition column and the medians are per condition.
#          51 has no channel-level "affected" rule (its labels score arms by the TNF gap), so no flag is set.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/75_all_marker_screen_common.R [<run_dir>]
#            (repo root; run dir must not exist; without <run_dir> it is run_out("all-marker-screen-common"))
# Inputs:  registry key common-gated (datasets/<ds>_runs.csv): manifest.csv, the clean3 FCS of every arm [97c], and
#            gs/as188_joint, the common-gate GatingSet [batch 0, 97d]
#          registry key common-cuts, only when the dataset's pos_rule or neg_rule is not fixed: cuts.csv [97e]
#          <data_dir>/unmix_cells/<arm>/   raw unmixed arms, e.g. swap_var102, swap_af102 [37_swap_af_variants.R]
#          <data_dir>/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv   dye to marker map [user zip]
#          datasets/<ds>.json (marker_parent, channels) and datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <run_dir>/marker_parent.csv   channel, dye, marker, parent populations, ambiguity note
#          <run_dir>/scores.csv          per arm, donor, population, channel: readouts (all three populations)
#          <run_dir>/summary.csv         per arm, population, channel (designated parents; 24H: and condition): donor
#                                        medians of the difference vs as188_joint (squeezed, above the cut) and of the
#                                        rSD ratio
#          <run_dir>/fig_screen_squeezed.{pdf,png,csv}   squeezed-share difference, channel x arm
#          <run_dir>/fig_screen_rsd.{pdf,png,csv}        log2 rSD ratio, channel x arm
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
  RUN_DIR <- run_out("all-marker-screen-common")
}
ARMS <- c(
  "as188_joint",
  "as102_var",
  "s102_c188",
  "s188_c102",
  "swap_var102",
  "swap_af102"
)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
# negative cut (rSD, squeeze cut; G13) and positive cut (share above; G7a): cut_of(arm, sample, channel, kind)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
# config marker_parent: default = populations scored for every marker, by_marker = markers with their own parent
# (3D: Viability, CD56 in live), no_marker = parent of a channel without a marker (3D: live), also_non_t = markers
# assigned to T cells but also expressed outside them (listed in the report as ambiguous)
MARKER_PARENT <- marker_parent_cfg()
POPS <- unique(c(
  MARKER_PARENT$default,
  unlist(MARKER_PARENT$by_marker),
  MARKER_PARENT$no_marker
))
COND_BY <- if (is.null(DS$conditions)) character() else "condition" # 24H only: medians per condition
CONTROL_FILE <- data_path(
  "autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"
)

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()

# Channels and parents from the control file ----------------------------------------------------------
marker_parent <- read_csv(CONTROL_FILE, show_col_types = FALSE) |>
  filter(!fluorophore %in% c("AF", "Negative")) |>
  transmute(
    channel = paste0(fluorophore, "-A"),
    dye = fluorophore,
    marker,
    parents = map_chr(marker, \(m) {
      paste(
        if (is.na(m)) {
          MARKER_PARENT$no_marker
        } else {
          MARKER_PARENT$by_marker[[m]] %||% MARKER_PARENT$default
        },
        collapse = ";"
      )
    }),
    note = case_when(
      is.na(marker) ~ "no marker in control file",
      marker %in%
        MARKER_PARENT$also_non_t ~ "assigned T cells; also expressed outside T cells",
      .default = ""
    )
  )
write_csv(marker_parent, file.path(RUN_DIR, "marker_parent.csv"))

# Scores per donor: every arm on as188_joint's population rows ---------------------------------------
pooled <- read_csv(POOLED_MANIFEST, show_col_types = FALSE)
gs <- load_gs(file.path(run_in("common-gated"), "gs/as188_joint"))
donors <- sort(pooled$sample[pooled$arm == "as188_joint"])
stopifnot(
  "scored samples are the sample table" = setequal(donors, samples_tbl()$sample)
)
scores <- donors |>
  set_names() |>
  map(\(s) {
    mats <- map(set_names(ARMS), \(a) {
      arm_matrix(a, s, marker_parent$channel, pooled)
    })
    neg_cuts <- map(set_names(ARMS), \(a) {
      cuts_by_channel(a, s, marker_parent$channel, "neg")
    })
    pos_cuts <- map(set_names(ARMS), \(a) {
      cuts_by_channel(a, s, marker_parent$channel, "pos")
    })
    pop_rows <- map(set_names(POPS), \(p) which(gh_pop_get_indices(gs[[s]], p)))
    stopifnot(
      "gate rows match cleaned rows" = length(gh_pop_get_indices(
        gs[[s]],
        "root"
      )) ==
        nrow(mats[[1]])
    )
    imap(mats, \(x, a) {
      map(pop_rows, \(r) {
        score_channels(
          x,
          mats[["as188_joint"]],
          r,
          neg_cut = neg_cuts[[a]],
          pos_cut = pos_cuts[[a]],
          ref_neg_cut = neg_cuts[["as188_joint"]]
        )
      }) |>
        list_rbind(names_to = "pop")
    }) |>
      list_rbind(names_to = "arm")
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

# Donor medians on each channel's designated parents ---------------------------------------------------
designated <- marker_parent |>
  separate_longer_delim(parents, ";") |>
  select(channel, marker, pop = parents)
summary_tbl <- scores |>
  inner_join(
    designated,
    by = c("channel", "pop"),
    relationship = "many-to-one"
  ) |>
  mutate(
    d_squeezed = squeezed_share - squeezed_share[arm == "as188_joint"],
    d_pos = pos_share - pos_share[arm == "as188_joint"],
    .by = c(donor, pop, channel)
  ) |>
  summarise(
    median_d_squeezed = median(d_squeezed),
    min_d_squeezed = min(d_squeezed),
    median_rsd_ratio = median(rsd_ratio),
    max_rsd_ratio = max(rsd_ratio),
    median_d_pos = median(d_pos),
    .by = c(arm, pop, channel, marker, all_of(COND_BY))
  )
write_csv(summary_tbl, file.path(RUN_DIR, "summary.csv"))

# Figures: channel x arm heatmaps, one panel per parent (24H: one row of panels per condition) ----------------
facet_pop <- if (is.null(DS$conditions)) {
  facet_wrap(~pop, nrow = 1, scales = "free_y")
} else {
  facet_grid(rows = vars(condition), cols = vars(pop), scales = "free_y")
}
plot_df <- summary_tbl |>
  filter(arm != "as188_joint") |>
  mutate(
    label = paste0(marker, " (", sub("-A$", "", channel), ")"),
    arm = factor(arm, levels = ARMS[-1]),
    pop = factor(pop, levels = POPS)
  )
p_sq <- ggplot(
  plot_df,
  aes(arm, reorder(label, median_d_squeezed), fill = median_d_squeezed)
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.2f", median_d_squeezed)),
    size = 5,
    size.unit = "pt"
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-1, 1) * max(abs(plot_df$median_d_squeezed)),
    name = "\u0394 squeezed\nshare"
  ) +
  facet_pop +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(5, 5, 16, 5) # room for the rotated labels
  )
save_panel(p_sq, plot_df, RUN_DIR, "fig_screen_squeezed", 183, 183)
p_rsd <- ggplot(
  plot_df,
  aes(arm, reorder(label, -median_rsd_ratio), fill = log2(median_rsd_ratio))
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.2f", median_rsd_ratio)),
    size = 5,
    size.unit = "pt"
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    direction = 1,
    limits = c(-1, 1) * max(abs(log2(plot_df$median_rsd_ratio))),
    name = "log2 rSD\nratio"
  ) +
  facet_pop +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(5, 5, 16, 5) # room for the rotated labels
  )
save_panel(p_rsd, plot_df, RUN_DIR, "fig_screen_rsd", 183, 183)

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
