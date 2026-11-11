# Purpose: Common-gates copy (plans/common_gates_plan.md, batch 3): gates and rows are the common SpectroFlo labels
#          (clean3), the squeeze cut is refitted as 0.5 x the negative rSD of as188_joint on those cells, and the
#          above-cut readout uses arcsinh 1.5 (G7a). Negatives for the rSD are below arcsinh 1.5 (G13).
#          Test 1 (README decision table, user 2026-11-10): check every fluorescence channel for the TNF-type
#          squeeze. For each arm, donor and channel: squeezed share near 0, negative rSD ratio vs as188_joint and
#          share above arcsinh 1.5 (39b / 51 metrics, cofactor 6000 in memory), on as188_joint's gates (rows paired
#          across arms). Parents: CD4T and CD8T for T-cell markers, live cells otherwise (marker_parent.csv).
#          51 has no channel-level "affected" rule (its labels score arms by the TNF gap), so no flag is set.
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/75_all_marker_screen_common.R \
#            experiments/runs/2026-11-11_bmv-all-marker-screen-common      (repo root; run dir must not exist)
# Inputs:  experiments/runs/2026-11-11_bmv-common-gated/manifest.csv   clean3 FCS of every arm [97c]
#          experiments/runs/2026-11-11_bmv-common-gated/gs/as188_joint   common-gate GatingSet [batch 0, 97d]
#          experiments/bmv_test/compare_3D/data/unmix_cells/{swap_var102,swap_af102}  [37_swap_af_variants.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv
#            dye to marker map [user zip]
# Outputs: <run_dir>/marker_parent.csv   channel, dye, marker, parent populations, ambiguity note
#          <run_dir>/scores.csv          per arm, donor, population, channel: readouts (all three populations)
#          <run_dir>/summary.csv         per arm, population, channel (designated parents): donor medians of the
#                                        difference vs as188_joint (squeezed, above 1.5) and of the rSD ratio
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

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
ARMS <- c(
  "as188_joint",
  "as102_var",
  "s102_c188",
  "s188_c102",
  "swap_var102",
  "swap_af102"
)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NEG_CUT <- 1.5 # negatives (rSD, squeeze cut) are events below arcsinh 1.5, as the positivity cut (G13, user 2026-11-11; was 1)
POS_CUT <- 1.5 # positives are above arcsinh 1.5 (G7a, user 2026-11-11; was 1)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
POPS <- c("CD4T", "CD8T", "live")
# markers scored in live cells, not T cells: the viability dye defines the live gate; CD56 marks NK cells;
# PE has no marker in the control file. Every other marker is scored in CD4T and CD8T (T-cell panel).
LIVE_ONLY <- c("Viability", "CD56", NA)
# assigned to T cells but also expressed outside them (listed in the report as ambiguous)
ALSO_NON_T <- c(
  "CD45RA",
  "CD7",
  "CD161",
  "CD57",
  "CD38",
  "CD71",
  "CCR7",
  "CD127",
  "CD25"
)
CONTROL_FILE <- "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()
source("experiments/benchmarks/05_bmv-unmix-compare/helper98_marker_scores.R")

# Channels and parents from the control file ----------------------------------------------------------
marker_parent <- read_csv(CONTROL_FILE, show_col_types = FALSE) |>
  filter(!fluorophore %in% c("AF", "Negative")) |>
  transmute(
    channel = paste0(fluorophore, "-A"),
    dye = fluorophore,
    marker,
    parents = if_else(marker %in% LIVE_ONLY, "live", "CD4T;CD8T"),
    note = case_when(
      is.na(marker) ~ "no marker in control file",
      marker %in%
        ALSO_NON_T ~ "assigned T cells; also expressed outside T cells",
      .default = ""
    )
  )
write_csv(marker_parent, file.path(RUN_DIR, "marker_parent.csv"))

# Scores per donor: every arm on as188_joint's population rows ---------------------------------------
pooled <- read_csv(POOLED_MANIFEST, show_col_types = FALSE)
gs <- load_gs("experiments/runs/2026-11-11_bmv-common-gated/gs/as188_joint")
donors <- sort(pooled$sample[pooled$arm == "as188_joint"])
scores <- donors |>
  set_names() |>
  map(\(s) {
    mats <- map(set_names(ARMS), \(a) {
      arm_matrix(a, s, marker_parent$channel, pooled)
    })
    pop_rows <- map(set_names(POPS), \(p) which(gh_pop_get_indices(gs[[s]], p)))
    stopifnot(
      "gate rows match cleaned rows" = length(gh_pop_get_indices(
        gs[[s]],
        "root"
      )) ==
        nrow(mats[[1]])
    )
    map(mats, \(x) {
      map(pop_rows, \(r) score_channels(x, mats[["as188_joint"]], r)) |>
        list_rbind(names_to = "pop")
    }) |>
      list_rbind(names_to = "arm")
  }) |>
  list_rbind(names_to = "donor")
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
    .by = c(arm, pop, channel, marker)
  )
write_csv(summary_tbl, file.path(RUN_DIR, "summary.csv"))

# Figures: channel x arm heatmaps, one panel per parent ---------------------------------------------------
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
  facet_wrap(~pop, nrow = 1, scales = "free_y") +
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
  facet_wrap(~pop, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(5, 5, 16, 5) # room for the rotated labels
  )
save_panel(p_rsd, plot_df, RUN_DIR, "fig_screen_rsd", 183, 183)

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
