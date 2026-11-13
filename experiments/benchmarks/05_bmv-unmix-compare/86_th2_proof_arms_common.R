# Purpose: Common-gates copy (plans/common_gates_plan.md, batch 3): gates and rows are the common SpectroFlo labels
#          (clean3), the squeeze cut is refitted as 0.5 x the negative rSD of as188_joint on those cells, and the
#          above-cut readout uses the positive cut (G7a). Negatives for the rSD are below the negative cut (G13); both
#          cuts come from cut_of() per arm, sample and channel (3D: arcsinh 1.5). The Q2 flag rule is unchanged.
#          Needs the Th2 and IL-10 dyes (config targets th2_dye, il10_dye); stops when the dataset lacks them.
#          Original purpose: grill Q3 (user, 2026-11-10): does the AF405 (GLUT1) proof chain move Th2 (eFluor 450) and IL-10
#          (BB700) the way it moved TNF (BV510)? Scores the existing proof arms with the 39b / 51 metrics on
#          as188_joint's CD4T and CD8T rows (helper75), applies the grill Q2 flag rule, and draws Th2 per arm
#          with donors joined across arms (the same events, so paired). No new unmixing.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/86_th2_proof_arms_common.R   (repo root)
# Inputs:  <data_dir>/unmix_cells/<arm>/   proof arms (24H: one subfolder per condition) [60, 61 scripts]; s102_c188
#            and var102_only_<glut1 dye> [06, 45]
#          registry key common-gated: manifest.csv (clean3 rows of every arm) [97c] and gs/as188_joint (CD4T, CD8T
#            gates) [batch 0, 97d]
#          registry key common-cuts, only when the dataset's pos_rule or neg_rule is not fixed: cuts.csv [97e]
#          datasets/<ds>.json (targets, channels, marker_parent) and datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <run_out("th2-proof-common")>/scores.csv   per donor, arm, pop, channel: rSD, squeezed share (24H: + condition)
#          <run_out("th2-proof-common")>/flags.csv    per arm, pop, channel (24H: and condition): Q2 rule components, flag
#          <run_out("th2-proof-common")>/fig_th2_proof_<pop>.{pdf,png,csv}   Th2 squeezed share per
#            arm, donors joined, median crossbar
#          <run_out("th2-proof-common")>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), run_out(), samples_tbl()
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()
source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lib/plots.R")
# shared scoring functions; sourcing it also sources dataset.R again, which is harmless
source("experiments/benchmarks/05_bmv-unmix-compare/helper98_marker_scores.R")

# Parameters ---------------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
# negative cut (rSD, squeeze cut; G13) and positive cut (share above; G7a): cut_of(arm, sample, channel, kind)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
D_SQUEEZE_MIN <- 0.10 # grill Q2 flag rule
RSD_RATIO_MAX <- 0.8 # grill Q2 flag rule
# the Th2 and IL-10 dyes exist in some panels only; a dataset without them has no Th2 proof to score
for (needed in c("th2_dye", "il10_dye")) {
  if (is.null(cfg("targets")[[needed]])) {
    stop(
      sprintf(
        "86_th2_proof_arms_common.R: needs target %s, missing in dataset %s",
        needed,
        DS$dataset
      ),
      call. = FALSE
    )
  }
}
CHANNELS <- c(
  paste0(cfg("targets")$th2_dye, "-A"),
  paste0(cfg("targets")$il10_dye, "-A"),
  cfg("channels")$TNF_slot
) |>
  setNames(c("Th2", "IL-10", cfg("targets")$tnf_label))
POPS <- marker_parent_cfg()$default # T-cell parents of the proof (3D: CD4T, CD8T)
COND_BY <- if (is.null(DS$conditions)) character() else "condition" # 24H only: flags per condition
GLUT1_ARM <- paste0("var102_only_", gsub("[ -]", "_", cfg("targets")$glut1_dye)) # the GLUT1 dye with its 1.0.2 variants
# arms in the order of fig2_interventions (63), with short labels for the x axis
ARMS <- c(
  as188_joint = "1.8.8",
  af405proof_ki_bv510 = "knock-in",
  af405proof_dose_000 = "dose 0",
  af405proof_dose_025 = "dose 0.25",
  af405proof_dose_050 = "dose 0.5",
  setNames("dose 1", GLUT1_ARM),
  af405proof_recentre = "1.0.2 spread on 1.8.8 ref",
  af405proof_ko_bv510 = "knock-out",
  s102_c188 = "all 1.0.2 (target)"
)
OUT_DIR <- run_out("th2-proof-common")

dir.create(OUT_DIR, recursive = TRUE)

# Scores per donor: every arm on as188_joint's population rows ----------------------------------------
pooled <- read_csv(POOLED_MANIFEST, show_col_types = FALSE)
gs <- load_gs(file.path(run_in("common-gated"), "gs/as188_joint"))
donors <- sort(pooled$sample[pooled$arm == "as188_joint"])
stopifnot(
  "scored samples are the sample table" = setequal(donors, samples_tbl()$sample)
)
scores <- donors |>
  set_names() |>
  map(\(s) {
    mats <- map(set_names(names(ARMS)), \(a) {
      arm_matrix(a, s, unname(CHANNELS), pooled)
    })
    neg_cuts <- map(set_names(names(ARMS)), \(a) {
      cuts_by_channel(a, s, unname(CHANNELS), "neg")
    })
    pos_cuts <- map(set_names(names(ARMS)), \(a) {
      cuts_by_channel(a, s, unname(CHANNELS), "pos")
    })
    pop_rows <- map(set_names(POPS), \(p) which(gh_pop_get_indices(gs[[s]], p)))
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
  list_rbind(names_to = "donor") |> # the sample id (3D: the well)
  mutate(marker = names(CHANNELS)[match(channel, CHANNELS)])
if (!is.null(DS$conditions)) {
  scores <- scores |>
    left_join(
      select(samples_tbl(), donor = sample, condition),
      by = "donor",
      relationship = "many-to-one",
      unmatched = "error"
    )
}
write_csv(scores, file.path(OUT_DIR, "scores.csv"))

# Q2 flag rule per arm, population and channel ---------------------------------------------------------
ref <- scores |>
  filter(arm == "as188_joint") |>
  select(donor, pop, channel, ref_squeezed = squeezed_share)
per_donor <- scores |>
  left_join(
    ref,
    by = join_by(donor, pop, channel),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(d_squeezed = squeezed_share - ref_squeezed)
flags <- per_donor |>
  filter(arm != "as188_joint") |>
  summarise(
    median_squeezed = median(squeezed_share),
    median_d_squeezed = median(d_squeezed),
    median_rsd_ratio = median(rsd_ratio),
    all_donors_up = all(d_squeezed > 0) && all(rsd_ratio < 1),
    .by = c(arm, pop, marker, all_of(COND_BY))
  ) |>
  mutate(
    flag = median_d_squeezed >= D_SQUEEZE_MIN &
      median_rsd_ratio <= RSD_RATIO_MAX &
      all_donors_up
  )
write_csv(flags, file.path(OUT_DIR, "flags.csv"))

# Figure: Th2 squeezed share per arm, donors joined across arms, one panel per population ---------------
arm_pal <- categorical_palette(unname(ARMS))
for (which_pop in POPS) {
  d <- per_donor |>
    filter(marker == "Th2", pop == which_pop) |>
    mutate(arm_label = unname(ARMS[arm])) |>
    select(donor, arm, arm_label, squeezed_share, rsd_ratio, all_of(COND_BY))
  p <- plot_box(
    d,
    group = "arm_label",
    value = "squeezed_share",
    pal = arm_pal,
    jitter_width = 0, # unused: paired donors are joined by lines, not jittered
    seed = 1,
    value_label = "Th2 cells within 0.5 rSD of 0\n(fraction of parent)",
    pair = "donor",
    title = sprintf("Th2 squeeze, AF405 proof arms, %s", which_pop)
  ) +
    rotate_x()
  # 24H: one panel per condition
  if (!is.null(DS$conditions)) {
    p <- p + facet_wrap(vars(condition))
  }
  save_panel(p, d, OUT_DIR, sprintf("fig_th2_proof_%s", which_pop), 89, 89)
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
