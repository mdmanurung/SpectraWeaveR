# Purpose: copy of 89b_pop_pct_fixed_cut.R on the common SpectroFlo gates (plans/common_gates_plan.md, G2, G7a,
#          G12). The common GatingSets already hold TNFpos and IFNgpos as the fixed arcsinh 1.5 cut, so 89
#          and this script agree; this script recomputes TNF+ and IFNg+ % of CD4T and CD8T from the events at
#          the positive cut and stops if they differ from the node counts. Lineage gates (nonDebris to DP/DN) are the common
#          SpectroFlo labels in both arms. Descriptive, no winner picked. Figure in the lucid-coder figure style.
#          The two cytokine slots, their labels and the positive cut come from the dataset config.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> \
#          sbatch experiments/slurm/run_script_lib.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/89b_pop_pct_fixed_cut_common.R \
#            [<run_dir>]   (repo root; default run_out("common-downstream-opt188-v2"); run dir must not exist)
# Inputs:  run_in("common-gated")/gs/{as188_joint,opt188_pregate_joint_p1af1}
#            common-gate GatingSets, channels arcsinh (cofactor 6000) [batch 0], README.md
#          cut_of(arm, sample, channel, "pos")  positive cut [config pos_rule; cuts.csv from 97e when not fixed]
#          experiments/benchmarks/05_bmv-unmix-compare/datasets/  config, run registry, sample table
# Outputs: <run_dir>/pop_pct.csv              arm, donor (= sample id), population, count, parent count, % of
#                                             parent; plus condition when the dataset has conditions
#          <run_dir>/fig_pop_pct.{pdf,png,csv} population % of parent, donors joined by lines (faceted by
#                                             condition when the dataset has conditions)
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
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), run_out(), cut_of()

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(RUN_DIR)) {
  RUN_DIR <- run_out("common-downstream-opt188-v2")
}
COMMON <- run_in("common-gated") # batch 0; gs/<arm> for all arms
ARMS <- c(as188_joint = "as188_joint", opt188 = "opt188_pregate_joint_p1af1") # grill Q1 (user, 2026-11-10)
# These T-cell nodes score TNF and IFNg; a dataset that labels the slots otherwise has no such target
if (
  !identical(
    c(cfg("targets")$tnf_label, cfg("targets")$ifng_label),
    c("TNF", "IFNg")
  )
) {
  stop(
    "89b_pop_pct_fixed_cut_common.R: needs TNF and IFNg, missing in dataset ",
    DS$dataset,
    " (slots are labelled ",
    cfg("targets")$tnf_label,
    " and ",
    cfg("targets")$ifng_label,
    ")",
    call. = FALSE
  )
}
# Cytokine+ above the dataset's positive cut, cofactor 6000 (3D: fixed arcsinh 1.5, G7a, user 2026-11-11).
# Slot names are the dataset's marker labels (3D: TNF, IFNg).
CYTOKINE_CHANNELS <- set_names(
  c(cfg("channels")$TNF_slot, cfg("channels")$IFNg_slot),
  c(cfg("targets")$tnf_label, cfg("targets")$ifng_label)
)
CUT_NOTE <- if (identical(cfg("pos_rule")$type, "fixed")) {
  "fixed cut"
} else {
  "sample cut"
} # label suffix
T_SUBSETS <- c("CD4T", "CD8T")
LIB <- Sys.getenv("SW_LIB", "experiments/benchmarks/05_bmv-unmix-compare/lib")

stopifnot(
  "batch 0 finished (README.md in the common dir)" = file.exists(file.path(
    COMMON,
    "README.md"
  )),
  "run dir must not exist" = !dir.exists(RUN_DIR)
)
dir.create(RUN_DIR, recursive = TRUE)
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))
source(file.path(LIB, "plots.R"))
# same palette and shapes as script 89, keyed by label
pal_arm <- categorical_palette(c(names(ARMS), "as188_joint, other seed"))
shape_arm <- c(as188_joint = 16, opt188 = 17, "as188_joint, other seed" = 1)
arm_n_labels <- \(donor) {
  \(g) {
    sprintf(
      "%s\n(n = %d)",
      c(as188_joint = "as188", opt188 = "opt188")[g],
      length(unique(donor))
    )
  }
}

# Lineage gates: % of parent from the common labels, cytokine nodes dropped -----------------------
gs_list <- map(ARMS, \(a) load_gs(file.path(COMMON, "gs", a)))
stopifnot(
  "both arms have the same gate tree" = identical(
    gs_get_pop_paths(gs_list$as188_joint),
    gs_get_pop_paths(gs_list$opt188)
  )
)
lineage_pct <- gs_list |>
  map(\(gs) {
    gs_pop_get_count_fast(gs) |>
      as_tibble() |>
      filter(!basename(Population) %in% c("TNFpos", "IFNgpos")) |>
      transmute(
        donor = name,
        population = basename(Population),
        count = Count,
        parent_count = ParentCount
      )
  }) |>
  list_rbind(names_to = "arm")

# Cytokine+ at the positive cut: events above it in the common CD4T and CD8T -------------------
cytokine_pct <- gs_list |>
  imap(\(gs, arm_name) {
    expand_grid(
      donor = sampleNames(gs),
      subset = T_SUBSETS,
      cytokine = names(CYTOKINE_CHANNELS)
    ) |>
      mutate(
        x = pmap(list(donor, subset, cytokine), \(s, p, k) {
          exprs(cytoframe_to_flowFrame(gh_pop_get_data(
            gs[[s]],
            p
          )))[, CYTOKINE_CHANNELS[[k]]]
        }),
        count = pmap_int(list(x, donor, cytokine), \(v, s, k) {
          sum(v > cut_of(ARMS[[arm_name]], s, CYTOKINE_CHANNELS[[k]], "pos"))
        }),
        parent_count = lengths(x),
        population = paste(subset, paste0(cytokine, "+ (", CUT_NOTE, ")"))
      ) |>
      select(donor, population, count, parent_count)
  }) |>
  list_rbind(names_to = "arm")

# the cytokine nodes of the common GatingSets are the same cut, so event counts must equal node counts
node_counts <- gs_list |>
  map(\(gs) {
    gs_pop_get_count_fast(gs) |>
      as_tibble() |>
      filter(basename(Population) %in% c("TNFpos", "IFNgpos")) |>
      transmute(
        donor = name,
        population = paste(
          basename(Parent),
          paste0(
            unname(
              c(
                TNF = names(CYTOKINE_CHANNELS)[[1]],
                IFNg = names(CYTOKINE_CHANNELS)[[2]]
              )[
                sub("pos$", "", basename(Population))
              ]
            ),
            "+ (",
            CUT_NOTE,
            ")"
          )
        ),
        node_count = Count
      )
  }) |>
  list_rbind(names_to = "arm")
stopifnot(
  "event counts at the positive cut equal the common TNFpos and IFNgpos node counts" = nrow(
    inner_join(
      cytokine_pct,
      node_counts,
      by = join_by(arm, donor, population),
      relationship = "one-to-one",
      unmatched = "error"
    ) |>
      filter(count != node_count)
  ) ==
    0
)

pop_pct <- bind_rows(lineage_pct, cytokine_pct) |>
  mutate(pct_of_parent = 100 * count / parent_count)
stopifnot(
  "gated samples equal the sample table" = setequal(
    pop_pct$donor,
    samples_tbl()$sample
  )
)
# Condition of each sample (24H only); 3D and EXV have none, so their tables stay as they are
if (!is.null(DS$conditions)) {
  pop_pct <- pop_pct |>
    left_join(
      select(samples_tbl(), donor = sample, condition),
      by = "donor",
      relationship = "many-to-one",
      unmatched = "error"
    )
}
write_csv(pop_pct, file.path(RUN_DIR, "pop_pct.csv"))

# Figure: population % of parent, donors joined by lines -------------------------------------------
pop_df <- pop_pct |>
  mutate(
    arm = factor(arm, levels = names(ARMS)),
    population = factor(population, levels = unique(population)) # gate-tree order, cytokines last
  )
p_pop <- ggplot(pop_df, aes(arm, pct_of_parent)) +
  geom_line(aes(group = donor), colour = "grey60") +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.5,
    colour = "black"
  ) +
  geom_point(aes(colour = arm, shape = arm), size = 1.5) +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_shape_manual(values = shape_arm, guide = "none") +
  scale_x_discrete(labels = arm_n_labels(pop_df$donor)) +
  facet_wrap(
    if (is.null(DS$conditions)) {
      vars(population)
    } else {
      vars(population, condition)
    },
    scales = "free_y",
    ncol = 4,
    strip.position = "left"
  ) +
  labs(x = NULL, y = NULL, title = "Gated populations, % of parent") +
  theme(strip.placement = "outside", strip.background = element_blank())
save_panel(
  p_pop,
  pop_df,
  RUN_DIR,
  "fig_pop_pct",
  width_mm = 183,
  height_mm = 130
)

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
