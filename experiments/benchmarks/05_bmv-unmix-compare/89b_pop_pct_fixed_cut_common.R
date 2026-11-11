# Purpose: copy of 89b_pop_pct_fixed_cut.R on the common SpectroFlo gates (plans/common_gates_plan.md, G2, G7a,
#          G12). The common GatingSets already hold TNFpos and IFNgpos as the fixed arcsinh 1.5 cut, so 89
#          and this script agree; this script recomputes TNF+ and IFNg+ % of CD4T and CD8T from the events at
#          CUT and stops if they differ from the node counts. Lineage gates (nonDebris to DP/DN) are the common
#          SpectroFlo labels in both arms. Descriptive, no winner picked. Figure in the lucid-coder figure style.
# Run:     sbatch experiments/slurm/run_script_lib.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/89b_pop_pct_fixed_cut_common.R \
#            experiments/runs/2026-11-11_bmv-common-downstream-opt188-v2   (repo root; run dir must not exist)
# Inputs:  experiments/runs/2026-11-11_bmv-common-gated/gs/{as188_joint,opt188_pregate_joint_p1af1}
#            common-gate GatingSets, channels arcsinh (cofactor 6000) [batch 0], README.md
# Outputs: <run_dir>/pop_pct.csv              arm, donor, population, count, parent count, % of parent
#          <run_dir>/fig_pop_pct.{pdf,png,csv} population % of parent, donors joined by lines
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
COMMON <- "experiments/runs/2026-11-11_bmv-common-gated" # batch 0; gs/<arm> for all arms
ARMS <- c(as188_joint = "as188_joint", opt188 = "opt188_pregate_joint_p1af1") # grill Q1 (user, 2026-11-10)
CUT <- 1.5 # cytokine+ above arcsinh 1.5, cofactor 6000 (G7a, user 2026-11-11; was 1)
CYTOKINE_CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A")
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

# Cytokine+ at the fixed cut: events above arcsinh 1.5 in the common CD4T and CD8T -------------------
cytokine_pct <- gs_list |>
  map(\(gs) {
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
        count = map_int(x, \(v) sum(v > CUT)),
        parent_count = lengths(x),
        population = paste(subset, paste0(cytokine, "+ (fixed cut)"))
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
          paste0(sub("pos$", "", basename(Population)), "+ (fixed cut)")
        ),
        node_count = Count
      )
  }) |>
  list_rbind(names_to = "arm")
stopifnot(
  "event counts at CUT equal the common TNFpos and IFNgpos node counts" = nrow(
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
    ~population,
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
