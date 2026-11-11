# Purpose: fix to 89_downstream_as188_vs_opt188.R part (1) (coordinator, 2026-11-10). The TNFpos and
#          IFNgpos nodes in the gated GatingSets are data-driven gates fitted per arm (19a_gate_arm.R:
#          mindensity for TNF, tailgate for IFNg), a recorded dead end (HANDOFF "Dead ends": TNF has
#          no valley). This script keeps the lineage gates (nonDebris to DP/DN) as gated, drops the
#          cytokine nodes, and computes TNF+ and IFNg+ % of CD4T and CD8T at the settled fixed cut
#          (arcsinh 1, cofactor 6000; README decision table), per donor, each arm on its own gates.
#          Descriptive, no winner picked. Figure in the lucid-coder figure style.
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/89b_pop_pct_fixed_cut.R \
#            experiments/runs/2026-11-10_bmv-downstream-opt188-v2   (repo root; run dir must not exist)
# Inputs:  experiments/runs/2026-11-10_bmv-ablation-gated/gs/{as188_joint,opt188_pregate_joint_p1af1}
#            GatingSets, channels already arcsinh (cofactor 6000) [19a_gate_arm.R via 51_score_ablation.R]
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
GATED <- "experiments/runs/2026-11-10_bmv-ablation-gated/gs"
ARMS <- c(as188_joint = "as188_joint", opt188 = "opt188_pregate_joint_p1af1") # grill Q1 (user, 2026-11-10)
CUT <- 1 # cytokine+ above arcsinh 1, cofactor 6000 (README decision table, user, 2026-11-10)
CYTOKINE_CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A")
T_SUBSETS <- c("CD4T", "CD8T")
STYLE <- "experiments/benchmarks/05_bmv-unmix-compare/style"

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)
source(file.path(STYLE, "theme.R"))
source(file.path(STYLE, "palettes.R"))
source(file.path(STYLE, "plots.R"))
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

# Lineage gates: % of parent from each arm's own gates, cytokine nodes dropped -----------------------
gs_list <- map(ARMS, \(a) load_gs(file.path(GATED, a)))
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

# Cytokine+ at the fixed cut: events above arcsinh 1 in each arm's own CD4T and CD8T -------------------
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
