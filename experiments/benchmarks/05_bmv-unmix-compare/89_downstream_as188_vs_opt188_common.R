# Purpose: copy of 89_downstream_as188_vs_opt188.R on the common SpectroFlo gates (plans/common_gates_plan.md,
#          G2, G3, G7a, G12). grill Q8/Q11 (README rows "Downstream comparison" and "Affected-marker flag",
#          2026-11-10): compare as188_joint and opt188 (opt188_pregate_joint_p1af1) downstream, descriptively,
#          no winner picked. (1) gated population % of parent per donor; both arms carry the same SpectroFlo
#          labels, and TNFpos and IFNgpos are the fixed arcsinh 1.5 cut (the old run used data-driven gates);
#          (2) the grill Q2 flag for every marker (opt188 vs as188_joint; CD4T and CD8T, live cells for non-T
#          markers), on the common labels with rows paired across arms (helper75); (3) SOM clustering with the
#          P10 settings and lineage markers on the common live events (as 03_downstream.R): ARI per donor, an ARI noise
#          floor from as188_joint reclustered with seeds 1-3, and metacluster abundance with each as188_joint
#          metacluster matched to the opt188 metacluster it overlaps most.
#          Figures follow the lucid-coder figure style (style/theme.R, palettes.R, plots.R).
#          Dataset: SW_DATASET (3d, 24h, exv) picks the config; with conditions (24H) the marker flags and the
#          metacluster abundance are summarised per condition.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> \
#          sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/89_downstream_as188_vs_opt188_common.R \
#            [<run_dir>]      (repo root; default run_out("common-downstream-opt188"); run dir must not exist)
# Inputs:  run_in("common-gated")/gs/{as188_joint,opt188_pregate_joint_p1af1}
#            common-gate GatingSets [batch 0], README.md
#          run_in("clean3")/manifest.csv   clean3 FCS of every arm [97c_reclean_pooled.R]
#          data_path("autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv")
#            dye to marker map [user zip]
#          experiments/benchmarks/05_bmv-unmix-compare/datasets/  config (channels, targets, neg_rule, pos_rule,
#            cluster_markers_other, live_only_markers), run registry, sample table
# Outputs: <run_dir>/pop_pct.csv              arm, donor (= sample id), population, count, parent count, % of
#                                             parent (+ condition when the dataset has conditions)
#          <run_dir>/marker_scores.csv        arm, donor, population, channel: 39b/51 readouts
#          <run_dir>/marker_flags.csv         per population and marker (and condition): donor medians, Q2 flag
#          <run_dir>/cluster_assignments.rds  donor, event row, metacluster per arm and noise seed
#          <run_dir>/clusters_ari.csv         comparison (opt188 or noise seed), donor, ARI vs as188_joint seed 42
#          <run_dir>/cluster_match.csv        as188_joint metacluster, matched opt188 metacluster, Jaccard
#          <run_dir>/cluster_abundance.csv    arm, donor, matched cluster, % of clustered events
#          <run_dir>/fig_pop_pct.{pdf,png,csv}           population % of parent, donors joined by lines
#          <run_dir>/fig_marker_flags.{pdf,png,csv}      donor-median squeezed-share change, flags marked
#          <run_dir>/fig_cluster_abundance.{pdf,png,csv} matched metacluster %, donors joined by lines
#          <run_dir>/fig_ari.{pdf,png,csv}               ARI per donor, opt188 and the seed noise floor
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
pkgload::load_all(".", export_all = FALSE, quiet = TRUE) # sw_cluster_run; run from the repo root

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(RUN_DIR)) {
  RUN_DIR <- run_out("common-downstream-opt188")
}
COMMON <- run_in("common-gated") # batch 0; gs/<arm> for all arms
CLEAN3_MANIFEST <- file.path(
  run_in(if (identical(cfg("dataset"), "3d")) "clean3" else "common-gated"),
  "manifest.csv"
) # the rows the common gs hold (3D: clean3 files; 24H and ex vivo: the uncleaned files of 97e)
ARMS <- c(as188_joint = "as188_joint", opt188 = "opt188_pregate_joint_p1af1") # grill Q1 (user, 2026-11-10)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
# helper75_marker-scores.R reads one global scalar CUT (arcsinh) for the negative rSD and the positive share,
# so only a dataset with the same fixed negative and positive cut can run here (3D: 1.5, G7a, user
# 2026-11-11; was 1). A per-sample cut needs helper75 to call cut_of() per arm, sample and channel.
if (
  !(identical(cfg("neg_rule")$type, "fixed") &&
    identical(cfg("pos_rule")$type, "fixed") &&
    identical(cfg("neg_rule")$value, cfg("pos_rule")$value))
) {
  stop(
    "89_downstream_as188_vs_opt188_common.R: helper75_marker-scores.R takes one scalar CUT; dataset ",
    DS$dataset,
    " has no single fixed cut for negatives and positives",
    call. = FALSE
  )
}
CUT <- cfg("neg_rule")$value
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
FLAG_D_SQUEEZED <- 0.10 # grill Q2, revised 2026-11-10 (README): donor-median squeezed-share rise at or above this, same sign in all donors; rSD criterion dropped
LIVE_ONLY <- cfg("live_only_markers") # scored in live cells, as in 75_all_marker_screen.R (3D: Viability, CD56, NA)
XDIM <- 10 # README "P10 clustering": SOM grid
YDIM <- 10
N_META <- 20 # P10: metaclusters
SOM_SEED <- 42 # P10
NOISE_SEEDS <- c(1, 2, 3) # ARI noise floor, as in 03_downstream.R (user, grill Q2 and Q4, 2026-11-09)
# README "Cluster markers": lineage only (user, 2026-11-09), channels as in 03_downstream.R
CLUSTER_MARKERS <- c(
  CD3 = cfg("channels")$CD3,
  CD4 = cfg("channels")$CD4,
  CD8 = cfg("channels")$CD8,
  unlist(cfg("cluster_markers_other")) # the other lineage markers, channel by marker (3D: CD45RA to CXCR5)
)
FLAG_LIMIT <- 0.25 # display only: symmetric limit of the diverging squeezed-share scale
CONTROL_FILE <- data_path(
  "autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"
)
COND_KEY <- if (is.null(DS$conditions)) character() else "condition" # extra summary key (24H)
LIB <- Sys.getenv("SW_LIB", "experiments/benchmarks/05_bmv-unmix-compare/lib")

stopifnot(
  "batch 0 finished (README.md in the common dir)" = file.exists(file.path(
    COMMON,
    "README.md"
  )),
  "run dir must not exist" = !dir.exists(RUN_DIR)
)
dir.create(RUN_DIR, recursive = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/helper75_marker-scores.R")
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))
source(file.path(LIB, "plots.R"))
# Condition of each sample (24H); 3D and EXV have none, so their tables stay as they are
with_condition <- function(d) {
  if (is.null(DS$conditions)) {
    return(d)
  }
  left_join(
    d,
    select(samples_tbl(), donor = sample, condition),
    by = "donor",
    relationship = "many-to-one",
    unmatched = "error"
  )
}
# one palette and shape set for every figure, keyed by label: the two arms and the seed noise floor
pal_arm <- categorical_palette(c(names(ARMS), "as188_joint, other seed"))
shape_arm <- c(as188_joint = 16, opt188 = 17, "as188_joint, other seed" = 1) # encode twice
# short x labels over n donors, so the two arm labels do not collide in narrow panels
arm_n_labels <- \(donor) {
  \(g) {
    sprintf(
      "%s\n(n = %d)",
      c(as188_joint = "as188", opt188 = "opt188")[g],
      length(unique(donor))
    )
  }
}

# (1) Population % of parent, common labels in both arms --------------------------------------------
gs_list <- map(ARMS, \(a) load_gs(file.path(COMMON, "gs", a)))
stopifnot(
  "both arms have the same gate tree" = identical(
    gs_get_pop_paths(gs_list$as188_joint),
    gs_get_pop_paths(gs_list$opt188)
  )
)
pop_pct <- gs_list |>
  map(\(gs) {
    gs_pop_get_count_fast(gs) |>
      as_tibble() |>
      transmute(
        donor = name,
        # cytokine gates repeat under CD4T and CD8T, so they carry their parent's name
        population = if_else(
          basename(Population) %in% c("TNFpos", "IFNgpos"),
          paste(basename(Parent), basename(Population)),
          basename(Population)
        ),
        parent = basename(Parent),
        count = Count,
        parent_count = ParentCount,
        pct_of_parent = 100 * Count / ParentCount
      )
  }) |>
  list_rbind(names_to = "arm") |>
  with_condition()
write_csv(pop_pct, file.path(RUN_DIR, "pop_pct.csv"))

# (2) Q2 flags: every marker, opt188 vs as188_joint on the common labels (paired rows) ------------
marker_parent <- read_csv(CONTROL_FILE, show_col_types = FALSE) |>
  filter(!fluorophore %in% c("AF", "Negative")) |>
  transmute(
    channel = paste0(fluorophore, "-A"),
    marker,
    parents = if_else(marker %in% LIVE_ONLY, "live", "CD4T;CD8T")
  )
pooled <- read_csv(CLEAN3_MANIFEST, show_col_types = FALSE)
donors <- sort(pooled$sample[pooled$arm == "as188_joint"])
stopifnot(
  "clean3 samples equal the sample table" = setequal(
    donors,
    samples_tbl()$sample
  )
)
gs_ref <- gs_list$as188_joint
mats <- donors |>
  set_names() |>
  map(\(s) map(ARMS, \(a) arm_matrix(a, s, marker_parent$channel, pooled)))
marker_scores <- donors |>
  set_names() |>
  map(\(s) {
    stopifnot(
      "gate rows match cleaned rows" = length(gh_pop_get_indices(
        gs_ref[[s]],
        "root"
      )) ==
        nrow(mats[[s]]$as188_joint)
    )
    pop_rows <- map(set_names(c("CD4T", "CD8T", "live")), \(p) {
      which(gh_pop_get_indices(gs_ref[[s]], p))
    })
    map(mats[[s]], \(x) {
      map(pop_rows, \(r) score_channels(x, mats[[s]]$as188_joint, r)) |>
        list_rbind(names_to = "pop")
    }) |>
      list_rbind(names_to = "arm")
  }) |>
  list_rbind(names_to = "donor")
write_csv(marker_scores, file.path(RUN_DIR, "marker_scores.csv"))

# Per donor: opt188 minus as188_joint on each marker's designated parents, then the Q2 rule
designated <- marker_parent |>
  separate_longer_delim(parents, ";") |>
  select(channel, marker, pop = parents)
per_donor <- marker_scores |>
  inner_join(
    designated,
    by = join_by(channel, pop),
    relationship = "many-to-one"
  ) |>
  select(donor, pop, channel, marker, arm, squeezed_share, rsd_ratio) |>
  pivot_wider(names_from = arm, values_from = c(squeezed_share, rsd_ratio)) |>
  mutate(
    d_squeezed = squeezed_share_opt188 - squeezed_share_as188_joint,
    rsd_ratio = rsd_ratio_opt188
  ) |>
  with_condition()
marker_flags <- per_donor |>
  summarise(
    median_d_squeezed = median(d_squeezed),
    median_rsd_ratio = median(rsd_ratio),
    all_donors_up = all(d_squeezed > 0),
    all_donors_narrower = all(rsd_ratio < 1),
    .by = all_of(c("pop", "channel", "marker", COND_KEY))
  ) |>
  mutate(
    flag = median_d_squeezed >= FLAG_D_SQUEEZED & all_donors_up
  )
write_csv(marker_flags, file.path(RUN_DIR, "marker_flags.csv"))

# (3) Clustering on the same paired live events: the common live label --------------------------
live_rows <- map(set_names(donors), \(s) {
  which(gh_pop_get_indices(gs_ref[[s]], "live"))
})
events <- tibble(
  donor = rep(donors, lengths(live_rows)),
  event = unlist(live_rows, use.names = FALSE)
)
arm_mats <- ARMS |>
  names() |>
  set_names() |>
  map(\(a) {
    map(donors, \(s) {
      m <- arm_matrix(ARMS[[a]], s, unname(CLUSTER_MARKERS), pooled)
      m[live_rows[[s]], , drop = FALSE]
    }) |>
      do.call(what = rbind) |>
      (\(m) asinh(m / COFACTOR))()
  })
cluster_som <- \(m, seed) {
  sw_cluster_run(
    m,
    lineage_markers = CLUSTER_MARKERS,
    method = "som",
    xdim = XDIM,
    ydim = YDIM,
    n_metaclusters = N_META,
    seed = seed
  )$assignments
}
labels <- c(
  map(arm_mats, \(m) cluster_som(m, SOM_SEED)),
  map(set_names(NOISE_SEEDS, paste0("as188_joint_seed", NOISE_SEEDS)), \(sd) {
    cluster_som(arm_mats$as188_joint, sd)
  })
)
assignments <- events |>
  bind_cols(as_tibble(labels))
saveRDS(assignments, file.path(RUN_DIR, "cluster_assignments.rds"))

# Adjusted Rand index (contingency formula), as in 03_downstream.R
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  choose2 <- \(x) x * (x - 1) / 2
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  expected <- sum_rows * sum_cols / choose2(length(a))
  (sum_cells - expected) / ((sum_rows + sum_cols) / 2 - expected)
}
compare_cols <- setdiff(names(labels), "as188_joint")
clusters_ari <- assignments |>
  pivot_longer(
    all_of(compare_cols),
    names_to = "comparison",
    values_to = "metacluster"
  ) |>
  summarise(
    ari = adjusted_rand_index(metacluster, as188_joint),
    .by = c(comparison, donor)
  )
write_csv(clusters_ari, file.path(RUN_DIR, "clusters_ari.csv"))

# Match each as188_joint metacluster to the opt188 metacluster it overlaps most (Jaccard)
overlap <- assignments |>
  count(ref = as188_joint, opt = opt188, name = "n_both") |>
  left_join(
    count(assignments, ref = as188_joint, name = "n_ref"),
    by = "ref",
    relationship = "many-to-one"
  ) |>
  left_join(
    count(assignments, opt = opt188, name = "n_opt"),
    by = "opt",
    relationship = "many-to-one"
  ) |>
  mutate(jaccard = n_both / (n_ref + n_opt - n_both))
cluster_match <- overlap |>
  slice_max(jaccard, n = 1, with_ties = FALSE, by = ref) |>
  select(ref, opt, jaccard)
write_csv(cluster_match, file.path(RUN_DIR, "cluster_match.csv"))
n_live <- tibble(donor = donors, n_live = lengths(live_rows))
cluster_abundance <- bind_rows(
  assignments |>
    count(donor, ref = as188_joint) |>
    mutate(arm = "as188_joint"),
  assignments |>
    count(donor, opt = opt188) |>
    inner_join(cluster_match, by = "opt", relationship = "many-to-many") |>
    select(donor, ref, n) |>
    mutate(arm = "opt188")
) |>
  left_join(
    n_live,
    by = "donor",
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(pct = 100 * n / n_live) |>
  with_condition()
write_csv(cluster_abundance, file.path(RUN_DIR, "cluster_abundance.csv"))

# Figure 1: population % of parent, donors joined by lines ------------------------------------------
pop_df <- pop_pct |>
  mutate(
    arm = factor(arm, levels = names(ARMS)),
    population = factor(population, levels = unique(population)) # gate-tree order
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
    vars(population, !!!rlang::syms(COND_KEY)),
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

# Figure 2: marker flags, donor-median squeezed-share change (opt188 - as188_joint) -------------------
flag_df <- marker_flags |>
  mutate(
    label = paste0(marker, " (", sub("-A$", "", channel), ")"),
    text = if_else(flag, "flag", sprintf("%.2f", median_d_squeezed))
  )
HEAT_COL <- "pop" # one column per population; per population and condition in 24H
if (!is.null(DS$conditions)) {
  flag_df <- mutate(flag_df, pop_condition = paste(pop, condition))
  HEAT_COL <- "pop_condition"
}
p_flag <- plot_heatmap(
  flag_df,
  row = "label",
  col = HEAT_COL,
  value = "median_d_squeezed",
  limit = FLAG_LIMIT,
  title = "Squeezed share change",
  x_label = NULL,
  y_label = NULL,
  fill_label = "opt188 minus\nas188, donor\nmedian"
) +
  geom_text(aes(label = text), size = 6, size.unit = "pt") +
  theme(legend.position = "right")
save_panel(
  p_flag,
  flag_df,
  RUN_DIR,
  "fig_marker_flags",
  width_mm = 89,
  height_mm = 140
)

# Figure 3: metacluster abundance, as188_joint metacluster vs its matched opt188 metacluster ----------
abund_df <- cluster_abundance |>
  left_join(cluster_match, by = "ref", relationship = "many-to-one") |>
  mutate(
    arm = factor(arm, levels = names(ARMS)),
    cluster = factor(sprintf("MC%02d (J = %.2f)", ref, jaccard))
  )
p_abund <- ggplot(abund_df, aes(arm, pct)) +
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
  scale_x_discrete(labels = arm_n_labels(abund_df$donor)) +
  facet_wrap(
    vars(cluster, !!!rlang::syms(COND_KEY)),
    scales = "free_y",
    ncol = 4,
    strip.position = "left"
  ) +
  labs(
    x = NULL,
    y = "% of live events",
    title = "Metacluster abundance, matched clusters"
  ) +
  theme(strip.placement = "outside", strip.background = element_blank())
save_panel(
  p_abund,
  abund_df,
  RUN_DIR,
  "fig_cluster_abundance",
  width_mm = 183,
  height_mm = 210
)

# Figure 4: ARI per donor, opt188 and the as188_joint seed noise floor ---------------------------------
ari_df <- clusters_ari |>
  mutate(
    comparison = if_else(
      comparison == "opt188",
      "opt188",
      "as188_joint, other seed"
    )
  )
ari_df <- ari_df |>
  mutate(
    comparison = factor(
      comparison,
      levels = c("opt188", "as188_joint, other seed")
    )
  )
p_ari <- ggplot(ari_df, aes(comparison, ari)) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.5,
    colour = "black"
  ) +
  geom_point(
    aes(colour = comparison, shape = comparison),
    position = position_jitter(width = 0.1, height = 0, seed = 1),
    size = 1.5
  ) +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_shape_manual(values = shape_arm, guide = "none") +
  scale_x_discrete(labels = n_labels(ari_df$comparison)) +
  labs(
    x = NULL,
    y = "ARI vs as188_joint (seed 42)",
    title = "Cluster agreement per donor"
  )
save_panel(p_ari, ari_df, RUN_DIR, "fig_ari", width_mm = 89, height_mm = 70)

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
