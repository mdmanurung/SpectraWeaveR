# Purpose: G17 copy of 94 (plans/common_gates_plan.md, G17): a cluster is shown only if it shifts the same way in
#          all donors AND its donor-median share of live cells in the as188_joint reference is at least
#          MIN_REF_PCT. The rule is stated in each figure caption and in cluster_pct.csv / the figure CSVs.
#          Original header of 94 follows. Copy of 94_cluster_raw_scatter.R on the common SpectroFlo gates (plans/common_gates_plan.md, G2, G3).
#          Raw per-cell view behind the opt188 vs as188_joint cluster-abundance figure (grill Q4, README row
#          "Raw companions", user 2026-11-11). Clusters shown are the as188_joint metaclusters whose matched
#          opt188 abundance shifts the same way in all donors (all five in 3D; within each condition when there are conditions; the rule that picked MC02, MC06, MC09 and
#          MC12 in the old 89; metacluster numbers belong to one clustering, so they are re-derived from the
#          new 89 run). For each, the two lineage markers whose cluster median differs most from the live-cell
#          median (as188_joint) are plotted against each other: ggcyto hexbin of all live cells, the cluster's
#          cells overlaid as black dots (as188_joint's cluster in its column, the matched opt188 cluster in the
#          opt188 column). The script stops unless FSC-A and SSC-A of every root row pair across the two arms.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> \
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/94b_cluster_raw_scatter_common_g17.R   (repo root)
#          The out dir is run_out("common-cluster-raw-scatter-g17").
# Inputs:  datasets/<SW_DATASET>.json, <SW_DATASET>_runs.csv, samples_tbl(); run folders come from these registry keys
#          run_in("common-downstream-opt188")/{cluster_assignments.rds,cluster_match.csv,
#            cluster_abundance.csv}   live rows, labels, abundance [89_downstream_as188_vs_opt188_common.R]
#          run_in("common-gated")/gs/{as188_joint,opt188_pregate_joint_p1af1}  [batch 0]
#          cfg("channels"), cfg("cluster_markers_other")   cluster markers and channels, built as in 89
# Outputs: with conditions, shifted_clusters.csv and cluster_pct.csv gain a condition column, the same-way
#          rule is applied within each condition, and each figure is one cluster in one condition (suffix _<condition>)
#          <out_dir>/cluster_selection.csv   per as188_joint cluster: same_way, median ref %, median shift, selected
#          <out_dir>/fig_cluster_MC<nn>.{png,pdf,csv}   one figure per shifted metacluster, 89 mm wide, 300 dpi PNG
#          <out_dir>/shifted_clusters.csv          as188_joint cluster, donor, opt188 minus as188_joint % of live
#          <out_dir>/defining_markers.csv          cluster, marker, cluster median minus live median
#          <out_dir>/cluster_pct.csv               % of live cells in the cluster, per arm and donor
#          <out_dir>/sessionInfo.txt

# dataset.R loads only dplyr and readr, so it can run before .libPaths puts ggplot2 3.5.2 (for ggcyto) first
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
LIB_GGPLOT2 <- grep("lib_ggplot2", cfg("lib_dirs"), value = TRUE, fixed = TRUE)
stopifnot("one ggplot2 lib dir in lib_dirs" = length(LIB_GGPLOT2) == 1)
.libPaths(c(LIB_GGPLOT2, .libPaths()))
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
})

# dataset.R calls a bare filter(); here flowCore's filter generic comes first on the search path, so give
# its two helpers a scope where filter is dplyr's (dataset.R should say dplyr::filter itself)
.dataset_scope <- new.env(parent = globalenv())
.dataset_scope$filter <- dplyr::filter
environment(run_in) <- .dataset_scope
environment(cut_of) <- .dataset_scope

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- run_out("common-cluster-raw-scatter-g17")
RUN_89 <- run_in("common-downstream-opt188")
COMMON <- run_in("common-gated") # the GatingSets 89 clustered on
LIB <- Sys.getenv("SW_LIB", "experiments/benchmarks/05_bmv-unmix-compare/lib")
ARMS <- c(as188_joint = "as188_joint", opt188 = "opt188_pregate_joint_p1af1")
PAIR_COLS <- c("FSC-A", "SSC-A") # equal per row means the arms hold the same cells (G8)
COFACTOR <- 6000 # README "Transform"
# 19a_gate_arm.R:113 writes these channels as arcsinh already; all other fluorescence channels are linear
STORED_ARCSINH <- unlist(
  cfg("channels")[c("CD3", "CD4", "CD8", "TNF_slot", "IFNg_slot")],
  use.names = FALSE
)
# marker -> channel in the order 89 clustered on (CLUSTER_MARKERS of 89; README "Cluster markers")
LINEAGE <- c(
  CD3 = cfg("channels")$CD3,
  CD4 = cfg("channels")$CD4,
  CD8 = cfg("channels")$CD8,
  unlist(cfg("cluster_markers_other")) # the other lineage markers, channel by marker, as in 89
)
MIN_REF_PCT <- 1 # minimum donor-median % of live cells in the as188_joint cluster (G17)
HEIGHT_PER_ROW_MM <- 30 # figure height per sample row; 5 rows give the 150 mm of the 3D figure
N_OVERLAY <- 1500 # cluster cells drawn as dots per panel; display only
OVERLAY_SEED <- 42 # which cluster cells are drawn
AXIS <- c(-1.5, 4.5) # display window only
BINS <- 60 # hexagons per panel; display only

stopifnot(
  "89 run finished" = file.exists(file.path(RUN_89, "sessionInfo.txt"))
)
dir.create(OUT_DIR)
source(file.path(LIB, "theme.R")) # save_panel, FIG_FONT

# Live cells, lineage markers on arcsinh, rows as in 89's assignments --------------------------------------
assignments <- readRDS(file.path(RUN_89, "cluster_assignments.rds"))
match_opt <- read_csv(
  file.path(RUN_89, "cluster_match.csv"),
  show_col_types = FALSE
)
gs_list <- map(ARMS, \(a) load_gs(file.path(COMMON, "gs", a)))
donors <- unique(assignments$donor)
smp <- samples_tbl()
stopifnot(
  "assignment samples are the sample table" = setequal(donors, smp$sample)
)
# with conditions, a condition column joins the sample-level tables and the rules run within each condition
COND_KEY <- if (is.null(DS$conditions)) character() else "condition"
add_condition <- \(d) {
  if (is.null(DS$conditions)) {
    return(d)
  }
  left_join(
    d,
    select(smp, donor = sample, condition),
    by = join_by(donor),
    relationship = "many-to-one",
    unmatched = "error"
  )
}
# rows per figure: all samples without conditions, else the samples of each condition
sample_sets <- if (is.null(DS$conditions)) {
  list(all = donors)
} else {
  split(
    smp$sample[smp$sample %in% donors],
    smp$condition[smp$sample %in% donors]
  )
}
stopifnot(
  "FSC-A and SSC-A of every root row pair across the two arms" = all(map_lgl(
    donors,
    \(d) {
      identical(
        exprs(gh_pop_get_data(gs_list$as188_joint[[d]], "root"))[, PAIR_COLS],
        exprs(gh_pop_get_data(gs_list$opt188[[d]], "root"))[, PAIR_COLS]
      )
    }
  ))
)

# Clusters whose matched abundance shifts the same way in all donors (opt188 minus as188_joint) -----------
shift <- read_csv(
  file.path(RUN_89, "cluster_abundance.csv"),
  show_col_types = FALSE
) |>
  select(donor, ref, arm, pct) |>
  tidyr::pivot_wider(names_from = arm, values_from = pct) |>
  mutate(diff_pct = opt188 - as188_joint) |>
  add_condition()
write_csv(shift, file.path(OUT_DIR, "shifted_clusters.csv"))
selection <- shift |>
  summarise(
    same_way = all(diff_pct > 0) | all(diff_pct < 0),
    median_ref_pct = median(as188_joint),
    median_shift_pp = median(diff_pct),
    .by = all_of(c("ref", COND_KEY))
  ) |>
  mutate(
    selected = same_way & median_ref_pct >= MIN_REF_PCT,
    set = if (length(COND_KEY)) .data$condition else "all"
  ) |>
  arrange(ref, set)
write_csv(select(selection, -set), file.path(OUT_DIR, "cluster_selection.csv"))
selected <- dplyr::filter(selection, selected)
CLUSTERS <- sort(unique(selected$ref))
# the donors counted in the rule are the samples of the figure's set
rule_of <- \(n_donors) {
  sprintf(
    "Shown: shifts the same way in all %d donors and holds >= %g%% of live cells in 1.8.8 (donor median)",
    n_donors,
    MIN_REF_PCT
  )
}
stopifnot(
  "at least one cluster shifts the same way in all donors" = length(CLUSTERS) >
    0
)
lineage <- map(gs_list, \(gs) {
  map(donors, \(d) {
    rows <- assignments$event[assignments$donor == d]
    m <- exprs(gh_pop_get_data(gs[[d]], "root"))[
      rows,
      unname(LINEAGE),
      drop = FALSE
    ]
    linear <- setdiff(colnames(m), STORED_ARCSINH)
    m[, linear] <- asinh(m[, linear] / COFACTOR)
    colnames(m) <- names(LINEAGE)
    m
  }) |>
    reduce(rbind)
})
stopifnot(
  "matrix rows match assignments" = nrow(lineage$as188_joint) ==
    nrow(assignments)
)

# Two defining markers per cluster: largest |cluster median - live median| in as188_joint ------------------
live_median <- apply(lineage$as188_joint, 2, median)
defining <- map(set_names(CLUSTERS), \(k) {
  in_k <- assignments$as188_joint == k
  diff <- apply(lineage$as188_joint[in_k, , drop = FALSE], 2, median) -
    live_median
  tibble(marker = names(diff), median_diff = unname(diff)) |>
    arrange(desc(abs(median_diff))) |>
    mutate(rank = row_number())
}) |>
  list_rbind(names_to = "cluster")
write_csv(defining, file.path(OUT_DIR, "defining_markers.csv"))

# % of live cells in each shown cluster, per arm and donor -------------------------------------------------
shown <- tibble(cluster = CLUSTERS) |>
  left_join(
    match_opt,
    by = join_by(cluster == ref),
    relationship = "one-to-one",
    unmatched = "drop"
  ) |>
  rename(opt_cluster = opt)
cluster_pct <- shown |>
  pmap(\(cluster, opt_cluster, jaccard) {
    assignments |>
      summarise(
        as188_joint = 100 * mean(as188_joint == cluster),
        opt188 = 100 * mean(opt188 == opt_cluster),
        .by = donor
      ) |>
      tidyr::pivot_longer(
        c(as188_joint, opt188),
        names_to = "arm",
        values_to = "pct_live"
      ) |>
      mutate(cluster = cluster, opt_cluster = opt_cluster)
  }) |>
  list_rbind()
cluster_pct <- add_condition(cluster_pct)
write_csv(cluster_pct, file.path(OUT_DIR, "cluster_pct.csv"))

# One ggcyto figure per cluster ---------------------------------------------------------------------------
cluster_figure <- \(which_cluster, opt_cluster, jaccard, set) {
  set_donors <- sample_sets[[set]]
  suffix <- if (is.null(DS$conditions)) "" else paste0("_", set)
  pair <- defining |>
    dplyr::filter(cluster == as.character(which_cluster), rank <= 2) |>
    pull(marker)
  label_of <- c(as188_joint = which_cluster, opt188 = opt_cluster)
  frames <- tidyr::expand_grid(arm = names(ARMS), donor = set_donors) |>
    mutate(name = paste(arm, donor, sep = "|")) |>
    with(set_names(
      map2(arm, donor, \(a, d) {
        flowFrame(lineage[[a]][assignments$donor == d, pair, drop = FALSE])
      }),
      name
    ))
  fs <- flowSet(frames)
  pd <- tibble(name = sampleNames(fs)) |>
    tidyr::separate_wider_delim(
      name,
      "|",
      names = c("arm", "donor"),
      cols_remove = FALSE
    )
  pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)
  # cluster cells per arm and donor, a seeded subsample drawn as dots
  set.seed(OVERLAY_SEED)
  overlay <- imap(label_of, \(k, a) {
    in_k <- assignments[[a]] == k
    tibble(
      arm = a,
      donor = assignments$donor[in_k],
      x = lineage[[a]][in_k, pair[1]],
      y = lineage[[a]][in_k, pair[2]]
    )
  }) |>
    list_rbind() |>
    slice_sample(n = N_OVERLAY, by = c(arm, donor)) |>
    dplyr::filter(donor %in% set_donors)
  rule <- rule_of(length(set_donors))
  labels <- dplyr::filter(
    cluster_pct,
    cluster == which_cluster,
    donor %in% set_donors
  ) |>
    mutate(rule = rule)
  sel <- dplyr::filter(selection, ref == which_cluster, set == .env$set)
  p <- ggcyto(fs, aes(x = !!rlang::sym(pair[1]), y = !!rlang::sym(pair[2]))) +
    geom_hex(bins = BINS) +
    geom_point(
      data = overlay,
      aes(x, y),
      inherit.aes = FALSE,
      size = 0.15,
      alpha = 0.4
    ) +
    geom_label(
      data = labels,
      aes(
        x = -Inf,
        y = Inf,
        label = sprintf("cluster %.1f%% of live", pct_live)
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 6,
      size.unit = "pt",
      fontface = "bold",
      fill = "white",
      alpha = 0.85,
      label.size = 0
    ) +
    facet_grid(donor ~ arm, switch = "y") +
    ggcyto_par_set(limits = list(x = AXIS, y = AXIS)) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      trans = "sqrt",
      guide = "none"
    ) +
    labs(
      x = sprintf("%s (arcsinh, cofactor 6000)", pair[1]),
      y = sprintf("%s (arcsinh, cofactor 6000)", pair[2]),
      title = sprintf(
        "MC%02d (1.8.8) vs MC%02d (opt188), J = %.2f",
        which_cluster,
        opt_cluster,
        jaccard
      ),
      caption = sprintf(
        "%s.\nDonor median: %.1f%% of live in 1.8.8, opt188 minus 1.8.8 = %+.1f pp.",
        stringr::str_wrap(rule, 70),
        sel$median_ref_pct,
        sel$median_shift_pp
      )
    ) +
    ggpubr::theme_pubr(base_size = 7, base_family = FIG_FONT, border = TRUE) +
    theme(
      plot.title = element_text(size = 7, face = "bold"),
      plot.caption = element_text(size = 5, hjust = 0),
      strip.text.x = element_text(size = 7),
      strip.text.y.left = element_text(size = 7, face = "bold", angle = 0),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(1, "mm")
    )
  save_panel(
    as.ggplot(p),
    labels,
    OUT_DIR,
    sprintf("fig_cluster_MC%02d%s", which_cluster, suffix),
    width_mm = 89,
    height_mm = HEIGHT_PER_ROW_MM * length(set_donors)
  )
}
# one figure per selected cluster and set (one set without conditions)
shown |>
  inner_join(
    select(selected, cluster = ref, set),
    by = join_by(cluster),
    relationship = "one-to-many",
    unmatched = "drop"
  ) |>
  select(which_cluster = cluster, opt_cluster, jaccard, set) |>
  pwalk(cluster_figure)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
