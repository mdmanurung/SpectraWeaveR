# Purpose: raw per-cell view behind the opt188 vs as188_joint cluster-abundance figure (grill Q4, README
#          row "Raw companions", user 2026-11-11). For each metacluster that shifts the same way in all five
#          donors (MC02, MC06, MC09, MC12 of as188_joint; 89's cluster_abundance.csv), the two lineage
#          markers whose cluster median differs most from the live-cell median (as188_joint) are plotted
#          against each other: ggcyto hexbin of all live cells, the cluster's cells overlaid as black dots
#          (as188_joint's cluster in its column, the matched opt188 cluster in the opt188 column).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/94_cluster_raw_scatter.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-downstream-opt188/cluster_assignments.rds  live rows, labels [89]
#          experiments/runs/2026-11-10_bmv-downstream-opt188/cluster_match.csv        matched clusters [89]
#          experiments/runs/2026-11-10_bmv-ablation-gated/gs/{as188_joint,opt188_pregate_joint_p1af1}  [19a]
# Outputs: <out_dir>/fig_cluster_MC<nn>.png, .pdf   one figure per shifted metacluster
#          <out_dir>/defining_markers.csv           cluster, marker, cluster median minus live median
#          <out_dir>/cluster_pct.csv                % of live cells in the cluster, per arm and donor
#          <out_dir>/sessionInfo.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
})

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- "experiments/runs/2026-11-11_bmv-cluster-raw-scatter"
RUN_89 <- "experiments/runs/2026-11-10_bmv-downstream-opt188"
GATED <- "experiments/runs/2026-11-10_bmv-ablation-gated/gs" # the GatingSets 89 clustered on
ARMS <- c(as188_joint = "as188_joint", opt188 = "opt188_pregate_joint_p1af1")
CLUSTERS <- c(2, 6, 9, 12) # shift the same way in all 5 donors (89 cluster_abundance.csv; directive)
COFACTOR <- 6000 # README "Transform"
# 19a_gate_arm.R:113 writes these channels as arcsinh already; all other fluorescence channels are linear
STORED_ARCSINH <- c(
  "SBUV510-A",
  "Spark Blue 574-A",
  "SBUV445-A",
  "BV510-A",
  "BV570-A"
)
# README "Cluster markers", channels as in 89
LINEAGE <- c(
  CD3 = "SBUV510-A",
  CD4 = "Spark Blue 574-A",
  CD8 = "SBUV445-A",
  CD45RA = "Spark UV 387-A",
  CCR7 = "PE-Fire 810-A",
  CD127 = "APC-Fire 810-A",
  CD25 = "SBUV605-A",
  CD56 = "RB744-A",
  TCRgd = "BUV395-A",
  CD161 = "BUV563-A",
  CD57 = "BV785-A",
  CD38 = "SBUV740-A",
  CD7 = "BUV805-A",
  CXCR5 = "BV750-A"
)
N_OVERLAY <- 1500 # cluster cells drawn as dots per panel; display only
OVERLAY_SEED <- 42 # which cluster cells are drawn
AXIS <- c(-1.5, 4.5) # display window only
BINS <- 60 # hexagons per panel; display only

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Live cells, lineage markers on arcsinh, rows as in 89's assignments --------------------------------------
assignments <- readRDS(file.path(RUN_89, "cluster_assignments.rds"))
match_opt <- read_csv(
  file.path(RUN_89, "cluster_match.csv"),
  show_col_types = FALSE
)
gs_list <- map(ARMS, \(a) load_gs(file.path(GATED, a)))
donors <- unique(assignments$donor)
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
write_csv(cluster_pct, file.path(OUT_DIR, "cluster_pct.csv"))

# One ggcyto figure per cluster ---------------------------------------------------------------------------
cluster_figure <- \(which_cluster, opt_cluster, jaccard) {
  pair <- defining |>
    dplyr::filter(cluster == as.character(which_cluster), rank <= 2) |>
    pull(marker)
  label_of <- c(as188_joint = which_cluster, opt188 = opt_cluster)
  frames <- tidyr::expand_grid(arm = names(ARMS), donor = donors) |>
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
    slice_sample(n = N_OVERLAY, by = c(arm, donor))
  labels <- dplyr::filter(cluster_pct, cluster == which_cluster)
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
      size = 7,
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
        "MC%02d (1.8.8 zip) vs MC%02d (opt188), J = %.2f",
        which_cluster,
        opt_cluster,
        jaccard
      )
    ) +
    ggpubr::theme_pubr(base_size = 8, border = TRUE) +
    theme(
      strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(1, "mm")
    )
  stem <- file.path(OUT_DIR, sprintf("fig_cluster_MC%02d", which_cluster))
  ggsave(
    paste0(stem, ".pdf"),
    as.ggplot(p),
    width = 100,
    height = 175,
    units = "mm",
    device = cairo_pdf
  )
  ggsave(
    paste0(stem, ".png"),
    as.ggplot(p),
    width = 100,
    height = 175,
    units = "mm",
    dpi = 200,
    bg = "white"
  )
}
shown |>
  select(which_cluster = cluster, opt_cluster, jaccard) |>
  pwalk(cluster_figure)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
