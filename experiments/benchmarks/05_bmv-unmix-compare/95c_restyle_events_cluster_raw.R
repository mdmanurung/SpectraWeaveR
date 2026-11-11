# Purpose: redraw the cluster raw companions from the common SpectroFlo gates (95_restyle_events_cluster_raw.R
#          repointed; plans/common_gates_plan.md G2, G8, G11). The cluster's two defining lineage markers, all
#          live cells as hexbin, the cluster's cells as dots, % of live cells per panel. Common-gate clustering
#          renumbers the shifted clusters: MC02, MC06 and MC12 are still shifted and are redrawn as the old
#          figures; MC09 is not shifted any more and has no new source; the other shifted clusters (MC01, MC04,
#          MC07, MC10, MC15, MC17, MC19) are drawn as extras. 89 mm wide, square panels.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_events_cluster_raw.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-11_bmv-common-downstream-opt188/{cluster_assignments.rds,cluster_match.csv}  [89 common]
#          experiments/runs/2026-11-11_bmv-common-cluster-raw-scatter/{defining_markers,cluster_pct}.csv   [94 common]
#          experiments/runs/2026-11-11_bmv-common-gated/gs/{as188_joint,opt188_pregate_joint_p1af1}   [97_common_gs.R]
# Outputs: experiments/runs/2026-11-11_bmv-restyle-common/21_raw_companions/cluster_raw/fig_cluster_MC<nn>.
#            {pdf,png,csv}   the .csv holds the % of live cells in the cluster per arm and donor

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(stringr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")
SCRIPT <- "95c_restyle_events_cluster_raw.R"

# Parameters ---------------------------------------------------------------------------------------
RUN_89 <- "experiments/runs/2026-11-11_bmv-common-downstream-opt188"
RUN_94 <- "experiments/runs/2026-11-11_bmv-common-cluster-raw-scatter-g17" # G17: same-sign shift and >= 1% of live cells (1.8.8 zip)
GATED <- "experiments/runs/2026-11-11_bmv-common-gated/gs" # the GatingSets 89 clustered on
OLD_CLUSTERS <- c(2, 6, 9, 12) # the certified figures
ARMS <- c(as188_joint = "as188_joint", opt188 = "opt188_pregate_joint_p1af1")
ARM_LABELS <- c(
  as188_joint = str_wrap(arm_label("as188_joint"), 16),
  opt188 = str_wrap(arm_label("opt188_pregate_joint_p1af1"), 16)
) # display names shared with the other 95c figures
COFACTOR <- 6000 # README "Transform"
STORED_ARCSINH <- c(
  "SBUV510-A",
  "Spark Blue 574-A",
  "SBUV445-A",
  "BV510-A",
  "BV570-A"
) # as 94
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
N_OVERLAY <- 1500 # cluster cells drawn as dots per panel; display only (94)
OVERLAY_SEED <- 42 # which cluster cells are drawn (94)
OVERLAY_ALPHA <- 0.2 # dot opacity, low so the hexes stay visible
AXIS <- c(-1.5, 4.5) # display window only
BINS <- 50 # hexagons per panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
WIDTH_MM <- 89 # Nature single column, two arms

defining <- read_csv(
  file.path(RUN_94, "defining_markers.csv"),
  show_col_types = FALSE
)
CLUSTERS <- sort(unique(defining$cluster)) # the shifted clusters of the common-gate 94 run
stems <- sprintf("fig_cluster_MC%02d", CLUSTERS)
OUT_DIR <- fresh_dir("21_raw_companions/cluster_raw")
# The G17 set replaces the earlier one: delete the old figures and record their names in the manifest.
old_files <- list.files(OUT_DIR, pattern = "^fig_cluster_MC", full.names = TRUE)
old_stems <- unique(tools::file_path_sans_ext(basename(old_files)))
unlink(old_files)
add_manifest(
  SCRIPT,
  "2026-11-11_bmv-restyle/21_raw_companions/cluster_raw/fig_cluster_MC09.png",
  NA_character_,
  RUN_94,
  status = "no new source",
  note = "MC09 is not shifted in all five donors under the common gates; 94 common has no defining markers for it"
)

# Live cells, lineage markers on arcsinh, rows as in 89's assignments (as 94) ------------------------------
assignments <- readRDS(file.path(RUN_89, "cluster_assignments.rds"))
match_opt <- read_csv(
  file.path(RUN_89, "cluster_match.csv"),
  show_col_types = FALSE
)
cluster_pct <- read_csv(
  file.path(RUN_94, "cluster_pct.csv"),
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
shown <- tibble(cluster = CLUSTERS) |>
  left_join(
    match_opt,
    by = join_by(cluster == ref),
    relationship = "one-to-one",
    unmatched = "drop"
  ) |>
  rename(opt_cluster = opt)

# One ggcyto figure per cluster ---------------------------------------------------------------------------
pwalk(
  list(shown$cluster, shown$opt_cluster, shown$jaccard, stems),
  \(which_cluster, opt_cluster, jaccard, stem) {
    pair <- defining |>
      dplyr::filter(cluster == which_cluster, rank <= 2) |>
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
      ) |>
      mutate(arm = factor(ARM_LABELS[arm], levels = ARM_LABELS))
    pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)
    # cluster cells per arm and donor, the same seeded subsample as 94
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
      mutate(arm = factor(ARM_LABELS[arm], levels = ARM_LABELS))
    labels <- dplyr::filter(cluster_pct, cluster == which_cluster) |>
      mutate(arm = factor(ARM_LABELS[arm], levels = ARM_LABELS))
    p <- ggcyto(fs, aes(x = !!rlang::sym(pair[1]), y = !!rlang::sym(pair[2]))) +
      geom_hex(bins = BINS) +
      geom_point(
        data = overlay,
        aes(x, y),
        inherit.aes = FALSE,
        size = 0.1,
        alpha = OVERLAY_ALPHA
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
        size = LABEL_PT,
        size.unit = "pt",
        family = FIG_FONT,
        fill = scales::alpha("white", 0.75),
        label.size = 0,
        label.padding = unit(0.1, "lines")
      ) +
      facet_grid(
        donor ~ arm,
        switch = "y",
        axes = "all",
        axis.labels = "margins"
      ) +
      ggcyto_par_set(limits = list(x = AXIS, y = AXIS)) +
      scale_fill_gradientn(
        colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
        trans = "sqrt",
        guide = "none"
      ) +
      labs(
        subtitle = "hexes = all live cells, black = cluster cells",
        caption = str_wrap(
          "hex colour = event density. Shown: clusters shifting the same way in all 5 donors and holding at least 1% of live cells (1.8.8 zip).",
          70
        ),
        x = sprintf("%s (arcsinh, cofactor 6000)", pair[1]),
        y = sprintf("%s (arcsinh, cofactor 6000)", pair[2]),
        title = sprintf(
          "MC%02d (1.8.8 zip) vs MC%02d (opt188), J = %.2f",
          which_cluster,
          opt_cluster,
          jaccard
        )
      ) +
      theme(
        strip.text.y.left = element_text(face = "bold", angle = 0),
        plot.caption = element_text(hjust = 0),
        strip.placement = "outside",
        panel.spacing = unit(1.2, "mm"),
        # both axes arcsinh, so the panels are square (r-0006)
        aspect.ratio = 1
      )
    # build on a cairo device, as save_panel() draws, so Nimbus Sans resolves through fontconfig
    cairo_pdf(tempfile(fileext = ".pdf"))
    g <- ggplotGrob(as.ggplot(p))
    dev.off()
    height_mm <- round((WIDTH_MM - 18) / length(ARMS) * length(donors) + 28)
    save_panel(g, labels, OUT_DIR, stem, WIDTH_MM, height_mm)
    add_manifest(
      SCRIPT,
      if (which_cluster %in% OLD_CLUSTERS) {
        paste0(
          "2026-11-11_bmv-restyle/21_raw_companions/cluster_raw/",
          stem,
          ".png"
        )
      } else {
        NA_character_
      },
      file.path(OUT_DIR, paste0(stem, ".png")),
      paste(RUN_89, RUN_94, sep = "; "),
      note = if (which_cluster %in% OLD_CLUSTERS) {
        ""
      } else {
        "extra: shifted under common gates, not in the old report"
      }
    )
  }
)

# Figures of the earlier set that the G17 rule dropped (MC02, MC19): the part row keeps the old file name,
# gets no new figure and the status "superseded".
part_file <- file.path(OUT_ROOT, "manifest_parts", paste0(SCRIPT, ".csv"))
part <- read_csv(part_file, col_types = cols(.default = "c"))
dropped <- !is.na(part$new_figure) & !file.exists(part$new_figure)
part$status[dropped] <- "superseded"
part$note[dropped] <- paste0(
  "earlier common-gate set replaced by the G17 selection (",
  basename(RUN_94),
  "); this cluster is not in it. Old file: ",
  part$new_figure[dropped]
)
part$new_figure[dropped] <- NA_character_
write_csv(part, part_file)
