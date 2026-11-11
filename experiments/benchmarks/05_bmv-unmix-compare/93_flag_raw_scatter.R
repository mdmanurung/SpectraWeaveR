# Purpose: raw per-cell view behind the all-marker flag figures (grill Q4, README row "Raw companions",
#          user 2026-11-11). For every marker flagged (revised Q2 rule) in as102_var, s188_c102 or opt188:
#          a ggcyto hexbin of that marker against GLUT1 (Alexa Fluor 405-A) in the population where it is
#          flagged (CD8T, else CD4T, else live), on as188_joint's gates with rows paired across arms, as the
#          flag scores were. Columns are arms (SpectroFlo and 1.8.8 zip as references), rows donors; the
#          fixed cut at arcsinh 1 drawn with geom_gate and the % above it printed per panel.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/93_flag_raw_scatter.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-marker-flags-v2/flags.csv         as102_var, s188_c102 flags [86b]
#          experiments/runs/2026-11-10_bmv-downstream-opt188/marker_flags.csv opt188 flags [89]
#          experiments/runs/2026-11-10_bmv-ablation-gated/gs/{as188_joint,opt188_pregate_joint_p1af1,as102_var}
#          experiments/runs/2026-11-10_bmv-gated-pooled/gs/{spectroflo,s188_c102}   GatingSets [19a]
# Outputs: <out_dir>/fig_flag_<marker>_<pop>.png, .pdf   one figure per flagged marker
#          <out_dir>/flagged_markers.csv                 marker, population shown, arms that flag it
#          <out_dir>/pct_above_cut.csv                   % above the cut per marker, arm and donor
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
OUT_DIR <- "experiments/runs/2026-11-11_bmv-flag-raw-scatter"
GATED <- "experiments/runs/2026-11-10_bmv-ablation-gated/gs"
GATED_POOLED <- "experiments/runs/2026-11-10_bmv-gated-pooled/gs"
# column label = GatingSet; references first, then the flagged arms
ARMS <- c(
  SpectroFlo = file.path(GATED_POOLED, "spectroflo"),
  "1.8.8 zip" = file.path(GATED, "as188_joint"),
  opt188 = file.path(GATED, "opt188_pregate_joint_p1af1"),
  "1.0.2 code,\n1.8.8 spectra" = file.path(GATED_POOLED, "s188_c102"),
  "1.0.2 zip" = file.path(GATED, "as102_var")
)
FLAG_ARM <- c(
  as102_var = "1.0.2 zip",
  s188_c102 = "1.0.2 code,\n1.8.8 spectra",
  opt188 = "opt188"
)
ARM_FILL <- c("#252525", "#d9d9d9", "#9ecae1", "#08519c", "#7f0000") # as the 72 scatters
DARK_FILLS <- c("#252525", "#08519c", "#7f0000") # strips that take white text
POP_PRIORITY <- c("CD8T", "CD4T", "live") # population shown when a marker is flagged in several
CUT <- 1 # fixed cut, arcsinh (user, 2026-11-10)
COFACTOR <- 6000 # README "Transform"
# 19a_gate_arm.R:113 writes these channels as arcsinh already; all other fluorescence channels are linear
STORED_ARCSINH <- c(
  "SBUV510-A",
  "Spark Blue 574-A",
  "SBUV445-A",
  "BV510-A",
  "BV570-A"
)
GLUT1 <- "Alexa Fluor 405-A"
AXIS <- c(-1.5, 4) # display window only; the % above the cut uses every cell
BINS <- 60 # hexagons per panel; display only

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Flagged markers and the population each is shown in ----------------------------------------------------
flags <- bind_rows(
  read_csv(
    "experiments/runs/2026-11-10_bmv-marker-flags-v2/flags.csv",
    show_col_types = FALSE
  ) |>
    dplyr::filter(source == "screen", arm %in% names(FLAG_ARM)),
  read_csv(
    "experiments/runs/2026-11-10_bmv-downstream-opt188/marker_flags.csv",
    show_col_types = FALSE
  ) |>
    mutate(arm = "opt188")
) |>
  dplyr::filter(flag)
flagged <- flags |>
  mutate(pop_rank = match(pop, POP_PRIORITY)) |>
  summarise(
    pop = POP_PRIORITY[min(pop_rank)],
    flagged_by = paste(sort(unique(arm)), collapse = "; "),
    .by = c(marker, channel)
  )
write_csv(flagged, file.path(OUT_DIR, "flagged_markers.csv"))

# Live cells of each arm and donor, rows paired on as188_joint's gates ------------------------------------
gs_list <- map(ARMS, load_gs)
gs_ref <- gs_list[["1.8.8 zip"]]
donors <- sampleNames(gs_ref)
channels <- unique(c(GLUT1, flagged$channel))
live <- map(set_names(donors), \(d) {
  which(gh_pop_get_indices(gs_ref[[d]], "live"))
})
# population membership within the live rows, for each population shown
in_pop <- map(set_names(donors), \(d) {
  map(set_names(POP_PRIORITY), \(p) {
    gh_pop_get_indices(gs_ref[[d]], p)[live[[d]]]
  })
})
cells <- imap(gs_list, \(gs, arm) {
  map(set_names(donors), \(d) {
    m <- exprs(gh_pop_get_data(gs[[d]], "root"))[
      live[[d]],
      c("FSC-A", channels),
      drop = FALSE
    ]
    linear <- setdiff(channels, STORED_ARCSINH)
    m[, linear] <- asinh(m[, linear] / COFACTOR)
    m
  })
})
# scatter does not depend on unmixing, so equal FSC-A per row shows every arm holds the same cells
stopifnot(
  "rows pair across arms" = all(map_lgl(cells, \(arm_cells) {
    all(map2_lgl(arm_cells, cells[["1.8.8 zip"]], \(a, r) {
      identical(a[, "FSC-A"], r[, "FSC-A"])
    }))
  }))
)

# % above the cut per flagged marker, arm and donor --------------------------------------------------------
pct_above <- flagged |>
  select(which_marker = marker, which_channel = channel, which_pop = pop) |>
  pmap(\(which_marker, which_channel, which_pop) {
    tidyr::expand_grid(arm = names(ARMS), donor = donors) |>
      mutate(
        marker = which_marker,
        pop = which_pop,
        pct_above_cut = map2_dbl(arm, donor, \(a, d) {
          100 *
            mean(
              cells[[a]][[d]][in_pop[[d]][[which_pop]], which_channel] >= CUT
            )
        })
      )
  }) |>
  list_rbind()
write_csv(
  mutate(pct_above, arm = gsub("\n", " ", arm)),
  file.path(OUT_DIR, "pct_above_cut.csv")
)

# One ggcyto figure per flagged marker ------------------------------------------------------------------
flag_figure <- \(which_marker, channel, pop, flagged_by) {
  frames <- tidyr::expand_grid(arm = names(ARMS), donor = donors) |>
    mutate(name = paste(arm, donor, sep = "|")) |>
    with(set_names(
      map2(arm, donor, \(a, d) {
        flowFrame(cells[[a]][[d]][in_pop[[d]][[pop]], c(GLUT1, channel)])
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
    mutate(arm = factor(arm, levels = names(ARMS)))
  pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)
  gate <- rectangleGate(
    filterId = "above cut",
    .gate = setNames(list(c(-Inf, Inf), c(CUT, Inf)), c(GLUT1, channel))
  )
  gate <- filterList(set_names(rep(list(gate), length(fs)), sampleNames(fs)))
  labels <- dplyr::filter(pct_above, marker == which_marker) |>
    mutate(arm = factor(arm, levels = names(ARMS)))
  flagged_labels <- FLAG_ARM[strsplit(flagged_by, "; ")[[1]]]
  p <- ggcyto(fs, aes(x = !!rlang::sym(GLUT1), y = !!rlang::sym(channel))) +
    geom_hex(bins = BINS) +
    geom_gate(gate, colour = "red") +
    geom_label(
      data = labels,
      aes(
        x = -Inf,
        y = Inf,
        label = sprintf("%s+ %.1f%%", which_marker, pct_above_cut)
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
      x = sprintf("GLUT1, %s (arcsinh, cofactor 6000)", GLUT1),
      y = sprintf("%s, %s (arcsinh, cofactor 6000)", which_marker, channel),
      title = sprintf(
        "%s: %s vs GLUT1 (flagged in %s)",
        pop,
        which_marker,
        paste(gsub("\n", " ", flagged_labels), collapse = ", ")
      )
    ) +
    ggpubr::theme_pubr(base_size = 8, border = TRUE) +
    theme(
      strip.text.x = element_text(size = 7),
      strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(1, "mm")
    )
  g <- ggplotGrob(as.ggplot(p))
  # colour each column strip by its arm; strip-t grobs run left to right in facet order
  strips <- grep("^strip-t", g$layout$name)
  strips <- strips[order(g$layout$l[strips])]
  for (i in seq_along(strips)) {
    rect <- g$grobs[[strips[i]]]$grobs[[1]]$children[[1]]
    rect$gp$fill <- ARM_FILL[i]
    g$grobs[[strips[i]]]$grobs[[1]]$children[[1]] <- rect
    txt <- g$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]]
    txt$gp$col <- if (ARM_FILL[i] %in% DARK_FILLS) "white" else "black"
    g$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]] <- txt
  }
  stem <- file.path(
    OUT_DIR,
    sprintf("fig_flag_%s_%s", gsub("[^A-Za-z0-9]", "", which_marker), pop)
  )
  ggsave(
    paste0(stem, ".pdf"),
    g,
    width = 160,
    height = 175,
    units = "mm",
    device = cairo_pdf
  )
  ggsave(
    paste0(stem, ".png"),
    g,
    width = 160,
    height = 175,
    units = "mm",
    dpi = 200,
    bg = "white"
  )
}
flagged |>
  select(which_marker = marker, channel, pop, flagged_by) |>
  pwalk(flag_figure)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
