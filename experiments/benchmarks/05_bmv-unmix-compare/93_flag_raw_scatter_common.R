# Purpose: copy of 93_flag_raw_scatter.R on the common SpectroFlo gates (plans/common_gates_plan.md, G2, G3, G7a).
#          Raw per-cell view behind the all-marker flag figures (grill Q4, README row "Raw companions", user
#          2026-11-11). For every marker flagged (revised Q2 rule) in as102_var, s188_c102 or opt188: a ggcyto
#          hexbin of that marker against GLUT1 (Alexa Fluor 405-A) in the population where it is flagged (CD8T,
#          else CD4T, else live), on the common SpectroFlo labels, as the flag scores were. Columns are arms
#          (SpectroFlo and 1.8.8 zip as references), rows donors; the fixed cut at arcsinh 1.5 (G7a) is drawn
#          with geom_gate and the % above it printed per panel. The script stops unless FSC-A and SSC-A of
#          every root row are identical in every arm (G8).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/93_flag_raw_scatter_common.R <flags.csv>
#            (repo root); <flags.csv> is the new 86b flags.csv from batch 3
# Inputs:  <flags.csv>                                                        as102_var, s188_c102 flags [86b, batch 3]
#          experiments/runs/2026-11-11_bmv-common-downstream-opt188/marker_flags.csv  opt188 flags [89 common]
#          experiments/runs/2026-11-11_bmv-common-gated/gs/{as188_joint,opt188_pregate_joint_p1af1,as102_var,
#            spectroflo,s188_c102}   common-gate GatingSets [batch 0], README.md
# Outputs: <out_dir>/fig_flag_<marker>_<pop>.{png,pdf,csv}   one figure per flagged marker, 183 mm wide, 300 dpi PNG
#          <out_dir>/flagged_markers.csv                 marker, population shown, arms that flag it
#          <out_dir>/pct_above_cut.csv                   % above the cut per marker, arm and donor
#          <out_dir>/pairing.csv                         arm, donor, root rows, FSC-A and SSC-A equal to 1.8.8 zip
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
OUT_DIR <- "experiments/runs/2026-11-11_bmv-common-flag-raw-scatter"
FLAGS_86B <- commandArgs(trailingOnly = TRUE)[1] # new 86b flags.csv, batch 3
FLAGS_89 <- "experiments/runs/2026-11-11_bmv-common-downstream-opt188/marker_flags.csv"
COMMON <- "experiments/runs/2026-11-11_bmv-common-gated" # batch 0; gs/<arm> for all arms
FIG_WIDTH_MM <- 183 # double column, the width save_panel takes besides 89
PAIR_COLS <- c("FSC-A", "SSC-A") # equal per row means every arm holds the same cells (G8)
# column label = GatingSet; references first, then the flagged arms
ARMS <- file.path(
  COMMON,
  "gs",
  c(
    SpectroFlo = "spectroflo",
    "1.8.8 zip" = "as188_joint",
    opt188 = "opt188_pregate_joint_p1af1",
    "1.0.2 code,\n1.8.8 spectra" = "s188_c102",
    "1.0.2 zip" = "as102_var"
  )
) |>
  set_names(c(
    "SpectroFlo",
    "1.8.8 zip",
    "opt188",
    "1.0.2 code,\n1.8.8 spectra",
    "1.0.2 zip"
  ))
FLAG_ARM <- c(
  as102_var = "1.0.2 zip",
  s188_c102 = "1.0.2 code,\n1.8.8 spectra",
  opt188 = "opt188"
)
ARM_FILL <- c("#252525", "#d9d9d9", "#9ecae1", "#08519c", "#7f0000") # as the 72 scatters
DARK_FILLS <- c("#252525", "#08519c", "#7f0000") # strips that take white text
POP_PRIORITY <- c("CD8T", "CD4T", "live") # population shown when a marker is flagged in several
CUT <- 1.5 # fixed cut, arcsinh, cofactor 6000 (G7a, user 2026-11-11; was 1)
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

stopifnot(
  "batch 0 finished (README.md in the common dir)" = file.exists(file.path(
    COMMON,
    "README.md"
  )),
  "pass the new 86b flags.csv as the first argument" = !is.na(FLAGS_86B) &&
    file.exists(FLAGS_86B),
  "89 common run finished" = file.exists(FLAGS_89),
  "output dir must not exist" = !dir.exists(OUT_DIR)
)
dir.create(OUT_DIR)
source(file.path(
  Sys.getenv("SW_LIB", "experiments/benchmarks/05_bmv-unmix-compare/lib"),
  "theme.R"
)) # save_panel, FIG_FONT

# Flagged markers and the population each is shown in ----------------------------------------------------
flags <- bind_rows(
  read_csv(FLAGS_86B, show_col_types = FALSE) |>
    dplyr::filter(source == "screen", arm %in% names(FLAG_ARM)),
  read_csv(FLAGS_89, show_col_types = FALSE) |>
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

# Live cells of each arm and donor, rows paired on the common labels ------------------------------------
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
# scatter does not depend on unmixing, so equal FSC-A and SSC-A on every root row show every arm holds the same cells
root_scatter <- \(gs, d) exprs(gh_pop_get_data(gs[[d]], "root"))[, PAIR_COLS]
pairing <- gs_list |>
  imap(\(gs, arm) {
    tibble(
      donor = donors,
      n_rows = map_int(donors, \(d) nrow(root_scatter(gs, d))),
      paired = map_lgl(donors, \(d) {
        identical(root_scatter(gs, d), root_scatter(gs_ref, d))
      })
    )
  }) |>
  list_rbind(names_to = "arm")
write_csv(
  mutate(pairing, arm = gsub("\n", " ", arm)),
  file.path(OUT_DIR, "pairing.csv")
)
stopifnot("rows pair across arms" = all(pairing$paired))

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
    ggpubr::theme_pubr(base_size = 7, base_family = FIG_FONT, border = TRUE) +
    theme(
      plot.title = element_text(size = 7, face = "bold"),
      strip.text.x = element_text(size = 7),
      strip.text.y.left = element_text(size = 7, face = "bold", angle = 0),
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
  save_panel(
    g,
    mutate(labels, arm = gsub("\n", " ", arm)),
    OUT_DIR,
    sprintf("fig_flag_%s_%s", gsub("[^A-Za-z0-9]", "", which_marker), pop),
    width_mm = FIG_WIDTH_MM,
    height_mm = 175
  )
}
flagged |>
  select(which_marker = marker, channel, pop, flagged_by) |>
  pwalk(flag_figure)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
