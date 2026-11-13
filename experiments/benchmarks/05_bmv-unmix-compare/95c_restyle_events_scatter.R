# Purpose: Nature-spec redraw on the common gates (plans/common_gates_plan.md G7a, G11) of the CD8T cytokine
#          hexbin grids of 95_restyle_events_scatter.R: same layout, but GatingSets come from the common-gate
#          dir, so every arm holds the same SpectroFlo-labelled CD8T cells, and the cut is the dataset's positive
#          cut (cut_of(); arcsinh 1.5 in 3D). The only line drawn is the cut (gate lines belong to SpectroFlo
#          panels of the gate figures, G11). Arm strips are filled with the arm colour key; one 183 mm figure per
#          x marker. The th2 figure stops on a dataset without Th2 (EXV).
# Run:     SW_DATASET=3d SW_RUN_DATE=YYYY-MM-DD pixi run Rscript
#            experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_events_scatter.R <figure>
#            from the repo root; <figure> is one of names(FIGURES): per_dye_v8, af_swap_v3, per_dye,
#            tnf_new, th2.
# Inputs:  run_in("common-gated")/gs/<arm>/   common-gate GatingSets [97e_common_gs.R]
#          experiments/benchmarks/05_bmv-unmix-compare/95c_common.R   theme, donor key, add_manifest()
#          datasets/<SW_DATASET>.json   lib_dirs (ggplot2 3.5 entry), channels, targets
# Outputs: run_out("restyle-common") (shared by the 95c scripts)/<folder>/<stem>[_<x>].pdf, .png, .csv
#            the .csv holds the % above the cut printed in each panel

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
# ggplot2 3.5.2 for ggcyto: only the ggplot2_35 entry of lib_dirs, as before
.libPaths(c(
  grep("lib_ggplot2_35$", cfg("lib_dirs"), value = TRUE),
  .libPaths()
))
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")
SCRIPT <- "95c_restyle_events_scatter.R"

# Parameters ---------------------------------------------------------------------------------------
FIGURE <- commandArgs(trailingOnly = TRUE)[1]
COMMON_GS <- file.path(run_in("common-gated"), "gs")
REFERENCE <- c("SpectroFlo\n(reference)" = "spectroflo")
REFERENCE_FILL <- "#252525"
# arms (label = arm folder), strip fill per arm, y marker, x markers, summary folder and file stem;
# arms and fills as in 72 (per_dye, tnf_new, th2) and its earlier versions (per_dye_v3, af_swap_v1)
FIGURES <- list(
  per_dye_v3 = list(
    arms = c(
      "1.8.8\nvariants" = "as188_joint",
      "+SBUV510\n1.0.2" = "var102_only_SBUV510",
      "+ AF405\n1.0.2" = "var102_only_Alexa_Fluor_405",
      "all 1.0.2\nvariants" = "swap_var102"
    ),
    fill = c("#d9d9d9", "#fdbb84", "#e34a33", "#7f0000"),
    y = "TNF",
    blocks = c("IFNg", "GLUT1", "SSC"),
    folder = "05_tnf_cause/per_dye",
    stem = "fig_per_dye_ggcyto",
    reference = FALSE
  ),
  af_swap_v1 = list(
    arms = c(
      "1.8.8\nall" = "as188_joint",
      "+AF\n1.0.2" = "swap_af102",
      "+variants\n1.0.2" = "swap_var102",
      "+AF+var\n1.0.2" = "swap_af102_var102",
      "1.0.2, AF\n1.8.8" = "rev_af188",
      "1.0.2\nall" = "s102_c188"
    ),
    fill = c("#d9d9d9", "#fdbb84", "#e34a33", "#7f0000", "#9ecae1", "#08519c"),
    y = "TNF",
    blocks = c("IFNg", "AF", "SSC"),
    folder = "05_tnf_cause/af_swap",
    stem = "fig_af_swap_ggcyto",
    reference = FALSE
  ),
  # latest 72 sets for section 05 (freshness audit, 2026-11-11): per-dye-ggcyto-v8 and af-swap-ggcyto-v3
  per_dye_v8 = list(
    arms = c(
      "1.8.8\nvariants" = "as188_joint",
      "opt188\n(optimized)" = "opt188_pregate_joint_p1af1",
      "+SBUV510\n1.0.2" = "var102_only_SBUV510",
      "+ AF405\n1.0.2" = "var102_only_Alexa_Fluor_405",
      "+SBUV510\n+AF405 1.0.2" = "var102_only_SBUV510__Alexa_Fluor_405",
      "all 1.0.2\nvariants" = "swap_var102"
    ),
    fill = c("#d9d9d9", "#9ecae1", "#fdbb84", "#e34a33", "#b30000", "#7f0000"),
    y = "TNF",
    blocks = c("IFNg", "GLUT1", "SSC"),
    folder = "05_tnf_cause/per_dye_v8",
    stem = "fig_per_dye_ggcyto",
    reference = TRUE
  ),
  af_swap_v3 = list(
    arms = c(
      "1.8.8\nall" = "as188_joint",
      "+AF\n1.0.2" = "swap_af102",
      "+variants\n1.0.2" = "swap_var102",
      "+AF+var\n1.0.2" = "swap_af102_var102",
      "1.0.2, AF\n1.8.8" = "rev_af188",
      "1.0.2\nall" = "s102_c188"
    ),
    fill = c("#d9d9d9", "#fdbb84", "#e34a33", "#7f0000", "#9ecae1", "#08519c"),
    y = "TNF",
    blocks = c("IFNg", "AF", "SSC"),
    folder = "05_tnf_cause/af_swap_v3",
    stem = "fig_af_swap_ggcyto",
    reference = TRUE
  ),
  per_dye = list(
    arms = c(
      "1.8.8\nvariants" = "as188_joint",
      "opt188\n(optimized)" = "opt188_pregate_joint_p1af1",
      "+SBUV510\n1.0.2" = "var102_only_SBUV510",
      "+ AF405\n1.0.2" = "var102_only_Alexa_Fluor_405",
      "+SBUV510\n+AF405 1.0.2" = "var102_only_SBUV510__Alexa_Fluor_405",
      "all 1.0.2\nvariants" = "swap_var102"
    ),
    fill = c("#d9d9d9", "#9ecae1", "#fdbb84", "#e34a33", "#b30000", "#7f0000"),
    y = "TNF",
    blocks = "GLUT1",
    folder = "16_scatter_spectroflo",
    stem = "per_dye_TNF_vs_GLUT1",
    reference = TRUE
  ),
  tnf_new = list(
    arms = c(
      "1.8.8\nzip" = "as188_joint",
      "opt188" = "opt188_pregate_joint_p1af1",
      "refined\nAF" = "af405_afrefine",
      "AF405 1.0.2\ncode, clean" = "var102code_clean_af405",
      "1.0.2 code\n1.8.8 spectra" = "s188_c102",
      "1.0.2\nzip" = "as102_var"
    ),
    fill = c("#d9d9d9", "#9ecae1", "#08519c", "#fdbb84", "#e34a33", "#7f0000"),
    y = "TNF",
    blocks = "GLUT1",
    folder = "16_scatter_spectroflo",
    stem = "tnf_new_TNF_vs_GLUT1",
    reference = TRUE
  ),
  th2 = list(
    arms = c(
      "1.8.8\nzip" = "as188_joint",
      "opt188" = "opt188_pregate_joint_p1af1",
      "+ AF405\n1.0.2" = "var102_only_Alexa_Fluor_405",
      "AF405 KO\n(no BV510)" = "af405proof_ko_bv510",
      "1.0.2 code\n1.8.8 spectra" = "s188_c102",
      "1.0.2\nzip" = "as102_var"
    ),
    fill = c("#d9d9d9", "#9ecae1", "#e34a33", "#fdbb84", "#08519c", "#7f0000"),
    y = "Th2",
    blocks = "GLUT1",
    folder = "16_scatter_spectroflo",
    stem = "th2_Th2_vs_GLUT1",
    reference = TRUE
  )
)
DARK_FILLS <- c("#e34a33", "#b30000", "#7f0000", "#08519c", REFERENCE_FILL) # strips that take white text
POP <- "CD8T"
COFACTOR <- 6000 # arcsinh cofactor for linear channels (CLAUDE.md)
AXIS <- c(-1.5, 4) # display window only; the % above the cut uses every event
BINS <- 50 # hexagon bins per panel; display only, for ~22 mm panels
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
# role (key used in FIGURES) -> channel; Th2 is added below when the figure needs it
CHANNELS <- c(
  TNF = cfg("channels")$TNF_slot,
  IFNg = cfg("channels")$IFNg_slot,
  GLUT1 = paste0(cfg("targets")$glut1_dye, "-A"),
  AF = "AF-A",
  SSC = "SSC-A"
)
# role -> label shown on axes and in titles (the dataset's marker names for the BV510 and BV570 slots)
LABELS <- c(
  TNF = cfg("targets")$tnf_label,
  IFNg = cfg("targets")$ifng_label,
  GLUT1 = "GLUT1",
  AF = "AF",
  SSC = "SSC"
)
LINEAR_TO_ARCSINH <- c(paste0(cfg("targets")$glut1_dye, "-A"), "AF-A")
SCALE_LABEL <- c(
  TNF = "arcsinh",
  Th2 = "arcsinh, cofactor 6000",
  IFNg = "arcsinh, cofactor 6000",
  GLUT1 = "arcsinh, cofactor 6000",
  AF = "arcsinh, cofactor 6000",
  SSC = "linear, x 10^6"
)
X_LIMITS <- list(IFNg = AXIS, GLUT1 = AXIS, AF = AXIS, SSC = c(0, 2.5e6))

spec <- FIGURES[[FIGURE]]
stopifnot("FIGURE is one of names(FIGURES)" = !is.null(spec))
if ("Th2" %in% c(spec$y, spec$blocks)) {
  CHANNELS[["Th2"]] <- need_target(SCRIPT, "Th2", "channels")
  LABELS[["Th2"]] <- need_target(SCRIPT, "th2_label")
  LINEAR_TO_ARCSINH <- c(LINEAR_TO_ARCSINH, CHANNELS[["Th2"]])
}
# only the channels this figure plots are cut from each frame
CH_USED <- CHANNELS[unique(c(spec$y, spec$blocks))]
LINEAR_TO_ARCSINH <- intersect(LINEAR_TO_ARCSINH, CH_USED)
ARMS <- if (spec$reference) c(REFERENCE, spec$arms) else spec$arms
ARM_FILL <- set_names(
  c(if (spec$reference) REFERENCE_FILL, spec$fill),
  names(ARMS)
)
Y_NAME <- LABELS[[spec$y]]
Y_CH <- CHANNELS[[spec$y]]
OUT_DIR <- file.path(OUT_ROOT, spec$folder)
stems <- if (length(spec$blocks) == 1) {
  spec$stem
} else {
  paste(spec$stem, spec$blocks, sep = "_")
}
stopifnot(
  "no output file exists yet" = !any(file.exists(file.path(
    OUT_DIR,
    paste0(stems, ".png")
  )))
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# One flowSet of CD8T events, arm x donor, linear channels on arcsinh --------------------------------
gs_path <- \(a) file.path(COMMON_GS, a)
frames <- ARMS |>
  imap(\(a, label) {
    gs <- load_gs(gs_path(a))
    sampleNames(gs) |>
      set_names(\(s) paste(label, s, sep = "|")) |>
      map(\(s) {
        cytoframe_to_flowFrame(gh_pop_get_data(gs[[s]], POP))[, unname(
          CH_USED
        )]
      })
  }) |>
  list_flatten(name_spec = "{inner}")
fs <- flowSet(frames)
fs <- transform(
  fs,
  transformList(
    LINEAR_TO_ARCSINH,
    arcsinhTransform(a = 0, b = 1 / COFACTOR, c = 0)
  )
)
pd <- tibble(name = sampleNames(fs)) |>
  tidyr::separate_wider_delim(
    name,
    "|",
    names = c("arm", "donor"),
    cols_remove = FALSE
  ) |>
  mutate(arm = factor(arm, levels = names(ARMS)))
pData(fs) <- as.data.frame(pd) |> `rownames<-`(pd$name)

# % of CD8T above the cut per arm and donor, printed in each panel ----------------------------------
# cut of one frame: its arm folder (ARMS value) and sample (donor) pick the cell of cuts.csv
cut_df <- pd |>
  mutate(
    cut = map2_dbl(arm, donor, \(a, d) {
      cut_of(ARMS[[as.character(a)]], d, Y_CH, "pos")
    })
  )
pct_pos <- cut_df |>
  mutate(
    pct_above_cut = pmap_dbl(
      list(name, cut),
      \(n, k) 100 * mean(exprs(fs[[n]])[, Y_CH] >= k)
    )
  ) |>
  select(name, pct_above_cut, arm, donor)
pct_table <- pct_pos |>
  mutate(
    arm_id = ARMS[as.character(arm)],
    arm = gsub("\n", " ", arm),
    marker = Y_NAME
  ) |>
  select(arm_id, arm, donor, marker, pct_above_cut)

# One figure per x marker: hexbin, fixed cut, %, donor x arm facets ----------------------------------
# Nature single column (89 mm) holds up to 4 arms, double column (183 mm) the rest; square panels,
# panel width = (figure width - ~18 mm for the y axis) / arms
fig_width_mm <- if (length(ARMS) <= 4) 89 else 183
fig_height_mm <- round(
  (fig_width_mm - 18) / length(ARMS) * n_distinct(pd$donor) + 28
)
walk2(spec$blocks, stems, \(x_name, stem) {
  x_ch <- CHANNELS[[x_name]]
  # the cut as a line: an open-ended gate would draw box edges at the display limits
  p <- ggcyto(fs, aes(x = !!rlang::sym(x_ch), y = !!rlang::sym(Y_CH))) +
    geom_hex(bins = BINS) +
    geom_hline(
      data = cut_df,
      aes(yintercept = cut),
      colour = "red",
      linewidth = FIG_LINE
    ) +
    geom_label(
      data = pct_pos,
      aes(
        x = -Inf,
        y = Inf,
        label = sprintf("%s+ %.1f%%", Y_NAME, pct_above_cut)
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = LABEL_PT,
      size.unit = "pt",
      family = FIG_FONT,
      fill = "white",
      alpha = 0.85,
      label.size = 0,
      label.padding = unit(0.1, "lines")
    ) +
    facet_grid(
      donor ~ arm,
      switch = "y",
      axes = "all",
      axis.labels = "margins"
    ) +
    ggcyto_par_set(limits = list(x = X_LIMITS[[x_name]], y = AXIS)) +
    # standard flow pseudocolor, as ggcyto's own default: blue (few events) to red (many)
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      trans = "sqrt",
      guide = "none"
    ) +
    labs(
      x = sprintf("%s, %s (%s)", LABELS[[x_name]], x_ch, SCALE_LABEL[[x_name]]),
      y = sprintf("%s, %s (%s)", Y_NAME, Y_CH, SCALE_LABEL[[spec$y]]),
      title = sprintf(
        "CD8 T cells: %s vs %s, red line = cut at %s",
        Y_NAME,
        LABELS[[x_name]],
        CUT_TXT
      )
    ) +
    theme(
      strip.text.y.left = element_text(face = "bold", angle = 0),
      strip.background.x = element_rect(colour = NA),
      strip.placement = "outside",
      panel.spacing = unit(1.2, "mm"),
      # both axes arcsinh (SSC aside), so the panels are square (r-0006)
      aspect.ratio = 1
    )
  # SSC ticks in millions, so the panels' labels do not overlap
  if (x_name == "SSC") {
    p <- p +
      scale_x_continuous(breaks = c(0, 1e6, 2e6), labels = c("0", "1", "2"))
  }
  # build on a cairo device, as save_panel() draws, so Nimbus Sans resolves through fontconfig
  cairo_pdf(tempfile(fileext = ".pdf"))
  g <- ggplotGrob(as.ggplot(p))
  dev.off()
  # fill each arm strip with the arm colour key; strip-t grobs run left to right in facet order
  strips <- grep("^strip-t", g$layout$name)
  strips <- strips[order(g$layout$l[strips])]
  for (i in seq_along(strips)) {
    fill <- ARM_FILL[[levels(pd$arm)[i]]]
    rect <- g$grobs[[strips[i]]]$grobs[[1]]$children[[1]]
    rect$gp$fill <- fill
    g$grobs[[strips[i]]]$grobs[[1]]$children[[1]] <- rect
    txt <- g$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]]
    txt$gp$col <- if (fill %in% DARK_FILLS) "white" else "black"
    g$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]] <- txt
  }
  save_panel(g, pct_table, OUT_DIR, stem, fig_width_mm, fig_height_mm)
  add_manifest(
    "95c_restyle_events_scatter",
    file.path("2026-11-11_bmv-restyle", spec$folder, paste0(stem, ".png")),
    file.path(OUT_DIR, paste0(stem, ".png")),
    COMMON_GS
  )
})
