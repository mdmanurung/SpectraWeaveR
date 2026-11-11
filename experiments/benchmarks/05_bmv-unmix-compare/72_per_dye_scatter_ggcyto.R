# Purpose: cytokine scatter per arm set with ggcyto: hexbin density (ggcyto's default pseudocolor, reversed
#          Brewer Spectral, sqrt counts), the fixed cut (arcsinh 1) drawn with geom_gate and the % above it as
#          text, column strips coloured by arm (user, 2026-11-10). CD8T, each arm on its own gates. Rows are
#          donors and columns arms, so one donor reads across methods. Blocks side by side, the y marker
#          against each x marker. Arm sets:
#            per_dye  y TNF; 1.0.2 variants swapped in one dye at a time, then SBUV510 + AF405 (45); x: IFNg, GLUT1, SSC
#            af_swap  y TNF; AF spectra and variants swapped between versions (37, 48); x: IFNg, AF, SSC
#            tnf_new  y TNF; arms from tests 3, Q4/Q9 and the downstream comparison; x: IFNg, GLUT1, SSC
#            th2      y Th2 (eFluor 450); arms that squeeze Th2 and the AF405 knock-out; x: GLUT1, IFNg, SSC
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/72_per_dye_scatter_ggcyto.R <set>
#            from the repo root.
# Inputs:  experiments/runs/2026-11-10_bmv-ablation-gated/gs/<arm>/   GatingSets [19a via 51]
#          experiments/runs/2026-11-10_bmv-gated-pooled/gs/{s188_c102,spectroflo}/  GatingSets [19a]
#          BV510/BV570 stored arcsinh; Alexa Fluor 405, AF-A and eFluor 450 stored linear
# Outputs: <out_dir>/fig_<y>_vs_<x>.pdf, .png              one figure per x marker, sized for print
#          <out_dir>/pct_above_cut.csv                  % of CD8T above the cut, per arm and donor
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
SET <- commandArgs(trailingOnly = TRUE)[1]
GATED <- "experiments/runs/2026-11-10_bmv-ablation-gated/gs"
GATED_POOLED <- "experiments/runs/2026-11-10_bmv-gated-pooled/gs" # holds s188_c102 and spectroflo
# arms per set (label = arm folder), strip colour per arm, y marker, x blocks, output folder
SETS <- list(
  per_dye = list(
    arms = c(
      "1.8.8\nvariants" = "as188_joint",
      "opt188\n(optimized)" = "opt188_pregate_joint_p1af1",
      "+SBUV510\n1.0.2" = "var102_only_SBUV510",
      "+ AF405\n1.0.2" = "var102_only_Alexa_Fluor_405",
      "+SBUV510\n+AF405 1.0.2" = "var102_only_SBUV510__Alexa_Fluor_405",
      "all 1.0.2\nvariants" = "swap_var102"
    ),
    # grey baseline, light to dark red as more of the TNF squeeze is swapped in
    fill = c("#d9d9d9", "#9ecae1", "#fdbb84", "#e34a33", "#b30000", "#7f0000"),
    y = "TNF",
    blocks = c("IFNg", "GLUT1", "SSC"),
    out_dir = "experiments/runs/2026-11-10_bmv-per-dye-ggcyto-v8"
  ),
  af_swap = list(
    arms = c(
      "1.8.8\nall" = "as188_joint",
      "+AF\n1.0.2" = "swap_af102",
      "+variants\n1.0.2" = "swap_var102",
      "+AF+var\n1.0.2" = "swap_af102_var102",
      "1.0.2, AF\n1.8.8" = "rev_af188",
      "1.0.2\nall" = "s102_c188"
    ),
    # 1.8.8 spectra: grey to dark red as 1.0.2 parts go in; 1.0.2 spectra: blues
    fill = c("#d9d9d9", "#fdbb84", "#e34a33", "#7f0000", "#9ecae1", "#08519c"),
    y = "TNF",
    blocks = c("IFNg", "AF", "SSC"),
    out_dir = "experiments/runs/2026-11-10_bmv-af-swap-ggcyto-v3"
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
    # 1.8.8 arms grey/blue, 1.0.2-code or 1.0.2-variant arms red
    fill = c("#d9d9d9", "#9ecae1", "#08519c", "#fdbb84", "#e34a33", "#7f0000"),
    y = "TNF",
    blocks = c("IFNg", "GLUT1", "SSC"),
    out_dir = "experiments/runs/2026-11-10_bmv-tnf-new-ggcyto-v3"
  ),
  th2 = list(
    arms = c(
      "1.8.8\nzip" = "as188_joint",
      "opt188" = "opt188_pregate_joint_p1af1",
      "+ AF405\n1.0.2" = "var102_only_Alexa_Fluor_405",
      "AF405 knock-out\n(no BV510 part)" = "af405proof_ko_bv510",
      "1.0.2 code\n1.8.8 spectra" = "s188_c102",
      "1.0.2\nzip" = "as102_var"
    ),
    fill = c("#d9d9d9", "#9ecae1", "#e34a33", "#fdbb84", "#08519c", "#7f0000"),
    y = "Th2",
    blocks = c("GLUT1", "IFNg", "SSC"),
    out_dir = "experiments/runs/2026-11-10_bmv-th2-ggcyto-v3"
  )
)
# SpectroFlo 3.3.0 export, first column of every set as the vendor reference (user, 2026-11-10)
REFERENCE <- c("SpectroFlo\n(reference)" = "spectroflo")
REFERENCE_FILL <- "#252525"
DARK_FILLS <- c("#e34a33", "#b30000", "#7f0000", "#08519c", REFERENCE_FILL) # strips that take white text
POP <- "CD8T"
CUT <- 1 # fixed cytokine cut, arcsinh (user, 2026-11-10)
COFACTOR <- 6000 # arcsinh cofactor for linear channels, as for the cytokines (CLAUDE.md)
AXIS <- c(-1.5, 4) # display window only; the % above the cut uses every event
BINS <- 60 # hexagon bins per panel; display only, fewer for the ~27 mm panels
PANEL_MM <- 23 # panel column width; seven columns fill a 183 mm double column
LABEL_PT <- 7 # in-panel % and strip text; house style floor is 6 pt at final width
CHANNELS <- c(
  TNF = "BV510-A",
  Th2 = "eFluor 450-A",
  IFNg = "BV570-A",
  GLUT1 = "Alexa Fluor 405-A",
  AF = "AF-A",
  SSC = "SSC-A"
)
LINEAR_TO_ARCSINH <- c("Alexa Fluor 405-A", "AF-A", "eFluor 450-A")
SCALE_LABEL <- c(
  TNF = "arcsinh",
  Th2 = "arcsinh, cofactor 6000",
  IFNg = "arcsinh, cofactor 6000",
  GLUT1 = "arcsinh, cofactor 6000",
  AF = "arcsinh, cofactor 6000",
  SSC = "linear, x 10^6"
)
# x display window per block; SSC 0-2.5e6 holds the CD8T cloud (~1e6) without the few outliers to 4e6
X_LIMITS <- list(IFNg = AXIS, GLUT1 = AXIS, AF = AXIS, SSC = c(0, 2.5e6))

stopifnot("set is one of the arm sets" = SET %in% names(SETS))
ARMS <- c(REFERENCE, SETS[[SET]]$arms)
ARM_FILL <- set_names(c(REFERENCE_FILL, SETS[[SET]]$fill), names(ARMS))
Y_NAME <- SETS[[SET]]$y
Y_CH <- CHANNELS[[Y_NAME]]
BLOCKS <- SETS[[SET]]$blocks
OUT_DIR <- SETS[[SET]]$out_dir

stopifnot(
  "output dir must not exist" = !dir.exists(OUT_DIR),
  "every arm has a strip colour" = setequal(names(ARMS), names(ARM_FILL))
)
dir.create(OUT_DIR)

# One flowSet of CD8T events, arm x donor, linear channels on arcsinh --------------------------------
gs_path <- \(a) {
  if (a %in% c("s188_c102", "spectroflo")) {
    file.path(GATED_POOLED, a)
  } else {
    file.path(GATED, a)
  }
}
frames <- ARMS |>
  imap(\(a, label) {
    gs <- load_gs(gs_path(a))
    sampleNames(gs) |>
      set_names(\(s) paste(label, s, sep = "|")) |>
      map(\(s) {
        cytoframe_to_flowFrame(gh_pop_get_data(gs[[s]], POP))[, unname(
          CHANNELS
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
pct_pos <- fsApply(fs, \(ff) 100 * mean(exprs(ff)[, Y_CH] >= CUT)) |>
  as.data.frame() |>
  tibble::rownames_to_column("name") |>
  rename(pct_above_cut = V1) |>
  left_join(pd, by = "name", relationship = "one-to-one", unmatched = "error")
write_csv(
  pct_pos |>
    mutate(
      arm_id = ARMS[as.character(arm)],
      arm = gsub("\n", " ", arm),
      marker = Y_NAME
    ) |>
    select(arm_id, arm, donor, marker, pct_above_cut),
  file.path(OUT_DIR, "pct_above_cut.csv")
)

# One ggcyto block per x marker: hexbin, fixed cut, %, donor x arm facets ----------------------------
block <- function(x_name) {
  x_ch <- CHANNELS[[x_name]]
  gate <- rectangleGate(
    filterId = "above cut",
    .gate = setNames(list(c(-Inf, Inf), c(CUT, Inf)), c(x_ch, Y_CH))
  )
  # ggcyto on a flowSet takes one gate per sample, named by sample
  gate <- filterList(set_names(rep(list(gate), length(fs)), sampleNames(fs)))
  # % from pct_pos, every event (geom_stats fails on the SSC block, so all blocks use this label)
  p <- ggcyto(fs, aes(x = !!rlang::sym(x_ch), y = !!rlang::sym(Y_CH))) +
    geom_hex(bins = BINS) +
    geom_gate(gate, colour = "red") +
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
      fontface = "bold",
      fill = "white",
      alpha = 0.85,
      label.size = 0,
      label.padding = unit(0.12, "lines")
    ) +
    facet_grid(donor ~ arm, switch = "y") +
    ggcyto_par_set(limits = list(x = X_LIMITS[[x_name]], y = AXIS)) +
    # standard flow pseudocolor, as ggcyto's own default: blue (few events) to red (many)
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      trans = "sqrt",
      guide = "none"
    ) +
    labs(
      x = sprintf("%s, %s (%s)", x_name, x_ch, SCALE_LABEL[[x_name]]),
      y = sprintf("%s, %s (%s)", Y_NAME, Y_CH, SCALE_LABEL[[Y_NAME]]),
      title = sprintf("CD8 T cells: %s vs %s", Y_NAME, x_name)
    ) +
    ggpubr::theme_pubr(base_size = 8, border = TRUE) +
    # donor strips left, outside, no box (house style); arm strips keep their fill
    theme(
      strip.text.x = element_text(size = LABEL_PT),
      strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      axis.text = element_text(size = LABEL_PT),
      panel.spacing = unit(1, "mm")
    )
  # SSC ticks in millions, so the panels' labels do not overlap
  if (x_name == "SSC") {
    p <- p +
      scale_x_continuous(breaks = c(0, 1e6, 2e6), labels = c("0", "1", "2"))
  }
  g <- ggplotGrob(as.ggplot(p))
  # colour each column strip by its arm; strip-t grobs run left to right in facet order
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
  g
}

# Figure ---------------------------------------------------------------------------------------------
# one figure per x marker, sized so text reads at print width (6 arms = 183 mm double column)
fig_width_mm <- PANEL_MM * length(ARMS) + 20
fig_height_mm <- PANEL_MM * n_distinct(pd$donor) + 35
walk(BLOCKS, \(x_name) {
  g <- block(x_name)
  stem <- file.path(OUT_DIR, sprintf("fig_%s_vs_%s", Y_NAME, x_name))
  ggsave(
    paste0(stem, ".pdf"),
    g,
    width = fig_width_mm,
    height = fig_height_mm,
    units = "mm",
    device = cairo_pdf
  )
  ggsave(
    paste0(stem, ".png"),
    g,
    width = fig_width_mm,
    height = fig_height_mm,
    units = "mm",
    dpi = 300,
    bg = "white"
  )
})

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
