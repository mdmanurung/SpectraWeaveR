# Purpose: Nature-spec redraw (user, 2026-11-11; lucid-coder r-0008, r-0015..r-0017) of the dashboard's
#          per-sample viewer: for each donor and parent (CD4T, CD8T), CD8T/CD4T events of all 11 arms as
#          ggcyto hexbins in three views, TNF x IFNg, SSC-A x TNF and SSC-A x IFNg, with the fixed cut at
#          arcsinh 1 and the % above it. Same GatingSets, cells and cut as 38_fixed_cut.R and
#          40_ssc_fixed_cut.R; the % printed must match 38's table, else the script stops.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_events_sample_viewer.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-10_bmv-gated-pooled/gs/<arm>/   GatingSets, 11 arms [19a_gate_arm.R]
#          experiments/runs/2026-11-10_bmv-figures-v2/gated-pooled/figures_fixed_cut_1/pct_fixed_source.csv
#            % above cut 1 per arm, donor, parent and cytokine [38_fixed_cut.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R   ARM_LABELS, VERSION_ROWS
# Outputs: experiments/runs/2026-11-11_bmv-restyle/08_dashboard_sample_viewer/
#            tnf_vs_ifng_fixed_<parent>_<donor>.{pdf,png,csv}
#            ssc_vs_<TNF|IFNg>_<parent>_<donor>.{pdf,png,csv}   the .csv holds the % printed per panel

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
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R")

# Parameters ---------------------------------------------------------------------------------------
GATED_POOLED <- "experiments/runs/2026-11-10_bmv-gated-pooled/gs"
PCT_38 <- "experiments/runs/2026-11-10_bmv-figures-v2/gated-pooled/figures_fixed_cut_1/pct_fixed_source.csv"
OUT_DIR <- "experiments/runs/2026-11-11_bmv-restyle/08_dashboard_sample_viewer"
CUT <- 1 # fixed cytokine cut, arcsinh (user, 2026-11-10)
CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A", SSC = "SSC-A") # TNF and IFNg stored as arcsinh (19a)
POPULATIONS <- c("CD4T", "CD8T")
AXIS <- c(-0.5, 5) # arcsinh display window, as 38; the % uses every event
SSC_LIMITS <- c(0, 2.5e6)
BINS <- 60 # hexagon bins per ~40 mm panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
N_COL <- 4 # 11 arms in 3 rows at 183 mm
VIEWS <- list(
  fixed = list(x = "TNF", y = "IFNg", stem = "tnf_vs_ifng_fixed"),
  ssc_TNF = list(x = "SSC", y = "TNF", stem = "ssc_vs_TNF"),
  ssc_IFNg = list(x = "SSC", y = "IFNg", stem = "ssc_vs_IFNg")
)

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE)
arm_order <- unlist(VERSION_ROWS, use.names = FALSE)
gating_sets <- map(set_names(arm_order), \(a) {
  load_gs(file.path(GATED_POOLED, a))
})
donors <- sampleNames(gating_sets[[1]])
pct_38 <- read_csv(PCT_38, show_col_types = FALSE)

# One figure per parent, donor and view ---------------------------------------------------------------
for (pop in POPULATIONS) {
  for (donor in donors) {
    # one flowFrame per arm: this donor's parent events
    fs <- flowSet(map(gating_sets, \(gs) {
      cytoframe_to_flowFrame(gh_pop_get_data(gs[[donor]], pop))[, unname(
        CHANNELS
      )]
    }))
    pData(fs)$arm <- factor(
      ARM_LABELS[sampleNames(fs)],
      levels = ARM_LABELS[arm_order]
    )
    # % above the cut per arm, checked against 38's table
    pct <- tibble(
      arm_id = sampleNames(fs),
      arm = pData(fs)$arm,
      TNF = fsApply(fs, \(ff) {
        100 * mean(exprs(ff)[, CHANNELS[["TNF"]]] >= CUT)
      })[, 1],
      IFNg = fsApply(fs, \(ff) {
        100 * mean(exprs(ff)[, CHANNELS[["IFNg"]]] >= CUT)
      })[, 1]
    )
    check <- pct |>
      tidyr::pivot_longer(
        c(TNF, IFNg),
        names_to = "cytokine",
        values_to = "pct_here"
      ) |>
      inner_join(
        dplyr::filter(pct_38, .data$donor == .env$donor, .data$parent == .env$pop),
        join_by(arm_id == arm, cytokine),
        relationship = "one-to-one",
        unmatched = "error"
      )
    stopifnot("% matches 38" = max(abs(check$pct_here - check$pct_pos)) < 1e-9)

    iwalk(VIEWS, \(view, view_name) {
      x_ch <- CHANNELS[[view$x]]
      y_ch <- CHANNELS[[view$y]]
      labels <- pct |>
        mutate(
          label = if (view_name == "fixed") {
            sprintf("TNF+ %.1f%%\nIFNg+ %.1f%%", .data$TNF, .data$IFNg)
          } else {
            sprintf("%s+ %.1f%%", view$y, .data[[view$y]])
          }
        )
      p <- ggcyto(fs, aes(x = !!rlang::sym(x_ch), y = !!rlang::sym(y_ch))) +
        geom_hex(bins = BINS) +
        geom_hline(yintercept = CUT, colour = "red", linewidth = FIG_LINE) +
        geom_label(
          data = labels,
          aes(x = -Inf, y = Inf, label = label),
          inherit.aes = FALSE,
          hjust = 0,
          vjust = 1,
          lineheight = 0.9,
          size = LABEL_PT,
          size.unit = "pt",
          family = FIG_FONT,
          fill = "white",
          alpha = 0.85,
          label.size = 0,
          label.padding = unit(0.1, "lines")
        ) +
        facet_wrap(~arm, ncol = N_COL, axes = "all", axis.labels = "margins") +
        ggcyto_par_set(
          limits = list(x = if (view$x == "SSC") SSC_LIMITS else AXIS, y = AXIS)
        ) +
        scale_fill_gradientn(
          colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
          trans = "sqrt",
          guide = "none"
        ) +
        labs(
          x = sprintf(
            "%s, %s (%s)",
            view$x,
            x_ch,
            if (view$x == "SSC") "linear, x 10^6" else "arcsinh"
          ),
          y = sprintf("%s, %s (arcsinh)", view$y, y_ch),
          title = sprintf("%s, donor %s: %s vs %s", pop, donor, view$y, view$x)
        ) +
        theme(panel.spacing = unit(1.5, "mm"), aspect.ratio = 1)
      # vertical cut for the TNF x IFNg view: IFNg is on y, TNF on x
      if (view_name == "fixed") {
        p <- p +
          geom_vline(xintercept = CUT, colour = "red", linewidth = FIG_LINE)
      }
      if (view$x == "SSC") {
        p <- p +
          scale_x_continuous(breaks = c(0, 1e6, 2e6), labels = c("0", "1", "2"))
      }
      save_panel(
        as.ggplot(p),
        select(labels, arm_id, arm, TNF, IFNg),
        OUT_DIR,
        sprintf("%s_%s_%s", view$stem, pop, donor),
        183,
        170
      )
    })
  }
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
