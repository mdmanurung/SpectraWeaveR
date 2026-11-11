# Purpose: Nature-spec redraw on the common gates (plans/common_gates_plan.md G7a, G11) of the dashboard's
#          per-sample viewer: for each donor and parent (CD4T, CD8T), the SpectroFlo-labelled events of all 11
#          arms as ggcyto hexbins in three views, TNF x IFNg, SSC-A x TNF and SSC-A x IFNg, with the fixed
#          cut at arcsinh 1.5 and the % above it. Same cells in every arm; the % must match 38c's table, else
#          the script stops. Copy of 95_restyle_events_sample_viewer.R with new inputs.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_events_sample_viewer.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-11_bmv-common-gated/gs/<arm>/   common-gate GatingSets, 11 arms [97e]
#          experiments/runs/2026-11-11_bmv-common-cut15-38c/pct_fixed_source.csv   % above 1.5 [38c]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R   ARM_LABELS, VERSION_ROWS
# Outputs: experiments/runs/2026-11-11_bmv-restyle-common/08_dashboard_sample_viewer/
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
source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R")

# Parameters ---------------------------------------------------------------------------------------
COMMON_GS <- "experiments/runs/2026-11-11_bmv-common-gated/gs"
PCT_38 <- "experiments/runs/2026-11-11_bmv-common-cut15-38c/pct_fixed_source.csv"
CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A", SSC = "SSC-A") # TNF and IFNg stored as arcsinh (19a)
POPULATIONS <- c("CD4T", "CD8T")
AXIS <- c(-1.5, 5) # arcsinh display window, same on every panel and donor so SpectroFlo's negative cloud is not clipped; the % uses every event
SSC_LIMITS <- c(0, 2.5e6)
BINS <- 60 # hexagon bins per ~40 mm panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
LABEL_INSET <- 0.08 # % label starts this fraction of the axis range in from the y axis
N_COL <- 4 # 11 arms in 3 rows at 183 mm
VIEWS <- list(
  fixed = list(x = "TNF", y = "IFNg", stem = "tnf_vs_ifng_fixed"),
  ssc_TNF = list(x = "SSC", y = "TNF", stem = "ssc_vs_TNF"),
  ssc_IFNg = list(x = "SSC", y = "IFNg", stem = "ssc_vs_IFNg")
)

OUT_DIR <- fresh_dir("08_dashboard_sample_viewer")
arm_order <- unlist(VERSION_ROWS, use.names = FALSE)
gating_sets <- map(set_names(arm_order), \(a) {
  load_gs(file.path(COMMON_GS, a))
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
        dplyr::filter(
          pct_38,
          .data$donor == .env$donor,
          .data$parent == .env$pop
        ),
        join_by(arm_id == arm, cytokine),
        relationship = "one-to-one",
        unmatched = "error"
      )
    stopifnot("% matches 38" = max(abs(check$pct_here - check$pct_pos)) < 1e-9)

    iwalk(VIEWS, \(view, view_name) {
      x_ch <- CHANNELS[[view$x]]
      y_ch <- CHANNELS[[view$y]]
      x_lim <- if (view$x == "SSC") SSC_LIMITS else AXIS
      labels <- pct |>
        mutate(
          label_x = x_lim[1] + LABEL_INSET * diff(x_lim),
          label = if (view_name == "fixed") {
            sprintf("TNF+ %.1f%%\nIFNg+ %.1f%%", .data$TNF, .data$IFNg)
          } else {
            sprintf("%s+ %.1f%%", view$y, .data[[view$y]])
          }
        )
      p <- ggcyto(fs, aes(x = !!rlang::sym(x_ch), y = !!rlang::sym(y_ch))) +
        geom_hex(bins = BINS, show.legend = TRUE) +
        geom_hline(yintercept = CUT, colour = "red", linewidth = FIG_LINE) +
        geom_label(
          data = labels,
          aes(x = label_x, y = Inf, label = label),
          inherit.aes = FALSE,
          hjust = 0,
          vjust = 1,
          lineheight = 0.9,
          size = LABEL_PT,
          size.unit = "pt",
          family = FIG_FONT,
          fill = scales::alpha("white", 0.75),
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
          name = "events per hex"
        ) +
        labs(
          x = if (view$x == "SSC") {
            "SSC-A (\u00d710\u2076)"
          } else {
            sprintf("%s, %s (arcsinh)", view$x, x_ch)
          },
          y = sprintf("%s, %s (arcsinh)", view$y, y_ch),
          title = sprintf(
            "%s, donor %s: %s vs %s, red line = cut at %s",
            pop,
            donor,
            view$y,
            view$x,
            CUT_TXT
          )
        ) +
        theme(
          panel.spacing = unit(1.5, "mm"),
          aspect.ratio = 1,
          legend.position = "right"
        )
      # vertical cut for the TNF x IFNg view: IFNg is on y, TNF on x
      if (view_name == "fixed") {
        p <- p +
          geom_vline(xintercept = CUT, colour = "red", linewidth = FIG_LINE)
      }
      if (view$x == "SSC") {
        p <- p +
          scale_x_continuous(breaks = c(0, 1e6, 2e6), labels = c("0", "1", "2"))
      }
      stem <- sprintf("%s_%s_%s", view$stem, pop, donor)
      save_panel(
        as.ggplot(p),
        select(labels, arm_id, arm, TNF, IFNg),
        OUT_DIR,
        stem,
        183,
        170
      )
      add_manifest(
        "95c_restyle_events_sample_viewer",
        file.path(
          "2026-11-11_bmv-restyle/08_dashboard_sample_viewer",
          paste0(stem, ".png")
        ),
        file.path(OUT_DIR, paste0(stem, ".png")),
        paste(COMMON_GS, "; % = 38c pct_fixed_source.csv")
      )
    })
  }
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
