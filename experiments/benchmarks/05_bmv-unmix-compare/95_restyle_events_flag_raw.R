# Purpose: Nature-spec redraw (user, 2026-11-11; lucid-coder r-0006, r-0008, r-0010, r-0016, r-0017) of the
#          flag raw companions from 93_flag_raw_scatter.R that the report cites (IFNg in CD4T, Th2 in CD8T):
#          flagged marker against GLUT1, live rows paired on as188_joint's gates, cut at arcsinh 1, % above
#          the cut per panel. Same cells and numbers as 93; only the layout changes: theme_prism 7 pt Nimbus
#          Sans, axes on every panel, labels on the outer edges, square panels, arm strips filled with the
#          arm colour key (as the 72 scatters), 183 mm wide.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_events_flag_raw.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-11_bmv-flag-raw-scatter/flagged_markers.csv   marker, population [93]
#          experiments/runs/2026-11-10_bmv-ablation-gated/gs/{as188_joint,opt188_pregate_joint_p1af1,as102_var}
#          experiments/runs/2026-11-10_bmv-gated-pooled/gs/{spectroflo,s188_c102}   GatingSets [19a]
# Outputs: experiments/runs/2026-11-11_bmv-restyle/21_raw_companions/flag_raw/fig_flag_<marker>_<pop>.pdf,
#            .png, .csv   the .csv holds the % above the cut per arm and donor

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

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- "experiments/runs/2026-11-11_bmv-restyle/21_raw_companions/flag_raw"
GATED <- "experiments/runs/2026-11-10_bmv-ablation-gated/gs"
GATED_POOLED <- "experiments/runs/2026-11-10_bmv-gated-pooled/gs"
CITED <- c("IFNg", "Th2") # flagged markers the report shows (90_report.qmd section 21)
# column label = GatingSet; references first, then the flagged arms (as 93)
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
POP_PRIORITY <- c("CD8T", "CD4T", "live")
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
BINS <- 50 # hexagons per panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
WIDTH_MM <- 183 # Nature double column, five arms

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Cited flagged markers and the population each is shown in -----------------------------------------------
flagged <- read_csv(
  "experiments/runs/2026-11-11_bmv-flag-raw-scatter/flagged_markers.csv",
  show_col_types = FALSE
) |>
  dplyr::filter(marker %in% CITED)
stems <- sprintf(
  "fig_flag_%s_%s",
  gsub("[^A-Za-z0-9]", "", flagged$marker),
  flagged$pop
)
stopifnot(
  "no output file exists yet" = !any(file.exists(file.path(
    OUT_DIR,
    paste0(stems, ".png")
  )))
)

# Live cells of each arm and donor, rows paired on as188_joint's gates (as 93) ----------------------------
gs_list <- map(ARMS, load_gs)
gs_ref <- gs_list[["1.8.8 zip"]]
donors <- sampleNames(gs_ref)
channels <- unique(c(GLUT1, flagged$channel))
live <- map(set_names(donors), \(d) {
  which(gh_pop_get_indices(gs_ref[[d]], "live"))
})
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

# One ggcyto figure per cited flagged marker -------------------------------------------------------------
pwalk(
  list(flagged$marker, flagged$channel, flagged$pop, flagged$flagged_by, stems),
  \(which_marker, which_channel, which_pop, flagged_by, stem) {
    pct_above <- tidyr::expand_grid(arm = names(ARMS), donor = donors) |>
      mutate(
        marker = which_marker,
        pop = which_pop,
        pct_above_cut = map2_dbl(arm, donor, \(a, d) {
          100 *
            mean(
              cells[[a]][[d]][in_pop[[d]][[which_pop]], which_channel] >= CUT
            )
        }),
        arm = factor(arm, levels = names(ARMS))
      )
    frames <- tidyr::expand_grid(arm = names(ARMS), donor = donors) |>
      mutate(name = paste(arm, donor, sep = "|")) |>
      with(set_names(
        map2(arm, donor, \(a, d) {
          flowFrame(cells[[a]][[d]][
            in_pop[[d]][[which_pop]],
            c(GLUT1, which_channel)
          ])
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
    flagged_labels <- FLAG_ARM[strsplit(flagged_by, "; ")[[1]]]
    p <- ggcyto(
      fs,
      aes(x = !!rlang::sym(GLUT1), y = !!rlang::sym(which_channel))
    ) +
      geom_hex(bins = BINS) +
      # the cut as a line: an open-ended gate would draw box edges at the display limits
      geom_hline(yintercept = CUT, colour = "red", linewidth = FIG_LINE) +
      geom_label(
        data = pct_above,
        aes(
          x = -Inf,
          y = Inf,
          label = sprintf("%s+ %.1f%%", which_marker, pct_above_cut)
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
      ggcyto_par_set(limits = list(x = AXIS, y = AXIS)) +
      scale_fill_gradientn(
        colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
        trans = "sqrt",
        guide = "none"
      ) +
      labs(
        x = sprintf("GLUT1, %s (arcsinh, cofactor 6000)", GLUT1),
        y = sprintf(
          "%s, %s (arcsinh, cofactor 6000)",
          which_marker,
          which_channel
        ),
        title = sprintf(
          "%s: %s vs GLUT1 (flagged in %s)",
          which_pop,
          which_marker,
          paste(gsub("\n", " ", flagged_labels), collapse = ", ")
        )
      ) +
      theme(
        strip.text.y.left = element_text(face = "bold", angle = 0),
        strip.background.x = element_rect(colour = NA),
        strip.placement = "outside",
        panel.spacing = unit(1.2, "mm"),
        # both axes arcsinh, so the panels are square (r-0006)
        aspect.ratio = 1
      )
    # build on a cairo device, as save_panel() draws, so Nimbus Sans resolves through fontconfig
    cairo_pdf(tempfile(fileext = ".pdf"))
    g <- ggplotGrob(as.ggplot(p))
    dev.off()
    # fill each arm strip with the arm colour key; strip-t grobs run left to right in facet order
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
    height_mm <- round((WIDTH_MM - 18) / length(ARMS) * length(donors) + 28)
    save_panel(
      g,
      mutate(pct_above, arm = gsub("\n", " ", arm)),
      OUT_DIR,
      stem,
      WIDTH_MM,
      height_mm
    )
  }
)
