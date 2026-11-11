# Purpose: Nature-spec redraw (user, 2026-11-11; lucid-coder r-0006, r-0008, r-0016, r-0017) of the CD4 vs
#          CD8 gate check from 19_fsc_vs_cytokines.R for donor A3 (the one the report cites): CD3+ cells per
#          arm, the CD4T and CD8T gates each arm's GatingSet holds (quadGate.tmix, pooled per arm), and the %
#          of CD3+ in each gate from the gate counts. Same cells, gates and % as 19; the layout changes: one
#          facet grid instead of a cowplot page, theme_prism 7 pt Nimbus Sans, axes on every panel, labels
#          on the outer edges, one shared display window and square panels (both axes arcsinh), 183 mm wide.
#          The amber tint for the zip arms is dropped; "(zip)" in the arm label still marks them.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95_restyle_events_gate_cd4_cd8.R [donor, default A3]
#            from the repo root.
# Inputs:  experiments/runs/2026-11-10_bmv-gated-pooled/gs/<arm>/   GatingSets [19a_gate_arm.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  ARM_LABELS, VERSION_ROWS
# Outputs: experiments/runs/2026-11-11_bmv-restyle/02_gates/cd4_cd8_tmix_pooled/gate_cd4_cd8_A3.pdf, .png,
#            .csv   the .csv holds the gate counts and % of CD3+ per arm

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
OUT_DIR <- "experiments/runs/2026-11-11_bmv-restyle/02_gates/cd4_cd8_tmix_pooled"
DONOR <- if (length(commandArgs(TRUE)) > 0) commandArgs(TRUE)[1] else "A3" # A3 is the donor the report shows; the dashboard viewer needs all five
CD4 <- "Spark Blue 574-A" # channels as in 19a_gate_arm.R
CD8 <- "SBUV445-A"
AXIS <- c(-3, 6) # shared display window, arcsinh; the % use every CD3+ cell
BINS <- 60 # hexagons per panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
WIDTH_MM <- 183 # Nature double column, four panels per row
N_COL <- 4

STEM <- paste0("gate_cd4_cd8_", DONOR)
stopifnot(
  "no output file exists yet" = !file.exists(file.path(
    OUT_DIR,
    paste0(STEM, ".png")
  ))
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# CD3+ events, gates and gate counts per arm, in version-row order ---------------------------------------
arms <- unlist(VERSION_ROWS, use.names = FALSE)
arm_levels <- unname(ARM_LABELS[arms])
gating_sets <- map(set_names(arms), \(a) load_gs(file.path(GATED_POOLED, a)))
frames <- map(gating_sets, \(gs) {
  cytoframe_to_flowFrame(gh_pop_get_data(gs[[DONOR]], "CD3pos"))[, c(CD4, CD8)]
})
fs <- flowSet(frames)
pData(fs)$arm <- factor(ARM_LABELS[sampleNames(fs)], levels = arm_levels)
# one gate per arm, each from its own GatingSet
gates_of <- \(pop) {
  filterList(map(gating_sets, \(gs) gh_pop_get_gate(gs[[DONOR]], pop)))
}
pct <- imap(gating_sets, \(gs, a) {
  gs_pop_get_count_fast(gs[DONOR]) |>
    as_tibble() |>
    mutate(arm_id = a)
}) |>
  list_rbind() |>
  dplyr::filter(
    Population %in%
      c(
        "/nonDebris/singlets/live/CD3pos/CD4T",
        "/nonDebris/singlets/live/CD3pos/CD8T"
      )
  ) |>
  mutate(
    gate = basename(Population),
    pct_of_cd3 = 100 * Count / ParentCount,
    arm = factor(ARM_LABELS[arm_id], levels = arm_levels),
    # CD8T top left, CD4T bottom right, as the quadrants sit
    x = if_else(gate == "CD8T", -Inf, Inf),
    y = if_else(gate == "CD8T", Inf, -Inf),
    hjust = if_else(gate == "CD8T", 0, 1),
    vjust = if_else(gate == "CD8T", 1, 0)
  )

# Figure ---------------------------------------------------------------------------------------------------
p <- ggcyto(fs, aes(x = !!rlang::sym(CD4), y = !!rlang::sym(CD8))) +
  geom_hex(bins = BINS) +
  geom_gate(gates_of("CD4T"), colour = "red", linewidth = FIG_LINE) +
  geom_gate(gates_of("CD8T"), colour = "red", linewidth = FIG_LINE) +
  geom_label(
    data = pct,
    aes(
      x = x,
      y = y,
      hjust = hjust,
      vjust = vjust,
      label = sprintf("%s %.1f%%", gate, pct_of_cd3)
    ),
    inherit.aes = FALSE,
    size = LABEL_PT,
    size.unit = "pt",
    family = FIG_FONT,
    fill = "white",
    alpha = 0.85,
    label.size = 0,
    label.padding = unit(0.1, "lines")
  ) +
  facet_wrap(~arm, ncol = N_COL, axes = "all", axis.labels = "margins") +
  ggcyto_par_set(limits = list(x = AXIS, y = AXIS)) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    trans = "sqrt",
    guide = "none"
  ) +
  labs(
    x = sprintf("CD4, %s (arcsinh)", CD4),
    y = sprintf("CD8, %s (arcsinh)", CD8),
    title = sprintf(
      "CD3+ T cells, donor %s: CD4T and CD8T gates per arm",
      DONOR
    )
  ) +
  theme(
    panel.spacing = unit(1.2, "mm"),
    # both axes arcsinh, so the panels are square (r-0006)
    aspect.ratio = 1
  )
# build on a cairo device, as save_panel() draws, so Nimbus Sans resolves through fontconfig
cairo_pdf(tempfile(fileext = ".pdf"))
g <- ggplotGrob(as.ggplot(p))
dev.off()
n_rows <- ceiling(length(arms) / N_COL)
height_mm <- round((WIDTH_MM - 18) / N_COL * n_rows + 32)
save_panel(
  g,
  select(pct, arm_id, arm, gate, Count, ParentCount, pct_of_cd3),
  OUT_DIR,
  STEM,
  WIDTH_MM,
  height_mm
)
