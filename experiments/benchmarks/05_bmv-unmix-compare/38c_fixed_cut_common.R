# Purpose: % TNF+ and % IFNg+ of CD4T and CD8T per arm and donor at ONE fixed cut (arcsinh 1.5, cofactor 6000)
#          on the common SpectroFlo gates (plan common_gates_plan.md: G1a, G7a, G10b). Every arm uses the same
#          cell labels, so differences come from unmixing only. Copy of 38_fixed_cut.R: new gated dir, cut
#          1.5, new run dir; figures redrawn in the Nature style with donors as points and a median bar.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/38c_fixed_cut_common.R <run_dir>
#          from the repo root, after batch 0 has written GATED_DIR/README.md. <run_dir> must not exist.
# Inputs:  GATED_DIR/manifest.csv, GATED_DIR/gs/<arm>/  common-gate GatingSets [97_common_gs.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R, report_figures.R  labels, builders
# Outputs: <run_dir>/pct_fixed_source.csv                 % above the cut per arm, donor, parent, cytokine
#          <run_dir>/fig_pct_fixed.{pdf,png,csv}          donors as points, median bar
#          <run_dir>/tnf_vs_ifng_fixed_<population>_<donor>.{pdf,png,csv}  TNF against IFNg hexbins, cut drawn
#          <run_dir>/sessionInfo_38c.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})
source("experiments/benchmarks/05_bmv-unmix-compare/report_figures.R") # theme, save_panel, builders, ARM_LABELS

# Parameters ------------------------------------------------------------------------------
GATED_DIR <- "experiments/runs/2026-11-11_bmv-common-gated" # batch 0 output (G9: read-only here)
CUT <- 1.5 # arcsinh cut for both cytokines (G7a, user 2026-11-11; replaces 1 in 38_fixed_cut.R)
CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A") # as in 19a_gate_arm.R; stored as arcsinh
NODES <- c(
  CD4T = "/nonDebris/singlets/live/CD3pos/CD4T",
  CD8T = "/nonDebris/singlets/live/CD3pos/CD8T"
)
BINS <- 60 # hexagon bins per axis; display only
AXIS_LIMITS <- c(-0.5, 5) # arcsinh display range (user, 2026-11-09: trims are plots only)
ARM_ORDER <- unlist(VERSION_ROWS, use.names = FALSE) # the 11 figure arms, in figure order

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("38c run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)
manifest <- read_csv(
  file.path(GATED_DIR, "manifest.csv"),
  show_col_types = FALSE
)
stopifnot(
  "all 11 figure arms are in the gated manifest" = all(
    ARM_ORDER %in% manifest$arm
  )
)

# Events of each arm, donor and population: the two cytokines --------------------------------
events <- ARM_ORDER |>
  set_names() |>
  map(\(a) {
    gs <- load_gs(file.path(GATED_DIR, "gs", a))
    sampleNames(gs) |>
      set_names() |>
      map(\(s) {
        NODES |>
          map(\(node) {
            x <- exprs(gh_pop_get_data(gs[[s]], node))[, CHANNELS, drop = FALSE]
            colnames(x) <- names(CHANNELS)
            as_tibble(x)
          }) |>
          list_rbind(names_to = "parent")
      }) |>
      list_rbind(names_to = "donor")
  }) |>
  list_rbind(names_to = "arm")

# % above the cut per arm, donor, parent and cytokine -------------------------------------------
pct <- events |>
  summarise(
    n_parent = n(),
    TNF = 100 * mean(TNF >= CUT),
    IFNg = 100 * mean(IFNg >= CUT),
    .by = c(arm, donor, parent)
  ) |>
  tidyr::pivot_longer(
    c(TNF, IFNg),
    names_to = "cytokine",
    values_to = "pct_pos"
  ) |>
  mutate(
    arm = factor(arm, levels = ARM_ORDER),
    cytokine = factor(cytokine, levels = names(CHANNELS))
  ) |>
  arrange(arm, donor, parent, cytokine) |>
  mutate(arm = as.character(arm), cytokine = as.character(cytokine)) |>
  select(arm, donor, parent, cytokine, n_parent, pct_pos)
stopifnot(
  "one value per arm, donor, parent, cytokine" = nrow(pct) ==
    length(ARM_ORDER) * 5 * 2 * 2
)
# Identity check: my counts above the cut equal the saved TNFpos / IFNgpos gates (fixed 1.5, README)
node_counts <- ARM_ORDER |>
  set_names() |>
  map(\(a) {
    gs_pop_get_count_fast(load_gs(file.path(GATED_DIR, "gs", a))) |>
      as_tibble() |>
      dplyr::filter(str_detect(Population, "/(CD4T|CD8T)/(TNF|IFNg)pos$"))
  }) |>
  list_rbind(names_to = "arm") |>
  transmute(
    arm,
    donor = name,
    parent = str_extract(Population, "CD[48]T"),
    cytokine = str_extract(Population, "(TNF|IFNg)(?=pos$)"),
    n_node = Count
  )
stopifnot(
  "counts above the cut equal the TNFpos / IFNgpos node counts" = nrow(
    inner_join(
      mutate(pct, n_above = round(pct_pos / 100 * n_parent)),
      node_counts,
      by = join_by(arm, donor, parent, cytokine),
      relationship = "one-to-one"
    ) |>
      dplyr::filter(n_above != n_node)
  ) ==
    0
)
write_csv(pct, file.path(run_dir, "pct_fixed_source.csv"))

# Figure 1: % positive per arm, donors as points over a median bar ------------------------------
redraw_panel(
  file.path(run_dir, "pct_fixed_source.csv"),
  "pct_pos",
  paste0("% positive at arcsinh ", CUT, ", common SpectroFlo gates"),
  run_dir,
  "fig_pct_fixed"
)

# Figure 2: TNF vs IFNg hexbins, one file per donor and parent, cut drawn ------------------------
pct_label <- pct |>
  pivot_wider(
    names_from = cytokine,
    values_from = pct_pos,
    names_prefix = "pct_"
  ) |>
  mutate(
    panel = sprintf(
      "%s\nTNF+ %.1f%%  IFNg+ %.1f%%",
      ARM_LABELS[arm],
      pct_TNF,
      pct_IFNg
    ),
    panel = factor(
      as.character(panel),
      levels = unique(as.character(panel))
    )
  )
for (s in sort(unique(events$donor))) {
  for (pop in names(NODES)) {
    plot_df <- events |>
      dplyr::filter(donor == s, parent == pop) |>
      left_join(
        dplyr::filter(pct_label, donor == s, parent == pop),
        by = join_by(arm, donor, parent),
        relationship = "many-to-one"
      ) |>
      mutate(
        panel = factor(
          as.character(panel),
          levels = levels(pct_label$panel)
        ) |>
          droplevels()
      )
    p <- ggplot(plot_df, aes(TNF, IFNg)) +
      geom_hex(bins = BINS) +
      geom_vline(xintercept = CUT, colour = "red", linewidth = 0.3) +
      geom_hline(yintercept = CUT, colour = "red", linewidth = 0.3) +
      scale_fill_viridis_c(
        option = "magma",
        transform = "log10",
        name = "Cells"
      ) +
      coord_cartesian(xlim = AXIS_LIMITS, ylim = AXIS_LIMITS) +
      facet_wrap(vars(panel), ncol = 4, axes = "all", axis.labels = "margins") +
      labs(
        x = "TNF (arcsinh)",
        y = "IFNg (arcsinh)",
        title = paste(pop, s, "- cut", CUT)
      )
    hex_counts <- layer_data(p, 1) |>
      transmute(panel = levels(plot_df$panel)[PANEL], x, y, cells = count)
    save_panel(
      p,
      hex_counts,
      run_dir,
      paste0("tnf_vs_ifng_fixed_", pop, "_", s),
      183,
      130
    )
  }
}
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_38c.txt")
)
