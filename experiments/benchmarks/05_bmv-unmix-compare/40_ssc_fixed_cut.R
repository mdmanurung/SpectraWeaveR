# Purpose: hexbin of TNF and IFNg (y) against SSC-A (x) on CD4 T and CD8 T cells with the fixed cytokine
#          cut drawn, one panel per arm in version rows, one figure per donor, population and cytokine;
#          plus % positive per arm with donors joined by lines. Same fixed gate as 38_fixed_cut.R, added to
#          the saved GatingSets in memory only.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/40_ssc_fixed_cut.R <run_dir> <cut> [out_root]
#          from the repo root, after 19a_gate_arm.R and 38_fixed_cut.R <run_dir> <cut> [out_root].
#          [out_root] holds figures_fixed_cut_<cut>/ from 38; default <run_dir>.
# Inputs:  <run_dir>/manifest.csv [01_manifest-gate.R]
#          <run_dir>/gs/<arm>/ saved GatingSets with CD4T and CD8T [19a_gate_arm.R]
#          <out_root>/figures_fixed_cut_<cut>/pct_fixed_source.csv  % above the cut [38_fixed_cut.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  labels, VERSION_ROWS, version_grid()
# Outputs: <out_root>/figures_fixed_cut_<cut>/ssc/ssc_vs_<cytokine>_<population>_<sample>.png  hexbins
#          <out_root>/figures_fixed_cut_<cut>/ssc/pct_pos_<population>.png  % positive, donors as lines
#          <out_root>/figures_fixed_cut_<cut>/ssc/sessionInfo_40.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto

# Bioconductor first, tidyverse parts last; dplyr::filter() is called by name.
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS, VERSION_ROWS, version_grid()

# Parameters ------------------------------------------------------------------------------
# Main cut is arcsinh 1 for both cytokines (user, 2026-11-10); the cut comes from the command line.
cli_args <- commandArgs(trailingOnly = TRUE)
CUT <- as.numeric(cli_args[[2]])
CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A") # as in 19a_gate_arm.R
POPULATIONS <- c("CD4T", "CD8T")
BINS <- 120 # hexagon bins; display only
Y_LIMITS <- c(-0.5, 6) # arcsinh display range (user, 2026-11-09: trims are plots only)
ARM_ORDER <- unlist(VERSION_ROWS, use.names = FALSE)
GROUP_BREAKS <- c(1.5, 5.5, 9.5) # dashed lines between the version groups

run_dir <- cli_args[[1]]
out_root <- if (length(cli_args) > 2) cli_args[[3]] else run_dir
cut_dir <- file.path(out_root, paste0("figures_fixed_cut_", CUT))
fig_dir <- file.path(cut_dir, "ssc")
stopifnot("40 outputs must not exist" = !dir.exists(fig_dir))
dir.create(fig_dir)
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
pct <- read_csv(
  file.path(cut_dir, "pct_fixed_source.csv"),
  show_col_types = FALSE
) |>
  mutate(
    arm = factor(arm, levels = ARM_ORDER),
    cytokine = factor(cytokine, levels = names(CHANNELS))
  )

# Load the GatingSets and add the fixed cytokine gates (in memory) -----------------------------
gating_sets <- unique(manifest$arm) |>
  set_names() |>
  map(\(a) {
    gs <- load_gs(file.path(run_dir, "gs", a))
    for (pop in POPULATIONS) {
      iwalk(CHANNELS, \(channel, cy) {
        bounds <- list(c(CUT, Inf)) |> set_names(channel)
        gs_pop_add(
          gs,
          rectangleGate(filterId = paste0(cy, "fix"), .gate = bounds),
          parent = pop
        )
      })
    }
    recompute(gs)
    gs
  })
pct_key <- set_names(
  pct$pct_pos,
  paste(pct$arm, pct$donor, pct$parent, pct$cytokine)
)

# SSC-A vs cytokine hexbins, fixed gate drawn ---------------------------------------------------
for (s in sort(unique(manifest$sample))) {
  for (pop in POPULATIONS) {
    for (cy in names(CHANNELS)) {
      panels <- gating_sets |>
        imap(\(gs, arm) {
          ggcyto(
            gs[s],
            aes(x = `SSC-A`, y = !!rlang::sym(CHANNELS[[cy]])),
            subset = pop
          ) +
            geom_hex(bins = BINS) +
            geom_hline(yintercept = CUT, colour = "red") +
            ggcyto_par_set(limits = list(y = Y_LIMITS)) +
            labs(
              title = ARM_LABELS[[arm]],
              subtitle = sprintf(
                "%s+ %.1f%%",
                cy,
                pct_key[[paste(arm, s, pop, cy)]]
              ),
              x = "SSC-A",
              y = paste0(cy, " (arcsinh)")
            ) +
            PANEL_THEME
        })
      ggsave(
        file.path(fig_dir, paste0("ssc_vs_", cy, "_", pop, "_", s, ".png")),
        version_grid(panels),
        width = 14,
        height = 10,
        dpi = 120,
        bg = "white"
      )
    }
  }
}

# % positive per arm, one figure per population, donors joined by lines ------------------------
for (pop in POPULATIONS) {
  p <- pct |>
    dplyr::filter(parent == pop) |>
    ggplot(aes(arm, pct_pos)) +
    annotate(
      "rect",
      xmin = match(ZIP_ARMS, ARM_ORDER) - 0.5,
      xmax = match(ZIP_ARMS, ARM_ORDER) + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = ZIP_FILL
    ) + # zip arms shaded
    geom_vline(
      xintercept = GROUP_BREAKS,
      linetype = "dashed",
      colour = "grey70"
    ) +
    geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey40") +
    geom_line(
      aes(group = donor, colour = donor),
      alpha = 0.6,
      linewidth = 0.4
    ) +
    geom_point(aes(colour = donor), size = 1.8) +
    scale_colour_brewer(palette = "Dark2", name = "Donor") +
    scale_x_discrete(labels = ARM_LABELS) +
    facet_wrap(~cytokine, scales = "free_y") +
    labs(
      x = NULL,
      y = paste0("% positive of ", pop),
      title = paste0(pop, ", cut at arcsinh ", CUT)
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(
    file.path(fig_dir, paste0("pct_pos_", pop, ".png")),
    p,
    width = 12,
    height = 6.5,
    dpi = 130,
    bg = "white"
  )
}
writeLines(
  capture.output(sessionInfo()),
  file.path(fig_dir, "sessionInfo_40.txt")
)
