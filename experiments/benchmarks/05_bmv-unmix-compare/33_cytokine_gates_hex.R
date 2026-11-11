# Purpose: show the cytokine gates with ggcyto: TNF (x) against IFNg (y) hexbin on CD4 T and CD8 T cells,
#          with each arm's own TNF+ and IFNg+ gates drawn by geom_gate, one panel per arm in version rows,
#          one page per donor and population. Loads the GatingSets saved by 19a_gate_arm.R.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/33_cytokine_gates_hex.R <run_dir> [out_root]
#          from the repo root, after 19a_gate_arm.R has run for every arm in <run_dir>/manifest.csv.
#          [out_root] is where the figure folder goes; default <run_dir>.
# Inputs:  <run_dir>/manifest.csv [01_manifest-gate.R]
#          <run_dir>/gs/<arm>/ saved GatingSets [19a_gate_arm.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  labels, VERSION_ROWS, version_grid()
#          experiments/bmv_test/compare_3D/data/lib_ggplot2_35: ggplot2 3.5.2 (ggcyto 1.34 fails on 4.0.3)
# Outputs: <out_root>/figures_cytokine_gates/tnf_vs_ifng_<population>_<sample>.png  hexbin with both gates
#          <out_root>/figures_cytokine_gates/sessionInfo_33.txt

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
TNF <- "BV510-A" # channels as in 19a_gate_arm.R
IFNG <- "BV570-A"
POPULATIONS <- c("CD4T", "CD8T")
BINS <- 100 # hexagon bins per panel; display only
AXIS_LIMITS <- c(-0.5, 5) # arcsinh display range on both axes (user, 2026-11-09: trims are plots only)

cli_args <- commandArgs(trailingOnly = TRUE)
run_dir <- cli_args[[1]]
out_root <- if (length(cli_args) > 1) cli_args[[2]] else run_dir
fig_dir <- file.path(out_root, "figures_cytokine_gates")
stopifnot("33 outputs must not exist" = !dir.exists(fig_dir))
dir.create(fig_dir, recursive = TRUE)
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)

# Load the saved GatingSets ------------------------------------------------------------------
arms <- unique(manifest$arm)
stopifnot(
  "every arm sits in one version row" = setequal(unlist(VERSION_ROWS), arms),
  "every arm has a label" = all(arms %in% names(ARM_LABELS))
)
gating_sets <- arms |>
  set_names() |>
  map(\(a) load_gs(file.path(run_dir, "gs", a)))

# % positive per arm, donor, population and cytokine, for the subtitles ------------------------
pct_of <- function(gs, sample, node) {
  stats <- gs_pop_get_stats(gs[sample], node, type = "percent")
  100 * stats$percent
}

# One page per donor and population: TNF vs IFNg hexbin with both cytokine gates ---------------------
for (s in sort(unique(manifest$sample))) {
  for (pop in POPULATIONS) {
    panels <- gating_sets |>
      imap(\(gs, arm) {
        tnf_node <- paste0(pop, "/TNFpos")
        ifng_node <- paste0(pop, "/IFNgpos")
        ggcyto(
          gs[s],
          aes(x = !!rlang::sym(TNF), y = !!rlang::sym(IFNG)),
          subset = pop
        ) +
          geom_hex(bins = BINS) +
          geom_gate(c(tnf_node, ifng_node)) +
          ggcyto_par_set(limits = list(x = AXIS_LIMITS, y = AXIS_LIMITS)) +
          labs(
            title = ARM_LABELS[[arm]],
            subtitle = sprintf(
              "TNF+ %.1f%%  IFNg+ %.1f%%",
              pct_of(gs, s, tnf_node),
              pct_of(gs, s, ifng_node)
            ),
            x = "TNF (arcsinh)",
            y = "IFNg (arcsinh)"
          ) +
          PANEL_THEME
      })
    ggsave(
      file.path(fig_dir, paste0("tnf_vs_ifng_", pop, "_", s, ".png")),
      version_grid(panels),
      width = 14,
      height = 10,
      dpi = 120,
      bg = "white"
    )
  }
}
writeLines(
  capture.output(sessionInfo()),
  file.path(fig_dir, "sessionInfo_33.txt")
)
