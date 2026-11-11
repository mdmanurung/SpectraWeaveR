# Purpose: draw the openCyto gating strategy (debris, singlets, live, CD3+, CD4 T and CD8 T from the pooled
#          t-mixture CD4 vs CD8 gate, TNF+ and IFNg+ on CD4 T) with ggcyto, each gate on its parent with the % of
#          parent, plus the gating tree. Loads the GatingSets saved by 19a_gate_arm.R (same gates as 19).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/10_ggcyto_strategy.R <run_dir> [out_root]
#          from the repo root, after 19a_gate_arm.R has run for the STRATEGY_ARMS.
#          [out_root] is where figures_ggcyto/ goes; default <run_dir>.
# Inputs:  <run_dir>/manifest.csv [01_manifest-gate.R]
#          <run_dir>/gs/<arm>/ saved GatingSets [19a_gate_arm.R]
#          experiments/bmv_test/compare_3D/data/lib_ggplot2_35: ggplot2 3.5.2, installed there because
#          ggcyto 1.34 fails with the project's ggplot2 4.0.3 (add_ggcyto error)
# Outputs: <out_root>/figures_ggcyto/strategy_<arm>_<sample>.png  one page per arm and donor
#          <out_root>/figures_ggcyto/tree_<arm>.png  gating hierarchy
#          <out_root>/figures_ggcyto/sessionInfo_10.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths()))

# Bioconductor first, tidyverse parts last. dplyr is attached without filter(): openCyto's quadrant gates
# calls filter() unqualified and must find flowCore::filter(), so dplyr::filter() is called by name.
# rlang is not attached:
# its exprs() would mask flowCore::exprs(), which openCyto calls unqualified.
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
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS: arm ID -> plot title
stopifnot(
  "ggplot2 3.5.x is needed for ggcyto" = packageVersion("ggplot2") < "4.0.0"
)

# Parameters ------------------------------------------------------------------------------
CD3 <- "SBUV510-A" # channels as in 19a_gate_arm.R
CD4 <- "Spark Blue 574-A"
CD8 <- "SBUV445-A"
TNF <- "BV510-A"
IFNG <- "BV570-A"
STRATEGY_ARMS <- c("spectroflo", "as102_var", "as188_joint") # one arm per method family; edit to add more
BINS <- 100 # hexagon bins for 2D panels; display only

cli_args <- commandArgs(trailingOnly = TRUE)
run_dir <- cli_args[[1]]
out_root <- if (length(cli_args) > 1) cli_args[[2]] else run_dir
fig_dir <- file.path(out_root, "figures_ggcyto")
stopifnot("10 outputs must not exist" = !dir.exists(fig_dir))
dir.create(fig_dir, recursive = TRUE)
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)

# Load one arm's saved GatingSet and its viability channel ---------------------------------------
load_arm <- function(arm_name) {
  list(
    gs = load_gs(file.path(run_dir, "gs", arm_name)),
    viability = manifest |>
      dplyr::filter(arm == arm_name) |>
      pull(viability) |>
      unique()
  )
}

# One ggcyto panel: the gate drawn on its parent, with % of parent -------------------------
# 1D gates get a density, 2D gates a hexbin pseudocolor; geom_gate() draws the openCyto gate.
panel <- function(gs, sample, parent, node, x, y = NULL) {
  base <- if (is.null(y)) {
    ggcyto(gs[sample], aes(x = !!rlang::sym(x)), subset = parent) +
      geom_histogram(bins = 150)
  } else {
    ggcyto(
      gs[sample],
      aes(x = !!rlang::sym(x), y = !!rlang::sym(y)),
      subset = parent
    ) +
      geom_hex(bins = BINS) +
      scale_fill_viridis_c(trans = "sqrt", guide = "none")
  }
  base +
    geom_gate(node) +
    geom_stats(adjust = 0.8, size = 3.5) +
    ggtitle(paste(basename(node), collapse = " + ")) +
    theme(
      strip.text = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
}

# One page per arm and donor ----------------------------------------------------------------
for (arm_name in STRATEGY_ARMS) {
  g <- load_arm(arm_name)
  # gating tree (same for every sample of this arm)
  png(
    file.path(fig_dir, paste0("tree_", arm_name, ".png")),
    width = 900,
    height = 700
  )
  plot(g$gs[[1]])
  dev.off()
  for (s in sampleNames(g$gs)) {
    panels <- list(
      panel(g$gs, s, "root", "nonDebris", "FSC-A"),
      panel(g$gs, s, "nonDebris", "singlets", "FSC-A", "FSC-H"),
      panel(g$gs, s, "singlets", "live", g$viability),
      panel(g$gs, s, "live", "CD3pos", CD3),
      panel(g$gs, s, "CD3pos", c("CD4T", "CD8T"), CD4, CD8)
    ) # cytokine gates are the fixed cut, drawn by 38_fixed_cut.R (user, 2026-11-10)
    ggsave(
      file.path(fig_dir, paste0("strategy_", arm_name, "_", s, ".png")),
      cowplot::plot_grid(
        cowplot::ggdraw() +
          cowplot::draw_label(
            paste0(ARM_LABELS[[arm_name]], ", donor ", s),
            fontface = "bold",
            size = 14
          ),
        cowplot::plot_grid(plotlist = map(panels, as.ggplot), ncol = 3),
        ncol = 1,
        rel_heights = c(0.06, 1)
      ),
      width = 16,
      height = 8,
      dpi = 120,
      bg = "white"
    )
  }
}
writeLines(
  capture.output(sessionInfo()),
  file.path(fig_dir, "sessionInfo_10.txt")
)
