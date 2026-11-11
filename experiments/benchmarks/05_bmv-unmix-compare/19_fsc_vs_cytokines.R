# Purpose: hexbin of TNF and IFNg (y) against FSC-A or SSC-A (x) on CD4 T and CD8 T cells, one panel per
#          unmixing arm, one figure per donor, with the cytokine gate and its % positive. Loads the GatingSets
#          saved by 19a_gate_arm.R, so every figure uses the same gates.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/19_fsc_vs_cytokines.R <run_dir> [scatter] [out_root]
#            scatter = FSC-A (default) or SSC-A, the x axis of the cytokine panels
#            out_root = where outputs go; default <run_dir>
#          from the repo root, after 19a_gate_arm.R has run for every arm in <run_dir>/manifest.csv.
# Inputs:  <run_dir>/manifest.csv [01_manifest-gate.R]
#          <run_dir>/gs/<arm>/ saved GatingSets [19a_gate_arm.R]
# Outputs: <out_root>/figures_<fsc|ssc>_cyto/<fsc|ssc>_vs_<cytokine>_<population>_<sample>.png  hexbin panels
#          <out_root>/<fsc|ssc>_cyto_gate_counts.csv  events per gate, arm and sample
#          <out_root>/figures_<fsc|ssc>_cyto/gate_cd4_cd8_<sample>.png  CD4 vs CD8 gate check per arm
#          <out_root>/figures_<fsc|ssc>_cyto/sessionInfo_19.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto

# Bioconductor first, tidyverse parts last; dplyr::filter() is called by name (as in 19a).
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS: arm ID -> plot title

# Parameters ------------------------------------------------------------------------------
CD4 <- "Spark Blue 574-A" # channels as in 19a_gate_arm.R
CD8 <- "SBUV445-A"
CYTOKINES <- c(TNF = "BV510-A", IFNg = "BV570-A")
POPULATIONS <- c("CD4T", "CD8T")
BINS <- 120 # hexagon bins per panel; display only
cli_args <- commandArgs(trailingOnly = TRUE)
run_dir <- cli_args[[1]]
scatter <- if (length(cli_args) > 1) cli_args[[2]] else "FSC-A"
out_root <- if (length(cli_args) > 2) cli_args[[3]] else run_dir
stopifnot("scatter is FSC-A or SSC-A" = scatter %in% c("FSC-A", "SSC-A"))
prefix <- tolower(sub("-A$", "", scatter)) # fsc or ssc, for file names
fig_dir <- file.path(out_root, paste0("figures_", prefix, "_cyto"))
stopifnot("19 outputs must not exist" = !dir.exists(fig_dir))
dir.create(fig_dir, recursive = TRUE)
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)

# Load the saved GatingSets ------------------------------------------------------------------
arms <- unique(manifest$arm)
gating_sets <- arms |>
  set_names() |>
  map(\(a) load_gs(file.path(run_dir, "gs", a)))

gate_counts <- gating_sets |>
  imap(\(gs, arm) {
    gs_pop_get_count_fast(gs) |> as_tibble() |> mutate(arm = arm)
  }) |>
  list_rbind()
write_csv(
  gate_counts,
  file.path(out_root, paste0(prefix, "_cyto_gate_counts.csv"))
) # before plots, so a plot failure keeps the counts
# % positive per arm, donor, population and cytokine, for the panel subtitles
pct_pos <- gate_counts |>
  dplyr::filter(str_detect(Population, "/(CD4T|CD8T)/(TNF|IFNg)pos$")) |>
  mutate(
    key = paste(
      arm,
      name,
      str_extract(Population, "CD[48]T"),
      str_extract(Population, "(TNF|IFNg)(?=pos$)")
    )
  ) |>
  (\(d) set_names(100 * d$Count / d$ParentCount, d$key))()


# Version grid: one labelled row per code version, panels named by arm --------------------------
stopifnot(
  "every arm sits in one version row" = setequal(unlist(VERSION_ROWS), arms),
  "every arm has a label" = all(arms %in% names(ARM_LABELS))
)
# One figure per donor, population and cytokine: ggcyto hexbin of the cytokine (y) against FSC-A (x) ---
# ggcyto needs ggplot2 3.5.x (private library above). Cytokines are arcsinh in the GatingSet; FSC-A stays
# linear. The y-scale lower limit drops events below DISPLAY_FLOOR from the plot only.
DISPLAY_FLOOR <- -0.5 # display trim on the cytokine axis, arcsinh (user, 2026-11-09: plots only)
DISPLAY_CEIL <- 6 # upper axis limit, arcsinh; above the brightest cytokine event, display only
for (s in sort(unique(manifest$sample))) {
  for (pop in POPULATIONS) {
    for (cy in names(CYTOKINES)) {
      panels <- gating_sets |>
        imap(\(gs, arm) {
          ggcyto(
            gs[s],
            aes(x = !!rlang::sym(scatter), y = !!rlang::sym(CYTOKINES[[cy]])),
            subset = pop
          ) +
            geom_hex(bins = BINS) +
            geom_gate(paste0(pop, "/", cy, "pos")) +
            scale_y_continuous(limits = c(DISPLAY_FLOOR, DISPLAY_CEIL)) +
            ggcyto_par_set(limits = list(y = c(DISPLAY_FLOOR, DISPLAY_CEIL))) +
            labs(
              title = ARM_LABELS[[arm]],
              subtitle = sprintf(
                "%s+ %.1f%%",
                cy,
                pct_pos[[paste(arm, s, pop, cy)]]
              ),
              x = scatter,
              y = paste0(cy, " (arcsinh)")
            ) +
            PANEL_THEME
        })
      ggsave(
        file.path(
          fig_dir,
          paste0(prefix, "_vs_", cy, "_", pop, "_", s, ".png")
        ),
        version_grid(panels),
        width = 14,
        height = 9,
        dpi = 120,
        bg = "white"
      )
    }
  }
}
# Gate counts and gate-check pages ----------------------------------------------------------
# Counts per gate, arm and sample, so a failed gate shows as a number. Gate-check pages: the
# CD4 vs CD8 quadrants per arm, with the % of parent.
for (s in sort(unique(manifest$sample))) {
  panels <- gating_sets |>
    imap(\(gs, arm) {
      ggcyto(
        gs[s],
        aes(x = !!rlang::sym(CD4), y = !!rlang::sym(CD8)),
        subset = "CD3pos"
      ) +
        geom_hex(bins = BINS) +
        geom_gate(c("CD4T", "CD8T")) +
        geom_stats() +
        labs(
          title = ARM_LABELS[[arm]],
          x = "CD4 (arcsinh)",
          y = "CD8 (arcsinh)"
        ) +
        PANEL_THEME
    })
  ggsave(
    file.path(fig_dir, paste0("gate_cd4_cd8_", s, ".png")),
    version_grid(panels),
    width = 14,
    height = 9,
    dpi = 120,
    bg = "white"
  )
}
writeLines(
  capture.output(sessionInfo()),
  file.path(fig_dir, "sessionInfo_19.txt")
)
