# Purpose: pseudocolor density of CD3 against TNF on CD4 T cells (CD4 T cells: CD3+ CD4+ CD8- (2D CD4 vs CD8 gate)), one panel per
#          unmixing arm, one figure per donor. Gates come from openCyto via sw_gate_run().
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/13_tnf_vs_markers.R <run_dir>
#          from the repo root.
# Inputs:  <run_dir>/manifest.csv [01_manifest-gate.R]
#          the FCS files listed in manifest.csv
# Outputs: <run_dir>/figures_tnf_vs_ggcyto/tnf_vs_<marker>_<sample>.png  pseudocolor panels, one per arm
#          <run_dir>/sessionInfo_13.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto

# Bioconductor first, tidyverse parts last, so dplyr::filter() is not masked.
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(ggcyto)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
DEBRIS_RANGE <- "c(2e5, 1.5e6)" # same debris window as 01 (user, 2026-11-09)
CYTOKINE_RANGE <- "c(1, 4)" # mindensity window for the TNF cut, same as 03 (user, grill Q1, 2026-11-09)
LINEAGE_RANGE <- "c(0, 3)" # mindensity window for CD3, CD4, CD8 cuts, same as 03 (user, 2026-11-09)
CD3 <- "SBUV510-A" # channels as in 03_downstream.R
CD4 <- "Spark Blue 574-A"
CD8 <- "SBUV445-A"
TNF <- "BV510-A"
IFNG <- "BV570-A" # IFNg channel as in 03_downstream.R
BINS <- 120 # hexagon bins per panel; display only

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
fig_dir <- file.path(run_dir, "figures_tnf_vs_ggcyto")
stopifnot("13 outputs must not exist" = !dir.exists(fig_dir))
dir.create(fig_dir)
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)

# Gate one arm: debris, singlets, live, CD3+, CD4+, CD8- -----------------------------------
gate_arm <- function(arm_name) {
  files <- manifest |>
    filter(arm == arm_name) |>
    arrange(sample)
  lineage_args <- paste0("gate_range=", LINEAGE_RANGE)
  cytokine_args <- paste0("gate_range=", CYTOKINE_RANGE)
  template <- tibble(
    alias = c("nonDebris", "singlets", "live", "CD3pos", "CD4T", "TNFpos"),
    pop = c("+", "+", "-", "+", "+-", "+"),
    parent = c("root", "nonDebris", "singlets", "live", "CD3pos", "CD4T"),
    dims = c("FSC-A", "FSC-A,FSC-H", files$viability[[1]], CD3, paste(CD4, CD8, sep = ","), TNF),
    gating_method = "mindensity",
    gating_args = c(paste0("gate_range=", DEBRIS_RANGE), "wider_gate=TRUE", NA, lineage_args, lineage_args, cytokine_args),
    collapseDataForGating = FALSE,
    groupBy = NA,
    preprocessing_method = NA,
    preprocessing_args = NA
  )
  template$gating_method[[2]] <- "singletGate"
  template_file <- file.path(run_dir, paste0("template_tnf_vs_", arm_name, ".csv"))
  write_csv(template, template_file)
  gs <- sw_gate_run(
    fcs_files = files$file,
    gating_template = template_file,
    transform_channels = c(files$viability[[1]], CD3, CD4, CD8, TNF, IFNG),
    cofactor = COFACTOR
  )
  sampleNames(gs) <- files$sample
  gs
}

arms <- unique(manifest$arm)
gating_sets <- arms |>
  set_names() |>
  map(gate_arm)

gate_counts <- gating_sets |>
  imap(\(gs, arm) gs_pop_get_count_fast(gs) |> as_tibble() |> mutate(arm = arm)) |>
  list_rbind()


# One figure per marker and donor: ggcyto hexbin of TNF (y) against the marker (x) on CD4 T cells ---
# ggcyto needs ggplot2 3.5.x (private library above). CD3, CD4, CD8, TNF, IFNg are arcsinh in the
# GatingSet; FSC-A and SSC-A stay linear.
X_MARKERS <- c(`FSC-A` = "FSC-A", `SSC-A` = "SSC-A", CD8 = CD8, CD4 = CD4, IFNg = IFNG) # display choice (user: "FSC or CD8", "IFNg as well")
for (s in sort(unique(manifest$sample))) {
  for (m in names(X_MARKERS)) {
    panels <- gating_sets |>
      imap(\(gs, arm) {
        ggcyto(gs[s], aes(x = !!rlang::sym(X_MARKERS[[m]]), y = !!rlang::sym(TNF)), subset = "CD4T") +
          geom_hex(bins = BINS) +
          ggcyto_par_set(limits = "data") +
          labs(title = arm, x = m, y = "TNF") +
          theme(strip.text = element_blank(), legend.position = "none")
      })
    ggsave(
      file.path(fig_dir, paste0("tnf_vs_", m, "_", s, ".png")),
      cowplot::plot_grid(plotlist = map(panels, as.ggplot), ncol = 4),
      width = 12,
      height = 9,
      dpi = 120
    )
  }
}
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo_13.txt"))
