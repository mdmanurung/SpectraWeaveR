# Purpose: pseudocolor density of CD3 against TNF on CD4 T cells (CD4 T cells: CD3+ CD4+ CD8- (2D CD4 vs CD8 gate)), one panel per
#          unmixing arm, one figure per donor. Gates come from openCyto via sw_gate_run().
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/09_cd3_tnf_pseudocolor.R <run_dir>
#          from the repo root.
# Inputs:  <run_dir>/manifest.csv [01_manifest-gate.R]
#          the FCS files listed in manifest.csv
# Outputs: <run_dir>/figures_cd3_tnf/cd3_tnf_<sample>.png  pseudocolor panels, one per arm
#          <run_dir>/cd3_tnf_gate_counts.csv               events per gate, arm and sample
#          <run_dir>/sessionInfo_09.txt

# Bioconductor first, tidyverse parts last, so dplyr::filter() is not masked.
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
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
BINS <- 120 # hexagon bins per panel; display only

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
fig_dir <- file.path(run_dir, "figures_cd3_tnf")
stopifnot("09 outputs must not exist" = !dir.exists(fig_dir))
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
  template_file <- file.path(run_dir, paste0("template_cd3_tnf_", arm_name, ".csv"))
  write_csv(template, template_file)
  gs <- sw_gate_run(
    fcs_files = files$file,
    gating_template = template_file,
    transform_channels = c(files$viability[[1]], CD3, CD4, CD8, TNF),
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
write_csv(gate_counts, file.path(run_dir, "cd3_tnf_gate_counts.csv"))

# One figure per donor: pseudocolor hexbin, one facet per arm -----------------------------
# Plain ggplot2: ggcyto fails with the installed ggplot2 (add_ggcyto error), as in 01.
# Channels in the GatingSet are already arcsinh-transformed for CD3 and TNF.
cd4_events <- function(gs, arm, s) {
  gh <- gs[[s]]
  flowCore::exprs(gh_pop_get_data(gh, "root"))[gh_pop_get_indices(gh, "CD4T"), c(CD3, TNF), drop = FALSE] |>
    as_tibble() |>
    mutate(arm = arm)
}
for (s in sort(unique(manifest$sample))) {
  d <- gating_sets |>
    imap(\(gs, arm) cd4_events(gs, arm, s)) |>
    list_rbind() |>
    mutate(arm = factor(arm, levels = arms))
  p <- ggplot(d, aes(.data[[CD3]], .data[[TNF]])) +
    geom_hex(bins = BINS) +
    scale_fill_viridis_c(trans = "sqrt", name = "events") +
    facet_wrap(vars(arm), ncol = 4) +
    labs(
      title = paste(s, "- CD4 T cells (CD4 T cells: CD3+ CD4+ CD8- (2D CD4 vs CD8 gate))"),
      x = "CD3 (arcsinh)",
      y = "TNF (arcsinh)"
    ) +
    theme_bw(base_size = 8)
  ggsave(file.path(fig_dir, paste0("cd3_tnf_", s, ".png")), p, width = 12, height = 9, dpi = 120)
}

# Gating strategy: each gate drawn on its parent, one arm, one figure per donor -------------
# Gate limits come from the openCyto gates themselves; 1D gates are drawn as vertical lines.
STRATEGY_ARM <- "as188_joint" # reference arm (user, 2026-11-09)
gate_limit <- function(gh, node) {
  g <- gh_pop_get_gate(gh, node)
  v <- c(g@min, g@max)
  unname(v[is.finite(v)])
}
strategy_panel <- function(gh, parent, node, x, y = NULL) {
  d <- flowCore::exprs(gh_pop_get_data(gh, parent)) |>
    as_tibble() |>
    slice_sample(n = 20000)
  p <- if (is.null(y)) {
    ggplot(d, aes(.data[[x]])) + geom_histogram(bins = 150)
  } else {
    ggplot(d, aes(.data[[x]], .data[[y]])) + geom_hex(bins = 100) + scale_fill_viridis_c(trans = "sqrt", guide = "none")
  }
  if (is.null(y)) p <- p + geom_vline(xintercept = gate_limit(gh, node), colour = "red")
  p + labs(title = paste0(node, " (on ", parent, ")")) + theme_bw(base_size = 8)
}
set.seed(1) # plot subsample only
for (s in sort(unique(manifest$sample))) {
  gh <- gating_sets[[STRATEGY_ARM]][[s]]
  viab <- manifest |> filter(arm == STRATEGY_ARM) |> pull(viability) |> first()
  panels <- list(
    strategy_panel(gh, "root", "nonDebris", "FSC-A"),
    strategy_panel(gh, "nonDebris", "singlets", "FSC-A", "FSC-H"),
    strategy_panel(gh, "singlets", "live", viab),
    strategy_panel(gh, "live", "CD3pos", CD3),
    strategy_panel(gh, "CD3pos", "CD4T", CD4, CD8),
    strategy_panel(gh, "CD4T", "TNFpos", TNF)
  )
  ggsave(
    file.path(fig_dir, paste0("strategy_", STRATEGY_ARM, "_", s, ".png")),
    cowplot::plot_grid(plotlist = panels, ncol = 3),
    width = 12,
    height = 7,
    dpi = 120
  )
}
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo_09.txt"))
