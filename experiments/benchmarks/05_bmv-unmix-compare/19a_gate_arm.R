# Purpose: gate one unmixing arm with openCyto and save its GatingSet, so the figure scripts (19, 10) load
#          the same gates without re-gating. One process per arm; run arms in parallel from the shell.
#          Gates: debris, singlets, live, CD3+, CD4 T / CD8 T (t-mixture, pooled over donors), TNF+ and IFNg+
#          under CD4 T and CD8 T.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/19a_gate_arm.R <run_dir> <arm>
#          from the repo root.
# Inputs:  <run_dir>/manifest.csv [01_manifest-gate.R]
#          the FCS files listed in manifest.csv
# Outputs: <run_dir>/gs/<arm>/              saved GatingSet (flowWorkspace::save_gs)
#          <run_dir>/gate_counts/<arm>.csv  events per gate and sample
#          <run_dir>/templates/<arm>.csv    openCyto template used

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_cytoutils", .libPaths())) # cytoUtils 0.1.1 (RGLab GitHub) for tailgate

# Bioconductor first, tidyverse parts last. dplyr is attached without filter(): openCyto's quadrant gates
# call filter() unqualified and must find flowCore::filter(), so dplyr::filter() is called by name.
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(cytoUtils)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
register_plugins(fun = cytoUtils:::.tailgate, methodName = "tailgate") # not registered by openCyto 2.18

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
DEBRIS_RANGE <- "c(2e5, 1.5e6)" # same debris window as 01 (user, 2026-11-09)
LINEAGE_RANGE <- "c(0, 3)" # mindensity window for the CD3 cut, same as 03 (user, 2026-11-09)
CD3 <- "SBUV510-A" # channels as in 03_downstream.R
CD4 <- "Spark Blue 574-A"
CD8 <- "SBUV445-A"
TNF <- "BV510-A"
IFNG <- "BV570-A"
# CD4 T and CD8 T from one 2D t-mixture quadrant gate (openCyto quadGate.tmix) on CD4 vs CD8 within CD3+,
# fitted once per arm on all donors pooled (user, 2026-11-10): K = 3 as in the openCyto vignette, ellipse
# quantile 0.8 (openCyto default). Per-donor tmix misplaced 3 of 55 cuts; mindensity and quadGate.seq pinned.
CD4_CD8_TMIX_ARGS <- "K=3, usePrior='no', quantile1=0.8, quantile3=0.8"
QUADRANTS <- c(
  CD4T = paste0(CD4, "+", CD8, "-"),
  CD8T = paste0(CD4, "-", CD8, "+"),
  DP = paste0(CD4, "+", CD8, "+"),
  DN = paste0(CD4, "-", CD8, "-")
) # openCyto's quadrant node names, renamed to these aliases after gating
T_PARENT <- "/nonDebris/singlets/live/CD3pos" # CD3+ node path, parent of the quadrant nodes
# Cytokine+ gates under CD4T and CD8T (user, 2026-11-10). No unstimulated control exists, so cuts are
# data-driven per arm: one gate fitted on all 5 donors pooled, applied to each donor (user, 2026-11-10). TNF: mindensity in the TNF window approved for 03 (grill Q1, 2026-11-09).
# IFNg: cytoUtils tailgate, tol 0.01 (cytoUtils default), 2 peaks with the lower one as reference.
CYTOKINE_GATES <- list(
  TNF = c(method = "mindensity", args = "gate_range=c(1, 4)"),
  IFNg = c(method = "tailgate", args = "tol=0.01, num_peaks=2, ref_peak=1")
)
CYTOKINES <- c(TNF = TNF, IFNg = IFNG)
POPULATIONS <- c("CD4T", "CD8T")

cli_args <- commandArgs(trailingOnly = TRUE)
run_dir <- cli_args[[1]]
arm_name <- cli_args[[2]]
gs_dir <- file.path(run_dir, "gs", arm_name)
stopifnot("GatingSet for this arm must not exist" = !dir.exists(gs_dir))
walk(file.path(run_dir, c("gs", "gate_counts", "templates")), \(d) {
  dir.create(d, showWarnings = FALSE)
})
files <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE) |>
  dplyr::filter(arm == arm_name) |>
  arrange(sample)
stopifnot("arm has 5 files" = nrow(files) == 5)

# Template up to the CD4/CD8 quadrants --------------------------------------------------------
# quadGate.tmix with pop "+/-+/-" makes four quadrant nodes named from the channels (openCyto ignores the
# alias). Pooling: collapseDataForGating with groupBy = number of donors puts all donors in one fit.
template <- tibble(
  alias = c("nonDebris", "singlets", "live", "CD3pos", "*"),
  pop = c("+", "+", "-", "+", "+/-+/-"),
  parent = c("root", "nonDebris", "singlets", "live", "CD3pos"),
  dims = c(
    "FSC-A",
    "FSC-A,FSC-H",
    files$viability[[1]],
    CD3,
    paste(CD4, CD8, sep = ",")
  ),
  gating_method = c(
    "mindensity",
    "singletGate",
    "mindensity",
    "mindensity",
    "quadGate.tmix"
  ),
  gating_args = c(
    paste0("gate_range=", DEBRIS_RANGE),
    "wider_gate=TRUE",
    NA,
    paste0("gate_range=", LINEAGE_RANGE),
    CD4_CD8_TMIX_ARGS
  ),
  collapseDataForGating = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  groupBy = c(NA, NA, NA, NA, nrow(files)),
  preprocessing_method = NA,
  preprocessing_args = NA
)
template_file <- file.path(run_dir, "templates", paste0(arm_name, ".csv"))
write_csv(template, template_file)

# Gate, rename the quadrants, add the cytokine gates -----------------------------------------
gs <- sw_gate_run(
  fcs_files = files$file,
  gating_template = template_file,
  transform_channels = c(files$viability[[1]], CD3, CD4, CD8, TNF, IFNG),
  cofactor = COFACTOR
)
sampleNames(gs) <- files$sample
iwalk(QUADRANTS, \(node, alias) {
  gs_pop_set_name(gs, paste0(T_PARENT, "/", node), alias)
})
# cytokine+ gates under each T cell population, named <cytokine>pos (TNFpos, IFNgpos); the template parser
# cannot name a quadrant node as a parent, so these are added after the rename
for (pop in POPULATIONS) {
  iwalk(CYTOKINES, \(channel, cy) {
    gs_add_gating_method(
      gs,
      alias = paste0(cy, "pos"),
      pop = "+",
      parent = pop,
      dims = channel,
      gating_method = CYTOKINE_GATES[[cy]][["method"]],
      gating_args = CYTOKINE_GATES[[cy]][["args"]],
      collapseDataForGating = TRUE,
      groupBy = nrow(files)
    )
  })
}

# Save -------------------------------------------------------------------------------------
gs_pop_get_count_fast(gs) |>
  as_tibble() |>
  mutate(arm = arm_name) |>
  write_csv(file.path(run_dir, "gate_counts", paste0(arm_name, ".csv")))
save_gs(gs, gs_dir)
