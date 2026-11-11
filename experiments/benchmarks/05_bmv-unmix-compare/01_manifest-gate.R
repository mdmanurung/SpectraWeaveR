# Purpose: gate nonDebris, singlets and live cells with openCyto in each unmixing arm, then
#          keep the events that are live in every arm (shared event set).
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/01_manifest-gate.R \
#            experiments/runs/2026-11-09_bmv-unmix-compare
# Inputs:  experiments/bmv_test/compare_3D/data/ (unzipped compare_3Dstim.zip) [user]
# Outputs: <run_dir>/manifest.csv           arm x sample x file
#          <run_dir>/gate_counts.csv        events per node, arm and sample
#          <run_dir>/shared_events.rds      per sample: logical keep vector (live in all arms)
#          <run_dir>/figures/gates_<sample>.png   scatter, singlet and viability gates per arm
#          <run_dir>/sessionInfo.txt

# tidyverse and here are not in the pixi env: load the parts used; run from the repo root
# Bioconductor first, tidyverse parts last, so dplyr::filter() is not masked.
suppressPackageStartupMessages({
  library(flowCore) # attached first: openCyto gate code calls exprs() unqualified
  library(flowWorkspace)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor for spectral flow (CLAUDE.md, confirmed by user)
DEBRIS_RANGE <- "c(2e5, 1.5e6)" # P1: FSC-A window searched for the debris/cell valley
DATA_DIR <- "experiments/bmv_test/compare_3D/data"
ARMS <- tribble(
  ~arm             , ~dir                                              , ~viability        ,
  "spectroflo"     , "manual"                                          , "LIVE DEAD NIR-A" ,
  "as102_basic"    , "autoSpectral_v1.0.2/AutospectralUnmix"           , "Zombie NIR-A"    ,
  "as102_var"      , "autoSpectral_v1.0.2/AutospectralUnmix_optimized" , "Zombie NIR-A"    ,
  "as188_joint"    , "autoSpectral_v1.8.8/AutoSpectral_unmixed"        , "Zombie NIR-A"    ,
  "s102_c188"      , "unmix_cells/s102_c188"                           , "Zombie NIR-A"    ,
  "s188_c102"      , "unmix_cells/s188_143803_c102"                    , "Zombie NIR-A"    ,
  "as188_fast_csv" , "unmix_cells/repro_s188_c188_fast"                , "Zombie NIR-A"    ,
  "best102_neg"    , "unmix_cells/best102_neg"                         , "Zombie NIR-A"    ,
  "best102_af"     , "unmix_cells/best102_af"                          , "Zombie NIR-A"    ,
  "best188_p1"     , "unmix_cells/best188_p1"                          , "Zombie NIR-A"    ,
  "best188_p5"     , "unmix_cells/best188_p5"                          , "Zombie NIR-A"
)

# Run dir ---------------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run dir must not exist" = !dir.exists(run_dir))
dir.create(file.path(run_dir, "figures"), recursive = TRUE)

# Manifest --------------------------------------------------------------------------------
manifest <- ARMS |>
  reframe(
    file = list.files(
      file.path(DATA_DIR, dir),
      pattern = "\\.fcs$",
      full.names = TRUE
    ),
    .by = c(arm, dir, viability)
  ) |>
  mutate(sample = str_extract(basename(file), "^[A-Z][0-9]+"))
write_csv(manifest, file.path(run_dir, "manifest.csv"))

# Gate one arm ----------------------------------------------------------------------------
# One template per arm because the viability channel is named differently in SpectroFlo.
gate_arm <- function(arm_row) {
  files <- manifest |>
    filter(arm == arm_row$arm) |>
    arrange(sample)
  template <- tibble(
    alias = c("nonDebris", "singlets", "live"),
    pop = c("+", "+", "-"),
    parent = c("root", "nonDebris", "singlets"),
    dims = c("FSC-A", "FSC-A,FSC-H", arm_row$viability),
    gating_method = c("mindensity", "singletGate", "mindensity"),
    gating_args = c(paste0("gate_range=", DEBRIS_RANGE), "wider_gate=TRUE", NA),
    collapseDataForGating = FALSE,
    groupBy = NA,
    preprocessing_method = NA,
    preprocessing_args = NA
  )
  template_file <- file.path(run_dir, paste0("template_", arm_row$arm, ".csv"))
  write_csv(template, template_file)
  gs <- sw_gate_run(
    fcs_files = files$file,
    gating_template = template_file,
    transform_channels = arm_row$viability,
    cofactor = COFACTOR
  )
  sampleNames(gs) <- files$sample
  gs
}

gating_sets <- ARMS |>
  split(~arm) |>
  map(gate_arm)

# Gate counts and shared live set ---------------------------------------------------------
gate_counts <- gating_sets |>
  imap(\(gs, arm) {
    gs_pop_get_count_fast(gs) |> as_tibble() |> mutate(arm = arm)
  }) |>
  list_rbind()
write_csv(gate_counts, file.path(run_dir, "gate_counts.csv"))

# Per sample: a logical vector over all events, TRUE when live in every arm.
live_mask <- function(gs, sample) {
  idx <- gh_pop_get_indices(gs[[sample]], "live")
  idx
}
samples <- sort(unique(manifest$sample))
shared <- samples |>
  set_names() |>
  map(\(s) gating_sets |> map(live_mask, sample = s) |> reduce(`&`))
saveRDS(shared, file.path(run_dir, "shared_events.rds"))

shared_summary <- tibble(
  sample = samples,
  n_shared = map_int(shared, sum),
  n_total = map_int(shared, length)
)
write_csv(shared_summary, file.path(run_dir, "shared_live_counts.csv"))

# Gate plots (inspect before any metric) --------------------------------------------------
# Plain ggplot2: ggcyto does not work with the installed ggplot2. Per arm: FSC-A with the
# debris cut, FSC-A vs FSC-H coloured by singlet call, viability with the live cut.
PLOT_EVENTS <- 20000 # events drawn per panel
set.seed(1) # plot subsample only; no metric uses it
gate_limit <- function(gh, node) {
  # a "-" population is stored as the "+" gate negated, so take whichever bound is finite
  g <- gh_pop_get_gate(gh, node)
  v <- c(g@min, g@max)
  unname(v[is.finite(v)])
}
plot_arm <- function(gs, arm, s, viability) {
  gh <- gs[[s]]
  d <- flowCore::exprs(gh_pop_get_data(gh, "root")) |>
    as_tibble() |>
    mutate(
      in_debris = gh_pop_get_indices(gh, "nonDebris"),
      in_singlet = gh_pop_get_indices(gh, "singlets"),
      in_live = gh_pop_get_indices(gh, "live")
    ) |>
    slice_sample(n = PLOT_EVENTS) # every sample has more events than this
  p_debris <- ggplot(d, aes(.data[["FSC-A"]])) +
    geom_histogram(bins = 200) +
    geom_vline(
      xintercept = gate_limit(gh, "nonDebris"),
      colour = "red"
    ) +
    labs(title = paste(arm, s, "nonDebris"), x = "FSC-A")
  p_singlet <- d |>
    filter(in_debris) |>
    ggplot(aes(.data[["FSC-A"]], .data[["FSC-H"]], colour = in_singlet)) +
    geom_point(size = 0.1, alpha = 0.3) +
    labs(title = paste(arm, s, "singlets"))
  p_live <- d |>
    filter(in_singlet) |>
    ggplot(aes(.data[[viability]])) +
    geom_histogram(bins = 200) +
    geom_vline(xintercept = gate_limit(gh, "live"), colour = "red") +
    labs(title = paste(arm, s, "live"), x = paste("arcsinh", viability))
  list(p_debris, p_singlet, p_live)
}
for (s in samples) {
  panels <- ARMS |>
    pmap(\(arm, dir, viability) {
      plot_arm(gating_sets[[arm]], arm, s, viability)
    }) |>
    list_flatten()
  p <- cowplot::plot_grid(plotlist = panels, ncol = 3)
  ggsave(
    file.path(run_dir, "figures", paste0("gates_", s, ".png")),
    p,
    width = 270,
    height = 360,
    units = "mm",
    dpi = 120
  )
}

writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
