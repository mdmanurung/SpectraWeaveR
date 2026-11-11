# Purpose: rate of impossible double positives per arm and donor, with positivity from openCyto gates
#          (user, 2026-11-09; replaces the unstained-cut version in 18, whose cuts were far too low):
#          CD4+CD8+ in CD3+ (the "++" quadrant of the 2D CD4 vs CD8 gate), TCRgd+ in CD3-, and
#          TCRgd+CD4+ in CD3+. False double positives come from spillover errors, so lower is better.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/21_exclusivity_gates.R <run_dir>
#          from the repo root. Uses the files in <run_dir>/manifest.csv (uncleaned for ranking).
# Inputs:  <run_dir>/manifest.csv [01_manifest-gate.R]
# Outputs: <run_dir>/exclusivity_gates.csv        per arm, sample: parent counts and % double positive
#          <run_dir>/exclusivity_gate_counts.csv  events per gate, arm and sample
#          <run_dir>/figures_exclusivity/gates_<sample>.png  CD4 vs CD8 and TCRgd vs CD3 gates per arm
#          <run_dir>/sessionInfo_21.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto

# Bioconductor first, tidyverse parts last. rlang is not attached: its exprs() would mask flowCore's.
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(ggcyto)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
DEBRIS_RANGE <- "c(2e5, 1.5e6)" # as 01 (user, 2026-11-09)
LINEAGE_RANGE <- "c(0, 3)" # mindensity window for CD3, CD4, CD8, TCRgd, as 03 (user, 2026-11-09)
CD3 <- "SBUV510-A" # channels as in 03_downstream.R
CD4 <- "Spark Blue 574-A"
CD8 <- "SBUV445-A"
TCRGD <- "BUV395-A"
T_PARENT <- "/nonDebris/singlets/live/CD3pos" # CD3+ node path, for the refGates
BINS <- 100 # display only

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
fig_dir <- file.path(run_dir, "figures_exclusivity")
out_files <- file.path(
  run_dir,
  c("exclusivity_gates.csv", "exclusivity_gate_counts.csv")
)
stopifnot(
  "21 outputs must not exist" = !dir.exists(fig_dir) &&
    !any(file.exists(out_files))
)
dir.create(fig_dir)
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
arms <- unique(manifest$arm)

# Gate one arm -------------------------------------------------------------------------------
# CD4T creates the hidden 1D gates "<CD4>+" and "<CD8>+" under CD3pos; DP and TCRgd_CD4 reuse them
# through refGate, so every quadrant of one arm shares one CD4 cut and one CD8 cut.
gate_arm <- function(arm_name) {
  files <- manifest |>
    filter(arm == arm_name) |>
    arrange(sample)
  lin <- paste0("gate_range=", LINEAGE_RANGE)
  template <- tribble(
    ~alias                                                    , ~pop , ~parent     , ~dims                        , ~gating_method , ~gating_args                        ,
    "nonDebris"                                               , "+"  , "root"      , "FSC-A"                      , "mindensity"   , paste0("gate_range=", DEBRIS_RANGE) ,
    "singlets"                                                , "+"  , "nonDebris" , "FSC-A,FSC-H"                , "singletGate"  , "wider_gate=TRUE"                   ,
    "live"                                                    , "-"  , "singlets"  , files$viability[[1]]         , "mindensity"   , NA                                  ,
    "CD3pos"                                                  , "+"  , "live"      , CD3                          , "mindensity"   , lin                                 ,
    "CD3neg"                                                  , "-"  , "live"      , CD3                          , "mindensity"   , lin                                 ,
    "CD4T"                                                    , "+-" , "CD3pos"    , paste(CD4, CD8, sep = ",")   , "mindensity"   , lin                                 ,
    "DP"                                                      , "++" , "CD3pos"    , paste(CD4, CD8, sep = ",")   , "refGate"      ,
    paste0(T_PARENT, "/", CD4, "+:", T_PARENT, "/", CD8, "+") ,
    "TCRgd_T"                                                 , "+"  , "CD3pos"    , TCRGD                        , "mindensity"   , lin                                 ,
    "TCRgd_CD4"                                               , "++" , "CD3pos"    , paste(TCRGD, CD4, sep = ",") , "refGate"      ,
    paste0(T_PARENT, "/TCRgd_T:", T_PARENT, "/", CD4, "+")    ,
    "TCRgd_CD3neg"                                            , "+"  , "CD3neg"    , TCRGD                        , "mindensity"   , lin
  ) |>
    mutate(
      collapseDataForGating = FALSE,
      groupBy = NA,
      preprocessing_method = NA,
      preprocessing_args = NA
    )
  template_file <- file.path(
    run_dir,
    paste0("template_exclusivity_", arm_name, ".csv")
  )
  write_csv(template, template_file)
  gs <- sw_gate_run(
    fcs_files = files$file,
    gating_template = template_file,
    transform_channels = c(files$viability[[1]], CD3, CD4, CD8, TCRGD),
    cofactor = COFACTOR
  )
  sampleNames(gs) <- files$sample
  gs
}

gating_sets <- arms |>
  set_names() |>
  map(gate_arm)

# Counts and rates ----------------------------------------------------------------------------
counts <- gating_sets |>
  imap(\(gs, arm) {
    gs_pop_get_count_fast(gs) |> as_tibble() |> mutate(arm = arm)
  }) |>
  list_rbind() |>
  mutate(pop = sub(".*/", "", Population))
write_csv(counts, out_files[[2]])

rates <- counts |>
  select(arm, sample = name, pop, Count) |>
  pivot_wider(names_from = pop, values_from = Count) |>
  transmute(
    arm,
    sample,
    n_cd3pos = CD3pos,
    n_cd3neg = CD3neg,
    pct_cd4cd8_dp = 100 * DP / CD3pos,
    pct_tcrgd_in_cd3neg = 100 * TCRgd_CD3neg / CD3neg,
    pct_tcrgd_cd4_dp = 100 * TCRgd_CD4 / CD3pos
  )
write_csv(rates, out_files[[1]])

# Gate-check pages: CD4 vs CD8 (CD4T, DP) and TCRgd vs CD3 on live, per arm --------------------
for (s in sort(unique(manifest$sample))) {
  panels <- gating_sets |>
    imap(\(gs, arm) {
      list(
        ggcyto(
          gs[s],
          aes(x = !!rlang::sym(CD4), y = !!rlang::sym(CD8)),
          subset = "CD3pos"
        ) +
          geom_hex(bins = BINS) +
          geom_gate(c("CD4T", "DP")) +
          geom_stats() +
          labs(title = paste(arm, "CD3+"), x = "CD4", y = "CD8"),
        ggcyto(
          gs[s],
          aes(x = !!rlang::sym(TCRGD), y = !!rlang::sym(CD3)),
          subset = "CD3neg"
        ) +
          geom_hex(bins = BINS) +
          geom_gate("TCRgd_CD3neg") +
          geom_stats() +
          labs(title = paste(arm, "CD3-"), x = "TCRgd", y = "CD3")
      ) |>
        map(\(p) {
          as.ggplot(
            p + theme(strip.text = element_blank(), legend.position = "none")
          )
        })
    }) |>
    list_flatten()
  ggsave(
    file.path(fig_dir, paste0("gates_", s, ".png")),
    cowplot::plot_grid(plotlist = panels, ncol = 6),
    width = 22,
    height = 9,
    dpi = 100,
    bg = "white"
  )
}
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_21.txt")
)
