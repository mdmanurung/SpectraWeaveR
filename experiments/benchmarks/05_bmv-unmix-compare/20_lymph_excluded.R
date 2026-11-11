# Purpose: check what the lymphocyte gate (flowClust K = 2 on FSC-A,SSC-A) leaves out. For live events
#          inside and outside the lymph gate, report % CD3+, and among CD3+ % CD4+, % CD8+ and % TNF+,
#          using the same openCyto cuts as the gated hierarchy (cuts are read from the gates).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/20_lymph_excluded.R <run_dir>
#          from the repo root.
# Inputs:  <run_dir>/manifest.csv [12_cleanup.R]; the cleaned FCS files it lists
# Outputs: <run_dir>/lymph_excluded.csv                    per arm, sample, inside/outside: counts and %
#          <run_dir>/figures_lymph_excluded/cd3_cd8_<arm>_<sample>.png  CD3 vs CD8, inside vs outside
#          <run_dir>/sessionInfo_20.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2

# Bioconductor first, tidyverse parts last. rlang is not attached: its exprs() would mask flowCore's.
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
DEBRIS_RANGE <- "c(2e5, 1.5e6)" # as 01 (user, 2026-11-09)
LINEAGE_RANGE <- "c(0, 3)" # as 03 (user, 2026-11-09)
CYTOKINE_RANGE <- "c(1, 4)" # as 03 (user, grill Q1)
LYMPH_ARGS <- "K=2" # flowClust, largest cluster (user, 2026-11-09)
CD3 <- "SBUV510-A" # channels as in 03_downstream.R
CD4 <- "Spark Blue 574-A"
CD8 <- "SBUV445-A"
TNF <- "BV510-A"
CHECK_ARMS <- c("as188_joint", "as102_var") # one arm per AutoSpectral version; scatter is identical
BINS <- 100 # display only
T_PARENT <- "/nonDebris/singlets/live/lymph/CD3pos"

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
fig_dir <- file.path(run_dir, "figures_lymph_excluded")
out_file <- file.path(run_dir, "lymph_excluded.csv")
stopifnot(
  "20 outputs must not exist" = !dir.exists(fig_dir) && !file.exists(out_file)
)
dir.create(fig_dir)
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)

# Gate one arm: same hierarchy as 10_ggcyto_strategy.R -------------------------------------------
gate_arm <- function(arm_name) {
  files <- manifest |>
    filter(arm == arm_name) |>
    arrange(sample)
  lin <- paste0("gate_range=", LINEAGE_RANGE)
  template <- tibble(
    alias = c(
      "nonDebris",
      "singlets",
      "live",
      "lymph",
      "CD3pos",
      "CD4T",
      "CD8T",
      "TNFpos"
    ),
    pop = c("+", "+", "-", "+", "+", "+-", "-+", "+"),
    parent = c(
      "root",
      "nonDebris",
      "singlets",
      "live",
      "lymph",
      "CD3pos",
      "CD3pos",
      "CD4T"
    ),
    dims = c(
      "FSC-A",
      "FSC-A,FSC-H",
      files$viability[[1]],
      "FSC-A,SSC-A",
      CD3,
      paste(CD4, CD8, sep = ","),
      paste(CD4, CD8, sep = ","),
      TNF
    ),
    gating_method = c(
      "mindensity",
      "singletGate",
      "mindensity",
      "flowClust",
      "mindensity",
      "mindensity",
      "refGate",
      "mindensity"
    ),
    gating_args = c(
      paste0("gate_range=", DEBRIS_RANGE),
      "wider_gate=TRUE",
      NA,
      LYMPH_ARGS,
      lin,
      lin,
      paste0(T_PARENT, "/", CD4, "+:", T_PARENT, "/", CD8, "+"),
      paste0("gate_range=", CYTOKINE_RANGE)
    ),
    collapseDataForGating = FALSE,
    groupBy = NA,
    preprocessing_method = NA,
    preprocessing_args = NA
  )
  template_file <- file.path(
    run_dir,
    paste0("template_lymph_excluded_", arm_name, ".csv")
  )
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

# 1D cut of a mindensity gate: whichever bound is finite --------------------------------------
gate_cut <- function(gh, node) {
  g <- gh_pop_get_gate(gh, node)
  v <- c(g@min, g@max)
  unname(v[is.finite(v)])
}

# Inside vs outside the lymph gate, per arm and sample --------------------------------------
res <- list()
for (a in CHECK_ARMS) {
  gs <- gate_arm(a)
  for (s in sampleNames(gs)) {
    gh <- gs[[s]]
    e <- flowCore::exprs(gh_pop_get_data(gh, "root"))
    cut <- c(
      CD3 = gate_cut(gh, "CD3pos"),
      CD4 = gate_cut(gh, paste0(T_PARENT, "/", CD4, "+")),
      CD8 = gate_cut(gh, paste0(T_PARENT, "/", CD8, "+")),
      TNF = gate_cut(gh, "TNFpos")
    )
    live <- gh_pop_get_indices(gh, "live")
    lymph <- gh_pop_get_indices(gh, "lymph")
    d <- tibble(
      CD3 = e[, CD3],
      CD4 = e[, CD4],
      CD8 = e[, CD8],
      TNF = e[, TNF],
      group = if_else(lymph, "inside lymph", "outside lymph")
    )[live, ]
    res[[length(res) + 1]] <- d |>
      summarise(
        n_live = n(),
        pct_cd3pos = 100 * mean(CD3 > cut[["CD3"]]),
        pct_cd4_of_cd3 = 100 *
          mean(
            CD4[CD3 > cut[["CD3"]]] > cut[["CD4"]] &
              CD8[CD3 > cut[["CD3"]]] <= cut[["CD8"]]
          ),
        pct_cd8_of_cd3 = 100 *
          mean(
            CD8[CD3 > cut[["CD3"]]] > cut[["CD8"]] &
              CD4[CD3 > cut[["CD3"]]] <= cut[["CD4"]]
          ),
        pct_tnf_of_cd3 = 100 * mean(TNF[CD3 > cut[["CD3"]]] > cut[["TNF"]]),
        .by = group
      ) |>
      mutate(arm = a, sample = s, .before = 1)
    p <- ggplot(d, aes(CD3, CD8)) +
      geom_hex(bins = BINS) +
      scale_fill_viridis_c(trans = "sqrt") +
      geom_vline(xintercept = cut[["CD3"]], colour = "red") +
      geom_hline(yintercept = cut[["CD8"]], colour = "red") +
      facet_wrap(vars(group)) +
      labs(
        title = paste(a, s, "- live events, CD3 vs CD8 (red = openCyto cuts)"),
        x = "CD3 (arcsinh)",
        y = "CD8 (arcsinh)"
      ) +
      theme_bw(base_size = 9)
    ggsave(
      file.path(fig_dir, paste0("cd3_cd8_", a, "_", s, ".png")),
      p,
      width = 9,
      height = 4.5,
      dpi = 120,
      bg = "white"
    )
  }
}
write_csv(list_rbind(res), out_file)
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_20.txt")
)
