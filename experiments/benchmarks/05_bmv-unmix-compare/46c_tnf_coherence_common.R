#!/usr/bin/env Rscript
# Purpose: check whether each arm's TNF signal hangs together with the cell's other activation readouts,
#          with no unstimulated control, on the common SpectroFlo gates (G1a, G7a). Two readouts per arm,
#          donor and parent (CD4T, CD8T):
#          (1) % TNF+ among IFNg+ cells, at the fixed cut (arcsinh 1.5 for both cytokines);
#          (2) Spearman rho of TNF, and of IFNg as the reference cytokine, with each activation or effector
#          marker in the panel. Rank correlation needs no further cut.
#          Copy of 46_tnf_coherence.R: new gated dir, cut 1.5, new run dir, Nature-style figures.
# Run:     SW_DATASET=<3d|24h> SW_RUN_DATE=<YYYY-MM-DD> \
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/46c_tnf_coherence_common.R [<run_dir>]
#          from the repo root, after batch 0 has written GATED_DIR/README.md. <run_dir> overrides the registry
#          output run_out("common-cut15-46c"); it must not exist. Stops on a dataset whose slots are not
#          TNF and IFNg (EXV: BV510 is CD1c, BV570 is HLA-DR).
# Inputs:  run_in("common-gated")/manifest.csv, gs/<arm>/  common-gate GatingSets [97_common_gs.R]
#          cut_of(arm, sample, channel, "pos")  positive cut [config pos_rule; cuts.csv from 97e when not fixed]
#          experiments/benchmarks/05_bmv-unmix-compare/datasets/  config, run registry, sample table
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R, report_figures.R  labels, builders
# Outputs: <run_dir>/tnf_among_ifng.csv           % TNF+ of IFNg+ cells per arm, sample, parent
#          <run_dir>/rho_activation.csv            Spearman rho per arm, sample, parent, cytokine, marker
#          <run_dir>/fig_a_tnf_among_ifng.{pdf,png,csv}   panel a, donors as points, median bar
#          <run_dir>/fig_b_rho_activation.{pdf,png,csv}   panel b
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), run_out(), cut_of()
source("experiments/benchmarks/05_bmv-unmix-compare/report_figures.R") # theme, save_panel, builders, ARM_LABELS

# Parameters ------------------------------------------------------------------------------
GATED_DIR <- run_in("common-gated") # batch 0 output (G9: read-only here)
# This script needs the TNF and IFNg slots; a dataset that labels them otherwise has no such target
if (
  !identical(
    c(cfg("targets")$tnf_label, cfg("targets")$ifng_label),
    c("TNF", "IFNg")
  )
) {
  stop(
    "46c_tnf_coherence_common.R: needs TNF and IFNg, missing in dataset ",
    DS$dataset,
    " (slots are labelled ",
    cfg("targets")$tnf_label,
    " and ",
    cfg("targets")$ifng_label,
    ")",
    call. = FALSE
  )
}
# Cytokine slots, as in 19a_gate_arm.R; arcsinh in the GatingSet. The cut comes from cut_of() per cell.
CYTOKINES <- c(TNF = cfg("channels")$TNF_slot, IFNg = cfg("channels")$IFNg_slot)
# Cut text for titles: the fixed value when the rule is fixed (3D: 1.5), else the per-sample cut
CUT_VALUE_TXT <- if (identical(cfg("pos_rule")$type, "fixed")) {
  cfg("pos_rule")$value
} else {
  "per sample"
}
# Every activation or effector marker in the panel (CD69, CD154, CD137 and GzmB are not in it). Rank
# correlation is invariant to the monotone transform, so these channels are used as stored (linear).
ACTIVATION <- c(
  `IL-2` = "BV605-A",
  CD25 = "SBUV605-A",
  CD71 = "Spark YG 593-A",
  CD38 = "SBUV740-A",
  `CTLA-4` = "BV650-A",
  `PD-1` = "PE-Fire 640-A"
)
NODES <- c(
  CD4T = "/nonDebris/singlets/live/CD3pos/CD4T",
  CD8T = "/nonDebris/singlets/live/CD3pos/CD8T"
)
POPULATIONS <- names(NODES)
ARM_ORDER <- unlist(VERSION_ROWS, use.names = FALSE)

run_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(run_dir)) {
  run_dir <- run_out("common-cut15-46c")
}
stopifnot("46c run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)
manifest <- read_csv(
  file.path(GATED_DIR, "manifest.csv"),
  show_col_types = FALSE
)
stopifnot(
  "all 11 figure arms are in the gated manifest" = all(
    ARM_ORDER %in% manifest$arm
  ),
  "gated samples equal the sample table" = setequal(
    manifest$sample,
    samples_tbl()$sample
  )
)

# Per arm, donor and parent: % TNF+ of IFNg+ and Spearman rho -----------------------------------
results <- ARM_ORDER |>
  set_names() |>
  map(\(a) {
    gs <- load_gs(file.path(GATED_DIR, "gs", a))
    sampleNames(gs) |>
      set_names() |>
      map(\(s) {
        POPULATIONS |>
          set_names() |>
          map(\(pop) {
            x <- exprs(gh_pop_get_data(gs[[s]], NODES[[pop]]))
            cut_tnf <- cut_of(a, s, CYTOKINES[["TNF"]], "pos")
            cut_ifng <- cut_of(a, s, CYTOKINES[["IFNg"]], "pos")
            ifng_pos <- x[, CYTOKINES[["IFNg"]]] > cut_ifng
            rho <- CYTOKINES |>
              map(\(cy_ch) {
                map_dbl(ACTIVATION, \(m) {
                  cor(x[, cy_ch], x[, m], method = "spearman")
                })
              })
            list(
              pct = tibble(
                n_ifng_pos = sum(ifng_pos),
                pct_tnf_of_ifng = 100 *
                  mean(x[ifng_pos, CYTOKINES[["TNF"]]] > cut_tnf),
                pct_tnf_of_ifng_neg = 100 *
                  mean(x[!ifng_pos, CYTOKINES[["TNF"]]] > cut_tnf)
              ),
              rho = rho |>
                map(\(r) enframe(r, name = "marker", value = "rho")) |>
                list_rbind(names_to = "cytokine")
            )
          })
      })
  })

# Long tables ------------------------------------------------------------------------------
# flatten arm -> donor -> parent into rows
tnf_among_ifng <- results |>
  map(\(by_s) {
    map(by_s, \(by_p) list_rbind(map(by_p, "pct"), names_to = "parent")) |>
      list_rbind(names_to = "donor")
  }) |>
  list_rbind(names_to = "arm")
rho_activation <- results |>
  map(\(by_s) {
    map(by_s, \(by_p) list_rbind(map(by_p, "rho"), names_to = "parent")) |>
      list_rbind(names_to = "donor")
  }) |>
  list_rbind(names_to = "arm")
write_csv(tnf_among_ifng, file.path(run_dir, "tnf_among_ifng.csv"))
write_csv(rho_activation, file.path(run_dir, "rho_activation.csv"))

# Figure a: % TNF+ among IFNg+ and IFNg- cells, donors as points over a median bar ------------
redraw_panel(
  file.path(run_dir, "tnf_among_ifng.csv"),
  "tnf_of_ifng",
  paste0("TNF+ among IFNg+ and IFNg- cells, cut ", CUT_VALUE_TXT),
  run_dir,
  "fig_a_tnf_among_ifng"
)

# Figure b: Spearman rho with activation markers, TNF and IFNg as rows ---------------------------
redraw_panel(
  file.path(run_dir, "rho_activation.csv"),
  "rho_marker",
  "Cytokine-marker Spearman rho",
  run_dir,
  "fig_b_rho_activation"
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
