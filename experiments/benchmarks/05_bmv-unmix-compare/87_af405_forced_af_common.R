# Purpose: report step of 87_af405_forced_af.R on the common SpectroFlo gates (plans/common_gates_plan.md,
#          batch 2). Copy with the build step cut out (no unmixing): it reads the unmixed arm folder and the
#          51c rescore, and redraws the panels. TNF and Th2 readouts and the gap-closed panel are rescored at
#          the fixed cut 1.5 (G7a) on the common labels (rows paired across arms, G8); the squeeze cut is
#          refitted on them. The AF-set and variant panels do not use gates and are redrawn unchanged.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/87_af405_forced_af_common.R   (repo root)
#            -> experiments/runs/2026-11-11_bmv-common-af405-forced-af (must not exist)
#            after 51c_score_ablation_common.R wrote 2026-11-11_bmv-common-af405-forced-af-score
# Inputs:  experiments/bmv_test/compare_3D/data/unmix_cells/af405_afrefine/   unmixed arm and its tables [87]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra, figure_spectral_variants  [user zip]
#          experiments/runs/2026-11-11_bmv-common-af405-forced-af-score/gap_closed.csv   [51c_score_ablation_common.R]
#          experiments/runs/2026-11-11_bmv-common-gated/gs/<arm>   common-gate GatingSets [batch 0]
#          helper87c_marker-scores-common.R; lib/theme.R, lib/palettes.R, lib/plots.R
# Outputs: <report_dir>/af_bv510_coef.{pdf,png,csv}, marker_scores.csv, q2_flags.csv, gap_closed.{pdf,png,csv},
#          sessionInfo.txt (and offset_spread.csv for 87b, af_bv510_coef.csv for 87)

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ---------------------------------------------------------------------------------------
STEP <- "report" # build and score ran before; this copy only reports on the common gates
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
ARM <- "af405_afrefine"
PIPELINE <- "joint" # zip script
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under joint
REPORT_DIR <- "experiments/runs/2026-11-11_bmv-common-af405-forced-af"
SCORE_DIR <- "experiments/runs/2026-11-11_bmv-common-af405-forced-af-score"
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
CUT <- 1.5 # negatives below / positives above arcsinh 1.5 (G7a, user, 2026-11-11)
COMMON <- "experiments/runs/2026-11-11_bmv-common-gated" # batch 0 output, read only
STORED_ARCSINH <- c("BV510-A", "BV570-A") # stored arcsinh in the GatingSets, as 19a
POP_PATH <- "/nonDebris/singlets/live/CD3pos/"
DONORS <- c("A3", "A6", "B11", "B2", "D1")
DONOR_SHAPES <- setNames(c(16, 17, 15, 18, 8), DONORS)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| < this fraction of as188_joint's rSD (39b; user, 2026-11-10)
Q2_DELTA <- 0.10 # grill Q2: donor-median rise in squeezed share that flags a marker (user, 2026-11-10)
Q2_RSD <- 0.8 # grill Q2: and donor-median negative rSD ratio at or below this
MARKERS <- c(TNF = "BV510-A", Th2 = "eFluor 450-A")
POPS <- c("CD4T", "CD8T")
JITTER_WIDTH <- 0.15 # display only
JITTER_SEED <- 1 # display only

d188 <- file.path(DATA, "autoSpectral_v1.8.8")
out_dir <- file.path(DATA, "unmix_cells", ARM)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra <- read_raw(file.path(
  d188,
  "table_spectra",
  "Automated_autospectral_spectra_20261105_143803.csv"
))
af188 <- read_raw(file.path(
  d188,
  "table_spectra",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
))

# Report: AF set BV510 part, TNF and Th2 readouts with the Q2 flag, gap-closed panel -------------------
if (STEP == "report") {
  stopifnot("report dir must not exist" = !dir.exists(REPORT_DIR))
  dir.create(REPORT_DIR, recursive = TRUE)
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")
  source("experiments/benchmarks/05_bmv-unmix-compare/lib/plots.R")
  source(
    "experiments/benchmarks/05_bmv-unmix-compare/helper87c_marker-scores-common.R"
  )

  # BV510 coefficient of each AF spectrum (script 74 metric), zip AF set vs refined AF set
  b510 <- spectra["BV510", ]
  af_refined <- read_raw(file.path(out_dir, "af_refined.csv"))
  af_coef <- list(
    `zip AF` = af188,
    `refined AF` = af_refined
  ) |>
    map(\(m) {
      tibble(
        af_spectrum = rownames(m),
        bv510_coef = as.vector(m[, names(b510)] %*% b510) / sum(b510^2)
      )
    }) |>
    list_rbind(names_to = "af_set")
  write_csv(af_coef, file.path(REPORT_DIR, "af_bv510_coef.csv"))
  af_pal <- categorical_palette(unique(af_coef$af_set))
  p_af <- plot_box(
    af_coef,
    group = "af_set",
    value = "bv510_coef",
    pal = af_pal,
    jitter_width = JITTER_WIDTH,
    seed = JITTER_SEED,
    value_label = "BV510 coefficient of AF spectrum",
    title = "AF spectra: BV510 part"
  )
  save_panel(
    p_af,
    af_coef,
    REPORT_DIR,
    "af_bv510_coef",
    width_mm = 89,
    height_mm = 70
  )

  # TNF and Th2 readouts per donor on as188_joint's CD4T and CD8T rows (helper87c, 39b metrics)
  gs <- load_gs(file.path(COMMON, "gs", "as188_joint")) # common labels are the same in every arm
  donors <- DONORS
  marker_scores <- donors |>
    set_names() |>
    map(\(s) {
      ref <- gs_matrix("as188_joint", s, unname(MARKERS))
      x <- gs_matrix(ARM, s, unname(MARKERS))
      map(set_names(POPS), \(p) {
        score_channels(
          x,
          ref,
          which(gh_pop_get_indices(gs[[s]], paste0(POP_PATH, p)))
        )
      }) |>
        list_rbind(names_to = "pop")
    }) |>
    list_rbind(names_to = "donor")
  ref_scores <- donors |>
    set_names() |>
    map(\(s) {
      ref <- gs_matrix("as188_joint", s, unname(MARKERS))
      map(set_names(POPS), \(p) {
        score_channels(
          ref,
          ref,
          which(gh_pop_get_indices(gs[[s]], paste0(POP_PATH, p)))
        )
      }) |>
        list_rbind(names_to = "pop")
    }) |>
    list_rbind(names_to = "donor") |>
    select(donor, pop, channel, ref_squeezed = squeezed_share)
  marker_scores <- marker_scores |>
    left_join(
      ref_scores,
      by = join_by(donor, pop, channel),
      relationship = "one-to-one",
      unmatched = "error"
    ) |>
    mutate(
      marker = names(MARKERS)[match(channel, MARKERS)],
      d_squeezed = squeezed_share - ref_squeezed
    )
  write_csv(marker_scores, file.path(REPORT_DIR, "marker_scores.csv"))

  # grill Q2 flag: median rise >= Q2_DELTA, median rSD ratio <= Q2_RSD, same sign of the rise in all donors
  q2_flags <- marker_scores |>
    summarise(
      median_d_squeezed = median(d_squeezed),
      median_rsd_ratio = median(rsd_ratio),
      same_sign = all(sign(d_squeezed) == sign(median(d_squeezed))),
      .by = c(marker, pop)
    ) |>
    mutate(
      flagged = median_d_squeezed >= Q2_DELTA &
        median_rsd_ratio <= Q2_RSD &
        same_sign
    )
  write_csv(q2_flags, file.path(REPORT_DIR, "q2_flags.csv"))

  # CD8T TNF gap closed per donor (51 scorer output), donors as points over the median
  gap <- read_csv(
    file.path(SCORE_DIR, "gap_closed.csv"),
    show_col_types = FALSE
  ) |>
    filter(arm == ARM) |>
    mutate(
      readout = if_else(
        readout == "squeezed_share",
        "TNF squeezed near 0",
        "TNF above arcsinh 1.5"
      )
    )
  gap_pal <- categorical_palette(c(
    "TNF squeezed near 0",
    "TNF above arcsinh 1.5"
  ))
  p_gap <- plot_box(
    gap,
    group = "readout",
    value = "frac_closed",
    pal = gap_pal,
    jitter_width = JITTER_WIDTH,
    seed = JITTER_SEED,
    value_label = "Fraction of 1.8.8 to 1.0.2 gap closed",
    title = "Refined AF: CD8T TNF gap closed"
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")
  save_panel(
    p_gap,
    gap,
    REPORT_DIR,
    "gap_closed",
    width_mm = 89,
    height_mm = 70
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(REPORT_DIR, "sessionInfo.txt")
  )
}
