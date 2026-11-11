# Purpose: report step of 87b_af405_afremove.R on the common SpectroFlo gates (plans/common_gates_plan.md,
#          batch 2). Copy with the build step cut out (no unmixing): it reads the unmixed arm folder and the
#          51c rescore, and redraws the panels. TNF and Th2 readouts and the gap-closed panel are rescored at
#          the fixed cut 1.5 (G7a) on the common labels (rows paired across arms, G8); the squeeze cut is
#          refitted on them. The AF-set and variant panels do not use gates and are redrawn unchanged.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/87b_af405_afremove_common.R   (repo root)
#            -> experiments/runs/2026-11-11_bmv-common-af405-afremove (must not exist)
#            after 51c_score_ablation_common.R wrote 2026-11-11_bmv-common-af405-afremove-score
# Inputs:  experiments/bmv_test/compare_3D/data/unmix_cells/af405_afremove/   unmixed arm and its tables [87b]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra, figure_spectral_variants  [user zip]
#          experiments/runs/2026-11-11_bmv-common-af405-afremove-score/gap_closed.csv   [51c_score_ablation_common.R]
#          experiments/runs/2026-11-11_bmv-common-gated/gs/<arm>   common-gate GatingSets [batch 0]
#          helper87c_marker-scores-common.R; lib/theme.R, lib/palettes.R, lib/plots.R
# Outputs: <report_dir>/variant_spread.{pdf,png,csv}, marker_scores.csv, q2_flags.csv, gap_closed.{pdf,png,csv},
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
ARM <- "af405_afremove"
DYE <- "Alexa Fluor 405"
LANDMARK_FILE <- normalizePath(
  "experiments/runs/2026-11-10_bmv-cell-gating/fcs_control_file_landmark.csv"
)
LANDMARK_N_CELLS <- 2000 # define.gate.landmarks, as chain7 (user, 2026-11-10)
LANDMARK_PERCENTILE <- 70 # same
POSITIVE_N <- 1e6 # above the AF405 positive count, so every cleaned positive is kept (user, 2026-11-10)
SEED <- 42 # repo seed (README: bird.seed + 42)
VARIANT_SLOTS <- c("variants", "delta.list", "delta.norms") # per-dye slots swapped together, as in 45
PIPELINE <- "joint" # zip script
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under joint
REPORT_DIR <- "experiments/runs/2026-11-11_bmv-common-af405-afremove"
SCORE_DIR <- "experiments/runs/2026-11-11_bmv-common-af405-afremove-score"
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
CUT <- 1.5 # negatives below / positives above arcsinh 1.5 (G7a, user, 2026-11-11)
COMMON <- "experiments/runs/2026-11-11_bmv-common-gated" # batch 0 output, read only
STORED_ARCSINH <- c("BV510-A", "BV570-A") # stored arcsinh in the GatingSets, as 19a
POP_PATH <- "/nonDebris/singlets/live/CD3pos/"
DONORS <- c("A3", "A6", "B11", "B2", "D1")
DONOR_SHAPES <- setNames(c(16, 17, 15, 18, 8), DONORS)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| < this fraction of as188_joint's rSD (39b; user, 2026-11-10)
Q2_DELTA <- 0.10 # revised Q2 flag: donor-median rise in squeezed share, same sign in all donors (user)
MARKERS <- c(TNF = "BV510-A", Th2 = "eFluor 450-A")
POPS <- c("CD4T", "CD8T")
JITTER_WIDTH <- 0.15 # display only
JITTER_SEED <- 1 # display only

d188 <- file.path(DATA, "autoSpectral_v1.8.8")
clean_dir <- file.path(d188, "Reference Group Clean")
out_dir <- file.path(DATA, "unmix_cells", ARM)
obj_dir <- file.path(out_dir, "objects")
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra_zip <- read_raw(file.path(
  d188,
  "table_spectra",
  "Automated_autospectral_spectra_20261105_143803.csv"
))
af_zip <- read_raw(file.path(
  d188,
  "table_spectra",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
))
variants_zip <- readRDS(file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
))

# Report: offset and spread, TNF and Th2 with the revised Q2 flag, gap closed ---------------------------
if (STEP == "report") {
  stopifnot("report dir must not exist" = !dir.exists(REPORT_DIR))
  dir.create(REPORT_DIR, recursive = TRUE)
  lib <- Sys.getenv("SW_LIB", "experiments/benchmarks/05_bmv-unmix-compare/lib")
  source(file.path(lib, "theme.R"))
  source(file.path(lib, "palettes.R"))
  source(file.path(lib, "plots.R"))
  source(
    "experiments/benchmarks/05_bmv-unmix-compare/helper87c_marker-scores-common.R"
  )

  # BV510 offset (reference - zip 1.8.8 reference) and variant spread (variant - own reference), script 74 metric
  b510 <- spectra_zip["BV510", ]
  proj_coef <- \(d) as.vector(d %*% b510) / sum(b510^2)
  spectra_used <- read_raw(file.path(out_dir, "spectra_used.csv"))
  variants_used <- readRDS(file.path(out_dir, "variants_used.rds"))
  spread <- list(
    `zip 1.8.8` = list(
      ref = spectra_zip[DYE, ],
      v = variants_zip$variants[[DYE]]
    ),
    `AF removed` = list(
      ref = spectra_used[DYE, ],
      v = variants_used$variants[[DYE]]
    )
  ) |>
    map(\(s) {
      tibble(
        spread_coef = proj_coef(sweep(
          s$v[, names(b510)],
          2,
          s$ref[names(b510)]
        )),
        offset_coef = sum(
          (s$ref[names(b510)] - spectra_zip[DYE, names(b510)]) * b510
        ) /
          sum(b510^2)
      )
    }) |>
    list_rbind(names_to = "set")
  write_csv(
    summarise(
      spread,
      n_variants = n(),
      median_spread = median(spread_coef),
      p90_spread = quantile(spread_coef, 0.9),
      offset = first(offset_coef),
      .by = set
    ),
    file.path(REPORT_DIR, "offset_spread.csv")
  )
  set_pal <- categorical_palette(unique(spread$set))
  p_spread <- plot_box(
    spread,
    group = "set",
    value = "spread_coef",
    pal = set_pal,
    jitter_width = JITTER_WIDTH,
    seed = JITTER_SEED,
    value_label = "BV510 coefficient of variant spread",
    title = "AF405 variants: BV510 part"
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")
  save_panel(
    p_spread,
    spread,
    REPORT_DIR,
    "variant_spread",
    width_mm = 89,
    height_mm = 70
  )

  # TNF and Th2 readouts per donor on as188_joint's CD4T and CD8T rows (helper87c, 39b metrics)
  gs <- load_gs(file.path(COMMON, "gs", "as188_joint")) # common labels are the same in every arm
  donors <- DONORS
  score_arm <- \(arm_name) {
    donors |>
      set_names() |>
      map(\(s) {
        ref <- gs_matrix("as188_joint", s, unname(MARKERS))
        x <- gs_matrix(arm_name, s, unname(MARKERS))
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
  }
  ref_scores <- score_arm("as188_joint") |>
    select(donor, pop, channel, ref_squeezed = squeezed_share)
  marker_scores <- score_arm(ARM) |>
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

  # revised Q2 flag: donor-median rise >= Q2_DELTA and the same sign in all donors (no rSD criterion)
  q2_flags <- marker_scores |>
    summarise(
      median_d_squeezed = median(d_squeezed),
      median_rsd_ratio = median(rsd_ratio),
      same_sign = all(sign(d_squeezed) == sign(median(d_squeezed))),
      .by = c(marker, pop)
    ) |>
    mutate(flagged = median_d_squeezed >= Q2_DELTA & same_sign)
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
    title = "AF removed: CD8T TNF gap closed"
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
