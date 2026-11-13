# Purpose: report step of 87b_af405_afremove.R on the common SpectroFlo gates (plans/common_gates_plan.md,
#          batch 2). Copy with the build step cut out (no unmixing): it reads the unmixed arm folder and the
#          51c rescore, and redraws the panels. TNF and Th2 readouts and the gap-closed panel are rescored at
#          the fixed cut 1.5 (G7a) on the common labels (rows paired across arms, G8); the squeeze cut is
#          refitted on them. The AF-set and variant panels do not use gates and are redrawn unchanged.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/87b_af405_afremove_common.R   (repo root)
#            -> run_out("common-af405-afremove") (must not exist)
#            after 51c_score_ablation_common.R wrote the run registered as common-af405-afremove-score
#          Stops on 24H (one spectra_used/variants_used per condition: which one the panel shows is not decided)
#          and on EXV (no Th2).
# Inputs:  data_path("unmix_cells", "af405_afremove")   unmixed arm and its tables [87b]
#          data_path("autoSpectral_v1.8.8", "table_spectra"), figure_spectral_variants  [user zip]; the spectra and AF
#            file names are cfg("stamps")$as188_spectra and $as188_af
#          run_in("common-af405-afremove-score")/gap_closed.csv   [51c_score_ablation_common.R]
#          run_in("common-gated")/gs/<arm>   common-gate GatingSets [batch 0]
#          run_in("cell-gating")/fcs_control_file_landmark.csv   (named LANDMARK_FILE; not read in the report step)
#          datasets/<SW_DATASET>.json (stamps, channels, targets, pos_rule, neg_rule), <ds>_samples.csv
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
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Config lookups ---------------------------------------------------------------------------
# a nested config value that the dataset may lack: stop naming this script and the value
need <- \(group, name) {
  v <- cfg(group)[[name]]
  if (is.null(v)) {
    stop(
      "87b_af405_afremove_common.R: needs ",
      name,
      ", missing in dataset ",
      DS$dataset,
      call. = FALSE
    )
  }
  v
}
if (!is.null(DS$conditions)) {
  stop(
    "87b_af405_afremove_common.R: unmix_cells/<arm>/<condition>/ holds one spectra_used and variants_used per ",
    "condition; which one the panel shows is not decided for dataset ",
    DS$dataset,
    call. = FALSE
  )
}

# Parameters ---------------------------------------------------------------------------------------
STEP <- "report" # build and score ran before; this copy only reports on the common gates
DATA <- normalizePath(data_path())
ARM <- "af405_afremove"
DYE <- cfg("targets")$glut1_dye
LANDMARK_FILE <- normalizePath(
  file.path(run_in("cell-gating"), "fcs_control_file_landmark.csv")
)
LANDMARK_N_CELLS <- 2000 # define.gate.landmarks, as chain7 (user, 2026-11-10)
LANDMARK_PERCENTILE <- 70 # same
POSITIVE_N <- 1e6 # above the AF405 positive count, so every cleaned positive is kept (user, 2026-11-10)
SEED <- 42 # repo seed (README: bird.seed + 42)
VARIANT_SLOTS <- c("variants", "delta.list", "delta.norms") # per-dye slots swapped together, as in 45
PIPELINE <- "joint" # zip script
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under joint
REPORT_DIR <- run_out("common-af405-afremove")
SCORE_DIR <- run_in("common-af405-afremove-score")
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
COMMON <- run_in("common-gated") # batch 0 output, read only
TNF_LABEL <- cfg("targets")$tnf_label
TNF_DYE <- cfg("targets")$tnf_dye
STORED_ARCSINH <- c(need("channels", "TNF_slot"), need("channels", "IFNg_slot")) # stored arcsinh in the GatingSets, as 19a
POP_PATH <- "/nonDebris/singlets/live/CD3pos/"
SAMPLE_IDS <- sort(samples_tbl()$sample) # 3D and EXV the well
# label of the positive cut in the figure: the fixed value, or the per-donor rule
POS_LABEL <- if (identical(cfg("pos_rule")$type, "fixed")) {
  sprintf("arcsinh %s", cfg("pos_rule")$value)
} else {
  "the positive cut"
}
SQUEEZE_FRAC <- 0.5 # squeezed = |x| < this fraction of as188_joint's rSD (39b; user, 2026-11-10)
Q2_DELTA <- 0.10 # revised Q2 flag: donor-median rise in squeezed share, same sign in all donors (user)
MARKERS <- set_names(
  c(need("channels", "TNF_slot"), need("channels", "Th2_slot")),
  c(TNF_LABEL, "Th2")
)
# helper87c_marker-scores-common.R takes one CUT (negatives below, positives above) for every sample, arm and channel
CUT <- tidyr::expand_grid(
  arm = c("as188_joint", ARM),
  sample = SAMPLE_IDS,
  channel = unname(MARKERS),
  kind = c("pos", "neg")
) |>
  mutate(cut = pmap_dbl(list(arm, sample, channel, kind), cut_of)) |>
  pull(cut) |>
  unique()
if (length(CUT) != 1) {
  stop(
    "87b_af405_afremove_common.R: the cut differs between samples, arms or kinds; ",
    "helper87c_marker-scores-common.R takes one CUT",
    call. = FALSE
  )
}
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
  cfg("stamps")$as188_spectra
))
af_zip <- read_raw(file.path(
  d188,
  "table_spectra",
  cfg("stamps")$as188_af
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

  # TNF-slot dye offset (reference - zip 1.8.8 reference) and variant spread (variant - own reference), script 74 metric
  b510 <- spectra_zip[TNF_DYE, ]
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
    value_label = sprintf("%s coefficient of variant spread", TNF_DYE),
    title = sprintf("AF405 variants: %s part", TNF_DYE)
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
  stopifnot(
    "gating set holds the samples of samples_tbl()" = setequal(
      sampleNames(gs),
      SAMPLE_IDS
    )
  )
  donors <- SAMPLE_IDS
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
        sprintf("%s squeezed near 0", TNF_LABEL),
        sprintf("%s above %s", TNF_LABEL, POS_LABEL)
      )
    )
  gap_pal <- categorical_palette(c(
    sprintf("%s squeezed near 0", TNF_LABEL),
    sprintf("%s above %s", TNF_LABEL, POS_LABEL)
  ))
  p_gap <- plot_box(
    gap,
    group = "readout",
    value = "frac_closed",
    pal = gap_pal,
    jitter_width = JITTER_WIDTH,
    seed = JITTER_SEED,
    value_label = "Fraction of 1.8.8 to 1.0.2 gap closed",
    title = sprintf("AF removed: CD8T %s gap closed", TNF_LABEL)
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
