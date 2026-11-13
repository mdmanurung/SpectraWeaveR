# Purpose: report step of 87_af405_forced_af.R on the common SpectroFlo gates (plans/common_gates_plan.md,
#          batch 2). Copy with the build step cut out (no unmixing): it reads the unmixed arm folder and the
#          51c rescore, and redraws the panels. TNF and Th2 readouts and the gap-closed panel are rescored at
#          the fixed cut 1.5 (G7a) on the common labels (rows paired across arms, G8); the squeeze cut is
#          refitted on them. The AF-set and variant panels do not use gates and are redrawn unchanged.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/87_af405_forced_af_common.R   (repo root)
#            -> run_out("common-af405-forced-af") (must not exist)
#            after 51c_score_ablation_common.R wrote the run registered as common-af405-forced-af-score
#          Stops on 24H (one AF refit per condition: which one the AF panel shows is not decided) and on EXV (no Th2).
# Inputs:  data_path("unmix_cells", "af405_afrefine")   unmixed arm and its tables [87]
#          data_path("autoSpectral_v1.8.8", "table_spectra"), figure_spectral_variants  [user zip]; the spectra and AF
#            file names are cfg("stamps")$as188_spectra and $as188_af
#          run_in("common-af405-forced-af-score")/gap_closed.csv   [51c_score_ablation_common.R]
#          run_in("common-gated")/gs/<arm>   common-gate GatingSets [batch 0]
#          datasets/<SW_DATASET>.json (stamps, channels, targets, pos_rule, neg_rule), <ds>_samples.csv
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
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Config lookups ---------------------------------------------------------------------------
# a nested config value that the dataset may lack: stop naming this script and the value
need <- \(group, name) {
  v <- cfg(group)[[name]]
  if (is.null(v)) {
    stop(
      "87_af405_forced_af_common.R: needs ",
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
    "87_af405_forced_af_common.R: unmix_cells/<arm>/<condition>/ holds one AF refit per condition; ",
    "which refit the AF panel shows is not decided for dataset ",
    DS$dataset,
    call. = FALSE
  )
}

# Parameters ---------------------------------------------------------------------------------------
STEP <- "report" # build and score ran before; this copy only reports on the common gates
DATA <- normalizePath(data_path())
ARM <- "af405_afrefine"
PIPELINE <- "joint" # zip script
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under joint
REPORT_DIR <- run_out("common-af405-forced-af")
SCORE_DIR <- run_in("common-af405-forced-af-score")
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
Q2_DELTA <- 0.10 # grill Q2: donor-median rise in squeezed share that flags a marker (user, 2026-11-10)
Q2_RSD <- 0.8 # grill Q2: and donor-median negative rSD ratio at or below this
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
    "87_af405_forced_af_common.R: the cut differs between samples, arms or kinds; ",
    "helper87c_marker-scores-common.R takes one CUT",
    call. = FALSE
  )
}
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
  cfg("stamps")$as188_spectra
))
af188 <- read_raw(file.path(
  d188,
  "table_spectra",
  cfg("stamps")$as188_af
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

  # TNF-slot dye coefficient of each AF spectrum (script 74 metric), zip AF set vs refined AF set
  b510 <- spectra[TNF_DYE, ]
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
    value_label = sprintf("%s coefficient of AF spectrum", TNF_DYE),
    title = sprintf("AF spectra: %s part", TNF_DYE)
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
  stopifnot(
    "gating set holds the samples of samples_tbl()" = setequal(
      sampleNames(gs),
      SAMPLE_IDS
    )
  )
  donors <- SAMPLE_IDS
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
    title = sprintf("Refined AF: CD8T %s gap closed", TNF_LABEL)
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
