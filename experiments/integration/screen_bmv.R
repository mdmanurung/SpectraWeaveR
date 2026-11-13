# Purpose: real-data regression of sw_unmix_negative_shape() and sw_unmix_marker_screen() (plan v2
#          step 22) against the common-gate benchmark 05 scripts 75 (all-marker screen, 75_all_marker_screen_common.R)
#          and 86b (revised flag rule, 86b_apply_flag_rule_v2_common.R) of one dataset (SW_DATASET: 3d, 24h or exv).
#          The package recomputes the negative-shape metrics from the same cleaned files and common gates as 75,
#          then applies the screen with 86b's rule values. Flags are compared only where 86b did not
#          coerce NA to FALSE; there the package gives NA by design (grill R11).
#          Samples: the package keys everything by a free `sample` label, but the screen counts samples as
#          "donors" and asks that all of them rise. With conditions (24H, three samples per donor) the shape and the
#          screen run once per condition, as 75 and 86b score per condition; without conditions (3D, EXV) there is
#          one group and nothing changes. min_donors is the size of the group (86b: every donor must rise).
#          The package takes one negative cut for every arm, sample and channel, so the script stops when the
#          dataset's neg_rule is not fixed (24H: per-donor control cuts).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_integration.sbatch \
#            experiments/integration/screen_bmv.R <pkg_dir> <repo_root> <out_dir> <dataset_json>
#            or, from the repo root, SW_DATASET=<ds> pixi run Rscript experiments/integration/screen_bmv.R \
#            <pkg_dir> <repo_root> <out_dir> <dataset_json>
#            pkg_dir       package checkout to load
#            repo_root     checkout that holds experiments/bmv_test and experiments/runs (the data); must be the
#                          working directory, because the registry paths are relative to it
#            out_dir       new run folder; must not exist
#            dataset_json  required, no default: experiments/benchmarks/05_bmv-unmix-compare/datasets/<ds>.json;
#                          its file name must be $SW_DATASET.json
# Inputs:  registry run keys (datasets/<ds>_runs.csv), all from the common-gate chain:
#            common-gated                 manifest.csv (cleaned FCS of every arm) and gs/as188_joint, the common-gate
#                                         GatingSet [97c, 97d]
#            all-marker-screen-common     scores.csv, marker_parent.csv [75_all_marker_screen_common.R]
#            marker-flags-v2-common       flags.csv [86b_apply_flag_rule_v2_common.R]
#          <data_dir>/{autoSpectral_v1.8.8/AutoSpectral_unmixed, unmix_cells/<arm>}   raw unmixed arms [user zip, 37_swap_af_variants.R]
#          datasets/<ds>_samples.csv    sample and file stem [00_samples.R]
#          config marker_parent (populations) and neg_rule (must be fixed)
#          The arm and sample reading is helper98_marker_scores.R (the helper of 75_all_marker_screen_common.R).
# Outputs: <out_dir>/metrics_compare.csv  arm, donor, pop, channel: package vs 75 metrics
#          <out_dir>/flags_compare.csv    arm, pop, channel: package flag vs 86b flag, and why they differ (24H adds condition)
#          <out_dir>/agreement.csv        counts: metric rows, largest differences, flags equal, expected NA
#          <out_dir>/screen.rds           the sw_marker_screen (input to the dashboards; 24H: screen_<condition>.rds)
#          <out_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  pkgload::load_all(
    commandArgs(trailingOnly = TRUE)[1],
    quiet = TRUE,
    export_all = FALSE
  )
  library(flowCore) # helper98 reads FCS with exprs() and read.FCS(); dplyr below masks filter()
  library(flowWorkspace)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Dataset config -----------------------------------------------------------------------------
# <dataset_json> (4th argument) names the dataset; dataset.R reads its folder (config, registry, sample table).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "screen_bmv.R: needs 4 arguments: <pkg_dir> <repo_root> <out_dir> <dataset_json>",
    call. = FALSE
  )
}
DATASET_JSON <- normalizePath(args[[4]], mustWork = TRUE)
if (basename(DATASET_JSON) != paste0(Sys.getenv("SW_DATASET"), ".json")) {
  stop(
    "screen_bmv.R: <dataset_json> must be $SW_DATASET.json; got ",
    DATASET_JSON,
    " with SW_DATASET='",
    Sys.getenv("SW_DATASET"),
    "'",
    call. = FALSE
  )
}
Sys.setenv(SW_DATASET_DIR = dirname(DATASET_JSON))
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), samples_tbl(), cut_of()
# arm_matrix(), cuts_by_channel(), marker_parent_cfg(), POOLED_MANIFEST, KEY_COLS (as 75_all_marker_screen_common.R)
source("experiments/benchmarks/05_bmv-unmix-compare/helper98_marker_scores.R")

# Parameters ---------------------------------------------------------------------------------
REPO <- args[2]
OUT_DIR <- args[3]
ARMS <- c(
  "as188_joint",
  "as102_var",
  "s102_c188",
  "s188_c102",
  "swap_var102",
  "swap_af102"
) # 75
REF_ARM <- "as188_joint"
MARKER_PARENT <- marker_parent_cfg() # config marker_parent; 75 scores the populations it lists (3D: CD4T, CD8T, live)
POPS <- unique(c(
  MARKER_PARENT$default,
  unlist(MARKER_PARENT$by_marker),
  MARKER_PARENT$no_marker
))
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
SQUEEZE_FRAC <- 0.5 # 75 (39b; user, 2026-11-10)
MIN_NEGATIVES <- 1 # helper75 had no floor; 1 reproduces it, so only 86b's NA coercion differs
MIN_DELTA <- 0.10 # 86b revised rule (user, 2026-11-10)
EXACT_TOL <- 1e-9
COND_BY <- if (is.null(DS$conditions)) character() else "condition" # 24H only: screen and flags per condition

stopifnot(
  "repo_root is the working directory (registry paths are relative to it)" = identical(
    normalizePath(REPO),
    normalizePath(getwd())
  ),
  "out dir is new" = !dir.exists(OUT_DIR)
)
dir.create(OUT_DIR, recursive = TRUE)

# Samples, populations and the negative cut ------------------------------------------------------------
# Cleaned files of every arm, the common gates and the channels of 75 (rows pair with as188_joint's cleaned rows)
pooled <- read_csv(POOLED_MANIFEST, show_col_types = FALSE)
marker_parent <- read_csv(
  file.path(run_in("all-marker-screen-common"), "marker_parent.csv"),
  show_col_types = FALSE
)
channels <- marker_parent$channel
gs <- load_gs(file.path(run_in("common-gated"), "gs/as188_joint"))
donors <- sort(pooled$sample[pooled$arm == REF_ARM])
stopifnot(
  "scored samples are the sample table" = setequal(donors, samples_tbl()$sample)
)
# sw_unmix_negative_shape() takes one cut; a fixed rule gives the same number in every arm, sample and channel
if (!identical(cfg("neg_rule")$type, "fixed")) {
  stop(
    "screen_bmv.R: the package takes one negative cut for every arm, sample and channel, but neg_rule of dataset ",
    DS$dataset,
    " is '",
    cfg("neg_rule")$type,
    "'; a per-sample cut needs a package change",
    call. = FALSE
  )
}
CUT <- cut_of(REF_ARM, donors[[1]], channels[[1]], "neg") # arcsinh cut of the negatives (75 via cut_of; G13)
SAMPLE_GROUPS <- if (is.null(DS$conditions)) {
  list(all = donors)
} else {
  map(
    split(samples_tbl()$sample, samples_tbl()$condition),
    \(g) donors[donors %in% g]
  )
}

# Negative shape and screen with the package, one group at a time ---------------------------------------
screens <- map(SAMPLE_GROUPS, \(group) {
  x <- map(set_names(ARMS), \(a) {
    map(set_names(group), \(s) arm_matrix(a, s, c(KEY_COLS, channels), pooled))
  })
  populations <- map(set_names(group), \(s) {
    map(set_names(POPS), \(p) {
      which(gh_pop_get_indices(gs[[s]], p))
    })
  })
  shape <- sw_unmix_negative_shape(
    x,
    reference = REF_ARM,
    populations = populations,
    channels = channels,
    cofactor = COFACTOR,
    cut = CUT,
    squeeze_frac = SQUEEZE_FRAC,
    min_negatives = MIN_NEGATIVES,
    key_cols = KEY_COLS
  )
  rm(x)
  list(
    shape = shape,
    screen = sw_unmix_marker_screen(
      shape,
      min_delta = MIN_DELTA,
      min_donors = length(group) # 86b: every donor of the group must rise
    )
  )
})
iwalk(screens, \(g, group_name) {
  saveRDS(
    g$screen,
    file.path(
      OUT_DIR,
      if (is.null(DS$conditions)) {
        "screen.rds"
      } else {
        paste0("screen_", group_name, ".rds")
      }
    )
  )
})
shape_metrics <- map(screens, \(g) g$shape$metrics) |> list_rbind()
screen_flags <- map(screens, \(g) g$screen$flags) |>
  list_rbind(
    names_to = if (is.null(DS$conditions)) rlang::zap() else "condition"
  )

# Metrics vs 75 -----------------------------------------------------------------------------------
scores75 <- read_csv(
  file.path(run_in("all-marker-screen-common"), "scores.csv"),
  show_col_types = FALSE
)
metrics_compare <- shape_metrics |>
  select(
    arm = setup,
    donor = sample,
    pop,
    channel,
    n_pop,
    rsd_neg,
    squeezed_share,
    rsd_ratio
  ) |>
  full_join(
    scores75 |>
      select(
        arm,
        donor,
        pop,
        channel,
        n,
        rsd_neg_75 = rsd_neg,
        squeezed_75 = squeezed_share,
        rsd_ratio_75 = rsd_ratio
      ),
    by = join_by(arm, donor, pop, channel),
    relationship = "one-to-one"
  ) |>
  mutate(
    d_squeezed = abs(squeezed_share - squeezed_75),
    d_rsd = abs(rsd_neg - rsd_neg_75),
    d_ratio = abs(rsd_ratio - rsd_ratio_75)
  )
write_csv(metrics_compare, file.path(OUT_DIR, "metrics_compare.csv"))

# Flags vs 86b (screen source) --------------------------------------------------------------------
flags86b <- read_csv(
  file.path(run_in("marker-flags-v2-common"), "flags.csv"),
  show_col_types = FALSE
) |>
  filter(source == "screen") |>
  mutate(
    # 86b wrote FALSE for a missing component (NA median or NA all_donors_up, coalesced at
    # 86b:121); only rows with every component present are comparable
    comparable = !is.na(median_d_squeezed) & !is.na(all_donors_up)
  )
flags_compare <- screen_flags |>
  select(
    arm = setup,
    pop,
    channel,
    any_of(COND_BY),
    n_donors,
    median_d_squeezed,
    all_donors_up,
    flag
  ) |>
  full_join(
    flags86b |>
      select(
        arm,
        pop,
        channel,
        any_of(COND_BY),
        n_donors_86b = n_donors,
        median_86b = median_d_squeezed,
        flag_86b = flag,
        comparable
      ),
    by = c("arm", "pop", "channel", COND_BY),
    relationship = "one-to-one"
  ) |>
  mutate(
    agree = if_else(comparable, flag %in% flag_86b & !is.na(flag), NA),
    why = case_when(
      is.na(n_donors_86b) | is.na(n_donors) ~ "row in one table only",
      !comparable &
        is.na(flag) ~ "86b coerced NA to FALSE; package NA (expected, R11)",
      !comparable ~ "86b coerced NA to FALSE; package not NA",
      is.na(flag) ~ "package NA, 86b complete",
      agree ~ "equal",
      .default = "differs"
    )
  )
write_csv(flags_compare, file.path(OUT_DIR, "flags_compare.csv"))

agreement <- tibble(
  metric_rows = nrow(metrics_compare),
  metric_unmatched = sum(
    is.na(metrics_compare$n_pop) | is.na(metrics_compare$n)
  ),
  n_differ = sum(metrics_compare$n_pop != metrics_compare$n, na.rm = TRUE),
  max_d_squeezed = max(metrics_compare$d_squeezed, na.rm = TRUE),
  max_d_rsd = max(metrics_compare$d_rsd, na.rm = TRUE),
  max_d_ratio = max(metrics_compare$d_ratio, na.rm = TRUE),
  flag_rows = nrow(flags_compare),
  flags_comparable = sum(flags_compare$comparable, na.rm = TRUE),
  flags_equal = sum(flags_compare$agree, na.rm = TRUE),
  flags_true = sum(flags_compare$flag, na.rm = TRUE),
  expected_na = sum(
    flags_compare$why == "86b coerced NA to FALSE; package NA (expected, R11)"
  ),
  other = sum(
    !flags_compare$why %in%
      c("equal", "86b coerced NA to FALSE; package NA (expected, R11)")
  )
)
write_csv(agreement, file.path(OUT_DIR, "agreement.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
print(agreement)
print(count(flags_compare, why))
stopifnot(
  "metrics match 75" = agreement$metric_unmatched == 0 &&
    agreement$max_d_squeezed <= EXACT_TOL,
  "comparable flags all agree with 86b" = agreement$flags_equal ==
    agreement$flags_comparable,
  "no unexplained flag differences" = agreement$other == 0
)
