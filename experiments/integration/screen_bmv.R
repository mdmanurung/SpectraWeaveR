# Purpose: real-data regression of sw_unmix_negative_shape() and sw_unmix_marker_screen() (plan v2
#          step 22) against benchmark 05 scripts 75 (all-marker screen) and 86b (revised flag rule).
#          The package recomputes the negative-shape metrics from the same cleaned files and gates as 75,
#          then applies the screen with 86b's rule values. Flags are compared only where 86b did not
#          coerce NA to FALSE (86b:121); there the package gives NA by design (grill R11).
# Run:     pixi run Rscript experiments/integration/screen_bmv.R <pkg_dir> <repo_root> <out_dir>
#            pkg_dir    package checkout to load
#            repo_root  checkout that holds experiments/bmv_test and experiments/runs (the data)
#            out_dir    new run folder; must not exist
# Inputs:  <repo_root>/experiments/runs/2026-11-10_bmv-gated-pooled/manifest.csv     cleaned FCS [12_cleanup.R]
#          <repo_root>/experiments/runs/2026-11-10_bmv-ablation-gated/gs/as188_joint  gates [19a_gate_arm.R]
#          <repo_root>/experiments/bmv_test/compare_3D/data/{autoSpectral_v1.8.8/AutoSpectral_unmixed,
#            unmix_cells/<arm>}                                                       raw unmixed arms
#          <repo_root>/experiments/runs/2026-11-10_bmv-all-marker-screen/{scores,marker_parent}.csv  [75]
#          <repo_root>/experiments/runs/2026-11-10_bmv-marker-flags-v2/flags.csv      [86b]
# Outputs: <out_dir>/metrics_compare.csv  arm, donor, pop, channel: package vs 75 metrics
#          <out_dir>/flags_compare.csv    arm, pop, channel: package flag vs 86b flag, and why they differ
#          <out_dir>/agreement.csv        counts: metric rows, largest differences, flags equal, expected NA
#          <out_dir>/screen.rds           the sw_marker_screen (input to the dashboards)
#          <out_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  pkgload::load_all(
    commandArgs(trailingOnly = TRUE)[1],
    quiet = TRUE,
    export_all = FALSE
  )
  library(dplyr)
  library(purrr)
  library(readr)
})

# Parameters ---------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
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
POPS <- c("CD4T", "CD8T", "live") # 75
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
CUT <- 1 # negatives below arcsinh 1 (75; user, 2026-11-10)
SQUEEZE_FRAC <- 0.5 # 75 (39b; user, 2026-11-10)
MIN_NEGATIVES <- 1 # helper75 had no floor; 1 reproduces it, so only 86b's NA coercion differs
MIN_DELTA <- 0.10 # 86b revised rule (user, 2026-11-10)
MIN_DONORS <- 5 # 86b: every donor must rise
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A")
EXACT_TOL <- 1e-9

RUNS <- file.path(REPO, "experiments/runs")
RAW_REF_DIR <- file.path(
  REPO,
  "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed"
)
UNMIX_DIR <- file.path(REPO, "experiments/bmv_test/compare_3D/data/unmix_cells")

stopifnot("out dir is new" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE)

# Events per arm and donor in as188_joint's cleaned row order (helper75 arm_matrix) -------------
read_x <- \(path) {
  flowCore::exprs(flowCore::read.FCS(
    path,
    transformation = FALSE,
    truncate_max_range = FALSE
  ))
}
key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))
pooled <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-gated-pooled/manifest.csv"),
  show_col_types = FALSE
) |>
  mutate(file = file.path(REPO, file))
arm_matrix <- function(arm, sample, channels) {
  ref_clean <- pooled$file[pooled$arm == REF_ARM & pooled$sample == sample]
  ref_key <- key_of(read_x(ref_clean))
  if (arm %in% pooled$arm) {
    e <- read_x(pooled$file[pooled$arm == arm & pooled$sample == sample])
  } else {
    raw <- read_x(file.path(UNMIX_DIR, arm, basename(ref_clean)))
    raw_ref_key <- key_of(read_x(file.path(RAW_REF_DIR, basename(ref_clean))))
    stopifnot(
      "raw rows pair with as188_joint raw" = nrow(raw) == length(raw_ref_key)
    )
    e <- raw[raw_ref_key %in% ref_key, , drop = FALSE]
  }
  stopifnot("rows pair with as188_joint" = identical(key_of(e), ref_key))
  e[, c(KEY_COLS, channels), drop = FALSE]
}

marker_parent <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-all-marker-screen/marker_parent.csv"),
  show_col_types = FALSE
)
channels <- marker_parent$channel
gs <- flowWorkspace::load_gs(file.path(
  RUNS,
  "2026-11-10_bmv-ablation-gated/gs/as188_joint"
))
donors <- sort(pooled$sample[pooled$arm == REF_ARM])

# Negative shape and screen with the package --------------------------------------------------
x <- map(set_names(ARMS), \(a) {
  map(set_names(donors), \(s) arm_matrix(a, s, channels))
})
populations <- map(set_names(donors), \(s) {
  map(set_names(POPS), \(p) {
    which(flowWorkspace::gh_pop_get_indices(gs[[s]], p))
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
screen <- sw_unmix_marker_screen(
  shape,
  min_delta = MIN_DELTA,
  min_donors = MIN_DONORS
)
saveRDS(screen, file.path(OUT_DIR, "screen.rds"))

# Metrics vs 75 -----------------------------------------------------------------------------------
scores75 <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-all-marker-screen/scores.csv"),
  show_col_types = FALSE
)
metrics_compare <- shape$metrics |>
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
  file.path(RUNS, "2026-11-10_bmv-marker-flags-v2/flags.csv"),
  show_col_types = FALSE
) |>
  filter(source == "screen") |>
  mutate(
    # 86b wrote FALSE for a missing component (NA median or NA all_donors_up, coalesced at
    # 86b:121); only rows with every component present are comparable
    comparable = !is.na(median_d_squeezed) & !is.na(all_donors_up)
  )
flags_compare <- screen$flags |>
  select(
    arm = setup,
    pop,
    channel,
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
        n_donors_86b = n_donors,
        median_86b = median_d_squeezed,
        flag_86b = flag,
        comparable
      ),
    by = join_by(arm, pop, channel),
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
