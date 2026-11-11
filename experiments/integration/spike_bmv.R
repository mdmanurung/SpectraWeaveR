# Purpose: real-data regression of the packaged spike-in functions (plan v2 step 22) against the
#          benchmark 05 runs 88 (T cells), 88b (panel 2), 88e (CD3-), 88d (co-spike) and 88h
#          (passes). The saved spike keys enter sw_unmix_spike_build() through `assignment`; the
#          existing unmixed files are scored with sw_unmix_spike_score(). No unmixing runs here.
# Run:     pixi run Rscript experiments/integration/spike_bmv.R <pkg_dir> <repo_root> <out_dir>
#            pkg_dir    package checkout to load
#            repo_root  checkout that holds experiments/bmv_test and experiments/runs (the data)
#            out_dir    new run folder; must not exist
# Inputs:  <repo_root>/experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/{Raw, table_spectra}
#          <repo_root>/experiments/bmv_test/compare_3D/data/spiked/                  88 spiked raw [88 build]
#          <repo_root>/experiments/bmv_test/compare_3D/data/unmix_cells/spk*_*/      unmixed files [88-88h]
#          <repo_root>/experiments/runs/2026-11-10_bmv-spike-in/                     88 keys, summary, baseline
#          <repo_root>/experiments/runs/2026-11-10_bmv-spike-in-panel2{,-build}/      88b
#          <repo_root>/experiments/runs/2026-11-11_bmv-spike-in-cd3neg{,-build}/      88e
#          <repo_root>/experiments/runs/2026-11-11_bmv-cospike-tnf-glut1/            88d summary, baseline
#          <repo_root>/experiments/runs/2026-11-11_bmv-spike-passes/                 88h summary, drift
# Outputs: <out_dir>/build_check.csv      case, marker, sample: max |package spiked - 88 spiked raw|
#          <out_dir>/summary_compare.csv  case, setup, set, pop, sample, level: package vs benchmark
#                                         n, median recovery, share lost, and their differences
#          <out_dir>/baseline_compare.csv case, setup, set: unspiked-cell counts and largest shifts,
#                                         package (arcsinh, all channels) vs benchmark
#          <out_dir>/file_status.csv      case, setup, sample, set: independent or coupled
#          <out_dir>/agreement.csv        one row per case: rows compared, unmatched, largest differences
#          <out_dir>/score_88h_TNF_T.rds  sw_spike_score of 88h TNF in T cells (input to the dashboards)
#          <out_dir>/sessionInfo.txt
# Expected differences (plan v2 "Expected regression differences"): the baseline metric is new
# (arcsinh, all fluorescence channels, share above BASELINE_TOL), so only counts and the
# zero/non-zero pattern are compared with 88's raw marker-channel max; summaries gain quantile
# and file_status columns. Recovery, share lost and n must match exactly.

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
LEVELS <- c(0, 1, 2, 3) # benchmark 05 grid (88, user-approved)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
LOSS_CUT <- 0.5 # "losing most of the spike" (64 readout, inherited by 88-88h)
BASELINE_TOL <- 0.05 # arcsinh shift that counts an unspiked cell as moved; inherited from 88h DRIFT_BINS[1]
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across unmixes of the same raw rows
EXACT_TOL <- 1e-9 # recovery and share differences above this count as disagreement

DATA <- file.path(REPO, "experiments/bmv_test/compare_3D/data")
RUNS <- file.path(REPO, "experiments/runs")
UNMIX <- file.path(DATA, "unmix_cells")
SPECTRA_FILE <- file.path(
  DATA,
  "autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
)
T_MARKERS <- c(TNF = "BV510", Th2 = "eFluor 450", IL10 = "BB700")
PASS_SETTINGS <- c("p5_af1", "p1_af5", "p5_af5")

stopifnot("out dir is new" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE)

# Raw files and spectra, read once -------------------------------------------------------------
read_x <- \(path) {
  flowCore::exprs(flowCore::read.FCS(
    path,
    transformation = FALSE,
    truncate_max_range = FALSE
  ))
}
spectra <- as.matrix(read.csv(SPECTRA_FILE, row.names = 1, check.names = FALSE))
key_tnf <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-spike-in/spike_key_TNF.csv"),
  show_col_types = FALSE
)
raw_files <- distinct(key_tnf, sample, raw_file)
raw <- map(set_names(raw_files$raw_file, raw_files$sample), \(f) {
  read_x(file.path(DATA, "autoSpectral_v1.8.8/Raw", f))
})
unmixed_name <- \(raw_file) sub("\\.fcs$", " AutoSpectral.fcs", raw_file)

# One case: build from a saved key, check against 88's spiked raw, score ----------------------
# `dirs(setup, marker)` gives c(spiked = <folder>, base = <folder>) of the unmixed files.
run_case <- function(
  case,
  key_file,
  marker,
  dye,
  setups,
  dirs,
  spiked_dir = NULL
) {
  key <- read_csv(key_file, show_col_types = FALSE)
  populations <- map(set_names(names(raw)), \(s) {
    k <- filter(key, sample == s)
    split(k$raw_row, k$pop)
  })
  design <- sw_unmix_spike_build(
    raw = raw,
    populations = populations,
    spectra = spectra,
    dyes = set_names(dye, marker),
    levels = LEVELS,
    cofactor = COFACTOR,
    spike_fraction = NULL,
    seed = NULL,
    key_cols = KEY_COLS,
    output_dir = NULL,
    assignment = select(key, sample, raw_row, pop, level)
  )
  build_check <- NULL
  if (!is.null(spiked_dir)) {
    # Package spiked matrix vs the benchmark's spiked raw file, first donor
    s <- names(raw)[1]
    theirs <- read_x(file.path(
      spiked_dir,
      raw_files$raw_file[raw_files$sample == s]
    ))
    ours <- design$data[[marker]][[s]]
    det <- colnames(spectra)
    build_check <- tibble(
      case = case,
      marker = marker,
      sample = s,
      max_abs_diff = max(abs(ours[, det] - theirs[, det])),
      max_rel_diff = max(
        abs(ours[, det] - theirs[, det]) / pmax(abs(theirs[, det]), 1)
      )
    )
  }
  design$data <- NULL
  unmixed <- tidyr::expand_grid(
    setup = setups,
    sample = names(raw),
    set = c(marker, "base")
  ) |>
    mutate(
      path = pmap_chr(list(setup, sample, set), \(setup, sample, set) {
        d <- dirs(setup, marker)
        folder <- if (set == "base") d[["base"]] else d[["spiked"]]
        file.path(
          UNMIX,
          folder,
          unmixed_name(raw_files$raw_file[raw_files$sample == sample])
        )
      })
    )
  stopifnot("every unmixed file exists" = all(file.exists(unmixed$path)))
  score <- sw_unmix_spike_score(
    design,
    unmixed,
    loss_cut = LOSS_CUT,
    baseline_tol = BASELINE_TOL
  )
  list(case = case, score = score, build_check = build_check)
}

# Benchmark tables in the package's column names
compare_summary <- function(case, score, ref) {
  ours <- score$summary |>
    select(
      setup,
      set,
      marker,
      pop,
      sample,
      level,
      n,
      median_recovery,
      share_lost
    )
  ours |>
    full_join(
      ref |>
        rename(
          n_ref = n_cells,
          median_ref = median_recovery,
          lost_ref = share_lost
        ),
      by = join_by(setup, marker, pop, sample, level),
      relationship = "one-to-one"
    ) |>
    mutate(
      case = case,
      d_median = abs(median_recovery - median_ref),
      d_lost = abs(share_lost - lost_ref),
      .before = 1
    )
}

# Cases 88, 88b, 88e and 88h -------------------------------------------------------------------
spk_dirs <- \(setup, marker) {
  c(
    spiked = paste0("spk_", setup, "_", marker),
    base = paste0("spk_", setup, "_base")
  )
}
neg_dirs <- \(setup, marker) {
  c(
    spiked = paste0("spkneg_", setup, "_", marker),
    base = paste0("spk_", setup, "_base")
  )
}
pass_dirs <- \(neg) {
  \(setup, marker) {
    if (setup == "p1_af1") {
      return(if (neg) neg_dirs("c188", marker) else spk_dirs("c188", marker))
    }
    set <- if (neg) paste0("neg_", marker) else marker
    c(
      spiked = paste0("spkpass_", setup, "_", set),
      base = paste0("spkpass_", setup, "_base")
    )
  }
}
panel2_markers <- c(
  ATGL = NA,
  ATP5A = NA,
  CD71 = NA,
  CPT1A = NA,
  CS = NA,
  G6PD = NA,
  GLUT1 = NA,
  IFNg = NA,
  IL2 = NA
)
neg_markers <- c(names(T_MARKERS), names(panel2_markers))
dye_of <- \(key_file) {
  unique(read_csv(key_file, show_col_types = FALSE, col_select = "dye")$dye)
}

cases <- c(
  imap(T_MARKERS, \(dye, m) {
    list(
      case = "88",
      key = file.path(
        RUNS,
        "2026-11-10_bmv-spike-in",
        paste0("spike_key_", m, ".csv")
      ),
      marker = m,
      setups = c("c188", "s102_c102", "s188_c102"),
      dirs = spk_dirs,
      cohort = "t",
      spiked_dir = file.path(DATA, "spiked", paste0("2026-11-10_spike_", m))
    )
  }),
  map(names(panel2_markers), \(m) {
    list(
      case = "88b",
      key = file.path(
        RUNS,
        "2026-11-10_bmv-spike-in-panel2-build",
        paste0("spike_key_", m, ".csv")
      ),
      marker = m,
      setups = c("c188", "afrefine"),
      dirs = spk_dirs,
      cohort = "t",
      spiked_dir = file.path(DATA, "spiked", paste0("2026-11-10_spike_", m))
    )
  }),
  map(neg_markers, \(m) {
    list(
      case = "88e",
      key = file.path(
        RUNS,
        "2026-11-11_bmv-spike-in-cd3neg-build",
        paste0("spike_key_", m, ".csv")
      ),
      marker = m,
      setups = c("c188", "afrefine"),
      dirs = neg_dirs,
      cohort = "neg",
      spiked_dir = file.path(
        DATA,
        "spiked",
        paste0("2026-11-11_spike_cd3neg_", m)
      )
    )
  }),
  map(names(T_MARKERS), \(m) {
    list(
      case = "88h",
      key = file.path(
        RUNS,
        "2026-11-10_bmv-spike-in",
        paste0("spike_key_", m, ".csv")
      ),
      marker = m,
      setups = c("p1_af1", PASS_SETTINGS),
      dirs = pass_dirs(FALSE),
      cohort = "t",
      spiked_dir = NULL
    )
  }),
  map(names(T_MARKERS), \(m) {
    list(
      case = "88h",
      key = file.path(
        RUNS,
        "2026-11-11_bmv-spike-in-cd3neg-build",
        paste0("spike_key_", m, ".csv")
      ),
      marker = m,
      setups = c("p1_af1", PASS_SETTINGS),
      dirs = pass_dirs(TRUE),
      cohort = "neg",
      spiked_dir = NULL
    )
  })
)

refs <- list(
  `88` = file.path(RUNS, "2026-11-10_bmv-spike-in"),
  `88b` = file.path(RUNS, "2026-11-10_bmv-spike-in-panel2"),
  `88e` = file.path(RUNS, "2026-11-11_bmv-spike-in-cd3neg"),
  `88h` = file.path(RUNS, "2026-11-11_bmv-spike-passes")
)
ref_summary <- map(refs, \(d) {
  read_csv(file.path(d, "spike_summary.csv"), show_col_types = FALSE) |>
    rename(any_of(c(setup = "setting")))
})

results <- map(cases, \(cs) {
  message(sprintf("case %s, marker %s", cs$case, cs$marker))
  res <- run_case(
    cs$case,
    cs$key,
    cs$marker,
    dye_of(cs$key),
    cs$setups,
    cs$dirs,
    cs$spiked_dir
  )
  ref <- ref_summary[[cs$case]] |>
    filter(
      marker == cs$marker,
      setup %in% cs$setups,
      pop %in% unique(res$score$summary$pop)
    )
  if (cs$case == "88h" && cs$marker == "TNF" && cs$cohort == "t") {
    saveRDS(res$score, file.path(OUT_DIR, "score_88h_TNF_T.rds"))
  }
  list(
    summary = compare_summary(cs$case, res$score, ref) |>
      mutate(cohort = cs$cohort, .after = 1),
    baseline = res$score$baseline_check |>
      mutate(case = cs$case, cohort = cs$cohort, .before = 1),
    file_status = res$score$file_status |>
      mutate(case = cs$case, cohort = cs$cohort, .before = 1),
    build = res$build_check
  )
})

summary_compare <- map(results, "summary") |> list_rbind()
baseline_all <- map(results, "baseline") |> list_rbind()
file_status <- map(results, "file_status") |> list_rbind()
build_check <- map(results, "build") |> list_rbind()

# Case 88d: co-spike on the TNF design ---------------------------------------------------------
key88 <- read_csv(
  file.path(RUNS, "2026-11-10_bmv-spike-in/spike_key_TNF.csv"),
  show_col_types = FALSE
)
design_tnf <- sw_unmix_spike_build(
  raw = raw,
  populations = map(set_names(names(raw)), \(s) {
    k <- filter(key88, sample == s)
    split(k$raw_row, k$pop)
  }),
  spectra = spectra,
  dyes = c(TNF = "BV510"),
  levels = LEVELS,
  cofactor = COFACTOR,
  spike_fraction = NULL,
  seed = NULL,
  key_cols = KEY_COLS,
  output_dir = NULL,
  assignment = select(key88, sample, raw_row, pop, level)
)
design_tnf$data <- NULL
cospike <- sw_unmix_spike_cospike(
  design_tnf,
  marker = "TNF",
  dye2 = c(GLUT1 = "Alexa Fluor 405"),
  levels2 = c(0, 1, 2),
  level1 = 1,
  output_dir = NULL
)
cospike$data <- NULL
co_sets <- cospike$sets$set
co_dir <- set_names(
  paste0("cospike_g", sub("^cospike_GLUT1_", "", co_sets)),
  co_sets
)
co_unmixed <- tidyr::expand_grid(
  setup = c("c188", "afrefine", "s102_c102"),
  sample = names(raw),
  set = c(co_sets, "base")
) |>
  mutate(
    path = pmap_chr(list(setup, sample, set), \(setup, sample, set) {
      folder <- if (set == "base") "base" else co_dir[[set]]
      file.path(
        UNMIX,
        paste0("spk_", setup, "_", folder),
        unmixed_name(raw_files$raw_file[raw_files$sample == sample])
      )
    })
  )
stopifnot(
  "every co-spike unmixed file exists" = all(file.exists(co_unmixed$path))
)
co_score <- sw_unmix_spike_score(
  cospike,
  co_unmixed,
  loss_cut = LOSS_CUT,
  baseline_tol = BASELINE_TOL
)
co_ref <- read_csv(
  file.path(RUNS, "2026-11-11_bmv-cospike-tnf-glut1/cospike_summary.csv"),
  show_col_types = FALSE
)
co_compare <- co_score$summary |>
  filter(level == 1) |>
  mutate(set_ref = unname(co_dir[set])) |>
  full_join(
    rename(co_ref, n_ref = n),
    by = join_by(setup, set_ref == set, pop, sample),
    relationship = "one-to-one"
  ) |>
  transmute(
    case = "88d",
    cohort = "t",
    setup,
    set,
    marker,
    pop,
    sample,
    level,
    n,
    n_ref,
    median_recovery,
    median_ref = median_tnf_recovery,
    share_lost,
    lost_ref = pct_tnf_lost / 100,
    d_median = abs(median_recovery - median_ref),
    d_lost = abs(share_lost - lost_ref),
    d_median_2 = abs(median_recovery_2 - median_glut1_recovery)
  )
summary_compare <- bind_rows(summary_compare, co_compare)
baseline_all <- bind_rows(
  baseline_all,
  co_score$baseline_check |> mutate(case = "88d", cohort = "t", .before = 1)
)
file_status <- bind_rows(
  file_status,
  co_score$file_status |> mutate(case = "88d", cohort = "t", .before = 1)
)

# Baseline: package (arcsinh, all channels) vs benchmark ----------------------------------------
ours_baseline <- baseline_all |>
  filter(group == "level0") |>
  summarise(
    n_unspiked = sum(n[channel == channel[1]]),
    max_shift_any = max(max_abs_diff, na.rm = TRUE),
    share_moved_max = max(share_moved, na.rm = TRUE),
    .by = c(case, cohort, setup, set)
  )
# 88, 88b, 88e: raw max |spiked - base| in the marker channel of level-0 cells, per pop
ref_baseline <- map(refs[c("88", "88b", "88e")], \(d) {
  read_csv(file.path(d, "baseline_check.csv"), show_col_types = FALSE)
}) |>
  list_rbind(names_to = "case") |>
  summarise(
    n_unspiked_ref = sum(n),
    max_marker_raw_ref = max(max_abs_diff),
    .by = c(case, setup, marker)
  ) |>
  rename(set = marker)
# 88h: largest arcsinh shift over all fluorescence channels (drift mode), per pop
ref_drift <- read_csv(
  file.path(refs[["88h"]], "drift_summary.csv"),
  show_col_types = FALSE
) |>
  summarise(
    n_unspiked_ref = sum(n_unspiked),
    max_shift_any_ref = max(max_any),
    .by = c(setting, marker, pop)
  )
drift_by_set <- ref_drift |>
  mutate(group = if_else(pop == "CD3neg", "neg", "t")) |>
  summarise(
    n_unspiked_ref = sum(n_unspiked_ref),
    max_shift_any_ref = max(max_shift_any_ref),
    .by = c(setting, marker, group)
  )
baseline_compare <- ours_baseline |>
  left_join(
    ref_baseline,
    by = join_by(case, setup, set),
    relationship = "many-to-one"
  ) |>
  left_join(
    drift_by_set |> mutate(case = "88h"),
    by = join_by(case, setup == setting, set == marker, cohort == group),
    relationship = "many-to-one"
  )
write_csv(baseline_compare, file.path(OUT_DIR, "baseline_compare.csv"))

# 88h coupling: one p1_af5 case must be coupled, p5_af1 independent -----------------------------
fs88h <- filter(file_status, case == "88h")
coupling <- fs88h |>
  summarise(
    n_files = n(),
    n_coupled = sum(file_status == "coupled"),
    .by = setup
  )

# Agreement per case ------------------------------------------------------------------------------
agreement <- summary_compare |>
  summarise(
    rows = n(),
    unmatched_ours = sum(is.na(n_ref)),
    unmatched_ref = sum(is.na(n)),
    n_differ = sum(n != n_ref, na.rm = TRUE),
    max_d_median = max(d_median, na.rm = TRUE),
    max_d_lost = max(d_lost, na.rm = TRUE),
    within_tol = all(c(d_median, d_lost) <= EXACT_TOL, na.rm = TRUE),
    .by = case
  )
write_csv(summary_compare, file.path(OUT_DIR, "summary_compare.csv"))
write_csv(file_status, file.path(OUT_DIR, "file_status.csv"))
write_csv(build_check, file.path(OUT_DIR, "build_check.csv"))
write_csv(agreement, file.path(OUT_DIR, "agreement.csv"))
write_csv(coupling, file.path(OUT_DIR, "coupling_88h.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

print(agreement)
print(coupling)
print(build_check)
# Never fail silently: the stated regression conditions stop the run when they break.
stopifnot(
  "summaries match the benchmark" = all(agreement$within_tol),
  "no unmatched summary rows" = all(
    agreement$unmatched_ours == 0 & agreement$unmatched_ref == 0
  ),
  "n matches" = all(agreement$n_differ == 0),
  "a p1_af5 file is coupled" = coupling$n_coupled[coupling$setup == "p1_af5"] >
    0,
  "p5_af1 is independent" = coupling$n_coupled[coupling$setup == "p5_af1"] == 0
)
