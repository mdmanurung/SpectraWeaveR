# Purpose: real-data regression of the packaged spike-in functions (plan v2 step 22) against the
#          common-gate benchmark 05 runs 88cg (T cells), 88bcg (panel 2), 88ecg (CD3-), 88dcg (co-spike) and
#          88hcg (passes) of one dataset (SW_DATASET: 3d, 24h or exv). The saved spike keys enter
#          sw_unmix_spike_build() through `assignment`; the existing unmixed files are scored with
#          sw_unmix_spike_score(). No unmixing runs here.
#          Samples: sw_unmix_spike_score() keys every table by a free `sample` label and never reads a donor,
#          but its by_setup table counts and pools those labels as "donors". With conditions (24H, three
#          samples per donor) the build and score run once per condition, so by_setup stays one row per
#          setting and condition; without conditions (3D, EXV) there is one group and nothing changes.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_integration.sbatch \
#            experiments/integration/spike_bmv.R <pkg_dir> <repo_root> <out_dir> <dataset_json> <spiked_stamp>
#            or, from the repo root, SW_DATASET=<ds> pixi run Rscript experiments/integration/spike_bmv.R \
#            <pkg_dir> <repo_root> <out_dir> <dataset_json> <spiked_stamp>
#            pkg_dir       package checkout to load
#            repo_root     checkout that holds experiments/bmv_test and experiments/runs (the data); must be the
#                          working directory, because the registry paths are relative to it
#            out_dir       new run folder; must not exist
#            dataset_json  required, no default: experiments/benchmarks/05_bmv-unmix-compare/datasets/<ds>.json;
#                          its file name must be $SW_DATASET.json
#            spiked_stamp  required, no default: build date (YYYY-MM-DD) of the folders <data_dir>/spiked/<stamp>_cg_spike_*
# Inputs:  <data_dir>/autoSpectral_v1.8.8/{Raw, table_spectra}   raw files; spectra file from config stamps
#          <data_dir>/spiked/<spiked_stamp>_cg_spike_*/              spiked raw [88cg, 88bcg, 88ecg builds]
#          <data_dir>/unmix_cells/{spkcg_,spkcgneg_,spkcgpass_}*  unmixed spiked sets [88cg-88hcg]; spk_<setup>_base and
#                                                                 spkpass_<setting>_base are the old base unmixes
#          datasets/<ds>_samples.csv                              sample and file stem [00_samples.R]
#          config spike_panels (panel1, panel2, cd3neg, passes) and targets (tnf_label, tnf_dye, glut1_dye)
#          registry run keys (datasets/<ds>_runs.csv), all from the common-gate chain:
#            spike-in-cg-build, spike-in-panel2-cg-build, spike-in-cd3neg-cg-build   spike keys [88cg, 88bcg, 88ecg]
#            spike-in-cg, spike-in-panel2-cg, spike-in-cd3neg-cg   spike_summary.csv, baseline_check.csv [88cg, 88bcg, 88ecg]
#            spike-passes-cg                                       spike_summary.csv, drift_summary.csv [88hcg]
#            cospike-tnf-glut1-cg                                  cospike_summary.csv [88dcg]
# Outputs: <out_dir>/build_check.csv      case, marker, sample: max |package spiked - 88 spiked raw|
#          <out_dir>/summary_compare.csv  case, setup, set, pop, sample, level: package vs benchmark
#                                         n, median recovery, share lost, and their differences (24H adds condition)
#          <out_dir>/baseline_compare.csv case, setup, set: unspiked-cell counts and largest shifts,
#                                         package (arcsinh, all channels) vs benchmark
#          <out_dir>/file_status.csv      case, setup, sample, set: independent or coupled (24H adds condition)
#          <out_dir>/agreement.csv        one row per case: rows compared, unmatched, largest differences
#          <out_dir>/coupling_88h.csv     per pass setting: files and coupled files
#          <out_dir>/score_88h_<tnf_label>_T.rds  sw_spike_score of 88h for the TNF-slot marker in T cells
#                                         (input to the dashboards; 24H: one file per condition, suffix _<condition>)
#          <out_dir>/sessionInfo.txt
# Expected differences (plan v2 "Expected regression differences"): the baseline metric is new
# (arcsinh, all fluorescence channels, share above BASELINE_TOL), so only counts and the
# zero/non-zero pattern are compared with the benchmark's raw marker-channel max; summaries gain quantile
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
  library(tibble)
  library(tidyr)
})

# Dataset config -----------------------------------------------------------------------------
# <dataset_json> (4th argument) names the dataset; dataset.R reads its folder (config, registry, sample table).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "spike_bmv.R: needs 5 arguments: <pkg_dir> <repo_root> <out_dir> <dataset_json> <spiked_stamp>",
    call. = FALSE
  )
}
DATASET_JSON <- normalizePath(args[[4]], mustWork = TRUE)
SPIKED_STAMP <- args[[5]] # build date of the spiked raw folders in <data_dir>/spiked, YYYY-MM-DD; no default
stopifnot(
  "spiked_stamp is YYYY-MM-DD" = grepl(
    "^[0-9]{4}-[0-9]{2}-[0-9]{2}$",
    SPIKED_STAMP
  )
)
if (basename(DATASET_JSON) != paste0(Sys.getenv("SW_DATASET"), ".json")) {
  stop(
    "spike_bmv.R: <dataset_json> must be $SW_DATASET.json; got ",
    DATASET_JSON,
    " with SW_DATASET='",
    Sys.getenv("SW_DATASET"),
    "'",
    call. = FALSE
  )
}
Sys.setenv(SW_DATASET_DIR = dirname(DATASET_JSON))
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), samples_tbl()

# Parameters ---------------------------------------------------------------------------------
REPO <- args[2]
OUT_DIR <- args[3]
LEVELS <- c(0, 1, 2, 3) # benchmark 05 grid (88, user-approved)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
LOSS_CUT <- 0.5 # "losing most of the spike" (64 readout, inherited by 88-88h)
BASELINE_TOL <- 0.05 # arcsinh shift that counts an unspiked cell as moved; inherited from 88h DRIFT_BINS[1]
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across unmixes of the same raw rows
EXACT_TOL <- 1e-9 # recovery and share differences above this count as disagreement

DATA <- data_path()
UNMIX <- file.path(DATA, "unmix_cells")
SPECTRA_FILE <- file.path(
  DATA,
  "autoSpectral_v1.8.8/table_spectra",
  cfg("stamps")$as188_spectra
)
TNF_LABEL <- cfg("targets")$tnf_label # marker label of the TNF-slot dye in the spike keys
TNF_DYE <- cfg("targets")$tnf_dye
GLUT1_DYE <- cfg("targets")$glut1_dye
# marker = dye panels of the common-gate spike scripts (config spike_panels; 3D panel1 = TNF BV510, Th2 eFluor 450, IL10 BB700)
SPIKE_PANELS <- cfg("spike_panels")
T_MARKERS <- unlist(SPIKE_PANELS$panel1) # 88cg
PANEL2_MARKERS <- unlist(SPIKE_PANELS$panel2) # 88bcg
NEG_MARKERS <- unlist(SPIKE_PANELS$cd3neg) # 88ecg
PASS_MARKERS <- unlist(SPIKE_PANELS$passes) # 88hcg
PASS_SETTINGS <- c("p5_af1", "p1_af5", "p5_af5")

stopifnot(
  "repo_root is the working directory (registry paths are relative to it)" = identical(
    normalizePath(REPO),
    normalizePath(getwd())
  ),
  "out dir is new" = !dir.exists(OUT_DIR)
)
dir.create(OUT_DIR, recursive = TRUE)

# Raw files and spectra, read once -------------------------------------------------------------
read_x <- \(path) {
  flowCore::exprs(flowCore::read.FCS(
    path,
    transformation = FALSE,
    truncate_max_range = FALSE
  ))
}
# One FCS file per sample in a folder (recursive: 24H keeps one subfolder per condition). The file stem
# (AutoSpectral appends " AutoSpectral" or " AutoSpectral_v1-8-8", cut by the pattern) joins one-to-one to the sample table, so a missing, doubled or
# stray file stops the run.
files_by_stem <- \(dir) {
  tibble(
    file = list.files(
      dir,
      pattern = "[.]fcs$",
      recursive = TRUE,
      full.names = TRUE
    )
  ) |>
    mutate(
      stem = sub(
        " AutoSpectral.*$",
        "",
        tools::file_path_sans_ext(basename(file))
      )
    ) |>
    inner_join(
      select(samples_tbl(), sample, stem),
      by = "stem",
      relationship = "one-to-one",
      unmatched = "error"
    ) |>
    select(sample, file)
}
spectra <- as.matrix(read.csv(SPECTRA_FILE, row.names = 1, check.names = FALSE))
key_tnf <- read_csv(
  file.path(
    run_in("spike-in-cg-build"),
    paste0("spike_key_", TNF_LABEL, ".csv")
  ),
  show_col_types = FALSE
)
raw_samples <- unique(key_tnf$sample) # the keys' sample order
stopifnot(
  "spike key samples are the sample table" = setequal(
    raw_samples,
    samples_tbl()$sample
  )
)
raw_paths <- files_by_stem(file.path(DATA, "autoSpectral_v1.8.8/Raw"))
raw <- map(set_names(raw_samples), \(s) {
  read_x(raw_paths$file[raw_paths$sample == s])
})
# Scoring groups: one per condition (24H), so by_setup counts samples of one condition; otherwise all samples.
SAMPLE_GROUPS <- if (is.null(DS$conditions)) {
  list(all = raw_samples)
} else {
  map(
    split(samples_tbl()$sample, samples_tbl()$condition),
    \(g) raw_samples[raw_samples %in% g]
  )
}

# One case: build from a saved key, check against the benchmark's spiked raw, score -----------
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
  spiked_paths <- if (!is.null(spiked_dir)) files_by_stem(spiked_dir)
  by_group <- map(SAMPLE_GROUPS, \(group) {
    populations <- map(set_names(group), \(s) {
      k <- filter(key, sample == s)
      split(k$raw_row, k$pop)
    })
    design <- sw_unmix_spike_build(
      raw = raw[group],
      populations = populations,
      spectra = spectra,
      dyes = set_names(dye, marker),
      levels = LEVELS,
      cofactor = COFACTOR,
      spike_fraction = NULL,
      seed = NULL,
      key_cols = KEY_COLS,
      output_dir = NULL,
      assignment = select(
        filter(key, sample %in% group),
        sample,
        raw_row,
        pop,
        level
      )
    )
    build_check <- NULL
    if (!is.null(spiked_dir)) {
      # Package spiked matrix vs the benchmark's spiked raw file, first sample of the group
      s <- group[1]
      theirs <- read_x(spiked_paths$file[spiked_paths$sample == s])
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
    # unmixed file of every setup, set and sample of the group, in setup, sample, set order
    unmixed <- expand_grid(setup = setups, set = c(marker, "base")) |>
      mutate(
        dir = map2_chr(setup, set, \(which_setup, which_set) {
          d <- dirs(which_setup, marker)
          file.path(
            UNMIX,
            if (which_set == "base") d[["base"]] else d[["spiked"]]
          )
        }),
        files = map(dir, files_by_stem)
      ) |>
      unnest(files) |>
      filter(sample %in% group) |>
      arrange(
        match(setup, setups),
        match(sample, group),
        match(set, c(marker, "base"))
      ) |>
      select(setup, sample, set, path = file)
    score <- sw_unmix_spike_score(
      design,
      unmixed,
      loss_cut = LOSS_CUT,
      baseline_tol = BASELINE_TOL
    )
    list(score = score, build_check = build_check)
  })
  list(
    case = case,
    scores = map(by_group, "score"),
    summary = map(by_group, \(g) g$score$summary) |> list_rbind(),
    baseline = map(by_group, \(g) g$score$baseline_check) |> list_rbind(),
    file_status = map(by_group, \(g) g$score$file_status) |> list_rbind(),
    build_check = map(by_group, "build_check") |> list_rbind()
  )
}

# Benchmark tables in the package's column names
compare_summary <- function(case, res, ref) {
  ours <- res$summary |>
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
# Unmixed folders of the common-gate scripts: spkcg_<setup>_<marker> (88cg, 88bcg), spkcgneg_ (88ecg), spkcgpass_
# (88hcg); every setup reads the old base unmix spk_<setup>_base, and 88hcg's other settings spkpass_<setting>_base.
spk_dirs <- \(setup, marker) {
  c(
    spiked = paste0("spkcg_", setup, "_", marker),
    base = paste0("spk_", setup, "_base")
  )
}
neg_dirs <- \(setup, marker) {
  c(
    spiked = paste0("spkcgneg_", setup, "_", marker),
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
      spiked = paste0("spkcgpass_", setup, "_", set),
      base = paste0("spkpass_", setup, "_base")
    )
  }
}
neg_markers <- names(NEG_MARKERS)
dye_of <- \(key_file) {
  unique(read_csv(key_file, show_col_types = FALSE, col_select = "dye")$dye)
}
key_of_run <- \(run_key, marker) {
  file.path(run_in(run_key), paste0("spike_key_", marker, ".csv"))
}
# Spiked raw folder <SPIKED_STAMP>_<prefix><marker>; the stamp is the build date the caller names (5th argument).
spiked_of <- \(prefix, marker) {
  out <- file.path(DATA, "spiked", paste0(SPIKED_STAMP, "_", prefix, marker))
  if (!dir.exists(out)) {
    stop("spike_bmv.R: spiked folder not found: ", out, call. = FALSE)
  }
  out
}

cases <- c(
  imap(T_MARKERS, \(dye, m) {
    list(
      case = "88",
      key = key_of_run("spike-in-cg-build", m),
      marker = m,
      setups = c("c188", "s102_c102", "s188_c102"),
      dirs = spk_dirs,
      cohort = "t",
      spiked_dir = spiked_of("cg_spike_", m)
    )
  }),
  map(names(PANEL2_MARKERS), \(m) {
    list(
      case = "88b",
      key = key_of_run("spike-in-panel2-cg-build", m),
      marker = m,
      setups = c("c188", "afrefine"),
      dirs = spk_dirs,
      cohort = "t",
      spiked_dir = spiked_of("cg_spike_", m)
    )
  }),
  map(neg_markers, \(m) {
    list(
      case = "88e",
      key = key_of_run("spike-in-cd3neg-cg-build", m),
      marker = m,
      setups = c("c188", "afrefine"),
      dirs = neg_dirs,
      cohort = "neg",
      spiked_dir = spiked_of("cg_spike_cd3neg_", m)
    )
  }),
  map(names(PASS_MARKERS), \(m) {
    list(
      case = "88h",
      key = key_of_run("spike-in-cg-build", m),
      marker = m,
      setups = c("p1_af1", PASS_SETTINGS),
      dirs = pass_dirs(FALSE),
      cohort = "t",
      spiked_dir = NULL
    )
  }),
  map(names(PASS_MARKERS), \(m) {
    list(
      case = "88h",
      key = key_of_run("spike-in-cd3neg-cg-build", m),
      marker = m,
      setups = c("p1_af1", PASS_SETTINGS),
      dirs = pass_dirs(TRUE),
      cohort = "neg",
      spiked_dir = NULL
    )
  })
)

refs <- list(
  `88` = run_in("spike-in-cg"),
  `88b` = run_in("spike-in-panel2-cg"),
  `88e` = run_in("spike-in-cd3neg-cg"),
  `88h` = run_in("spike-passes-cg")
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
      pop %in% unique(res$summary$pop)
    )
  list(
    # the score objects the dashboards read: 88h, TNF-slot marker, T cells (one per scoring group)
    score_88h = if (
      cs$case == "88h" && cs$marker == TNF_LABEL && cs$cohort == "t"
    ) {
      res$scores
    },
    summary = compare_summary(cs$case, res, ref) |>
      mutate(cohort = cs$cohort, .after = 1),
    baseline = res$baseline |>
      mutate(case = cs$case, cohort = cs$cohort, .before = 1),
    file_status = res$file_status |>
      mutate(case = cs$case, cohort = cs$cohort, .before = 1),
    build = res$build_check
  )
})
score_88h <- results |>
  map("score_88h") |>
  compact()
if (length(score_88h) != 1) {
  stop(
    "spike_bmv.R: expected one 88h score for marker ",
    TNF_LABEL,
    " in T cells, found ",
    length(score_88h),
    "; check config spike_panels$passes in dataset ",
    DS$dataset,
    call. = FALSE
  )
}
iwalk(score_88h[[1]], \(group_score, group_name) {
  saveRDS(
    group_score,
    file.path(
      OUT_DIR,
      paste0(
        "score_88h_",
        TNF_LABEL,
        "_T",
        if (is.null(DS$conditions)) "" else paste0("_", group_name),
        ".rds"
      )
    )
  )
})

summary_compare <- map(results, "summary") |> list_rbind()
baseline_all <- map(results, "baseline") |> list_rbind()
file_status <- map(results, "file_status") |> list_rbind()
build_check <- map(results, "build") |> list_rbind()

# Case 88d: co-spike on the TNF design ---------------------------------------------------------
# one design per scoring group, as above; the sets and their folders are the same in every group
co_groups <- map(SAMPLE_GROUPS, \(group) {
  design_tnf <- sw_unmix_spike_build(
    raw = raw[group],
    populations = map(set_names(group), \(s) {
      k <- filter(key_tnf, sample == s)
      split(k$raw_row, k$pop)
    }),
    spectra = spectra,
    dyes = set_names(TNF_DYE, TNF_LABEL),
    levels = LEVELS,
    cofactor = COFACTOR,
    spike_fraction = NULL,
    seed = NULL,
    key_cols = KEY_COLS,
    output_dir = NULL,
    assignment = select(
      filter(key_tnf, sample %in% group),
      sample,
      raw_row,
      pop,
      level
    )
  )
  design_tnf$data <- NULL
  cospike <- sw_unmix_spike_cospike(
    design_tnf,
    marker = TNF_LABEL,
    dye2 = c(GLUT1 = GLUT1_DYE),
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
  co_unmixed <- expand_grid(
    setup = c("c188", "afrefine", "s102_c102"),
    set = c(co_sets, "base")
  ) |>
    mutate(
      dir = map2_chr(setup, set, \(which_setup, which_set) {
        folder <- if (which_set == "base") "base" else co_dir[[which_set]]
        file.path(
          UNMIX,
          paste0(
            if (which_set == "base") "spk_" else "spkcg_",
            which_setup,
            "_",
            folder
          )
        )
      }),
      files = map(dir, files_by_stem)
    ) |>
    unnest(files) |>
    filter(sample %in% group) |>
    arrange(
      match(setup, c("c188", "afrefine", "s102_c102")),
      match(sample, group),
      match(set, c(co_sets, "base"))
    ) |>
    select(setup, sample, set, path = file)
  co_score <- sw_unmix_spike_score(
    cospike,
    co_unmixed,
    loss_cut = LOSS_CUT,
    baseline_tol = BASELINE_TOL
  )
  list(score = co_score, co_dir = co_dir)
})
co_dir <- co_groups[[1]]$co_dir
co_summary <- map(co_groups, \(g) g$score$summary) |> list_rbind()
co_ref <- read_csv(
  file.path(run_in("cospike-tnf-glut1-cg"), "cospike_summary.csv"),
  show_col_types = FALSE
)
co_compare <- co_summary |>
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
  map(co_groups, \(g) g$score$baseline_check) |>
    list_rbind() |>
    mutate(case = "88d", cohort = "t", .before = 1)
)
file_status <- bind_rows(
  file_status,
  map(co_groups, \(g) g$score$file_status) |>
    list_rbind() |>
    mutate(case = "88d", cohort = "t", .before = 1)
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
# 24H: say which condition each sample row belongs to (the tables stay one row per sample)
if (!is.null(DS$conditions)) {
  sample_condition <- select(samples_tbl(), sample, condition)
  summary_compare <- left_join(
    summary_compare,
    sample_condition,
    by = "sample",
    relationship = "many-to-one",
    unmatched = "error"
  )
  file_status <- left_join(
    file_status,
    sample_condition,
    by = "sample",
    relationship = "many-to-one",
    unmatched = "error"
  )
}
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
  "coupling holds rows for p1_af5 and p5_af1" = all(
    c("p1_af5", "p5_af1") %in% coupling$setup
  )
)
# The coupled / independent result is the 3D finding (benchmark 88h); other datasets report it only.
if (DS$dataset == "3d") {
  stopifnot(
    "a p1_af5 file is coupled" = all(
      coupling$n_coupled[coupling$setup == "p1_af5"] > 0
    ),
    "p5_af1 is independent" = all(
      coupling$n_coupled[coupling$setup == "p5_af1"] == 0
    )
  )
} else {
  message(
    "spike_bmv.R: coupling_88h.csv reported for ",
    DS$dataset,
    "; no result expectation checked"
  )
}
