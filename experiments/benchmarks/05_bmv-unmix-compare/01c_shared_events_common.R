# Purpose: shared_events for 34 on the common SpectroFlo gates (plans/common_gates_plan.md, audit note and G10b).
#          01 kept the events that are live in every arm's own gating. Here shared_events is the SpectroFlo
#          `live` label, mapped to the raw rows 34 reads (34 indexes each manifest file's raw rows with
#          shared[[sample]]). A raw row is TRUE when its clean3 twin is live in the labels. Stops unless every
#          raw file of the 11 arms pairs with the Raw SpectroFlo file on Time, FSC-A, FSC-H and SSC-A.
# Run:     SW_DATASET=3d SW_RUN_DATE=2026-11-11 pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/01c_shared_events_common.R
#          from the repo root; on SLURM: sbatch experiments/slurm/run_common_gs.sbatch 01c_shared_events_common.R
#          SW_DATASET is 3d, 24h or exv and SW_RUN_DATE is YYYY-MM-DD (neither has a default).
# Inputs:  registry keys of datasets/<SW_DATASET>_runs.csv, files under data_path():
#          run_in("unmix-compare-final")/manifest.csv   raw files of the arms [01]
#          run_in("unmix-compare-final")/shared_live_counts.csv   old counts, for comparison [01]
#          run_in("common-gates-clean3")/labels.csv.gz   sample, row, raw_row, live [97d]
#          data_path("autoSpectral_v1.8.8", "Raw")/**/*.fcs   raw SpectroFlo files, which raw_row indexes;
#                                                              joined to the samples on the file stem
#          datasets/<SW_DATASET>_samples.csv   samples
# Outputs: <out_dir> = run_out("common-final")
#          <out_dir>/manifest.csv            the raw files of every arm per sample, copied from the final run
#          <out_dir>/shared_events.rds       per sample: logical vector over raw rows, TRUE = SpectroFlo live
#          <out_dir>/shared_live_counts.csv  sample, n_shared, n_total, n_shared_old (01 per-arm gating)
#          <out_dir>/pairing_raw.csv         arm, sample, n_rows_equal, time/fsca/fsch/ssca_equal vs the Raw file
#          <out_dir>/sessionInfo_01c.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- run_out("common-final") # new dated dir; must not exist
SRC_RUN <- run_in("unmix-compare-final") # the dir 34 was run on
LABELS <- file.path(run_in(cfg("label_key")), "labels.csv.gz") # 3D: 97d labels; 24H and ex vivo: 99a labels (read_labels)
RAW_DIR <- data_path("autoSpectral_v1.8.8", "Raw") # as 97d_spectroflo_labels_clean3.R
PAIR_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # must be equal on every raw row
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE)
manifest <- read_csv(file.path(SRC_RUN, "manifest.csv"), show_col_types = FALSE)
write_csv(manifest, file.path(OUT_DIR, "manifest.csv"))
samples <- sort(unique(manifest$sample), method = "radix") # C order, as arrange()

# Pairing of every raw arm file with the Raw SpectroFlo file ---------------------------------------
read_pair <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))[,
    PAIR_COLS,
    drop = FALSE
  ]
}
# raw file of each sample: stem join, one file per sample (24H Raw has one subfolder per condition)
raw_file_of <- tibble(
  file = list.files(
    RAW_DIR,
    pattern = "\\.fcs$",
    full.names = TRUE,
    recursive = TRUE
  )
) |>
  mutate(
    stem = basename(file) |>
      str_remove("\\.fcs$") |>
      str_remove("[ _]*AutoSpectral.*$")
  ) |>
  inner_join(
    select(samples_tbl(), sample, stem),
    by = join_by(stem),
    relationship = "one-to-one",
    unmatched = "error"
  )
stopifnot(
  "manifest samples are the sample table" = setequal(
    samples,
    raw_file_of$sample
  )
)
raw_ref <- map(
  set_names(raw_file_of$file, raw_file_of$sample)[samples],
  read_pair
)
paired <- parallel::mclapply(
  seq_len(nrow(manifest)),
  \(i) {
    m <- read_pair(manifest$file[[i]])
    ref <- raw_ref[[manifest$sample[[i]]]]
    same_n <- nrow(m) == nrow(ref)
    same_col <- \(col) same_n && identical(unname(m[, col]), unname(ref[, col]))
    tibble(
      arm = manifest$arm[[i]],
      sample = manifest$sample[[i]],
      n_rows_equal = same_n,
      time_equal = same_col("Time"),
      fsca_equal = same_col("FSC-A"),
      fsch_equal = same_col("FSC-H"),
      ssca_equal = same_col("SSC-A")
    )
  },
  mc.cores = N_CORES
)
stopifnot(
  "a pairing worker failed" = !any(map_lgl(paired, \(o) {
    inherits(o, "try-error") || is.null(o)
  }))
)
pairing <- list_rbind(paired)
write_csv(pairing, file.path(OUT_DIR, "pairing_raw.csv"))
stopifnot(
  "every raw arm file pairs with the Raw SpectroFlo file" = all(unlist(select(
    pairing,
    -arm,
    -sample
  )))
)

# SpectroFlo live on raw rows ----------------------------------------------------------------------
labels <- read_labels(c("sample", "raw_row", "live")) |>
  split(~sample)
stopifnot("labels cover the donors" = setequal(names(labels), samples))
shared <- map(set_names(samples), \(s) {
  n_raw <- nrow(raw_ref[[s]])
  live_raw <- labels[[s]]$raw_row[labels[[s]]$live]
  stopifnot(
    "raw_row is unique and inside the raw file" = !anyDuplicated(live_raw) &&
      all(live_raw %in% seq_len(n_raw))
  )
  seq_len(n_raw) %in% live_raw
})
saveRDS(shared, file.path(OUT_DIR, "shared_events.rds"))

shared_counts <- tibble(
  sample = samples,
  n_shared = map_int(shared, sum),
  n_total = map_int(shared, length)
) |>
  left_join(
    read_csv(
      file.path(SRC_RUN, "shared_live_counts.csv"),
      show_col_types = FALSE
    ) |>
      select(sample, n_shared_old = n_shared),
    by = join_by(sample),
    relationship = "one-to-one",
    unmatched = "error"
  )
write_csv(shared_counts, file.path(OUT_DIR, "shared_live_counts.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo_01c.txt")
)
