# Copy of 62_af405_spike_build.R for the common SpectroFlo labels (plans/common_gates_plan.md G5, G10b): CD4T and CD8T
#          cells come from experiments/runs/2026-11-11_bmv-common-gates-clean3/labels.csv.gz (raw_row key) instead of the
#          as188_joint gates. Levels, seed 42 and the BV510 spike are unchanged. Writes a new run dir and a new spiked
#          dir (the old 62 had no must-not-exist guard on the run dir; this copy has one).
# Purpose: proof chain step 5 (user, 2026-11-10): spike a known TNF-slot dye signal (BV510 in 3D) into real raw events.
#          For CD4T and CD8T cells (SpectroFlo labels), add k x the 1.8.8 dye spectrum to the raw detector
#          values, with k = 6000 x sinh(level) for levels 1, 2, 3 (arcsinh units, cofactor 6000). Spectra rows
#          have max 1, so an exact unmix raises the dye by k. Per sample, T cells are shuffled (set.seed(42),
#          user-approved) and split in quarters: unspiked, level 1, level 2, level 3. Other events untouched.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_spike_cg.sbatch 62cg_af405_spike_build.R
#          SW_SMOKE_DONOR: a donor or sample id (its samples only); empty = full run.
# Inputs:  <data_dir>/autoSpectral_v1.8.8/Raw/**/*.fcs  raw files, joined to datasets/<ds>_samples.csv on the file stem [user zip]
#          <data_dir>/autoSpectral_v1.8.8/table_spectra/<stamps.as188_spectra>  1.8.8 spectra [user zip; dataset config]
#          run key common-gates-clean3: labels.csv.gz  [97d_spectroflo_labels_clean3.R]
#          config targets$tnf_dye: the dye to spike [dataset config]
# Outputs: <data_dir>/spiked/2026-11-11_cg_af405_proof/[<condition>/]*.fcs  spiked raw files
#          run_out("af405-proof-cg")/spike_key.csv   sample, raw_row, pop, level, k (raw_file is the path under Raw/)
#          run_out("af405-proof-cg")/sessionInfo_62.txt
#          The run folder is experiments/runs/<SW_RUN_DATE>_bmv<run_tag>-af405-proof-cg (smoke: -smoke).

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ------------------------------------------------------------------------------
DATA <- data_path()
SMOKE_DONOR <- Sys.getenv("SW_SMOKE_DONOR", "") # donor or sample id for a smoke run: all its samples; empty = full run
SMOKE_TAG <- if (nzchar(SMOKE_DONOR)) "smoke_" else "" # smoke runs write to their own dirs
SMOKE_SUFFIX <- if (nzchar(SMOKE_DONOR)) "-smoke" else ""
RUN_DIR <- run_out(paste0("af405-proof-cg", SMOKE_SUFFIX)) # stops when the folder exists
LABELS <- file.path(run_in(cfg("label_key")), "labels.csv.gz") # 3D: SpectroFlo labels of clean3 (97d); 24H and ex vivo: pooled AS 1.8.8 (99a); read with read_labels()
TNF_DYE <- cfg("targets")$tnf_dye # dye whose spectrum is spiked (BV510 in 3D)
OUT_FCS <- file.path(
  DATA,
  "spiked",
  paste0(SMOKE_TAG, "2026-11-11_cg_af405_proof")
)
LEVELS <- c(0, 1, 2, 3) # arcsinh units added; 0 = unspiked control quarter (user-approved grid)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
SEED <- 42 # user-approved, which cells get which level
POPS <- c("CD4T", "CD8T") # label columns of labels.csv.gz

stopifnot("spiked fcs dir must not exist" = !dir.exists(OUT_FCS))
dir.create(RUN_DIR, recursive = TRUE)
dir.create(OUT_FCS, recursive = TRUE)

# TNF-slot dye spectrum on the detector channels -------------------------------------------
s188 <- as.matrix(read.csv(
  file.path(
    DATA,
    "autoSpectral_v1.8.8/table_spectra",
    cfg("stamps")$as188_spectra
  ),
  row.names = 1,
  check.names = FALSE
))
tnf_spectrum <- s188[TNF_DYE, ]

# CD4T and CD8T raw rows per donor, from the SpectroFlo labels -------------------------------------
labels <- read_labels(c("sample", "raw_row", "FSC_A", "CD4T", "CD8T"))
raw_dir <- file.path(DATA, "autoSpectral_v1.8.8/Raw")

# Raw files of the 1.8.8 arm, joined to the sample table on the full file stem (one file per sample)
raw_files <- tibble(
  raw_name = list.files(raw_dir, pattern = "\\.fcs$", recursive = TRUE)
) |>
  mutate(stem = sub("\\.fcs$", "", basename(raw_name))) |>
  inner_join(
    samples_tbl(),
    by = "stem",
    relationship = "one-to-one",
    unmatched = "error"
  )
if (nzchar(SMOKE_DONOR)) {
  raw_files <- filter(raw_files, donor == SMOKE_DONOR | sample == SMOKE_DONOR)
  if (nrow(raw_files) == 0) {
    stop(
      "62cg_af405_spike_build.R: SW_SMOKE_DONOR ",
      SMOKE_DONOR,
      " matches no donor or sample of dataset ",
      DS$dataset,
      call. = FALSE
    )
  }
  labels <- filter(labels, sample %in% raw_files$sample)
}

# Spike each donor's raw file and record the key ------------------------------------------------
set.seed(SEED)
spike_key <- map(set_names(unique(labels$sample)), \(s) {
  raw_name <- raw_files$raw_name[raw_files$sample == s]
  stopifnot("one raw file per sample" = length(raw_name) == 1)
  donor <- labels[labels$sample == s, ]
  ff <- read.FCS(
    file.path(raw_dir, raw_name),
    transformation = FALSE,
    truncate_max_range = FALSE
  )
  key <- imap(set_names(POPS), \(label_col, p) {
    tibble(
      raw_row = donor$raw_row[donor[[label_col]]],
      FSC_A = donor$FSC_A[donor[[label_col]]],
      pop = p
    )
  }) |>
    list_rbind() |>
    distinct(raw_row, .keep_all = TRUE) |> # CD4T and CD8T are disjoint; keep first if not
    mutate(
      level = sample(rep_len(LEVELS, n())),
      k = COFACTOR * sinh(level),
      sample = s,
      raw_file = raw_name
    )
  stopifnot(
    "label raw_row pairs with the raw file (FSC-A)" = max(abs(
      exprs(ff)[key$raw_row, "FSC-A"] - key$FSC_A
    )) <
      0.5 # float32 spacing at 1e6 is 0.06 to 0.12
  )
  x <- exprs(ff)
  x[key$raw_row, colnames(s188)] <- x[key$raw_row, colnames(s188)] +
    outer(key$k, tnf_spectrum)
  exprs(ff) <- x
  out_file <- file.path(OUT_FCS, raw_name)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE) # condition subfolder in 24H
  write.FCS(ff, out_file)
  key
}) |>
  list_rbind()
write_csv(spike_key, file.path(RUN_DIR, "spike_key.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_62.txt")
)
