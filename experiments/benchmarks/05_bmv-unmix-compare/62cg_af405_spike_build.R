# Copy of 62_af405_spike_build.R for the common SpectroFlo labels (plans/common_gates_plan.md G5, G10b): CD4T and CD8T
#          cells come from experiments/runs/2026-11-11_bmv-common-gates-clean3/labels.csv.gz (raw_row key) instead of the
#          as188_joint gates. Levels, seed 42 and the BV510 spike are unchanged. Writes a new run dir and a new spiked
#          dir (the old 62 had no must-not-exist guard on the run dir; this copy has one).
# Purpose: proof chain step 5 (user, 2026-11-10): spike a known BV510 (TNF) signal into real raw events.
#          For CD4T and CD8T cells (as188_joint gates), add k x the 1.8.8 BV510 spectrum to the raw detector
#          values, with k = 6000 x sinh(level) for levels 1, 2, 3 (arcsinh units, cofactor 6000). Spectra rows
#          have max 1, so an exact unmix raises BV510 by k. Per donor, T cells are shuffled (set.seed(42),
#          user-approved) and split in quarters: unspiked, level 1, level 2, level 3. Other events untouched.
# Run:     sbatch experiments/slurm/run_spike_cg.sbatch 62cg_af405_spike_build.R
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Raw/*.fcs            raw files [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed  row keys [user zip]
#          experiments/runs/2026-11-11_bmv-common-gates-clean3/labels.csv.gz              [97d_spectroflo_labels_clean3.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/..._143803.csv
# Outputs: experiments/bmv_test/compare_3D/data/spiked/2026-11-11_cg_af405_proof/*.fcs  spiked raw files
#          experiments/runs/2026-11-11_bmv-af405-proof-cg/spike_key.csv   sample, raw_row, pop, level, k
#          experiments/runs/2026-11-11_bmv-af405-proof-cg/sessionInfo_62.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
SMOKE_DONOR <- Sys.getenv("SW_SMOKE_DONOR", "") # one whole donor (all its rows) for a smoke run; empty = full run
SMOKE_TAG <- if (nzchar(SMOKE_DONOR)) "smoke_" else "" # smoke runs write to their own dirs
SMOKE_SUFFIX <- if (nzchar(SMOKE_DONOR)) "-smoke" else ""
RUN_DIR <- paste0(
  "experiments/runs/2026-11-11_bmv-af405-proof-cg",
  SMOKE_SUFFIX
)
LABELS <- "experiments/runs/2026-11-11_bmv-common-gates-clean3/labels.csv.gz" # SpectroFlo labels, clean3 pooled gates (G10b)
OUT_FCS <- file.path(
  DATA,
  "spiked",
  paste0(SMOKE_TAG, "2026-11-11_cg_af405_proof")
)
LEVELS <- c(0, 1, 2, 3) # arcsinh units added; 0 = unspiked control quarter (user-approved grid)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
SEED <- 42 # user-approved, which cells get which level
POPS <- c("CD4T", "CD8T") # label columns of labels.csv.gz

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
stopifnot("spiked fcs dir must not exist" = !dir.exists(OUT_FCS))
dir.create(RUN_DIR, recursive = TRUE)
dir.create(OUT_FCS, recursive = TRUE)

# BV510 spectrum on the detector channels -------------------------------------------------
s188 <- as.matrix(read.csv(
  file.path(
    DATA,
    "autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
  ),
  row.names = 1,
  check.names = FALSE
))
b510 <- s188["BV510", ]

# CD4T and CD8T raw rows per donor, from the SpectroFlo labels -------------------------------------
labels <- read_csv(
  LABELS,
  col_types = cols_only(
    sample = "c",
    raw_row = "d",
    FSC_A = "d",
    CD4T = "l",
    CD8T = "l"
  )
)
if (nzchar(SMOKE_DONOR)) {
  labels <- filter(labels, sample == SMOKE_DONOR)
}
raw_dir <- file.path(DATA, "autoSpectral_v1.8.8/Raw")
raw_names <- list.files(raw_dir, pattern = "\\.fcs$")

# Spike each donor's raw file and record the key ------------------------------------------------
set.seed(SEED)
spike_key <- map(set_names(unique(labels$sample)), \(s) {
  raw_name <- raw_names[startsWith(raw_names, paste0(s, " "))]
  stopifnot("one raw file per donor" = length(raw_name) == 1)
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
    outer(key$k, b510)
  exprs(ff) <- x
  write.FCS(ff, file.path(OUT_FCS, raw_name))
  key
}) |>
  list_rbind()
write_csv(spike_key, file.path(RUN_DIR, "spike_key.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_62.txt")
)
