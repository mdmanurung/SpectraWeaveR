# Purpose: build the sample table of one dataset (3D, 24H or EXV) from the delivered files, so that every
#          later script joins arm files to samples on the file stem and never parses a well code. One sample
#          is one acquisition file. Stops unless every delivered arm folder holds exactly one file per sample.
# Run:     SW_DATASET=3d pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/00_samples.R
#          from the repo root. SW_DATASET is 3d, 24h or exv (no default). It stops when the csv exists, unless
#          SW_OVERWRITE_SAMPLES=1. SW_RUN_DATE is not used.
# Inputs:  <data_dir>/{manual,autoSpectral_v1.0.2/AutospectralUnmix,autoSpectral_v1.0.2/AutospectralUnmix_optimized,
#          autoSpectral_v1.8.8/AutoSpectral_unmixed,autoSpectral_v1.8.8/Raw}/**/*.fcs   delivered files [user, zip]
#          datasets/<SW_DATASET>.json   data_dir, raw_layout, conditions, donor_regex
# Outputs: datasets/<SW_DATASET>_samples.csv   sample, donor, condition, well, acq_batch, stem
#            sample     3D and EXV: the well code; 24H: <donor>_<condition> (a well code repeats across conditions)
#            stem       file name without ".fcs" and without the trailing " AutoSpectral..." arm suffix, so the
#                       manual, 1.0.2 and 1.8.8 files of one sample share it

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ---------------------------------------------------------------------------------------
# delivered arm folders under data_dir; "raw" is the reference that lists the samples
ARM_DIRS <- c(
  spectroflo = "manual",
  as102_basic = "autoSpectral_v1.0.2/AutospectralUnmix",
  as102_var = "autoSpectral_v1.0.2/AutospectralUnmix_optimized",
  as188_joint = "autoSpectral_v1.8.8/AutoSpectral_unmixed",
  raw = "autoSpectral_v1.8.8/Raw"
)
OUT_CSV <- file.path(.sw05_dir, paste0(DS$dataset, "_samples.csv"))

if (file.exists(OUT_CSV) && Sys.getenv("SW_OVERWRITE_SAMPLES") != "1") {
  stop(
    "00_samples.R: ",
    OUT_CSV,
    " exists; set SW_OVERWRITE_SAMPLES=1 to overwrite it",
    call. = FALSE
  )
}
if (!dir.exists(data_path())) {
  stop(
    "00_samples.R: data dir missing for dataset ",
    DS$dataset,
    ": ",
    data_path(),
    call. = FALSE
  )
}

# Files of every arm -------------------------------------------------------------------------------
# recursive: the 24H Raw, manual and 1.8.8 folders have one subfolder per condition
files <- imap(ARM_DIRS, \(dir, arm) {
  tibble(
    arm = arm,
    file = list.files(
      data_path(dir),
      pattern = "\\.fcs$",
      recursive = TRUE,
      full.names = TRUE
    )
  )
}) |>
  list_rbind()

stem_of <- \(path) {
  path |>
    basename() |>
    str_remove("\\.fcs$") |>
    str_remove("[ _]*AutoSpectral.*$")
}

# Samples from the Raw files -----------------------------------------------------------------------
cond_regex <- if (is.null(DS$conditions)) {
  NA_character_
} else {
  paste0("_24h(", paste(DS$conditions, collapse = "|"), ")_")
}
samples <- files |>
  filter(arm == "raw") |>
  mutate(
    name = basename(file),
    donor = str_extract(name, cfg("donor_regex")),
    well = str_extract(name, "^[A-Z][0-9]+"),
    acq_batch = str_match(name, "_BMV_Aurora_(.+)_Batch1")[, 2],
    condition = if (is.na(cond_regex)) {
      NA_character_
    } else {
      str_match(name, cond_regex)[, 2]
    },
    stem = stem_of(file)
  )
if (anyNA(select(samples, donor, well, acq_batch, stem))) {
  stop(
    "00_samples.R: donor, well or batch not found in a Raw file name (donor_regex ",
    cfg("donor_regex"),
    ")",
    call. = FALSE
  )
}
if (!is.na(cond_regex)) {
  # the condition in the file name must be the subfolder it sits in
  folder_cond <- basename(dirname(samples$file))
  if (anyNA(samples$condition) || cfg("raw_layout") != "by_condition") {
    stop(
      "00_samples.R: condition missing in a Raw file name, or raw_layout is not by_condition",
      call. = FALSE
    )
  }
  if (!identical(samples$condition, folder_cond)) {
    stop(
      "00_samples.R: file-name condition differs from the subfolder for: ",
      paste(samples$name[samples$condition != folder_cond], collapse = "; "),
      call. = FALSE
    )
  }
}
samples <- samples |>
  mutate(
    sample = if (cfg("raw_layout") == "flat") {
      well
    } else {
      paste(donor, condition, sep = "_")
    }
  ) |>
  select(sample, donor, condition, well, acq_batch, stem) |>
  arrange(sample)
if (anyDuplicated(samples$sample) || anyDuplicated(samples$stem)) {
  stop(
    "00_samples.R: sample or stem is not unique in the Raw files",
    call. = FALSE
  )
}

# One file per sample in every arm -----------------------------------------------------------------
arm_check <- files |>
  mutate(stem = stem_of(file)) |>
  summarise(
    n_files = n(),
    n_stems = n_distinct(stem),
    n_known = sum(stem %in% samples$stem),
    .by = arm
  ) |>
  filter(n_files != nrow(samples) | n_stems != n_files | n_known != n_files)
missing_arm <- setdiff(names(ARM_DIRS), unique(files$arm))
if (nrow(arm_check) > 0 || length(missing_arm) > 0) {
  stop(
    "00_samples.R: arms without exactly one file per sample (",
    nrow(samples),
    " samples): ",
    paste(c(arm_check$arm, missing_arm), collapse = ", "),
    call. = FALSE
  )
}

write_csv(samples, OUT_CSV)
