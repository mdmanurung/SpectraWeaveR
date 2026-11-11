# Purpose: compare the AutoSpectral 1.0.2 and 1.8.8 spectra (P13), check how much the assumed
#          1.8.8 file choice matters, and run the spectra checks on both sets. Nothing is graded.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/04_spectra.R \
#            experiments/runs/2026-11-09_bmv-unmix-compare
#          (start from the repository root; the run dir already exists from 01_manifest-gate.R)
# Inputs:  experiments/bmv_test/compare_3D/data/ (unzipped compare_3Dstim.zip) [user]
#            autoSpectral_v1.0.2/Reference Group/table_spectra/   spectra and AF CSVs [user, ReRunUnmix.R]
#            autoSpectral_v1.8.8/table_spectra/                   Automated, original, Clean, AF CSVs
#          SpectroFlo exports no spectra, so it is not in this script.
# Outputs: <run_dir>/spectra_compare.csv     per fluorophore, 1.0.2 vs 1.8.8: angle_deg, max_abs_delta
#          <run_dir>/spectra_timestamps.csv  same metrics between the 1.8.8 files and timestamps
#          <run_dir>/spectra_checks.csv      sw_unmix_diagnose checks per set, plus condition number
#          <run_dir>/spectra_provenance.csv  file paths, MD5, matrix shape and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
DATA_DIR <- "experiments/bmv_test/compare_3D/data"
DIR_102 <- file.path(
  DATA_DIR,
  "autoSpectral_v1.0.2/Reference Group/table_spectra"
)
DIR_188 <- file.path(DATA_DIR, "autoSpectral_v1.8.8/table_spectra")
FILE_102 <- "Cells Matching Negative_autospectral_spectra.csv" # user-confirmed, from ReRunUnmix.R
FILE_102_AF <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv"
STAMP_188 <- "20261105_143803" # ASSUMED, pending user confirm: which 1.8.8 run is the reference
STAMPS_188 <- c("20261105_142441", "20261105_143455", "20261105_143803") # all three runs
# The Automated, Clean and AF files of one run have stamps a few minutes apart. The AF stamp
# is paired to its run by order (142610, 143604, 143933): an inference, not read from a file.
AF_STAMP_188 <- c(
  "20261105_142441" = "20261105_142610",
  "20261105_143455" = "20261105_143604",
  "20261105_143803" = "20261105_143933"
)
# Fluorophore names that differ between sets: list c("name in 1.0.2" = "name in 1.8.8") here
# after reading the only-in-one-set rows of spectra_compare.csv. None are needed for these files.
NAME_ALIASES <- character(0)

# Run dir ---------------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_files <- file.path(
  run_dir,
  c(
    "spectra_compare.csv",
    "spectra_timestamps.csv",
    "spectra_checks.csv",
    "spectra_provenance.csv"
  )
)
stopifnot("spectra outputs must not exist" = !any(file.exists(out_files)))
stopifnot(
  "run dir must exist (made by 01_manifest-gate.R)" = dir.exists(run_dir)
)

# Read spectra ----------------------------------------------------------------------------
# One CSV is fluorophores in rows, detectors in columns, as AutoSpectral writes it.
read_spectra <- function(path) {
  tab <- read_csv(
    path,
    col_types = cols(.default = col_double(), `...1` = col_character()),
    name_repair = "unique_quiet"
  )
  names(tab)[[1]] <- "fluorophore"
  m <- as.matrix(tab[-1])
  rownames(m) <- tab[["fluorophore"]]
  m
}

# Clean files carry their own stamp, a minute after the Automated one: match by order.
clean_stamps <- str_extract(
  sort(list.files(DIR_188, pattern = "^Clean_")),
  "[0-9]{8}_[0-9]{6}"
)
files_188 <- tribble(
  ~set        , ~prefix                       , ~stamp     , ~file_stamp  ,
  "automated" , "Automated"                   , STAMPS_188 , STAMPS_188   ,
  "original"  , "Automated_original_(unused)" , STAMPS_188 , STAMPS_188   ,
  "clean"     , "Clean"                       , STAMPS_188 , clean_stamps
) |>
  tidyr::unnest(c(stamp, file_stamp)) |>
  mutate(
    path = file.path(
      DIR_188,
      str_c(prefix, "_autospectral_spectra_", file_stamp, ".csv")
    ),
    label = str_c(set, "_", stamp)
  )
spectra_188 <- files_188 |>
  pull(path, name = label) |>
  map(read_spectra)

spectra_102 <- read_spectra(file.path(DIR_102, FILE_102))
spectra_188_assumed <- spectra_188[[str_c("automated_", STAMP_188)]]
af_102 <- read_spectra(file.path(DIR_102, FILE_102_AF))
af_188 <- read_spectra(
  file.path(
    DIR_188,
    str_c(
      "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra ",
      AF_STAMP_188[[STAMP_188]],
      ".csv"
    )
  )
)

# Compare two sets ------------------------------------------------------------------------
# Align by fluorophore name (x names mapped through NAME_ALIASES), keep names found in only
# one set as rows with status only_in_x / only_in_y, and use sw_unmix_compare_spectra() on the rest.
compare_sets <- function(x, y, aliases = NAME_ALIASES) {
  rownames(x) <- coalesce(unname(aliases[rownames(x)]), rownames(x))
  both <- intersect(rownames(x), rownames(y))
  stopifnot("detector names match" = setequal(colnames(x), colnames(y)))
  y <- y[, colnames(x), drop = FALSE]
  tibble(fluorophore = union(rownames(x), rownames(y))) |>
    left_join(
      sw_unmix_compare_spectra(
        x[both, , drop = FALSE],
        y[both, , drop = FALSE]
      ) |>
        select(fluorophore, angle_deg, max_abs_delta),
      by = join_by(fluorophore),
      relationship = "one-to-one"
    ) |>
    mutate(
      status = case_when(
        fluorophore %in% both ~ "both",
        fluorophore %in% rownames(x) ~ "only_in_x",
        .default = "only_in_y"
      )
    )
}

# 1.0.2 vs 1.8.8 --------------------------------------------------------------------------
# x is 1.0.2, y is each 1.8.8 run (user, 2026-11-09: report all three; which run unmixed the data is not proven).
spectra_compare <- STAMPS_188 |>
  set_names() |>
  map(\(stamp) {
    compare_sets(spectra_102, spectra_188[[str_c("automated_", stamp)]]) |>
      mutate(
        status = str_replace(status, "_x$", "_1.0.2") |>
          str_replace("_y$", "_1.8.8"),
        stamp_188 = stamp,
        file_102 = FILE_102,
        file_188 = basename(files_188$path[
          files_188$label == str_c("automated_", stamp)
        ])
      )
  }) |>
  list_rbind()
write_csv(spectra_compare, file.path(run_dir, "spectra_compare.csv"))

# Timestamp sensitivity -------------------------------------------------------------------
# Two questions: do the three Automated files differ, does the Automated file differ from its
# "original (unused)" and Clean siblings. The assumed file is one of the three.
auto_labels <- str_c("automated_", STAMPS_188)
timestamp_pairs <- bind_rows(
  as_tibble(t(combn(auto_labels, 2)), .name_repair = ~ c("x", "y")) |>
    mutate(question = "automated_between_stamps"),
  tibble(
    x = auto_labels,
    y = str_c("original_", STAMPS_188),
    question = "automated_vs_original"
  ),
  tibble(
    x = auto_labels,
    y = str_c("clean_", STAMPS_188),
    question = "automated_vs_clean"
  )
)
spectra_timestamps <- timestamp_pairs |>
  mutate(
    result = map2(x, y, \(a, b) {
      compare_sets(spectra_188[[a]], spectra_188[[b]], aliases = character(0))
    })
  ) |>
  tidyr::unnest(result)
write_csv(spectra_timestamps, file.path(run_dir, "spectra_timestamps.csv"))

# Spectra checks --------------------------------------------------------------------------
# Each set goes through sw_unmix_result() and sw_unmix_diagnose() with no limits, so every value
# is not_judged. diagnose computes max_cosine, max_sif and condition_number in base R.
result_102 <- sw_unmix_result(
  spectra = spectra_102,
  af_spectra = af_102,
  method = "AutoSpectral 1.0.2",
  run_id = "as102"
)
result_188 <- sw_unmix_result(
  spectra = spectra_188_assumed,
  af_spectra = af_188,
  method = "AutoSpectral 1.8.8",
  run_id = "as188"
)
checks <- list(as102 = result_102, as188 = result_188) |>
  map(\(r) sw_unmix_diagnose(r, limits = list())[["checks"]]) |>
  list_rbind(names_to = "set")

# Condition number from AutoSpectral itself, as a cross-check of the base R value. Needs
# AutoSpectral: when it is missing the row stays, with NA and the reason.
has_autospectral <- requireNamespace("AutoSpectral", quietly = TRUE)
condition_autospectral <- list(
  as102 = spectra_102,
  as188 = spectra_188_assumed
) |>
  imap(\(m, set) {
    tibble(
      set = set,
      check_id = "condition_number_autospectral",
      scope = "panel",
      item = "all",
      metric = "condition_number_autospectral",
      value = if (has_autospectral) sw_unmix_condition_number(m) else NA_real_,
      status = "not_judged",
      reason = if (has_autospectral) {
        "AutoSpectral::calculate.condition.number, rounded to 2 decimals"
      } else {
        "skipped: AutoSpectral not installed"
      }
    )
  }) |>
  list_rbind()
spectra_checks <- bind_rows(
  mutate(checks, limit = as.character(limit)),
  condition_autospectral
)
write_csv(spectra_checks, file.path(run_dir, "spectra_checks.csv"))

# Provenance ------------------------------------------------------------------------------
used_files <- c(
  as102_spectra = file.path(DIR_102, FILE_102),
  as102_af = file.path(DIR_102, FILE_102_AF),
  as188_spectra = files_188$path[
    files_188$label == str_c("automated_", STAMP_188)
  ],
  as188_af = file.path(
    DIR_188,
    str_c(
      "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra ",
      AF_STAMP_188[[STAMP_188]],
      ".csv"
    )
  )
)
spectra_provenance <- bind_rows(
  tibble(
    role = names(used_files),
    path = unname(used_files),
    status = c(
      "user-confirmed",
      "user-confirmed",
      "ASSUMED, pending user confirm",
      "inferred pairing by stamp order"
    )
  ) |>
    mutate(md5 = unname(tools::md5sum(path))),
  tibble(
    role = str_c("sensitivity_", files_188$label),
    path = files_188$path,
    status = "read for spectra_timestamps.csv",
    md5 = unname(tools::md5sum(files_188$path))
  )
) |>
  mutate(
    n_rows = map_int(path, \(p) nrow(read_spectra(p))),
    n_detectors = map_int(path, \(p) ncol(read_spectra(p))),
    SpectraWeaveR = as.character(packageVersion("SpectraWeaveR")),
    AutoSpectral = if (has_autospectral) {
      as.character(packageVersion("AutoSpectral"))
    } else {
      NA_character_
    },
    run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
write_csv(spectra_provenance, file.path(run_dir, "spectra_provenance.csv"))
