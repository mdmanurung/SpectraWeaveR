# Purpose: check which fluorophores of the BMV panel have BOTH a bead and a cell single-stain control, as the first
#          step of a bead-versus-cell spectral mismatch check (AutoSpectral spectral.mismatch.plot() or cosine angle).
#          The script stops with an error when no fluorophore has both, because there is nothing to compare.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/24_bead_cell_mismatch.R <run_dir>
#          from the repo root. <run_dir> is e.g. experiments/runs/2026-11-10_bmv-diagnostics (created if absent;
#          the subfolder 24_bead_cell_mismatch inside it must not exist).
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group/*.fcs                 control file names
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv   fluorophore and control.type
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/*.fcs                 same panel, other run
# Outputs: <run_dir>/24_bead_cell_mismatch/control_pairing.csv   one row per fluorophore: files, types, has_both
#          <run_dir>/24_bead_cell_mismatch/files_found.csv       every control file by folder
#          <run_dir>/24_bead_cell_mismatch/pairing_summary.txt   result in words
#          <run_dir>/24_bead_cell_mismatch/sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
DYES_ASKED <- c(
  "Alexa Fluor 405",
  "Alexa Fluor 488",
  "Alexa Fluor 594",
  "APC-Cy7",
  "PE-Cy5.5",
  "PE-Cy7"
) # the 6 dyes the user expected to have both control types (user, 2026-11-10)

# Setup -----------------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_dir <- file.path(run_dir, "24_bead_cell_mismatch")
stopifnot("24 output folder must not exist" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
})

# Files on disk ---------------------------------------------------------------------------
folders <- c(
  v188 = file.path(DATA, "autoSpectral_v1.8.8", "Reference Group"),
  v188_clean = file.path(DATA, "autoSpectral_v1.8.8", "Reference Group Clean"),
  v102 = file.path(DATA, "autoSpectral_v1.0.2", "Reference Group")
)
files_found <- bind_rows(lapply(names(folders), \(k) {
  tibble(
    folder = k,
    file = list.files(folders[[k]], pattern = "\\.fcs$", recursive = TRUE)
  )
})) |>
  mutate(
    type = case_when(
      grepl("\\(Beads\\)", file) ~ "beads",
      grepl("\\(Cells\\)", file) ~ "cells",
      TRUE ~ "other"
    ),
    fluorophore = sub(" \\((Beads|Cells)\\)\\.fcs$", "", basename(file))
  )
write_csv(files_found, file.path(out_dir, "files_found.csv"))

# Pairing ---------------------------------------------------------------------------------
ctrl <- read_csv(
  file.path(folders[["v188_clean"]], "fcs_control_file.csv"),
  show_col_types = FALSE
)
pairing <- files_found |>
  filter(type %in% c("beads", "cells")) |>
  distinct(fluorophore, type, folder) |>
  group_by(fluorophore) |>
  summarise(
    folders_with_beads = paste(folder[type == "beads"], collapse = ";"),
    folders_with_cells = paste(folder[type == "cells"], collapse = ";"),
    has_beads = any(type == "beads"),
    has_cells = any(type == "cells"),
    .groups = "drop"
  ) |>
  mutate(
    has_both = has_beads & has_cells,
    asked = fluorophore %in% DYES_ASKED
  ) |>
  left_join(
    ctrl |> select(fluorophore, control_file_type = control.type),
    by = "fluorophore"
  )
write_csv(pairing, file.path(out_dir, "control_pairing.csv"))

asked <- pairing |> filter(asked)
msg <- c(
  paste0(
    "dyes asked: ",
    length(DYES_ASKED),
    "; found in the folders: ",
    nrow(asked)
  ),
  paste0("dyes asked with beads: ", sum(asked$has_beads)),
  paste0("dyes asked with cells: ", sum(asked$has_cells)),
  paste0("dyes asked with both: ", sum(asked$has_both)),
  paste0("dyes in the whole panel with both: ", sum(pairing$has_both))
)
writeLines(msg, file.path(out_dir, "pairing_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(as.data.frame(pairing))
cat(msg, sep = "\n")

if (!any(pairing$has_both)) {
  stop(
    "No fluorophore has both a bead and a cell control in the Reference Group folders, ",
    "so no bead-versus-cell spectral mismatch can be computed. See ",
    file.path(out_dir, "control_pairing.csv"),
    call. = FALSE
  )
}
stop(
  "Some fluorophores have both control types, but the mismatch step is not implemented in this script. ",
  "Extend it with get.spectra.automated() per control type before use.",
  call. = FALSE
)
