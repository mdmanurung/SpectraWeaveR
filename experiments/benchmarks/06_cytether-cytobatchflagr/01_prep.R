#!/usr/bin/env Rscript
# Purpose: read the CytoBatchFlagR controls, subsample, transform, and assign frozen FlowSOM
#          populations per leave-one-donor-out fold (plan: plans/cytether_t6_cytobatchflagr_plan.md).
# Run:     pixi run Rscript experiments/benchmarks/06_cytether-cytobatchflagr/01_prep.R <run_dir>
# Inputs:  $SW_DATA_DIR/cytobatchflagr_zenodo/*.fcs   72 control files (Zenodo 15388817)
#          $SW_DATA_DIR/cytobatchflagr_zenodo/Table_2.csv   file, batch, control, sample_id
# Outputs (in <run_dir>, which must not exist):
#   cells.rds            one row per kept cell: sample_id, batch, control, then the 26 markers (arcsinh)
#   markers.csv          channel name and FCS desc label of the 26 markers
#   populations.rds      one integer column per fold (fold_tcu1 ...): FlowSOM metacluster of each cell
#   population_counts.csv  cells per sample, fold and metacluster
#   samples.csv          sample metadata with cells read and kept
#   sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(flowCore)
  library(FlowSOM)
})

# Parameters (approved, plan T6) -------------------------------------------------
CELLS_PER_FILE <- 50000 # N3; files with fewer cells contribute all of them
SUBSAMPLE_SEED <- 1 # N6
# N2b: the 20 channels with the same antibody and fluorophore in all 18 batches. Dropped:
# CD45, CD14 (gating), AF, APC (GITR -> CCR4 at batch05), Alexa 488 / 700 (HELIOS and CCR7
# swap at batch08), BV570 CD27, BV785 CD28, PerCP CD163/CD28 (absent in later batches).
MARKER_CHANNELS <- c(
  "APC-Fire 750-A", "Alexa Fluor 647-A", "BUV395-A", "BUV496-A", "BUV563-A", "BUV661-A",
  "BUV737-A", "BUV805-A", "BV421-A", "BV480-A", "BV605-A", "BV650-A", "BV711-A", "BV750-A",
  "PE-A", "PE-Cy5-A", "PE-Cy7-A", "PE-Dazzle594-A", "Pacific Blue-A", "PerCP-eFluor 710-A"
)
COFACTOR <- 6000 # Q4: spectral
ANCHOR <- "batch01" # Q8: all batches hold 4 references, tie goes to the earliest
FLOWSOM_XDIM <- 10 # N5
FLOWSOM_YDIM <- 10 # N5
FLOWSOM_NCLUS <- 20 # N5
FLOWSOM_SEED <- 1 # N6

data_dir <- file.path(Sys.getenv("SW_DATA_DIR", "/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata"), "cytobatchflagr_zenodo")
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
if (dir.exists(run_dir)) stop(sprintf("Run dir %s exists; use a new one.", run_dir), call. = FALSE)
dir.create(run_dir, recursive = TRUE)

# Markers ------------------------------------------------------------------------
# Labels from the last batch (its panel matches Table_1.pdf); channels are the identifiers.
samples <- read_csv(file.path(data_dir, "Table_2.csv"), show_col_types = FALSE)
par <- pData(parameters(read.FCS(file.path(data_dir, last(samples$file)), which.lines = 1, truncate_max_range = FALSE)))
markers <- tibble(channel = MARKER_CHANNELS, label = par$desc[match(MARKER_CHANNELS, par$name)])
write_csv(markers, file.path(run_dir, "markers.csv"))

# Read, subsample, transform -----------------------------------------------------
# One seeded draw per file, in Table_2 order, so the subsample reproduces.
set.seed(SUBSAMPLE_SEED)
read_one <- function(file) {
  ff <- read.FCS(file.path(data_dir, file), truncate_max_range = FALSE, transformation = FALSE)
  x <- exprs(ff)[, markers$channel, drop = FALSE]
  keep <- if (nrow(x) > CELLS_PER_FILE) sort(sample.int(nrow(x), CELLS_PER_FILE)) else seq_len(nrow(x))
  list(x = asinh(x[keep, , drop = FALSE] / COFACTOR), n_read = nrow(x))
}
read <- map(samples$file, read_one)
samples <- samples |>
  mutate(n_read = map_int(read, "n_read"), n_kept = map_int(read, \(r) nrow(r$x)))
cells <- map(read, "x") |>
  map(as_tibble) |>
  list_rbind(names_to = "i") |>
  mutate(sample_id = samples$sample_id[i], batch = samples$batch[i], control = samples$control[i], .before = 1) |>
  select(-i)
write_csv(samples, file.path(run_dir, "samples.csv"))
saveRDS(cells, file.path(run_dir, "cells.rds"))

# Frozen populations per fold ------------------------------------------------------
# FlowSOM on the anchor batch's three reference controls, then every cell mapped to it (Q5).
xm <- as.matrix(cells[, markers$channel])
fit_fold <- function(held_out) {
  anchor_ref <- cells$batch == ANCHOR & cells$control != held_out
  fsom <- FlowSOM(xm[anchor_ref, ], colsToUse = markers$channel, xdim = FLOWSOM_XDIM, ydim = FLOWSOM_YDIM,
                  nClus = FLOWSOM_NCLUS, seed = FLOWSOM_SEED, silent = TRUE)
  as.integer(GetMetaclusters(NewData(fsom, xm)))
}
controls <- sort(unique(samples$control))
populations <- map(controls, fit_fold) |>
  set_names(paste0("fold_", controls)) |>
  as_tibble()
saveRDS(populations, file.path(run_dir, "populations.rds"))

counts <- populations |>
  mutate(sample_id = cells$sample_id) |>
  tidyr::pivot_longer(-sample_id, names_to = "fold", values_to = "metacluster") |>
  count(fold, sample_id, metacluster)
write_csv(counts, file.path(run_dir, "population_counts.csv"))

writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
