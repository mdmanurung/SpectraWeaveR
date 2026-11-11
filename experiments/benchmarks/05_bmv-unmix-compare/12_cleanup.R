# Purpose: write cleaned FCS files with the same events removed in every arm (paired design kept).
#          An event is dropped when it is a shared live event and, in ANY arm, either
#          (1) PeacoQC RemoveMargins (sw_filter_margins_peacoqc) removes it: any fluorescence channel below
#              -21,761 or above 1e7, or a scatter channel at its instrument limit, or
#          (2) it is below -NEG_T (arcsinh) in at least NEG_N fluorescence channels.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/12_cleanup.R <src_run_dir> <new_run_dir>
#          from the repo root. <new_run_dir> must not exist.
# Inputs:  <src_run_dir>/manifest.csv, shared_events.rds [01_manifest-gate.R]
#          the FCS files listed in manifest.csv; event order is identical across arms (checked in 01)
# Outputs: <new_run_dir>/manifest.csv                 same columns as the source, file points to cleaned FCS
#          <new_run_dir>/cleanup_counts.csv           per sample: events, shared live, dropped by margins, by negatives, total
#          experiments/bmv_test/compare_3D/data/cleaned/<new run dir name>/<arm>/<file>.fcs   cleaned FCS (git-ignored, under data/)
#          <new_run_dir>/sessionInfo_12.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NEG_T <- 2 # below -T in arcsinh units (user, 2026-11-09, from the cleanup_grid.csv options)
NEG_N <- 3 # in at least N fluorescence channels (user, 2026-11-09)
MARGIN_SPECS <- list(AllFluoChannels = c(-21761, 1e7)) # PeacoQC RemoveMargins on every fluorescence
#   channel: minRange -21,761 (arcsinh -2), maxRange 1e7 (user, 2026-11-09); scatter keeps PeacoQC's
#   instrument range, so events at the scatter maximum are removed
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; finds the rows PeacoQC kept
MAX_MARGIN_FRAC <- 0.5 # stop if margins flag more than half the shared live events (defaults did 99.9%)
CLEAN_ROOT <- "experiments/bmv_test/compare_3D/data/cleaned" # files go to CLEAN_ROOT/<new run dir name>/<arm>/

args <- commandArgs(trailingOnly = TRUE)
src_dir <- args[[1]]
new_dir <- args[[2]]
stopifnot("new run dir must not exist" = !dir.exists(new_dir))
dir.create(new_dir, recursive = TRUE)
manifest <- read_csv(file.path(src_dir, "manifest.csv"), show_col_types = FALSE)
shared <- readRDS(file.path(src_dir, "shared_events.rds"))

# Per arm and sample: rows flagged at the scatter maximum and by negatives ------------------------------
flag_file <- function(f, s) {
  ff <- sw_io_read_fcs(f)
  e <- flowCore::exprs(ff)
  key <- do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))
  kept <- sw_filter_margins_peacoqc(
    ff,
    channel_specifications = MARGIN_SPECS
  ) |>
    flowCore::exprs()
  kept_key <- do.call(paste, c(as.data.frame(kept[, KEY_COLS]), sep = "|"))
  ch <- colnames(e)[
    grepl("-A$", colnames(e)) &
      !grepl("^(FSC|SSC)|Zombie|LIVE DEAD|^AF", colnames(e))
  ]
  list(
    margin = !(key %in% kept_key),
    neg = rowSums(asinh(e[, ch, drop = FALSE] / COFACTOR) < -NEG_T) >= NEG_N
  )
}
flags <- manifest |>
  mutate(flag = map2(file, sample, flag_file))

# Union over arms, restricted to shared live events ------------------------------------------
drops <- unique(manifest$sample) |>
  set_names() |>
  map(\(s) {
    f <- flags$flag[flags$sample == s]
    margin <- reduce(map(f, "margin"), `|`) & shared[[s]]
    neg <- reduce(map(f, "neg"), `|`) & shared[[s]]
    list(margin = margin, neg = neg, drop = margin | neg)
  })

counts <- imap(drops, \(d, s) {
  tibble(
    sample = s,
    n_events = length(d$drop),
    n_shared_live = sum(shared[[s]]),
    n_margin = sum(d$margin),
    n_neg = sum(d$neg),
    n_dropped = sum(d$drop)
  )
}) |>
  list_rbind()
write_csv(counts, file.path(new_dir, "cleanup_counts.csv"))
stopifnot(
  "margins flag more than MAX_MARGIN_FRAC of shared live events" = all(
    counts$n_margin / counts$n_shared_live < MAX_MARGIN_FRAC
  )
)

# Write cleaned FCS and the new manifest ---------------------------------------------------
manifest_clean <- manifest |>
  mutate(
    clean_file = file.path(CLEAN_ROOT, basename(new_dir), arm, basename(file)),
    wrote = pmap_lgl(list(file, sample, clean_file), \(f, s, out) {
      dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
      ff <- sw_io_read_fcs(f)
      flowCore::write.FCS(ff[!drops[[s]]$drop, ], out)
      TRUE
    })
  ) |>
  mutate(file = clean_file) |>
  select(arm, dir, viability, file, sample)
write_csv(manifest_clean, file.path(new_dir, "manifest.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(new_dir, "sessionInfo_12.txt")
)
