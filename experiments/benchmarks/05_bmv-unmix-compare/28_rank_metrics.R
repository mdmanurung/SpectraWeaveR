# Purpose: guardrails and tie-break metrics for the revised decision rule (user, 2026-11-10), on the
#          same paired events for every arm. Negatives and positives come from the reference arm's cut
#          only, so adding arms never changes earlier numbers.
#          Guardrails: pile-up index (MAD-based rSD / outer-quantile spread of the negatives; a spike at
#          zero shrinks MAD, not the outer quantiles), TNF-CD8 streak (Spearman rho among events negative
#          for both), hypernegative fraction (below -3 x reference rSD). Tie-breaks: outer-quantile spread
#          (2.5-97.5 / 3.92) and stain index built on it.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/28_rank_metrics.R <src_run_dir> <new_run_dir>
#          from the repo root. <new_run_dir> must not exist.
# Inputs:  <src_run_dir>/manifest.csv, shared_events.rds [01_manifest-gate.R]
#          <src_run_dir>/mindensity_cuts.csv  reference-arm cuts [02_metrics.R]
#          experiments/bmv_test/compare_3D/data/unmix_cells/<EXTRA_ARMS>/*.fcs [25_joint_sweep.R, 22_fix_my_unmix.R]
# Outputs: <new_run_dir>/manifest.csv       all arms scored
#          <new_run_dir>/reference_cuts.csv per sample, channel: reference cut and pinned flag
#          <new_run_dir>/rank_channel.csv   per arm, sample, channel: pinned, n_neg, rsd_mad, spread_q, pileup_index,
#                                           pct_hyperneg, stain_index_q (marker channels)
#          <new_run_dir>/rank_streak.csv    per arm, sample: TNF-CD8 rho among double negatives, n
#          <new_run_dir>/sessionInfo_28.txt

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
REFERENCE_ARM <- "as188_joint" # negatives and positives from this arm's cut only (user, 2026-11-10)
QUANTILES <- c(0.025, 0.975) # outer quantiles (user, 2026-11-09)
Z_WIDTH <- 3.92 # width of the 2.5-97.5 interval in SD units under normality
HYPERNEG_K <- 3 # hypernegative: below -K x reference rSD (user, grill Q19, 2026-11-09)
PIN_TOL <- 0.02 # a cut within this of the window edge counts as pinned (report rule, 2026-11-09)
CUT_WINDOW <- c(0, 3) # mindensity window used for the cuts in 02 (user, 2026-11-09)
STREAK <- c(TNF = "BV510-A", CD8 = "SBUV445-A") # pair from the stripe bisect (2026-11-10)
EXTRA_ARMS <- c(
  "repro188",
  "jpr_false",
  "coll03",
  "coll07",
  "alpha025",
  "alpha075"
) # 25_joint_sweep.R
EXTRA_DIR <- "experiments/bmv_test/compare_3D/data/unmix_cells"

args <- commandArgs(trailingOnly = TRUE)
src_dir <- args[[1]]
new_dir <- args[[2]]
stopifnot("new run dir must not exist" = !dir.exists(new_dir))
dir.create(new_dir, recursive = TRUE)

# Manifest: original arms plus extra arms (same raw files, same row order) ---------------------
manifest_src <- read_csv(
  file.path(src_dir, "manifest.csv"),
  show_col_types = FALSE
)
manifest_extra <- EXTRA_ARMS |>
  map(\(a) {
    tibble(
      arm = a,
      dir = file.path("unmix_cells", a),
      viability = "Zombie NIR-A",
      file = list.files(
        file.path(EXTRA_DIR, a),
        pattern = "\\.fcs$",
        full.names = TRUE
      )
    )
  }) |>
  list_rbind() |>
  mutate(sample = stringr::str_extract(basename(file), "^[A-Z][0-9]+"))
manifest <- bind_rows(manifest_src, manifest_extra)
stopifnot("5 files per arm" = all(count(manifest, arm)$n == 5))
write_csv(manifest, file.path(new_dir, "manifest.csv"))

shared <- readRDS(file.path(src_dir, "shared_events.rds"))
ref_cuts <- read_csv(
  file.path(src_dir, "mindensity_cuts.csv"),
  show_col_types = FALSE
) |>
  filter(arm == REFERENCE_ARM)
# Pinned or missing reference cuts: excluded per sample x channel, not whole channels (user, grill Q21,
# 2026-11-09: "exclude pinned, show both"); the streak check keeps raw cuts (it needs no valley).
ref_cuts <- ref_cuts |>
  mutate(
    pinned = is.na(cut) |
      abs(cut - CUT_WINDOW[[1]]) <= PIN_TOL |
      abs(cut - CUT_WINDOW[[2]]) <= PIN_TOL
  )
write_csv(ref_cuts, file.path(new_dir, "reference_cuts.csv"))
channels <- unique(ref_cuts$channel)
stopifnot("no missing reference cut" = !anyNA(ref_cuts$cut))
message(
  "Pinned sample x channel cells excluded from spread/SI: ",
  sum(ref_cuts$pinned),
  " of ",
  nrow(ref_cuts)
)

# One sample: reference negatives and positives, metrics for every arm ------------------------
process_sample <- function(s) {
  cut_s <- ref_cuts |>
    filter(sample == s) |>
    pull(cut, name = channel)
  pinned_s <- ref_cuts |>
    filter(sample == s) |>
    pull(pinned, name = channel)
  lin <- manifest |>
    filter(sample == s) |>
    (\(m) set_names(m$file, m$arm))() |>
    map(\(f) {
      flowCore::exprs(sw_io_read_fcs(f))[shared[[s]], channels, drop = FALSE]
    })
  ref_asinh <- asinh(lin[[REFERENCE_ARM]] / COFACTOR)
  neg <- sweep(ref_asinh, 2, cut_s[channels], "<")
  pos <- sweep(ref_asinh, 2, cut_s[channels], ">")
  ref_rsd <- map_dbl(set_names(channels), \(ch) {
    1.4826 * mad(lin[[REFERENCE_ARM]][neg[, ch], ch], constant = 1)
  })

  per_channel <- lin |>
    imap(\(m, a) {
      map(channels, \(ch) {
        x <- m[neg[, ch], ch]
        q <- quantile(x, QUANTILES, names = FALSE)
        spread_q <- (q[[2]] - q[[1]]) / Z_WIDTH
        rsd_mad <- 1.4826 * mad(x, constant = 1)
        tibble(
          arm = a,
          sample = s,
          channel = ch,
          pinned = pinned_s[[ch]],
          n_neg = length(x),
          n_pos = sum(pos[, ch]),
          rsd_mad = rsd_mad,
          spread_q = spread_q,
          pileup_index = rsd_mad / spread_q,
          pct_hyperneg = 100 * mean(x < -HYPERNEG_K * ref_rsd[[ch]]),
          stain_index_q = (median(m[pos[, ch], ch]) - median(x)) /
            (2 * spread_q)
        )
      }) |>
        list_rbind()
    }) |>
    list_rbind()

  both_neg <- neg[, STREAK[["TNF"]]] & neg[, STREAK[["CD8"]]]
  streak <- lin |>
    imap(\(m, a) {
      tibble(
        arm = a,
        sample = s,
        n_both_neg = sum(both_neg),
        rho_tnf_cd8 = cor(
          m[both_neg, STREAK[["TNF"]]],
          m[both_neg, STREAK[["CD8"]]],
          method = "spearman"
        )
      )
    }) |>
    list_rbind()
  list(per_channel = per_channel, streak = streak)
}

results <- sort(unique(manifest$sample)) |>
  set_names() |>
  map(process_sample)

write_csv(
  results |> map("per_channel") |> list_rbind(),
  file.path(new_dir, "rank_channel.csv")
)
write_csv(
  results |> map("streak") |> list_rbind(),
  file.path(new_dir, "rank_streak.csv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(new_dir, "sessionInfo_28.txt")
)
