# Purpose: negative spread and stain index from outer quantiles, which a pile-up of events near zero
#          barely moves. Spread = (97.5th - 2.5th percentile of the shared negatives) / 3.92 (SD under
#          normality); stain index = (median pos - median neg) / (2 x spread). Same events as 02.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/15_quantile_spread.R <run_dir>
#          from the repo root.
# Inputs:  <run_dir>/manifest.csv, shared_events.rds [01_manifest-gate.R]
#          <run_dir>/mindensity_cuts.csv  per arm, sample, channel cut [02_metrics.R]
#          <run_dir>/stain_index.csv      marker, channel and positive rule per arm and sample [02_metrics.R]
# Outputs: <run_dir>/quantile_spread.csv  per arm, sample, channel: n_neg, median, q_lo, q_hi, spread_q
#          <run_dir>/quantile_stain_index.csv  per arm, sample, marker: stain_index_q
#          <run_dir>/sessionInfo_15.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor; only used to apply the saved cuts (CLAUDE.md, confirmed by user)
QUANTILES <- c(0.025, 0.975) # outer quantiles (user, 2026-11-09)
Z_WIDTH <- 3.92 # width of the 2.5-97.5 interval in SD units under normality
REFERENCE_ARM <- "as188_joint" # positives fall back to this arm's cut, as in 02 (user, grill Q3)
MIN_SHARED_POS <- 100 # same fallback rule as 02_metrics.R

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_files <- file.path(
  run_dir,
  c("quantile_spread.csv", "quantile_stain_index.csv")
)
stopifnot("15 outputs must not exist" = !any(file.exists(out_files)))
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
shared <- readRDS(file.path(run_dir, "shared_events.rds"))
cuts <- read_csv(
  file.path(run_dir, "mindensity_cuts.csv"),
  show_col_types = FALSE
)
si_markers <- read_csv(
  file.path(run_dir, "stain_index.csv"),
  show_col_types = FALSE
) |>
  distinct(marker, channel)
arms <- unique(manifest$arm)
channels <- unique(cuts$channel)

# One sample: shared negatives and shared positives per channel, all arms ---------------------
process_sample <- function(s) {
  lin <- arms |>
    set_names() |>
    map(\(a) {
      f <- manifest |> filter(arm == a, sample == s) |> pull(file)
      flowCore::exprs(sw_io_read_fcs(f))[shared[[s]], channels, drop = FALSE]
    })
  above <- arms |>
    set_names() |>
    map(\(a) {
      cut_a <- cuts |>
        filter(arm == a, sample == s) |>
        pull(cut, name = channel)
      m <- sweep(asinh(lin[[a]] / COFACTOR), 2, cut_a[channels], ">")
      m[is.na(m)] <- FALSE
      m
    })
  neg <- above |> map(\(m) !m) |> reduce(`&`)
  pos_all <- above |> reduce(`&`)

  spread <- arms |>
    map(\(a) {
      map(channels, \(ch) {
        x <- lin[[a]][neg[, ch], ch]
        q <- quantile(x, QUANTILES, names = FALSE)
        tibble(
          arm = a,
          sample = s,
          channel = ch,
          n_neg = length(x),
          median_neg = median(x),
          q_lo = q[[1]],
          q_hi = q[[2]],
          spread_q = (q[[2]] - q[[1]]) / Z_WIDTH
        )
      }) |>
        list_rbind()
    }) |>
    list_rbind()

  # positives: above the cut in all arms; reference arm's cut when fewer than MIN_SHARED_POS
  si <- si_markers |>
    pmap(\(marker, channel) {
      pos <- if (sum(pos_all[, channel]) >= MIN_SHARED_POS) {
        pos_all[, channel]
      } else {
        above[[REFERENCE_ARM]][, channel]
      }
      map(arms, \(a) {
        sp <- spread |> filter(arm == a, channel == .env$channel)
        tibble(
          arm = a,
          sample = s,
          marker = marker,
          channel = channel,
          n_pos = sum(pos),
          stain_index_q = (median(lin[[a]][pos, channel]) - sp$median_neg) /
            (2 * sp$spread_q)
        )
      }) |>
        list_rbind()
    }) |>
    list_rbind()
  list(spread = spread, si = si)
}

results <- sort(unique(manifest$sample)) |>
  set_names() |>
  map(process_sample)

write_csv(results |> map("spread") |> list_rbind(), out_files[[1]])
write_csv(results |> map("si") |> list_rbind(), out_files[[2]])
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_15.txt")
)
