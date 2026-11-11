# Purpose: test for value clamps in the unmixed data. Per arm, sample and fluorescence channel, count
#          shared live events sitting exactly at the channel minimum (a floor clamp, e.g. the IFNg column
#          at about -4.4 arcsinh) or exactly at zero, and the excess of events near zero.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/14_clamp_check.R <run_dir>
#          from the repo root. Uses the uncleaned files listed in <run_dir>/manifest.csv.
# Inputs:  <run_dir>/manifest.csv, shared_events.rds [01_manifest-gate.R]
# Outputs: <run_dir>/clamp_check.csv  per arm, sample, channel: min value (arcsinh), n and % of shared live
#                                     events at the minimum, % exactly 0, zero-window density ratio
#          <run_dir>/sessionInfo_14.txt

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
ZERO_W <- 0.2 # zero window half-width, arcsinh; display only, same w as the zero pile-up check

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_file <- file.path(run_dir, "clamp_check.csv")
stopifnot("14 outputs must not exist" = !file.exists(out_file))
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
shared <- readRDS(file.path(run_dir, "shared_events.rds"))

# One file: per channel, events at the minimum, at zero, and the zero-window ratio ----------
clamp_file <- function(f, s) {
  e <- flowCore::exprs(sw_io_read_fcs(f))[shared[[s]], , drop = FALSE]
  ch <- colnames(e)[grepl("-A$", colnames(e)) & !grepl("^(FSC|SSC)|Zombie|LIVE DEAD|^AF", colnames(e))]
  map(ch, \(j) {
    x <- e[, j]
    a <- asinh(x / COFACTOR)
    near <- sum(abs(a) < ZERO_W)
    neighbours <- (sum(a >= -3 * ZERO_W & a < -ZERO_W) + sum(a > ZERO_W & a <= 3 * ZERO_W)) / 2
    tibble(
      channel = j,
      min_arcsinh = min(a),
      n_at_min = sum(x == min(x)),
      pct_at_min = 100 * mean(x == min(x)),
      pct_exact_zero = 100 * mean(x == 0),
      zero_ratio = near / neighbours
    )
  }) |>
    list_rbind()
}

result <- manifest |>
  mutate(res = map2(file, sample, clamp_file)) |>
  select(arm, sample, res) |>
  tidyr::unnest(res)
write_csv(result, out_file)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo_14.txt"))
