# Purpose: show how many shared live events a multi-channel negative filter would remove, for a grid of
#          thresholds, so the user can pick T and N. Decision support only: nothing is filtered here.
#          An event is flagged when it is below -T (arcsinh) in at least N fluorescence channels, in any arm.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/11_cleanup_grid.R <run_dir>
#          from the repo root.
# Inputs:  <run_dir>/manifest.csv, shared_events.rds [01_manifest-gate.R]
#          the FCS files listed in manifest.csv
# Outputs: <run_dir>/cleanup_grid.csv  per sample, T, N: % of shared live events flagged in any arm,
#                                      and per arm % flagged
#          <run_dir>/sessionInfo_11.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
GRID_T <- c(1, 2, 3) # thresholds shown, arcsinh units; display grid only, the user picks T
GRID_N <- c(3, 5, 8) # channel counts shown; display grid only, the user picks N

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_file <- file.path(run_dir, "cleanup_grid.csv")
stopifnot("11 outputs must not exist" = !file.exists(out_file))
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
shared <- readRDS(file.path(run_dir, "shared_events.rds"))

# Per arm and sample: number of channels below -T for each shared event ---------------------
n_below <- manifest |>
  mutate(
    counts = map2(file, sample, \(f, s) {
      e <- flowCore::exprs(sw_io_read_fcs(f))[shared[[s]], , drop = FALSE]
      ch <- colnames(e)[grepl("-A$", colnames(e)) & !grepl("^(FSC|SSC)|Zombie|LIVE DEAD|^AF", colnames(e))]
      a <- asinh(e[, ch, drop = FALSE] / COFACTOR)
      set_names(GRID_T) |> map(\(t) rowSums(a < -t))
    })
  )

grid <- expand_grid(sample = unique(manifest$sample), t = GRID_T, n = GRID_N) |>
  mutate(
    flags = pmap(list(sample, t, n), \(s, t, n) {
      rows <- which(n_below$sample == s)
      set_names(map(rows, \(i) n_below$counts[[i]][[as.character(t)]] >= n), n_below$arm[rows])
    }),
    pct_any_arm = map_dbl(flags, \(f) 100 * mean(reduce(f, `|`))),
    pct_per_arm = map(flags, \(f) round(100 * map_dbl(f, mean), 2))
  )

grid |>
  select(sample, t, n, pct_any_arm, pct_per_arm) |>
  unnest_wider(pct_per_arm) |>
  write_csv(out_file)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo_11.txt"))
