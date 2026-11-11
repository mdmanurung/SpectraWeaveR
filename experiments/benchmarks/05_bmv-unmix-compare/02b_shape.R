# Purpose: shape of each channel's negative population (centred on zero, symmetric, round) and the
#          fraction of strongly negative events, per arm, on the shared live singlet events.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/02b_shape.R <run_dir> [out_dir]
#          from the repo root. Reads the run dir written by 01 and 02; writes to out_dir (default run_dir).
# Inputs:  <run_dir>/manifest.csv, shared_events.rds [01_manifest-gate.R]
#          <run_dir>/mindensity_cuts.csv: per arm, sample, channel cut on arcsinh [02_metrics.R]
#          the FCS files listed in manifest.csv
# Outputs: <out_dir>/shape_negatives.csv  per arm, sample, channel: centring, symmetry on the shared negatives
#          <out_dir>/negative_cells.csv   per arm, sample, channel: fraction of events below -K x reference rSD
#          <out_dir>/roundness.csv        per arm, sample, channel pair: axis ratio of the standardised negative cloud
#          <out_dir>/sessionInfo_02b.txt

# Bioconductor first, tidyverse parts last, so dplyr::filter() is not masked.
suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor for spectral flow (CLAUDE.md, confirmed by user); only used to apply the saved cuts
REFERENCE_ARM <- "as188_joint" # its negative rSD sets the common negative-cell threshold (user, 2026-11-09)
NEG_CELL_K <- 3 # negative cell: below -K x reference rSD of that channel (user, grill Q19, 2026-11-09)
SYMMETRY_Q <- c(0.1587, 0.8413) # median -/+ one SD under normality; symmetry = log2(lower half width / upper half width)
WINSOR_RSD <- 4 # roundness: values clipped to median +/- this many rSD before the correlation (limits outlier leverage)

# Run dirs --------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[[1]]
out_dir <- if (length(args) >= 2) args[[2]] else run_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_files <- file.path(
  out_dir,
  c("shape_negatives.csv", "negative_cells.csv", "roundness.csv")
)
stopifnot("02b outputs must not exist" = !any(file.exists(out_files)))

# Inputs ----------------------------------------------------------------------------------
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
shared <- readRDS(file.path(run_dir, "shared_events.rds"))
cuts <- read_csv(
  file.path(run_dir, "mindensity_cuts.csv"),
  show_col_types = FALSE
)
arms <- unique(manifest$arm)
samples <- sort(unique(manifest$sample))
markers <- unique(cuts$channel)
rsd <- \(x) 1.4826 * mad(x, constant = 1) # mad() already scales by 1.4826 unless constant = 1

# One sample: all arms, shared events, shared negatives ------------------------------------
process_sample <- function(s) {
  lin <- arms |>
    set_names() |>
    map(\(a) {
      f <- manifest |> filter(arm == a, sample == s) |> pull(file)
      flowCore::exprs(sw_io_read_fcs(f))[shared[[s]], markers, drop = FALSE]
    })
  # shared negative set per channel: below the arm's own cut in every arm (as in 02_metrics.R)
  neg <- arms |>
    map(\(a) {
      cut_a <- cuts |>
        filter(arm == a, sample == s) |>
        pull(cut, name = channel)
      sweep(asinh(lin[[a]] / COFACTOR), 2, cut_a[markers], "<") |>
        (\(m) {
          m[is.na(m)] <- FALSE
          m
        })()
    }) |>
    reduce(`&`)
  ref_rsd <- map_dbl(markers |> set_names(), \(ch) {
    rsd(lin[[REFERENCE_ARM]][neg[, ch], ch])
  })

  shape <- arms |>
    map(\(a) {
      markers |>
        map(\(ch) {
          x <- lin[[a]][neg[, ch], ch]
          q <- quantile(
            x,
            c(SYMMETRY_Q[[1]], 0.5, SYMMETRY_Q[[2]]),
            names = FALSE
          )
          tibble(
            arm = a,
            sample = s,
            channel = ch,
            n_neg = length(x),
            centring = median(x) / rsd(x),
            symmetry_log2 = log2((q[[2]] - q[[1]]) / (q[[3]] - q[[2]]))
          )
        }) |>
        list_rbind()
    }) |>
    list_rbind()

  neg_cells <- arms |>
    map(\(a) {
      tibble(
        arm = a,
        sample = s,
        channel = markers,
        n_events = nrow(lin[[a]]),
        threshold = -NEG_CELL_K * ref_rsd[markers],
        n_below = map_int(seq_along(markers), \(i) {
          sum(lin[[a]][, markers[[i]]] < threshold[[i]])
        })
      ) |>
        mutate(frac_below = n_below / n_events)
    }) |>
    list_rbind()

  # roundness: standardised, winsorised negative values; pairs use events negative for both
  pair_idx <- combn(length(markers), 2)
  roundness <- arms |>
    map(\(a) {
      z <- lin[[a]]
      for (j in seq_along(markers)) {
        v <- z[, j]
        v[!neg[, j]] <- NA
        ctr <- median(v, na.rm = TRUE)
        sc <- rsd(v[!is.na(v)])
        z[, j] <- pmin(pmax((v - ctr) / sc, -WINSOR_RSD), WINSOR_RSD)
      }
      r <- suppressWarnings(cor(z, use = "pairwise.complete.obs"))
      n_both <- crossprod(!is.na(z))
      tibble(
        arm = a,
        sample = s,
        channel_a = markers[pair_idx[1, ]],
        channel_b = markers[pair_idx[2, ]],
        n_both = n_both[cbind(pair_idx[1, ], pair_idx[2, ])],
        r_wins = r[cbind(pair_idx[1, ], pair_idx[2, ])]
      ) |>
        mutate(axis_ratio = sqrt((1 - abs(r_wins)) / (1 + abs(r_wins))))
    }) |>
    list_rbind()
  list(shape = shape, neg_cells = neg_cells, roundness = roundness)
}

results <- samples |> set_names() |> map(process_sample)

# Write -----------------------------------------------------------------------------------
write_csv(results |> map("shape") |> list_rbind(), out_files[[1]])
write_csv(results |> map("neg_cells") |> list_rbind(), out_files[[2]])
write_csv(results |> map("roundness") |> list_rbind(), out_files[[3]])
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo_02b.txt")
)
