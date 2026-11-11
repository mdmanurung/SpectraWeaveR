# Purpose: rate of impossible double positives on the stained samples, per arm and donor:
#          CD4+CD8+ in CD3+, TCRgd+ in CD3-, TCRgd+CD4+ in CD3+. False double positives come from
#          spillover errors, so lower is better. Positive = above the 99.5th percentile of the same
#          channel in that arm's unmixed H12 unstained sample; SpectroFlo has no H12 export, so its cut is
#          the 99.5th percentile of its own shared negatives in each sample (not directly comparable).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/18_exclusivity.R <run_dir>
#          from the repo root.
# Inputs:  <run_dir>/manifest.csv, shared_events.rds [01_manifest-gate.R]
#          <run_dir>/mindensity_cuts.csv [02_metrics.R]  (shared negatives, SpectroFlo cut only)
#          experiments/bmv_test/compare_3D/data/unmix_controls/<arm>/H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs
#          [16_unmix_controls.R]
# Outputs: <run_dir>/exclusivity.csv       per arm, sample: n parent and % double positive per pair
#          <run_dir>/exclusivity_cuts.csv  per arm, sample, marker: cut (linear) and its source
#          <run_dir>/sessionInfo_18.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor; only used to apply the saved mindensity cuts
POS_QUANTILE <- 0.995 # positive = above this quantile of the unstained (user, 2026-11-09)
MARKERS <- c(
  CD3 = "SBUV510-A",
  CD4 = "Spark Blue 574-A",
  CD8 = "SBUV445-A",
  TCRgd = "BUV395-A"
) # as in 03
H12 <- "experiments/bmv_test/compare_3D/data/unmix_controls/%s/H12 Unstained_BMV_Aurora_3dCD328_Batch1 AutoSpectral.fcs"
NO_H12_ARM <- "spectroflo"

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_files <- file.path(run_dir, c("exclusivity.csv", "exclusivity_cuts.csv"))
stopifnot("18 outputs must not exist" = !any(file.exists(out_files)))
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
shared <- readRDS(file.path(run_dir, "shared_events.rds"))
mcuts <- read_csv(
  file.path(run_dir, "mindensity_cuts.csv"),
  show_col_types = FALSE
)
arms <- unique(manifest$arm)

# Cuts: H12 per arm; SpectroFlo from its shared negatives per sample ---------------------------
h12_cuts <- setdiff(arms, NO_H12_ARM) |>
  set_names() |>
  map(\(a) {
    e <- flowCore::exprs(sw_io_read_fcs(sprintf(H12, a)))
    map_dbl(MARKERS, \(ch) quantile(e[, ch], POS_QUANTILE, names = FALSE))
  })

# One arm x sample: double-positive rates ---------------------------------------------------
excl_one <- function(a, s, f) {
  e <- flowCore::exprs(sw_io_read_fcs(f))[shared[[s]], MARKERS, drop = FALSE]
  colnames(e) <- names(MARKERS)
  cut <- if (a == NO_H12_ARM) {
    # shared negatives of this sample: below the mindensity cut in every arm (as 02_metrics.R)
    below_all <- map(arms, \(b) {
      fb <- manifest |> filter(arm == b, sample == s) |> pull(file)
      cb <- mcuts |> filter(arm == b, sample == s) |> pull(cut, name = channel)
      sweep(
        asinh(
          flowCore::exprs(sw_io_read_fcs(fb))[
            shared[[s]],
            MARKERS,
            drop = FALSE
          ] /
            COFACTOR
        ),
        2,
        cb[MARKERS],
        "<"
      )
    }) |>
      reduce(`&`)
    map2_dbl(names(MARKERS), seq_along(MARKERS), \(m, j) {
      quantile(e[below_all[, j] %in% TRUE, m], POS_QUANTILE, names = FALSE)
    })
  } else {
    h12_cuts[[a]]
  }
  names(cut) <- names(MARKERS)
  pos <- sweep(e, 2, cut[colnames(e)], ">")
  cd3p <- pos[, "CD3"]
  list(
    rates = tibble(
      arm = a,
      sample = s,
      n_cd3pos = sum(cd3p),
      n_cd3neg = sum(!cd3p),
      pct_cd4cd8_dp = 100 * mean(pos[cd3p, "CD4"] & pos[cd3p, "CD8"]),
      pct_tcrgd_in_cd3neg = 100 * mean(pos[!cd3p, "TCRgd"]),
      pct_tcrgd_cd4_dp = 100 * mean(pos[cd3p, "TCRgd"] & pos[cd3p, "CD4"])
    ),
    cuts = tibble(
      arm = a,
      sample = s,
      marker = names(cut),
      cut = cut,
      source = if (a == NO_H12_ARM) "own shared negatives" else "H12 unstained"
    )
  )
}

res <- manifest |>
  select(arm, sample, file) |>
  pmap(\(arm, sample, file) excl_one(arm, sample, file))

write_csv(res |> map("rates") |> list_rbind(), out_files[[1]])
write_csv(res |> map("cuts") |> list_rbind(), out_files[[2]])
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_18.txt")
)
