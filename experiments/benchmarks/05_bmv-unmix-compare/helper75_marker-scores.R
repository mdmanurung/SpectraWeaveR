# Purpose: shared by 75_all_marker_screen.R and 76_per_dye_channel_map.R. Loads an arm's events for one donor
#          in as188_joint's cleaned row order and scores every channel's negatives with the 39b / 51 metrics.
# Run:     source("experiments/benchmarks/05_bmv-unmix-compare/helper75_marker-scores.R") after the caller
#          has set COFACTOR, CUT and SQUEEZE_FRAC.
# Inputs:  experiments/runs/2026-11-10_bmv-gated-pooled/manifest.csv   cleaned FCS per pooled arm [12_cleanup.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed   raw as188_joint
#          experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/   raw unmixed arms [37, 45, 48 scripts]
# Outputs: none (functions only)

POOLED_MANIFEST <- "experiments/runs/2026-11-10_bmv-gated-pooled/manifest.csv"
RAW_REF_DIR <- "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed"
UNMIX_DIR <- "experiments/bmv_test/compare_3D/data/unmix_cells"
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; maps raw rows to cleaned rows

read_x <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
}
key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))

# One arm and donor, all channels, linear, rows = as188_joint's cleaned rows. Pooled arms are read from
# their cleaned files; other arms are raw unmixed files masked to as188_joint's cleaned rows (as 51 does).
arm_matrix <- function(arm, sample, channels, pooled) {
  ref_clean <- pooled$file[
    pooled$arm == "as188_joint" & pooled$sample == sample
  ]
  ref_key <- key_of(read_x(ref_clean))
  if (arm %in% pooled$arm) {
    e <- read_x(pooled$file[pooled$arm == arm & pooled$sample == sample])
  } else {
    raw <- read_x(file.path(UNMIX_DIR, arm, basename(ref_clean)))
    raw_ref_key <- key_of(read_x(file.path(RAW_REF_DIR, basename(ref_clean))))
    if (nrow(raw) != length(raw_ref_key)) {
      stop(
        sprintf("%s %s: raw rows differ from as188_joint raw", arm, sample),
        call. = FALSE
      )
    }
    e <- raw[raw_ref_key %in% ref_key, , drop = FALSE]
  }
  if (!identical(key_of(e), ref_key)) {
    stop(
      sprintf("%s %s: rows do not pair with as188_joint", arm, sample),
      call. = FALSE
    )
  }
  e[, channels, drop = FALSE]
}

# 39b / 51 readouts per channel on given rows: Nguyen rSD (q84 - q50) of events below arcsinh CUT, its ratio
# to as188_joint on the same rows, share with |x| < SQUEEZE_FRAC x as188_joint's rSD, share above arcsinh CUT
nguyen <- \(v) unname(diff(quantile(v, c(0.5, 0.84))))
score_channels <- function(x, x_ref, rows) {
  colnames(x) |>
    set_names() |>
    map(\(ch) {
      lin <- x[rows, ch]
      lin_ref <- x_ref[rows, ch]
      ref_rsd <- nguyen(lin_ref[asinh(lin_ref / COFACTOR) < CUT])
      tibble(
        n = length(rows),
        rsd_neg = nguyen(lin[asinh(lin / COFACTOR) < CUT]),
        ref_rsd = ref_rsd,
        rsd_ratio = rsd_neg / ref_rsd,
        squeezed_share = mean(abs(lin) < SQUEEZE_FRAC * ref_rsd),
        pos_share = mean(asinh(lin / COFACTOR) > CUT)
      )
    }) |>
    list_rbind(names_to = "channel")
}
