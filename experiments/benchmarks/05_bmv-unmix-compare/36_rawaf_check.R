# Purpose: does loading the AF csv raw (not row-rescaled by read.spectra) change the s188_c102 arm?
#          Verifies s188_c102_rawaf against s188_143803_c102, compares per donor and channel,
#          checks the AUTOFLUORESCENCE keyword against the raw csv, and plots BV510-A.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/36_rawaf_check.R
#          from the repo root, after 06b_unmix_c102_rawaf.R has written s188_c102_rawaf.
# Inputs:  data/unmix_cells/s188_143803_c102          [06_unmix_c102.R, read.spectra AF]
#          data/unmix_cells/s188_c102_rawaf           [06b_unmix_c102_rawaf.R, raw AF]
#          experiments/runs/2026-11-10_bmv-gated/manifest.csv   [01_manifest-gate.R; as188_joint]
#          data/autoSpectral_v1.8.8/table_spectra/<AF csv>      [user zip]
# Outputs: data/unmix_cells/s188_c102_rawaf/verify.csv       event count and channel names per file
#          <RUN>/comparison.csv                              per donor and channel differences
#          <RUN>/bv_tail.csv                                 BV510-A and BV570-A noise and positive share
#          <RUN>/af_keyword.csv                              keyword vs raw csv vs read.spectra csv
#          <RUN>/bv510_density.png, bv510_density_source.csv, sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
RUN <- "experiments/runs/2026-11-10_bmv-rawaf-check" # new dir, never overwritten
COFACTOR <- 6000 # arcsinh cofactor for spectral flow (CLAUDE.md)
NEG_CUT <- 1 # arcsinh units: below = negative population for the Nguyen SD
DIFF_TOL <- 1 # linear units: an event "differs" when |diff| exceeds this (float noise is far below)
BV_CHANNELS <- c("BV510-A", "BV570-A")

# Setup -----------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(ggplot2)
})
stopifnot("run dir must not exist" = !dir.exists(RUN))
dir.create(RUN, recursive = TRUE)
theme_set(ggpubr::theme_pubr(base_size = 8))
data_dir <- "experiments/bmv_test/compare_3D/data"
ref_dir <- file.path(data_dir, "unmix_cells/s188_143803_c102")
raw_dir <- file.path(data_dir, "unmix_cells/s188_c102_rawaf")
read_fcs <- \(f) read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
sample_of <- \(f) str_extract(basename(f), "^[A-Z]+[0-9]+")

files <- tibble(
  ref = list.files(ref_dir, pattern = "\\.fcs$", full.names = TRUE)
) |>
  mutate(sample = sample_of(ref), raw = file.path(raw_dir, basename(ref)))

# Verify ----------------------------------------------------------------------------------
# Event count and channel names must equal the reference per file, or stop.
frames <- files |>
  mutate(
    ref_ff = map(ref, read_fcs),
    raw_ff = map(raw, read_fcs)
  )
verify <- frames |>
  reframe(
    sample,
    n_ref = map_int(ref_ff, nrow),
    n_raw = map_int(raw_ff, nrow),
    channels_equal = map2_lgl(ref_ff, raw_ff, \(a, b) {
      identical(colnames(a), colnames(b))
    })
  ) |>
  mutate(ok = n_ref == n_raw & channels_equal)
write_csv(verify, file.path(raw_dir, "verify.csv"))
if (!all(verify$ok)) {
  stop("verify.csv has mismatches; see ", file.path(raw_dir, "verify.csv"))
}

# Per donor and channel differences -------------------------------------------------------
comparison <- frames |>
  reframe(
    sample,
    ch = map2(ref_ff, raw_ff, \(a, b) {
      x <- exprs(a)
      y <- exprs(b)[, colnames(x)]
      d <- abs(x - y)
      tibble(
        channel = colnames(x),
        frac_ne = colMeans(d != 0),
        frac_gt_tol = colMeans(d > DIFF_TOL),
        median_abs = apply(d, 2, median),
        max_abs = apply(d, 2, max)
      )
    })
  ) |>
  unnest(ch)
write_csv(comparison, file.path(RUN, "comparison.csv"))

# BV510 and BV570 tails -------------------------------------------------------------------
# Nguyen robust SD (q84 - q50) of arcsinh values below NEG_CUT; share above NEG_CUT.
nguyen_sd <- \(a) {
  unname(quantile(a[a < NEG_CUT], 0.84) - quantile(a[a < NEG_CUT], 0.50))
}
arm_values <- \(ff, ch) asinh(exprs(ff)[, ch] / COFACTOR)
bv_tail <- frames |>
  select(sample, ref_ff, raw_ff) |>
  pivot_longer(c(ref_ff, raw_ff), names_to = "arm", values_to = "ff") |>
  mutate(
    arm = recode(arm, ref_ff = "s188_c102", raw_ff = "s188_c102_rawaf")
  ) |>
  cross_join(tibble(channel = BV_CHANNELS)) |>
  mutate(a = map2(ff, channel, arm_values)) |>
  reframe(
    nguyen_sd = map_dbl(a, nguyen_sd),
    share_above = map_dbl(a, \(v) mean(v >= NEG_CUT)),
    .by = c(sample, arm, channel)
  )
write_csv(bv_tail, file.path(RUN, "bv_tail.csv"))

# AUTOFLUORESCENCE keyword ----------------------------------------------------------------
# Keyword layout: nrow, ncol, row names, col names, then the values row-major.
parse_kw <- \(s) {
  tok <- str_split_1(s, ",")
  nr <- as.integer(tok[[1]])
  nc <- as.integer(tok[[2]])
  vals <- as.numeric(tok[-seq_len(2 + nr + nc)])
  m <- matrix(
    vals,
    nr,
    nc,
    byrow = TRUE,
    dimnames = list(tok[3:(2 + nr)], tok[(3 + nr):(2 + nr + nc)])
  )
  m
}
af_csv <- file.path(
  data_dir,
  "autoSpectral_v1.8.8/table_spectra/H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
)
af_raw <- as.matrix(read.csv(af_csv, row.names = 1, check.names = FALSE))
af_rescaled <- af_raw / apply(af_raw, 1, max)
af_keyword <- frames |>
  pivot_longer(c(ref_ff, raw_ff), names_to = "arm", values_to = "ff") |>
  reframe(
    sample,
    arm = recode(arm, ref_ff = "s188_c102", raw_ff = "s188_c102_rawaf"),
    kw = map(ff, \(f) {
      parse_kw(keyword(f)[["$AUTOFLUORESCENCE"]])[, colnames(af_raw)]
    }),
    max_abs_vs_raw = map_dbl(kw, \(m) max(abs(m - af_raw))),
    max_abs_vs_rescaled = map_dbl(kw, \(m) max(abs(m - af_rescaled)))
  ) |>
  select(!kw)
write_csv(af_keyword, file.path(RUN, "af_keyword.csv"))

# BV510-A density figure ------------------------------------------------------------------
manifest <- read_csv(
  "experiments/runs/2026-11-10_bmv-gated/manifest.csv",
  show_col_types = FALSE
)
joint <- manifest |>
  filter(arm == "as188_joint") |>
  reframe(sample, joint_ff = map(file, read_fcs))
dens_source <- frames |>
  select(sample, s188_c102 = ref_ff, s188_c102_rawaf = raw_ff) |>
  inner_join(
    joint |> rename(as188_joint = joint_ff),
    by = "sample",
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  pivot_longer(!sample, names_to = "arm", values_to = "ff") |>
  reframe(sample, arm, bv510_asinh = map(ff, \(f) arm_values(f, "BV510-A"))) |>
  unnest(bv510_asinh) |>
  mutate(
    arm = factor(arm, levels = c("as188_joint", "s188_c102", "s188_c102_rawaf"))
  )
write_csv(dens_source, file.path(RUN, "bv510_density_source.csv"))

arm_cols <- c(
  as188_joint = "black",
  s188_c102 = "#D55E00",
  s188_c102_rawaf = "#0072B2"
)
p <- ggplot(dens_source, aes(bv510_asinh, colour = arm, linetype = arm)) +
  geom_density(linewidth = 0.4) +
  scale_colour_manual(values = arm_cols) +
  scale_linetype_manual(
    values = c(
      as188_joint = "solid",
      s188_c102 = "dashed",
      s188_c102_rawaf = "dotted"
    )
  ) +
  facet_wrap(vars(sample), nrow = 1) +
  labs(
    x = "BV510-A, arcsinh(x / 6000)",
    y = "Event density",
    colour = NULL,
    linetype = NULL
  )
ggsave(
  file.path(RUN, "bv510_density.png"),
  plot = p,
  width = 183,
  height = 60,
  units = "mm",
  dpi = 300
)

writeLines(
  capture.output(sessionInfo()),
  file.path(RUN, "sessionInfo.txt")
)
