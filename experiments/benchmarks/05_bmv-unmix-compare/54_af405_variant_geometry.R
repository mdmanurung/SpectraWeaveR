# Purpose: describe the Alexa Fluor 405 (GLUT1) control and its spectral variants in AutoSpectral 1.0.2 and
#          1.8.8, to test whether the 1.0.2 variants that squeeze BV510 (TNF) look like autofluorescence (AF)
#          or like BV510. Also defines the variant subsets unmixed by 55_af405_subset_unmix.R.
#          Test 1: brightness of the control and of the matched negative in the AF405 peak and BV510 peak
#          detectors. Test 2: cosine of every variant, and of every variant minus the 1.8.8 AF405
#          spectrum (the delta the 1.8.8 joint unmixer uses), to the AF spectra, BV510 and BV570.
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/54_af405_variant_geometry.R \
#            experiments/runs/2026-11-10_bmv-af405-mechanism
#          (from the repo root; run dir must not exist)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/  controls, spectra,
#            AF, Spectral_variants.rds [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/  Reference Group (raw), Reference Group
#            Clean, spectra 143803, AF 143933, Spectral_variants.rds [user zip]
# Outputs: <run_dir>/control_brightness.csv  per version and file: events, quantiles in AF405 and BV510 peak
#                                            detectors
#          <run_dir>/variant_cosines.csv     per variant: cosines of variant and delta to AF, BV510, BV570
#          <run_dir>/variant_subsets.rds     named list of 1.0.2 AF405 row indices per subset arm
#          <run_dir>/variant_subsets.csv     the same, long format, with the split values
#          <run_dir>/thresholds.csv          AF405 positivity threshold in each version
#          <run_dir>/fig_variant_cosines.png distributions per version (+ _source.csv)
#          <run_dir>/fig_variant_spectra.png 1.0.2 and 1.8.8 AF405 variants with BV510 and mean AF spectra
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
DYE <- "Alexa Fluor 405"
TARGETS <- c("BV510", "BV570") # TNF and IFNg dyes, the channels that moved
SUBSET_SEED <- 42 # seed for the random-26 count control (user-approved seed for variant subsampling)
N_RANDOM <- 26 # 1.8.8 holds 26 variants per dye

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)
theme_set(ggpubr::theme_pubr(base_size = 11))

# Inputs ----------------------------------------------------------------------------------
d102 <- file.path(DATA, "autoSpectral_v1.0.2", "Reference Group")
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
# raw csv read: read.spectra() rescales rows (user, 2026-11-10)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra188 <- read_raw(file.path(
  d188,
  "table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
))
spectra102 <- read_raw(file.path(
  d102,
  "table_spectra/Cells Matching Negative_autospectral_spectra.csv"
))
af188 <- read_raw(file.path(
  d188,
  "table_spectra/H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
))
af102 <- read_raw(file.path(
  d102,
  "table_spectra/H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv"
))
v188 <- readRDS(file.path(
  d188,
  "figure_spectral_variants/Spectral_variants.rds"
))
v102 <- readRDS(file.path(
  d102,
  "figure_spectral_variants/Spectral_variants.rds"
))
detectors <- colnames(spectra188)
stopifnot(
  "detector columns identical" = identical(
    colnames(v102$variants[[DYE]]),
    detectors
  ) &&
    identical(colnames(af102), detectors) &&
    identical(colnames(af188), detectors)
)

# Test 1: control brightness in the AF405 and BV510 peak detectors -----------------------------
# peak detectors from the 1.8.8 spectra; quantiles of the raw (unmixed-free) control events
peak <- c(
  AF405 = detectors[which.max(spectra188[DYE, ])],
  BV510 = detectors[which.max(spectra188["BV510", ])]
)
control_files <- tribble(
  ~version , ~set                                      , ~file                                                                   ,
  "1.0.2"  , "Reference Group (raw, used by 1.0.2)"    , file.path(d102, "Alexa Fluor 405 (Cells).fcs")                          ,
  "1.0.2"  , "Reference Group (raw, used by 1.0.2)"    , file.path(d102, "Unstained (Cells).fcs")                                ,
  "1.8.8"  , "Reference Group Clean (rectangle-gated)" , file.path(d188, "Reference Group Clean", "Alexa Fluor 405 (Cells).fcs") ,
  "1.8.8"  , "Reference Group Clean (rectangle-gated)" , file.path(d188, "Reference Group Clean", "Unstained (Cells).fcs")
)
control_brightness <- control_files |>
  mutate(
    values = map(file, \(f) {
      e <- exprs(read.FCS(
        f,
        transformation = FALSE,
        truncate_max_range = FALSE
      ))
      e[, peak, drop = FALSE]
    }),
    n_events = map_int(values, nrow),
    summary = map(values, \(m) {
      imap(peak, \(det, nm) {
        q <- quantile(m[, det], c(0.5, 0.9, 0.99))
        tibble(
          peak_of = nm,
          detector = det,
          q50 = q[[1]],
          q90 = q[[2]],
          q99 = q[[3]]
        )
      }) |>
        list_rbind()
    })
  ) |>
  mutate(file = basename(file)) |>
  select(-values) |>
  unnest(summary)
write_csv(control_brightness, file.path(run_dir, "control_brightness.csv"))

# Test 2: cosines of variants and deltas to AF, BV510, BV570 ------------------------------------
# delta = variant - 1.8.8 AF405 spectrum: what the 1.8.8 joint unmixer subtracts in the swap arm
cos_rows <- \(m, ref) as.vector(m %*% ref) / sqrt(rowSums(m^2) * sum(ref^2))
max_abs_cos <- \(m, refs) {
  apply(
    abs(m %*% t(refs)) / outer(sqrt(rowSums(m^2)), sqrt(rowSums(refs^2))),
    1,
    max
  )
}
ref188 <- spectra188[DYE, ]
variant_sets <- list(
  `1.0.2` = v102$variants[[DYE]],
  `1.8.8` = v188$variants[[DYE]]
)
variant_cosines <- imap(variant_sets, \(m, set_version) {
  delta <- sweep(m, 2, ref188)
  delta_norm <- sqrt(rowSums(delta^2))
  # a variant identical to the reference has no delta direction
  delta_ok <- delta_norm > 1e-9
  ref_own <- if (set_version == "1.0.2") spectra102[DYE, ] else ref188
  tibble(
    version = set_version,
    row = seq_len(nrow(m)),
    variant = rownames(m),
    delta_norm = delta_norm,
    cos_variant_ref188 = cos_rows(m, ref188),
    cos_variant_ref_own = cos_rows(m, ref_own),
    maxcos_variant_af188 = max_abs_cos(m, af188),
    maxcos_variant_af102 = max_abs_cos(m, af102),
    maxcos_delta_af188 = if_else(delta_ok, max_abs_cos(delta, af188), NA_real_),
    maxcos_delta_af102 = if_else(delta_ok, max_abs_cos(delta, af102), NA_real_),
    cos_delta_bv510 = if_else(
      delta_ok,
      cos_rows(delta, spectra188["BV510", ]),
      NA_real_
    ),
    cos_delta_bv570 = if_else(
      delta_ok,
      cos_rows(delta, spectra188["BV570", ]),
      NA_real_
    ),
    delta_bv510_peak = delta[, peak[["BV510"]]]
  )
}) |>
  list_rbind()
write_csv(variant_cosines, file.path(run_dir, "variant_cosines.csv"))
write_csv(
  tibble(
    version = c("1.0.2", "1.8.8"),
    threshold = c(v102$thresholds[[DYE]], v188$thresholds[[DYE]]),
    n_variants = c(nrow(v102$variants[[DYE]]), nrow(v188$variants[[DYE]]))
  ),
  file.path(run_dir, "thresholds.csv")
)

# Subset definitions for 55: median splits of the 1.0.2 rows, plus a random 26 --------------------
c102 <- variant_cosines |>
  filter(version == "1.0.2", !is.na(maxcos_delta_af188))
med_af <- median(c102$maxcos_delta_af188)
med_bv <- median(c102$cos_delta_bv510)
set.seed(SUBSET_SEED)
subsets <- list(
  af405_102_lowAF = c102$row[c102$maxcos_delta_af188 <= med_af],
  af405_102_highAF = c102$row[c102$maxcos_delta_af188 > med_af],
  af405_102_lowBV510 = c102$row[c102$cos_delta_bv510 <= med_bv],
  af405_102_highBV510 = c102$row[c102$cos_delta_bv510 > med_bv],
  af405_102_rand26 = sort(sample(nrow(v102$variants[[DYE]]), N_RANDOM))
)
saveRDS(subsets, file.path(run_dir, "variant_subsets.rds"))
imap(subsets, \(rows, arm) tibble(arm = arm, row = rows)) |>
  list_rbind() |>
  left_join(
    c102 |> select(row, maxcos_delta_af188, cos_delta_bv510),
    by = "row",
    relationship = "many-to-one"
  ) |>
  mutate(median_af = med_af, median_bv510 = med_bv) |>
  write_csv(file.path(run_dir, "variant_subsets.csv"))

# Figure: cosine distributions per version -------------------------------------------------
cos_long <- variant_cosines |>
  select(
    version,
    row,
    maxcos_delta_af188,
    cos_delta_bv510,
    cos_delta_bv570,
    delta_norm
  ) |>
  pivot_longer(-c(version, row), names_to = "measure") |>
  mutate(
    measure = factor(
      measure,
      levels = c(
        "delta_norm",
        "maxcos_delta_af188",
        "cos_delta_bv510",
        "cos_delta_bv570"
      ),
      labels = c(
        "Delta size (norm)",
        "Delta vs AF (max |cos|)",
        "Delta vs BV510 (cos)",
        "Delta vs BV570 (cos)"
      )
    )
  )
write_csv(cos_long, file.path(run_dir, "fig_variant_cosines_source.csv"))
p_cos <- ggplot(cos_long, aes(version, value, colour = version)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0, seed = 1),
    size = 1.2,
    alpha = 0.7
  ) +
  facet_wrap(~measure, scales = "free_y", nrow = 1) +
  scale_colour_manual(
    values = c(`1.0.2` = "#D55E00", `1.8.8` = "#0072B2"),
    guide = "none"
  ) +
  labs(x = "AF405 variant set", y = NULL)
ggsave(
  file.path(run_dir, "fig_variant_cosines.png"),
  p_cos,
  width = 260,
  height = 90,
  units = "mm",
  dpi = 150,
  bg = "white"
)

# Figure: the variant spectra with BV510 and mean AF ------------------------------------------
spec_long <- bind_rows(
  imap(variant_sets, \(m, version) {
    as_tibble(m, rownames = "variant") |>
      mutate(version = version, kind = "AF405 variant")
  }) |>
    list_rbind(),
  tibble(variant = "BV510", version = "both", kind = "BV510") |>
    bind_cols(as_tibble(t(spectra188["BV510", ]))),
  tibble(
    variant = "AF mean (1.8.8)",
    version = "both",
    kind = "AF mean (1.8.8)"
  ) |>
    bind_cols(as_tibble(t(colMeans(af188) / max(colMeans(af188)))))
) |>
  pivot_longer(all_of(detectors), names_to = "detector") |>
  mutate(detector = factor(detector, levels = detectors))
p_spec <- ggplot(
  spec_long,
  aes(detector, value, group = interaction(version, variant))
) +
  geom_line(
    data = \(d) filter(d, kind == "AF405 variant"),
    aes(colour = version),
    alpha = 0.35,
    linewidth = 0.4
  ) +
  geom_line(
    data = \(d) filter(d, kind != "AF405 variant"),
    aes(linetype = kind),
    linewidth = 0.7
  ) +
  scale_colour_manual(
    values = c(`1.0.2` = "#D55E00", `1.8.8` = "#0072B2"),
    name = "AF405 variants"
  ) +
  scale_linetype_manual(
    values = c(BV510 = "solid", `AF mean (1.8.8)` = "dashed"),
    name = NULL
  ) +
  labs(x = "Detector", y = "Normalised intensity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
ggsave(
  file.path(run_dir, "fig_variant_spectra.png"),
  p_spec,
  width = 260,
  height = 110,
  units = "mm",
  dpi = 150,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
