# Purpose: proof chain step 1 to 4 (user, 2026-11-10): show how the AutoSpectral 1.0.2 Alexa Fluor 405
#          (GLUT1) variants overlap BV510 (TNF), and build the variants objects for the intervention arms.
#          The 1.8.8 joint unmixer (AutoSpectralRcpp 1.3.4, src/unmix_autospectral_joint_pipeline.cpp:223-227)
#          uses delta = variant - spectra[dye, ] with the spectra passed to unmix, so in the swap arm the
#          effective AF405 delta is d = v102 - s188 = e (spread around the 1.0.2 reference, v102 - s102) plus
#          o (reference offset, s102 - s188). Arms change one part of d at a time:
#            ko_bv510   1.0.2 variants with the BV510-parallel part of d removed
#            ki_bv510   1.8.8 variants plus a BV510 part matching the 1.0.2 coefficients by rank
#            dose_000   every variant = s102 (offset only, no spread); dose_025, dose_050 = s102 + a * e
#            recentre   s188 + e (1.0.2 spread, no offset)
#          alpha = 1 is the existing arm var102_only_Alexa_Fluor_405 (45_variant_per_dye.R).
# Run:     sbatch experiments/slurm/run_af405_proof.sbatch 60 (or pixi run Rscript <this script>)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/Automated_..._143803.csv
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/table_spectra/
#            Cells Matching Negative_autospectral_spectra.csv   1.0.2 reference spectra [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/figure_spectral_variants/Spectral_variants.rds
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/figure_spectral_variants/
#            Spectral_variants.rds
# Outputs: experiments/runs/2026-11-10_bmv-af405-proof/objects/<arm>.rds   variants object per arm
#          experiments/runs/2026-11-10_bmv-af405-proof/overlap_coef.csv      per variant: projection and cosine
#          experiments/runs/2026-11-10_bmv-af405-proof/fig1_delta_overlap.png  deltas with BV510 / BV570 shape
#          experiments/runs/2026-11-10_bmv-af405-proof/fig1b_coef.png        coefficient and cosine per version
#          experiments/runs/2026-11-10_bmv-af405-proof/fig1_delta_overlap_source.csv
#          experiments/runs/2026-11-10_bmv-af405-proof/objects_check.csv     per arm: variants, norms, coef
#          experiments/runs/2026-11-10_bmv-af405-proof/sessionInfo_60.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
RUN_DIR <- "experiments/runs/2026-11-10_bmv-af405-proof"
DYE <- "Alexa Fluor 405"
TARGETS <- c("BV510", "BV570") # BV510 = TNF; BV570 = IFNg, the contrast
DOSES <- c(dose_000 = 0, dose_025 = 0.25, dose_050 = 0.5) # user-approved grid (alpha 1 = existing arm)

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(file.path(RUN_DIR, "objects"), recursive = TRUE)

# Inputs ----------------------------------------------------------------------------------
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
s188 <- read_raw(file.path(
  DATA,
  "autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
))
s102 <- read_raw(file.path(
  DATA,
  "autoSpectral_v1.0.2/Reference Group/table_spectra/Cells Matching Negative_autospectral_spectra.csv"
))[rownames(s188), colnames(s188)]
var188 <- readRDS(file.path(
  DATA,
  "autoSpectral_v1.8.8/figure_spectral_variants/Spectral_variants.rds"
))
var102 <- readRDS(file.path(
  DATA,
  "autoSpectral_v1.0.2/Reference Group/figure_spectral_variants/Spectral_variants.rds"
))
v102 <- var102$variants[[DYE]]
v188 <- var188$variants[[DYE]]
stopifnot(
  "detector order equal" = identical(colnames(v102), colnames(s188)) &&
    identical(colnames(v188), colnames(s188))
)
ref188 <- s188[DYE, ]
ref102 <- s102[DYE, ]

# Step 1: effective deltas and their BV510 / BV570 parts ------------------------------------
proj_coef <- \(d, b) as.vector(d %*% b) / sum(b^2) # coefficient of b in each row of d
cosine <- \(d, b) as.vector(d %*% b) / (sqrt(rowSums(d^2)) * sqrt(sum(b^2)))
deltas <- list(
  `1.0.2 (as used in 1.8.8: v - s188)` = sweep(v102, 2, ref188),
  `1.0.2 spread only (v - s102)` = sweep(v102, 2, ref102),
  `1.8.8 (v - s188)` = sweep(v188, 2, ref188)
)
overlap <- imap(deltas, \(d, version) {
  map(set_names(TARGETS), \(t) {
    tibble(
      variant = rownames(d),
      coef = proj_coef(d, s188[t, ]),
      cosine = cosine(d, s188[t, ]),
      delta_norm = sqrt(rowSums(d^2))
    )
  }) |>
    list_rbind(names_to = "target")
}) |>
  list_rbind(names_to = "version")
offset <- ref102 - ref188
offset_row <- tibble(
  version = "offset s102 - s188",
  target = TARGETS,
  variant = "offset",
  coef = map_dbl(TARGETS, \(t) sum(offset * s188[t, ]) / sum(s188[t, ]^2)),
  cosine = map_dbl(TARGETS, \(t) {
    sum(offset * s188[t, ]) / sqrt(sum(offset^2) * sum(s188[t, ]^2))
  }),
  delta_norm = sqrt(sum(offset^2))
)
write_csv(
  bind_rows(overlap, offset_row),
  file.path(RUN_DIR, "overlap_coef.csv")
)

# Figure 1: deltas across detectors with the BV510 and BV570 shapes ----------------------------
detector_order <- colnames(s188) # instrument order: UV, V, B, YG, R lasers by wavelength
delta_long <- imap(deltas[c(1, 3)], \(d, version) {
  as_tibble(d, rownames = "variant") |>
    pivot_longer(-variant, names_to = "detector", values_to = "delta") |>
    mutate(version = version)
}) |>
  list_rbind()
shape_long <- map(set_names(TARGETS), \(t) {
  tibble(
    detector = detector_order,
    shape = s188[t, ] * max(abs(delta_long$delta)) # scaled to the delta range, display only
  )
}) |>
  list_rbind(names_to = "target")
plot_frame <- delta_long |>
  mutate(detector = factor(detector, levels = detector_order))
shape_frame <- shape_long |>
  mutate(detector = factor(detector, levels = detector_order))
write_csv(plot_frame, file.path(RUN_DIR, "fig1_delta_overlap_source.csv"))
p1 <- ggplot(plot_frame, aes(detector, delta, group = variant)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(alpha = 0.35, linewidth = 0.3, colour = "grey30") +
  geom_line(
    data = shape_frame,
    aes(detector, shape, colour = target, group = target),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  scale_colour_manual(
    values = c(BV510 = "#d95f02", BV570 = "#1b9e77"),
    name = "Dye shape"
  ) +
  facet_wrap(~version, ncol = 1) +
  labs(x = "Detector", y = "AF405 variant minus 1.8.8 AF405 spectrum") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
ggsave(
  file.path(RUN_DIR, "fig1_delta_overlap.png"),
  p1,
  width = 12,
  height = 7,
  dpi = 130,
  bg = "white"
)

p1b <- overlap |>
  pivot_longer(c(coef, cosine), names_to = "metric", values_to = "value") |>
  ggplot(aes(version, value)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40") +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 42),
    size = 1,
    alpha = 0.6
  ) +
  facet_grid(metric ~ target, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(
  file.path(RUN_DIR, "fig1b_coef.png"),
  p1b,
  width = 10,
  height = 6,
  dpi = 130,
  bg = "white"
)

# Steps 2 to 4: variants objects for the intervention arms -------------------------------------
# Rows are renormalised to max 1, as AutoSpectral stores variants. delta.list feeds only the
# covariance of the leakage term (cpp:243-246); it is kept on the swap arm's convention (v - s102)
# for 1.0.2-based arms and (v - s188) for 1.8.8-based arms. delta.norms only gate the "no variation"
# check (sanitize.optimization.inputs), so dose_000 keeps the 1.0.2 norms to stay optimised.
norm_rows <- \(m) m / apply(m, 1, max)
b510 <- s188["BV510", ]
d_eff <- sweep(v102, 2, ref188) # what the swap arm feeds the unmixer
c102 <- proj_coef(d_eff, b510)
c188 <- proj_coef(sweep(v188, 2, ref188), b510)
c_ki <- quantile(
  c102,
  probs = (rank(c188, ties.method = "first") - 0.5) / length(c188),
  names = FALSE
)
spread102 <- sweep(v102, 2, ref102)
new_variants <- list(
  ko_bv510 = norm_rows(v102 - outer(c102, b510)), # s188 + d - c b
  ki_bv510 = norm_rows(v188 + outer(c_ki, b510)),
  recentre = norm_rows(sweep(spread102, 2, ref188, `+`))
)
new_variants <- c(
  new_variants,
  map(DOSES, \(a) norm_rows(sweep(a * spread102, 2, ref102, `+`)))
)
delta_ref <- c(
  ko_bv510 = "102",
  ki_bv510 = "188",
  recentre = "102",
  dose_000 = "102",
  dose_025 = "102",
  dose_050 = "102"
)
build_object <- \(v, arm) {
  dimnames(v) <- dimnames(if (nrow(v) == nrow(v102)) v102 else v188)
  d <- sweep(v, 2, if (delta_ref[[arm]] == "102") ref102 else ref188)
  norms <- sqrt(rowSums(d^2))
  if (all(norms < 1e-12)) {
    norms <- var102$delta.norms[[DYE]]
  } # dose_000: keep AF405 optimised
  obj <- var188
  obj$variants[[DYE]] <- v
  obj$delta.list[[DYE]] <- d
  obj$delta.norms[[DYE]] <- norms
  obj
}
objects <- imap(new_variants, build_object)
iwalk(objects, \(obj, arm) {
  saveRDS(obj, file.path(RUN_DIR, "objects", paste0(arm, ".rds")))
})
check <- imap(objects, \(obj, arm) {
  v <- obj$variants[[DYE]]
  tibble(
    arm = arm,
    n_variants = nrow(v),
    row_max_min = min(apply(v, 1, max)),
    median_coef_bv510 = median(proj_coef(sweep(v, 2, ref188), b510)),
    median_eff_delta_norm = median(sqrt(rowSums(sweep(v, 2, ref188)^2))),
    threshold = obj$thresholds[[DYE]],
    other_dyes_equal_188 = identical(
      obj$variants[names(obj$variants) != DYE],
      var188$variants[names(var188$variants) != DYE]
    )
  )
}) |>
  list_rbind()
check <- bind_rows(
  tibble(
    arm = "var102_only_Alexa_Fluor_405 (alpha 1)",
    n_variants = nrow(v102),
    row_max_min = min(apply(v102, 1, max)),
    median_coef_bv510 = median(c102),
    median_eff_delta_norm = median(sqrt(rowSums(d_eff^2))),
    threshold = var188$thresholds[[DYE]],
    other_dyes_equal_188 = TRUE
  ),
  tibble(
    arm = "as188_joint",
    n_variants = nrow(v188),
    row_max_min = min(apply(v188, 1, max)),
    median_coef_bv510 = median(c188),
    median_eff_delta_norm = median(sqrt(rowSums(sweep(v188, 2, ref188)^2))),
    threshold = var188$thresholds[[DYE]],
    other_dyes_equal_188 = TRUE
  ),
  check
)
write_csv(check, file.path(RUN_DIR, "objects_check.csv"))
stopifnot("other dyes untouched in every arm" = all(check$other_dyes_equal_188))
writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_60.txt")
)
