# Purpose: measure the BV510 (TNF) part of the Alexa Fluor 405 (GLUT1) and SBUV510 (CD3) spectral variants in
#          every saved spectra set: the 1.0.2 and 1.8.8 zips and the rebuild-chain arms (49_ablation_chain.R),
#          to see at which preprocessing step the BV510 part appears (user, 2026-11-10). Metric as in
#          60_af405_proof_objects.R: projection coefficient of BV510 (1.8.8 spectrum) in a difference spectrum.
#          spread = variant - the set's own reference spectrum for that dye;
#          offset = the set's reference - the 1.8.8 zip reference (one value per set and dye).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/74_variant_bv510_by_set.R   (repo root)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/  table_spectra (run 143803),
#            figure_spectral_variants/Spectral_variants.rds  [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/  table_spectra (Cells
#            Matching Negative), figure_spectral_variants/Spectral_variants.rds  [user zip]
#          experiments/bmv_test/compare_3D/data/unmix_cells/<chain arm>/objects/spectra.rds,
#            spectral_variants.rds  [49_ablation_chain.R]
# Outputs: experiments/runs/2026-11-10_bmv-variant-bv510-by-set/variant_coef.csv   per set, dye, variant: spread coef
#          experiments/runs/2026-11-10_bmv-variant-bv510-by-set/set_summary.csv     per set and dye: n, median and
#                                                                                  90th percentile spread, offset
#          experiments/runs/2026-11-10_bmv-variant-bv510-by-set/fig_variant_bv510_by_set.png (+ _source.csv)
#          experiments/runs/2026-11-10_bmv-variant-bv510-by-set/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ---------------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
DYES <- c("Alexa Fluor 405", "SBUV510")
TARGET <- "BV510"
# sets in preprocessing order: zips, then the chain from the 1.8.8 setup towards the 1.0.2 route
CHAIN <- c(
  "chain1_var",
  "chain2_filt",
  "chain3_af",
  "chain4_spec",
  "chain5_bead",
  "chain6a_neg",
  "chain7a_pregate",
  "chain8_legacy"
)
OUT_DIR <- "experiments/runs/2026-11-10_bmv-variant-bv510-by-set"

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Spectra and variants per set ---------------------------------------------------------------------
read_csv_matrix <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
s188 <- read_csv_matrix(file.path(
  DATA,
  "autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
))
sets <- c(
  list(
    zip_188 = list(
      spectra = s188,
      variants = readRDS(file.path(
        DATA,
        "autoSpectral_v1.8.8/figure_spectral_variants/Spectral_variants.rds"
      ))
    ),
    zip_102 = list(
      spectra = read_csv_matrix(file.path(
        DATA,
        "autoSpectral_v1.0.2/Reference Group/table_spectra/Cells Matching Negative_autospectral_spectra.csv"
      )),
      variants = readRDS(file.path(
        DATA,
        "autoSpectral_v1.0.2/Reference Group/figure_spectral_variants/Spectral_variants.rds"
      ))
    )
  ),
  map(set_names(CHAIN), \(a) {
    d <- file.path(DATA, "unmix_cells", a, "objects")
    list(
      spectra = readRDS(file.path(d, "spectra.rds")),
      variants = readRDS(file.path(d, "spectral_variants.rds"))
    )
  })
)

# BV510 coefficient of each variant's spread, and of each set's reference offset ---------------------
b510 <- s188[TARGET, ]
proj_coef <- \(d, b) as.vector(d %*% b) / sum(b^2) # coefficient of b in each row of d
variant_coef <- imap(sets, \(x, set) {
  map(set_names(DYES), \(dye) {
    v <- x$variants$variants[[dye]][, names(b510)]
    ref <- x$spectra[dye, names(b510)]
    tibble(
      variant = seq_len(nrow(v)),
      spread_coef = proj_coef(sweep(v, 2, ref), b510),
      offset_coef = sum((ref - s188[dye, names(b510)]) * b510) / sum(b510^2)
    )
  }) |>
    list_rbind(names_to = "dye")
}) |>
  list_rbind(names_to = "set") |>
  mutate(set = factor(set, levels = names(sets)))
write_csv(variant_coef, file.path(OUT_DIR, "variant_coef.csv"))

set_summary <- variant_coef |>
  summarise(
    n_variants = n(),
    median_spread = median(spread_coef),
    p90_spread = quantile(spread_coef, 0.9),
    offset = first(offset_coef),
    .by = c(set, dye)
  )
write_csv(set_summary, file.path(OUT_DIR, "set_summary.csv"))

# Figure: every variant as a point, the set's offset as a cross ---------------------------------------
plot_df <- variant_coef
offset_df <- set_summary
p <- ggplot(plot_df, aes(set, spread_coef)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(
    position = position_jitter(width = 0.2, height = 0, seed = 1),
    size = 0.6,
    alpha = 0.5
  ) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.5) +
  geom_point(
    data = offset_df,
    aes(y = offset),
    shape = 4,
    size = 3,
    colour = "red"
  ) +
  facet_wrap(~dye, ncol = 1, scales = "free_y") +
  labs(
    x = NULL,
    y = "BV510 coefficient (projection on the 1.8.8 BV510 spectrum)",
    title = "BV510 part of each variant (points: variant - own reference; red x: reference - 1.8.8 reference)"
  ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(
  file.path(OUT_DIR, "fig_variant_bv510_by_set.png"),
  p,
  width = 7.5,
  height = 6,
  dpi = 150,
  bg = "white"
)
write_csv(plot_df, file.path(OUT_DIR, "fig_variant_bv510_by_set_source.csv"))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
