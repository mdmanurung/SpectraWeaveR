# Purpose: generalise the script 60 / 74 overlap metric to every dye pair (user, 2026-11-10): how much of
#          channel X's 1.8.8 dye spectrum sits in dye D's variants. Projection coefficient of s188[X, ] in
#            spread  = D's 1.0.2 variant - D's 1.0.2 reference   (median over D's variants)
#            offset  = D's 1.0.2 reference - D's 1.8.8 reference (one value)
#            swapped = D's 1.0.2 variant - D's 1.8.8 reference   (median; what a per-dye arm of 45 uses)
#            zip_188 = D's 1.8.8 variant - D's 1.8.8 reference   (median; contrast)
#          The predicted driver of channel X is the dye D != X with the largest median `swapped` coefficient.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/80_variant_overlap_matrix.R   (repo root)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra (run 143803),
#            figure_spectral_variants/Spectral_variants.rds  [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/table_spectra (Cells Matching
#            Negative), figure_spectral_variants/Spectral_variants.rds  [user zip]
# Outputs: experiments/runs/2026-11-10_bmv-variant-overlap-matrix/overlap.csv            per dye D, channel X, part
#          experiments/runs/2026-11-10_bmv-variant-overlap-matrix/predicted_drivers.csv  per channel X: top dye D
#          experiments/runs/2026-11-10_bmv-variant-overlap-matrix/fig_variant_overlap_matrix.png (+ _source.csv)
#          experiments/runs/2026-11-10_bmv-variant-overlap-matrix/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ---------------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
OUT_DIR <- "experiments/runs/2026-11-10_bmv-variant-overlap-matrix"
PARTS <- c("swapped", "spread", "offset", "zip_188") # panel order

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Spectra and variants ------------------------------------------------------------------------------
read_csv_matrix <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
s188 <- read_csv_matrix(file.path(
  DATA,
  "autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
))
s102 <- read_csv_matrix(file.path(
  DATA,
  "autoSpectral_v1.0.2/Reference Group/table_spectra/Cells Matching Negative_autospectral_spectra.csv"
))[rownames(s188), colnames(s188)]
var188 <- readRDS(file.path(
  DATA,
  "autoSpectral_v1.8.8/figure_spectral_variants/Spectral_variants.rds"
))$variants
var102 <- readRDS(file.path(
  DATA,
  "autoSpectral_v1.0.2/Reference Group/figure_spectral_variants/Spectral_variants.rds"
))$variants
dyes <- rownames(s188)

# Coefficient of every channel's spectrum in every dye's difference spectra ----------------------------
# coef[i, X] = (d_i . s188[X, ]) / |s188[X, ]|^2, one column per channel X
coef_matrix <- \(d) {
  (d %*% t(s188)) / matrix(rowSums(s188^2), nrow(d), length(dyes), byrow = TRUE)
}
overlap <- dyes |>
  set_names() |>
  map(\(dye) {
    parts <- list(
      swapped = coef_matrix(sweep(
        var102[[dye]][, colnames(s188)],
        2,
        s188[dye, ]
      )),
      spread = coef_matrix(sweep(
        var102[[dye]][, colnames(s188)],
        2,
        s102[dye, ]
      )),
      offset = coef_matrix(t(s102[dye, ] - s188[dye, ])),
      zip_188 = coef_matrix(sweep(
        var188[[dye]][, colnames(s188)],
        2,
        s188[dye, ]
      ))
    )
    imap(parts, \(m, part) {
      tibble(
        channel_dye = dyes,
        coef = apply(m, 2, median),
        n_variants = nrow(m)
      )
    }) |>
      list_rbind(names_to = "part")
  }) |>
  list_rbind(names_to = "variant_dye")
write_csv(overlap, file.path(OUT_DIR, "overlap.csv"))

predicted <- overlap |>
  filter(part == "swapped", variant_dye != channel_dye) |>
  slice_max(coef, n = 1, by = channel_dye, with_ties = FALSE) |>
  arrange(desc(coef))
write_csv(predicted, file.path(OUT_DIR, "predicted_drivers.csv"))

# Figure: variant dye x channel dye, one panel per part ---------------------------------------------------
plot_df <- overlap |>
  mutate(
    part = factor(part, levels = PARTS),
    variant_dye = factor(variant_dye, levels = rev(dyes)),
    channel_dye = factor(channel_dye, levels = dyes),
    coef_shown = if_else(variant_dye == channel_dye, NA_real_, coef) # own channel not shown
  )
lim <- max(abs(plot_df$coef_shown), na.rm = TRUE)
p <- ggplot(plot_df, aes(channel_dye, variant_dye, fill = coef_shown)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-lim, lim),
    na.value = "grey85",
    name = "coefficient"
  ) +
  facet_wrap(~part, nrow = 2) +
  labs(
    x = "Channel X (its 1.8.8 dye spectrum)",
    y = "Dye D whose variants or reference are tested",
    title = "Coefficient of channel X's spectrum in dye D's difference spectra (median over variants; grey = own channel)"
  ) +
  theme_bw(base_size = 6) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(
  file.path(OUT_DIR, "fig_variant_overlap_matrix.png"),
  p,
  width = 11,
  height = 10,
  dpi = 150,
  bg = "white"
)
write_csv(plot_df, file.path(OUT_DIR, "fig_variant_overlap_matrix_source.csv"))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
