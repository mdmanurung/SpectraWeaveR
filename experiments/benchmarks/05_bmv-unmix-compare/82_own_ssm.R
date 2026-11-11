# Purpose: compute the spillover spread matrix (SSM) from the AutoSpectral 1.8.8 single-stain controls with a
#          stated direction, and use it to find the orientation of the zip's raw SSM (user, 2026-11-10; README
#          row "Test 6, own SSM (82)"). Direction here: row = stained dye i whose positive events spread,
#          column = dye channel j receiving the spread.
#          Unmixing: ordinary least squares with the 33 zip 1.8.8 dye spectra (run 143803), no AF row (the
#          usual SSM setting: spread from fluorochromes only).
#          Populations: positives = events of dye i's stained control above the zip 1.8.8 positivity threshold
#          of dye i (Spectral_variants.rds $thresholds); negatives = that control's universal negative from the
#          zip control file.
#          ss(i -> j) = sqrt(max(0, rSD_j(pos_i)^2 - rSD_j(neg_i)^2)) / sqrt(median_i(pos_i) - median_i(neg_i)),
#          rSD = Nguyen q84 - q50, NA below 100 events (as 34_rank_symmetric.R).
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/82_own_ssm.R experiments/runs/2026-11-10_bmv-own-ssm
#          (from the repo root; run dir must not exist)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/*.fcs,
#            fcs_control_file.csv  [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/
#            Automated_autospectral_spectra_20261105_143803.csv  [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/figure_spectral_variants/Spectral_variants.rds
#            ($thresholds)  [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/figure_similarity_heatmap/
#            autospectral_spillover_spread_matrix.csv  [user zip]
# Outputs: <run_dir>/own_ssm.csv           long: source (row) dye, target (column) dye, ss, population sizes
#          <run_dir>/own_ssm_matrix.csv    the same as a matrix, rows = spreading dye
#          <run_dir>/orientation.csv       correlation of ours with the zip SSM as stored and transposed
#          <run_dir>/fig_own_ssm.png       heatmap, values printed
#          <run_dir>/fig_ssm_vs_zip.png    ours vs zip, as stored and transposed (+ _source.csv)
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})

# Parameters ---------------------------------------------------------------------------------------
ZIP <- "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8"
NGUYEN_P <- c(0.5, 0.84) # Nguyen 2013 robust SD = q84 - q50 (as 34_rank_symmetric.R)
MIN_N <- 100 # NA below this many events (as 34_rank_symmetric.R)
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)

# Spectra, thresholds, control file ------------------------------------------------------------------
spectra <- as.matrix(read.csv(
  file.path(
    ZIP,
    "table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
  ),
  row.names = 1,
  check.names = FALSE
))
thresholds <- readRDS(file.path(
  ZIP,
  "figure_spectral_variants/Spectral_variants.rds"
))$thresholds
controls <- read_csv(
  file.path(ZIP, "Reference Group Clean/fcs_control_file.csv"),
  show_col_types = FALSE
) |>
  filter(fluorophore %in% rownames(spectra))
stopifnot(
  "one control and one threshold per dye" = setequal(
    controls$fluorophore,
    rownames(spectra)
  ) &&
    setequal(names(thresholds), rownames(spectra))
)
dyes <- rownames(spectra)
# OLS projection: abundances = raw %*% t(S) %*% (S t(S))^-1
ols <- t(spectra) %*% solve(spectra %*% t(spectra))

# OLS-unmixed events per control file ---------------------------------------------------------------
unmix_file <- \(f) {
  ff <- read.FCS(
    file.path(ZIP, "Reference Group Clean", f),
    truncate_max_range = FALSE,
    transformation = FALSE
  )
  exprs(ff)[, colnames(spectra)] %*% ols
}
files <- unique(c(controls$filename, controls$universal.negative))
unmixed <- map(set_names(files), unmix_file)

# SSM: spread of every other channel on dye i's positives vs its negatives ----------------------------
rsd <- \(x) {
  if (length(x) < MIN_N) {
    NA_real_
  } else {
    diff(quantile(x, NGUYEN_P, names = FALSE))
  }
}
med <- \(x) if (length(x) < MIN_N) NA_real_ else median(x)
own_ssm <- controls |>
  select(source = fluorophore, filename, universal.negative) |>
  pmap(\(source, filename, universal.negative) {
    pos <- unmixed[[filename]]
    pos <- pos[pos[, source] > thresholds[[source]], , drop = FALSE]
    neg <- unmixed[[universal.negative]]
    den <- med(pos[, source]) - med(neg[, source])
    tibble(
      source = source,
      target = setdiff(dyes, source),
      n_pos = nrow(pos),
      n_neg = nrow(neg),
      ss = map_dbl(target, \(j) {
        if (is.na(den) || den <= 0) {
          return(NA_real_)
        }
        sqrt(max(0, rsd(pos[, j])^2 - rsd(neg[, j])^2)) / sqrt(den)
      })
    )
  }) |>
  list_rbind()
write_csv(own_ssm, file.path(RUN_DIR, "own_ssm.csv"))

own_mat <- matrix(
  NA_real_,
  length(dyes),
  length(dyes),
  dimnames = list(dyes, dyes)
)
own_mat[cbind(own_ssm$source, own_ssm$target)] <- own_ssm$ss
write.csv(own_mat, file.path(RUN_DIR, "own_ssm_matrix.csv"))

# Orientation of the zip SSM: correlate ours with it as stored and transposed --------------------------
zip_mat <- as.matrix(read.csv(
  file.path(
    ZIP,
    "figure_similarity_heatmap/autospectral_spillover_spread_matrix.csv"
  ),
  row.names = 1,
  check.names = FALSE
))[dyes, dyes]
compare_df <- bind_rows(
  as_tibble(
    as.table(own_mat),
    .name_repair = ~ c("source", "target", "ours")
  ) |>
    mutate(
      zip = zip_mat[cbind(source, target)],
      zip_reading = "zip as stored (row = spreading)"
    ),
  as_tibble(
    as.table(own_mat),
    .name_repair = ~ c("source", "target", "ours")
  ) |>
    mutate(
      zip = zip_mat[cbind(target, source)],
      zip_reading = "zip transposed (column = spreading)"
    )
) |>
  filter(source != target, !is.na(ours), !is.na(zip))
orientation <- compare_df |>
  summarise(
    n_pairs = n(),
    spearman = cor(ours, zip, method = "spearman"),
    pearson_log1p = cor(log1p(ours), log1p(zip)),
    .by = zip_reading
  )
write_csv(orientation, file.path(RUN_DIR, "orientation.csv"))

# Figure 1: our SSM heatmap ----------------------------------------------------------------------------
heat_df <- own_ssm |>
  mutate(
    source = factor(source, levels = rev(dyes)),
    target = factor(target, levels = dyes),
    label = if_else(is.na(ss) | ss < 1, "", sprintf("%.0f", ss))
  )
p_heat <- ggplot(heat_df, aes(target, source, fill = ss)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = label), size = 1.8) +
  scale_fill_distiller(
    palette = "YlOrRd",
    direction = 1,
    trans = "sqrt",
    na.value = "grey85",
    name = "SSM"
  ) +
  coord_fixed() +
  labs(
    x = "receiving dye channel (column)",
    y = "stained dye whose positives spread (row)",
    title = "SSM from 1.8.8 single-stain controls (OLS, zip spectra); values >= 1 printed"
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave(
  file.path(RUN_DIR, "fig_own_ssm.png"),
  p_heat,
  width = 8,
  height = 7.5,
  dpi = 170,
  bg = "white"
)

# Figure 2: ours vs zip, both readings of the zip matrix -------------------------------------------------
scatter_labels <- orientation |>
  mutate(
    text = sprintf(
      "Spearman %.2f\nPearson (log1p) %.2f\nn = %d",
      spearman,
      pearson_log1p,
      n_pairs
    )
  )
p_scatter <- ggplot(compare_df, aes(log1p(zip), log1p(ours))) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_text(
    data = scatter_labels,
    aes(x = -Inf, y = Inf, label = text),
    hjust = -0.05,
    vjust = 1.1,
    size = 2.5
  ) +
  facet_wrap(~zip_reading) +
  labs(
    x = "zip SSM, log1p",
    y = "our SSM (row = spreading), log1p",
    title = "Which way does the zip SSM run?"
  ) +
  theme_bw(base_size = 8)
ggsave(
  file.path(RUN_DIR, "fig_ssm_vs_zip.png"),
  p_scatter,
  width = 7,
  height = 3.6,
  dpi = 170,
  bg = "white"
)
write_csv(compare_df, file.path(RUN_DIR, "fig_ssm_vs_zip_source.csv"))

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
