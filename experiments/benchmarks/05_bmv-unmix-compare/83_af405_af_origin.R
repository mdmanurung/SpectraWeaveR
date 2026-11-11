# Purpose: test 7 (README decision table, user 2026-11-10): is the BV510-like spread in the AutoSpectral
#          1.0.2 Alexa Fluor 405 (GLUT1) variants cell autofluorescence (AF)? Split the AF405 positive control
#          events (clean controls) at the median of their unmixed AF, run the 1.0.2 get.spectral.variants on each
#          half (test 3 settings, 77_af405_variants_102code_clean.R), seeds 42 and 43, and measure the BV510
#          spread (script 74 metric) and the shape of the mean spread against the H12 AF spectra and BV510.
#          Unmixed AF: the zip 1.8.8 setup (spectra 143803, AF set 143933, zip variants, joint, speed slow)
#          applied to the AF405 control file with unmix.folder; AF abundance = its AF-A channel. Positives:
#          unmixed AF405 above the zip 1.8.8 AF405 threshold (variants object slot thresholds).
#          Steps: split (1.8.8 code), variants (1.0.2 code; one half and seed per call), measure (no AutoSpectral).
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/83_af405_af_origin.R <run_dir> split
#          then with ALLOW_EXISTING_RUN=1: <run_dir> variants <low|high> <42|43> (4 calls), then <run_dir> measure.
#          run_dir = experiments/runs/2026-11-10_bmv-af405-af-origin
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/  clean controls +
#            fcs_control_file.csv [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/ spectra 143803, AF 143933, variants [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/table_spectra/
#            'Cells Matching Negative' spectra and H12 AF csv; figure_spectral_variants/Spectral_variants.rds
#          experiments/bmv_test/compare_3D/data/lib_as102  AutoSpectral 1.0.2 library
#          experiments/runs/2026-11-10_bmv-af405-102code-clean/variants_102code_clean.rds  test 3 variants [77]
# Outputs: <run_dir>/af405_unmixed/                 the AF405 control unmixed with the zip 1.8.8 setup (split)
#          <run_dir>/split.csv                      per event kept as positive: AF-A, AF405, half (split)
#          <run_dir>/split_summary.csv              events per half, AF-A median used for the split (split)
#          <run_dir>/controls_low/, controls_high/  control dirs: clean controls with AF405 = that half (split)
#          <run_dir>/variants_<half>_s<seed>.rds    1.0.2-code variants object (variants)
#          <run_dir>/figure_spectral_variants_<half>_s<seed>/  1.0.2 variant figures (variants)
#          <run_dir>/spread.csv                     per set and variant: BV510 coefficient of variant - ref (measure)
#          <run_dir>/spread_summary.csv             per set: n, median, 90th percentile (measure)
#          <run_dir>/shape_cosines.csv              per set: cosine of mean spread to BV510, to the AF mean,
#                                                   max and median to the H12 AF rows (measure)
#          <run_dir>/fig_af_origin_spread.png (+ _source.csv)   BV510 coefficient per set (measure)
#          <run_dir>/fig_af_origin_shape.png (+ _source.csv)    mean spread shape vs AF spectra and BV510
#          <run_dir>/sessionInfo_<step>.txt

args <- commandArgs(trailingOnly = TRUE)
RUN_DIR <- args[[1]]
STEP <- args[[2]]
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
if (STEP == "variants") {
  .libPaths(c(file.path(DATA, "lib_as102"), .libPaths()))
}

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})
if (STEP != "measure") {
  suppressPackageStartupMessages(library(AutoSpectral))
}

# Parameters ---------------------------------------------------------------------------------------
DYE <- "Alexa Fluor 405"
DYE_FILE <- "Alexa Fluor 405 (Cells).fcs"
PIPELINE <- "joint" # zip 1.8.8 setup (inherited: user zip)
SPEED <- "slow" # always slow (user, 2026-11-10)
HALVES <- c("low", "high")
SEEDS <- c(42, 43) # README seed rule plus replicate (user, 2026-11-10, test 7)

# Inputs -------------------------------------------------------------------------------------------
d102 <- file.path(DATA, "autoSpectral_v1.0.2", "Reference Group")
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
clean_dir <- file.path(d188, "Reference Group Clean")
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
} # no rescaling
s102_file <- file.path(
  d102,
  "table_spectra",
  "Cells Matching Negative_autospectral_spectra.csv"
)
af102_file <- file.path(
  d102,
  "table_spectra",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs_H12 Unstained_BMV_Aurora_3dCD328_Batch1_AF.csv"
)
s188_file <- file.path(
  d188,
  "table_spectra",
  "Automated_autospectral_spectra_20261105_143803.csv"
)
af188_file <- file.path(
  d188,
  "table_spectra",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
)
var188_file <- file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)

# Step split: unmixed AF of the AF405 control, positives split at the AF-A median ---------------------
if (STEP == "split") {
  stopifnot(
    "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) ==
      "1.8.8",
    "run dir must not exist" = !dir.exists(RUN_DIR)
  )
  dir.create(RUN_DIR, recursive = TRUE)
  RUN_DIR <- normalizePath(RUN_DIR)
  var188 <- readRDS(var188_file)

  # unmix the AF405 control file alone with the zip 1.8.8 setup
  in_dir <- file.path(RUN_DIR, "af405_raw")
  dir.create(in_dir)
  file.symlink(file.path(clean_dir, DYE_FILE), file.path(in_dir, DYE_FILE))
  unmixed_dir <- file.path(RUN_DIR, "af405_unmixed")
  dir.create(unmixed_dir)
  asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
  flow.control <- reload.flow.control(
    clean_dir,
    file.path(clean_dir, "fcs_control_file.csv"),
    asp
  )
  unmix.folder(
    fcs.dir = in_dir,
    spectra = read_raw(s188_file),
    asp = asp,
    flow.control = flow.control,
    method = "AutoSpectral",
    af.spectra = read_raw(af188_file),
    spectra.variants = var188,
    pipeline = PIPELINE,
    speed = SPEED,
    output.dir = unmixed_dir,
    parallel = TRUE,
    verbose = TRUE
  )

  # raw and unmixed events are the same events in the same order
  raw <- read.FCS(
    file.path(clean_dir, DYE_FILE),
    transformation = FALSE,
    truncate_max_range = FALSE
  )
  # unmix.folder names its output "<input name> AutoSpectral.fcs"
  unm <- exprs(read.FCS(
    file.path(unmixed_dir, sub("\\.fcs$", " AutoSpectral.fcs", DYE_FILE)),
    transformation = FALSE,
    truncate_max_range = FALSE
  ))
  stopifnot("unmixed file keeps every raw event" = nrow(unm) == nrow(raw))
  threshold <- var188$thresholds[[DYE]]
  events <- tibble(
    event = seq_len(nrow(unm)),
    af = unm[, "AF-A"],
    af405 = unm[, paste0(DYE, "-A")]
  )
  positives <- events |>
    filter(af405 > threshold)
  af_median <- median(positives$af)
  positives <- positives |>
    mutate(half = if_else(af <= af_median, "low", "high"))
  write_csv(positives, file.path(RUN_DIR, "split.csv"))
  write_csv(
    positives |>
      summarise(
        n = n(),
        af_median_all_positives = af_median,
        af405_threshold = threshold,
        .by = half
      ) |>
      bind_rows(tibble(
        half = "negative (dropped)",
        n = nrow(events) - nrow(positives)
      )),
    file.path(RUN_DIR, "split_summary.csv")
  )

  # one control dir per half: every clean control linked, AF405 replaced by that half's positive events
  for (h in HALVES) {
    ctrl <- file.path(RUN_DIR, paste0("controls_", h))
    dir.create(ctrl)
    others <- setdiff(list.files(clean_dir), DYE_FILE)
    file.symlink(file.path(clean_dir, others), file.path(ctrl, others))
    keep <- positives$event[positives$half == h]
    write.FCS(raw[keep, ], file.path(ctrl, DYE_FILE))
  }
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_split.txt")
  )
}

# Step variants: 1.0.2 get.spectral.variants on one half, one seed (test 3 settings) -------------------
if (STEP == "variants") {
  half <- args[[3]]
  seed <- as.integer(args[[4]])
  stopifnot(
    "1.0.2 code loaded" = as.character(packageVersion("AutoSpectral")) ==
      "1.0.2",
    "half and seed as planned" = half %in% HALVES && seed %in% SEEDS
  )
  RUN_DIR <- normalizePath(RUN_DIR)
  ctrl <- file.path(RUN_DIR, paste0("controls_", half))
  tag <- paste0(half, "_s", seed)
  out_rds <- file.path(RUN_DIR, paste0("variants_", tag, ".rds"))
  stopifnot(
    "variants for this half and seed not yet made" = !file.exists(out_rds)
  )
  asp <- get.autospectral.param(cytometer = "aurora", figures = TRUE)
  set.seed(seed)
  # the zip call used relative paths from inside control.dir (as in 77)
  old_wd <- setwd(ctrl)
  variants <- get.spectral.variants(
    ctrl,
    "fcs_control_file.csv",
    asp,
    read_raw(s102_file),
    read_raw(af102_file),
    output.dir = file.path(RUN_DIR, paste0("figure_spectral_variants_", tag)),
    parallel = FALSE,
    verbose = TRUE
  )
  setwd(old_wd)
  saveRDS(variants, out_rds)
  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, paste0("sessionInfo_variants_", tag, ".txt"))
  )
}

# Step measure: BV510 spread per set, and the mean spread shape against AF and BV510 ---------------------
if (STEP == "measure") {
  s188 <- read_raw(s188_file)
  s102 <- read_raw(s102_file)
  af102 <- read_raw(af102_file)
  det <- colnames(s188)
  b510 <- s188["BV510", det]
  ref102 <- s102[DYE, det]

  # variant matrices per set; every 1.0.2-code set is measured against the 1.0.2 reference, as in 77
  var_sets <- c(
    list(
      zip_188 = list(
        v = readRDS(var188_file)$variants[[DYE]],
        ref = s188[DYE, det]
      ),
      zip_102 = list(
        v = readRDS(file.path(
          d102,
          "figure_spectral_variants",
          "Spectral_variants.rds"
        ))$variants[[DYE]],
        ref = ref102
      ),
      clean_all_s42 = list(
        v = readRDS(
          "experiments/runs/2026-11-10_bmv-af405-102code-clean/variants_102code_clean.rds"
        )$variants[[DYE]],
        ref = ref102
      )
    ),
    map(
      set_names(paste0(rep(HALVES, each = 2), "_s", SEEDS)),
      \(tag) {
        list(
          v = readRDS(file.path(
            RUN_DIR,
            paste0("variants_", tag, ".rds")
          ))$variants[[DYE]],
          ref = ref102
        )
      }
    )
  )
  set_levels <- names(var_sets)
  deltas <- map(var_sets, \(x) sweep(x$v[, det, drop = FALSE], 2, x$ref))

  # BV510 coefficient of each variant's spread (script 74 metric)
  spread <- imap(deltas, \(d, set) {
    tibble(
      variant = seq_len(nrow(d)),
      spread_coef = as.vector(d %*% b510) / sum(b510^2)
    )
  }) |>
    list_rbind(names_to = "set") |>
    mutate(set = factor(set, levels = set_levels))
  write_csv(spread, file.path(RUN_DIR, "spread.csv"))
  write_csv(
    summarise(
      spread,
      n_variants = n(),
      median = median(spread_coef),
      p90 = quantile(spread_coef, 0.9),
      .by = set
    ),
    file.path(RUN_DIR, "spread_summary.csv")
  )

  # mean spread shape per set, and its cosine to BV510, to the AF mean and to each H12 AF row
  cos <- \(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))
  af_rows <- af102[, det, drop = FALSE]
  af_mean <- colMeans(af_rows)
  mean_shape <- map(deltas, colMeans)
  shape_cosines <- imap(mean_shape, \(m, set) {
    af_cos <- apply(af_rows, 1, \(a) cos(m, a))
    tibble(
      set = set,
      cos_bv510 = cos(m, b510),
      cos_af_mean = cos(m, af_mean),
      cos_af_max = max(af_cos),
      cos_af_median = median(af_cos)
    )
  }) |>
    list_rbind()
  write_csv(shape_cosines, file.path(RUN_DIR, "shape_cosines.csv"))

  # figure 1: BV510 coefficient per set, each variant a point
  p1 <- ggplot(spread, aes(set, spread_coef)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_point(
      position = position_jitter(width = 0.2, height = 0, seed = 1),
      size = 0.7,
      alpha = 0.6
    ) +
    geom_boxplot(outlier.shape = NA, fill = NA, width = 0.5) +
    labs(
      x = NULL,
      y = "BV510 coefficient of (variant - own reference)",
      title = "AF405 variants: 1.0.2 code on the low-AF and high-AF halves of the positives"
    ) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(
    file.path(RUN_DIR, "fig_af_origin_spread.png"),
    p1,
    width = 5.5,
    height = 3.8,
    dpi = 150,
    bg = "white"
  )
  write_csv(spread, file.path(RUN_DIR, "fig_af_origin_spread_source.csv"))

  # figure 2: mean spread shapes (each scaled to max |value| 1) over the AF rows and BV510 (scaled to max 1)
  scale_max <- \(x) x / max(abs(x))
  shape_df <- imap(mean_shape, \(m, set) {
    tibble(detector = det, value = scale_max(m), set = set)
  }) |>
    list_rbind() |>
    mutate(set = factor(set, levels = set_levels))
  af_df <- as_tibble(t(apply(af_rows, 1, scale_max)), rownames = "af") |>
    mutate(af = as.character(seq_len(n()))) |>
    pivot_longer(-af, names_to = "detector", values_to = "value")
  ref_df <- bind_rows(
    tibble(
      detector = det,
      value = scale_max(af_mean),
      line = "AF mean (H12, 1.0.2 set)"
    ),
    tibble(detector = det, value = scale_max(b510), line = "BV510 (1.8.8)")
  )
  plot_det <- \(df) mutate(df, detector = factor(detector, levels = det))
  p2 <- ggplot(plot_det(shape_df), aes(detector, value, group = 1)) +
    geom_line(
      data = plot_det(af_df),
      aes(group = af),
      colour = "grey80",
      linewidth = 0.2
    ) +
    geom_line(
      data = plot_det(ref_df),
      aes(group = line, linetype = line),
      colour = "black",
      linewidth = 0.5
    ) +
    geom_line(colour = "#d95f02", linewidth = 0.5) +
    facet_wrap(~set, ncol = 2) +
    labs(
      x = "Detector",
      y = "Scaled intensity",
      linetype = NULL,
      title = "Mean AF405 variant spread (orange) vs H12 AF spectra (grey) and BV510"
    ) +
    theme_bw(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4),
      legend.position = "top"
    )
  ggsave(
    file.path(RUN_DIR, "fig_af_origin_shape.png"),
    p2,
    width = 9,
    height = 8,
    dpi = 150,
    bg = "white"
  )
  write_csv(
    bind_rows(shape_df, ref_df),
    file.path(RUN_DIR, "fig_af_origin_shape_source.csv")
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(RUN_DIR, "sessionInfo_measure.txt")
  )
}
