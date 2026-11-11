# Purpose: score the held-out controls (G1). Same scorer as 17_score_controls.R, applied to the unmixed TEST
#          halves and the full H12 file from 26_heldout_controls.R. In a single-stain control every channel
#          except the target should read 0, and in a negative control every channel should. Per arm and
#          control, score each off-target channel on the positive events (or all events of a negative
#          control): spread = (97.5th - 2.5th percentile) / 3.92 and bias = |median|, linear scale.
#          Also compare the train-half spectra with the full-data spectra of each code line.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/27_heldout_score.R <run_dir> [held_root]
#          from the repo root. run_dir must not exist. held_root is only given for a smoke test.
# Inputs:  <held_root>/unmixed/<arm>/*.fcs and verify.csv [26_heldout_controls.R]
#          <held_root>/spectra_c188/spectra.rds, spectra_c102/spectra.rds [26_heldout_controls.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/table_spectra/Cells Matching Negative_autospectral_spectra.csv
# Outputs: <run_dir>/heldout_scores.csv         per arm, control, off-target channel: n events, spread_q, bias
#          <run_dir>/heldout_summary.csv        per arm, control: median spread_q and bias, log2 ratio vs REF_ARM
#          <run_dir>/heldout_groups.csv         per arm and control group: median log2 ratio over controls
#          <run_dir>/heldout_spectra_angle.csv  per code line and fluorophore: angle (deg) train vs full-data spectrum
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
POS_QUANTILE <- 0.995 # positive = target above this quantile of the matching negative (user, 2026-11-09)
QUANTILES <- c(0.025, 0.975) # off-target spread from outer quantiles (user, 2026-11-09)
Z_WIDTH <- 3.92 # width of the 2.5-97.5 interval in SD units under normality
NEGATIVES <- c(Beads = "Negative_1 (Beads)", Cells = "Unstained (Cells)") # matching universal negatives
NEGATIVE_CONTROLS <- c(
  "Negative_1 (Beads)",
  "Unstained (Cells)",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1"
)
REF_ARM <- "as188_joint" # reference arm (README, approved choices)
FULL_SPECTRA <- c(
  # full-data spectra of each code line, as used by the full-data arms in 16
  c188 = file.path(
    DATA,
    "autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
  ),
  c102 = file.path(
    DATA,
    "autoSpectral_v1.0.2/Reference Group/table_spectra/Cells Matching Negative_autospectral_spectra.csv"
  )
)

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
held_root <- if (length(commandArgs(trailingOnly = TRUE)) >= 2) {
  commandArgs(trailingOnly = TRUE)[[2]]
} else {
  file.path(DATA, "heldout_controls")
}
stopifnot("run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)

# Control table: arm, file, control name, type, target channel -------------------------------
# control names come from the raw file each output was made from (verify.csv), because 1.0.2 names
# outputs from $FIL (e.g. "LIVE DEAD NIR (Cells) AutoSpectral.fcs" for "Zombie NIR (Cells).fcs")
controls <- list.files(
  file.path(held_root, "unmixed"),
  pattern = "^verify\\.csv$",
  recursive = TRUE,
  full.names = TRUE
) |>
  map(\(v) {
    read_csv(v, show_col_types = FALSE) |> mutate(arm = basename(dirname(v)))
  }) |>
  list_rbind() |>
  transmute(
    file = out,
    arm = arm,
    control = str_remove(basename(raw), "\\.fcs$"),
    type = str_extract(control, "(?<=\\()(Beads|Cells)(?=\\))"),
    fluor = str_remove(control, " \\((Beads|Cells)\\)$"),
    target = if_else(
      control %in% NEGATIVE_CONTROLS,
      NA_character_,
      paste0(fluor, "-A")
    )
  )
stopifnot(
  "same number of controls in every arm" = n_distinct(count(controls, arm)$n) ==
    1
)

# Score one arm -------------------------------------------------------------------------------
score_arm <- function(arm_files) {
  read <- \(f) flowCore::exprs(sw_io_read_fcs(f))
  neg <- map(NEGATIVES, \(n) read(arm_files$file[arm_files$control == n]))
  arm_files |>
    pmap(\(file, arm, control, type, fluor, target) {
      e <- read(file)
      ch <- colnames(e)[
        grepl("-A$", colnames(e)) & !grepl("^(FSC|SSC)|^AF", colnames(e))
      ]
      rows <- if (is.na(target)) {
        rep(TRUE, nrow(e))
      } else {
        e[, target] >
          quantile(neg[[type]][, target], POS_QUANTILE, names = FALSE)
      }
      map(setdiff(ch, target), \(j) {
        x <- e[rows, j]
        q <- quantile(x, QUANTILES, names = FALSE)
        tibble(
          arm = arm,
          control = control,
          type = type,
          target = target,
          channel = j,
          n = length(x),
          spread_q = (q[[2]] - q[[1]]) / Z_WIDTH,
          bias = abs(median(x))
        )
      }) |>
        list_rbind()
    }) |>
    list_rbind()
}

scores <- controls |>
  split(~arm) |>
  map(\(d) score_arm(select(d, file, arm, control, type, fluor, target))) |>
  list_rbind()

# Summary per control, log2 ratio vs the reference arm, then median per control group ---------
summary <- scores |>
  summarise(
    n = first(n),
    spread_q = median(spread_q),
    bias = median(bias),
    .by = c(arm, control, type, target)
  )
ref <- summary |>
  filter(arm == REF_ARM) |>
  select(control, ref_spread_q = spread_q, ref_bias = bias)
summary <- summary |>
  left_join(ref, by = "control", relationship = "many-to-one") |>
  mutate(
    group = case_when(
      control %in% NEGATIVE_CONTROLS ~ control,
      type == "Cells" ~ "cell single-stain",
      type == "Beads" ~ "bead single-stain"
    ),
    log2_spread = log2(spread_q / ref_spread_q),
    log2_bias = log2(bias / ref_bias)
  )
groups <- summary |>
  summarise(
    n_controls = n(),
    median_log2_spread = median(log2_spread),
    median_log2_bias = median(log2_bias),
    n_nonfinite_bias = sum(!is.finite(log2_bias)), # bias exactly 0 or NA gives an infinite ratio
    .by = c(arm, group)
  )

# Train-half vs full-data spectra: angle per fluorophore -------------------------------------
angle_deg <- function(a, b) {
  acos(pmin(1, sum(a * b) / sqrt(sum(a^2) * sum(b^2)))) * 180 / pi
}
angles <- imap(FULL_SPECTRA, \(full_file, line) {
  full <- as.matrix(read.csv(full_file, row.names = 1, check.names = FALSE))
  train <- readRDS(file.path(
    held_root,
    paste0("spectra_", line),
    "spectra.rds"
  ))
  fl <- setdiff(intersect(rownames(full), rownames(train)), "AF")
  tibble(
    line = line,
    fluorophore = fl,
    angle = map_dbl(fl, \(f) {
      angle_deg(full[f, colnames(train)], train[f, colnames(train)])
    })
  )
}) |>
  list_rbind()

write_csv(scores, file.path(run_dir, "heldout_scores.csv"))
write_csv(summary, file.path(run_dir, "heldout_summary.csv"))
write_csv(groups, file.path(run_dir, "heldout_groups.csv"))
write_csv(angles, file.path(run_dir, "heldout_spectra_angle.csv"))
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
