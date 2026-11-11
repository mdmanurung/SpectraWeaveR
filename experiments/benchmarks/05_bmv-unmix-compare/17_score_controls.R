# Purpose: ground-truth unmixing error on the controls. In a single-stain control every channel except
#          the target should read 0, and in an unstained or negative control every channel should.
#          Per arm and control, score each off-target channel on the positive events (or all events of a
#          negative control): spread = (97.5th - 2.5th percentile) / 3.92 and bias = |median|, linear scale.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/17_score_controls.R <run_dir>
#          from the repo root.
# Inputs:  experiments/bmv_test/compare_3D/data/unmix_controls/<arm>/*.fcs [16_unmix_controls.R]
# Outputs: <run_dir>/control_scores.csv   per arm, control, off-target channel: n events, spread_q, bias
#          <run_dir>/control_summary.csv  per arm, control: median spread_q and bias over off-target channels
#          <run_dir>/sessionInfo_17.txt

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
CONTROL_DIR <- "experiments/bmv_test/compare_3D/data/unmix_controls"
POS_QUANTILE <- 0.995 # positive = target above this quantile of the matching negative (user, 2026-11-09)
QUANTILES <- c(0.025, 0.975) # off-target spread from outer quantiles (user, 2026-11-09)
Z_WIDTH <- 3.92 # width of the 2.5-97.5 interval in SD units under normality
NEGATIVES <- c(Beads = "Negative_1 (Beads)", Cells = "Unstained (Cells)") # matching universal negatives
NEGATIVE_CONTROLS <- c(
  "Negative_1 (Beads)",
  "Unstained (Cells)",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1"
)

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_files <- file.path(run_dir, c("control_scores.csv", "control_summary.csv"))
stopifnot("17 outputs must not exist" = !any(file.exists(out_files)))

# Control table: arm, file, control name, type, target channel -------------------------------
# control names come from the raw file each output was made from (verify.csv), because 1.0.2 names
# outputs from $FIL (e.g. "LIVE DEAD NIR (Cells) AutoSpectral.fcs" for "Zombie NIR (Cells).fcs")
controls <- list.files(
  CONTROL_DIR,
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
stopifnot("36 controls per arm" = all(count(controls, arm)$n == 36))

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

summary <- scores |>
  summarise(
    n = first(n),
    spread_q = median(spread_q),
    bias = median(bias),
    .by = c(arm, control, type, target)
  )

write_csv(scores, out_files[[1]])
write_csv(summary, out_files[[2]])
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_17.txt")
)
