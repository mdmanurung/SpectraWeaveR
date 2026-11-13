# Purpose: refit mindensity_cuts.csv for 34 on the common shared cells (plans/common_gates_plan.md, audit note).
#          Same fit as 02_metrics.R: per arm, donor and marker channel, gate_mindensity on arcsinh values of
#          the shared events, window CUT_RANGE; a cut at the data minimum means no valley and is NA. Only the
#          cuts are refitted, because 34 reads nothing else from 02. The 02 quality and spread tables are not
#          rebuilt here.
# Run:     SW_DATASET=3d pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/02c_metrics_common.R
#          from the repo root, after 01c_shared_events_common.R; on SLURM:
#          sbatch experiments/slurm/run_common_gs.sbatch 02c_metrics_common.R
#          SW_DATASET is 3d, 24h or exv (no default). SW_RUN_DATE is not used: the cuts go into the 01c dir.
# Inputs:  registry keys of datasets/<SW_DATASET>_runs.csv:
#          run_in("common-final")/manifest.csv, shared_events.rds   raw files and SpectroFlo live on raw rows [01c]
#          run_in("unmix-compare-final")/mindensity_cuts.csv   old cuts, for comparison [02]; the comparison
#                                  file is left out, with a message, when 02 has not run on the dataset
# Outputs: <RUN_DIR> = run_in("common-final"), the 01c dir (it must be registered before this script runs)
#          <RUN_DIR>/mindensity_cuts.csv         arm, sample, channel, cut, found (format of 02)
#          <RUN_DIR>/mindensity_cuts_vs_old.csv  per arm, sample, channel: new and old cut
#          <RUN_DIR>/sessionInfo_02c.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(openCyto)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE) # sw_io_read_fcs, sw_channel_is_fluor
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- run_in("common-final") # written by 01c
OLD_CUTS <- file.path(run_in("unmix-compare-final"), "mindensity_cuts.csv")
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
CUT_RANGE <- c(0, 3) # mindensity window on arcsinh (user, 2026-11-09), as 02_metrics.R
REFERENCE_ARM <- "as188_joint" # as 02: its channels define the marker set
N_CORES <- 5 # one worker per sample

stopifnot(
  "cuts do not exist" = !file.exists(file.path(RUN_DIR, "mindensity_cuts.csv"))
)
manifest <- read_csv(file.path(RUN_DIR, "manifest.csv"), show_col_types = FALSE)
shared <- readRDS(file.path(RUN_DIR, "shared_events.rds"))
samples <- names(shared)
stopifnot(
  "shared events cover the samples" = setequal(samples, samples_tbl()$sample)
)
arms <- unique(manifest$arm)

# One donor: cuts for every arm and marker channel -------------------------------------------------
fit_cuts <- function(s) {
  arm_art <- map(set_names(arms), \(a) {
    row <- dplyr::filter(manifest, arm == a, sample == s)
    ff <- sw_io_read_fcs(row$file)
    stopifnot(
      "shared mask length matches event count" = length(shared[[s]]) == nrow(ff)
    )
    ff <- ff[shared[[s]], ]
    # marker set as 02: fluorescence channels minus AF-A, AF Index and this arm's viability channel
    chans <- flowCore::colnames(ff)[sw_channel_is_fluor(ff)] |>
      setdiff(c("AF-A", "AF Index", row$viability))
    asinh(exprs(ff)[, chans] / COFACTOR)
  })
  markers <- colnames(arm_art[[REFERENCE_ARM]])
  stopifnot(
    "same marker channels in every arm" = all(map_lgl(arm_art, \(m) {
      setequal(colnames(m), markers)
    }))
  )
  imap(arm_art, \(m, a) {
    fr <- flowFrame(m[, markers])
    cut <- map_dbl(set_names(markers), \(ch) {
      cut <- getExportedValue(
        "openCyto",
        paste0("gate_", cfg("mindensity_method"))
      )(fr, channel = ch, gate_range = CUT_RANGE)@min[[1]]
      if (cut <= min(m[, ch])) NA_real_ else cut # a cut at the data minimum is no valley
    })
    tibble(
      arm = a,
      sample = s,
      channel = markers,
      cut = unname(cut),
      found = !is.na(cut)
    )
  }) |>
    list_rbind()
}
out <- parallel::mclapply(samples, fit_cuts, mc.cores = N_CORES)
stopifnot(
  "a cut worker failed" = !any(map_lgl(out, \(o) {
    inherits(o, "try-error") || is.null(o)
  }))
)
cuts <- list_rbind(out)
write_csv(cuts, file.path(RUN_DIR, "mindensity_cuts.csv"))

if (file.exists(OLD_CUTS)) {
  old <- read_csv(OLD_CUTS, show_col_types = FALSE)
  cuts |>
    select(arm, sample, channel, cut_new = cut) |>
    left_join(
      select(old, arm, sample, channel, cut_old = cut),
      by = join_by(arm, sample, channel),
      relationship = "one-to-one"
    ) |>
    write_csv(file.path(RUN_DIR, "mindensity_cuts_vs_old.csv"))
} else {
  message(
    "02c_metrics_common.R: no old cuts at ",
    OLD_CUTS,
    " (02_metrics.R has not run on dataset ",
    DS$dataset,
    "): mindensity_cuts_vs_old.csv is not written"
  )
}
writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_02c.txt")
)
