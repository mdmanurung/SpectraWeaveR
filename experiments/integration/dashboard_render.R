# Purpose: render the packaged dashboards once from the real 88h spike score and the 75/86b screen
#          (plan v2 step 20 and decision R3.6): a run folder with the report and F7/F8, the single-file
#          dashboard with its Spike & screen tab, and the Quarto dashboard. Local files only.
# Run:     pixi run Rscript experiments/integration/dashboard_render.R <pkg_dir> <repo_root> <out_dir>
#            pkg_dir    package checkout to load
#            repo_root  checkout that holds experiments/bmv_test and experiments/runs (the data)
#            out_dir    new run folder; must not exist
# Inputs:  <repo_root>/experiments/runs/2026-11-11_pkg-spike-regression/score_88h_TNF_T.rds  [spike_bmv.R]
#          <repo_root>/experiments/runs/2026-11-11_pkg-screen-regression/screen.rds          [screen_bmv.R]
#          <repo_root>/experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/  1.8.8 spectra
# Outputs: <out_dir>/run_88h/report/          sw_unmix_report() with spike and screen sections and CSVs
#          <out_dir>/run_88h/figures/panels/  F7, F8 (PNG 300 dpi, PDF, source_data CSV)
#          <out_dir>/dashboard.html           sw_unmix_dashboard() of run_88h
#          <out_dir>/quarto/                  sw_unmix_dashboard_quarto(): dashboard.html, img/, img_log.csv
#          <out_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  pkgload::load_all(
    commandArgs(trailingOnly = TRUE)[1],
    quiet = TRUE,
    export_all = FALSE
  )
  library(dplyr)
})

# Parameters ---------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
REPO <- args[2]
OUT_DIR <- args[3]
READOUT_LEVEL <- 1 # +1 arcsinh, as the 88h figures (inherited, Q6)
PSEUDO_LOG_SIGMA <- 1000 # dashboard event binning (user, 2026-11-09, OMIP dashboards); no event tables here
BINS <- 200 # as above
RUNS <- file.path(REPO, "experiments/runs")
SPECTRA_FILE <- file.path(
  REPO,
  "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
)

stopifnot("out dir is new" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE)
score <- readRDS(file.path(
  RUNS,
  "2026-11-11_pkg-spike-regression/score_88h_TNF_T.rds"
))
screen <- readRDS(file.path(
  RUNS,
  "2026-11-11_pkg-screen-regression/screen.rds"
))

# Run folder: report and figures --------------------------------------------------------------
run_dir <- file.path(OUT_DIR, "run_88h")
spectra <- as.matrix(read.csv(SPECTRA_FILE, row.names = 1, check.names = FALSE))
result <- sw_unmix_result(spectra = spectra, run_id = "bmv-88h-passes")
dx <- sw_unmix_diagnose(result, spike = score, screen = screen)
sw_unmix_report(
  result,
  dx,
  path = file.path(run_dir, "report"),
  notes = "Real BMV data: 88h pass settings (TNF spike, CD4T and CD8T) and the 75/86b marker screen. No limits were passed, so the new rows are not_judged."
)
figs <- list(
  F7 = sw_unmix_figure_spike_summary(
    score,
    READOUT_LEVEL,
    stem = "spike-summary_TNF_T"
  ),
  F8 = sw_unmix_figure_spike_ecdf(
    score,
    READOUT_LEVEL,
    stem = "spike-ecdf_TNF_T"
  )
)
saved <- sw_unmix_save_figures(
  figs,
  dir = file.path(run_dir, "figures/panels"),
  check = "nature"
)

# Single-file dashboard -------------------------------------------------------------------------
sw_unmix_dashboard(
  c(`88h passes` = run_dir),
  path = file.path(OUT_DIR, "dashboard.html"),
  pseudo_log_sigma = PSEUDO_LOG_SIGMA,
  bins = BINS,
  title = "BMV 88h spike passes"
)

# Quarto dashboard: numbers from the score, figures from the panels ----------------------------
# by_setup carries n_coupled: donors whose file failed the independence check
lvl <- score$by_setup |>
  filter(level == READOUT_LEVEL)
verdicts <- lvl |>
  filter(pop == "CD8T") |>
  transmute(
    label = sprintf("%s, %s %s: cells losing the spike", setup, marker, pop),
    value = sprintf("%.1f%% (donor median)", 100 * share_lost_median),
    note = sprintf(
      "%d donors, range %.1f-%.1f%%; %d files coupled",
      n_donors,
      100 * share_lost_min,
      100 * share_lost_max,
      n_coupled
    ),
    # colour = the independence check, not a judgement of the loss
    status = if_else(n_coupled > 0, "bad", "good")
  )
panel_png <- \(stem) file.path(run_dir, "figures/panels", paste0(stem, ".png"))
finding_tbl <- lvl |>
  transmute(
    text = sprintf(
      "- %s, %s: %.1f%% of cells lose the spike (donor median; range %.1f-%.1f%%); median recovery %.3f; %d of %d files coupled.",
      setup,
      pop,
      100 * share_lost_median,
      100 * share_lost_min,
      100 * share_lost_max,
      recovery_median,
      n_coupled,
      n_donors
    )
  )
flag_tbl <- screen$flags |>
  summarise(
    flagged = sum(flag, na.rm = TRUE),
    na = sum(is.na(flag)),
    n = n(),
    .by = setup
  )
sections <- list(
  list(
    heading = "TNF spike-in recovery by pass setting",
    finding = paste(
      c(
        sprintf(
          "TNF (BV510) spiked at +%s arcsinh into CD4T and CD8T cells of 5 donors; the score compares each setting's spiked and unspiked unmixes cell by cell. Colours of the cards follow the independence check: a coupled file has unspiked cells that moved by more than %s arcsinh.",
          READOUT_LEVEL,
          score$baseline_tol
        ),
        "",
        finding_tbl$text
      ),
      collapse = "\n"
    ),
    figures = data.frame(
      path = c(panel_png("spike-summary_TNF_T"), panel_png("spike-ecdf_TNF_T")),
      caption = c(
        "F7: share of cells losing the spike and donor median recovery per setting.",
        "F8: cumulative distribution of per-cell recovery per setting."
      )
    ),
    notes = score$notes
  ),
  list(
    heading = "Marker screen against as188_joint",
    finding = paste(
      c(
        sprintf("Rule: %s.", screen$rule$text),
        "",
        sprintf(
          "- %s: %d of %d channel-population pairs flagged; %d NA.",
          flag_tbl$setup,
          flag_tbl$flagged,
          flag_tbl$n,
          flag_tbl$na
        )
      ),
      collapse = "\n"
    ),
    figures = NULL,
    notes = character()
  )
)
sw_unmix_dashboard_quarto(
  sections,
  path = file.path(OUT_DIR, "quarto"),
  title = "BMV unmix diagnostics (package render)",
  verdicts = verdicts,
  subtitle = "88h pass settings and the 75/86b marker screen, rendered with sw_unmix_dashboard_quarto()",
  footer = "SpectraWeaveR package render, 2026-11-11; local file, not published."
)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
