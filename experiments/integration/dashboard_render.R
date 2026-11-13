# Purpose: render the packaged dashboards once from the real 88h spike score and the 75/86b screen
#          (plan v2 step 20 and decision R3.6) of one dataset (SW_DATASET: 3d, 24h or exv): a run folder with the
#          report and F7/F8, the single-file dashboard with its Spike & screen tab, and the Quarto dashboard. Local
#          files only. With conditions (24H) the score and the screen exist once per condition (spike_bmv.R and
#          screen_bmv.R score each condition on its own), so there is one run folder per condition, all in the one
#          dashboard, and one section pair per condition in the Quarto dashboard.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_integration.sbatch \
#            experiments/integration/dashboard_render.R <pkg_dir> <repo_root> <out_dir> <dataset_json>
#            or, from the repo root, SW_DATASET=<ds> pixi run Rscript experiments/integration/dashboard_render.R \
#            <pkg_dir> <repo_root> <out_dir> <dataset_json>
#            pkg_dir       package checkout to load
#            repo_root     checkout that holds experiments/bmv_test and experiments/runs (the data); must be the
#                          working directory, because the registry paths are relative to it
#            out_dir       new run folder; must not exist
#            dataset_json  required, no default: experiments/benchmarks/05_bmv-unmix-compare/datasets/<ds>.json;
#                          its file name must be $SW_DATASET.json
# Inputs:  registry run keys (datasets/<ds>_runs.csv); the two keys name the out_dir you gave spike_bmv.R and
#          screen_bmv.R, and you register them after approving those runs:
#            pkg-spike-regression   score_88h_<tnf_label>_T.rds (24H: score_88h_<tnf_label>_T_<condition>.rds) [spike_bmv.R]
#            pkg-screen-regression  screen.rds (24H: screen_<condition>.rds)                                    [screen_bmv.R]
#          <data_dir>/autoSpectral_v1.8.8/table_spectra/<stamps.as188_spectra>   1.8.8 spectra
#          config targets (tnf_label, tnf_dye) and run_tag
# Outputs: <out_dir>/run_88h[_<condition>]/report/          sw_unmix_report() with spike and screen sections and CSVs
#          <out_dir>/run_88h[_<condition>]/figures/panels/  F7, F8 (PNG 300 dpi, PDF, source_data CSV)
#          <out_dir>/dashboard.html           sw_unmix_dashboard() of the run folders
#          <out_dir>/quarto/                  sw_unmix_dashboard_quarto(): dashboard.html, img/, img_log.csv
#          <out_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  pkgload::load_all(
    commandArgs(trailingOnly = TRUE)[1],
    quiet = TRUE,
    export_all = FALSE
  )
  library(dplyr)
  library(purrr)
})

# Dataset config -----------------------------------------------------------------------------
# <dataset_json> (4th argument) names the dataset; dataset.R reads its folder (config, registry, sample table).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "dashboard_render.R: needs 4 arguments: <pkg_dir> <repo_root> <out_dir> <dataset_json>",
    call. = FALSE
  )
}
DATASET_JSON <- normalizePath(args[[4]], mustWork = TRUE)
if (basename(DATASET_JSON) != paste0(Sys.getenv("SW_DATASET"), ".json")) {
  stop(
    "dashboard_render.R: <dataset_json> must be $SW_DATASET.json; got ",
    DATASET_JSON,
    " with SW_DATASET='",
    Sys.getenv("SW_DATASET"),
    "'",
    call. = FALSE
  )
}
Sys.setenv(SW_DATASET_DIR = dirname(DATASET_JSON))
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in()

# Parameters ---------------------------------------------------------------------------------
REPO <- args[2]
OUT_DIR <- args[3]
READOUT_LEVEL <- 1 # +1 arcsinh, as the 88h figures (inherited, Q6)
PSEUDO_LOG_SIGMA <- 1000 # dashboard event binning (user, 2026-11-09, OMIP dashboards); no event tables here
BINS <- 200 # as above
TNF_LABEL <- cfg("targets")$tnf_label # marker label of the TNF-slot dye in the spike score file names
TNF_DYE <- cfg("targets")$tnf_dye
TAG <- paste0("BMV", cfg("run_tag")) # title prefix: BMV (3D), BMV24h, BMVexv
SPECTRA_FILE <- data_path(
  "autoSpectral_v1.8.8/table_spectra",
  cfg("stamps")$as188_spectra
)
# scoring groups of spike_bmv.R and screen_bmv.R: one per condition (24H); "" = all samples, no suffix in file names
GROUPS <- if (is.null(DS$conditions)) "" else DS$conditions
group_suffix <- \(g) if (nzchar(g)) paste0("_", g) else ""

stopifnot(
  "repo_root is the working directory (registry paths are relative to it)" = identical(
    normalizePath(REPO),
    normalizePath(getwd())
  ),
  "out dir is new" = !dir.exists(OUT_DIR)
)
dir.create(OUT_DIR, recursive = TRUE)
spectra <- as.matrix(read.csv(SPECTRA_FILE, row.names = 1, check.names = FALSE))
result <- sw_unmix_result(
  spectra = spectra,
  run_id = paste0("bmv", cfg("run_tag"), "-88h-passes")
)

# Run folder: report and figures, one per scoring group -----------------------------------------
# panel_png(): the PNG of one panel stem inside a run folder
panel_png <- \(run_dir, stem) {
  file.path(run_dir, "figures/panels", paste0(stem, ".png"))
}
stem_summary <- paste0("spike-summary_", TNF_LABEL, "_T")
stem_ecdf <- paste0("spike-ecdf_", TNF_LABEL, "_T")
by_group <- map(set_names(GROUPS), \(g) {
  score <- readRDS(file.path(
    run_in("pkg-spike-regression"),
    paste0("score_88h_", TNF_LABEL, "_T", group_suffix(g), ".rds")
  ))
  screen <- readRDS(file.path(
    run_in("pkg-screen-regression"),
    paste0("screen", group_suffix(g), ".rds")
  ))
  run_dir <- file.path(OUT_DIR, paste0("run_88h", group_suffix(g)))
  dx <- sw_unmix_diagnose(result, spike = score, screen = screen)
  sw_unmix_report(
    result,
    dx,
    path = file.path(run_dir, "report"),
    notes = paste0(
      "Real BMV data",
      if (nzchar(g)) paste0(" (", DS$dataset, ", ", g, ")") else "",
      ": 88h pass settings (",
      TNF_LABEL,
      " spike, CD4T and CD8T) and the 75/86b marker screen. No limits were passed, so the new rows are not_judged."
    )
  )
  figs <- list(
    F7 = sw_unmix_figure_spike_summary(
      score,
      READOUT_LEVEL,
      stem = stem_summary
    ),
    F8 = sw_unmix_figure_spike_ecdf(score, READOUT_LEVEL, stem = stem_ecdf)
  )
  saved <- sw_unmix_save_figures(
    figs,
    dir = file.path(run_dir, "figures/panels"),
    check = "nature"
  )
  list(score = score, screen = screen, run_dir = run_dir)
})

# Single-file dashboard -------------------------------------------------------------------------
dashboard_runs <- set_names(
  map_chr(by_group, "run_dir"),
  paste0(
    "88h passes",
    map_chr(GROUPS, \(g) if (nzchar(g)) paste0(" ", g) else "")
  )
)
sw_unmix_dashboard(
  dashboard_runs,
  path = file.path(OUT_DIR, "dashboard.html"),
  pseudo_log_sigma = PSEUDO_LOG_SIGMA,
  bins = BINS,
  title = paste0(TAG, " 88h spike passes")
)

# Quarto dashboard: numbers from the score, figures from the panels ----------------------------
# by_setup carries n_coupled: donors whose file failed the independence check
verdicts <- imap(by_group, \(grp, g) {
  grp$score$by_setup |>
    filter(level == READOUT_LEVEL, pop == "CD8T") |>
    transmute(
      label = sprintf(
        "%s%s, %s %s: cells losing the spike",
        if (nzchar(g)) paste0(g, ", ") else "",
        setup,
        marker,
        pop
      ),
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
}) |>
  list_rbind()
sections <- imap(by_group, \(grp, g) {
  lvl <- filter(grp$score$by_setup, level == READOUT_LEVEL)
  finding_lines <- lvl |>
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
  flag_tbl <- grp$screen$flags |>
    summarise(
      flagged = sum(flag, na.rm = TRUE),
      na = sum(is.na(flag)),
      n = n(),
      .by = setup
    )
  where <- if (nzchar(g)) paste0(" (", g, ")") else ""
  list(
    list(
      heading = paste0(TNF_LABEL, " spike-in recovery by pass setting", where),
      finding = paste(
        c(
          sprintf(
            "%s (%s) spiked at +%s arcsinh into CD4T and CD8T cells of %d donors; the score compares each setting's spiked and unspiked unmixes cell by cell. Colours of the cards follow the independence check: a coupled file has unspiked cells that moved by more than %s arcsinh.",
            TNF_LABEL,
            TNF_DYE,
            READOUT_LEVEL,
            max(lvl$n_donors),
            grp$score$baseline_tol
          ),
          "",
          finding_lines$text
        ),
        collapse = "\n"
      ),
      figures = data.frame(
        path = c(
          panel_png(grp$run_dir, stem_summary),
          panel_png(grp$run_dir, stem_ecdf)
        ),
        caption = c(
          "F7: share of cells losing the spike and donor median recovery per setting.",
          "F8: cumulative distribution of per-cell recovery per setting."
        )
      ),
      notes = grp$score$notes
    ),
    list(
      heading = paste0("Marker screen against as188_joint", where),
      finding = paste(
        c(
          sprintf("Rule: %s.", grp$screen$rule$text),
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
}) |>
  list_flatten(name_spec = "{inner}") |>
  unname()
sw_unmix_dashboard_quarto(
  sections,
  path = file.path(OUT_DIR, "quarto"),
  title = paste0(TAG, " unmix diagnostics (package render)"),
  verdicts = verdicts,
  subtitle = "88h pass settings and the 75/86b marker screen, rendered with sw_unmix_dashboard_quarto()",
  footer = paste0(
    "SpectraWeaveR package render, ",
    Sys.getenv("SW_RUN_DATE"),
    "; local file, not published."
  )
)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
