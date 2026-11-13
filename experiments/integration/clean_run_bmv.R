# Purpose: real-data smoke test of sw_clean_run() (margins, arcsinh, PeacoQC, nonDebris, singlets, live)
#          on the samples of one dataset (SW_DATASET: 3d, 24h or exv) in one unmixed arm (as188_joint,
#          AutoSpectral 1.8.8 joint), pooled into one gate fit per gate group (config gate_group: 24H =
#          condition, so sw_clean_run() is called once per condition; 3D and EXV have no group and one call).
#          Uses the unmixed, UNCLEANED FCS files, the same inputs that fed benchmark 05's 12_cleanup.R.
#          Nothing is unmixed here. Every sw_clean_run() value below is a TEST VALUE the user approved for
#          the 3D smoke test; none is a recommendation, and none was re-approved for 24H or EXV.
#          No seed: PeacoQC and openCyto mindensity do not sample (checked: PeacoQC draws random numbers
#          only in two plot builders, which are not called).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_integration.sbatch \
#            experiments/integration/clean_run_bmv.R <pkg_dir> <repo_root> \
#            <repo_root>/experiments/runs/2026-11-11_clean-run-smoke2 <dataset_json>
#            or, from the repo root, SW_DATASET=<ds> pixi run Rscript experiments/integration/clean_run_bmv.R \
#            <pkg_dir> <repo_root> <out_dir> <dataset_json>
#            pkg_dir       package checkout to load (holds R/clean_gate.R)
#            repo_root     checkout that holds experiments/bmv_test (the data) and experiments/runs; must be the
#                          working directory, because the registry paths are relative to it
#            out_dir       new run folder; must not exist
#            dataset_json  required, no default: experiments/benchmarks/05_bmv-unmix-compare/datasets/<ds>.json;
#                          its file name must be $SW_DATASET.json
# Inputs:  <repo_root>/<data_dir>/autoSpectral_v1.8.8/AutoSpectral_unmixed/**.fcs  [user zip; data_dir from <dataset_json>]
#            one file per sample of datasets/<ds>_samples.csv (24H: in one subfolder per condition),
#            joined on the file stem one-to-one; a missing or extra file stops the run [00_samples.R]
#          <dataset_json> channels (viability_autospectral) and palette (colour per plotted unit)
#          No registry run key is read (no run_in call).
# Outputs: <out_dir>/audit.csv                    sample x step: events in and out (sw_clean_run $audit)
#          <out_dir>/audit_total.csv              the same over all donors
#          <out_dir>/gates.csv                    sample x gate: shape, keeps, cut_at_window_floor
#          <out_dir>/qc.csv                       PeacoQC summary per donor with IT_limit and MAD
#          <out_dir>/template.csv                 openCyto template that was run
#          <out_dir>/parameters.csv               the values sw_clean_run() received
#          <out_dir>/warnings.txt                 every warning raised in the script, one per line
#          <out_dir>/figures/fig_<name>.{png,pdf} 300 dpi PNG and vector PDF, three figures
#          <out_dir>/source_data/fig_<name>_*.csv the plotted data of each figure
#          <out_dir>/sessionInfo.txt              written last

suppressPackageStartupMessages({
  pkgload::load_all(
    commandArgs(trailingOnly = TRUE)[[1]],
    quiet = TRUE,
    export_all = FALSE
  )
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

# Dataset config ---------------------------------------------------------------
# <dataset_json> (4th argument) names the dataset; dataset.R reads its folder (config, registry, sample table).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "clean_run_bmv.R: needs 4 arguments: <pkg_dir> <repo_root> <out_dir> <dataset_json>",
    call. = FALSE
  )
}
DATASET_JSON <- normalizePath(args[[4]], mustWork = TRUE)
if (basename(DATASET_JSON) != paste0(Sys.getenv("SW_DATASET"), ".json")) {
  stop(
    "clean_run_bmv.R: <dataset_json> must be $SW_DATASET.json; got ",
    DATASET_JSON,
    " with SW_DATASET='",
    Sys.getenv("SW_DATASET"),
    "'",
    call. = FALSE
  )
}
Sys.setenv(SW_DATASET_DIR = dirname(DATASET_JSON))
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), samples_tbl()

# Parameters -------------------------------------------------------------------
# TEST VALUES, user-approved for this smoke test only; not recommendations.
MARGIN_RANGE <- c(-21761, 1e7) # TEST VALUE; the margin range benchmark 05 passed (12_cleanup.R:28)
DEBRIS_WINDOW <- c(2e5, 1.5e6) # TEST VALUE; inherited from 19a_gate_arm.R:32 (raw FSC-A)
SINGLET_WIDER_GATE <- TRUE # TEST VALUE
POOL_SAMPLES <- TRUE # TEST VALUE; exercises the pooled fit over the five donors
VIABILITY_CHANNEL <- cfg("channels")$viability_autospectral # TEST VALUE (3D: Zombie NIR-A); the dye channel of the as188_joint arm
VIABILITY_METHOD <- "mindensity" # TEST VALUE
VIABILITY_CUT <- NULL # TEST VALUE; NULL because the method is mindensity
IT_LIMIT <- 0.55 # TEST VALUE; PeacoQC isolation tree limit
MAD_LIMIT <- 6 # TEST VALUE; PeacoQC MAD limit
COFACTOR <- 6000 # TEST VALUE; arcsinh cofactor
# Not sw_clean_run() values
ARM_DIR <- "autoSpectral_v1.8.8/AutoSpectral_unmixed" # arm as188_joint in 01_manifest-gate.R; same folder in every dataset
DATA_DIR <- data_path() # config data_dir
DENSITY_POINTS <- 512 # grid points of each plotted density
FIG_WIDTH_MM <- 183 # Nature double column
FIG_HEIGHT_MM <- 65
DONOR_COLOURS <- unlist(cfg("palette")) # config palette: one colour per plotted unit, same in every figure (3D: Tol bright, per well)
GATE_GROUP <- DS$gate_group # NULL: one gate fit over all samples; else the sample-table column that splits the fits
COND_COL <- if (is.null(DS$conditions)) character() else "condition" # condition column in tables and facets, only when the dataset has conditions
UNIT_COL <- if (is.null(DS$conditions)) "sample" else "donor" # plotted unit: the sample (3D: the well); with conditions the donor, one panel per condition
STEPS <- c("raw", "margins", "qc", "non_debris", "singlets", "live") # run order of the audit

# Run dir and warning log ------------------------------------------------------
REPO <- args[[2]]
OUT_DIR <- args[[3]]
stopifnot(
  "repo_root is the working directory (registry paths are relative to it)" = identical(
    normalizePath(REPO),
    normalizePath(getwd())
  ),
  "run dir must not exist" = !dir.exists(OUT_DIR)
)
walk(
  file.path(OUT_DIR, c("figures", "source_data")),
  dir.create,
  recursive = TRUE
)
warning_file <- file.path(OUT_DIR, "warnings.txt")
file.create(warning_file)
# Each warning is appended at once (kept even if the script later fails) and echoed to the job log.
globalCallingHandlers(
  warning = \(w) {
    line <- str_c(
      conditionMessage(w),
      " [call: ",
      deparse(conditionCall(w))[[1]],
      "]"
    )
    cat(line, "\n", file = warning_file, append = TRUE)
    message("WARNING: ", line)
    invokeRestart("muffleWarning")
  }
)

# Figure theme -----------------------------------------------------------------
theme_set(
  ggprism::theme_prism(
    base_size = 7,
    base_family = "Nimbus Sans",
    base_fontface = "plain",
    base_line_size = 0.4
  ) +
    theme(
      plot.title = element_text(face = "bold", size = 7),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      strip.clip = "off"
    )
)

# Manifest and files -----------------------------------------------------------
# One file per sample: the file stem (AutoSpectral appends " AutoSpectral" or " AutoSpectral_v1-8-8", cut by the pattern) joins one-to-one to the sample table,
# so a missing, doubled or stray .fcs file stops the run. recursive: 24H arm folders hold one subfolder per condition.
# sample_units: the plotted unit of each sample (column `donor` of the plot frames; 3D: the well).
sample_units <- samples_tbl() |>
  transmute(
    sample,
    stem,
    donor = .data[[UNIT_COL]],
    condition,
    gate_by = .data[[GATE_GROUP %||% "sample"]]
  )
manifest <- tibble(
  file = list.files(
    file.path(REPO, DATA_DIR, ARM_DIR),
    pattern = "[.]fcs$",
    recursive = TRUE,
    full.names = TRUE
  )
) |>
  mutate(
    stem = str_remove(
      tools::file_path_sans_ext(basename(file)),
      " AutoSpectral.*$"
    )
  ) |>
  inner_join(
    sample_units,
    by = "stem",
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  arrange(sample)
# every plotted unit needs a palette colour (the scale would otherwise drop it silently)
stopifnot(
  "palette covers every plotted unit" = all(
    manifest$donor %in% names(DONOR_COLOURS)
  )
)
ffs <- manifest$file |>
  map(sw_io_read_fcs) |>
  set_names(manifest$sample)

# Clean ------------------------------------------------------------------------
# sw_gate_run() attaches flowCore, so filter() below is called as dplyr::filter().
# One call per gate group (config gate_group; 24H: one per condition). Without a group it is one call on all samples.
gate_groups <- split(manifest$sample, manifest$gate_by)
if (is.null(GATE_GROUP)) {
  gate_groups <- list(manifest$sample)
}
res_by_group <- map(gate_groups, \(group_samples) {
  sw_clean_run(
    ffs[group_samples],
    margin_range = MARGIN_RANGE,
    debris_window = DEBRIS_WINDOW,
    singlet_wider_gate = SINGLET_WIDER_GATE,
    pool_samples = POOL_SAMPLES,
    viability_channel = VIABILITY_CHANNEL,
    viability_method = VIABILITY_METHOD,
    viability_cut = VIABILITY_CUT,
    IT_limit = IT_LIMIT,
    MAD = MAD_LIMIT,
    cofactor = COFACTOR
  )
})
# one table per result name; with a gate group the group is added as column gate_group
res_table <- \(name) {
  tables <- map(res_by_group, name)
  if (is.null(GATE_GROUP)) {
    tables[[1]]
  } else {
    list_rbind(tables, names_to = "gate_group")
  }
}
res <- list(
  audit = res_table("audit"),
  audit_total = res_table("audit_total"),
  gates = res_table("gates"),
  qc = res_table("qc"),
  template = res_table("template"),
  parameters = res_by_group[[1]]$parameters # the same values went into every call
)
# Tables, written first so they survive a figure-stage error -------------------
write_csv(res$audit, file.path(OUT_DIR, "audit.csv"))
write_csv(res$audit_total, file.path(OUT_DIR, "audit_total.csv"))
write_csv(res$gates, file.path(OUT_DIR, "gates.csv"))
write_csv(res$qc, file.path(OUT_DIR, "qc.csv"))
write_csv(res$template, file.path(OUT_DIR, "template.csv"))
res$parameters |>
  imap(\(value, name) tibble(parameter = name, value = deparse1(value))) |>
  list_rbind() |>
  write_csv(file.path(OUT_DIR, "parameters.csv"))

message("audit_total:")
print(res$audit_total)
message("cut_at_window_floor:")
print(
  dplyr::filter(res$gates, gate == "nonDebris") |>
    select(sample, definition, cut_at_window_floor)
)

# Plot frames ------------------------------------------------------------------
# Densities use every event as read (before margins), one curve per donor on a common grid.
# FRAME_BY: columns that identify one curve (3D: the donor, i.e. the well; with conditions: donor and condition)
FRAME_BY <- c("donor", COND_COL)
raw_values <- ffs |>
  imap(\(ff, which_sample) {
    e <- flowCore::exprs(ff)
    tibble(
      sample = which_sample,
      fsc_a = e[, "FSC-A"],
      zombie = asinh(e[, VIABILITY_CHANNEL] / COFACTOR)
    )
  }) |>
  list_rbind() |>
  left_join(
    select(sample_units, sample, donor, condition),
    by = "sample",
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  select(all_of(FRAME_BY), fsc_a, zombie)
density_frame <- \(values, grid_from, grid_to) {
  with(
    density(values, n = DENSITY_POINTS, from = grid_from, to = grid_to),
    tibble(x = x, density = y)
  )
}
# Cuts: the one finite limit of each rectangle gate (events at or beyond it are cut).
cuts <- res$gates |>
  dplyr::filter(gate %in% c("nonDebris", "live")) |>
  separate_wider_regex(
    definition,
    c(".*: \\[", lower = "[^,]+", ", ", upper = "[^\\]]+", "\\]"),
    cols_remove = FALSE
  ) |>
  mutate(
    cut = if_else(
      is.finite(as.numeric(lower)),
      as.numeric(lower),
      as.numeric(upper)
    )
  ) |>
  left_join(
    select(sample_units, sample, donor, condition),
    by = "sample",
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  select(all_of(FRAME_BY), gate, cut, cut_at_window_floor)

fsc_density <- raw_values |>
  reframe(
    density_frame(fsc_a, min(raw_values$fsc_a), max(raw_values$fsc_a)),
    .by = all_of(FRAME_BY)
  ) |>
  mutate(
    density = density * 1e6,
    donor = factor(donor, levels = names(DONOR_COLOURS))
  )
fsc_cuts <- cuts |>
  dplyr::filter(gate == "nonDebris")
fsc_window <- tibble(edge = c("min", "max"), fsc_a = DEBRIS_WINDOW)
fsc_caption <- "All events as read. Dashed: debris_window. Black: fitted nonDebris cut, one per donor (one line when pooled)."

zombie_density <- raw_values |>
  reframe(
    density_frame(zombie, min(raw_values$zombie), max(raw_values$zombie)),
    .by = all_of(FRAME_BY)
  ) |>
  mutate(donor = factor(donor, levels = names(DONOR_COLOURS)))
zombie_cuts <- cuts |>
  dplyr::filter(gate == "live")
zombie_caption <- "All events as read, arcsinh(x / cofactor). Black: fitted live cut, one per donor (one line when pooled)."

kept <- res$audit |>
  left_join(
    select(sample_units, sample, donor, condition),
    by = "sample",
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(
    donor = factor(donor, levels = names(DONOR_COLOURS)),
    step = factor(step, levels = STEPS)
  ) |>
  pivot_longer(
    c(n_events, pct_of_initial),
    names_to = "measure",
    values_to = "value"
  ) |>
  mutate(
    measure = factor(
      measure,
      levels = c("n_events", "pct_of_initial"),
      labels = c("Events kept", "% of raw events")
    )
  ) |>
  select(all_of(FRAME_BY), step, measure, value)

# Figures ----------------------------------------------------------------------
# One call per figure: PNG at 300 dpi, vector PDF, and the plotted data as CSV.
save_fig <- \(plot, name, frames) {
  path <- file.path(OUT_DIR, "figures", str_c("fig_", name))
  ggsave(
    str_c(path, ".png"),
    plot,
    width = FIG_WIDTH_MM,
    height = FIG_HEIGHT_MM,
    units = "mm",
    dpi = 300
  )
  ggsave(
    str_c(path, ".pdf"),
    plot,
    width = FIG_WIDTH_MM,
    height = FIG_HEIGHT_MM,
    units = "mm",
    device = cairo_pdf
  )
  iwalk(frames, \(frame, frame_name) {
    write_csv(
      frame,
      file.path(
        OUT_DIR,
        "source_data",
        str_c("fig_", name, "_", frame_name, ".csv")
      )
    )
  })
}

# 24H: one panel (row) per condition, axes on every panel; without conditions the figures keep their single panel
facet_condition <- if (is.null(DS$conditions)) {
  NULL
} else {
  facet_wrap(vars(condition), axes = "all", axis.labels = "margins")
}
facet_kept <- if (is.null(DS$conditions)) {
  facet_wrap(vars(measure), scales = "free_y")
} else {
  facet_grid(
    rows = vars(condition),
    cols = vars(measure),
    scales = "free_y",
    axes = "all",
    axis.labels = "margins"
  )
}

p_fsc <- ggplot(fsc_density, aes(x = x, y = density)) +
  geom_line(aes(colour = donor), linewidth = 0.4) +
  geom_vline(data = fsc_cuts, aes(xintercept = cut), linewidth = 0.3) +
  geom_vline(
    data = fsc_window,
    aes(xintercept = fsc_a),
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  scale_colour_manual(values = DONOR_COLOURS) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6)) +
  facet_condition +
  labs(
    title = "Raw FSC-A and debris cut",
    x = "FSC-A (x 1e6, raw)",
    y = "Density (per 1e6 FSC-A)",
    colour = "Donor",
    caption = fsc_caption
  )
save_fig(
  p_fsc,
  "fsc_debris",
  list(density = fsc_density, cuts = fsc_cuts, window = fsc_window)
)

p_zombie <- ggplot(zombie_density, aes(x = x, y = density)) +
  geom_line(aes(colour = donor), linewidth = 0.4) +
  geom_vline(data = zombie_cuts, aes(xintercept = cut), linewidth = 0.3) +
  scale_colour_manual(values = DONOR_COLOURS) +
  facet_condition +
  labs(
    title = paste(VIABILITY_CHANNEL, "and live cut"),
    x = paste0(VIABILITY_CHANNEL, " (arcsinh, cofactor ", COFACTOR, ")"),
    y = "Density",
    colour = "Donor",
    caption = zombie_caption
  )
save_fig(
  p_zombie,
  "live_zombie",
  list(density = zombie_density, cuts = zombie_cuts)
)

p_kept <- ggplot(
  kept,
  aes(x = step, y = value, colour = donor, group = donor) # one line per donor in each panel
) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.2) +
  facet_kept +
  scale_colour_manual(values = DONOR_COLOURS) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title = "Events kept per step",
    x = "Step (run order)",
    y = NULL,
    colour = "Donor",
    caption = "Counts from the sw_clean_run() audit table; each point is one donor."
  )
save_fig(p_kept, "events_kept", list(audit = kept))

# Session info (last) ----------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
