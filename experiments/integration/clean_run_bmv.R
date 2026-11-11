# Purpose: real-data smoke test of sw_clean_run() (margins, arcsinh, PeacoQC, nonDebris, singlets, live)
#          on the five BMV donors of one unmixed arm (as188_joint, AutoSpectral 1.8.8 joint), pooled
#          into one gate fit. Uses the unmixed, UNCLEANED FCS files, the same inputs that fed benchmark
#          05's 12_cleanup.R. Nothing is unmixed here. Every sw_clean_run() value below is a TEST VALUE
#          the user approved for this smoke test; none is a recommendation.
#          No seed: PeacoQC and openCyto mindensity do not sample (checked: PeacoQC draws random numbers
#          only in two plot builders, which are not called).
# Run:     sbatch experiments/slurm/run_integration.sbatch \
#            experiments/integration/clean_run_bmv.R <pkg_dir> <repo_root> \
#            <repo_root>/experiments/runs/2026-11-11_clean-run-smoke
#            pkg_dir    package checkout to load (holds R/clean_gate.R)
#            repo_root  checkout that holds experiments/bmv_test (the data) and experiments/runs
#            out_dir    new run folder; must not exist
# Inputs:  <repo_root>/experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed/*.fcs
#            one file per donor A3, A6, B11, B2, D1 (arm as188_joint in 01_manifest-gate.R) [user]
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

# Parameters -------------------------------------------------------------------
# TEST VALUES, user-approved for this smoke test only; not recommendations.
DEBRIS_WINDOW <- c(2e5, 1.5e6) # TEST VALUE; inherited from 19a_gate_arm.R:32 (raw FSC-A)
SINGLET_WIDER_GATE <- TRUE # TEST VALUE
POOL_SAMPLES <- TRUE # TEST VALUE; exercises the pooled fit over the five donors
VIABILITY_CHANNEL <- "Zombie NIR-A" # TEST VALUE; the dye channel of the as188_joint arm
VIABILITY_METHOD <- "mindensity" # TEST VALUE
VIABILITY_CUT <- NULL # TEST VALUE; NULL because the method is mindensity
IT_LIMIT <- 0.55 # TEST VALUE; PeacoQC isolation tree limit
MAD_LIMIT <- 6 # TEST VALUE; PeacoQC MAD limit
COFACTOR <- 6000 # TEST VALUE; arcsinh cofactor
# Not sw_clean_run() values
ARM_DIR <- "autoSpectral_v1.8.8/AutoSpectral_unmixed" # arm as188_joint in 01_manifest-gate.R
DATA_DIR <- "experiments/bmv_test/compare_3D/data"
DENSITY_POINTS <- 512 # grid points of each plotted density
FIG_WIDTH_MM <- 183 # Nature double column
FIG_HEIGHT_MM <- 65
DONOR_COLOURS <- c(
  A3 = "#4477AA",
  A6 = "#EE6677",
  B11 = "#228833",
  B2 = "#CCBB44",
  D1 = "#66CCEE"
) # Tol bright, one colour per donor, same in every figure
STEPS <- c("raw", "margins", "qc", "non_debris", "singlets", "live") # run order of the audit

# Run dir and warning log ------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
REPO <- args[[2]]
OUT_DIR <- args[[3]]
stopifnot("run dir must not exist" = !dir.exists(OUT_DIR))
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
manifest <- tibble(
  file = list.files(
    file.path(REPO, DATA_DIR, ARM_DIR),
    pattern = "\\.fcs$",
    full.names = TRUE
  )
) |>
  mutate(donor = str_extract(basename(file), "^[A-Z][0-9]+")) |>
  arrange(donor)
# Donor labels in the file names must be the five donors of the palette (also guards the plots).
stopifnot(
  "five donors found" = identical(manifest$donor, sort(names(DONOR_COLOURS)))
)
ffs <- manifest$file |>
  map(sw_io_read_fcs) |>
  set_names(manifest$donor)

# Clean ------------------------------------------------------------------------
# sw_gate_run() attaches flowCore, so filter() below is called as dplyr::filter().
res <- sw_clean_run(
  ffs,
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
raw_values <- ffs |>
  imap(\(ff, donor) {
    e <- flowCore::exprs(ff)
    tibble(
      donor = donor,
      fsc_a = e[, "FSC-A"],
      zombie = asinh(e[, VIABILITY_CHANNEL] / COFACTOR)
    )
  }) |>
  list_rbind()
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
    donor = sample,
    cut = if_else(
      is.finite(as.numeric(lower)),
      as.numeric(lower),
      as.numeric(upper)
    )
  ) |>
  select(donor, gate, cut, cut_at_window_floor)

fsc_density <- raw_values |>
  reframe(
    density_frame(fsc_a, min(raw_values$fsc_a), max(raw_values$fsc_a)),
    .by = donor
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
    .by = donor
  ) |>
  mutate(donor = factor(donor, levels = names(DONOR_COLOURS)))
zombie_cuts <- cuts |>
  dplyr::filter(gate == "live")
zombie_caption <- "All events as read, arcsinh(x / cofactor). Black: fitted live cut, one per donor (one line when pooled)."

kept <- res$audit |>
  mutate(
    donor = factor(sample, levels = names(DONOR_COLOURS)),
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
  select(donor, step, measure, value)

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
  labs(
    title = "Zombie NIR-A and live cut",
    x = "Zombie NIR-A (arcsinh, cofactor 6000)",
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
  aes(x = step, y = value, colour = donor, group = donor)
) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.2) +
  facet_wrap(vars(measure), scales = "free_y") +
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
