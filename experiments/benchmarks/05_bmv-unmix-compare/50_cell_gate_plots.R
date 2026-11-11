# Purpose: plots for two gating choices the user makes by eye, no unmixing. (1) FSC-A x SSC-A of the cell
#          universal negative and three cell controls with candidate cell rectangles (current FSC floor
#          6e5, and 4e5 and 3e5), with events kept per candidate. (2) AutoSpectral 1.8.8 tune.gate() grid
#          (landmark gating, docs article 08) for the stained cell controls.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/50_cell_gate_plots.R
#          from the repo root.
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group/  raw controls [zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/
#            Zombie NIR (Cells).fcs, 84,528-event file (user, 2026-11-10, D13); fcs_control_file.csv
# Outputs: experiments/runs/2026-11-10_bmv-cell-gating/cell_rectangles.png   FSC x SSC with candidates
#          experiments/runs/2026-11-10_bmv-cell-gating/cell_rectangles_kept.csv  events kept per candidate
#          experiments/runs/2026-11-10_bmv-cell-gating/fcs_control_file_landmark.csv  control file copy
#            with gate.name and gate.define for landmark gating (the original is not touched)
#          experiments/runs/2026-11-10_bmv-cell-gating/figure_gate_tuning/  tune.gate() grid, PDF and PNG
#          experiments/runs/2026-11-10_bmv-cell-gating/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(AutoSpectral)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
DATA_DIR <- "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8"
RAW_DIR <- file.path(DATA_DIR, "Reference Group")
CLEAN_DIR <- file.path(DATA_DIR, "Reference Group Clean")
# Files to show: the cell universal negative, two stained cell controls, and the viability control
# (no CD3/CD4/CD8 cell controls exist in this panel; those markers were run on beads)
PLOT_FILES <- c(
  `Unstained (Cells)` = file.path(RAW_DIR, "Unstained (Cells).fcs"),
  `Alexa Fluor 488 (ATP5A)` = file.path(RAW_DIR, "Alexa Fluor 488 (Cells).fcs"),
  `PE-Cy7 (ATGL)` = file.path(RAW_DIR, "PE-Cy7 (Cells).fcs"),
  `Zombie NIR (84,528-event file)` = file.path(
    CLEAN_DIR,
    "Zombie NIR (Cells).fcs"
  )
)
# Candidate cell rectangles: the zip gate (AutoSpectral_unmix_v1.8.8.R:39-42) and two lower FSC floors
CANDIDATES <- tribble(
  ~candidate         , ~fsc_min , ~fsc_max , ~ssc_min , ~ssc_max ,
  "FSC >= 6e5 (zip)" ,      6e5 ,      4e6 ,      1e5 , 1.6e6    ,
  "FSC >= 4e5"       ,      4e5 ,      4e6 ,      1e5 , 1.6e6    ,
  "FSC >= 3e5"       ,      3e5 ,      4e6 ,      1e5 , 1.6e6
)
CANDIDATE_COLOURS <- c(
  "FSC >= 6e5 (zip)" = "black",
  "FSC >= 4e5" = "#D55E00",
  "FSC >= 3e5" = "#0072B2"
)
# Stained cell controls that share one landmark gate and define it (user approved landmark gating,
# 2026-11-10, O2); viability keeps its own gate, the unstained uses the gate without defining it
LANDMARK_DEFINE <- c(
  "Alexa Fluor 405",
  "Alexa Fluor 488",
  "Alexa Fluor 594",
  "PE-Cy5.5",
  "PE-Cy7",
  "APC-Cy7"
)
TUNE_N_CELLS <- c(100, 500, 2000) # user, 2026-11-10 (O2), docs article 08 grid
TUNE_PERCENTILES <- c(30, 50, 70)

run_dir <- "experiments/runs/2026-11-10_bmv-cell-gating"
stopifnot("run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir)

# Read scatter of each file -------------------------------------------------------------------
scatter <- PLOT_FILES |>
  map(\(f) {
    exprs(read.FCS(
      f,
      truncate_max_range = FALSE,
      column.pattern = "^(FSC|SSC)-A$"
    )) |>
      as_tibble()
  }) |>
  list_rbind(names_to = "file") |>
  mutate(file = factor(file, levels = names(PLOT_FILES)))

# Events kept per file and candidate -----------------------------------------------------------
kept <- scatter |>
  cross_join(CANDIDATES) |>
  summarise(
    n_total = n(),
    n_kept = sum(
      `FSC-A` >= fsc_min &
        `FSC-A` <= fsc_max &
        `SSC-A` >= ssc_min &
        `SSC-A` <= ssc_max
    ),
    .by = c(file, candidate)
  ) |>
  mutate(pct_kept = 100 * n_kept / n_total)
write_csv(kept, file.path(run_dir, "cell_rectangles_kept.csv"))

# Figure: FSC x SSC with the candidate rectangles ------------------------------------------------
p_rect <- ggplot(scatter, aes(`FSC-A`, `SSC-A`)) +
  geom_hex(bins = 120) +
  scale_fill_viridis_c(trans = "log10", name = "Events") +
  geom_rect(
    data = CANDIDATES,
    aes(
      xmin = fsc_min,
      xmax = fsc_max,
      ymin = ssc_min,
      ymax = ssc_max,
      colour = candidate
    ),
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.6
  ) +
  scale_colour_manual(values = CANDIDATE_COLOURS, name = "Cell gate") +
  coord_cartesian(xlim = c(0, 4.2e6), ylim = c(0, 2.5e6)) +
  facet_wrap(~file, ncol = 2) +
  labs(x = "FSC-A", y = "SSC-A") +
  theme_bw(base_size = 11)
ggsave(
  file.path(run_dir, "cell_rectangles.png"),
  p_rect,
  width = 10,
  height = 8,
  dpi = 130,
  bg = "white"
)

# Landmark gating grid ------------------------------------------------------------------------
# control file copy: one shared landmark gate for the stained cell controls; originals untouched
control_copy <- read_csv(
  file.path(CLEAN_DIR, "fcs_control_file.csv"),
  show_col_types = FALSE
) |>
  mutate(
    gate.name = if_else(
      control.type == "cells" & !is.viability,
      "lymphocyte",
      gate.name
    ),
    gate.define = if_else(
      control.type == "cells" & !is.viability,
      fluorophore %in% LANDMARK_DEFINE,
      gate.define
    )
  )
control_copy_file <- file.path(run_dir, "fcs_control_file_landmark.csv")
write_csv(control_copy, control_copy_file)

asp <- get.autospectral.param(cytometer = "aurora", figures = TRUE)
tune.gate(
  control.file = control_copy_file,
  control.dir = RAW_DIR,
  asp = asp,
  gate.name = "lymphocyte",
  n.cells = TUNE_N_CELLS,
  percentiles = TUNE_PERCENTILES,
  output.dir = file.path(run_dir, "figure_gate_tuning")
)

# PNG copy of the tuning grid for viewing (tune.gate writes only a PDF)
tuning_pdf <- list.files(
  file.path(run_dir, "figure_gate_tuning"),
  "\\.pdf$",
  full.names = TRUE
)
system2(
  "pdftoppm",
  c(
    "-r",
    "110",
    "-png",
    shQuote(tuning_pdf),
    shQuote(sub("\\.pdf$", "", tuning_pdf))
  )
)

writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
