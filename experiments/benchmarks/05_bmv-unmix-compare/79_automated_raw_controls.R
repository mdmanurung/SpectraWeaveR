# Purpose: test 5 (README decision table, user 2026-11-10): does the 1.8.8 automated spectra route keep the
#          AF405 (GLUT1) and SBUV510 (CD3) references free of a BV510 part even on the raw, ungated controls
#          that the 1.0.2 route used? Runs get.spectra.automated() with the zip 1.8.8 arguments and the zip
#          control file, on the raw 'Reference Group' FCS files instead of 'Reference Group Clean'. Then
#          measures each dye's reference offset vs the zip 1.8.8 spectra with the script 74 metric (projection
#          coefficient on the 1.8.8 BV510 spectrum). No unmixing.
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/79_automated_raw_controls.R \
#            experiments/runs/2026-11-10_bmv-automated-raw-controls        (repo root; run dir must not exist)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/*.fcs  raw controls [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv
#            zip control file (same file names as the raw folder) [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/ run 143803 spectra [user zip]
#          experiments/runs/2026-11-10_bmv-variant-bv510-by-set/set_summary.csv  zip_102, chain offsets [74]
# Outputs: <run_dir>/spectra_automated_raw.csv   dye spectra from the raw controls
#          <run_dir>/offsets.csv                 per dye: BV510 coefficient of (reference - 1.8.8 zip reference)
#          <run_dir>/fig_offsets.png             AF405 and SBUV510 offsets: zip_102, chain4_spec, chain5_bead,
#                                                this run (+ fig_offsets_source.csv)
#          <run_dir>/table_spectra/, figure_*/   AutoSpectral's own outputs
#          <run_dir>/args.txt, sessionInfo.txt

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# Parameters ---------------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
RUN_DIR <- "experiments/runs/2026-11-10_bmv-automated-raw-controls"
N_CANDIDATES <- 10000 # zip get.spectra.automated() call (AutoSpectral_unmix_v1.8.8.R:217)
N_SPECTRAL <- 500 # same
COSINE_THRESHOLD <- 0.95 # same
SEED <- 42 # settled seed (README); the zip run set none
DYES <- c("Alexa Fluor 405", "SBUV510")
TARGET <- "BV510"

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)

# Spectra from the raw controls, zip arguments, outputs written into the run dir ------------------------
asp <- get.autospectral.param(cytometer = "aurora", figures = TRUE) # zip: figures = TRUE
out_dirs <- grep("dir", names(asp), value = TRUE)
out_dirs <- out_dirs[vapply(asp[out_dirs], is.character, logical(1))]
asp[out_dirs] <- lapply(asp[out_dirs], \(d) file.path(RUN_DIR, d))
control_dir <- file.path(DATA, "autoSpectral_v1.0.2", "Reference Group")
control_file <- file.path(
  DATA,
  "autoSpectral_v1.8.8",
  "Reference Group Clean",
  "fcs_control_file.csv"
)
writeLines(
  c(
    paste("controls (raw):", control_dir),
    paste("control file (zip 1.8.8):", control_file),
    paste("n.candidates:", N_CANDIDATES),
    paste("n.spectral:", N_SPECTRAL),
    paste("cosine.threshold:", COSINE_THRESHOLD),
    "other arguments: package defaults (as the zip)",
    paste("seed:", SEED)
  ),
  file.path(RUN_DIR, "args.txt")
)
set.seed(SEED)
spectra_raw <- get.spectra.automated(
  control.dir = control_dir,
  control.def.file = control_file,
  asp = asp,
  n.candidates = N_CANDIDATES,
  n.spectral = N_SPECTRAL,
  cosine.threshold = COSINE_THRESHOLD
)
write.csv(spectra_raw, file.path(RUN_DIR, "spectra_automated_raw.csv"))

# Offsets vs the zip 1.8.8 references, script 74 metric ------------------------------------------------
s188 <- as.matrix(read.csv(
  file.path(
    DATA,
    "autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
  ),
  row.names = 1,
  check.names = FALSE
))
stopifnot("detectors align" = identical(colnames(spectra_raw), colnames(s188)))
b510 <- s188[TARGET, ]
offsets <- tibble(
  dye = intersect(rownames(spectra_raw), rownames(s188)),
  offset = vapply(
    dye,
    \(d) sum((spectra_raw[d, ] - s188[d, ]) * b510) / sum(b510^2),
    numeric(1)
  )
)
write_csv(offsets, file.path(RUN_DIR, "offsets.csv"))

# Figure: AF405 and SBUV510 offset per set ---------------------------------------------------------
plot_df <- read_csv(
  "experiments/runs/2026-11-10_bmv-variant-bv510-by-set/set_summary.csv",
  show_col_types = FALSE
) |>
  filter(set %in% c("zip_102", "chain4_spec", "chain5_bead"), dye %in% DYES) |>
  select(set, dye, offset) |>
  bind_rows(
    offsets |> filter(dye %in% DYES) |> mutate(set = "automated_raw")
  ) |>
  mutate(
    set = factor(
      set,
      levels = c("zip_102", "chain4_spec", "chain5_bead", "automated_raw")
    )
  )
p <- ggplot(plot_df, aes(set, offset)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_col(width = 0.6) +
  facet_wrap(~dye, nrow = 1) +
  labs(
    x = NULL,
    y = "BV510 coefficient of (reference - 1.8.8 zip reference)",
    title = "BV510 part of each set's AF405 and SBUV510 reference"
  ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(
  file.path(RUN_DIR, "fig_offsets.png"),
  p,
  width = 6,
  height = 3.2,
  dpi = 150,
  bg = "white"
)
write_csv(plot_df, file.path(RUN_DIR, "fig_offsets_source.csv"))

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
