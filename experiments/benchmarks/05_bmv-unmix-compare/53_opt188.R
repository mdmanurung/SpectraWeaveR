# Purpose: build the optimised AutoSpectral 1.8.8 arm (opt188) from the raw controls with the settings the user
#          approved on 2026-11-10 (README rows: bead pre-gate, Zombie NIR, cell landmark gate, Unstained in the
#          cell gate group, AF refine, af.remove FALSE, n.passes, stained.sample, speed slow, seeds), and unmix.
#          Two gate routes (user, 2026-11-10), because get.spectra.automated() takes no gate list:
#            pregate  the landmark polygon is applied to the stained cell controls (not the cell negative), gated copies
#                     are written, then get.spectra.automated(10000, 500, 0.95) and get.spectral.variants() run
#                     on them unchanged (option A)
#            legacy   define.flow.control(gate.list = list(lymphocyte = <landmark gate>)) -> clean.controls(
#                     af.remove = FALSE) -> get.fluorophore.spectra(); variants from the files without the
#                     landmark gate (option B)
#          Two cell-negative options, pending the chain6a/6b result (user, 2026-11-10):
#            neg_same_fsc4e5  Unstained (Cells) re-gated from raw with the wider cell rectangle
#            neg_h12          the H12 unstained sample, inside the same rectangle, as the cell negative
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/53_opt188.R <route> <neg> <n_passes> <stage>
#            from the repo root, on SLURM (experiments/slurm/run_opt188.sbatch); route = pregate or legacy;
#            neg = neg_same_fsc4e5 or neg_h12; stage = extract (controls, spectra, AF, variants, object check
#            only; the smoke test) or full (also unmix and verify).
#            extract writes to experiments/runs/2026-11-10_bmv-opt188-smoke/<arm>/, full to
#            experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/; either must not exist.
#            arm = opt188_pregate or opt188_legacyspec, then _<neg>_p<n_passes>.
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group/  raw controls and
#            unstained_AF/H12 [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/Zombie NIR (Cells).fcs
#            84,528-event file (user, 2026-11-10, D13)
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Raw/  the 5 samples [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/  143803 spectra, for the
#            object check [user zip]
#          experiments/runs/2026-11-10_bmv-cell-gating/fcs_control_file_landmark.csv  control-file copy with
#            gate.name lymphocyte for the cell controls and Unstained (Cells) [50_cell_gate_plots.R]
# Outputs: <out_dir>/controls/          control files used (re-gated) and the control-file copy
#          <out_dir>/work/              AutoSpectral tables and figures, incl. figure_gate/
#          <out_dir>/control_gate_log.csv  events kept per control and gate
#          <out_dir>/spectra.csv, af_spectra.csv, variants.rds, landmark_gate.rds  built objects
#          <out_dir>/object_check.csv   built objects vs the zip objects, per dye
#          <out_dir>/args.txt, sessionInfo.txt; full stage also *.fcs and verify.csv

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
SMOKE_ROOT <- "experiments/runs/2026-11-10_bmv-opt188-smoke"
PIPELINE <- "joint" # zip 1.8.8 script (inherited: user zip)
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint"
N_AF_PASSES <- 1L # default; the passes grid (52) may change n.passes and n.af.passes later
EXTRACT_SEED <- 42 # set.seed() around extraction (user, 2026-11-10, D5); inside, asp$bird.seed
BEAD_FSC <- c(2e5, 9e5) # widened bead rectangle (user, 2026-11-10, O1)
BEAD_SSC <- c(1.3e6, 2.25e6) # same
CELL_FSC <- c(4e5, 4e6) # wider cell rectangle (user, 2026-11-10, chain6, after 50's plots)
CELL_SSC <- c(1e5, 1.6e6) # zip SSC bounds kept (same)
LANDMARK_N_CELLS <- 2000 # define.gate.landmarks n.cells (user, 2026-11-10, from the tune.gate grid)
LANDMARK_PERCENTILE <- 70 # define.gate.landmarks percentile (same)
STAINED_SAMPLE <- "A3 BMV_009C_3dCD328_BMV_Aurora_3dCD328_Batch1.fcs" # get.spectral.variants stained.sample: donor A3 raw (user, 2026-11-10)
NECESSITY_THRESHOLD <- 0.01 # optimize.necessity.threshold, package default kept (user, 2026-11-10)
LANDMARK_FILE <- normalizePath(
  "experiments/runs/2026-11-10_bmv-cell-gating/fcs_control_file_landmark.csv"
)

# Setup -----------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
route <- cli_args[[1]]
neg <- cli_args[[2]]
n_passes <- as.integer(cli_args[[3]])
stage <- cli_args[[4]]
stopifnot(
  "route is pregate or legacy" = route %in% c("pregate", "legacy"),
  # neg_h12 dropped (user, 2026-11-10): the H12 file holds the same events as Unstained (Cells)
  "neg is neg_same_fsc4e5" = neg == "neg_same_fsc4e5",
  "stage is extract or full" = stage %in% c("extract", "full"),
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)
arm <- paste0(
  if (route == "pregate") "opt188_pregate" else "opt188_legacyspec",
  "_",
  neg,
  "_p",
  n_passes
)
out_dir <- if (stage == "extract") {
  file.path(normalizePath(SMOKE_ROOT, mustWork = FALSE), arm)
} else {
  file.path(DATA, "unmix_cells", arm)
}
stopifnot("out dir must not exist" = !dir.exists(out_dir))
ctrl_dir <- file.path(out_dir, "controls")
work_dir <- file.path(out_dir, "work")
dir.create(ctrl_dir, recursive = TRUE)
dir.create(work_dir)
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
raw_ctrl_dir <- file.path(d188, "Reference Group")
clean_dir <- file.path(d188, "Reference Group Clean")
h12_file <- file.path(
  raw_ctrl_dir,
  "unstained_AF",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs"
)
stained_file <- file.path(d188, "Raw", STAINED_SAMPLE)
stopifnot("stained sample exists" = file.exists(stained_file))
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
dir_fields <- names(asp)[grepl("dir", names(asp)) & map_lgl(asp, is.character)]
asp[dir_fields] <- map(asp[dir_fields], \(d) file.path(work_dir, d))
t0 <- Sys.time()

# Controls: rectangles on the raw files, control-file copy with the chosen cell negative --------------
# Beads: widened bead rectangle. Cells and the cell negative: wider cell rectangle. Zombie NIR: the
# 84,528-event clean file, copied as is. neg_h12 renames the H12 unstained into the Unstained (Cells) row,
# so H12 is the AF control and the universal negative for every cell control.
ctrl <- read.csv(LANDMARK_FILE, check.names = FALSE)
cell_negative <- "Unstained (Cells).fcs"
if (neg == "neg_h12") {
  h12_name <- "H12 Unstained (Cells).fcs"
  ctrl$filename[ctrl$filename == cell_negative] <- h12_name
  ctrl$universal.negative[ctrl$universal.negative == cell_negative] <- h12_name
  cell_negative <- h12_name
}
regate <- \(src, dest, fsc, ssc, kind) {
  ff <- read.FCS(src, transformation = FALSE, truncate_max_range = FALSE)
  kept <- Subset(
    ff,
    rectangleGate(filterId = kind, list("FSC-A" = fsc, "SSC-A" = ssc))
  )
  write.FCS(kept, file.path(ctrl_dir, dest))
  tibble(file = dest, gate = kind, n_in = nrow(ff), n_kept = nrow(kept))
}
bead_files <- ctrl$filename[ctrl$control.type == "beads"]
cell_files <- setdiff(
  ctrl$filename[ctrl$control.type == "cells"],
  c("Zombie NIR (Cells).fcs", cell_negative)
)
negative_src <- if (neg == "neg_h12") {
  h12_file
} else {
  file.path(raw_ctrl_dir, cell_negative)
}
gate_log <- bind_rows(
  map(bead_files, \(f) {
    regate(file.path(raw_ctrl_dir, f), f, BEAD_FSC, BEAD_SSC, "bead_rect")
  }),
  map(cell_files, \(f) {
    regate(file.path(raw_ctrl_dir, f), f, CELL_FSC, CELL_SSC, "cell_rect")
  }),
  regate(negative_src, cell_negative, CELL_FSC, CELL_SSC, "cell_rect")
)
file.copy(file.path(clean_dir, "Zombie NIR (Cells).fcs"), ctrl_dir)
control_file <- file.path(ctrl_dir, "fcs_control_file.csv")
write.csv(ctrl, control_file, row.names = FALSE)

# Landmark gate on the cell controls (gate.define TRUE for the six metabolic controls) ------------------
landmark_gate <- define.gate.landmarks(
  control.file = control_file,
  control.dir = ctrl_dir,
  asp = asp,
  gate.name = "lymphocyte",
  n.cells = LANDMARK_N_CELLS,
  percentile = LANDMARK_PERCENTILE,
  output.dir = file.path(work_dir, "figure_gate")
)
saveRDS(landmark_gate, file.path(out_dir, "landmark_gate.rds"))

# pregate: apply the landmark polygon to the lymphocyte-group files and keep the gated copies ----------
# legacy: keep a copy of the rectangle-only folder for the variants, gate inside define.flow.control
lymph_files <- ctrl$filename[ctrl$gate.name == "lymphocyte"]
variant_dir <- ctrl_dir
if (route == "legacy") {
  variant_dir <- file.path(out_dir, "controls_rect_only")
  dir.create(variant_dir)
  file.copy(list.files(ctrl_dir, full.names = TRUE), variant_dir)
}
if (route == "pregate") {
  # stained cell controls only; the cell negative keeps its rectangle gate (user, 2026-11-10): the polygon
  # left it 578 of 6,374 events, below AutoSpectral's 1,000-event floor
  landmark_log <- setdiff(lymph_files, cell_negative) |>
    map(\(f) {
      ff <- read.FCS(
        file.path(ctrl_dir, f),
        transformation = FALSE,
        truncate_max_range = FALSE
      )
      # AutoSpectral's own point-in-polygon test, carried through a row-index column
      indexed <- cbind(exprs(ff), .row = seq_len(nrow(ff)))
      kept_rows <- apply.gate(indexed, landmark_gate, asp = asp)[, ".row"]
      write.FCS(ff[kept_rows, ], file.path(ctrl_dir, f))
      tibble(
        file = f,
        gate = "landmark",
        n_in = nrow(ff),
        n_kept = length(kept_rows)
      )
    }) |>
    list_rbind()
  gate_log <- bind_rows(gate_log, landmark_log)
}
write_csv(gate_log, file.path(out_dir, "control_gate_log.csv"))

# Dye spectra --------------------------------------------------------------------------------------
set.seed(EXTRACT_SEED)
if (route == "pregate") {
  spectra <- get.spectra.automated(
    control.dir = ctrl_dir,
    control.def.file = control_file,
    asp = asp,
    n.candidates = 10000,
    n.spectral = 500,
    cosine.threshold = 0.95
  )
} else {
  # beads and the viability control keep density gates; the cell group uses the landmark gate
  flow_control_spec <- define.flow.control(
    ctrl_dir,
    control_file,
    asp,
    gate.list = list(lymphocyte = landmark_gate)
  )
  clean_control <- clean.controls(
    flow_control_spec,
    asp,
    universal.negative = TRUE,
    af.remove = FALSE,
    main.figures = FALSE
  )
  spectra <- get.fluorophore.spectra(
    clean_control,
    asp,
    use.clean.expr = TRUE,
    title = "opt188 legacy spectra",
    figures = FALSE
  )
  spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in the 06/07/16/26/49 scripts
}

# AF spectra and variants --------------------------------------------------------------------------
af_spectra <- get.af.spectra(
  unstained.sample = h12_file,
  asp = asp,
  spectra = spectra,
  refine = TRUE
)
spectra_variants <- get.spectral.variants(
  control.dir = variant_dir,
  control.def.file = file.path(variant_dir, "fcs_control_file.csv"),
  asp = asp,
  spectra = spectra,
  output.dir = work_dir,
  stained.sample = stained_file,
  optimize.necessity.threshold = NECESSITY_THRESHOLD
)
write.csv(spectra, file.path(out_dir, "spectra.csv"))
write.csv(af_spectra, file.path(out_dir, "af_spectra.csv"))
saveRDS(spectra_variants, file.path(out_dir, "variants.rds"))
extract_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

# Object check: built objects vs the zip objects ---------------------------------------------------
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra_zip <- read_raw(file.path(
  d188,
  "table_spectra",
  "Automated_autospectral_spectra_20261105_143803.csv"
))
variants_zip <- readRDS(file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
))
angle_deg <- \(a, b) {
  acos(min(1, sum(a * b) / sqrt(sum(a^2) * sum(b^2)))) * 180 / pi
}
stopifnot(
  "built fluors equal zip fluors" = setequal(
    rownames(spectra),
    rownames(spectra_zip)
  )
)
object_check <- tibble(dye = rownames(spectra_zip)) |>
  mutate(
    spectra_angle_vs_zip = map_dbl(dye, \(d) {
      angle_deg(spectra[d, colnames(spectra_zip)], spectra_zip[d, ])
    }),
    n_variants = map_int(dye, \(d) nrow(spectra_variants$variants[[d]])),
    n_variants_zip = map_int(dye, \(d) nrow(variants_zip$variants[[d]])),
    threshold = unname(spectra_variants$thresholds[dye]),
    threshold_zip = unname(variants_zip$thresholds[dye]),
    recommended = unname(spectra_variants$optimize.recommended[dye]),
    af_rows = nrow(af_spectra),
    extract_min = extract_min
  )
write_csv(object_check, file.path(out_dir, "object_check.csv"))

# Arguments record ---------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste("route:", route, " neg:", neg, " stage:", stage),
    paste(
      "bead rectangle FSC-A:",
      paste(BEAD_FSC, collapse = "-"),
      " SSC-A:",
      paste(BEAD_SSC, collapse = "-")
    ),
    paste(
      "cell rectangle FSC-A:",
      paste(CELL_FSC, collapse = "-"),
      " SSC-A:",
      paste(CELL_SSC, collapse = "-")
    ),
    paste("cell negative:", cell_negative, "from", negative_src),
    "Zombie NIR: 84,528-event Reference Group Clean file (user, 2026-11-10, D13)",
    "control file: fcs_control_file_landmark.csv copy (Unstained in gate.name lymphocyte)",
    paste(
      "landmark gate: define.gate.landmarks(n.cells =",
      LANDMARK_N_CELLS,
      ", percentile =",
      LANDMARK_PERCENTILE,
      ")"
    ),
    if (route == "pregate") {
      "spectra: landmark polygon applied to lymphocyte-group files; get.spectra.automated(10000, 500, 0.95)"
    } else {
      "spectra: define.flow.control(gate.list = lymphocyte landmark; density for other groups) -> clean.controls(universal.negative = TRUE, af.remove = FALSE) -> get.fluorophore.spectra(use.clean.expr = TRUE); variants from files without the landmark gate"
    },
    "get.af.spectra(H12, refine = TRUE)",
    paste(
      "get.spectral.variants(stained.sample =",
      STAINED_SAMPLE,
      ", optimize.necessity.threshold =",
      NECESSITY_THRESHOLD,
      ")"
    ),
    paste(
      "extraction seed (set.seed):",
      EXTRACT_SEED,
      " asp$bird.seed:",
      asp$bird.seed
    ),
    paste(
      "pipeline:",
      PIPELINE,
      " speed:",
      SPEED,
      " n.passes:",
      n_passes,
      " n.af.passes:",
      N_AF_PASSES
    ),
    paste("AF rows:", nrow(af_spectra), " extraction minutes:", extract_min)
  ),
  file.path(out_dir, "args.txt")
)

# Unmix and verify (full stage only) -------------------------------------------------------------
if (stage == "full") {
  flow.control <- reload.flow.control(ctrl_dir, control_file, asp)
  t1 <- Sys.time()
  unmix.folder(
    fcs.dir = file.path(d188, "Raw"),
    spectra = spectra,
    asp = asp,
    flow.control = flow.control,
    method = "AutoSpectral",
    af.spectra = af_spectra,
    spectra.variants = spectra_variants,
    pipeline = PIPELINE,
    speed = SPEED,
    n.passes = n_passes,
    n.af.passes = N_AF_PASSES,
    output.dir = out_dir,
    parallel = TRUE,
    verbose = TRUE
  )
  message(
    "unmix minutes: ",
    round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1)
  )
  ref_dir <- file.path(d188, "AutoSpectral_unmixed")
  read_hdr <- \(path) {
    ff <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
    list(n = nrow(ff), channels = colnames(ff))
  }
  verify <- tibble(file = list.files(ref_dir, pattern = "\\.fcs$")) |>
    mutate(
      ref = map(file.path(ref_dir, file), read_hdr),
      new = map(file.path(out_dir, file), read_hdr),
      n_ref = map_int(ref, "n"),
      n_new = map_int(new, "n"),
      channels_equal = map2_lgl(ref, new, \(r, n) {
        identical(r$channels, n$channels)
      })
    ) |>
    select(file, n_ref, n_new, channels_equal)
  write_csv(verify, file.path(out_dir, "verify.csv"))
  stopifnot(
    "5 files unmixed" = nrow(verify) == 5,
    "event counts equal as188_joint" = all(verify$n_ref == verify$n_new),
    "channel names equal as188_joint" = all(verify$channels_equal)
  )
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
