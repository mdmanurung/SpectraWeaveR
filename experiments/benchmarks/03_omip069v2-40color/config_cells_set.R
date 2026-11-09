# Settings for benchmark 03 (OMIP-069v2, 40 colours), cells spectra set from one combined control file.
# The control file holds cell and bead controls (AutoSpectral article 02b naming). Both spectra sets are unmixed, one run each.
list(
  dataset = "omip069v2",
  control_dir = "data/Raw/96 U-Bottom/Reference Group",
  control_file_repo = "experiments/benchmarks/03_omip069v2-40color/control_combined.csv", # one control file: cell rows plus bead rows named "<dye> beads" (AutoSpectral article 02b)
  spectra_set = "cells", # cells: cell spectra only; beads: bead spectrum where one exists, else the cell spectrum
  unstained = "data/Raw/96 U-Bottom/Reference Group/D5 Unstained (Cells)_96 U-Bottom.fcs",
  stained = "data/Raw/Multicolor/40-Color Revised Panel Donor 2.fcs",
  spectroflo = "data/Unmixed/Multicolor/40-Color Revised Panel Donor 2.fcs",
  cytometer = "aurora",
  spectracle_tissue_type = "Cells",
  # spectracle's names -> ours, matched by name from the first run's spectra_unmatched.csv (to be confirmed by the user)
  spectracle_rename = c("DEAD Blue" = "LIVE DEAD Blue", "PE-Dazzle594" = "PE-Dazzle 594", "Spark NIR685" = "Spark NIR 685"),
  methods = c("ols", "autospectral_fast", "autospectral_slow"),
  limits = sw_unmix_limits_autospectral(), # SIF warn 4, fail 6; user choice 2026-11-09
  seed = 1,
  n_events = 50000,
  pseudo_log_sigma = 1000,
  biplot_bins = 200,
  distribution_channels = c("BUV395-A", "BV421-A", "BV480-A", "Super Bright 436-A"), # F3 subset for a 40+ colour panel; not yet confirmed by the user
  biplot_pair = c("BUV395-A", "BV421-A"),
  title = "OMIP-069v2 40-colour PBMC, cells spectra set: SpectraWeaveR beta run",
  notes = "Real-data run on one stained file. Comparators: spectracle spectra and SpectroFlo unmixed output. Limits: SIF warn 4, fail 6 (AutoSpectral article 01, chosen by the user 2026-11-09); every other metric is not judged."
)
