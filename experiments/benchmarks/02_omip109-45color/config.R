# Settings for benchmark 02 (OMIP-109, 45 colours). Methods, seed, sample size,
# sigma and bins are carried over from benchmark 01 and are NOT yet confirmed by
# the user for this dataset; change them here.
list(
  dataset = "omip109",
  control_dir = "data/Raw/Reference Group",
  control_file = NULL, # AutoSpectral builds the control file from the file names (46 controls found)
  unstained = "data/Raw/Reference Group/Unstained (Cells).fcs",
  stained = "data/Raw/Multicolor/45-Color Donor 1.fcs",
  spectroflo = "data/Unmixed/Multicolor/45-Color Donor 1.fcs",
  cytometer = "aurora",
  spectracle_tissue_type = "Cells",
  # spectracle's names -> ours, matched by name from the first run's spectra_unmatched.csv (to be confirmed by the user)
  spectracle_rename = c(
    "Spark UV387" = "Spark UV 387", "DEAD Blue" = "LIVE DEAD Blue", "BUV805 +BSB" = "BUV805",
    "RealBlue 780" = "RB780", "PE-Dazzle594" = "PE-Dazzle 594", "PE-Fire 640 250ng" = "PE-Fire 640",
    "APC 500ng" = "APC", "Spark NIR685" = "Spark NIR 685"
  ),
  methods = c("ols", "autospectral_fast", "autospectral_slow"),
  limits = sw_unmix_limits_autospectral(), # SIF warn 4, fail 6; user choice 2026-11-09
  seed = 1,
  n_events = 50000,
  pseudo_log_sigma = 1000,
  biplot_bins = 200,
  distribution_channels = c("BUV395-A", "BV421-A", "BV480-A", "Super Bright 436-A"), # F3 subset for a 40+ colour panel; not yet confirmed by the user
  biplot_pair = c("BUV395-A", "BV421-A"), # to be checked against this panel's detectors
  title = "OMIP-109 45-colour PBMC: SpectraWeaveR beta run",
  notes = "Real-data run on one stained file. Comparators: spectracle spectra and SpectroFlo unmixed output. Limits: SIF warn 4, fail 6 (AutoSpectral article 01, chosen by the user 2026-11-09); every other metric is not judged."
)
