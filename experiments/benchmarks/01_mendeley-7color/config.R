# Settings for benchmark 01. Choices marked "user" were set by the user on 2026-11-09.
list(
  dataset = "autospectral_mendeley", # key in experiments/datasets/registry.csv
  control_dir = "data/Raw/Set1/Reference Group",
  control_file = "fcs_control_file.csv", # shipped with the data
  unstained = "data/Raw/Set1/Unstained/D2 Lung_Set1.fcs", # matched to the stained lung file, as in the article
  stained = "data/Raw/Set1/Stained/D5 Lung_Set1.fcs",
  spectroflo = "data/Unmixed/Set1/Stained/D5 Lung_Set1.fcs", # SpectroFlo comparator
  cytometer = "aurora", # Cytek Aurora, 5 laser
  spectracle_tissue_type = "Cells",
  spectracle_rename = c("F PE" = "PE", "e780" = "eFluor 780"), # spectracle names parsed from file names
  methods = c("ols", "autospectral_fast", "autospectral_slow"), # user
  limits = sw_unmix_limits_autospectral(), # SIF warn 4, fail 6; user choice 2026-11-09
  seed = 1, # user
  n_events = 50000, # per file, F3 and F4 (user)
  pseudo_log_sigma = 1000, # F3 and F4 (user)
  biplot_bins = 200, # F4 (user)
  biplot_pair = c("BUV395-A", "BV421-A"), # pair used in AutoSpectral article 01
  title = "Mendeley 7-colour mouse lung: SpectraWeaveR beta run",
  notes = "Real-data run on one stained lung file. Comparators: spectracle spectra and SpectroFlo unmixed output. Limits: SIF warn 4, fail 6 (AutoSpectral article 01, chosen by the user 2026-11-09); every other metric is not judged."
)
