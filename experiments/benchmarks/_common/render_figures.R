#!/usr/bin/env Rscript
# Purpose: build the comparison table and figures F1 to F4 from a finished run's saved outputs.
#          Sourced at the end of run_unmix.R; also runs alone to redo figures without recomputing.
# Run:     Rscript experiments/benchmarks/_common/render_figures.R <run_dir> <config.R> [output_folder_name]
#          output_folder_name defaults to "figures"; an existing folder is refused, so a redo uses a new name (for example figures_v2)
# Inputs:  <run_dir>/report/spectra.csv, <run_dir>/unmixed/<method>/*.fcs, the config, the dataset
# Outputs: <run_dir>/<output_folder_name>/ holding F1 to F4 (with source_data/), spectra_compare.csv,
#          spectra_unmatched.csv and sessionInfo.txt

if (!exists("cfg")) { # standalone use
  suppressPackageStartupMessages({ library(dplyr); library(tibble) })
  pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
  args <- commandArgs(trailingOnly = TRUE)
  run_dir <- args[[1]]
  CONFIG_PATH <- args[[2]]
  out_name <- if (length(args) >= 3) args[[3]] else "figures"
  cfg <- source(CONFIG_PATH, local = TRUE)$value
  data_dir <- Sys.getenv("SW_DATA_DIR", "/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata")
  control_dir <- file.path(data_dir, cfg$dataset, cfg$control_dir)
  spectroflo <- file.path(data_dir, cfg$dataset, cfg$spectroflo)
}
SEED <- cfg$seed
N_EVENTS <- cfg$n_events
PSEUDO_LOG_SIGMA <- cfg$pseudo_log_sigma
BIPLOT_BINS <- cfg$biplot_bins
BIPLOT_PAIR <- cfg$biplot_pair
if (!exists("out_name")) out_name <- "figures"
out_dir <- file.path(run_dir, out_name)
stopifnot("output folder already exists in this run" = !dir.exists(out_dir))
dir.create(out_dir, recursive = TRUE)

# Our reference spectra, as saved by the report (long form: fluorophore, detector, value)
spec_long <- read.csv(file.path(run_dir, "report", "spectra.csv"), check.names = FALSE)
fluors <- unique(spec_long$fluorophore)
dets <- unique(spec_long$detector)
our_spectra <- matrix(
  spec_long$value[match(paste(rep(fluors, times = length(dets)), rep(dets, each = length(fluors))),
                        paste(spec_long$fluorophore, spec_long$detector))],
  nrow = length(fluors), dimnames = list(fluors, dets)
)
unmixed_files <- list.files(file.path(run_dir, "unmixed"), pattern = "\\.fcs$", recursive = TRUE, full.names = TRUE)
names(unmixed_files) <- basename(dirname(unmixed_files))

# Comparator spectra from spectracle, same controls ----------------------------
sp <- sw_unmix_spectracle(control_dir, af = "separate", tissue_type = cfg$spectracle_tissue_type)
sp_spectra <- sp$spectra
# spectracle names detectors without the "-A" suffix AutoSpectral uses
to_as_name <- setNames(paste0(colnames(sp_spectra), "-A"), colnames(sp_spectra))
colnames(sp_spectra) <- unname(to_as_name[colnames(sp_spectra)])
# spectracle parses fluorophore names from file names; the config maps its names to ours
rn <- rownames(sp_spectra)
hit <- rn %in% names(cfg$spectracle_rename)
rn[hit] <- unname(cfg$spectracle_rename[rn[hit]])
rownames(sp_spectra) <- rn
unmatched <- tibble(
  fluorophore = union(
    setdiff(rownames(our_spectra), rownames(sp_spectra)),
    setdiff(rownames(sp_spectra), rownames(our_spectra))
  )
) |>
  mutate(
    in_ours = fluorophore %in% rownames(our_spectra),
    in_spectracle = fluorophore %in% rownames(sp_spectra)
  )
write.csv(unmatched, file.path(out_dir, "spectra_unmatched.csv"), row.names = FALSE)
shared_fluors <- intersect(rownames(our_spectra), rownames(sp_spectra))
shared_dets <- intersect(colnames(our_spectra), colnames(sp_spectra))
spectra_list <- list(
  ours = our_spectra[shared_fluors, shared_dets, drop = FALSE],
  spectracle = sp_spectra[shared_fluors, shared_dets, drop = FALSE]
)
spectra_compare <- sw_unmix_compare_spectra(spectra_list$ours, spectra_list$spectracle)
write.csv(spectra_compare, file.path(out_dir, "spectra_compare.csv"), row.names = FALSE)

# Figures F1 to F4 -------------------------------------------------------------
read_exprs <- function(path) {
  flowCore::exprs(flowCore::read.FCS(path, truncate_max_range = FALSE))
}
event_list <- c(lapply(unmixed_files, read_exprs), list(spectroflo = read_exprs(spectroflo)))
names(event_list) <- c(names(unmixed_files), "spectroflo")
common <- Reduce(intersect, lapply(event_list, colnames))
common <- common[grepl("-A$", common) & !grepl("^(FSC|SSC|Time)", common)]
events <- sw_unmix_sample_events(
  lapply(event_list, \(m) m[, common, drop = FALSE]),
  n_events = N_EVENTS,
  seed = SEED
)
# F3 shows the channels the config names (all common channels when it names none)
dist_channels <- if (is.null(cfg$distribution_channels)) common else intersect(cfg$distribution_channels, common)
stopifnot("F3 channels found in the sample" = length(dist_channels) > 0)
panels <- list(
  F1 = sw_unmix_figure_spectra(spectra_list),
  F2 = sw_unmix_figure_similarity(our_spectra),
  F3 = sw_unmix_figure_distributions(events, channels = dist_channels, pseudo_log_sigma = PSEUDO_LOG_SIGMA),
  F4 = sw_unmix_figure_biplot(
    events,
    x_channel = BIPLOT_PAIR[[1]],
    y_channel = BIPLOT_PAIR[[2]],
    pseudo_log_sigma = PSEUDO_LOG_SIGMA,
    bins = BIPLOT_BINS
  )
)
# F1 gets taller with the number of fluorophores (one facet each)
n_fluors <- nrow(spectra_list$ours)
sw_unmix_save_figures(
  panels,
  file.path(out_dir, "panels"),
  sizes = list(
    F1 = c(180, max(140, 20 * ceiling(sqrt(n_fluors)) + 40)),
    F2 = rep(max(120, 5 * nrow(our_spectra) + 50), 2) # one tile per fluorophore; labels need room
  )
)


writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
