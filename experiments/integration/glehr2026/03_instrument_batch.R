# experiments/integration/glehr2026/03_instrument_batch.R
# Purpose: real batch test. The same 9 BAD donors were stained once and acquired on an Aurora and a
#          Fortessa, so the instrument is a known batch with identical biology. Score each marker's
#          batch effect with sw_batch_scores() before CD3 relativisation, after the authors'
#          relativisation, and after ours (sw_relativise_fit + sw_relativise_apply, refitted per
#          instrument), and read the instrument gain spread from sw_relativise_diagnose().
# Run:     Rscript experiments/integration/glehr2026/03_instrument_batch.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*.fcs   CD3 controls, panel, authors' relativised panel
# Outputs: <run_dir>/batch_scores_long.csv   stage, check, marker, batch pair, value
#          <run_dir>/relativise_diagnose.csv  gain checks of our model with instrument as batch
#          <run_dir>/fig_instrument_eta2.png, figure_source_data_instrument_eta2.csv
#
# Choices approved by the user on 2026-11-10:
#   Q1 to Q4 as in 02_reproduce_cd3.R (gate box 1-99%, asinh(x / 6000) split, seed 1, nstart 10,
#   min_events 200); Q5 linear scale, ranked by eta2; Q6 2000 events per file, seed 1.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) stop("run dir exists: ", run_dir, call. = FALSE)
dir.create(run_dir, recursive = TRUE)

INSTRUMENTS <- c("aurora", "fortessa")
GATE_PROBS <- c(0.01, 0.99) # Q1
COFACTOR <- 6000 # Q2
SEED <- 1 # Q4, Q6
NSTART <- 10 # Q4
MIN_EVENTS <- 200 # Q4
N_EVENTS <- 2000 # Q6
LOW <- 0 # paper target for CD3- median
HIGH <- 1e4 # paper target for CD3+ median
KEY <- c("TIME", "FSC-A", "FSC-H", "SSC-A", "SSC-H")

read_exprs <- function(path) {
  exprs(read.FCS(path, truncate_max_range = FALSE, transformation = FALSE, emptyValue = FALSE))
}

# Read controls and panels per donor and instrument --------------------------------------------

donor_ids <- sort(basename(list.dirs(file.path(data_dir, "BAD", "aurora"), recursive = FALSE)))
units <- expand.grid(donor = seq_along(donor_ids), instrument = INSTRUMENTS, stringsAsFactors = FALSE)
units$unit <- paste0("d", units$donor, "_", units$instrument)

loaded <- lapply(seq_len(nrow(units)), function(j) {
  dir_j <- file.path(data_dir, "BAD", units$instrument[j], donor_ids[units$donor[j]])
  fs <- list.files(dir_j, pattern = "\\.fcs$", full.names = TRUE)
  rel <- read_exprs(fs[grepl("relativised", fs)])
  orig <- read_exprs(fs[grepl("_ungated_.*12-panel", fs)])
  # The authors' CD3+ events, before relativisation: matched on unchanged scatter and time.
  k_orig <- do.call(paste, as.data.frame(orig[, KEY]))
  k_rel <- do.call(paste, as.data.frame(rel[, KEY]))
  uniq <- !(k_orig %in% k_orig[duplicated(k_orig)])
  hit <- match(k_rel, k_orig[uniq])
  before <- orig[uniq, , drop = FALSE][hit[!is.na(hit)], ]
  box <- apply(rel[, c("FSC-A", "SSC-A")], 2, stats::quantile, probs = GATE_PROBS)
  ctrl_files <- fs[grepl("_ungated_.*_[0-9]{2}-CD3-", fs)]
  fluor <- sub(".*-CD3-(.*)\\.fcs$", "\\1", ctrl_files)
  controls <- lapply(ctrl_files, function(f) {
    e <- read_exprs(f)
    keep <- e[, "FSC-A"] >= box[1, 1] & e[, "FSC-A"] <= box[2, 1] &
      e[, "SSC-A"] >= box[1, 2] & e[, "SSC-A"] <= box[2, 2]
    e[keep, , drop = FALSE]
  })
  names(controls) <- paste0(units$unit[j], "_", fluor)
  list(
    before = before,
    authors = rel[!is.na(hit), ],
    controls = controls,
    map = data.frame(unit = units$unit[j], channel = paste0(fluor, "-A"), control = names(controls))
  )
})
names(loaded) <- units$unit
markers <- unique(unlist(lapply(loaded, \(x) x$map$channel)))

# Fit and apply our relativisation -------------------------------------------------------------

controls <- unlist(lapply(loaded, \(x) x$controls), recursive = FALSE)
names(controls) <- unlist(lapply(loaded, \(x) names(x$controls)))
control_map <- bind_rows(lapply(loaded, \(x) x$map))
positive <- sw_relativise_split_kmeans(controls, control_map, \(x) asinh(x / COFACTOR), SEED, NSTART)
model <- sw_relativise_fit(controls, control_map, positive, min_events = MIN_EVENTS)
before <- lapply(loaded, \(x) x$before)
ours <- sw_relativise_apply(
  before,
  model,
  unit = setNames(names(before), names(before)),
  channels = markers,
  low = LOW,
  high = HIGH,
  input_scale = "linear"
)
diag <- sw_relativise_diagnose(model, batch = setNames(units$instrument, units$unit), order = NULL)
write.csv(diag, file.path(run_dir, "relativise_diagnose.csv"), row.names = FALSE)

# Batch scores at each stage -------------------------------------------------------------------

batch <- setNames(units$instrument, units$unit)
stages <- list(before = before, authors = lapply(loaded, \(x) x$authors), ours = ours)
scores <- lapply(names(stages), function(s) {
  sw_batch_scores(stages[[s]], batch, markers, n_events = N_EVENTS, seed = SEED, order = NULL, reference_group = NULL, population = NULL) |>
    mutate(stage = s)
}) |>
  bind_rows()
write.csv(scores, file.path(run_dir, "batch_scores_long.csv"), row.names = FALSE)

eta <- scores |>
  filter(check == "eta2") |>
  mutate(stage = factor(stage, levels = names(stages)))
print(as.data.frame(summarise(eta, median_eta2 = median(value), max_eta2 = max(value), .by = stage)))
print(as.data.frame(tidyr::pivot_wider(select(eta, marker, stage, value), names_from = stage, values_from = value)))
print(as.data.frame(filter(diag, grepl("batch", check))))

# Figure ---------------------------------------------------------------------------------------

write.csv(eta, file.path(run_dir, "figure_source_data_instrument_eta2.csv"), row.names = FALSE)
p <- ggplot(eta, aes(value, marker, colour = stage)) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = c(before = "#B2182B", authors = "#4D4D4D", ours = "#2166AC")) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "eta2: share of file-median variance explained by instrument", y = NULL, colour = NULL)
ggsave(file.path(run_dir, "fig_instrument_eta2.png"), p, width = 120, height = 80, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
