# experiments/integration/cytobatchflagr/02_batch_cd3.R
# Purpose: real batch series with constant biology. The same four control samples were acquired in
#          18 batches (dates). Score each marker's batch effect with sw_batch_scores(), fit an
#          internal CD3 reference per file (each file's own CD3 split, sw_relativise_internal_map),
#          read the CD3 gain across batches and dates with sw_relativise_diagnose(), and compare
#          the CD3 channel's batch eta2 before and after relativisation.
# Run:     Rscript experiments/integration/cytobatchflagr/02_batch_cd3.R <extract_run_dir> <run_dir>
# Inputs:  <extract_run_dir>/events.rds, files.csv, marker_names.rds   from 01_extract.R
# Outputs: <run_dir>/batch_scores_long.csv    stage, check, marker, batch, value
#          <run_dir>/cd3_model.csv            per file: m, M, gain, separation, batch, sample, date
#          <run_dir>/relativise_diagnose.csv  gain checks with batch and date order
#          <run_dir>/fig_cbf_eta2.png, figure_source_data_cbf_eta2.csv
#          <run_dir>/fig_cbf_cd3_gain.png, figure_source_data_cbf_cd3_gain.csv
#
# Choices approved by the user on 2026-11-10 (reuse of the Glehr choices): k-means split on
# asinh(x / 6000), seed 1, nstart 10, min_events 200, low 0, high 1e4, batch scores on the linear
# scale ranked by eta2, 2000 events per file, seed 1; the CD3 split uses all events (no gate);
# batch scores use the 30 markers shared by every file.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
extract_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) stop("run dir exists: ", run_dir, call. = FALSE)
dir.create(run_dir, recursive = TRUE)

CD3_CHANNEL <- "BUV563-A"
COFACTOR <- 6000
SEED <- 1
NSTART <- 10
MIN_EVENTS <- 200
N_EVENTS <- 2000
LOW <- 0
HIGH <- 1e4

ev <- readRDS(file.path(extract_dir, "events.rds"))
info <- read.csv(file.path(extract_dir, "files.csv"))
marker_names <- readRDS(file.path(extract_dir, "marker_names.rds"))
stopifnot("files match events" = identical(info$file, names(ev)))
ids <- paste0("f", seq_along(ev))
names(ev) <- ids
info$id <- ids
# Markers present in every file (user decision 2026-11-10): batches 5-18 lack some of CD27
# (BV570-A), CD28 (BV785-A) and CD163 (PerCP-A), so those three are left out.
shared <- Reduce(intersect, lapply(ev, colnames))
markers <- intersect(names(marker_names)[!is.na(marker_names)], shared)
batch <- setNames(info$batch, ids)
date_order <- setNames(as.numeric(as.Date(info$date, format = "%d-%b-%Y")), ids)

# Batch scores before relativisation ------------------------------------------------------------

before <- sw_batch_scores(ev, batch, markers, n_events = N_EVENTS, seed = SEED, order = NULL, reference_group = NULL, population = NULL)

# Internal CD3 reference per file ---------------------------------------------------------------

cd3_map <- sw_relativise_internal_map(ev, CD3_CHANNEL)
positive <- sw_relativise_split_kmeans(ev, cd3_map, \(x) asinh(x / COFACTOR), SEED, NSTART)
model <- sw_relativise_fit(ev, cd3_map, positive, min_events = MIN_EVENTS)
diag <- sw_relativise_diagnose(model, batch = batch, order = date_order)
write.csv(diag, file.path(run_dir, "relativise_diagnose.csv"), row.names = FALSE)
model_tbl <- model |>
  as_tibble() |>
  left_join(
    select(info, unit = id, batch, sample = control, date),
    by = join_by(unit),
    relationship = "one-to-one",
    unmatched = "error"
  )
write.csv(model_tbl, file.path(run_dir, "cd3_model.csv"), row.names = FALSE)

after_data <- sw_relativise_apply(
  ev,
  model,
  unit = setNames(ids, ids),
  channels = CD3_CHANNEL,
  low = LOW,
  high = HIGH,
  input_scale = "linear"
)
after <- sw_batch_scores(after_data, batch, CD3_CHANNEL, n_events = N_EVENTS, seed = SEED, order = NULL, reference_group = NULL, population = NULL)
scores <- bind_rows(mutate(before, stage = "before"), mutate(after, stage = "after CD3"))
write.csv(scores, file.path(run_dir, "batch_scores_long.csv"), row.names = FALSE)

# Summaries -------------------------------------------------------------------------------------

eta <- scores |>
  filter(check == "eta2") |>
  mutate(marker_name = unname(marker_names[marker]))
print(as.data.frame(arrange(select(eta, stage, marker_name, value), desc(value))))
print(as.data.frame(filter(diag, check %in% c("gain_order_rho"))))
print(as.data.frame(summarise(model_tbl, gain_cv = sd(gain) / mean(gain), sep_min = min(separation),
  sep_median = median(separation), .by = sample)))

# Figures ---------------------------------------------------------------------------------------

eta_plot <- eta |>
  mutate(marker_name = forcats::fct_reorder(marker_name, value, .fun = max))
write.csv(eta_plot, file.path(run_dir, "figure_source_data_cbf_eta2.csv"), row.names = FALSE)
p1 <- ggplot(eta_plot, aes(value, marker_name, colour = stage)) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c(before = "#B2182B", `after CD3` = "#2166AC")) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "eta2: share of file-median variance explained by batch (18 batches)", y = NULL, colour = NULL)
ggsave(file.path(run_dir, "fig_cbf_eta2.png"), p1, width = 120, height = 140, units = "mm", dpi = 200)

gain_plot <- model_tbl |>
  mutate(date = as.Date(date, format = "%d-%b-%Y"))
write.csv(gain_plot, file.path(run_dir, "figure_source_data_cbf_cd3_gain.csv"), row.names = FALSE)
p2 <- ggplot(gain_plot, aes(date, gain, colour = sample)) +
  geom_point(size = 1.2) +
  geom_line(alpha = 0.4) +
  scale_y_log10() +
  labs(x = "Acquisition date (one batch per date)", y = "CD3 gain M - m (linear units, log scale)", colour = NULL)
ggsave(file.path(run_dir, "fig_cbf_cd3_gain.png"), p2, width = 140, height = 75, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
