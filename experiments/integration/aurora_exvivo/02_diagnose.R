# experiments/integration/aurora_exvivo/02_diagnose.R
# Purpose: run the technical diagnostics (header, signal, anomaly, batch, cluster) on the
#          subsampled real files and draw four figures. Numbers only: nothing is called
#          bad. Every setting is an illustrative choice for a first check.
# Run:     Rscript experiments/integration/aurora_exvivo/02_diagnose.R <fcs_dir> <run_dir>
# Inputs:  <fcs_dir>/*.fcs                 headers are read again for the header audit
#          <run_dir>/events.rds, headers.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/diagnosis.csv         sw_diagnose() table, one row per subject and check
#          <run_dir>/top5.csv              five largest values per check
#          <run_dir>/header_differences.csv rows where a header item differs from the reference
#          <run_dir>/fig_batch_pca.png     PCA of per-file features, coloured by batch
#          <run_dir>/fig_batch_markers.png per-marker eta2 and largest pairwise EMD
#          <run_dir>/fig_outlier_scores.png forest score per file, by batch and acquisition order
#          <run_dir>/fig_header.png        header differences per file, by batch
#          <run_dir>/figure_source_data_*.csv  the numbers behind each figure

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
fcs_dir <- args[[1]]
run_dir <- args[[2]]

# Settings (illustrative) ----------------------------------------------------
N_EVENTS <- 2000 # events per file used by every module
SEED <- 1
QUANTILE_STEP <- 0.1
AGG_SIZE <- 20000
FOREST_TREES <- 300
FOREST_SAMPLE <- 64
CLUSTER_K <- 12
CLUSTER_NSTART <- 3
MIN_CLUSTER_EVENTS <- 30
TIME_BINS <- 50 # time bins per file; the width is the median file span / TIME_BINS
SCALE_NOTE <- "as exported by OMIQ: marker channels already transformed (cofactor not recorded in the file)"
# Reference file for the header audit: the first file of the first batch.

# Load ------------------------------------------------------------------------
events <- readRDS(file.path(run_dir, "events.rds"))
headers <- readRDS(file.path(run_dir, "headers.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(events)))
channels <- setdiff(colnames(events[[1]]), "Time")
batch <- setNames(info$batch, info$file)

# Acquisition order from $DATE and $BTIM (within and across batches).
stamp <- vapply(
  headers,
  \(h) paste(h[["$DATE"]], h[["$BTIM"]]),
  character(1)
)
acq <- as.numeric(as.POSIXct(stamp, format = "%d-%b-%Y %H:%M:%S", tz = "UTC"))
names(acq) <- names(headers)
stopifnot("acquisition time parsed" = !anyNA(acq))

span <- vapply(events, \(m) diff(range(m[, "Time"])), numeric(1))
time_bin <- as.numeric(median(span) / TIME_BINS)

# Diagnose --------------------------------------------------------------------
paths <- file.path(fcs_dir, info$file)
reference <- paths[info$batch == sort(unique(info$batch))[1]][1]
dx <- sw_diagnose(
  events,
  header = list(
    paths = paths,
    keywords = c("$CYT", "$CYTSN", "$INST", "CREATOR", "$PnV"),
    reference = reference
  ),
  signal = list(
    channels = channels,
    range_max = NULL, # instrument maximum on the transformed scale is not known
    range_min = NULL,
    time_channel = "Time",
    time_bin = time_bin
  ),
  anomaly = list(
    channels = channels,
    features = c("quantile", "emd"),
    quantile_step = QUANTILE_STEP,
    n_events = N_EVENTS,
    seed = SEED,
    scale = SCALE_NOTE,
    aggregate_from = "all",
    agg_size = AGG_SIZE,
    ntrees = FOREST_TREES,
    sample_size = FOREST_SAMPLE,
    ndim = 1
  ),
  batch = list(
    batch = batch,
    markers = channels,
    n_events = N_EVENTS,
    seed = SEED,
    order = acq
  ),
  cluster = list(
    markers = channels,
    k = CLUSTER_K,
    n_events = N_EVENTS,
    seed = SEED,
    nstart = CLUSTER_NSTART,
    min_cluster_events = MIN_CLUSTER_EVENTS
  ),
  relativise = NULL, # no reference controls in these files
  otd = NULL,
  negative = NULL,
  gate = NULL
)
write.csv(dx, file.path(run_dir, "diagnosis.csv"), row.names = FALSE)
write.csv(
  sw_diagnosis_top(dx, n = 5),
  file.path(run_dir, "top5.csv"),
  row.names = FALSE
)
hd <- dx[dx$module == "header" & dx$value == 1, ]
write.csv(hd, file.path(run_dir, "header_differences.csv"), row.names = FALSE)

# Figure 1: PCA of the per-file features, coloured by batch ------------------
feat <- sw_anomaly_features(
  events,
  channels,
  c("quantile", "emd"),
  QUANTILE_STEP,
  N_EVENTS,
  SEED,
  SCALE_NOTE,
  "all",
  AGG_SIZE
)
# PCA is for display only; it is run on the unscaled feature columns.
pc <- prcomp(as.matrix(feat[-1]), center = TRUE, scale. = FALSE)
pca_df <- tibble(
  file = feat$file,
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  batch = batch[feat$file]
)
var_expl <- 100 * pc$sdev^2 / sum(pc$sdev^2)
write.csv(
  pca_df,
  file.path(run_dir, "figure_source_data_pca.csv"),
  row.names = FALSE
)
p1 <- ggplot(pca_df, aes(PC1, PC2, colour = batch)) +
  geom_point(size = 1.4, alpha = 0.85) +
  labs(
    x = sprintf("PC1 (%.0f%% of feature variance)", var_expl[1]),
    y = sprintf("PC2 (%.0f%%)", var_expl[2]),
    colour = NULL
  ) +
  coord_fixed()
ggsave(
  file.path(run_dir, "fig_batch_pca.png"),
  p1,
  width = 89,
  height = 80,
  units = "mm",
  dpi = 200
)

# Figure 2: per-marker batch scores -------------------------------------------
bs <- dx[dx$module == "batch", ]
eta <- bs |> filter(check == "eta2") |> select(marker = item, eta2 = value)
emd <- bs |>
  filter(check == "pair_emd") |>
  summarise(max_pair_emd = max(value), .by = item) |>
  rename(marker = item)
spr <- bs |>
  filter(check == "spread_ratio") |>
  summarise(max_abs_log_spread = max(abs(log(value))), .by = item) |>
  rename(marker = item)
marker_df <- eta |>
  left_join(emd, by = "marker", relationship = "one-to-one") |>
  left_join(spr, by = "marker", relationship = "one-to-one") |>
  mutate(marker = reorder(marker, eta2))
write.csv(
  marker_df,
  file.path(run_dir, "figure_source_data_markers.csv"),
  row.names = FALSE
)
long2 <- marker_df |>
  pivot_longer(c(eta2, max_pair_emd, max_abs_log_spread), names_to = "score")
p2 <- ggplot(long2, aes(value, marker)) +
  geom_point(size = 1.2) +
  facet_wrap(vars(score), scales = "free_x", nrow = 1) +
  labs(x = NULL, y = NULL)
ggsave(
  file.path(run_dir, "fig_batch_markers.png"),
  p2,
  width = 140,
  height = 130,
  units = "mm",
  dpi = 200
)

# Figure 3: outlier score per file, by batch and acquisition order -----------
out_df <- dx |>
  filter(check == "outlier_score") |>
  transmute(
    file = subject,
    score = value,
    batch = batch[subject],
    acquired = acq[subject]
  ) |>
  arrange(acquired) |>
  mutate(order = row_number())
write.csv(
  out_df,
  file.path(run_dir, "figure_source_data_outlier.csv"),
  row.names = FALSE
)
p3 <- ggplot(out_df, aes(order, score, colour = batch)) +
  geom_point(size = 1.2) +
  labs(
    x = "Acquisition order (all files)",
    y = "Isolation-forest score",
    colour = NULL
  )
ggsave(
  file.path(run_dir, "fig_outlier_scores.png"),
  p3,
  width = 120,
  height = 70,
  units = "mm",
  dpi = 200
)

# Figure 4: header differences per file --------------------------------------
hd_df <- dx |>
  filter(module == "header") |>
  summarise(n_differing = sum(value == 1), .by = subject) |>
  mutate(batch = batch[subject]) |>
  arrange(batch, subject)
write.csv(
  hd_df,
  file.path(run_dir, "figure_source_data_header.csv"),
  row.names = FALSE
)
p4 <- ggplot(hd_df, aes(batch, n_differing)) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0, seed = 1),
    size = 1
  ) +
  labs(x = NULL, y = "Header items that differ from the reference file")
ggsave(
  file.path(run_dir, "fig_header.png"),
  p4,
  width = 89,
  height = 70,
  units = "mm",
  dpi = 200
)

writeLines(
  c(
    paste("time_bin used:", time_bin),
    paste("reference file:", basename(reference)),
    paste("channels:", length(channels), "; files:", length(events)),
    capture.output(sessionInfo())
  ),
  file.path(run_dir, "run_notes.txt")
)
cat("diagnosis rows:", nrow(dx), "; modules:", attr(dx, "modules"), "\n")
