# Purpose: compare the four unmixing arms downstream on the shared live singlet events: SOM
#          metaclusters (ARI per donor against the reference arm) and cytokine+ % in CD4+ and
#          CD8+ T cells. Events are paired across arms, so each arm is clustered on the same
#          events (all samples pooled) and its metaclusters are compared row by row.
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/03_downstream.R \
#            experiments/runs/2026-11-09_bmv-unmix-compare [smoke_n_events_per_sample]
#          The run dir already exists (made by 01); this script only adds files to it.
# Inputs:  <run_dir>/manifest.csv       arm, dir, viability, file, sample [01_manifest-gate.R]
#          <run_dir>/shared_events.rds  per sample: logical over ALL events, TRUE = shared live
#                                       singlet event [01_manifest-gate.R]
#          the FCS files listed in manifest.csv
# Outputs: <run_dir>/cluster_assignments.rds  arm, sample, event (row in the FCS file), metacluster
#          <run_dir>/clusters_ari.csv         arm, sample, ari against the reference arm
#          <run_dir>/cluster_freq.csv         arm, sample, metacluster, n events
#          <run_dir>/cytokine_freq.csv        arm, sample, parent, cytokine, pct_positive, threshold
#          <run_dir>/lineage_thresholds.csv   arm, sample, marker, threshold (CD3, CD4, CD8 cuts)
#          <run_dir>/sessionInfo_downstream.txt

suppressPackageStartupMessages({
  library(flowCore) # attached first: openCyto gates call exprs() unqualified; dplyr then wins masks
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE) # tidyverse and here are not in the pixi env; run from the repo root

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor for spectral flow (CLAUDE.md, confirmed by user)
REF_ARM <- "as188_joint" # reference arm (user, 2026-11-09)
XDIM <- 10 # P10: SOM grid, same in every arm
YDIM <- 10
N_META <- 20 # P10: metaclusters
SOM_SEED <- 42 # P10
SMOKE_SEED <- 42 # only used to pick the smoke-test subsample
NOT_CLUSTERED <- c("AF-A", "AF Index") # autofluorescence channels; viability and scatter drop out below
LINEAGE <- c(CD3 = "SBUV510-A", CD4 = "Spark Blue 574-A", CD8 = "SBUV445-A") # P11, checked vs FCS headers
MINDENSITY_RANGE <- c(0, 3) # mindensity window for the CD3, CD4, CD8 lineage cuts (user, 2026-11-09)
CYTOKINE_RANGE <- c(1, 4) # mindensity window for cytokine cuts, per arm (user, grill Q1, 2026-11-09): TNF pinned at 0 with c(0, 3)
NOISE_SEEDS <- c(1, 2, 3) # ARI noise floor: as188_joint reclustered with these seeds (user, grill Q2 and Q4, 2026-11-09)
# Clustering markers: lineage only (user, 2026-11-09), so clusters reflect cell type, not stimulation state.
CLUSTER_MARKERS <- c(
  CD3 = "SBUV510-A",
  CD4 = "Spark Blue 574-A",
  CD8 = "SBUV445-A",
  CD45RA = "Spark UV 387-A",
  CCR7 = "PE-Fire 810-A",
  CD127 = "APC-Fire 810-A",
  CD25 = "SBUV605-A",
  CD56 = "RB744-A",
  TCRgd = "BUV395-A",
  CD161 = "BUV563-A",
  CD57 = "BV785-A",
  CD38 = "SBUV740-A",
  CD7 = "BUV805-A",
  CXCR5 = "BV750-A"
)
CYTOKINES <- c(TNF = "BV510-A", IFNg = "BV570-A", `IL-2` = "BV605-A") # P11, checked vs FCS headers

# Run dir and smoke-test size -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[[1]]
SMOKE_N <- as.integer(args[2]) # NA = use every shared event; set only for the smoke test
stopifnot(
  "outputs must not exist" = !file.exists(file.path(
    run_dir,
    "cluster_assignments.rds"
  ))
)

# Inputs ----------------------------------------------------------------------------------
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
shared <- readRDS(file.path(run_dir, "shared_events.rds"))
samples <- sort(unique(manifest$sample))
arms <- unique(manifest$arm) |> set_names()

# Event indices to use: the same rows in every arm because events are paired.
set.seed(SMOKE_SEED)
keep <- samples |>
  set_names() |>
  map(\(s) {
    idx <- which(shared[[s]])
    if (is.na(SMOKE_N)) idx else sort(sample(idx, min(SMOKE_N, length(idx))))
  })
events <- tibble(
  sample = rep(samples, lengths(keep)),
  event = unlist(keep, use.names = FALSE)
)

# Channels: clustered detectors, taken from the reference arm, joined by detector name -----
ref_ff <- manifest |>
  filter(arm == REF_ARM) |>
  slice(1) |>
  pull(file) |>
  sw_io_read_fcs(which.lines = 10)
CHANNELS <- ref_ff |>
  sw_channel_is_fluor() |>
  keep(isTRUE) |>
  names() |>
  setdiff(c(
    NOT_CLUSTERED,
    manifest |> filter(arm == REF_ARM) |> pull(viability) |> unique()
  ))
stopifnot(
  "gate channels exist in the FCS files" = all(
    c(LINEAGE, CYTOKINES, CLUSTER_MARKERS) %in% CHANNELS
  )
)

# Arcsinh matrix per arm: rows = shared events (samples pooled), columns = CHANNELS -----------
arm_mats <- arms |>
  map(\(a) {
    files <- manifest |>
      filter(arm == a) |>
      arrange(sample)
    stopifnot(
      "manifest samples match shared events" = identical(files$sample, samples)
    )
    map2(files$file, keep, \(f, idx) {
      flowCore::exprs(sw_io_read_fcs(f))[idx, CHANNELS]
    }) |>
      do.call(what = rbind) |>
      (\(m) asinh(m / COFACTOR))()
  })

# Cluster each arm with identical settings ------------------------------------------------
clusters <- arm_mats |>
  map(\(m) {
    sw_cluster_run(
      m,
      lineage_markers = CLUSTER_MARKERS,
      method = "som",
      xdim = XDIM,
      ydim = YDIM,
      n_metaclusters = N_META,
      seed = SOM_SEED
    )
  })
assignments <- clusters |>
  imap(\(res, arm) mutate(events, arm = arm, metacluster = res$assignments)) |>
  list_rbind()
saveRDS(assignments, file.path(run_dir, "cluster_assignments.rds"))

# Adjusted Rand index (contingency formula) -----------------------------------------------
# aricode and mclust are not installed. Tested below on identical and shuffled labels.
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  choose2 <- \(x) x * (x - 1) / 2
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  expected <- sum_rows * sum_cols / choose2(length(a))
  (sum_cells - expected) / ((sum_rows + sum_cols) / 2 - expected)
}
set.seed(SMOKE_SEED)
labels_test <- sample(N_META, 5000, replace = TRUE)
stopifnot(
  "ari of identical labels is 1" = adjusted_rand_index(
    labels_test,
    labels_test
  ) ==
    1,
  "ari of shuffled labels is near 0" = abs(adjusted_rand_index(
    labels_test,
    sample(labels_test)
  )) <
    0.01
)

# ARI per arm and donor against the reference arm, on the same events ---------------------
ref_labels <- assignments |>
  filter(arm == REF_ARM) |>
  select(sample, event, ref_metacluster = metacluster)
clusters_ari <- assignments |>
  left_join(
    ref_labels,
    by = join_by(sample, event),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  summarise(
    ari = adjusted_rand_index(metacluster, ref_metacluster),
    .by = c(arm, sample)
  )
write_csv(clusters_ari, file.path(run_dir, "clusters_ari.csv"))

# ARI noise floor: the reference arm reclustered with other seeds ------------------------
# Same events and settings as the seed-42 run; only the SOM seed changes.
noise_labels <- NOISE_SEEDS |>
  set_names() |>
  map(\(seed) {
    sw_cluster_run(
      arm_mats[[REF_ARM]],
      lineage_markers = CLUSTER_MARKERS,
      method = "som",
      xdim = XDIM,
      ydim = YDIM,
      n_metaclusters = N_META,
      seed = seed
    )$assignments
  })
noise_labels[[as.character(SOM_SEED)]] <- clusters[[REF_ARM]]$assignments
seed_pairs <- combn(names(noise_labels), 2, simplify = FALSE)
ari_noise_floor <- seed_pairs |>
  map(\(pr) {
    events |>
      mutate(a = noise_labels[[pr[[1]]]], b = noise_labels[[pr[[2]]]]) |>
      summarise(ari = adjusted_rand_index(a, b), .by = sample) |>
      mutate(seed_a = pr[[1]], seed_b = pr[[2]])
  }) |>
  list_rbind()
write_csv(ari_noise_floor, file.path(run_dir, "ari_noise_floor.csv"))

# Events per metacluster, arm and donor ---------------------------------------------------
cluster_freq <- assignments |>
  count(arm, sample, metacluster)
write_csv(cluster_freq, file.path(run_dir, "cluster_freq.csv"))

# Cytokine+ % in CD4+ and CD8+ T cells ----------------------------------------------------
# openCyto mindensity cut on one channel; NA when openCyto returns no gate.
mindensity_cut <- possibly(
  \(x, range = MINDENSITY_RANGE) {
    gate <- openCyto::gate_mindensity(
      flowCore::flowFrame(cbind(v = x)),
      channel = "v",
      gate_range = range
    )
    unname(gate@min)
  },
  otherwise = NA_real_
)

# Gate one arm x donor matrix: CD3+, then CD4+CD8- / CD8+CD4-, then each cytokine per parent.
gate_cytokines <- function(m) {
  cut_cd3 <- mindensity_cut(m[, LINEAGE[["CD3"]]])
  t_cells <- m[which(m[, LINEAGE[["CD3"]]] > cut_cd3), , drop = FALSE]
  cut_cd4 <- mindensity_cut(t_cells[, LINEAGE[["CD4"]]])
  cut_cd8 <- mindensity_cut(t_cells[, LINEAGE[["CD8"]]])
  hi_cd4 <- t_cells[, LINEAGE[["CD4"]]] > cut_cd4
  hi_cd8 <- t_cells[, LINEAGE[["CD8"]]] > cut_cd8
  parents <- list(
    `CD4+` = t_cells[which(hi_cd4 & !hi_cd8), , drop = FALSE],
    `CD8+` = t_cells[which(hi_cd8 & !hi_cd4), , drop = FALSE]
  )
  cytokine <- parents |>
    imap(\(cells, parent) {
      tibble(
        parent = parent,
        cytokine = names(CYTOKINES),
        threshold = map_dbl(CYTOKINES, \(ch) {
          mindensity_cut(cells[, ch], CYTOKINE_RANGE)
        }),
        pct_positive = map2_dbl(CYTOKINES, threshold, \(ch, thr) {
          100 * mean(cells[, ch] > thr)
        })
      )
    }) |>
    list_rbind()
  # arcsinh density per parent and cytokine, for the cut-inspection figure
  density_tbl <- parents |>
    imap(\(cells, parent) {
      CYTOKINES |>
        imap(\(ch, cyt) {
          d <- stats::density(cells[, ch], from = -1, to = 6, n = 128)
          tibble(parent = parent, cytokine = cyt, x = d$x, y = d$y)
        }) |>
        list_rbind()
    }) |>
    list_rbind()
  lineage <- tibble(
    marker = names(LINEAGE),
    threshold = c(cut_cd3, cut_cd4, cut_cd8)
  )
  list(cytokine = cytokine, lineage = lineage, density = density_tbl)
}

gated <- arm_mats |>
  imap(\(m, arm) {
    samples |>
      set_names() |>
      map(\(s) gate_cytokines(m[events$sample == s, , drop = FALSE]))
  })

# Pull one table (cytokine or lineage) out of the nested arm x sample list.
bind_gated <- function(part) {
  gated |>
    map(\(by_sample) {
      by_sample |> map(\(g) g[[part]]) |> list_rbind(names_to = "sample")
    }) |>
    list_rbind(names_to = "arm")
}
cytokine_freq <- bind_gated("cytokine") |>
  select(arm, sample, parent, cytokine, pct_positive, threshold)
write_csv(cytokine_freq, file.path(run_dir, "cytokine_freq.csv"))
write_csv(bind_gated("lineage"), file.path(run_dir, "lineage_thresholds.csv"))
write_csv(bind_gated("density"), file.path(run_dir, "cytokine_density.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_downstream.txt")
)
