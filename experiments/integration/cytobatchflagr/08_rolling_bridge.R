# experiments/integration/cytobatchflagr/08_rolling_bridge.R
# Purpose: can a rolling reference or a bridge-sample adjustment bring the frozen-reference false
#          flags on clean later batches (9.4% at |z| >= 3, G36) back toward the reference in-sample
#          baseline, while a fault in 80% of files is still found? Three arms score the same files:
#          frozen (reference batch01-09), rolling (the 9 batches before each scored batch) and bridge
#          (each batch shifted and scaled per marker so its bridge file matches the bridge's
#          batch01-09 median and SD, then a frozen reference on the adjusted batch01-09). Each of
#          the 4 controls is the bridge in turn; the other 3 controls are reference and scored files.
# Run:     Rscript experiments/integration/cytobatchflagr/08_rolling_bridge.R <extract_run_dir> <run_dir> <n_seeds>
# Inputs:  <extract_run_dir>/events.rds, files.csv, marker_names.rds   from 01_extract.R
# Outputs: <run_dir>/null_z.csv              bridge, arm, file, batch, check, marker, z (clean files)
#          <run_dir>/fault_z.csv             bridge, seed, marker, file, faulty, arm, median_z (hit marker)
#          <run_dir>/summary.csv             per bridge and arm: false flags at |z| >= 3, AUROC, hits
#          <run_dir>/cutoff_share.csv        share of |z| >= 2-5 per bridge, arm and kind
#          <run_dir>/fig_rolling_bridge_cutoff.png, figure_source_data_rolling_bridge_cutoff.csv
#          <run_dir>/fig_rolling_bridge_null_z.png, figure_source_data_rolling_bridge_null_z.csv
#          <run_dir>/session_info.txt
#
# Decisions: G37 (this test), G36 (cutoff grid, baseline), G32 (batch01-09 reference), G33
# (median_z, sd_z), G34 (fault and seed rule as 06), G9 (25 named markers, asinh(x / 6000)).

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
extract_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) {
  stop("run dir exists: ", run_dir, call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# Settings --------------------------------------------------------------------------------------
SEEDS <- seq_len(as.integer(args[[3]])) # G34: 10 seeds
SHIFT_SD <- 1 # G34
FRACTION_BAD <- 0.8 # G34, of the non-bridge scored files
REF_BATCHES <- 1:9 # G32
SCORED_BATCHES <- 10:18 # G32
ROLL_WINDOW <- 9 # G37: batches before each scored batch
COFACTOR <- 6000 # G9
CUTOFFS <- c(2, 3, 4, 5) # G36
RULE_FLAGS <- 0.056 # G37: false flags at |z| >= 3 at most in-sample 4.1% + 1.5 points
RULE_AUROC <- 0.9 # G37 / D1
G <- 1 # seed rule as 06

# Data ------------------------------------------------------------------------------------------
ev <- readRDS(file.path(extract_dir, "events.rds"))
info <- read.csv(file.path(extract_dir, "files.csv"))
marker_names <- readRDS(file.path(extract_dir, "marker_names.rds"))
stopifnot("files match events" = identical(info$file, names(ev)))
ids <- paste0("f", seq_along(ev))
names(ev) <- ids
shared <- Reduce(intersect, lapply(ev, colnames))
markers <- intersect(names(marker_names)[!is.na(marker_names)], shared)
base <- lapply(ev, \(m) asinh(m[, markers] / COFACTOR))
file_info <- tibble(
  file = ids,
  batch = info$batch,
  batch_no = as.integer(sub("batch", "", info$batch, fixed = TRUE)),
  control = info$control
)
stopifnot(
  "one file per batch and control" = !anyDuplicated(file_info[c(
    "batch_no",
    "control"
  )])
)
controls <- sort(unique(file_info$control))
marker_sd <- apply(do.call(rbind, base), 2, sd)

audit_z <- function(files, reference) {
  sa <- sw_fcs_signal_audit(
    files,
    markers,
    NULL,
    NULL,
    time_channel = NULL,
    time_bin = NULL,
    area_height = NULL,
    area_height_quantile = NULL,
    reference_group = NULL,
    reference = reference,
    summary = NULL
  )
  sa <- sa[
    sa$check %in% c("median_z", "sd_z"),
    c("file", "check", "channel", "value")
  ]
  names(sa) <- c("file", "check", "marker", "z")
  sa
}

# Bridge adjustment: per batch and marker, map the bridge file onto its batch01-09 median and SD.
bridge_adjust <- function(data, bridge) {
  br <- file_info[file_info$control == bridge, ]
  br_med <- t(vapply(
    data[br$file],
    \(m) apply(m, 2, median),
    numeric(length(markers))
  ))
  br_sd <- t(vapply(
    data[br$file],
    \(m) apply(m, 2, sd),
    numeric(length(markers))
  ))
  in_ref <- br$batch_no %in% REF_BATCHES
  target_med <- apply(br_med[in_ref, , drop = FALSE], 2, median)
  target_sd <- apply(br_sd[in_ref, , drop = FALSE], 2, median)
  out <- data
  for (f in names(data)) {
    k <- match(file_info$batch_no[file_info$file == f], br$batch_no)
    gain <- target_sd / br_sd[k, ]
    out[[f]] <- sweep(sweep(data[[f]], 2, br_med[k, ], `-`), 2, gain, `*`) |>
      sweep(2, target_med, `+`)
  }
  out
}

# Score the scored files of one condition under the three arms.
score_arms <- function(data, bridge) {
  keep <- file_info$control != bridge
  ref_ids <- file_info$file[keep & file_info$batch_no %in% REF_BATCHES]
  scored_ids <- file_info$file[keep & file_info$batch_no %in% SCORED_BATCHES]
  frozen_ref <- sw_signal_reference_fit(data[ref_ids], markers)
  adjusted <- bridge_adjust(data, bridge)
  bridge_ref <- sw_signal_reference_fit(adjusted[ref_ids], markers)
  rolling <- lapply(SCORED_BATCHES, \(b) {
    win <- seq(b - ROLL_WINDOW, b - 1)
    win_ids <- file_info$file[keep & file_info$batch_no %in% win]
    b_ids <- file_info$file[keep & file_info$batch_no == b]
    audit_z(data[b_ids], sw_signal_reference_fit(data[win_ids], markers))
  }) |>
    bind_rows()
  bind_rows(
    mutate(audit_z(data[scored_ids], frozen_ref), arm = "frozen"),
    mutate(rolling, arm = "rolling"),
    mutate(audit_z(adjusted[scored_ids], bridge_ref), arm = "bridge"),
    mutate(audit_z(data[ref_ids], frozen_ref), arm = "reference in-sample")
  )
}

# Null and fault, per bridge --------------------------------------------------------------------
auroc <- function(score, pos) {
  r <- rank(score)
  n1 <- sum(pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * sum(!pos))
}
null_z <- list()
fault_z <- list()
for (bridge in controls) {
  null_z[[bridge]] <- mutate(score_arms(base, bridge), bridge = bridge)
  scored_ids <- file_info$file[
    file_info$control != bridge & file_info$batch_no %in% SCORED_BATCHES
  ]
  fault_z[[bridge]] <- lapply(SEEDS, \(seed) {
    set.seed(seed * 1000 + G)
    mk <- sample(markers, 1)
    bad <- scored_ids %in%
      sample(scored_ids, round(FRACTION_BAD * length(scored_ids)))
    shifted <- base
    shifted[scored_ids[bad]] <- lapply(shifted[scored_ids[bad]], \(x) {
      x[, mk] <- x[, mk] + SHIFT_SD * marker_sd[[mk]]
      x
    })
    score_arms(shifted, bridge) |>
      filter(arm != "reference in-sample", check == "median_z", marker == mk) |>
      transmute(
        bridge = bridge,
        seed = seed,
        marker = mk,
        file,
        faulty = file %in% scored_ids[bad],
        arm,
        median_z = z
      )
  }) |>
    bind_rows()
}
null_z <- bind_rows(null_z) |>
  left_join(
    file_info[c("file", "batch")],
    by = join_by(file),
    relationship = "many-to-one",
    unmatched = "error"
  )
fault_z <- bind_rows(fault_z)
write.csv(null_z, file.path(run_dir, "null_z.csv"), row.names = FALSE)
write.csv(fault_z, file.path(run_dir, "fault_z.csv"), row.names = FALSE)

# Summary ---------------------------------------------------------------------------------------
flags <- null_z |>
  filter(!is.na(z)) |>
  summarise(false_flags_3 = mean(abs(z) >= 3), .by = c(bridge, arm))
fault_sum <- fault_z |>
  filter(!is.na(median_z)) |>
  summarise(
    auroc = auroc(abs(median_z), faulty),
    hits_3 = mean(abs(median_z[faulty]) >= 3),
    .by = c(bridge, arm, seed)
  ) |>
  summarise(
    auroc_mean = mean(auroc),
    auroc_ci95_half = 1.96 * sd(auroc) / sqrt(n()),
    hits_3 = mean(hits_3),
    .by = c(bridge, arm)
  )
summary_tbl <- left_join(
  flags,
  fault_sum,
  by = join_by(bridge, arm),
  relationship = "one-to-one"
) |>
  mutate(improves = false_flags_3 <= RULE_FLAGS & auroc_mean >= RULE_AUROC)
write.csv(summary_tbl, file.path(run_dir, "summary.csv"), row.names = FALSE)
print(as.data.frame(mutate(
  summary_tbl,
  across(where(is.numeric), \(x) round(x, 3))
)))

# Share at or above each cutoff -------------------------------------------------------------------
false_share <- null_z |>
  filter(!is.na(z)) |>
  reframe(
    cutoff = CUTOFFS,
    share = vapply(CUTOFFS, \(k) mean(abs(z) >= k), numeric(1)),
    .by = c(bridge, arm)
  ) |>
  mutate(kind = "clean files, all markers (false flags)")
hit_share <- fault_z |>
  filter(faulty, !is.na(median_z)) |>
  reframe(
    cutoff = CUTOFFS,
    share = vapply(CUTOFFS, \(k) mean(abs(median_z) >= k), numeric(1)),
    .by = c(bridge, arm)
  ) |>
  mutate(kind = "faulty files, hit marker (hits)")
cutoff_share <- bind_rows(false_share, hit_share)
write.csv(
  cutoff_share,
  file.path(run_dir, "cutoff_share.csv"),
  row.names = FALSE
)

# Figures ---------------------------------------------------------------------------------------
arm_levels <- c("reference in-sample", "frozen", "rolling", "bridge")
arm_colours <- c(
  `reference in-sample` = "grey60",
  frozen = "#D55E00",
  rolling = "#0072B2",
  bridge = "#009E73"
)
cut_plot <- mutate(cutoff_share, arm = factor(arm, levels = arm_levels))
write.csv(
  cut_plot,
  file.path(run_dir, "figure_source_data_rolling_bridge_cutoff.csv"),
  row.names = FALSE
)
p_cut <- ggplot(cut_plot, aes(cutoff, share, colour = arm)) +
  geom_line(aes(group = interaction(arm, bridge)), alpha = 0.35) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.9) +
  stat_summary(fun = mean, geom = "point", size = 1.2) +
  facet_wrap(vars(kind), nrow = 1) +
  scale_colour_manual(values = arm_colours) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "|z| cutoff",
    y = "Share at or above the cutoff",
    colour = NULL,
    caption = sprintf(
      "Thin lines: each of the 4 bridge controls; thick: mean. Hits: %d seeds, +1 SD in 80%% of scored files",
      length(SEEDS)
    )
  )
ggsave(
  file.path(run_dir, "fig_rolling_bridge_cutoff.png"),
  p_cut,
  width = 170,
  height = 75,
  units = "mm",
  dpi = 200
)

null_plot <- null_z |>
  filter(!is.na(z), arm != "reference in-sample") |>
  mutate(abs_z = abs(z), arm = factor(arm, levels = arm_levels[-1]))
write.csv(
  null_plot,
  file.path(run_dir, "figure_source_data_rolling_bridge_null_z.csv"),
  row.names = FALSE
)
p_null <- ggplot(null_plot, aes(batch, abs_z, colour = arm)) +
  geom_hline(yintercept = 3, linetype = 3, colour = "grey60") +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.8),
    linewidth = 0.3
  ) +
  facet_wrap(vars(check), ncol = 1) +
  scale_colour_manual(values = arm_colours[-1]) +
  coord_cartesian(ylim = c(0, 5)) +
  labs(
    x = NULL,
    y = "|z| per file and marker (clean files)",
    colour = NULL,
    caption = "All 4 bridge rotations pooled; boxes without outliers, y cut at 5; dotted: |z| = 3"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(run_dir, "fig_rolling_bridge_null_z.png"),
  p_null,
  width = 170,
  height = 110,
  units = "mm",
  dpi = 200
)

writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
