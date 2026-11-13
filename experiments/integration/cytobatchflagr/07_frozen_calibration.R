# experiments/integration/cytobatchflagr/07_frozen_calibration.R
# Purpose: is the frozen-reference z calibrated on later batches that carry real batch offsets?
#          (1) Null: no fault. |median_z| and |sd_z| per file and marker for the clean scored files
#          (batch10-18) under the frozen reference and under the cohort z, and for the reference
#          files scored in-sample (baseline); trend of |z| with acquisition date. (2) The G34 fault
#          (+1 SD on one marker in 80% of scored files) rerun with per-file z kept, so hits and false
#          flags can be read at |z| >= 2, 3, 4, 5. No cutoff is chosen here.
# Run:     Rscript experiments/integration/cytobatchflagr/07_frozen_calibration.R <extract_run_dir> <run_dir> <n_seeds>
# Inputs:  <extract_run_dir>/events.rds, files.csv, marker_names.rds   from 01_extract.R
# Outputs: <run_dir>/null_z.csv                  file, batch, date, set, arm, check, marker, z
#          <run_dir>/null_trend.csv              Spearman of |z| with date, per arm and check
#          <run_dir>/fault_z.csv                 seed, marker, file, faulty, arm, median_z of hit marker
#          <run_dir>/fault_auroc.csv             seed, arm, AUROC (should match 06)
#          <run_dir>/cutoff_share.csv            share of |z| >= cutoff: false flags and hits
#          <run_dir>/fig_frozen_null_z.png, figure_source_data_frozen_null_z.csv
#          <run_dir>/fig_frozen_cutoff.png, figure_source_data_frozen_cutoff.csv
#          <run_dir>/session_info.txt
#
# Decisions: G36 (this test), G32 (reference batch01-09, scored batch10-18), G33 (median_z and
# sd_z), G34 (fault, seeds, seed rule as 06), G9 (25 named markers, asinh(x / 6000)).

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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
FRACTION_BAD <- 0.8 # G34
REF_BATCHES <- 1:9 # G32
SCORED_BATCHES <- 10:18 # G32
COFACTOR <- 6000 # G9
CUTOFFS <- c(2, 3, 4, 5) # G36: reported, none chosen
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
  date = as.Date(info$date, format = "%d-%b-%Y"),
  batch_no = as.integer(sub("batch", "", info$batch, fixed = TRUE))
)
ref_ids <- ids[file_info$batch_no %in% REF_BATCHES]
scored_ids <- ids[file_info$batch_no %in% SCORED_BATCHES]
marker_sd <- apply(do.call(rbind, base), 2, sd)
ref <- sw_signal_reference_fit(base[ref_ids], markers)

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
  sa[
    sa$check %in% c("median_z", "sd_z"),
    c("file", "check", "channel", "value")
  ]
}

# (1) Null: clean files ---------------------------------------------------------------------------
# Scored files under both arms; reference files in-sample under the frozen arm as the baseline.
null_z <- bind_rows(
  mutate(audit_z(base[scored_ids], ref), set = "scored", arm = "frozen"),
  mutate(audit_z(base[scored_ids], NULL), set = "scored", arm = "cohort"),
  mutate(
    audit_z(base[ref_ids], ref),
    set = "reference (in-sample)",
    arm = "frozen"
  )
) |>
  rename(marker = channel, z = value) |>
  left_join(
    file_info,
    by = join_by(file),
    relationship = "many-to-one",
    unmatched = "error"
  )
write.csv(null_z, file.path(run_dir, "null_z.csv"), row.names = FALSE)

# Trend of |z| with date over the scored files, per arm and check.
null_trend <- null_z |>
  filter(set == "scored", !is.na(z)) |>
  summarise(
    spearman_abs_z_date = cor(abs(z), as.numeric(date), method = "spearman"),
    median_abs_z = median(abs(z)),
    .by = c(arm, check)
  )
write.csv(null_trend, file.path(run_dir, "null_trend.csv"), row.names = FALSE)
print(as.data.frame(mutate(
  null_trend,
  across(where(is.numeric), \(x) round(x, 3))
)))

# (2) Fault rerun with per-file z -------------------------------------------------------------
fault_z <- lapply(SEEDS, \(seed) {
  set.seed(seed * 1000 + G)
  mk <- sample(markers, 1)
  bad <- scored_ids %in%
    sample(scored_ids, round(FRACTION_BAD * length(scored_ids)))
  shifted <- base[scored_ids]
  shifted[bad] <- lapply(shifted[bad], \(x) {
    x[, mk] <- x[, mk] + SHIFT_SD * marker_sd[[mk]]
    x
  })
  lapply(c("cohort", "frozen"), \(which_arm) {
    a <- audit_z(shifted, if (which_arm == "frozen") ref else NULL)
    mz <- a[a$check == "median_z" & a$channel == mk, ]
    tibble(
      seed = seed,
      marker = mk,
      file = scored_ids,
      faulty = bad,
      arm = which_arm,
      median_z = mz$value[match(scored_ids, mz$file)]
    )
  }) |>
    bind_rows()
}) |>
  bind_rows()
write.csv(fault_z, file.path(run_dir, "fault_z.csv"), row.names = FALSE)

auroc <- function(score, pos) {
  r <- rank(score)
  n1 <- sum(pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * sum(!pos))
}
fault_auroc <- fault_z |>
  summarise(auroc = auroc(abs(median_z), faulty), .by = c(seed, arm))
write.csv(fault_auroc, file.path(run_dir, "fault_auroc.csv"), row.names = FALSE)
print(as.data.frame(summarise(
  fault_auroc,
  auroc = round(mean(auroc), 3),
  .by = arm
)))

# Share of |z| at or above each cutoff ----------------------------------------------------------
# False flags: clean scored files, every marker and both checks (null). Hits: faulty files,
# hit marker, median_z. The reference in-sample baseline is shown for the frozen arm.
false_flags <- null_z |>
  filter(!is.na(z)) |>
  reframe(
    cutoff = CUTOFFS,
    share = vapply(CUTOFFS, \(k) mean(abs(z) >= k), numeric(1)),
    .by = c(set, arm)
  ) |>
  mutate(kind = paste0("clean files, all markers (", set, ")"))
hits <- fault_z |>
  filter(faulty, !is.na(median_z)) |>
  reframe(
    cutoff = CUTOFFS,
    share = vapply(CUTOFFS, \(k) mean(abs(median_z) >= k), numeric(1)),
    .by = arm
  ) |>
  mutate(set = "scored", kind = "faulty files, hit marker (+1 SD)")
cutoff_share <- bind_rows(false_flags, hits) |>
  select(kind, set, arm, cutoff, share)
write.csv(
  cutoff_share,
  file.path(run_dir, "cutoff_share.csv"),
  row.names = FALSE
)
print(as.data.frame(mutate(cutoff_share, share = round(share, 3))))

# Figures ---------------------------------------------------------------------------------------
null_plot <- null_z |>
  filter(!is.na(z)) |>
  mutate(
    abs_z = abs(z),
    panel = paste(arm, "z,", set)
  )
write.csv(
  null_plot,
  file.path(run_dir, "figure_source_data_frozen_null_z.csv"),
  row.names = FALSE
)
p_null <- ggplot(null_plot, aes(batch, abs_z)) +
  geom_hline(yintercept = CUTOFFS, linetype = 3, colour = "grey70") +
  geom_point(
    position = position_jitter(width = 0.25, height = 0, seed = 1),
    size = 0.3,
    alpha = 0.3
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 3,
    size = 2,
    colour = "#D55E00"
  ) +
  facet_grid(
    rows = vars(check),
    cols = vars(panel),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 2, 3, 5, 10, 20, 50)) +
  labs(
    x = NULL,
    y = "|z| per file and marker (log1p scale)",
    caption = "Clean CBF files, no fault; 25 markers x 4 files per batch; cross: batch median; dotted: |z| = 2, 3, 4, 5"
  ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(
  file.path(run_dir, "fig_frozen_null_z.png"),
  p_null,
  width = 180,
  height = 110,
  units = "mm",
  dpi = 200
)

write.csv(
  cutoff_share,
  file.path(run_dir, "figure_source_data_frozen_cutoff.csv"),
  row.names = FALSE
)
p_cut <- ggplot(cutoff_share, aes(cutoff, share, colour = arm)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(vars(kind), nrow = 1) +
  scale_colour_manual(values = c(cohort = "grey40", frozen = "#D55E00")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "|z| cutoff",
    y = "Share at or above the cutoff",
    colour = NULL,
    caption = sprintf(
      "Hits: %d seeds, 80%% of 36 scored files shifted by +1 SD",
      length(SEEDS)
    )
  )
ggsave(
  file.path(run_dir, "fig_frozen_cutoff.png"),
  p_cut,
  width = 180,
  height = 70,
  units = "mm",
  dpi = 200
)

writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
