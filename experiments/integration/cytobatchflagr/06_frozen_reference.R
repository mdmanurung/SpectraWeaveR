# experiments/integration/cytobatchflagr/06_frozen_reference.R
# Purpose: does a frozen reference (sw_signal_reference_fit, passed as `reference` to
#          sw_fcs_signal_audit) find a fault that hits most scored files, where the cohort z cannot?
#          Reference = the 36 files of batches B01-B09, clean. Scored = the 36 files of B10-B18.
#          One hit marker is shifted by +1 SD in 80% of the scored files. Two arms score the same
#          files: "cohort" (reference = NULL, z over the scored files) and "frozen" (z against the
#          reference). Score per file = |median_z| of the hit marker; AUROC ranks faulty above
#          clean files over the scored set. A no-fault null (same files unshifted, same 80% labels)
#          is scored for both arms.
#          Expected: with 80% of files faulty, the fault becomes the cohort norm, so the cohort
#          arm gives AUROC near or below 0.5. The frozen arm does not move its centre.
# Run:     Rscript experiments/integration/cytobatchflagr/06_frozen_reference.R <extract_run_dir> <run_dir> <n_seeds>
# Inputs:  <extract_run_dir>/events.rds, files.csv, marker_names.rds   from 01_extract.R
#          (batch labels in files.csv are "batch01" to "batch18"; four control files per batch)
# Outputs: <run_dir>/frozen_reference_long.csv                one row per seed, condition, arm
#          <run_dir>/frozen_reference_summary.csv             mean AUROC and 95% CI per condition, arm
#          <run_dir>/figure_source_data_frozen_reference.csv  the plotted summary
#          <run_dir>/fig_frozen_reference.png
#          <run_dir>/signal_reference.rds                     the frozen reference (36 files, 25 markers)
#          <run_dir>/session_info.txt
#
# Decisions: G32 (CBF is the real multi-batch dataset; reference B01-B09, scored B10-B18),
# G33 (per-marker median and MAD of file-level medians and SDs; median_z and sd_z only),
# G34 (1 SD shift on one marker in 80% of scored files; 10 seeds), G9 (25 shared named markers (G9 says 30; the code keeps the 25 named ones),
# asinh(x / 6000), as in 04_reference_group.R). Marker SD for the shift is taken over all files,
# as in 04. Seed rule set.seed(seed * 1000 + g) as in 04, with g = 1 (one fault row).
# The bar from the user (D1): mean AUROC >= 0.90.

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
SEEDS <- seq_len(as.integer(args[[3]])) # G34: 10 seeds for the full run
SHIFT_SD <- 1 # G34: shift in marker SD
FRACTION_BAD <- 0.8 # G34: share of scored files shifted
REF_BATCHES <- 1:9 # G32: reference batches B01-B09
SCORED_BATCHES <- 10:18 # G32: scored batches B10-B18
COFACTOR <- 6000 # arcsinh scale, G9
BAR <- 0.9 # user's bar D1
G <- 1 # fault row in the seed rule

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
batch_no <- as.integer(sub("batch", "", info$batch, fixed = TRUE))
ref_ids <- ids[batch_no %in% REF_BATCHES]
scored_ids <- ids[batch_no %in% SCORED_BATCHES]
marker_sd <- apply(do.call(rbind, base), 2, sd)

# Frozen reference from the clean reference files (G33), with a save and load round trip
ref <- sw_signal_reference_fit(base[ref_ids], markers)
sw_signal_reference_save(ref, file.path(run_dir, "signal_reference.rds"))
ref <- sw_signal_reference_load(file.path(run_dir, "signal_reference.rds"))

auroc <- function(score, pos) {
  r <- rank(score)
  n1 <- sum(pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * sum(!pos))
}

# Seeds -----------------------------------------------------------------------------------------

long <- lapply(SEEDS, \(seed) {
  set.seed(seed * 1000 + G)
  mk <- sample(markers, 1) # hit marker
  bad <- scored_ids %in%
    sample(scored_ids, round(FRACTION_BAD * length(scored_ids)))
  names(bad) <- scored_ids
  clean <- base[scored_ids]
  shifted <- clean
  shifted[bad] <- lapply(shifted[bad], \(x) {
    x[, mk] <- x[, mk] + SHIFT_SD * marker_sd[[mk]]
    x
  })
  cond <- list(fault = shifted, null = clean)
  lapply(names(cond), \(condition) {
    lapply(c("cohort", "frozen"), \(arm) {
      sa <- sw_fcs_signal_audit(
        cond[[condition]],
        markers,
        NULL,
        NULL,
        time_channel = NULL,
        time_bin = NULL,
        area_height = NULL,
        area_height_quantile = NULL,
        reference_group = NULL,
        reference = if (arm == "frozen") ref else NULL
      )
      mz <- sa[sa$check == "median_z" & sa$channel == mk, ]
      tibble(
        seed = seed,
        marker = mk,
        condition = condition,
        arm = arm,
        auroc = auroc(abs(mz$value[match(scored_ids, mz$file)]), bad),
        median_abs_z_faulty = median(abs(mz$value[match(
          scored_ids[bad],
          mz$file
        )]))
      )
    }) |>
      bind_rows()
  }) |>
    bind_rows()
}) |>
  bind_rows()
write.csv(
  long,
  file.path(run_dir, "frozen_reference_long.csv"),
  row.names = FALSE
)

# Summary ---------------------------------------------------------------------------------------

summary_tbl <- long |>
  summarise(
    n_seeds = n(),
    ci95_half = 1.96 * sd(auroc) / sqrt(n()),
    auroc_mean = mean(auroc),
    .by = c(condition, arm)
  ) |>
  mutate(meets_bar = auroc_mean >= BAR)
write.csv(
  summary_tbl,
  file.path(run_dir, "frozen_reference_summary.csv"),
  row.names = FALSE
)
write.csv(
  summary_tbl,
  file.path(run_dir, "figure_source_data_frozen_reference.csv"),
  row.names = FALSE
)
print(as.data.frame(mutate(
  summary_tbl,
  across(c(auroc_mean, ci95_half), \(x) round(x, 3))
)))

# Figure ----------------------------------------------------------------------------------------

arm_levels <- c("cohort", "frozen")
long_plot <- mutate(long, arm = factor(arm, levels = arm_levels))
summary_plot <- mutate(summary_tbl, arm = factor(arm, levels = arm_levels))
p <- ggplot(summary_plot, aes(auroc_mean, arm, colour = arm)) +
  geom_vline(xintercept = c(0.5, BAR), linetype = 2, colour = "grey60") +
  geom_point(
    data = long_plot,
    aes(auroc, arm),
    size = 0.6,
    alpha = 0.6,
    position = position_jitter(height = 0.15, width = 0, seed = 1)
  ) +
  geom_pointrange(
    aes(xmin = auroc_mean - ci95_half, xmax = pmin(auroc_mean + ci95_half, 1)),
    size = 0.3
  ) +
  facet_wrap(vars(condition), ncol = 1) +
  scale_colour_manual(values = c(cohort = "grey40", frozen = "#D55E00")) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    x = sprintf(
      "AUROC, faulty vs clean files (%d scored files; mean and 95%% CI over %d seeds; dots = seeds; dashed: 0.5 and bar %.1f)",
      length(scored_ids),
      length(SEEDS),
      BAR
    ),
    y = NULL,
    colour = "arm"
  ) +
  theme(legend.position = "none")
ggsave(
  file.path(run_dir, "fig_frozen_reference.png"),
  p,
  width = 120,
  height = 70,
  units = "mm",
  dpi = 200
)
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
