# experiments/integration/glehr2026/02_reproduce_cd3.R
# Purpose: fit m and M from the CD3 single-stain controls with sw_relativise_split_kmeans() and
#          sw_relativise_fit(), and compare them with the authors' own m and M recovered by
#          01_invert_authors.R (Glehr et al. 2026, BAD study, Aurora, 9 donors).
# Run:     Rscript experiments/integration/glehr2026/02_reproduce_cd3.R <data_dir> <invert_run_dir> <run_dir> <instrument>
# Inputs:  <data_dir>/BAD/aurora/<donor>/*_0N-CD3-<fluor>.fcs   CD3 control per channel
#          <data_dir>/BAD/aurora/<donor>/*relativised*12-panel.fcs   authors' gated events (scatter range)
#          <invert_run_dir>/authors_m_M.csv   authors' m and M per donor and channel
# Outputs: <run_dir>/compare_m_M.csv   donor, channel, ours and authors' m, M, log ratios
#          <run_dir>/fig_cd3_reproduce.png, figure_source_data_cd3_reproduce.csv
#
# Choices approved by the user on 2026-11-10 (Q1 to Q4):
#   Q1 gate:  FSC-A x SSC-A box from the 1-99% quantiles of the authors' gated panel events, per donor
#   Q2 split: k-means (k = 2) on asinh(x / 6000)
#   Q3 bar:   median |log(M_ours / M_authors)| <= 0.05 and worst donor <= 0.15
#   Q4 knobs: seed 1, nstart 10, min_events 200

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[[1]]
invert_dir <- args[[2]]
run_dir <- args[[3]]
instrument <- args[[4]] # "aurora" or "fortessa"
if (dir.exists(run_dir)) stop("run dir exists: ", run_dir, call. = FALSE)
dir.create(run_dir, recursive = TRUE)

GATE_PROBS <- c(0.01, 0.99) # Q1
COFACTOR <- 6000 # Q2
SEED <- 1 # Q4
NSTART <- 10 # Q4
MIN_EVENTS <- 200 # Q4
BAR_MEDIAN <- 0.05 # Q3
BAR_WORST <- 0.15 # Q3

read_exprs <- function(path) {
  exprs(read.FCS(path, truncate_max_range = FALSE, transformation = FALSE, emptyValue = FALSE))
}

# Fit per donor --------------------------------------------------------------------------------

donors <- sort(list.dirs(file.path(data_dir, "BAD", instrument), recursive = FALSE))
ours <- lapply(seq_along(donors), function(i) {
  unit <- paste0("d", i)
  fs <- list.files(donors[i], pattern = "\\.fcs$", full.names = TRUE)
  # Scatter box from the authors' gated events of the same donor.
  gated <- read_exprs(fs[grepl("relativised", fs)])
  box <- apply(gated[, c("FSC-A", "SSC-A")], 2, stats::quantile, probs = GATE_PROBS)
  ctrl_files <- fs[grepl("_ungated_.*_[0-9]{2}-CD3-", fs)]
  fluor <- sub(".*-CD3-(.*)\\.fcs$", "\\1", ctrl_files)
  controls <- lapply(ctrl_files, function(f) {
    e <- read_exprs(f)
    keep <- e[, "FSC-A"] >= box[1, 1] & e[, "FSC-A"] <= box[2, 1] &
      e[, "SSC-A"] >= box[1, 2] & e[, "SSC-A"] <= box[2, 2]
    e[keep, , drop = FALSE]
  })
  names(controls) <- paste0(unit, "_", fluor)
  control_map <- data.frame(unit = unit, channel = paste0(fluor, "-A"), control = names(controls))
  positive <- sw_relativise_split_kmeans(
    controls,
    control_map,
    split_transform = \(x) asinh(x / COFACTOR),
    seed = SEED,
    nstart = NSTART
  )
  sw_relativise_fit(controls, control_map, positive, min_events = MIN_EVENTS) |>
    mutate(n_gated = vapply(controls[control], nrow, integer(1)))
}) |>
  bind_rows()

# Compare with the authors ---------------------------------------------------------------------

authors <- read.csv(file.path(invert_dir, "authors_m_M.csv"))
cmp <- ours |>
  select(donor = unit, channel, n_gated, n_negative, n_positive, m, M, separation) |>
  inner_join(
    select(authors, donor, channel, m_authors = m, M_authors = M),
    by = join_by(donor, channel),
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  mutate(
    log_ratio_M = log(M / M_authors),
    log_ratio_gain = log((M - m) / (M_authors - m_authors)),
    m_diff = m - m_authors
  )
write.csv(cmp, file.path(run_dir, "compare_m_M.csv"), row.names = FALSE)

verdict <- cmp |>
  summarise(
    median_abs_log_M = median(abs(log_ratio_M)),
    worst_abs_log_M = max(abs(log_ratio_M)),
    median_abs_log_gain = median(abs(log_ratio_gain)),
    worst_abs_log_gain = max(abs(log_ratio_gain))
  ) |>
  mutate(pass = median_abs_log_M <= BAR_MEDIAN & worst_abs_log_M <= BAR_WORST)
write.csv(verdict, file.path(run_dir, "verdict.csv"), row.names = FALSE)
print(as.data.frame(verdict))
print(as.data.frame(summarise(cmp, med = median(log_ratio_M), worst = max(abs(log_ratio_M)), .by = channel)))

# Figure ---------------------------------------------------------------------------------------

plot_frame <- cmp |>
  tidyr::pivot_longer(c(log_ratio_M, log_ratio_gain), names_to = "quantity", values_to = "log_ratio")
write.csv(plot_frame, file.path(run_dir, "figure_source_data_cd3_reproduce.csv"), row.names = FALSE)
p <- ggplot(plot_frame, aes(log_ratio, channel)) +
  annotate("rect", xmin = -BAR_MEDIAN, xmax = BAR_MEDIAN, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = c(-BAR_WORST, BAR_WORST), linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0) +
  geom_point(size = 0.9, alpha = 0.7) +
  facet_wrap(vars(quantity)) +
  labs(x = "log(ours / authors), one point per donor; band 0.05, dashes 0.15", y = NULL)
ggsave(file.path(run_dir, "fig_cd3_reproduce.png"), p, width = 160, height = 70, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
