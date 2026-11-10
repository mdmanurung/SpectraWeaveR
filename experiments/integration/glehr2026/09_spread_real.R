# experiments/integration/glehr2026/09_spread_real.R
# Purpose: test whether spillover spread explains the wider AF700 and PB negatives on the Fortessa
#          (08: 5.6x and 2.8x wider). For every directed pair among the panel channels,
#          sw_negative_audit() reports spread_log_ratio: the target-negative SD where the source is
#          high over where it is low. A source that spreads into AF700 or PB on the Fortessa only
#          shows a higher ratio on the Fortessa than on the Aurora for the same donor.
# Run:     Rscript experiments/integration/glehr2026/09_spread_real.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*relativised*12-panel.fcs
# Outputs: <run_dir>/spread_audit.csv
#          <run_dir>/fig_spread_af700_pb.png, figure_source_data_spread_af700_pb.csv
#
# Scale, channels, events and hyper_k as in 08 (illustrative). The two quantiles are reported as
# variants; spread_quantile is a required argument the user approved as the one knob.

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

# All ten fluorescent channels of the 12-panel files, so every possible source is tested.
CHANNELS <- c(
  "AF700-A", "PB-A", "FITC-A", "APC-A", "PE-A", "PC7-A", "ECD-A", "PC5.5-A", "AA750-A", "KrO-A"
)
TARGETS <- c("AF700-A", "PB-A")
N_EVENTS <- 2000
HYPER_K <- 4
QUANTILES <- c(0.9, 0.95)
INSTRUMENTS <- c("aurora", "fortessa")

# Read files ------------------------------------------------------------------
donor_ids <- sort(basename(list.dirs(file.path(data_dir, "BAD", "aurora"), recursive = FALSE)))
grid <- expand.grid(donor = seq_along(donor_ids), instrument = INSTRUMENTS, stringsAsFactors = FALSE)
files <- lapply(seq_len(nrow(grid)), function(j) {
  d <- file.path(data_dir, "BAD", grid$instrument[j], donor_ids[grid$donor[j]])
  f <- list.files(d, pattern = "relativised.*12-panel\\.fcs$", full.names = TRUE)
  e <- exprs(read.FCS(f, truncate_max_range = FALSE, transformation = FALSE, emptyValue = FALSE))
  asinh(e[, CHANNELS] / 1000)
})
names(files) <- paste0("d", grid$donor, "_", grid$instrument)

# Spread per directed pair, two quantile variants -----------------------------
audit <- lapply(QUANTILES, function(q) {
  sw_negative_audit(files, CHANNELS, N_EVENTS, seed = 1, hyper_k = HYPER_K, pairs = "all", spread_quantile = q) |>
    filter(check == "spread_log_ratio") |>
    mutate(spread_quantile = q)
}) |>
  bind_rows()
write.csv(audit, file.path(run_dir, "spread_audit.csv"), row.names = FALSE)

# Plot frame: pairs into AF700 and PB, per donor and instrument ----------------
plot_frame <- audit |>
  mutate(
    instrument = sub(".*_", "", file),
    donor = sub("_.*", "", file),
    source = sub(" > .*", "", channel),
    target = sub(".* > ", "", channel)
  ) |>
  filter(target %in% TARGETS)
write.csv(plot_frame, file.path(run_dir, "figure_source_data_spread_af700_pb.csv"), row.names = FALSE)
print(as.data.frame(plot_frame |>
  summarise(ratio = round(exp(median(value, na.rm = TRUE)), 2), .by = c(spread_quantile, target, source, instrument)) |>
  tidyr::pivot_wider(names_from = instrument, values_from = ratio) |>
  arrange(spread_quantile, target, desc(fortessa))))

p <- ggplot(plot_frame, aes(instrument, exp(value), group = donor)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_line(alpha = 0.4) +
  geom_point(size = 0.8) +
  scale_y_log10() +
  facet_grid(rows = vars(target, spread_quantile), cols = vars(source), labeller = label_both) +
  labs(x = NULL, y = "Target-negative SD, source high / source low (lines: donors)")
ggsave(file.path(run_dir, "fig_spread_af700_pb.png"), p, width = 260, height = 150, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
