# experiments/integration/glehr2026/08_af700_pb_negative.R
# Purpose: test whether sw_negative_audit() flags the cause suggested by 07 (fig_af700_pb_deciles):
#          after CD3 relativisation, the negative population of AF700 and PB differs between the
#          Aurora and the Fortessa while the positive peak lines up. Audit the authors' relativised
#          BAD files (9 donors x 2 instruments) and compare neg_median and neg_sd by instrument.
# Run:     Rscript experiments/integration/glehr2026/08_af700_pb_negative.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*relativised*12-panel.fcs
# Outputs: <run_dir>/negative_audit.csv
#          <run_dir>/fig_af700_pb_negative.png, figure_source_data_af700_pb_negative.csv
#
# Settings as in aurora_exvivo/06_negative_real.R (illustrative); the scale asinh(x / 1000) of the
# relativised values matches 07's figure (CD3- at 0, CD3+ at 1e4).

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

CHANNELS <- c("AF700-A", "PB-A", "FITC-A", "APC-A", "PE-A", "PC7-A")
N_EVENTS <- 2000
HYPER_K <- 4
INSTRUMENTS <- c("aurora", "fortessa")

donor_ids <- sort(basename(list.dirs(file.path(data_dir, "BAD", "aurora"), recursive = FALSE)))
grid <- expand.grid(donor = seq_along(donor_ids), instrument = INSTRUMENTS, stringsAsFactors = FALSE)
files <- lapply(seq_len(nrow(grid)), function(j) {
  d <- file.path(data_dir, "BAD", grid$instrument[j], donor_ids[grid$donor[j]])
  f <- list.files(d, pattern = "relativised.*12-panel\\.fcs$", full.names = TRUE)
  e <- exprs(read.FCS(f, truncate_max_range = FALSE, transformation = FALSE, emptyValue = FALSE))
  asinh(e[, CHANNELS] / 1000)
})
names(files) <- paste0("d", grid$donor, "_", grid$instrument)

audit <- sw_negative_audit(files, CHANNELS, N_EVENTS, seed = 1, hyper_k = HYPER_K, pairs = NULL)
write.csv(audit, file.path(run_dir, "negative_audit.csv"), row.names = FALSE)

plot_frame <- audit |>
  filter(check %in% c("neg_median", "neg_sd", "separation")) |>
  mutate(instrument = sub(".*_", "", file), donor = sub("_.*", "", file))
write.csv(plot_frame, file.path(run_dir, "figure_source_data_af700_pb_negative.csv"), row.names = FALSE)
print(as.data.frame(plot_frame |>
  summarise(value = round(median(value, na.rm = TRUE), 3), .by = c(check, channel, instrument)) |>
  tidyr::pivot_wider(names_from = instrument, values_from = value)))

p <- ggplot(plot_frame, aes(instrument, value, group = donor)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 0.8) +
  facet_grid(rows = vars(check), cols = vars(channel), scales = "free_y") +
  labs(x = NULL, y = "Negative population on asinh(x / 1000) of relativised values (lines: donors)")
ggsave(file.path(run_dir, "fig_af700_pb_negative.png"), p, width = 180, height = 110, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
