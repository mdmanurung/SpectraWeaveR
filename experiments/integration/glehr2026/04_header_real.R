# experiments/integration/glehr2026/04_header_real.R
# Purpose: real test of sw_fcs_header_audit(). The 9 BAD donors' panel files from the Aurora and the
#          Fortessa are audited against donor 1's Aurora file. Files from the other instrument
#          should differ on the instrument keywords and channel ranges; same-instrument files
#          should match on everything but per-file values (event count, dates).
# Run:     Rscript experiments/integration/glehr2026/04_header_real.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*_ungated_*12-panel.fcs
# Outputs: <run_dir>/header_audit.csv     file label (d<i>_<instrument>), kind, key, channel, status
#          <run_dir>/header_summary.csv   status counts per file and key
#          <run_dir>/fig_header_real.png, figure_source_data_header_real.csv

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) stop("run dir exists: ", run_dir, call. = FALSE)
dir.create(run_dir, recursive = TRUE)

KEYWORDS <- c("$CYT", "$CYTSN", "$PnN", "$PnS", "$PnR", "$PnE", "$PnV", "$TOT", "$DATE")

donor_ids <- sort(basename(list.dirs(file.path(data_dir, "BAD", "aurora"), recursive = FALSE)))
files <- expand.grid(donor = seq_along(donor_ids), instrument = c("aurora", "fortessa"), stringsAsFactors = FALSE)
files$path <- vapply(seq_len(nrow(files)), function(j) {
  d <- file.path(data_dir, "BAD", files$instrument[j], donor_ids[files$donor[j]])
  list.files(d, pattern = "_ungated_.*12-panel\\.fcs$", full.names = TRUE)
}, character(1))
files$label <- paste0("d", files$donor, "_", files$instrument)

audit <- sw_fcs_header_audit(files$path, KEYWORDS, reference = files$path[1])
audit$file <- files$label[match(audit$file, basename(files$path))]
write.csv(audit, file.path(run_dir, "header_audit.csv"), row.names = FALSE)

# Rows that do not match, per file and key (channel rows collapsed to their key).
plot_frame <- audit |>
  filter(file != "d1_aurora") |>
  summarise(n_not_match = sum(status != "match"), .by = c(file, key)) |>
  mutate(
    instrument = sub(".*_", "", file),
    file = factor(file, levels = files$label)
  )
write.csv(plot_frame, file.path(run_dir, "figure_source_data_header_real.csv"), row.names = FALSE)
write.csv(sw_fcs_header_summary(audit), file.path(run_dir, "header_summary.csv"), row.names = FALSE)
print(as.data.frame(summarise(plot_frame, files_with_difference = sum(n_not_match > 0), .by = c(instrument, key))))

p <- ggplot(plot_frame, aes(key, file, fill = n_not_match)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "grey95", high = "#B2182B") +
  labs(x = "Header key", y = NULL, fill = "rows not\nmatching") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(run_dir, "fig_header_real.png"), p, width = 120, height = 100, units = "mm", dpi = 200)
