# experiments/integration/glehr2026/07_af700_pb_deciles.R
# Purpose: find why AF700 (CD8) and PB (CD57) keep an instrument effect after CD3 relativisation
#          (03 and 05: eta2 0.5-0.9 in both studies). Compare each donor's deciles of the authors'
#          relativised values between the Aurora and the Fortessa. A gain mismatch shifts every
#          decile by the same factor; a median that sits in the valley of a bimodal marker moves
#          only the middle deciles.
# Run:     Rscript experiments/integration/glehr2026/07_af700_pb_deciles.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*relativised*12-panel.fcs
# Outputs: <run_dir>/deciles.csv   donor, instrument, channel, prob, value
#          <run_dir>/fig_af700_pb_deciles.png, figure_source_data_af700_pb_deciles.csv

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

CHANNELS <- c("AF700-A", "PB-A", "FITC-A", "APC-A") # two residual channels, two clean ones
PROBS <- seq(0.1, 0.9, by = 0.1)
INSTRUMENTS <- c("aurora", "fortessa")

donor_ids <- sort(basename(list.dirs(file.path(data_dir, "BAD", "aurora"), recursive = FALSE)))
deciles <- lapply(seq_along(donor_ids), function(i) {
  lapply(INSTRUMENTS, function(inst) {
    d <- file.path(data_dir, "BAD", inst, donor_ids[i])
    f <- list.files(d, pattern = "relativised.*12-panel\\.fcs$", full.names = TRUE)
    e <- exprs(read.FCS(f, truncate_max_range = FALSE, transformation = FALSE, emptyValue = FALSE))
    lapply(CHANNELS, function(ch) {
      tibble(
        donor = paste0("d", i),
        instrument = inst,
        channel = ch,
        prob = PROBS,
        value = as.numeric(stats::quantile(e[, ch], PROBS))
      )
    }) |>
      bind_rows()
  }) |>
    bind_rows()
}) |>
  bind_rows()
write.csv(deciles, file.path(run_dir, "deciles.csv"), row.names = FALSE)

# Fortessa minus Aurora per donor, channel and decile, on asinh(x / 1000) of the relativised
# values (CD3- maps to 0 and CD3+ to 1e4, so 1000 is a tenth of the CD3+ median).
paired <- deciles |>
  tidyr::pivot_wider(names_from = instrument, values_from = value) |>
  mutate(diff = asinh(fortessa / 1000) - asinh(aurora / 1000))
write.csv(paired, file.path(run_dir, "figure_source_data_af700_pb_deciles.csv"), row.names = FALSE)
print(as.data.frame(summarise(paired, median_diff = median(diff), .by = c(channel, prob)) |>
  tidyr::pivot_wider(names_from = prob, values_from = median_diff) |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))))

p <- ggplot(paired, aes(prob, diff, group = donor)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(alpha = 0.5) +
  geom_point(size = 0.7) +
  facet_wrap(vars(channel), nrow = 1) +
  labs(x = "Decile", y = "Fortessa minus Aurora, asinh(x / 1000) of relativised value\n(one line per donor)")
ggsave(file.path(run_dir, "fig_af700_pb_deciles.png"), p, width = 180, height = 70, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
