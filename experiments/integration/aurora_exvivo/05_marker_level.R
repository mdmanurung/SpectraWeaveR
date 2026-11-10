# experiments/integration/aurora_exvivo/05_marker_level.R
# Purpose: split the real-data detection limit into two causes. (a) Natural file-to-file spread
#          of each marker's median and SD in the cohort. (b) AUROC of the injected marker's own
#          z-score (an oracle that knows which marker was hit) against the max over 40 markers.
# Run:     Rscript experiments/integration/aurora_exvivo/05_marker_level.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/marker_level_long.csv, marker_level_matrix.csv
#          <run_dir>/natural_spread.csv      per marker: SD of file medians and of file SDs
#          <run_dir>/fig_marker_level.png, figure_source_data_marker_level.csv
#
# Every value below is illustrative; none is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
BATCH <- "Batch1"
SEEDS <- 1:10
FRACTION_BAD <- 0.25
SHIFT_SD <- c(0.3, 0.5, 1)
STRETCH <- c(1.1, 1.3, 1.5)
TIME_BINS <- 50

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
base <- ev[info$file[info$batch == BATCH & grepl("BMV_[0-9]+[A-Z]", info$file)]]
ch <- setdiff(colnames(base[[1]]), "Time")
marker_sd <- apply(do.call(rbind, base)[, ch], 2, sd)
time_bin <- as.numeric(median(vapply(base, \(m) diff(range(m[, "Time"])), numeric(1))) / TIME_BINS)

# Natural spread, in units of the pooled marker SD --------------------------------
file_med <- t(vapply(base, \(m) apply(m[, ch], 2, median), numeric(length(ch))))
file_sd <- t(vapply(base, \(m) apply(m[, ch], 2, sd), numeric(length(ch))))
spread <- tibble(
  marker = ch,
  sd_of_file_medians = apply(file_med, 2, sd) / marker_sd[ch],
  cv_of_file_sds = apply(file_sd, 2, sd) / colMeans(file_sd)
)
write.csv(spread, file.path(run_dir, "natural_spread.csv"), row.names = FALSE)
cat("median over markers of SD(file medians)/pooled SD:", round(median(spread$sd_of_file_medians), 2), "\n")
cat("median over markers of CV(file SDs):", round(median(spread$cv_of_file_sds), 2), "\n")

# Injected faults: max over markers against the injected marker only ----------------
auroc <- function(score, pos) {
  r <- rank(score)
  n1 <- sum(pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * sum(!pos))
}
jobs <- c(
  lapply(SHIFT_SD, \(s) {
    list(fault = "shift", size = s, check = "median_z", fn = \(x, m) {
      x[, m] <- x[, m] + s * marker_sd[[m]]
      x
    })
  }),
  lapply(STRETCH, \(s) {
    list(fault = "stretch", size = s, check = "sd_z", fn = \(x, m) {
      x[, m] <- mean(x[, m]) + (x[, m] - mean(x[, m])) * s
      x
    })
  })
)
long <- lapply(SEEDS, \(sd) {
  lapply(jobs, \(j) {
    set.seed(sd * 1000 + round(j$size * 100))
    n_bad <- round(FRACTION_BAD * length(base))
    bad <- seq_along(base) %in% sample(seq_along(base), n_bad)
    marker <- sample(ch, 1)
    files <- base
    files[bad] <- lapply(files[bad], j$fn, m = marker)
    names(files) <- paste0("f", seq_along(files))
    aud <- sw_fcs_signal_audit(
      files,
      ch,
      range_max = NULL,
      range_min = NULL,
      time_channel = "Time",
      time_bin = time_bin, area_height = NULL, area_height_quantile = NULL, reference_group = NULL, reference = NULL
    )
    a <- aud[aud$check == j$check, ]
    mx <- as.numeric(tapply(abs(a$value), a$file, max, na.rm = TRUE)[names(files)])
    a_own <- a[a$channel == marker, ]
    own <- abs(a_own$value)[match(names(files), a_own$file)]
    tibble(
      seed = sd,
      fault = j$fault,
      size = j$size,
      score = c("max over 40 markers", "injected marker only (oracle)"),
      auroc = c(auroc(mx, bad), auroc(own, bad))
    )
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(long, file.path(run_dir, "marker_level_long.csv"), row.names = FALSE)
mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc) / sqrt(n()),
    auroc = mean(auroc),
    .by = c(fault, size, score)
  )
write.csv(mat, file.path(run_dir, "marker_level_matrix.csv"), row.names = FALSE)
write.csv(mat, file.path(run_dir, "figure_source_data_marker_level.csv"), row.names = FALSE)
p <- ggplot(mat, aes(size, auroc, colour = score)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = auroc - ci95_half, ymax = auroc + ci95_half), size = 0.25) +
  geom_line() +
  facet_wrap(vars(fault), scales = "free_x", nrow = 1) +
  labs(
    x = "Effect size (SD added, or width factor)",
    y = "AUROC (mean, 95% CI over seeds)",
    colour = NULL
  )
ggsave(
  file.path(run_dir, "fig_marker_level.png"),
  p,
  width = 140,
  height = 70,
  units = "mm",
  dpi = 200
)
