# experiments/integration/aurora_exvivo/12_composition_real.R
# Purpose: the last two coverage rows with no real-file evidence. (1) Composition change: a share of
#          the events of the smallest k-means cluster is removed. (2) All markers widened: every
#          channel stretched around the file median. Faults go into 25% of one batch's real donor
#          files; AUROC ranks faulty files above clean ones.
# Run:     Rscript experiments/integration/aurora_exvivo/12_composition_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/composition_real_long.csv
#          <run_dir>/fig_composition_real.png, figure_source_data_composition_real.csv
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
SEEDS <- 1:5
FRACTION_BAD <- 0.25
REMOVE_SHARE <- c(0.5, 0.8) # share of the smallest cluster's events removed
WIDEN <- c(1.1, 1.2) # width factor on every channel
CLUSTER_K <- 6 # as in 10_file_faults_real.R
CLUSTER_EVENTS <- 2000

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
base <- ev[info$file[info$batch == BATCH & grepl("BMV_[0-9]+[A-Z]", info$file)]]
names(base) <- paste0("f", seq_along(base))
ch <- setdiff(colnames(base[[1]]), "Time")
pooled <- do.call(rbind, base)
set.seed(1)
km <- stats::kmeans(pooled[sample.int(nrow(pooled), 20000), ch], CLUSTER_K, nstart = 3, iter.max = 100)
smallest <- which.min(km$size)
largest <- which.max(km$size)
print(round(km$size / sum(km$size), 3)) # cluster shares in the clean pool
nearest <- function(x) {
  d2 <- outer(rowSums(x^2), rowSums(km$centers^2), "+") - 2 * tcrossprod(x, km$centers)
  max.col(-d2, ties.method = "first")
}
auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}
top <- \(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)

grid <- bind_rows(
  tibble(fault = "composition_smallest", size = REMOVE_SHARE),
  tibble(fault = "composition_largest", size = REMOVE_SHARE),
  tibble(fault = "all_widened", size = WIDEN)
)
long <- lapply(SEEDS, function(seed) {
  lapply(seq_len(nrow(grid)), function(g) {
    fault <- grid$fault[g]
    size <- grid$size[g]
    set.seed(seed * 1000 + g)
    bad <- seq_along(base) %in% sample(seq_along(base), round(FRACTION_BAD * length(base)))
    files <- base
    files[bad] <- lapply(files[bad], function(x) {
      if (startsWith(fault, "composition")) {
        target <- if (fault == "composition_smallest") smallest else largest
        # Logical mask: x[-integer(0), ] would drop every row when nothing is removed.
        drop_event <- nearest(x[, ch]) == target & stats::runif(nrow(x)) < size
        x[!drop_event, , drop = FALSE]
      } else {
        med <- apply(x[, ch], 2, stats::median)
        x[, ch] <- sweep(sweep(x[, ch], 2, med) * size, 2, med, "+")
        x
      }
    })
    b <- sample(ch, 1) # one marker's row, as a user ranks per marker (Q14)
    sa <- sw_fcs_signal_audit(files, ch, NULL, NULL, time_channel = NULL, time_bin = NULL, area_height = NULL, area_height_quantile = NULL, reference_group = NULL, reference = NULL, summary = NULL)
    sd_z <- sa[sa$check == "sd_z" & sa$channel == b, ]
    cs <- sw_cluster_scores(files, ch, CLUSTER_K, CLUSTER_EVENTS, seed, 3, 30, reference = NULL, reference_group = NULL)
    fz <- cs[cs$check == "freq_z", ]
    gz <- sa[sa$check == "global_sd_z", ]
    scores <- list(
      sd_z_one_marker = abs(sd_z$value[match(names(files), sd_z$file)]),
      global_sd_z = abs(gz$value[match(names(files), gz$file)]),
      cluster_freq_z = as.numeric(tapply(abs(fz$value), fz$file, top)[names(files)])
    )
    if (startsWith(fault, "composition")) {
      # The row sw_diagnose ranks: the found cluster nearest the target cluster.
      target <- if (fault == "composition_smallest") smallest else largest
      hit_found <- which.min(rowSums(sweep(attr(cs, "centers"), 2, km$centers[target, ])^2))
      fz_hit <- fz[fz$cluster == hit_found, ]
      scores$cluster_freq_z_hitrow <- abs(fz_hit$value[match(names(files), fz_hit$file)])
    }
    tibble(seed = seed, fault = fault, size = size, check = names(scores),
      auroc = vapply(scores, \(v) auroc(v, bad), numeric(1)))
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(long, file.path(run_dir, "composition_real_long.csv"), row.names = FALSE)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(fault, size, check)
  ) |>
  mutate(case = paste(fault, size))
write.csv(mat, file.path(run_dir, "figure_source_data_composition_real.csv"), row.names = FALSE)
print(as.data.frame(mutate(mat, across(c(auroc, ci95_half), \(x) round(x, 2)))))
p <- ggplot(mat, aes(auroc, check)) +
  geom_vline(xintercept = c(0.5, 0.9), linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)), size = 0.2) +
  facet_wrap(vars(case), nrow = 1) +
  labs(x = "AUROC on real files (mean, 95% CI over 5 seeds)", y = NULL)
ggsave(file.path(run_dir, "fig_composition_real.png"), p, width = 180, height = 50, units = "mm", dpi = 200)
