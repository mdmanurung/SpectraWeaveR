# experiments/integration/aurora_exvivo/13_rare_oracle_real.R
# Purpose: test why cluster_median_z had wide intervals for a shift in the smallest cluster
#          (11_limits_real.R). Score the same faults two ways: with sw_cluster_scores() (clusters
#          found on the pooled, partly faulty files) and with the injection's own clusters fixed
#          (oracle: per-file median of the hit marker among events nearest the smallest reference
#          centre, robust z over files). If the oracle is near 1 where sw_cluster_scores is not,
#          the weak point is cluster matching, not the statistic.
# Run:     Rscript experiments/integration/aurora_exvivo/13_rare_oracle_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/rare_oracle_long.csv
#          <run_dir>/fig_rare_oracle.png, figure_source_data_rare_oracle.csv
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
RARE_SD <- c(1, 2, 3) # as in 11_limits_real.R
CLUSTER_K <- 6
CLUSTER_EVENTS <- 2000

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
base <- ev[info$file[info$batch == BATCH & grepl("BMV_[0-9]+[A-Z]", info$file)]]
names(base) <- paste0("f", seq_along(base))
ch <- setdiff(colnames(base[[1]]), "Time")
pooled <- do.call(rbind, base)
marker_sd <- apply(pooled[, ch], 2, sd)
set.seed(1)
km <- stats::kmeans(pooled[sample.int(nrow(pooled), 20000), ch], CLUSTER_K, nstart = 3, iter.max = 100)
smallest <- which.min(km$size)
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
robust_z <- function(x) (x - median(x)) / (1.4826 * mad(x, constant = 1))
top <- \(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
# Cluster labels of the clean files, fixed before any fault (the oracle's clusters).
clean_hit <- lapply(base, \(x) nearest(x[, ch]) == smallest)

long <- lapply(SEEDS, function(seed) {
  lapply(seq_along(RARE_SD), function(g) {
    size <- RARE_SD[g]
    set.seed(seed * 1000 + g + 3) # +3: same draws as the rare_pop rows of 11_limits_real.R
    bad <- seq_along(base) %in% sample(seq_along(base), round(FRACTION_BAD * length(base)))
    b <- sample(ch, 1)
    files <- base
    for (i in which(bad)) {
      files[[i]][clean_hit[[i]], b] <- files[[i]][clean_hit[[i]], b] + size * marker_sd[[b]]
    }
    cs <- sw_cluster_scores(files, ch, CLUSTER_K, CLUSTER_EVENTS, seed, 3, 30, reference = NULL, reference_group = NULL)
    d <- cs[cs$check == "median_z" & cs$marker == b, ]
    oracle <- vapply(seq_along(files), \(i) stats::median(files[[i]][clean_hit[[i]], b]), numeric(1))
    # Which found cluster best matches the injected one: share of its events that are injected.
    found_share <- attr(cs, "centers") |>
      (\(cent) {
        hit_centre <- km$centers[smallest, ]
        min(sqrt(rowSums(sweep(cent, 2, hit_centre)^2)))
      })()
    # The package fix: clusters fitted on known-good files only (here the clean files), the
    # rest assigned to the nearest centre.
    cs_ref <- sw_cluster_scores(files, ch, CLUSTER_K, CLUSTER_EVENTS, seed, 3, 30,
      reference = names(files)[!bad], reference_group = NULL)
    d_ref <- cs_ref[cs_ref$check == "median_z" & cs_ref$marker == b, ]
    # The row a user would read: the found cluster whose centre is nearest the injected one.
    hit_found <- which.min(rowSums(sweep(attr(cs, "centers"), 2, km$centers[smallest, ])^2))
    d_hit <- d[d$cluster == hit_found, ]
    scores <- list(
      cluster_median_z_hitrow = abs(d_hit$value[match(names(files), d_hit$file)]),
      cluster_median_z = as.numeric(tapply(abs(d$value), d$file, top)[names(files)]),
      cluster_median_z_ref = as.numeric(tapply(abs(d_ref$value), d_ref$file, top)[names(files)]),
      oracle_cluster_z = abs(robust_z(oracle))
    )
    tibble(seed = seed, size = size, check = names(scores),
      auroc = vapply(scores, \(v) auroc(v, bad), numeric(1)),
      centre_distance = found_share)
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(long, file.path(run_dir, "rare_oracle_long.csv"), row.names = FALSE)
print(as.data.frame(long))

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc) / sqrt(n()),
    auroc = mean(auroc),
    .by = c(size, check)
  )
write.csv(mat, file.path(run_dir, "figure_source_data_rare_oracle.csv"), row.names = FALSE)
p <- ggplot(long, aes(factor(size), auroc, colour = check)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = 2, colour = "grey60") +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5, seed = 1), size = 1) +
  geom_point(data = mat, position = position_dodge(width = 0.5), shape = 95, size = 6) +
  labs(x = "Shift in the smallest cluster (marker SD)", y = "AUROC on real files (points: seeds; bar: mean)", colour = NULL)
ggsave(file.path(run_dir, "fig_rare_oracle.png"), p, width = 120, height = 75, units = "mm", dpi = 200)
