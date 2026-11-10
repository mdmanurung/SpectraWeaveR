# experiments/integration/aurora_exvivo/11_limits_real.R
# Purpose: detection limits on real files for the two fault classes that missed the AUROC bar at
#          one size in 10_file_faults_real.R: a marker split into two modes, and a shift inside
#          the smallest cluster. Sweeps the fault size and scores sd_z, a per-channel EMD to the
#          pooled cohort (robust z over files) and, for the rare population, cluster_median_z.
# Run:     Rscript experiments/integration/aurora_exvivo/11_limits_real.R <run_dir> [<out_dir>]
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <out_dir>/limits_real_long.csv     (<out_dir> defaults to <run_dir>; must not exist if given)
#          <out_dir>/fig_limits_real.png, figure_source_data_limits_real.csv
#          bimodality_z arm added for G20-G23.
#
# Every value below is illustrative; none is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[[1]]
out_dir <- if (length(args) > 1) args[[2]] else run_dir
if (out_dir != run_dir) {
  if (dir.exists(out_dir)) stop("out dir exists: ", out_dir, call. = FALSE)
  dir.create(out_dir, recursive = TRUE)
}
BATCH <- "Batch1"
SEEDS <- 1:5
FRACTION_BAD <- 0.25
SPLIT_SD <- c(0.5, 1, 1.5) # half-distance between the two modes, in marker SD
RARE_SD <- c(1, 2, 3) # shift inside the smallest cluster, in marker SD
CLUSTER_K <- 6 # as in 10_file_faults_real.R
CLUSTER_EVENTS <- 2000
EMD_EVENTS <- as.integer(Sys.getenv("SW_EMD_EVENTS", "2000")) # files hold 5000 events
EMD_POOL <- 20000
ONLY_FAULT <- Sys.getenv("SW_ONLY_FAULT", "") # "" = all faults
SUFFIX <- if (EMD_EVENTS == 2000 && ONLY_FAULT == "") "" else paste0("_emd", EMD_EVENTS, ONLY_FAULT)

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

grid <- bind_rows(
  tibble(fault = "bimodal", size = SPLIT_SD),
  tibble(fault = "rare_pop", size = RARE_SD)
)
if (ONLY_FAULT != "") grid <- grid[grid$fault == ONLY_FAULT, ]
long <- lapply(SEEDS, function(seed) {
  lapply(seq_len(nrow(grid)), function(g) {
    fault <- grid$fault[g]
    size <- grid$size[g]
    set.seed(seed * 1000 + g)
    bad <- seq_along(base) %in% sample(seq_along(base), round(FRACTION_BAD * length(base)))
    b <- sample(ch, 1)
    files <- base
    files[bad] <- lapply(files[bad], function(x) {
      if (fault == "bimodal") {
        x[, b] <- x[, b] + sample(c(-size, size), nrow(x), replace = TRUE) * marker_sd[[b]]
      } else {
        hit <- nearest(x[, ch]) == smallest
        x[hit, b] <- x[hit, b] + size * marker_sd[[b]]
      }
      x
    })
    sa <- sw_fcs_signal_audit(files, ch, NULL, NULL, time_channel = NULL, time_bin = NULL, area_height = NULL, area_height_quantile = NULL, reference_group = NULL, reference = NULL)
    sd_z <- sa[sa$check == "sd_z" & sa$channel == b, ]
    bz <- sa[sa$check == "bimodality_z" & sa$channel == b, ]
    af <- sw_anomaly_features(
      lapply(files, \(x) x[, b, drop = FALSE]),
      channels = b,
      features = "emd",
      quantile_step = NULL,
      n_events = EMD_EVENTS,
      seed = seed,
      scale = "OMIQ-transformed",
      aggregate_from = "all",
      agg_size = EMD_POOL,
      reference = NULL, reference_group = NULL
    )
    scores <- list(
      sd_z = abs(sd_z$value[match(names(files), sd_z$file)]),
      bimodality_z = abs(bz$value[match(names(files), bz$file)]),
      # G20 fallback: Hartigan dip of the hit marker (diptest, GPL, called not copied), robust z.
      # diptest is not a package dependency; without it the arm is NA.
      dip_z = if (requireNamespace("diptest", quietly = TRUE)) {
        abs(robust_z(vapply(files, \(x) diptest::dip(x[, b]), numeric(1))))
      } else {
        rep(NA_real_, length(files))
      },
      emd_z = abs(robust_z(af[[paste0(b, "_emd")]][match(names(files), af$file)]))
    )
    if (fault == "rare_pop") {
      cs <- sw_cluster_scores(files, ch, CLUSTER_K, CLUSTER_EVENTS, seed, 3, 30, reference = NULL, reference_group = NULL)
      d <- cs[cs$check == "median_z" & cs$marker == b, ]
      top <- tapply(abs(d$value), d$file, \(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
      scores$cluster_median_z <- as.numeric(top[names(files)])
    }
    tibble(seed = seed, fault = fault, size = size, check = names(scores),
      auroc = vapply(scores, \(v) auroc(v, bad), numeric(1)))
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(long, file.path(out_dir, paste0("limits_real_long", SUFFIX, ".csv")), row.names = FALSE)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(fault, size, check)
  )
write.csv(mat, file.path(out_dir, paste0("figure_source_data_limits_real", SUFFIX, ".csv")), row.names = FALSE)
print(as.data.frame(mutate(mat, across(c(auroc, ci95_half), \(x) round(x, 2)))))
p <- ggplot(mat, aes(size, auroc, colour = check)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = auroc - ci95_half, ymax = pmin(auroc + ci95_half, 1)), size = 0.2) +
  geom_line() +
  facet_wrap(vars(fault), scales = "free_x") +
  labs(x = "Fault size (marker SD)", y = "AUROC on real files (mean, 95% CI over 5 seeds)", colour = NULL)
ggsave(file.path(out_dir, paste0("fig_limits_real", SUFFIX, ".png")), p, width = 150, height = 70, units = "mm", dpi = 200)
