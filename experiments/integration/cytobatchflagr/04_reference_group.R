# experiments/integration/cytobatchflagr/04_reference_group.R
# Purpose: does reference_group (robust z within the same donor, or on the residuals of a
#          two-way median polish over donor and batch) close the donor-limited gaps on real files? The same four control samples were acquired in 18 batches, so `control` is
#          the donor and each donor has 18 files. Faults go into 25% of files (batch faults into
#          all four files of one batch). Each fault is scored twice on the same injected files:
#          without reference_group (z over all 72 files), with reference_group = control (one-way,
#          arm "with") and with a data frame of file, control and batch (two-way, arm "two_way",
#          file-level faults only: a two-way design would remove the batch effect the batch fault is).
#          AUROC ranks faulty files above clean ones; the hit marker's own row is scored, as a
#          user ranks per marker. The batch fault is scored as the hit marker among all markers.
# Run:     Rscript experiments/integration/cytobatchflagr/04_reference_group.R <extract_run_dir> <run_dir>
# Inputs:  <extract_run_dir>/events.rds, files.csv, marker_names.rds   from 01_extract.R
# Outputs: <run_dir>/reference_group_long.csv                one row per seed, fault, check, arm
#          <run_dir>/figure_source_data_reference_group.csv  mean AUROC and 95% CI per fault, check, arm
#          <run_dir>/fig_reference_group.png
#          <run_dir>/session_info.txt
#
# Choices reused from earlier scripts (approved by the user): fault sizes and seeds of
# aurora_exvivo/11_limits_real.R and 12_composition_real.R; cluster settings of
# 12_composition_real.R; 25 shared named markers (G9 says 30; the code keeps the 25 named ones) and asinh(x / 6000) as in 02_batch_cd3.R;
# 2000 events per file for the batch scores as in 02_batch_cd3.R; batch fault as in
# aurora_exvivo/09_batch_real.R. Marker SD for the fault sizes is taken over all files.
# The bar from the user: mean AUROC >= 0.90 and lower CI >= 0.80.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
extract_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) stop("run dir exists: ", run_dir, call. = FALSE)
dir.create(run_dir, recursive = TRUE)

# Settings (illustrative, reused from the scripts named above) ----------------------------------
SEEDS <- 1:5
FRACTION_BAD <- 0.25
SPLIT_SD <- 0.5 # half-distance between the two modes, in marker SD
WIDEN <- c(1.1, 1.2) # width factor on every marker
REMOVE_SHARE <- 0.5 # share of the hit cluster's events removed
BATCH_SPLIT_SD <- 0.5 # two-mode split of one marker in every file of one batch
COFACTOR <- 6000 # arcsinh scale, as in 02_batch_cd3.R
CLUSTER_K <- 6
CLUSTER_EVENTS <- 2000
CLUSTER_NSTART <- 3
MIN_CLUSTER_EVENTS <- 30
EMD_EVENTS <- 2000
EMD_POOL <- 20000 # pooled events: over all files, or per group with reference_group
BATCH_EVENTS <- 2000

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
donor <- setNames(as.character(info$control), ids)
batch <- setNames(as.character(info$batch), ids)
TWO_WAY <- data.frame(file = ids, donor = as.character(info$control), batch = as.character(info$batch))
ARMS <- list(without = NULL, with = donor, two_way = TWO_WAY)

marker_sd <- apply(do.call(rbind, base), 2, sd)
pooled <- do.call(rbind, base)
set.seed(1)
km <- stats::kmeans(
  pooled[sample.int(nrow(pooled), 20000), ],
  CLUSTER_K,
  nstart = CLUSTER_NSTART,
  iter.max = 100
)
smallest <- which.min(km$size)
largest <- which.max(km$size)
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
# Robust z over all files (g = NULL), within group g, or (g a data frame) on the residuals of a
# two-way median polish over its two factor columns (maxiter 10; one file per cell here).
robust_z <- function(x, g) {
  z <- function(v) (v - median(v)) / (1.4826 * mad(v, constant = 1))
  if (is.null(g)) {
    return(z(x))
  }
  if (is.data.frame(g)) {
    tab <- tapply(x, list(g$donor, g$batch), mean)
    mp <- stats::medpolish(tab, maxiter = 10, trace.iter = FALSE, na.rm = TRUE)
    resid <- x - mp$overall - mp$row[g$donor] - mp$col[g$batch]
    return(z(resid))
  }
  out <- numeric(length(x))
  for (lv in unique(g)) out[g == lv] <- z(x[g == lv])
  out
}
top <- \(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)

# Faults ----------------------------------------------------------------------------------------

grid <- bind_rows(
  tibble(fault = "two_mode_split", size = SPLIT_SD),
  tibble(fault = "all_widened", size = WIDEN),
  tibble(fault = "composition_smallest", size = REMOVE_SHARE),
  tibble(fault = "composition_largest", size = REMOVE_SHARE),
  tibble(fault = "batch_split", size = BATCH_SPLIT_SD)
)
long <- lapply(SEEDS, function(seed) {
  lapply(seq_len(nrow(grid)), function(g) {
    fault <- grid$fault[g]
    size <- grid$size[g]
    set.seed(seed * 1000 + g)
    mk <- sample(markers, 1) # hit marker
    files <- base
    if (fault == "batch_split") {
      hit_batch <- sample(unique(batch), 1)
      bad <- batch == hit_batch
      files[bad] <- lapply(files[bad], function(x) {
        x[, mk] <- x[, mk] + sample(c(-size, size), nrow(x), replace = TRUE) * marker_sd[[mk]]
        x
      })
    } else {
      bad <- seq_along(base) %in% sample(seq_along(base), round(FRACTION_BAD * length(base)))
      files[bad] <- lapply(files[bad], function(x) {
        if (fault == "two_mode_split") {
          x[, mk] <- x[, mk] + sample(c(-size, size), nrow(x), replace = TRUE) * marker_sd[[mk]]
        } else if (fault == "all_widened") {
          med <- apply(x, 2, median)
          x <- sweep(sweep(x, 2, med) * size, 2, med, "+")
        } else {
          target <- if (fault == "composition_smallest") smallest else largest
          # Logical mask: x[-integer(0), ] would drop every row when nothing is removed.
          drop_event <- nearest(x) == target & stats::runif(nrow(x)) < size
          x <- x[!drop_event, , drop = FALSE]
        }
        x
      })
    }
    bad <- setNames(bad, names(files))
    arms <- if (fault == "batch_split") names(ARMS)[1:2] else names(ARMS)
    lapply(arms, function(arm) {
      grp <- ARMS[[arm]]
      if (fault == "batch_split") {
        sc <- sw_batch_scores(files, batch, markers, BATCH_EVENTS, seed, order = NULL, reference_group = grp, population = NULL)
        per_marker <- list(
          pair_emd = sc[sc$check == "pair_emd", ] |>
            summarise(v = max(value), .by = marker),
          spread_ratio = sc[sc$check == "spread_ratio", ] |>
            summarise(v = max(abs(log(value))), .by = marker)
        )
        res <- vapply(per_marker, \(d) auroc(d$v[match(markers, d$marker)], markers == mk), numeric(1))
      } else if (fault == "two_mode_split") {
        sa <- sw_fcs_signal_audit(files, markers, NULL, NULL, time_channel = NULL, time_bin = NULL,
          area_height = NULL, area_height_quantile = NULL, reference_group = grp, reference = NULL)
        sd_z <- sa[sa$check == "sd_z" & sa$channel == mk, ]
        af <- sw_anomaly_features(
          lapply(files, \(x) x[, mk, drop = FALSE]),
          channels = mk,
          features = "emd",
          quantile_step = NULL,
          n_events = EMD_EVENTS,
          seed = seed,
          scale = "asinh(x / 6000)",
          aggregate_from = "all",
          agg_size = EMD_POOL,
          reference = NULL,
          reference_group = grp
        )
        emd <- af[[paste0(mk, "_emd")]][match(names(files), af$file)]
        res <- c(
          sd_z = auroc(abs(sd_z$value[match(names(files), sd_z$file)]), bad),
          emd_z = auroc(abs(robust_z(emd, if (is.data.frame(grp)) grp[match(names(files), grp$file), ] else grp[names(files)])), bad)
        )
      } else if (fault == "all_widened") {
        sa <- sw_fcs_signal_audit(files, markers, NULL, NULL, time_channel = NULL, time_bin = NULL,
          area_height = NULL, area_height_quantile = NULL, reference_group = grp, reference = NULL)
        gz <- sa[sa$check == "global_sd_z", ]
        res <- c(global_sd_z = auroc(abs(gz$value[match(names(files), gz$file)]), bad))
      } else {
        cs <- sw_cluster_scores(files, markers, CLUSTER_K, CLUSTER_EVENTS, seed, CLUSTER_NSTART,
          MIN_CLUSTER_EVENTS, reference = NULL, reference_group = grp)
        # The row sw_diagnose ranks: the found cluster nearest the target cluster.
        target <- if (fault == "composition_smallest") smallest else largest
        hit_found <- which.min(rowSums(sweep(attr(cs, "centers"), 2, km$centers[target, ])^2))
        fz <- cs[cs$check == "freq_z" & cs$cluster == hit_found, ]
        res <- c(freq_z_hitrow = auroc(abs(fz$value[match(names(files), fz$file)]), bad))
      }
      tibble(seed = seed, fault = fault, size = size, check = names(res), arm = arm, auroc = unname(res))
    }) |>
      bind_rows()
  }) |>
    bind_rows()
}) |>
  bind_rows()
write.csv(long, file.path(run_dir, "reference_group_long.csv"), row.names = FALSE)

# Summary ---------------------------------------------------------------------------------------

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(fault, size, check, arm)
  ) |>
  mutate(
    case = paste(fault, size),
    arm = factor(arm, levels = c("without", "with", "two_way")),
    meets_bar = auroc >= 0.9 & auroc - ci95_half >= 0.8
  )
write.csv(mat, file.path(run_dir, "figure_source_data_reference_group.csv"), row.names = FALSE)
print(as.data.frame(mutate(mat, across(c(auroc, ci95_half), \(x) round(x, 3)))[c("case", "check", "arm", "auroc", "ci95_half", "meets_bar")]))

p <- ggplot(mat, aes(auroc, check, colour = arm)) +
  geom_vline(xintercept = c(0.5, 0.9), linetype = 2, colour = "grey60") +
  geom_point(
    data = mutate(long, case = paste(fault, size), arm = factor(arm, levels = c("without", "with", "two_way"))),
    position = position_dodge(width = 0.6), size = 0.4, alpha = 0.5
  ) +
  geom_pointrange(
    aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)),
    position = position_dodge(width = 0.6), size = 0.2
  ) +
  facet_wrap(vars(case), scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(without = "grey40", with = "#D55E00", two_way = "#0072B2")) +
  labs(
    x = "AUROC on 72 real files (mean and 95% CI over 5 seeds; dots = seeds)",
    y = NULL,
    colour = "reference_group"
  )
ggsave(file.path(run_dir, "fig_reference_group.png"), p, width = 150, height = 130, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
