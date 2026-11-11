# Purpose: three user-approved metrics (2026-11-10) added to the ranking, on the same paired events
#          for every arm. Negatives and positives come from the reference arm's cut only.
#          1. Primary spread: Nguyen 2013 robust SD of the reference negatives, linear scale,
#             q84 - q50. The 2.5-97.5 outer-quantile spread stays as a check column.
#          2. Quadrature spillover-spread SS_ij = sqrt(max(0, rSD_j(i+ j-)^2 - rSD_j(i- j-)^2)) /
#             sqrt(median_i(i+) - median_i(i-)), rSD = Nguyen, all populations from the reference cuts.
#          3. Tilt: Spearman rho for every channel pair among events negative for both channels.
#          Columns of 28 (MAD rSD, pile-up index, hypernegatives, stain index, TNF-CD8 rho) are
#          recomputed unchanged. One figure per metric is written to the run dir.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/30_rank_metrics2.R <src_run_dir> <new_run_dir>
#          from the repo root. <new_run_dir> must not exist.
# Inputs:  <src_run_dir>/manifest.csv, shared_events.rds [01_manifest-gate.R]
#          <src_run_dir>/mindensity_cuts.csv  reference-arm cuts [02_metrics.R]
#          experiments/bmv_test/compare_3D/data/unmix_cells/<EXTRA_ARMS>/*.fcs [25_joint_sweep.R, 22_fix_my_unmix.R]
# Outputs: <new_run_dir>/manifest.csv          all arms scored
#          <new_run_dir>/reference_cuts.csv    per sample, channel: reference cut and pinned flag
#          <new_run_dir>/rank_channel.csv      per arm, sample, channel: pinned, n_neg, n_pos, rsd_nguyen (primary),
#                                              spread_q (check), rsd_mad, pileup_index, pct_hyperneg,
#                                              stain_index (Nguyen), stain_index_q
#          <new_run_dir>/rank_ssm.csv          per arm, sample, source i, target j: ss, n of each population, low_n
#          <new_run_dir>/rank_tilt.csv         per arm, sample, channel pair: n_both_neg, rho
#          <new_run_dir>/rank_streak.csv       per arm, sample: TNF-CD8 rho among double negatives (as in 28)
#          <new_run_dir>/summary_spread.csv    per arm, sample: median rsd_nguyen over unpinned channels
#          <new_run_dir>/summary_ssm.csv       per arm, sample: median ss over unpinned i, j
#          <new_run_dir>/summary_tilt.csv      per arm, sample: median |rho| over pairs
#          <new_run_dir>/tilt_top10_as188_joint.csv  10 pairs with largest median |rho| over donors in the reference arm
#          <new_run_dir>/fig_spread.png, fig_ssm.png, fig_tilt.png and fig_*_source.csv (plotted data)
#          <new_run_dir>/sessionInfo_30.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
REFERENCE_ARM <- "as188_joint" # negatives and positives from this arm's cut only (user, 2026-11-10)
NGUYEN_P <- c(0.5, 0.84) # Nguyen 2013 robust SD = q84 - q50 (user, 2026-11-10)
QUANTILES <- c(0.025, 0.975) # outer quantiles, check column (user, 2026-11-09)
Z_WIDTH <- 3.92 # width of the 2.5-97.5 interval in SD units under normality
HYPERNEG_K <- 3 # hypernegative: below -K x reference rSD (user, grill Q19, 2026-11-09)
PIN_TOL <- 0.02 # a cut within this of the window edge counts as pinned (report rule, 2026-11-09)
CUT_WINDOW <- c(0, 3) # mindensity window used for the cuts in 02 (user, 2026-11-09)
MIN_N <- 100 # NA for SS or rho when a population has fewer events; 28 sets none, 100 is the README's
#               P9 minimum for shared positives. NOT user-set: needs confirmation (2026-11-10)
N_TOP_PAIRS <- 10 # tilt pairs listed for the reference arm (user, 2026-11-10)
N_CORES <- 5 # one worker per donor
STREAK <- c(TNF = "BV510-A", CD8 = "SBUV445-A") # pair from the stripe bisect (2026-11-10)
EXTRA_ARMS <- c(
  "repro188",
  "jpr_false",
  "coll03",
  "coll07",
  "alpha025",
  "alpha075",
  "sig188",
  "fix188"
) # 25_joint_sweep.R, 22_fix_my_unmix.R
EXTRA_DIR <- "experiments/bmv_test/compare_3D/data/unmix_cells"

args <- commandArgs(trailingOnly = TRUE)
src_dir <- args[[1]]
new_dir <- args[[2]]
stopifnot("new run dir must not exist" = !dir.exists(new_dir))
dir.create(new_dir, recursive = TRUE)

# Manifest: original arms plus extra arms (same raw files, same row order) ---------------------
manifest_src <- read_csv(
  file.path(src_dir, "manifest.csv"),
  show_col_types = FALSE
)
manifest_extra <- EXTRA_ARMS |>
  map(\(a) {
    tibble(
      arm = a,
      dir = file.path("unmix_cells", a),
      viability = "Zombie NIR-A",
      file = list.files(
        file.path(EXTRA_DIR, a),
        pattern = "\\.fcs$",
        full.names = TRUE
      )
    )
  }) |>
  list_rbind() |>
  mutate(sample = stringr::str_extract(basename(file), "^[A-Z][0-9]+"))
manifest <- bind_rows(manifest_src, manifest_extra)
stopifnot("5 files per arm" = all(count(manifest, arm)$n == 5))
stopifnot(
  "5 distinct donors per arm" = all(
    count(distinct(manifest, arm, sample), arm)$n == 5
  )
)
write_csv(manifest, file.path(new_dir, "manifest.csv"))

# Reference cuts, pinned flag -------------------------------------------------------------
shared <- readRDS(file.path(src_dir, "shared_events.rds"))
ref_cuts <- read_csv(
  file.path(src_dir, "mindensity_cuts.csv"),
  show_col_types = FALSE
) |>
  filter(arm == REFERENCE_ARM) |>
  # Pinned cuts are excluded per sample x channel (user, grill Q21, 2026-11-09). Tilt and the
  # TNF-CD8 streak keep raw cuts: they need no valley.
  mutate(
    pinned = is.na(cut) |
      abs(cut - CUT_WINDOW[[1]]) <= PIN_TOL |
      abs(cut - CUT_WINDOW[[2]]) <= PIN_TOL
  )
write_csv(ref_cuts, file.path(new_dir, "reference_cuts.csv"))
channels <- unique(ref_cuts$channel)
stopifnot("no missing reference cut" = !anyNA(ref_cuts$cut))
message(
  "Pinned sample x channel cells excluded from spread, SI and SS: ",
  sum(ref_cuts$pinned),
  " of ",
  nrow(ref_cuts)
)

# Nguyen robust SD, linear scale ----------------------------------------------------------
rsd_nguyen <- function(x) {
  q <- quantile(x, NGUYEN_P, names = FALSE)
  q[[2]] - q[[1]]
}

# One sample: reference populations, then every metric for every arm ----------------------
process_sample <- function(s) {
  cut_s <- ref_cuts |>
    filter(sample == s) |>
    pull(cut, name = channel)
  pinned_s <- ref_cuts |>
    filter(sample == s) |>
    pull(pinned, name = channel)
  lin <- manifest |>
    filter(sample == s) |>
    (\(m) set_names(m$file, m$arm))() |>
    map(\(f) {
      flowCore::exprs(sw_io_read_fcs(f))[shared[[s]], channels, drop = FALSE]
    })
  stopifnot(
    "same event count in every arm" = n_distinct(map_int(lin, nrow)) == 1
  )
  ref_asinh <- asinh(lin[[REFERENCE_ARM]] / COFACTOR)
  neg <- sweep(ref_asinh, 2, cut_s[channels], "<")
  pos <- sweep(ref_asinh, 2, cut_s[channels], ">")
  ref_rsd <- map_dbl(set_names(channels), \(ch) {
    1.4826 * mad(lin[[REFERENCE_ARM]][neg[, ch], ch], constant = 1)
  })
  usable <- channels[!pinned_s[channels]]
  pair_idx <- combn(channels, 2, simplify = FALSE)

  # Per channel: spreads, pile-up, hypernegatives, stain index
  per_channel <- lin |>
    imap(\(m, a) {
      map(channels, \(ch) {
        x <- m[neg[, ch], ch]
        q <- quantile(x, QUANTILES, names = FALSE)
        spread_q <- (q[[2]] - q[[1]]) / Z_WIDTH
        rsd_mad <- 1.4826 * mad(x, constant = 1)
        rsd_ng <- rsd_nguyen(x)
        delta <- median(m[pos[, ch], ch]) - median(x)
        tibble(
          arm = a,
          sample = s,
          channel = ch,
          pinned = pinned_s[[ch]],
          n_neg = length(x),
          n_pos = sum(pos[, ch]),
          rsd_nguyen = rsd_ng,
          spread_q = spread_q,
          rsd_mad = rsd_mad,
          pileup_index = rsd_mad / spread_q,
          pct_hyperneg = 100 * mean(x < -HYPERNEG_K * ref_rsd[[ch]]),
          stain_index = delta / (2 * rsd_ng),
          stain_index_q = delta / (2 * spread_q)
        )
      }) |>
        list_rbind()
    }) |>
    list_rbind()

  # Quadrature spillover-spread, i and j both unpinned, i != j
  ssm <- lin |>
    imap(\(m, a) {
      map(usable, \(i) {
        den <- median(m[pos[, i], i]) - median(m[neg[, i], i])
        map(setdiff(usable, i), \(j) {
          x_ipjn <- m[pos[, i] & neg[, j], j]
          x_injn <- m[neg[, i] & neg[, j], j]
          tibble(
            arm = a,
            sample = s,
            source = i,
            target = j,
            n_ip = sum(pos[, i]),
            n_in = sum(neg[, i]),
            n_ipjn = length(x_ipjn),
            n_injn = length(x_injn),
            low_n = min(
              sum(pos[, i]),
              sum(neg[, i]),
              length(x_ipjn),
              length(x_injn)
            ) <
              MIN_N,
            delta_i = den,
            ss = if (low_n || den <= 0) {
              NA_real_
            } else {
              sqrt(max(0, rsd_nguyen(x_ipjn)^2 - rsd_nguyen(x_injn)^2)) /
                sqrt(den)
            }
          )
        }) |>
          list_rbind()
      }) |>
        list_rbind()
    }) |>
    list_rbind()

  # Tilt: Spearman rho for every pair among events negative for both channels
  tilt <- lin |>
    imap(\(m, a) {
      map(pair_idx, \(p) {
        both <- neg[, p[[1]]] & neg[, p[[2]]]
        tibble(
          arm = a,
          sample = s,
          channel_a = p[[1]],
          channel_b = p[[2]],
          n_both_neg = sum(both),
          rho = if (sum(both) < MIN_N) {
            NA_real_
          } else {
            cor(m[both, p[[1]]], m[both, p[[2]]], method = "spearman")
          }
        )
      }) |>
        list_rbind()
    }) |>
    list_rbind()

  streak <- tilt |>
    filter(
      channel_a %in% STREAK,
      channel_b %in% STREAK
    ) |>
    transmute(arm, sample, n_both_neg, rho_tnf_cd8 = rho)
  list(per_channel = per_channel, ssm = ssm, tilt = tilt, streak = streak)
}

# Run all donors in parallel; any worker error stops the script ---------------------------
samples <- sort(unique(manifest$sample))
results <- parallel::mclapply(samples, process_sample, mc.cores = N_CORES) |>
  set_names(samples)
failed <- map_lgl(results, \(r) inherits(r, "try-error") || !is.list(r))
if (any(failed)) {
  stop(
    "worker failed for: ",
    paste(names(results)[failed], collapse = ", "),
    "\n",
    paste(map_chr(results[failed], as.character), collapse = "\n"),
    call. = FALSE
  )
}

rank_channel <- results |> map("per_channel") |> list_rbind()
rank_ssm <- results |> map("ssm") |> list_rbind()
rank_tilt <- results |> map("tilt") |> list_rbind()
rank_streak <- results |> map("streak") |> list_rbind()

# Summaries per arm and donor, on unpinned channels where the metric has a valley ------------
summary_spread <- rank_channel |>
  filter(!pinned) |>
  summarise(
    median_rsd_nguyen = median(rsd_nguyen),
    median_spread_q = median(spread_q),
    n_channels = n(),
    .by = c(arm, sample)
  )
summary_ssm <- rank_ssm |>
  summarise(
    median_ss = median(ss, na.rm = TRUE),
    n_cells = n(),
    n_na = sum(is.na(ss)),
    n_low_n = sum(low_n),
    .by = c(arm, sample)
  )
summary_tilt <- rank_tilt |>
  summarise(
    median_abs_rho = median(abs(rho), na.rm = TRUE),
    n_pairs = n(),
    n_na = sum(is.na(rho)),
    .by = c(arm, sample)
  )
# Top pairs in the reference arm: median over donors of |rho|
tilt_top <- rank_tilt |>
  filter(arm == REFERENCE_ARM) |>
  summarise(
    median_abs_rho = median(abs(rho), na.rm = TRUE),
    rho_by_donor = paste0(sample, "=", round(rho, 3), collapse = "; "),
    .by = c(channel_a, channel_b)
  ) |>
  slice_max(median_abs_rho, n = N_TOP_PAIRS, with_ties = FALSE)

# Never fail silently: no metric column may be all NA --------------------------------------
metric_cols <- list(
  rank_channel = c(
    "rsd_nguyen",
    "spread_q",
    "rsd_mad",
    "pileup_index",
    "pct_hyperneg",
    "stain_index",
    "stain_index_q"
  ),
  rank_ssm = c("ss"),
  rank_tilt = c("rho"),
  rank_streak = c("rho_tnf_cd8"),
  summary_spread = c("median_rsd_nguyen"),
  summary_ssm = c("median_ss"),
  summary_tilt = c("median_abs_rho")
)
iwalk(metric_cols, \(cols, tbl) {
  d <- get(tbl)
  walk(cols, \(col) {
    if (all(is.na(d[[col]]))) {
      stop(sprintf("%s$%s is all NA", tbl, col), call. = FALSE)
    }
  })
})
stopifnot("5 donors in every summary" = all(count(summary_tilt, arm)$n == 5))
stopifnot("10 top pairs" = nrow(tilt_top) == N_TOP_PAIRS)
message(
  "SS cells NA (low n or non-positive delta): ",
  sum(is.na(rank_ssm$ss)),
  " of ",
  nrow(rank_ssm)
)
message(
  "SS cells flagged low_n (< ",
  MIN_N,
  " events): ",
  sum(rank_ssm$low_n),
  " of ",
  nrow(rank_ssm)
)

write_csv(rank_channel, file.path(new_dir, "rank_channel.csv"))
write_csv(rank_ssm, file.path(new_dir, "rank_ssm.csv"))
write_csv(rank_tilt, file.path(new_dir, "rank_tilt.csv"))
write_csv(rank_streak, file.path(new_dir, "rank_streak.csv"))
write_csv(summary_spread, file.path(new_dir, "summary_spread.csv"))
write_csv(summary_ssm, file.path(new_dir, "summary_ssm.csv"))
write_csv(summary_tilt, file.path(new_dir, "summary_tilt.csv"))
write_csv(tilt_top, file.path(new_dir, "tilt_top10_as188_joint.csv"))

# Figures: donors as points, bar = median over the 5 donors, arms ordered by that median -------
theme_set(ggpubr::theme_pubr(base_size = 8))
donor_pal <- c(
  A3 = "#1b9e77",
  A6 = "#d95f02",
  B11 = "#7570b3",
  B2 = "#e7298a",
  D1 = "#66a61e"
)

# Plot frame: one row per arm and donor, arm as a factor ordered by its donor median
arm_frame <- function(d, value) {
  d |>
    rename(value = all_of(value)) |>
    mutate(arm = reorder(arm, value, FUN = median))
}
plot_arm_summary <- function(d, ylab) {
  med <- summarise(d, value = median(value), .by = arm)
  ggplot(d, aes(arm, value)) +
    geom_col(data = med, fill = "grey80", width = 0.7) +
    geom_point(
      aes(colour = sample),
      position = position_jitter(width = 0.15, height = 0, seed = 1),
      size = 1.2
    ) +
    scale_colour_manual(values = donor_pal, breaks = names(donor_pal)) +
    labs(x = NULL, y = ylab, colour = "Donor") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

fig_spread <- arm_frame(summary_spread, "median_rsd_nguyen")
write_csv(fig_spread, file.path(new_dir, "fig_spread_source.csv"))
p_spread <- plot_arm_summary(
  fig_spread,
  "Median over channels of Nguyen rSD\n(q84 - q50, reference negatives, linear)"
)
ggsave(
  file.path(new_dir, "fig_spread.png"),
  p_spread,
  width = 140,
  height = 90,
  units = "mm",
  dpi = 200
)

fig_ssm <- arm_frame(summary_ssm, "median_ss")
fig_ssm_heat <- rank_ssm |>
  filter(arm == REFERENCE_ARM) |>
  summarise(ss = median(ss, na.rm = TRUE), .by = c(source, target)) |>
  mutate(ss = if_else(is.nan(ss), NA_real_, ss))
write_csv(fig_ssm, file.path(new_dir, "fig_ssm_source.csv"))
write_csv(fig_ssm_heat, file.path(new_dir, "fig_ssm_heatmap_source.csv"))
p_ssm_a <- plot_arm_summary(fig_ssm, "Median SS over source i, target j")
p_ssm_b <- ggplot(fig_ssm_heat, aes(target, source, fill = ss)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "sqrt",
    na.value = "white",
    name = "SS (median\nover donors)"
  ) +
  labs(
    x = "Target j",
    y = "Source i",
    title = paste0(REFERENCE_ARM, ": unpinned channels, white = NA")
  ) +
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.text.y = element_text(size = 5)
  )
p_ssm <- ggpubr::ggarrange(
  p_ssm_a,
  p_ssm_b,
  ncol = 2,
  widths = c(1, 1.3),
  labels = c("a", "b")
)
ggsave(
  file.path(new_dir, "fig_ssm.png"),
  p_ssm,
  width = 250,
  height = 110,
  units = "mm",
  dpi = 200,
  bg = "white"
)

fig_tilt <- arm_frame(summary_tilt, "median_abs_rho")
write_csv(fig_tilt, file.path(new_dir, "fig_tilt_source.csv"))
p_tilt <- plot_arm_summary(
  fig_tilt,
  "Median |Spearman rho| over channel pairs\n(events negative for both)"
)
ggsave(
  file.path(new_dir, "fig_tilt.png"),
  p_tilt,
  width = 140,
  height = 90,
  units = "mm",
  dpi = 200,
  bg = "white"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(new_dir, "sessionInfo_30.txt")
)
