# Purpose: [34d: G13 variant of 34c, one fixed cut of arcsinh 1.5 for negatives and positives instead of per-arm mindensity cuts] rescore negative spread and the guardrail metrics on the common SpectroFlo gates (G1a), with a symmetric negative definition
#          (user decision, 2026-11-10): negatives = cells negative in EVERY arm, each by its own
#          mindensity cut, so every arm is scored on the same cells. 30 used the reference arm's cut,
#          which trims only the reference's upper tail (fable review). 30 is left unchanged.
#          "Every arm" means the 11 figure arms only (ARMS below). The sweep arms of 30 are not scored.
#          Positives (stain index, SSM source) = cells positive in every arm, same symmetry.
#          A donor x channel is excluded from spread, stain index, hypernegatives and SSM when ANY of
#          the 11 arms has a pinned cut there (within PIN_TOL of 0 or 3 arcsinh). Tilt and the
#          TNF-CD8 streak keep all channels (they need no valley), as in 30.
#          Metrics per arm, donor, channel (linear scale; cuts are on arcsinh):
#          - Nguyen robust SD (q84 - q50) primary; 2.5-97.5/3.92 spread as check
#          - stain index = (median pos - median neg) / (2 x Nguyen rSD)
#          - hypernegative % = shared negatives below -3 x the channel's median-over-arms Nguyen rSD
#            (30 used the reference arm's rSD; this is the symmetric replacement)
#          - TNF-CD8 streak: Spearman rho of BV510-A vs SBUV445-A among cells negative for both in every arm
#          - tilt: median |rho| over all channel pairs, cells negative for both in every arm
#          - quadrature SSM as in 30, symmetric populations
#          Every metric is NA when its population has fewer than MIN_N events (stain index and rSD too).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> \
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/34d_rank_symmetric_common_cut15.R [<new_run_dir>]
#          from the repo root, after batch 0 has written SRC_DIR/README.md. <new_run_dir> overrides the registry
#          output run_out("common-cut15-34d"); it must not exist.
#          Copy of 34_rank_symmetric.R: source dir is the common-gate run (SRC_DIR), figures are redrawn
#          in the Nature style (donors as points, median bar).
# Inputs:  run_in("common-final")/manifest.csv, shared_events.rds  [batch 0, common gates]
#          SRC_DIR/mindensity_cuts.csv  cuts for all 11 arms, all samples, 32 channels [batch 0, 02_metrics.R]
#          experiments/benchmarks/05_bmv-unmix-compare/datasets/  config (streak pair), run registry, sample table
#          cut_of(arm, sample, channel, "pos" and "neg")  fixed-cut channels [config pos_rule and neg_rule]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R, report_figures.R  arm order and labels
# Outputs: <new_run_dir>/shared_counts.csv   per donor, channel: pinned_any, n_neg_shared, n_pos_shared
#          <new_run_dir>/rank_channel.csv    per arm, donor, channel: n_neg, n_pos, rsd_nguyen, spread_q,
#                                            stain_index, pct_hyperneg
#          <new_run_dir>/rank_ssm.csv        per arm, donor, source, target: ss and population sizes
#          <new_run_dir>/rank_tilt.csv       per arm, donor, channel pair: n_both_neg, rho
#          <new_run_dir>/rank_streak.csv     per arm, donor: TNF-CD8 rho
#          <new_run_dir>/fig_{spread,stain_index,hyperneg,streak,tilt,ssm}_sym.{pdf,png,csv} and ..._sym_source.csv
#          <new_run_dir>/sessionInfo_34d.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), run_out(), cut_of()
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/report_figures.R") # theme, save_panel, builders, ARM_LABELS, VERSION_ROWS

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NGUYEN_P <- c(0.5, 0.84) # Nguyen 2013 robust SD = q84 - q50 (user, 2026-11-10)
QUANTILES <- c(0.025, 0.975) # outer quantiles, check column (user, 2026-11-09)
Z_WIDTH <- 3.92 # width of the 2.5-97.5 interval in SD units under normality
HYPERNEG_K <- 3 # hypernegative: below -K x median-over-arms Nguyen rSD (K from user grill Q19; rSD choice 2026-11-10)
PIN_TOL <- 0.02 # a cut within this of the window edge counts as pinned (as in 30)
CUT_WINDOW <- c(0, 3) # mindensity window used for the cuts in 02 (user, 2026-11-09)
MIN_N <- 100 # every metric is NA below this many events in its population (user, 2026-11-10)
N_CORES <- 5 # one worker per donor
# TNF-slot x CD8 channel pair from the stripe bisect (2026-11-10); the dataset config owns the pair
STREAK <- cfg("targets")$streak_pair
if (length(STREAK) != 2) {
  stop(
    "34d_rank_symmetric_common_cut15.R: config targets$streak_pair (two channel names) is missing for dataset ",
    DS$dataset,
    call. = FALSE
  )
}
ARMS <- unlist(VERSION_ROWS, use.names = FALSE) # the 11 figure arms, in version-group order

# Batch 0 writes the common-gate shared events and mindensity cuts here (README.md names the files).
SRC_DIR <- run_in("common-final")
new_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(new_dir)) {
  new_dir <- run_out("common-cut15-34d")
}
src_dir <- SRC_DIR
stopifnot("new run dir must not exist" = !dir.exists(new_dir))
dir.create(new_dir, recursive = TRUE)

# Manifest, shared events, cuts for the 11 arms ---------------------------------------------
manifest <- read_csv(file.path(src_dir, "manifest.csv"), show_col_types = FALSE)
stopifnot(
  "manifest arms equal the 11 figure arms" = setequal(manifest$arm, ARMS),
  "manifest samples equal the sample table" = setequal(
    manifest$sample,
    samples_tbl()$sample
  ),
  "one file per sample and arm" = all(
    count(manifest, arm)$n == nrow(samples_tbl())
  )
)
shared <- readRDS(file.path(src_dir, "shared_events.rds"))
# G13 (user, 2026-11-11): negatives are cells below the dataset's negative cut in every arm, for every
# channel; positives are cells above its positive cut (3D: both arcsinh 1.5). The per-arm mindensity cuts
# are replaced by cut_of(), so no sample x channel is pinned.
cuts <- read_csv(
  file.path(src_dir, "mindensity_cuts.csv"),
  show_col_types = FALSE
) |>
  mutate(
    cut_neg = pmap_dbl(list(arm, sample, channel), \(a, s, ch) {
      cut_of(a, s, ch, "neg")
    }),
    cut_pos = pmap_dbl(list(arm, sample, channel), \(a, s, ch) {
      cut_of(a, s, ch, "pos")
    }),
    pinned = FALSE
  )
stopifnot(
  "cuts exist for all 11 arms" = setequal(cuts$arm, ARMS),
  "no missing cut" = !anyNA(c(cuts$cut_neg, cuts$cut_pos))
)
channels <- unique(cuts$channel)

# Nguyen robust SD, linear scale; NA below MIN_N ------------------------------------------
rsd_min <- function(x) {
  if (length(x) < MIN_N) {
    return(NA_real_)
  }
  q <- quantile(x, NGUYEN_P, names = FALSE)
  q[[2]] - q[[1]]
}
median_min <- function(x) if (length(x) < MIN_N) NA_real_ else median(x)

# One donor: shared populations, then every metric for every arm ----------------------------
process_sample <- function(s) {
  lin <- manifest |>
    filter(sample == s) |>
    (\(m) set_names(m$file, m$arm))() |>
    map(\(f) {
      flowCore::exprs(sw_io_read_fcs(f))[shared[[s]], channels, drop = FALSE]
    })
  lin <- lin[ARMS]
  stopifnot(
    "same event count in every arm" = n_distinct(map_int(lin, nrow)) == 1
  )

  # Negative and positive cut and pinned flag per arm, aligned to channels by name
  cut_neg_arm <- map(ARMS, \(a) {
    filter(cuts, sample == s, arm == a) |>
      pull(cut_neg, name = channel) |>
      (\(v) v[channels])()
  }) |>
    set_names(ARMS)
  cut_pos_arm <- map(ARMS, \(a) {
    filter(cuts, sample == s, arm == a) |>
      pull(cut_pos, name = channel) |>
      (\(v) v[channels])()
  }) |>
    set_names(ARMS)
  pinned_any <- cuts |>
    filter(sample == s) |>
    summarise(p = any(pinned), .by = channel) |>
    pull(p, name = channel) |>
    (\(v) v[channels])()

  # A cell is negative (positive) in every arm: below (above) each arm's own cut
  below <- map(ARMS, \(a) {
    sweep(asinh(lin[[a]] / COFACTOR), 2, cut_neg_arm[[a]], "<")
  })
  above <- map(ARMS, \(a) {
    sweep(asinh(lin[[a]] / COFACTOR), 2, cut_pos_arm[[a]], ">")
  })
  neg_all <- reduce(below, `&`)
  pos_all <- reduce(above, `&`)
  # Pinned donor x channel: empty population, so every valley-based metric is NA via MIN_N
  neg_use <- sweep(neg_all, 2, !pinned_any, "&")
  pos_use <- sweep(pos_all, 2, !pinned_any, "&")
  usable <- channels[!pinned_any]

  counts <- tibble(
    sample = s,
    channel = channels,
    pinned_any = unname(pinned_any),
    n_neg_shared = colSums(neg_all),
    n_pos_shared = colSums(pos_all)
  )

  # Median over arms of the Nguyen rSD per channel: the hypernegative scale
  rsd_med <- map(lin, \(m) {
    map_dbl(set_names(channels), \(ch) rsd_min(m[neg_use[, ch], ch]))
  }) |>
    reduce(rbind) |>
    apply(2, median)

  # Per channel: rSD, stain index, hypernegatives
  per_channel <- imap(lin, \(m, a) {
    map(channels, \(ch) {
      x <- m[neg_use[, ch], ch]
      xp <- m[pos_use[, ch], ch]
      q <- quantile(x, QUANTILES, names = FALSE)
      rsd_ng <- rsd_min(x)
      tibble(
        arm = a,
        sample = s,
        channel = ch,
        n_neg = length(x),
        n_pos = length(xp),
        rsd_nguyen = rsd_ng,
        spread_q = if (length(x) < MIN_N) {
          NA_real_
        } else {
          (q[[2]] - q[[1]]) / Z_WIDTH
        },
        stain_index = (median_min(xp) - median_min(x)) / (2 * rsd_ng),
        pct_hyperneg = if (length(x) < MIN_N) {
          NA_real_
        } else {
          100 * mean(x < -HYPERNEG_K * rsd_med[[ch]])
        }
      )
    }) |>
      list_rbind()
  }) |>
    list_rbind()

  # Quadrature spillover-spread, i and j both usable, i != j
  ssm <- imap(lin, \(m, a) {
    map(usable, \(i) {
      den <- median_min(m[pos_use[, i], i]) - median_min(m[neg_use[, i], i])
      map(setdiff(usable, i), \(j) {
        x_ipjn <- m[pos_use[, i] & neg_use[, j], j]
        x_injn <- m[neg_use[, i] & neg_use[, j], j]
        tibble(
          arm = a,
          sample = s,
          source = i,
          target = j,
          n_ip = sum(pos_use[, i]),
          n_in = sum(neg_use[, i]),
          n_ipjn = length(x_ipjn),
          n_injn = length(x_injn),
          ss = if (is.na(den) || den <= 0) {
            NA_real_
          } else {
            sqrt(max(0, rsd_min(x_ipjn)^2 - rsd_min(x_injn)^2)) / sqrt(den)
          }
        )
      }) |>
        list_rbind()
    }) |>
      list_rbind()
  }) |>
    list_rbind()

  # Tilt: Spearman rho for every pair among cells negative for both channels in every arm
  tilt <- imap(lin, \(m, a) {
    combn(channels, 2, simplify = FALSE) |>
      map(\(p) {
        both <- neg_all[, p[[1]]] & neg_all[, p[[2]]]
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
    filter(channel_a %in% STREAK, channel_b %in% STREAK) |>
    transmute(arm, sample, n_both_neg, rho_tnf_cd8 = rho)
  list(
    counts = counts,
    per_channel = per_channel,
    ssm = ssm,
    tilt = tilt,
    streak = streak
  )
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
shared_counts <- results |> map("counts") |> list_rbind()
rank_channel <- results |> map("per_channel") |> list_rbind()
rank_ssm <- results |> map("ssm") |> list_rbind()
rank_tilt <- results |> map("tilt") |> list_rbind()
rank_streak <- results |> map("streak") |> list_rbind()

# Summaries per arm and donor ----------------------------------------------------------------
sum_channel <- rank_channel |>
  summarise(
    across(c(rsd_nguyen, spread_q, stain_index, pct_hyperneg), \(v) {
      median(v, na.rm = TRUE)
    }),
    .by = c(arm, sample)
  )
sum_ssm <- rank_ssm |>
  summarise(median_ss = median(ss, na.rm = TRUE), .by = c(arm, sample))
sum_tilt <- rank_tilt |>
  summarise(
    median_abs_rho = median(abs(rho), na.rm = TRUE),
    .by = c(arm, sample)
  )
summaries <- list(
  spread = sum_channel |> transmute(arm, sample, value = rsd_nguyen),
  stain_index = sum_channel |> transmute(arm, sample, value = stain_index),
  hyperneg = sum_channel |> transmute(arm, sample, value = pct_hyperneg),
  streak = rank_streak |> transmute(arm, sample, value = rho_tnf_cd8),
  tilt = sum_tilt |> transmute(arm, sample, value = median_abs_rho),
  ssm = sum_ssm |> transmute(arm, sample, value = median_ss)
)

# Never fail silently: no all-NA metric column, every sample x 11 arms in every summary ----------
iwalk(
  list(
    rank_channel = rank_channel[c(
      "rsd_nguyen",
      "spread_q",
      "stain_index",
      "pct_hyperneg"
    )],
    rank_ssm = rank_ssm["ss"],
    rank_tilt = rank_tilt["rho"],
    rank_streak = rank_streak["rho_tnf_cd8"]
  ),
  \(d, tbl) {
    iwalk(d, \(v, col) {
      if (all(is.na(v))) {
        stop(sprintf("%s$%s is all NA", tbl, col), call. = FALSE)
      }
    })
  }
)
iwalk(summaries, \(d, nm) {
  if (all(is.na(d$value))) {
    stop("summary ", nm, " is all NA", call. = FALSE)
  }
  stopifnot(
    "one row per arm and sample" = nrow(d) == length(ARMS) * length(samples)
  )
})
n_excl <- sum(shared_counts$pinned_any)
message(
  "Excluded donor x channel cells (pinned in any arm): ",
  n_excl,
  " of ",
  nrow(shared_counts)
)
message(
  "Metric cells NA: ",
  paste(
    imap_chr(summaries, \(d, nm) paste0(nm, "=", sum(is.na(d$value)))),
    collapse = ", "
  )
)

write_csv(shared_counts, file.path(new_dir, "shared_counts.csv"))
write_csv(rank_channel, file.path(new_dir, "rank_channel.csv"))
write_csv(rank_ssm, file.path(new_dir, "rank_ssm.csv"))
write_csv(rank_tilt, file.path(new_dir, "rank_tilt.csv"))
write_csv(rank_streak, file.path(new_dir, "rank_streak.csv"))

# Figures: arm down the y axis in version-group order, donors as points over a median bar -----
metric_titles <- c(
  spread = "fig_spread_sym",
  stain_index = "fig_stain_index_sym",
  hyperneg = "fig_hyperneg_sym",
  streak = "fig_streak_sym",
  tilt = "fig_tilt_sym",
  ssm = "fig_ssm_sym"
)
iwalk(metric_titles, \(fig_name, name) {
  src <- file.path(new_dir, sprintf("%s_source.csv", fig_name))
  summaries[[name]] |>
    mutate(arm = factor(arm, levels = ARMS)) |>
    arrange(arm, sample) |>
    mutate(arm = as.character(arm)) |>
    write_csv(src)
  redraw_panel(src, "rank", fig_name, new_dir, fig_name)
})

writeLines(
  capture.output(sessionInfo()),
  file.path(new_dir, "sessionInfo_34d.txt")
)
