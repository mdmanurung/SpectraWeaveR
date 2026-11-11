# Purpose: rescore negative spread and the guardrail metrics with a symmetric negative definition
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
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/34_rank_symmetric.R <src_run_dir> <new_run_dir>
#          from the repo root. <new_run_dir> must not exist.
# Inputs:  <src_run_dir>/manifest.csv, shared_events.rds [01_manifest-gate.R]
#          <src_run_dir>/mindensity_cuts.csv  cuts for all 11 arms, 5 donors, 32 channels [02_metrics.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  arm order and labels
# Outputs: <new_run_dir>/shared_counts.csv   per donor, channel: pinned_any, n_neg_shared, n_pos_shared
#          <new_run_dir>/rank_channel.csv    per arm, donor, channel: n_neg, n_pos, rsd_nguyen, spread_q,
#                                            stain_index, pct_hyperneg
#          <new_run_dir>/rank_ssm.csv        per arm, donor, source, target: ss and population sizes
#          <new_run_dir>/rank_tilt.csv       per arm, donor, channel pair: n_both_neg, rho
#          <new_run_dir>/rank_streak.csv     per arm, donor: TNF-CD8 rho
#          <new_run_dir>/fig_{spread,stain_index,hyperneg,streak,tilt,ssm}_sym.png and ..._source.csv
#          <new_run_dir>/sessionInfo_34.txt

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
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R")

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
STREAK <- c(TNF = "BV510-A", CD8 = "SBUV445-A") # pair from the stripe bisect (2026-11-10)
ARMS <- unlist(VERSION_ROWS, use.names = FALSE) # the 11 figure arms, in version-group order
DONOR_PAL <- c(
  A3 = "#1b9e77",
  A6 = "#d95f02",
  B11 = "#7570b3",
  B2 = "#e7298a",
  D1 = "#66a61e"
)

args <- commandArgs(trailingOnly = TRUE)
src_dir <- args[[1]]
new_dir <- args[[2]]
stopifnot("new run dir must not exist" = !dir.exists(new_dir))
dir.create(new_dir, recursive = TRUE)

# Manifest, shared events, cuts for the 11 arms ---------------------------------------------
manifest <- read_csv(file.path(src_dir, "manifest.csv"), show_col_types = FALSE)
stopifnot(
  "manifest arms equal the 11 figure arms" = setequal(manifest$arm, ARMS),
  "5 files per arm" = all(count(manifest, arm)$n == 5)
)
shared <- readRDS(file.path(src_dir, "shared_events.rds"))
cuts <- read_csv(
  file.path(src_dir, "mindensity_cuts.csv"),
  show_col_types = FALSE
) |>
  mutate(
    pinned = is.na(cut) |
      abs(cut - CUT_WINDOW[[1]]) <= PIN_TOL |
      abs(cut - CUT_WINDOW[[2]]) <= PIN_TOL
  )
stopifnot(
  "cuts exist for all 11 arms" = setequal(cuts$arm, ARMS),
  "no missing cut" = !anyNA(cuts$cut)
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

  # Cut and pinned flag per arm, aligned to channels by name
  cut_arm <- map(ARMS, \(a) {
    filter(cuts, sample == s, arm == a) |>
      pull(cut, name = channel) |>
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
    sweep(asinh(lin[[a]] / COFACTOR), 2, cut_arm[[a]], "<")
  })
  above <- map(ARMS, \(a) {
    sweep(asinh(lin[[a]] / COFACTOR), 2, cut_arm[[a]], ">")
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

# Never fail silently: no all-NA metric column, 5 donors x 11 arms in every summary ----------
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
  stopifnot("55 arm x donor rows" = nrow(d) == length(ARMS) * length(samples))
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

# Figures: arm on x in version-group order, donors paired by lines, boxplot per arm ----------
theme_set(ggpubr::theme_pubr(base_size = 8))
group_end <- cumsum(lengths(VERSION_ROWS))
group_mid <- group_end - lengths(VERSION_ROWS) / 2 + 0.5
group_names <- c("SpectroFlo", "AS 1.0.2", "AS 1.8.8", "crossed")
caption <- ARM_LABELS[ARMS] |>
  str_replace_all("\n", ": ") |>
  (\(l) paste0(ARMS, " = ", l))() |>
  paste(collapse = "; ") |>
  str_wrap(width = 150)

plot_paired <- function(name, ylab) {
  d <- summaries[[name]] |> mutate(arm = factor(arm, levels = ARMS))
  write_csv(d, file.path(new_dir, sprintf("fig_%s_sym_source.csv", name)))
  p <- ggplot(d, aes(arm, value)) +
    geom_boxplot(
      outlier.shape = NA,
      fill = "grey92",
      colour = "grey40",
      width = 0.55
    ) +
    geom_line(
      aes(group = sample, colour = sample),
      linewidth = 0.3,
      alpha = 0.7
    ) +
    geom_point(aes(colour = sample), size = 1.3) +
    geom_vline(
      xintercept = head(group_end, -1) + 0.5,
      linetype = "dashed",
      colour = "grey50"
    ) +
    annotate(
      "text",
      x = group_mid,
      y = Inf,
      label = group_names,
      vjust = 1.4,
      size = 2.6,
      fontface = "bold"
    ) +
    scale_colour_manual(
      values = DONOR_PAL,
      breaks = names(DONOR_PAL),
      limits = names(DONOR_PAL)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = ylab, colour = "Donor", caption = caption) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      plot.caption = element_text(size = 5, hjust = 0)
    )
  ggsave(
    file.path(new_dir, sprintf("fig_%s_sym.png", name)),
    p,
    width = 170,
    height = 105,
    units = "mm",
    dpi = 200,
    bg = "white"
  )
}
plot_paired(
  "spread",
  "Median over channels of Nguyen rSD\n(q84 - q50, cells negative in every arm, linear)"
)
plot_paired(
  "stain_index",
  "Median over channels of stain index\n(Nguyen rSD, shared positives and negatives)"
)
plot_paired(
  "hyperneg",
  "Median over channels of hypernegatives (%)\n(below -3 x median-over-arms Nguyen rSD)"
)
plot_paired(
  "streak",
  "TNF-CD8 Spearman rho\n(BV510-A vs SBUV445-A, double negative in every arm)"
)
plot_paired(
  "tilt",
  "Median |Spearman rho| over channel pairs\n(negative for both in every arm)"
)
plot_paired(
  "ssm",
  "Median quadrature SS over source i, target j\n(symmetric populations)"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(new_dir, "sessionInfo_34.txt")
)
