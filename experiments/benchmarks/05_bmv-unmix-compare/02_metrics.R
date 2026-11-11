# Purpose: spread, spillover-spread matrix, lower-tail artefacts, negative-population leaning and
#          stain index for the four unmixing arms, on the shared live singlet events only.
# Run:     sbatch experiments/slurm/run_script.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/02_metrics.R \
#            experiments/runs/2026-11-09_bmv-unmix-compare
# Inputs:  <run_dir>/manifest.csv        arm, dir, viability, file, sample [01_manifest-gate.R]
#          <run_dir>/shared_events.rds   per sample: logical vector over all events, TRUE = event
#                                        is nonDebris + singlet + live in every arm [01_manifest-gate.R]
# Outputs: <run_dir>/quality.csv          sw_unmix_quality() per arm and sample (no CV flag column)
#          <run_dir>/mindensity_cuts.csv  mindensity cut per arm, sample, channel (NA = no valley)
#          <run_dir>/negative_spread.csv  P5 median and rSD of the shared negatives per channel
#          <run_dir>/artefacts.csv        P7 fraction of negatives below median - 3 rSD
#          <run_dir>/spread_matrix.csv    P6 SS(i -> j), i = separation marker, j = marker channel
#          <run_dir>/neg_correlation.csv  P8 Spearman rho per channel pair, events negative for both
#          <run_dir>/stain_index.csv      P9 stain index per separation marker
#          <run_dir>/sessionInfo_02.txt   (not sessionInfo.txt: step 1 already wrote that name)

# Bioconductor first, tidyverse parts last, so dplyr::filter() and select() are not masked.
# tidyverse and here are not in the pixi env; run from the repo root.
suppressPackageStartupMessages({
  library(flowCore)
  library(openCyto)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

# Parameters ------------------------------------------------------------------------------
COFACTOR <- 6000 # arcsinh cofactor for spectral flow (CLAUDE.md, confirmed by user)
CUT_RANGE <- c(0, 3) # mindensity search window on arcsinh (user, 2026-11-09): default put cuts in the tails
MIN_SHARED_POS <- 100 # P9: fewest shared positives before the reference cut is used (user, 2026-11-09)
LOWER_TAIL_SD <- 3 # P7: lower tail is below median - 3 rSD
REFERENCE_ARM <- "as188_joint" # its $PnS descriptions define the marker labels
# Separation markers (README, user 2026-11-09) and their detectors; verified against the
# as188_joint $PnS descriptions in process_sample().
SEP_MARKERS <- tribble(
  ~marker  , ~channel           ,
  "CD3"    , "SBUV510-A"        ,
  "CD4"    , "Spark Blue 574-A" ,
  "CD8"    , "SBUV445-A"        ,
  "CD45RA" , "Spark UV 387-A"   ,
  "CCR7"   , "PE-Fire 810-A"    ,
  "CD127"  , "APC-Fire 810-A"   ,
  "CD25"   , "SBUV605-A"        ,
  "CD56"   , "RB744-A"          ,
  "TCRgd"  , "BUV395-A"         ,
  "CD161"  , "BUV563-A"         ,
  "TNF"    , "BV510-A"          ,
  "IFNg"   , "BV570-A"          ,
  "IL-2"   , "BV605-A"          ,
  "IL-10"  , "BB700-A"
)

# Run dir ---------------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot(
  "step 2 outputs must not exist" = !file.exists(file.path(
    run_dir,
    "quality.csv"
  ))
)
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
shared <- readRDS(file.path(run_dir, "shared_events.rds"))
samples <- names(shared)
arms <- unique(manifest$arm)

# Robust SD (P5): 1.4826 x MAD on the linear scale. An empty set gives NA.
rsd <- function(x) 1.4826 * mad(x, constant = 1)

# One sample ------------------------------------------------------------------------------
# Reads the four arms, keeps shared events, and returns one tibble per output file.
process_sample <- function(s) {
  mask <- shared[[s]]

  # Read each arm; linear values are the raw values, because sinh(asinh(x / 6000)) * 6000 = x.
  arm_data <- arms |>
    set_names() |>
    map(\(a) {
      row <- manifest |> filter(arm == a, sample == s)
      ff <- sw_io_read_fcs(row$file)
      stopifnot(
        "shared mask length matches event count" = length(mask) == nrow(ff)
      )
      ff <- ff[mask, ]
      # Marker set: fluorescence channels minus AF-A, AF Index and this arm's viability channel.
      chans <- flowCore::colnames(ff)[sw_channel_is_fluor(ff)] |>
        setdiff(c("AF-A", "AF Index", row$viability))
      list(
        lin = exprs(ff)[, chans],
        descriptions = markernames(ff),
        quality = suppressMessages(sw_unmix_quality(ff)) |>
          select(-flagged) |>
          mutate(arm = a, sample = s, .before = 1)
      )
    })

  # Channels are joined by detector name: every arm needs the same marker channels, and the
  # separation detectors must carry the expected label in the reference arm.
  markers <- colnames(arm_data[[REFERENCE_ARM]]$lin)
  ref_desc <- arm_data[[REFERENCE_ARM]]$descriptions
  stopifnot(
    "same marker channels in every arm" = all(map_lgl(arm_data, \(d) {
      setequal(colnames(d$lin), markers)
    })),
    "separation channels carry the expected labels" = identical(
      unname(ref_desc[SEP_MARKERS$channel]),
      SEP_MARKERS$marker
    ),
    "separation channels are marker channels" = all(
      SEP_MARKERS$channel %in% markers
    )
  )
  arm_lin <- arm_data |> map(\(d) d$lin[, markers])
  arm_art <- arm_lin |> map(\(m) asinh(m / COFACTOR))

  # Per arm and channel: mindensity cut on arcsinh. A cut at the data minimum means openCyto
  # found no valley (unimodal channel), so it is NA.
  cuts <- arm_art |>
    map(\(m) {
      fr <- flowFrame(m)
      markers |>
        set_names() |>
        map_dbl(\(ch) {
          cut <- gate_mindensity(fr, channel = ch, gate_range = CUT_RANGE)@min[[
            1
          ]]
          if (cut <= min(m[, ch])) NA_real_ else cut
        })
    })
  cuts_tbl <- cuts |>
    imap(\(v, a) {
      tibble(
        arm = a,
        sample = s,
        channel = markers,
        cut = unname(v),
        found = !is.na(v)
      )
    }) |>
    list_rbind()

  # P4 negatives: below the arm's own cut in every arm. Channels with an NA cut in any arm
  # get an empty negative set in all arms, so all downstream numbers are NA there.
  neg <- arm_art |>
    imap(\(m, a) {
      below <- sweep(m, 2, cuts[[a]][markers], "<")
      below[is.na(below)] <- FALSE
      below
    }) |>
    reduce(`&`)
  # Positives for the separation markers: at or above the arm's own cut (used by P6).
  pos <- arm_art |>
    imap(\(m, a) {
      above <- sweep(
        m[, SEP_MARKERS$channel],
        2,
        cuts[[a]][SEP_MARKERS$channel],
        ">="
      )
      above[is.na(above)] <- FALSE
      above
    })
  pos_shared <- pos |> reduce(`&`)

  # Metrics per arm.
  per_arm <- arms |>
    set_names() |>
    map(\(a) {
      lin <- arm_lin[[a]]
      above <- pos[[a]]

      # P5 spread and P7 lower tail of the shared negatives, per channel.
      neg_stats <- markers |>
        map(\(ch) {
          x <- lin[neg[, ch], ch]
          med <- median(x)
          sd_neg <- rsd(x)
          tibble(
            arm = a,
            sample = s,
            channel = ch,
            n_neg = length(x),
            median_neg = med,
            rsd_neg = sd_neg,
            lower_tail_frac = mean(x < med - LOWER_TAIL_SD * sd_neg)
          )
        }) |>
        list_rbind()

      # P6 spread matrix. Medians of i are taken over the whole i+ (arm cut) and i- (shared
      # negatives) populations; rSD of j is taken in the i+ j- and i- j- subsets.
      spread <- SEP_MARKERS |>
        pmap(\(marker, channel) {
          i_pos <- above[, channel]
          i_neg <- neg[, channel]
          med_pos <- median(lin[i_pos, channel])
          med_neg <- median(lin[i_neg, channel])
          gap <- med_pos - med_neg
          setdiff(markers, channel) |>
            map(\(ch) {
              set_pos <- i_pos & neg[, ch]
              set_neg <- i_neg & neg[, ch]
              rsd_pos <- rsd(lin[set_pos, ch])
              rsd_neg <- rsd(lin[set_neg, ch])
              tibble(
                arm = a,
                sample = s,
                marker_i = marker,
                channel_i = channel,
                channel_j = ch,
                n_pos_jneg = sum(set_pos),
                n_neg_jneg = sum(set_neg),
                median_i_pos = med_pos,
                median_i_neg = med_neg,
                rsd_j_pos = rsd_pos,
                rsd_j_neg = rsd_neg,
                ss = if (is.na(gap) || gap <= 0) {
                  NA_real_
                } else {
                  (rsd_pos - rsd_neg) / sqrt(gap)
                }
              )
            }) |>
            list_rbind()
        }) |>
        list_rbind()

      # P8 Spearman rho between every pair of marker channels, events negative for both.
      corr <- combn(markers, 2, simplify = FALSE) |>
        map(\(p) {
          both <- neg[, p[[1]]] & neg[, p[[2]]]
          n_both <- sum(both)
          tibble(
            arm = a,
            sample = s,
            channel_a = p[[1]],
            channel_b = p[[2]],
            n_both = n_both,
            rho = if (n_both < 3) {
              NA_real_
            } else {
              cor(lin[both, p[[1]]], lin[both, p[[2]]], method = "spearman")
            }
          )
        }) |>
        list_rbind()

      # P9 stain index (user, 2026-11-09): positives are events above the cut in ALL arms, so every arm
      # scores the same cells; negatives are the shared negative set. Where fewer than
      # MIN_SHARED_POS shared positives exist, positives fall back to the reference arm's cut.
      stain <- SEP_MARKERS |>
        pmap(\(marker, channel) {
          shared_pos <- pos_shared[, channel]
          use_shared <- sum(shared_pos) >= MIN_SHARED_POS
          x_pos <- lin[
            if (use_shared) shared_pos else pos[[REFERENCE_ARM]][, channel],
            channel
          ]
          x_neg <- lin[neg[, channel], channel]
          tibble(
            arm = a,
            sample = s,
            marker = marker,
            channel = channel,
            pos_rule = if (use_shared) "shared_all_arms" else "reference_cut",
            n_pos = length(x_pos),
            n_neg = length(x_neg),
            median_pos = median(x_pos),
            median_neg = median(x_neg),
            rsd_neg = rsd(x_neg),
            stain_index = (median(x_pos) - median(x_neg)) / (2 * rsd(x_neg))
          )
        }) |>
        list_rbind()

      list(neg_stats = neg_stats, spread = spread, corr = corr, stain = stain)
    })

  list(
    quality = arm_data |> map("quality") |> list_rbind(),
    cuts = cuts_tbl,
    neg_stats = per_arm |> map("neg_stats") |> list_rbind(),
    spread = per_arm |> map("spread") |> list_rbind(),
    corr = per_arm |> map("corr") |> list_rbind(),
    stain = per_arm |> map("stain") |> list_rbind()
  )
}

# Run all samples and write ---------------------------------------------------------------
results <- samples |>
  set_names() |>
  map(process_sample)

results |>
  map("quality") |>
  list_rbind() |>
  write_csv(file.path(run_dir, "quality.csv"))
results |>
  map("cuts") |>
  list_rbind() |>
  write_csv(file.path(run_dir, "mindensity_cuts.csv"))
results |>
  map("neg_stats") |>
  list_rbind() |>
  select(arm, sample, channel, n_neg, median_neg, rsd_neg) |>
  write_csv(file.path(run_dir, "negative_spread.csv"))
results |>
  map("neg_stats") |>
  list_rbind() |>
  select(arm, sample, channel, n_neg, lower_tail_frac) |>
  write_csv(file.path(run_dir, "artefacts.csv"))
results |>
  map("spread") |>
  list_rbind() |>
  write_csv(file.path(run_dir, "spread_matrix.csv"))
results |>
  map("corr") |>
  list_rbind() |>
  write_csv(file.path(run_dir, "neg_correlation.csv"))
results |>
  map("stain") |>
  list_rbind() |>
  write_csv(file.path(run_dir, "stain_index.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_02.txt")
)
