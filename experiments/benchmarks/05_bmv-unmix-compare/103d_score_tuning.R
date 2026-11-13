# Purpose: score the AS 1.8.8 tuning arms on the smoke donor (stage A, metrics from the unmixed files). Thresholds are
#          data-driven per arm and channel from the arm's own unmixed unstained control (103c): floor = q0.5, edge = q99.5.
#          Per arm, sample and channel, in live cells: % below the floor; spread_q of negatives (below the edge) =
#          (q97.5 - q2.5) / 3.92; stain index_q = (median pos - median neg) / (2 spread_q), positives above the edge; pile-up
#          index = MAD rSD / spread_q. Guardrails: % TNF+ in CD4 and CD8 T cells in SEB against the donor uRBC cut (relative
#          change vs the baseline arm, limit 10%); TNF-CD8 Spearman rho among CD3+ cells below both edges. Conditioning:
#          condition number of the dye + AF spectra matrix. Held-out control error and spike-in recovery are stage B.
#          SW_TUNE_SCOPE=smoke scores the smoke donor (tuning-score); SW_TUNE_SCOPE=all scores all 30 files (confirm-score).
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> SW_TUNE_SCOPE=<smoke|all> sbatch experiments/slurm/run_common_gs.sbatch 103d_score_tuning.R <tuning root>
# Inputs:  <tuning root>/<arm>/*.fcs, spectra.csv, af_spectra.csv, unstained/thresholds.csv   [103, 103c]
#          run_in("pop-gates")/labels.csv.gz   pooled AS 1.8.8 labels by raw row, with FSC_A and SSC_A for pairing [99a]
#          run_in("common-cuts")/cuts.csv       donor uRBC q99.5 TNF cut [97e]
# Outputs: <out_dir> = run_out("tuning-score")
#          <out_dir>/channel_metrics.csv    arm, sample, channel, marker, floor, edge, pct_below_floor, spread_q, stain_index_q, pileup
#          <out_dir>/guardrails.csv         arm, sample, population, % TNF+, relative change vs baseline, pass; TNF-CD8 rho
#          <out_dir>/arm_summary.csv        arm: medians over samples and channels, condition number, guardrail pass
#          <out_dir>/figures/arm_metrics.png      metric per arm, one point per sample (summary, not events)
#          <out_dir>/figures/arm_metrics_paired.png  same metrics, one line per sample across arms
#          <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters -----------------------------------------------------------------------------
ROOT <- commandArgs(trailingOnly = TRUE)[[1]]
SCOPE <- Sys.getenv("SW_TUNE_SCOPE") # smoke = smoke donor; all = every sample (confirmation)
stopifnot("SW_TUNE_SCOPE is smoke or all" = SCOPE %in% c("smoke", "all"))
BASELINE <- "a188_rFALSE_p1_t0.01" # zip settings
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
TNF <- "BV510-A"
CD8 <- cfg("channels")$CD8
TNF_TOLERANCE <- 0.10 # SEB TNF+ may not drop by more than 10% relative to baseline (user 2026-11-14)
MIN_POS <- 100 # stain index needs this many positives, as 3D (02_metrics.R)
Z_WIDTH <- 3.92 # width of the 2.5-97.5 interval in SD units under normality (15_quantile_spread.R)
OUT_DIR <- run_out(c(smoke = "tuning-score", all = "confirm-score")[[SCOPE]])
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Inputs --------------------------------------------------------------------------------------
smoke <- samples_tbl() |>
  dplyr::filter(SCOPE == "all" | .data$donor == cfg("smoke_donor"))
labels <- read_csv(file.path(run_in("pop-gates"), "labels.csv.gz"), show_col_types = FALSE) |>
  dplyr::filter(sample %in% smoke$sample)
tnf_cut <- read_csv(file.path(run_in("common-cuts"), "cuts.csv"), show_col_types = FALSE) |>
  dplyr::filter(arm == "as188_joint", sample %in% smoke$sample, channel == TNF) |>
  select(sample, pos_cut)
ctrl <- read_csv(data_path("autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"), show_col_types = FALSE) |>
  dplyr::filter(!fluorophore %in% c("AF", "Negative")) |>
  distinct(fluorophore, marker)
marker_of <- set_names(ctrl$marker, paste0(ctrl$fluorophore, "-A"))
arms <- basename(list.dirs(ROOT, recursive = FALSE))
arms <- arms[file.exists(file.path(ROOT, arms, "unstained", "thresholds.csv"))]
stopifnot("baseline arm scored" = BASELINE %in% arms)

# Per arm and sample ------------------------------------------------------------------------------
score_one <- \(arm, s) {
  st <- smoke$stem[smoke$sample == s]
  f <- file.path(ROOT, arm, paste0(st, " AutoSpectral.fcs"))
  m <- exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
  lb <- labels[labels$sample == s, ]
  stopifnot(
    "unmixed rows pair with the label rows" = nrow(m) == nrow(lb) &&
      isTRUE(all.equal(unname(m[, "FSC-A"]), lb$FSC_A, tolerance = 1e-6)) &&
      isTRUE(all.equal(unname(m[, "SSC-A"]), lb$SSC_A, tolerance = 1e-6))
  )
  th <- read_csv(file.path(ROOT, arm, "unstained", "thresholds.csv"), show_col_types = FALSE)
  live <- asinh(m[lb$live, th$channel, drop = FALSE] / COFACTOR)
  ch_tbl <- pmap(th, \(channel, n, floor, edge) {
    x <- live[, channel]
    neg <- x[x < edge]
    pos <- x[x >= edge]
    spread <- unname(diff(quantile(neg, c(0.025, 0.975)))) / Z_WIDTH
    tibble(arm, sample = s, channel, floor, edge, pct_below_floor = 100 * mean(x < floor), spread_q = spread,
           stain_index_q = if (length(pos) >= MIN_POS) (median(pos) - median(neg)) / (2 * spread) else NA_real_,
           pileup = mad(neg) / spread)
  }) |> list_rbind()
  cut <- tnf_cut$pos_cut[tnf_cut$sample == s]
  tnf <- asinh(m[, TNF] / COFACTOR)
  cd8 <- asinh(m[, CD8] / COFACTOR)
  e <- set_names(th$edge, th$channel)
  dn <- lb$CD3pos & tnf < e[[TNF]] & cd8 < e[[CD8]]
  gr <- tibble(arm, sample = s, population = c("CD4T", "CD8T"),
               pct_tnf = c(100 * mean(tnf[lb$CD4T] >= cut), 100 * mean(tnf[lb$CD8T] >= cut)),
               tnf_cd8_rho = suppressWarnings(cor(tnf[dn], cd8[dn], method = "spearman")))
  list(ch = ch_tbl, gr = gr)
}
grid <- expand.grid(arm = arms, s = smoke$sample, stringsAsFactors = FALSE)
res <- pmap(grid, \(arm, s) score_one(arm, s))
ch <- map(res, "ch") |> list_rbind() |> mutate(marker = unname(marker_of[channel]))
gr <- map(res, "gr") |>
  list_rbind() |>
  left_join(select(smoke, sample, condition), by = join_by(sample), relationship = "many-to-one")
base <- gr |> dplyr::filter(arm == BASELINE) |> select(sample, population, base_tnf = pct_tnf)
gr <- gr |>
  left_join(base, by = join_by(sample, population), relationship = "many-to-one") |>
  mutate(rel_change = pct_tnf / base_tnf - 1, pass = condition != "SEB" | rel_change >= -TNF_TOLERANCE)
write_csv(ch, file.path(OUT_DIR, "channel_metrics.csv"))
write_csv(gr, file.path(OUT_DIR, "guardrails.csv"))

# Conditioning and arm summary -----------------------------------------------------------------------
cond_num <- map(arms, \(a) {
  sp <- as.matrix(read.csv(file.path(ROOT, a, "spectra.csv"), row.names = 1, check.names = FALSE))
  af <- as.matrix(read.csv(file.path(ROOT, a, "af_spectra.csv"), row.names = 1, check.names = FALSE))
  tibble(arm = a, n_af = nrow(af), kappa_dyes = kappa(t(sp), exact = TRUE), kappa_dyes_af = kappa(t(rbind(sp, af)), exact = TRUE))
}) |> list_rbind()
summ <- ch |>
  summarise(tnf_below_floor = median(pct_below_floor[channel == TNF]), pct_below_floor = median(pct_below_floor),
            spread_q = median(spread_q), stain_index_q = median(stain_index_q, na.rm = TRUE), pileup = median(pileup), .by = arm) |>
  left_join(summarise(gr, guard_pass = all(pass), n_guard_fail = sum(!pass),
                      seb_cd4_tnf_change = median(rel_change[condition == "SEB" & population == "CD4T"]),
                      seb_cd8_tnf_change = median(rel_change[condition == "SEB" & population == "CD8T"]),
                      tnf_cd8_rho = median(tnf_cd8_rho), .by = arm), by = join_by(arm), relationship = "one-to-one") |>
  left_join(cond_num, by = join_by(arm), relationship = "one-to-one")
write_csv(summ, file.path(OUT_DIR, "arm_summary.csv"))

# Figure: per-arm metrics, one point per sample ------------------------------------------------------------
per_sample <- ch |>
  summarise(`% live below floor (all channels, median)` = median(pct_below_floor),
            `% live below floor (TNF)` = pct_below_floor[channel == TNF],
            `spread_q (median over channels)` = median(spread_q),
            `stain index_q (median)` = median(stain_index_q, na.rm = TRUE), .by = c(arm, sample)) |>
  left_join(gr |> dplyr::filter(population == "CD4T") |> select(arm, sample, `% CD4 TNF+` = pct_tnf, `TNF-CD8 rho` = tnf_cd8_rho),
            by = join_by(arm, sample), relationship = "one-to-one") |>
  pivot_longer(-c(arm, sample), names_to = "metric") |>
  left_join(select(smoke, sample, condition), by = join_by(sample), relationship = "many-to-one") |>
  mutate(arm = factor(arm, levels = sort(arms)))
p <- ggplot(per_sample, aes(value, arm)) +
  geom_point(aes(colour = condition), size = 1.2, alpha = 0.7) +
  facet_wrap(vars(metric), scales = "free_x", nrow = 1) +
  labs(x = NULL, y = NULL, colour = NULL, caption = sprintf("%d samples (%s); thresholds from each arm's unstained control.", n_distinct(per_sample$sample), SCOPE)) +
  theme(legend.position = "right")
ggsave(file.path(OUT_DIR, "figures", "arm_metrics.png"), p, width = 290, height = 110, units = "mm", dpi = 150, bg = "white")
p2 <- ggplot(per_sample, aes(arm, value, group = sample)) +
  geom_line(aes(colour = condition), linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(colour = condition), size = 0.9) +
  facet_wrap(vars(metric), scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, colour = NULL, caption = "One line per sample across arms.") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")
ggsave(file.path(OUT_DIR, "figures", "arm_metrics_paired.png"), p2, width = 290, height = 100, units = "mm", dpi = 150, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
