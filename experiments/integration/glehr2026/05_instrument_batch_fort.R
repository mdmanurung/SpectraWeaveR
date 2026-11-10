# experiments/integration/glehr2026/05_instrument_batch_fort.R
# Purpose: repeat 03_instrument_batch.R on the FORT study (Glehr et al. 2026): 14 samples stained once
#          and acquired on a Fortessa, an LX and a Navios, so the instrument is a known batch with
#          identical biology. Score each marker's batch effect with sw_batch_scores() before CD3
#          relativisation, after the authors' relativisation, and after ours (sw_relativise_fit +
#          sw_relativise_apply, refitted per instrument and sample), and read the instrument gain
#          spread and the control separation from sw_relativise_diagnose().
# Run:     Rscript experiments/integration/glehr2026/05_instrument_batch_fort.R <data_dir> <run_dir>
# Inputs:  <data_dir>/FORT/<instrument>/<sample>/*.fcs   CD3 controls (01 to 10), 12-panel, and the
#          authors' *.relativized-*12-panel.fcs (note the z; BAD files use an s)
# Outputs: <run_dir>/inventory.csv            per unit: events, share of authors' events matched, maps
#          <run_dir>/batch_scores_long.csv    stage, check, marker, batch pair, value
#          <run_dir>/relativise_fit.csv       our m, M, gain, separation per unit and channel
#          <run_dir>/relativise_diagnose.csv  gain checks of our model with instrument as batch
#          <run_dir>/separation_summary.csv   separation and gain per instrument and channel
#          <run_dir>/fig_instrument_eta2.png, figure_source_data_instrument_eta2.csv
#          <run_dir>/fig_separation.png, figure_source_data_separation.csv
#
# Channels: all three instruments use the same FCS parameter names (FITC-A ... KrO-A, FSC-A, FSC-W,
#   SSC-A, TIME), and $PnS holds the same marker-fluor labels, so no channel mapping is needed. The
#   fluor in the CD3 control file name is the channel name without "-A".
#
# Authors' rescale (approved 2026-11-10, "recover linear scatter map"): the authors' relativized files
#   rescale FSC-A, SSC-A, FSC-W and TIME to 0..10000 per file, as a min-max over their gated subset.
#   It is not a fixed map: not raw / $PnR * 1e4 and not full-file min-max (header $PnR of the authors'
#   files is 1048576 for all channels, and the subset min and max are not the raw file's). The map
#   is recovered per sample and instrument: the TIME span a and start t0 from the step structure of
#   the early events (rel diffs are integer multiples of a * the raw tick), then scatter by least
#   squares on events whose time cell holds one raw event, then all events are matched on all four
#   columns within EXACT_TOL = 0.01 (float32 noise is about 0.001). The same map applies to the
#   controls of that sample and instrument; the gate box is mapped back to raw units with it.
#   Navios TIME wraps once at 2^20; the map is affine on the wrapped values, so no unwrapping.
#
# Choices approved by the user on 2026-11-10 for the BAD dataset and reused unchanged:
#   gate box 1-99% of FSC-A x SSC-A of the authors' gated panel events per sample and instrument;
#   asinh(x / 6000) split, seed 1, nstart 10, min_events 200; low 0, high 1e4; batch scores on the
#   linear scale ranked by eta2; 2000 events per file, seed 1; order = NULL.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) stop("run dir exists: ", run_dir, call. = FALSE)
dir.create(run_dir, recursive = TRUE)

INSTRUMENTS <- c("fortessa", "lx", "navios")
GATE_PROBS <- c(0.01, 0.99)
COFACTOR <- 6000
SEED <- 1
NSTART <- 10
MIN_EVENTS <- 200
N_EVENTS <- 2000
LOW <- 0 # paper target for CD3- median
HIGH <- 1e4 # paper target for CD3+ median
MAP_COLS <- c("TIME", "FSC-A", "SSC-A", "FSC-W")
FIT_COLS <- c("FSC-A", "SSC-A", "FSC-W")
EXACT_TOL <- 0.01 # on the 0..10000 scale; float32 noise is about 0.001
PASS_TOL <- 0.02 # tolerance while the map is still being refined
EARLY_NOISE <- 6e-8 # float32 relative precision

read_exprs <- function(path) {
  exprs(read.FCS(path, truncate_max_range = FALSE, transformation = FALSE, emptyValue = FALSE))
}

# Recover the authors' linear rescale for one sample x instrument ------------------------------
# o: raw panel matrix. r: authors' relativised matrix (a subset of o, rescaled per column).

lin_fit <- function(x, y) {
  f <- lm.fit(cbind(1, x), y)
  c(intercept = f$coefficients[[1]], slope = f$coefficients[[2]])
}

# TIME span and offset from the step structure of the authors' early events.
recover_time <- function(ot, rt) {
  su <- sort(unique(ot))
  res <- signif(min(diff(su)), 2)
  a_c <- unname(diff(quantile(rt, c(0.1, 0.9))) / diff(quantile(ot, c(0.1, 0.9))))
  cand <- a_c * exp(seq(log(0.8), log(1.25), by = 2e-6))
  dl <- diff(rt[seq_len(min(61, length(rt)))])
  dl <- dl[dl > 0]
  sc <- vapply(cand, function(a) { z <- dl / (a * res); sum((z - round(z))^2) }, numeric(1))
  a <- cand[which.min(sc)]
  g <- round((su - su[1]) / res)
  occ <- logical(max(g) + 1)
  occ[g + 1] <- TRUE
  early <- which(rt <= min(1e4, 0.05 * a * res / EARLY_NOISE))
  k <- round((rt[early] - rt[early[1]]) / (a * res))
  surv <- g
  for (n in c(10, 50, 300, length(k))) {
    use <- seq_len(min(n, length(k)))
    idx <- outer(surv, k[use], "+")
    ok <- rowSums(idx <= max(g) & matrix(occ[pmin(idx, max(g)) + 1], nrow(idx))) == length(use)
    surv <- surv[ok]
    if (length(surv) <= 1) break
  }
  list(res = res, a = a, a_coarse = a_c, score = min(sc), n_start = length(surv),
       t0 = if (length(surv) == 1) su[1] + surv * res else NA_real_, first = early[1])
}

# Candidate panel events within a raw time window, scored on scatter; best candidate per event.
best_match <- function(o, r, ord, ots, traw, fits, W, tol) {
  lo <- findInterval(traw - W, ots) + 1
  hi <- findInterval(traw + W, ots)
  cnt <- pmax(hi - lo + 1, 0)
  ridx <- rep(seq_along(traw), cnt)
  cidx <- ord[sequence(cnt, from = lo)]
  d <- matrix(0, length(ridx), length(FIT_COLS))
  for (j in seq_along(FIT_COLS)) {
    ch <- FIT_COLS[j]
    d[, j] <- abs(fits[["intercept", ch]] + fits[["slope", ch]] * o[cidx, ch] - r[ridx, ch])
  }
  worst <- apply(d, 1, max)
  keep <- worst <= tol
  ridx <- ridx[keep]; cidx <- cidx[keep]; worst <- worst[keep]
  o2 <- order(ridx, worst)
  first <- !duplicated(ridx[o2])
  row <- rep(NA_integer_, nrow(r))
  row[ridx[o2][first]] <- cidx[o2][first]
  row
}

recover_map <- function(o, r) {
  ot <- o[, "TIME"]; rt <- r[, "TIME"]
  tm <- recover_time(ot, rt)
  out <- list(time = tm, row = rep(NA_integer_, nrow(r)), fit = NULL, share_exact = 0)
  if (tm$n_start != 1) return(out)
  ord <- order(ot); ots <- ot[ord]
  # Refine the TIME map in growing stages of the authors' events (sorted by rescaled time), so that
  # extrapolation stays within half a tick: a fit on the first events fixes the next stage's cells.
  tf <- c(intercept = rt[tm$first] - tm$a * tm$t0, slope = tm$a)
  ordr <- order(rt)
  for (n in unique(pmin(round(150 * 3^(0:14)), length(rt)))) {
    E <- ordr[seq_len(n)]
    tr <- (rt[E] - tf[[1]]) / tf[[2]]
    lo <- findInterval(tr - 0.51 * tm$res, ots) + 1
    hi <- findInterval(tr + 0.51 * tm$res, ots)
    occ <- hi >= lo
    if (mean(occ) < 0.5) return(out)
    tf <- lin_fit(ots[lo[occ]], rt[E][occ])
  }
  traw <- (rt - tf[[1]]) / tf[[2]]
  # Anchors: authors' events whose time cell holds exactly one panel event; fit on the earliest,
  # then keep anchors that agree with the fit and refit.
  lo <- findInterval(traw - 0.51 * tm$res, ots) + 1
  hi <- findInterval(traw + 0.51 * tm$res, ots)
  anchor <- which(hi - lo == 0)
  if (length(anchor) < 8) return(out)
  anchor <- anchor[order(rt[anchor])]
  arow <- ord[lo[anchor]]
  fits <- vapply(FIT_COLS, function(ch) lin_fit(o[arow[1:min(300, length(arow))], ch], r[anchor[1:min(300, length(arow))], ch]), numeric(2))
  good <- apply(vapply(FIT_COLS, function(ch) abs(fits[["intercept", ch]] + fits[["slope", ch]] * o[arow, ch] - r[anchor, ch]), numeric(length(arow))), 1, max) <= PASS_TOL
  if (sum(good) < 8) return(out)
  fits <- vapply(FIT_COLS, function(ch) lin_fit(o[arow[good], ch], r[anchor[good], ch]), numeric(2))
  span <- max(ot) - min(ot)
  for (W in c(max(5 * tm$res, 1e-5 * span), max(3 * tm$res, 2e-6 * span), 3 * tm$res)) {
    row <- best_match(o, r, ord, ots, traw, fits, W, PASS_TOL)
    ok <- !is.na(row)
    if (sum(ok) < 20) return(out)
    fits <- vapply(FIT_COLS, function(ch) lin_fit(o[row[ok], ch], r[ok, ch]), numeric(2))
    tf <- lin_fit(ot[row[ok]], rt[ok])
    traw <- (rt - tf[["intercept"]]) / tf[["slope"]]
  }
  # Final exact check on all four columns with the fitted maps.
  ok <- !is.na(row)
  allfit <- cbind(TIME = tf, fits)
  resid <- vapply(MAP_COLS, function(ch) {
    abs(allfit["intercept", ch] + allfit["slope", ch] * o[row[ok], ch] - r[ok, ch])
  }, numeric(sum(ok)))
  exact <- rep(FALSE, nrow(r))
  exact[ok] <- apply(resid, 1, max) <= EXACT_TOL
  row[!exact] <- NA_integer_
  list(time = tm, row = row, fit = allfit, share_exact = mean(exact),
       n_dup_rows = sum(duplicated(row[!is.na(row)])),
       resid_median = median(resid), resid_max = max(resid))
}

# Read controls and panels per sample and instrument -------------------------------------------

sample_dirs <- sort(unique(unlist(lapply(INSTRUMENTS, function(i) {
  basename(list.dirs(file.path(data_dir, "FORT", i), recursive = FALSE))
}))))
units <- expand.grid(sample = seq_along(sample_dirs), instrument = INSTRUMENTS, stringsAsFactors = FALSE)
units$unit <- paste0("s", units$sample, "_", units$instrument)

loaded <- lapply(seq_len(nrow(units)), function(j) {
  dir_j <- file.path(data_dir, "FORT", units$instrument[j], sample_dirs[units$sample[j]])
  fs <- list.files(dir_j, pattern = "\\.(fcs|lmd)$", full.names = TRUE, ignore.case = TRUE)
  rel_f <- fs[grepl("relativi[sz]ed", fs)]
  orig_f <- fs[grepl("_ungated_.*12-panel", fs)]
  ctrl_files <- fs[grepl("_ungated_.*_[0-9]{2}-CD3-", fs)]
  if (length(rel_f) != 1 || length(orig_f) != 1 || length(ctrl_files) == 0) {
    message("skip ", units$unit[j], ": relativised ", length(rel_f), ", panel ", length(orig_f), ", controls ", length(ctrl_files))
    return(NULL)
  }
  rel <- read_exprs(rel_f)
  orig <- read_exprs(orig_f)
  rec <- recover_map(orig, rel)
  matched <- !is.na(rec$row)
  message(units$unit[j], ": matched ", sum(matched), " of ", nrow(rel))
  inv <- data.frame(
    unit = units$unit[j], n_panel = nrow(orig), n_authors = nrow(rel), n_matched = sum(matched),
    share_exact = rec$share_exact, n_dup_rows = if (is.null(rec$fit)) NA else rec$n_dup_rows,
    time_a = rec$time$a, time_res = rec$time$res, time_t0 = rec$time$t0,
    resid_median = if (is.null(rec$fit)) NA else rec$resid_median,
    resid_max = if (is.null(rec$fit)) NA else rec$resid_max
  )
  if (is.null(rec$fit) || sum(matched) == 0) {
    message("skip ", units$unit[j], ": authors' map not recovered")
    return(list(failed = inv))
  }
  # Gate box: authors' FSC-A and SSC-A quantiles, mapped back to raw units with the recovered map.
  box <- apply(rel[, c("FSC-A", "SSC-A")], 2, stats::quantile, probs = GATE_PROBS)
  for (ch in c("FSC-A", "SSC-A")) box[, ch] <- (box[, ch] - rec$fit[["intercept", ch]]) / rec$fit[["slope", ch]]
  fluor <- sub(".*-CD3-(.*)\\.fcs$", "\\1", ctrl_files, ignore.case = TRUE)
  n_ctrl <- integer(length(ctrl_files))
  controls <- lapply(seq_along(ctrl_files), function(i) {
    e <- read_exprs(ctrl_files[i])
    keep <- e[, "FSC-A"] >= box[1, 1] & e[, "FSC-A"] <= box[2, 1] &
      e[, "SSC-A"] >= box[1, 2] & e[, "SSC-A"] <= box[2, 2]
    n_ctrl[i] <<- nrow(e)
    e[keep, , drop = FALSE]
  })
  names(controls) <- paste0(units$unit[j], "_", fluor)
  share <- vapply(controls, nrow, integer(1)) / n_ctrl
  inv$n_controls <- length(controls)
  inv$min_gated_control <- min(vapply(controls, nrow, integer(1)))
  inv$gated_share_min <- min(share)
  inv$gated_share_median <- median(share)
  list(
    before = orig[rec$row[matched], , drop = FALSE],
    authors = rel[matched, , drop = FALSE],
    controls = controls,
    map = data.frame(unit = units$unit[j], channel = paste0(fluor, "-A"), control = names(controls)),
    inv = inv
  )
})
names(loaded) <- units$unit
failed <- Filter(\(x) !is.null(x$failed), loaded)
loaded <- Filter(\(x) !is.null(x) && is.null(x$failed), loaded)
inventory <- bind_rows(c(lapply(loaded, \(x) x$inv), lapply(failed, \(x) x$failed)))
write.csv(inventory, file.path(run_dir, "inventory.csv"), row.names = FALSE)
print(as.data.frame(inventory))
if (length(failed) > 0) message("DROPPED units (map not recovered): ", paste(names(failed), collapse = ", "))
units <- units[units$unit %in% names(loaded), ]
markers <- unique(unlist(lapply(loaded, \(x) x$map$channel)))
message("markers: ", paste(markers, collapse = " "))

# Fit and apply our relativisation -------------------------------------------------------------

controls <- unlist(lapply(loaded, \(x) x$controls), recursive = FALSE)
names(controls) <- unlist(lapply(loaded, \(x) names(x$controls)))
control_map <- bind_rows(lapply(loaded, \(x) x$map))
positive <- sw_relativise_split_kmeans(controls, control_map, \(x) asinh(x / COFACTOR), SEED, NSTART)
model <- sw_relativise_fit(controls, control_map, positive, min_events = MIN_EVENTS)
write.csv(model, file.path(run_dir, "relativise_fit.csv"), row.names = FALSE)
before <- lapply(loaded, \(x) x$before)
ours <- sw_relativise_apply(
  before,
  model,
  unit = setNames(names(before), names(before)),
  channels = markers,
  low = LOW,
  high = HIGH,
  input_scale = "linear"
)
batch <- setNames(units$instrument, units$unit)
diag <- sw_relativise_diagnose(model, batch = batch, order = NULL)
write.csv(diag, file.path(run_dir, "relativise_diagnose.csv"), row.names = FALSE)

# Batch scores at each stage -------------------------------------------------------------------

stages <- list(before = before, authors = lapply(loaded, \(x) x$authors), ours = ours)
scores <- lapply(names(stages), function(s) {
  sw_batch_scores(stages[[s]], batch, markers, n_events = N_EVENTS, seed = SEED, order = NULL) |>
    mutate(stage = s)
}) |>
  bind_rows()
write.csv(scores, file.path(run_dir, "batch_scores_long.csv"), row.names = FALSE)

eta <- scores |>
  filter(check == "eta2") |>
  mutate(stage = factor(stage, levels = names(stages)))
print(as.data.frame(summarise(eta, median_eta2 = median(value), max_eta2 = max(value), .by = stage)))
print(as.data.frame(tidyr::pivot_wider(select(eta, marker, stage, value), names_from = stage, values_from = value)))
print(as.data.frame(filter(diag, check == "gain_batch_ratio") |> tidyr::pivot_wider(names_from = unit, values_from = value)))
sep <- model |>
  mutate(instrument = batch[unit]) |>
  summarise(
    median_separation = median(separation),
    min_separation = min(separation),
    median_gain = median(gain),
    .by = c(instrument, channel)
  )
print(as.data.frame(sep))
write.csv(sep, file.path(run_dir, "separation_summary.csv"), row.names = FALSE)

# Figures --------------------------------------------------------------------------------------

write.csv(eta, file.path(run_dir, "figure_source_data_instrument_eta2.csv"), row.names = FALSE)
p <- ggplot(eta, aes(value, marker, colour = stage)) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = c(before = "#B2182B", authors = "#4D4D4D", ours = "#2166AC")) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "eta2: share of file-median variance explained by instrument", y = NULL, colour = NULL)
ggsave(file.path(run_dir, "fig_instrument_eta2.png"), p, width = 120, height = 80, units = "mm", dpi = 200)

sep_units <- model |> mutate(instrument = batch[unit]) |> select(unit, instrument, channel, separation)
write.csv(sep_units, file.path(run_dir, "figure_source_data_separation.csv"), row.names = FALSE)
p2 <- ggplot(sep_units, aes(separation, channel, colour = instrument)) +
  geom_point(size = 1, alpha = 0.7) +
  labs(x = "separation = gain / (sd_negative + sd_positive), one point per sample", y = NULL, colour = NULL)
ggsave(file.path(run_dir, "fig_separation.png"), p2, width = 120, height = 80, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
