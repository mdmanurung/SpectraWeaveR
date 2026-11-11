# experiments/integration/glehr2026/06_gate_benchmark.R
# Purpose: real-data gating benchmark. The authors' relativised panel files (Glehr et al. 2026) hold
#          exactly their lymphocyte -> singlet -> CD3+ events, so they are an expert gate. For each
#          unit (BAD: donor x instrument; FORT: sample x instrument) the candidate events are the
#          original panel events inside the authors' scatter box, and an event is CD3+ if it matches
#          an authors' event. Two gates are scored against that truth on asinh(AA750-A / 6000):
#          A, sw_gate_density_cut (valley) per unit; B, the reference unit's cut moved to the other
#          units with sw_gate_transfer (same-instrument reference, and a cross-instrument reference).
# Run:     Rscript experiments/integration/glehr2026/06_gate_benchmark.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*.fcs   12-panel (original, compensated) and the
#                                                      authors' *relativised* panel
#          <data_dir>/FORT/<instrument>/<sample>/*.fcs 12-panel and the authors' *relativized* panel
# Outputs: <run_dir>/inventory.csv      per unit: events, box, matched authors' events, status
#          <run_dir>/cuts.csv           per unit and method: cut and the reason when it is NA
#          <run_dir>/metrics.csv        per unit and method: F1, precision, recall, cut offset, box share
#          <run_dir>/summary.csv        per dataset, instrument and method: median and range of F1
#          <run_dir>/fallback_count.csv units that used the transfer fallback (density_or_transfer)
#          <run_dir>/gate_qc_bad_cross.csv  sw_gate_qc numbers of the BAD cross transfer, with F1
#          <run_dir>/fig_gate_f1.png, figure_source_data_gate_f1.csv
#          <run_dir>/session_info.txt
#
# Matching of the authors' events (copied from 03_instrument_batch.R and 05_instrument_batch_fort.R):
#   BAD: exact match on TIME, FSC-A, FSC-H, SSC-A, SSC-H; original events whose key is not unique are
#        not candidates (their truth is ambiguous) and are counted in inventory.csv.
#   FORT: the authors rescaled TIME, FSC-A, SSC-A, FSC-W per file; a per-file linear map is recovered
#        (see 05) and events are matched on all four columns within EXACT_TOL.
#
# Caveat: the candidate box is not the authors' exact lymphocyte and singlet gate. CD3+ events inside
#   the box that are not in the authors' file may lie outside their singlet gate. metrics.csv holds
#   share_extra_above_thr: of the events above a method's cut, the share that is not in the authors'
#   file but lies at or above the authors' implicit threshold.

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

# Knobs (all user-approved 2026-11-10) ---------------------------------------------------------

COFACTOR <- 6000 # user-approved: asinh(x / 6000) of the CD3 channel
CD3_CHANNEL <- "AA750-A" # user-approved: CD3 channel in the panel, original compensated values
GATE_PROBS <- c(0.01, 0.99) # user-approved: box = 1-99% of FSC-A x SSC-A of the authors' events
DENSITY_METHOD <- "valley" # user-approved
MIN_PEAK_HEIGHT <- 0.05 # user-approved
PEAK_DISTANCE_BW <- 2 # user-approved: min_peak_distance = 2 * bandwidth
N_LANDMARKS <- 2 # user-approved
MAX_SHIFT <- Inf # user-approved
DATASETS <- list(BAD = c("aurora", "fortessa"), FORT = c("fortessa", "lx", "navios"))
CROSS_REFERENCE <- c(BAD = "aurora", FORT = "fortessa") # user-approved: first unit of this instrument
KEY <- c("TIME", "FSC-A", "FSC-H", "SSC-A", "SSC-H") # BAD matching
MAP_COLS <- c("TIME", "FSC-A", "SSC-A", "FSC-W") # FORT matching, as in 05
FIT_COLS <- c("FSC-A", "SSC-A", "FSC-W")
EXACT_TOL <- 0.01 # as in 05
PASS_TOL <- 0.02 # as in 05
EARLY_NOISE <- 6e-8 # as in 05

read_exprs <- function(path) {
  exprs(read.FCS(path, truncate_max_range = FALSE, transformation = FALSE, emptyValue = FALSE))
}

# Recover the authors' linear rescale for FORT (copied from 05) --------------------------------

lin_fit <- function(x, y) {
  f <- lm.fit(cbind(1, x), y)
  c(intercept = f$coefficients[[1]], slope = f$coefficients[[2]])
}

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
  ok <- !is.na(row)
  allfit <- cbind(TIME = tf, fits)
  resid <- vapply(MAP_COLS, function(ch) {
    abs(allfit["intercept", ch] + allfit["slope", ch] * o[row[ok], ch] - r[ok, ch])
  }, numeric(sum(ok)))
  exact <- rep(FALSE, nrow(r))
  exact[ok] <- apply(resid, 1, max) <= EXACT_TOL
  row[!exact] <- NA_integer_
  list(time = tm, row = row, fit = allfit, share_exact = mean(exact))
}

# Candidate events and truth per unit ----------------------------------------------------------
# Each unit holds x = asinh(CD3 / cofactor) of the candidates, truth (matches an authors' event),
# thr = the authors' implicit threshold (min asinh CD3 over their matched events), and counts.

dirs <- lapply(names(DATASETS), \(ds) {
  sort(unique(unlist(lapply(DATASETS[[ds]], \(i) {
    basename(list.dirs(file.path(data_dir, ds, i), recursive = FALSE))
  }))))
})
names(dirs) <- names(DATASETS)

units <- bind_rows(lapply(names(DATASETS), \(ds) {
  expand.grid(index = seq_along(dirs[[ds]]), instrument = DATASETS[[ds]], stringsAsFactors = FALSE) |>
    mutate(dataset = ds, unit = paste0(if_else(ds == "BAD", "d", "s"), index, "_", instrument))
}))

loaded <- lapply(seq_len(nrow(units)), function(j) {
  ds <- units$dataset[j]
  dir_j <- file.path(data_dir, ds, units$instrument[j], dirs[[ds]][units$index[j]])
  fs <- list.files(dir_j, pattern = "\\.(fcs|lmd)$", full.names = TRUE, ignore.case = TRUE)
  rel_f <- fs[grepl("relativi[sz]ed", fs)]
  orig_f <- fs[grepl("_ungated_.*12-panel", fs)]
  base <- data.frame(unit = units$unit[j], dataset = ds, instrument = units$instrument[j])
  if (length(rel_f) != 1 || length(orig_f) != 1) {
    return(list(inv = mutate(base, status = "panel or relativised file missing"), x = NULL))
  }
  rel <- read_exprs(rel_f)
  orig <- read_exprs(orig_f)
  box <- apply(rel[, c("FSC-A", "SSC-A")], 2, stats::quantile, probs = GATE_PROBS)
  if (ds == "BAD") {
    k_orig <- do.call(paste, as.data.frame(orig[, KEY]))
    k_rel <- do.call(paste, as.data.frame(rel[, KEY]))
    uniq <- !(k_orig %in% k_orig[duplicated(k_orig)])
    hit <- match(k_rel, k_orig[uniq])
    pool <- orig[uniq, , drop = FALSE]
    truth_rows <- hit[!is.na(hit)]
    n_ambiguous <- sum(!uniq)
  } else {
    rec <- recover_map(orig, rel)
    if (is.null(rec$fit) || all(is.na(rec$row))) {
      return(list(inv = mutate(base, n_authors = nrow(rel), status = "authors' map not recovered"), x = NULL))
    }
    for (ch in c("FSC-A", "SSC-A")) box[, ch] <- (box[, ch] - rec$fit[["intercept", ch]]) / rec$fit[["slope", ch]]
    pool <- orig
    truth_rows <- rec$row[!is.na(rec$row)]
    n_ambiguous <- 0L
  }
  in_box <- pool[, "FSC-A"] >= box[1, 1] & pool[, "FSC-A"] <= box[2, 1] &
    pool[, "SSC-A"] >= box[1, 2] & pool[, "SSC-A"] <= box[2, 2]
  cand <- which(in_box)
  x_all <- asinh(pool[, CD3_CHANNEL] / COFACTOR)
  truth_rows <- unique(truth_rows)
  inv <- mutate(
    base,
    n_panel = nrow(orig), n_authors = nrow(rel), n_authors_matched = length(truth_rows),
    n_ambiguous_dropped = n_ambiguous, n_candidates = length(cand),
    n_authors_in_box = sum(truth_rows %in% cand), status = "ok"
  )
  list(
    inv = inv,
    x = x_all[cand],
    truth = cand %in% truth_rows,
    thr = min(x_all[truth_rows])
  )
})
names(loaded) <- units$unit
inventory <- bind_rows(lapply(loaded, `[[`, "inv"))
write.csv(inventory, file.path(run_dir, "inventory.csv"), row.names = FALSE)
print(as.data.frame(inventory))

# Method A: per-unit density cut ---------------------------------------------------------------

density_cut <- function(x) {
  bw <- stats::bw.nrd0(x) # Silverman, per unit
  sw_gate_density_cut(
    x, method = DENSITY_METHOD, bandwidth = bw, min_peak_height = MIN_PEAK_HEIGHT,
    min_peak_distance = PEAK_DISTANCE_BW * bw, prob = NULL, k = NULL, mad_side = NULL
  )
}

cuts_a <- bind_rows(lapply(units$unit, \(u) {
  x <- loaded[[u]]$x
  if (is.null(x)) {
    return(data.frame(unit = u, method = "density", cut = NA_real_, reason = loaded[[u]]$inv$status))
  }
  cut <- density_cut(x)
  data.frame(unit = u, method = "density", cut = as.numeric(cut), reason = attr(cut, "reason"))
}))

# Method B: reference cut moved to target units ------------------------------------------------

transfer_cuts <- function(ref_unit, target_units, allow_scale, method) {
  na_rows <- \(why) data.frame(unit = target_units, method = method, cut = NA_real_, reason = why)
  ref_cut <- cuts_a$cut[cuts_a$unit == ref_unit]
  if (is.null(loaded[[ref_unit]]$x)) return(na_rows(paste("reference unit unavailable:", ref_unit)))
  if (is.na(ref_cut)) return(na_rows(paste("reference cut is NA:", cuts_a$reason[cuts_a$unit == ref_unit])))
  has_x <- !vapply(loaded[target_units], \(l) is.null(l$x), logical(1))
  res <- sw_gate_transfer(
    list(sw_gate_range("cd3pos", "CD3", min = ref_cut, max = Inf, parent = NA)),
    reference_data = cbind(CD3 = loaded[[ref_unit]]$x),
    data = lapply(loaded[target_units[has_x]], \(l) cbind(CD3 = l$x)),
    bandwidth = stats::bw.nrd0(loaded[[ref_unit]]$x), # Silverman of the reference unit
    n_landmarks = N_LANDMARKS,
    allow_scale = allow_scale,
    max_shift = MAX_SHIFT
  )
  out <- bind_rows(lapply(names(res$gates), \(u) {
    status <- res$qc$status[res$qc$sample == u]
    data.frame(
      unit = u, method = method,
      cut = if_else(status == "adapted", res$gates[[u]]$cd3pos$vertices["min", 1], NA_real_),
      reason = if_else(status == "adapted", NA_character_, status)
    )
  }))
  out <- bind_rows(out, na_rows("target unit unavailable")[!has_x, ]) |>
    mutate(reference = ref_unit)
  attr(out, "gates") <- res$gates
  out
}

# Within-instrument reference: the first unit (unit order) of that instrument with a valley cut
# (user-approved 2026-11-10, "transfer from a gated sample"); NA when no unit has one.
within_reference <- \(ds, ins) {
  u <- units$unit[units$dataset == ds & units$instrument == ins]
  ok <- u[!is.na(cuts_a$cut[match(u, cuts_a$unit)])]
  if (length(ok) > 0) ok[1] else NA_character_
}
first_unit <- \(ds, ins) units$unit[units$dataset == ds & units$instrument == ins & units$index == 1]
cuts_b <- bind_rows(lapply(names(DATASETS), \(ds) {
  within <- lapply(DATASETS[[ds]], \(ins) {
    ref <- within_reference(ds, ins)
    all_u <- units$unit[units$dataset == ds & units$instrument == ins]
    if (is.na(ref)) {
      return(bind_rows(lapply(c("noscale", "scale"), \(v) {
        data.frame(
          unit = all_u, method = paste0("transfer_within_", v), cut = NA_real_,
          reason = "no unit of this instrument has a valley cut", reference = NA_character_
        )
      })))
    }
    tg <- setdiff(all_u, ref)
    bind_rows(
      transfer_cuts(ref, tg, FALSE, "transfer_within_noscale"),
      transfer_cuts(ref, tg, TRUE, "transfer_within_scale")
    )
  })
  ref <- first_unit(ds, CROSS_REFERENCE[[ds]])
  tg <- units$unit[units$dataset == ds & units$instrument != CROSS_REFERENCE[[ds]]]
  cross <- bind_rows(
    transfer_cuts(ref, tg, FALSE, "transfer_cross_noscale"),
    transfer_cuts(ref, tg, TRUE, "transfer_cross_scale")
  )
  bind_rows(within, list(cross))
}))

# Method C: own valley cut when it exists, else the within-instrument transfer ------------------

cuts_c <- bind_rows(lapply(c("noscale", "scale"), \(v) {
  tr <- filter(cuts_b, method == paste0("transfer_within_", v))
  cuts_a |>
    select(unit, own_cut = cut) |>
    left_join(
      select(tr, unit, tr_cut = cut, tr_reason = reason, tr_reference = reference),
      by = "unit", relationship = "one-to-one"
    ) |>
    mutate(
      method = paste0("density_or_transfer_", v),
      used_fallback = is.na(own_cut),
      cut = if_else(used_fallback, tr_cut, own_cut),
      reason = if_else(used_fallback & is.na(cut), coalesce(tr_reason, "no transfer result"), NA_character_),
      reference = if_else(used_fallback, tr_reference, NA_character_)
    ) |>
    select(unit, method, cut, reason, reference, used_fallback)
}))
fallback_count <- cuts_c |>
  left_join(select(units, unit, dataset, instrument), by = "unit", relationship = "many-to-one") |>
  summarise(n_units = n(), n_fallback = sum(used_fallback), n_fallback_still_na = sum(used_fallback & is.na(cut)), .by = c(dataset, instrument, method))
write.csv(fallback_count, file.path(run_dir, "fallback_count.csv"), row.names = FALSE)
print(as.data.frame(fallback_count))

cuts <- bind_rows(cuts_a, cuts_b, cuts_c) |>
  left_join(select(units, unit, dataset, instrument), by = "unit", relationship = "many-to-one")
write.csv(cuts, file.path(run_dir, "cuts.csv"), row.names = FALSE)

# Metrics against the authors' gate ------------------------------------------------------------

metrics <- bind_rows(lapply(seq_len(nrow(cuts)), \(i) {
  l <- loaded[[cuts$unit[i]]]
  cut <- cuts$cut[i]
  if (is.null(l$x) || is.na(cut)) {
    return(data.frame(i = i))
  }
  above <- l$x >= cut
  tp <- sum(above & l$truth)
  data.frame(
    i = i,
    n_candidates = length(l$x), n_truth = sum(l$truth), n_above = sum(above),
    precision = if (sum(above) > 0) tp / sum(above) else NA_real_,
    recall = tp / sum(l$truth),
    f1 = 2 * tp / (sum(above) + sum(l$truth)),
    thr_authors = l$thr,
    cut_offset = cut - l$thr,
    share_extra_above_thr = if (sum(above) > 0) sum(above & !l$truth & l$x >= l$thr) / sum(above) else NA_real_
  )
}))
metrics <- cuts |>
  mutate(i = row_number()) |>
  left_join(metrics, by = "i", relationship = "one-to-one") |>
  select(-i)
write.csv(metrics, file.path(run_dir, "metrics.csv"), row.names = FALSE)

# Gate QC of the BAD cross-instrument transfer (shift-only vs scale) ---------------------------
# Reference = the Aurora unit the gate came from. edge_width = one Silverman bandwidth of that unit;
# this value is illustrative (not a user-approved choice), it only sets edge_share. The QC reads
# kept_share and kept_share_z; no cutoff is applied. One call per variant, so z is across the
# reference and the nine Fortessa units of that variant.

qc_ref <- first_unit("BAD", CROSS_REFERENCE[["BAD"]])
qc_targets <- units$unit[units$dataset == "BAD" & units$instrument != CROSS_REFERENCE[["BAD"]]]
qc <- bind_rows(lapply(c(noscale = FALSE, scale = TRUE), \(allow_scale) {
  tr <- transfer_cuts(qc_ref, qc_targets, allow_scale, "x")
  variant <- if (allow_scale) "scale" else "noscale"
  gates <- c(
    setNames(list(list(sw_gate_range("cd3pos", "CD3", min = cuts_a$cut[cuts_a$unit == qc_ref], max = Inf, parent = NA))), qc_ref),
    attr(tr, "gates")
  )
  sw_gate_qc(
    data = lapply(loaded[names(gates)], \(l) cbind(CD3 = l$x)),
    gates = gates,
    reference = qc_ref,
    edge_width = stats::bw.nrd0(loaded[[qc_ref]]$x)
  ) |>
    filter(check %in% c("kept_share", "kept_share_z", "edge_share", "centroid_shift")) |>
    mutate(variant = variant)
}))
qc_wide <- qc |>
  select(sample, variant, check, value) |>
  tidyr::pivot_wider(names_from = check, values_from = value) |>
  left_join(
    metrics |>
      filter(method %in% c("transfer_cross_noscale", "transfer_cross_scale"), dataset == "BAD") |>
      transmute(sample = unit, variant = sub("transfer_cross_", "", method), f1, true_share = n_truth / n_candidates),
    by = c("sample", "variant"), relationship = "one-to-one"
  )
write.csv(qc_wide, file.path(run_dir, "gate_qc_bad_cross.csv"), row.names = FALSE)
print(as.data.frame(qc_wide))

# Summary --------------------------------------------------------------------------------------

summary_tbl <- metrics |>
  summarise(
    n_units = n(), n_na = sum(is.na(cut)),
    n_fallback = sum(used_fallback %in% TRUE),
    na_reasons = paste(unique(reason[is.na(cut)]), collapse = "; "),
    f1_median = median(f1, na.rm = TRUE),
    f1_min = if (all(is.na(f1))) NA_real_ else min(f1, na.rm = TRUE),
    f1_max = if (all(is.na(f1))) NA_real_ else max(f1, na.rm = TRUE),
    cut_offset_median = median(cut_offset, na.rm = TRUE),
    share_extra_median = median(share_extra_above_thr, na.rm = TRUE),
    .by = c(dataset, instrument, method)
  )
write.csv(summary_tbl, file.path(run_dir, "summary.csv"), row.names = FALSE)
print(as.data.frame(summary_tbl))

# Figure: F1 per unit by method and instrument -------------------------------------------------

METHOD_LEVELS <- c(
  "density", "transfer_within_noscale", "transfer_within_scale",
  "transfer_cross_noscale", "transfer_cross_scale",
  "density_or_transfer_noscale", "density_or_transfer_scale"
)
METHOD_COLOURS <- c(
  density = "#4D4D4D", transfer_within_noscale = "#2166AC", transfer_within_scale = "#67A9CF",
  transfer_cross_noscale = "#B2182B", transfer_cross_scale = "#EF8A62",
  density_or_transfer_noscale = "#1B7837", density_or_transfer_scale = "#7FBF7B"
)
plot_df <- metrics |>
  filter(!is.na(f1)) |>
  mutate(
    method = forcats::fct(method, levels = METHOD_LEVELS),
    panel = paste(dataset, instrument)
  ) |>
  select(unit, panel, method, f1)
write.csv(plot_df, file.path(run_dir, "figure_source_data_gate_f1.csv"), row.names = FALSE)
p <- ggplot(plot_df, aes(method, f1, colour = method)) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 1.2) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.6, linewidth = 0.3, colour = "black") +
  facet_wrap(vars(panel), nrow = 1) +
  scale_colour_manual(values = METHOD_COLOURS, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "F1 vs authors' gate") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(run_dir, "fig_gate_f1.png"), p, width = 280, height = 110, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
