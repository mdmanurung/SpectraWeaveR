# Purpose: shared by 87_af405_forced_af_common.R and 87b_af405_afremove_common.R. Reads an arm's events for
#          one donor from the common gated dir (rows paired across arms, G8) and scores each channel's
#          negatives with the 39b / 51 metrics. Replaces helper75_marker-scores.R, which rebuilt the rows
#          from clean2 files; here the rows are the common GatingSet's.
# Run:     source("experiments/benchmarks/05_bmv-unmix-compare/helper87c_marker-scores-common.R") after the
#          caller has set COFACTOR, CUT, SQUEEZE_FRAC, COMMON, STORED_ARCSINH and POP_PATH.
# Inputs:  <COMMON>/gs/<arm>   common-gate GatingSets [batch 0]
# Outputs: none (functions only)

# One arm and donor, selected channels, linear scale (the gs stores STORED_ARCSINH channels as arcsinh).
gs_matrix <- function(arm, sample, channels) {
  gs <- load_gs(file.path(COMMON, "gs", arm))
  e <- exprs(gh_pop_get_data(gs[[sample]], "root"))[, channels, drop = FALSE]
  stored <- colnames(e) %in% STORED_ARCSINH
  e[, stored] <- COFACTOR * sinh(e[, stored])
  e
}

# 39b / 51 readouts per channel on given rows: Nguyen rSD (q84 - q50) of events below arcsinh CUT, its ratio
# to as188_joint on the same rows, share with |x| < SQUEEZE_FRAC x as188_joint's rSD, share above arcsinh CUT
nguyen <- \(v) unname(diff(quantile(v, c(0.5, 0.84))))
score_channels <- function(x, x_ref, rows) {
  colnames(x) |>
    set_names() |>
    map(\(ch) {
      lin <- x[rows, ch]
      lin_ref <- x_ref[rows, ch]
      ref_rsd <- nguyen(lin_ref[asinh(lin_ref / COFACTOR) < CUT])
      tibble(
        n = length(rows),
        rsd_neg = nguyen(lin[asinh(lin / COFACTOR) < CUT]),
        ref_rsd = ref_rsd,
        rsd_ratio = rsd_neg / ref_rsd,
        squeezed_share = mean(abs(lin) < SQUEEZE_FRAC * ref_rsd),
        pos_share = mean(asinh(lin / COFACTOR) > CUT)
      )
    }) |>
    list_rbind(names_to = "channel")
}
