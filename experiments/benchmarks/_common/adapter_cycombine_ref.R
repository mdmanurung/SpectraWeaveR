# Purpose: cyCombine fitted on calibrating reference cells only ("fair" arm), then transferred to every cell.
# Run:     sourced by a benchmark run script; nothing runs here.
# Inputs:  x (cells by markers), meta (batch, sample_id, is_reference, is_holdout_reference, arrival), markers, anchor
# Outputs: list(x, supported, abstained, reason, fits_on_heldout = FALSE, linear_max_resid)
# Pinned:  cyCombine 0.3.0 (fcd42db), wrapped through exported functions only: create_som, correct_data.
#
# cyCombine has no fit/apply split: batch_correct() refits on all cells it gets. Its steps are
#   (1) per-batch scale() of each marker, (2) kohonen SOM on the scaled data (8 x 8, rlen 10, seed 473),
#   (3) per SOM cluster, ComBat on the raw values, then a clip to that cluster's [min, max] per marker
#   (correct_data's data-frame path, which batch_correct uses for a data frame).
# Within one batch x cluster x marker, step (3) is clip(a + b * raw). So:
#   (1) uses means and SDs of calibrating cells; (2) trains on calibrating cells and maps every cell with
#   kohonen::map; (3) runs correct_data on calibrating cells only, recovers a and b from their unclipped
#   values, and applies clip(a + b * raw) to every cell of that batch x cluster with the same [min, max].
# Cells of a batch x cluster without calibrating cells (or without a recoverable map) keep their values.

CYCOMBINE_REF_SEED <- 473 # cyCombine package default (H1)
CYCOMBINE_REF_XDIM <- 8 # package default
CYCOMBINE_REF_YDIM <- 8 # package default
CYCOMBINE_REF_RLEN <- 10 # package default

run_cycombine_ref <- function(x, meta, markers, anchor) {
  cal <- meta$is_reference & !meta$is_holdout_reference
  if ("arrival" %in% names(meta)) cal <- cal & meta$arrival == "initial"
  known <- meta$batch %in% unique(meta$batch[cal])

  # (1) Per-batch scaling from calibrating cells.
  mu <- rowsum(x[cal, markers], meta$batch[cal]) / as.vector(table(meta$batch[cal]))
  sdv <- vapply(markers, \(m) tapply(x[cal, m], meta$batch[cal], stats::sd), numeric(nrow(mu)))
  z <- (x[known, markers] - mu[meta$batch[known], ]) / sdv[meta$batch[known], ]

  # (2) SOM on calibrating cells; every known cell is mapped to its nearest unit.
  cal_k <- cal[known]
  zc <- z[cal_k, , drop = FALSE]
  som <- cyCombine::create_som(
    tibble::as_tibble(zc), markers = markers, seed = CYCOMBINE_REF_SEED, rlen = CYCOMBINE_REF_RLEN,
    xdim = CYCOMBINE_REF_XDIM, ydim = CYCOMBINE_REF_YDIM, return_model = TRUE
  )
  label <- kohonen::map(som, newdata = z)$unit.classif

  # (3) ComBat on calibrating cells only (reference mode first, as in H1b).
  # Data-frame path, as batch_correct uses for a data frame.
  raw_cal <- x[known, markers][cal_k, , drop = FALSE]
  df <- tibble::as_tibble(raw_cal)
  df$batch <- meta$batch[known][cal_k]
  df$id <- seq_len(nrow(df))
  lab_cal <- label[cal_k]
  ref_batch <- anchor
  fit <- function(rb) {
    suppressMessages(cyCombine::correct_data(df, label = lab_cal, markers = markers, covar = NULL,
                                             ref.batch = rb, method = "ComBat", parametric = TRUE))
  }
  corr <- tryCatch(fit(ref_batch), error = function(e) {
    ref_batch <<- NULL
    fit(NULL)
  })
  corr <- as.matrix(corr[order(corr$id), markers])

  # Transfer: per batch x cluster x marker, a + b * raw from unclipped calibrating cells, then the clip.
  out <- x
  raw_k <- x[known, markers]
  xk <- raw_k
  bk <- meta$batch[known]
  grp_cal <- paste(bk[cal_k], label[cal_k], sep = "\r")
  grp_all <- paste(bk, label, sep = "\r")
  resid_max <- 0
  changed <- matrix(FALSE, nrow(xk), length(markers))
  for (g in unique(grp_cal)) {
    ic <- which(grp_cal == g)
    ia <- which(grp_all == g)
    lab <- lab_cal[ic[1]]
    for (j in seq_along(markers)) {
      lo <- min(raw_cal[lab_cal == lab, j])
      hi <- max(raw_cal[lab_cal == lab, j])
      r <- raw_cal[ic, j]
      y <- corr[ic, j]
      free <- !is.na(y) & y > lo & y < hi # NA: cyCombine could not correct this group
      if (sum(free) < 2 || stats::var(r[free]) == 0) next
      b <- stats::cov(r[free], y[free]) / stats::var(r[free])
      a <- mean(y[free]) - b * mean(r[free])
      resid_max <- max(resid_max, abs(y[free] - (a + b * r[free])))
      xk[ia, j] <- pmin(pmax(a + b * raw_k[ia, j], lo), hi)
      changed[ia, j] <- TRUE
    }
  }
  out[known, markers] <- xk
  # Check: on calibrating cells the transfer must reproduce cyCombine's own output.
  repro <- max(abs(xk[cal_k, ][changed[cal_k, ]] - corr[changed[cal_k, ]]), na.rm = TRUE)
  list(
    x = out, supported = TRUE, abstained = sum(!known) + sum(!changed) / length(markers),
    reason = paste0("cyCombine fit on calibrating reference cells only; ref.batch=",
                    if (is.null(ref_batch)) "NULL" else ref_batch,
                    "; affine transfer max residual ", signif(resid_max, 3),
                    "; max deviation from cyCombine on calibrating cells ", signif(repro, 3)),
    fits_on_heldout = FALSE, linear_max_resid = resid_max, repro_max = repro
  )
}
