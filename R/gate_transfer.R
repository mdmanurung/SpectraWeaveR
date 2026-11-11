# R/gate_transfer.R
# Adapt reference gates to new samples by density landmark registration on each
# gate axis. The idea follows landmark registration (Hahne et al. 2009,
# flowStats gaussNorm); the code is written from that idea.

#' Adapt reference gates to new samples by density landmarks
#'
#' @description
#' Move the vertices of reference gates onto each new sample. For every gate
#' axis the function estimates the 1D density of the parent-gated events in
#' the reference and in the sample, finds peaks (landmarks) in both, and
#' shifts the gate on that axis by how far the landmarks moved. With
#' `allow_scale = TRUE` and two landmarks it also rescales the axis about the
#' first landmark.
#'
#' Gates are adapted hierarchically. A parent is adapted first, and its
#' children take the events inside the adapted parent. The reference sample is
#' gated with the unchanged reference gates.
#'
#' @details
#' **Landmarks.** The density uses a Gaussian kernel with your `bandwidth` on
#' a common grid. A landmark is an interior local maximum, refined by a
#' parabola through the three grid points. The `n_landmarks` highest peaks are
#' kept in each sample and matched in order of position. Shift with one
#' landmark is the difference of the two peaks. With two landmarks and
#' `allow_scale = FALSE` it is the mean of the two differences. With two
#' landmarks and `allow_scale = TRUE`, the new coordinate is
#' `s1 + scale * (x - r1)` with `scale = (s2 - s1) / (r2 - r1)`, where `r` and
#' `s` are the reference and sample landmarks. `shift` in the output is
#' `s1 - r1` in that case.
#'
#' **Refusals.** A gate is refused, and returned unchanged, when any axis
#' fails. `refused_no_landmark`: fewer than `n_landmarks` peaks in the
#' reference or the sample, or fewer than 2 finite parent events.
#' `refused_max_shift`: the largest landmark displacement `|s - r|` on an axis
#' exceeds `max_shift` for that channel. An unchanged gate still passes its
#' events to its children. Read `qc` before you use `gates`.
#'
#' **Limits.** The highest peaks must be the same populations in the reference
#' and the sample. If two populations have similar size, a different peak can
#' win and the shift is wrong; the `max_shift` limit is the guard. Landmarks
#' sit on parent-gated events, so a poor parent gate carries into children.
#'
#' **Scale is yours.** `data` must be on the scale you gate on. This function
#' never transforms.
#'
#' **No defaults.** `bandwidth`, `n_landmarks`, `allow_scale` and `max_shift`
#' change the gates, so all must be passed.
#'
#' @param gates A list of `sw_gate` objects (the reference gates).
#' @param reference_data A numeric matrix of events of the reference sample.
#' @param data A named list of numeric matrices, one per sample to adapt.
#' @param bandwidth A positive number, or a named numeric vector with one
#'   value per gate channel, in channel units.
#' @param n_landmarks 1 or 2, how many peaks to match per axis.
#' @param allow_scale `TRUE` or `FALSE`. Allow a rescale between two
#'   landmarks. Needs `n_landmarks = 2` to have any effect.
#' @param max_shift A positive number, or a named numeric vector with one
#'   value per gate channel, in channel units. The largest landmark
#'   displacement that is accepted.
#'
#' @return A list with
#' * `gates`: a named list (by sample) of named lists (by gate id) of
#'   `sw_gate`. Refused gates are the reference gates.
#' * `qc`: a tibble with `sample`, `gate`, `channel`, `shift`, `scale`
#'   (both `NA` when no landmark pair was found), `n_landmarks_found` (the
#'   smaller peak count of reference and sample, capped at `n_landmarks`) and
#'   `status` (`"adapted"`, `"refused_no_landmark"` or `"refused_max_shift"`).
#'   The status is per gate and repeats on each of its channels.
#'
#' @references Hahne F, et al. (2009) Per-channel basis normalization methods
#'   for flow cytometry data. Cytometry A 77:121-131.
#'
#' @seealso [sw_gate_qc()], [sw_gate_plot()]
#'
#' @examples
#' set.seed(1)
#' mk <- function(shift) {
#'   x <- c(stats::rnorm(600, 0), stats::rnorm(300, 6))
#'   cbind(A = x + shift, B = stats::rnorm(900))
#' }
#' g <- sw_gate_range("pop", "A", -2, 2, parent = NA)
#' res <- sw_gate_transfer(
#'   list(g),
#'   reference_data = mk(0),
#'   data = list(s1 = mk(1)),
#'   bandwidth = 0.3, # illustrative
#'   n_landmarks = 1,
#'   allow_scale = FALSE,
#'   max_shift = 3 # illustrative
#' )
#' res$qc
#'
#' @export
sw_gate_transfer <- function(
  gates,
  reference_data,
  data,
  bandwidth,
  n_landmarks,
  allow_scale,
  max_shift
) {
  gl <- .sw_gq_order(gates)
  ref <- .sw_gq_matrix(reference_data)
  data <- .sw_gq_check_data(data)
  chans <- unique(unlist(lapply(gl, `[[`, "channels")))
  bw <- .sw_gq_chan_arg(bandwidth, chans, "bandwidth")
  mx <- .sw_gq_chan_arg(max_shift, chans, "max_shift")
  if (!is.numeric(n_landmarks) || length(n_landmarks) != 1 || !n_landmarks %in% c(1, 2)) {
    stop("`n_landmarks` must be 1 or 2.", call. = FALSE)
  }
  if (!is.logical(allow_scale) || length(allow_scale) != 1 || is.na(allow_scale)) {
    stop("`allow_scale` must be TRUE or FALSE.", call. = FALSE)
  }
  ref_masks <- .sw_gq_masks(gl, ref)

  adapted <- list()
  qc <- list()
  for (s in names(data)) {
    m <- data[[s]]
    cur <- gl
    inside <- list()
    for (id in names(gl)) {
      g <- gl[[id]]
      absent <- setdiff(g$channels, colnames(m))
      if (length(absent) > 0) {
        stop("Sample '", s, "' lacks channels of gate '", id, "': ", paste(absent, collapse = ", "), call. = FALSE)
      }
      pm <- if (is.na(g$parent)) rep(TRUE, nrow(m)) else inside[[g$parent]]
      ref_par <- ref[ref_masks[[id]]$parent, , drop = FALSE]
      smp_par <- m[pm, , drop = FALSE]
      fits <- lapply(g$channels, function(ch) {
        .sw_gq_register(
          ref_par[, ch], smp_par[, ch], bw[[ch]], n_landmarks, allow_scale, mx[[ch]]
        )
      })
      names(fits) <- g$channels
      reasons <- vapply(fits, function(f) f$status, character(1))
      status <- if (any(reasons == "refused_no_landmark")) {
        "refused_no_landmark"
      } else if (any(reasons == "refused_max_shift")) {
        "refused_max_shift"
      } else {
        "adapted"
      }
      if (status == "adapted") cur[[id]] <- .sw_gq_move_gate(g, fits)
      inside[[id]] <- pm & sw_gate_apply(cur[[id]], m)
      for (ch in g$channels) {
        f <- fits[[ch]]
        qc[[length(qc) + 1L]] <- data.frame(
          sample = s, gate = id, channel = ch, shift = f$shift, scale = f$scale,
          n_landmarks_found = f$found, status = status, stringsAsFactors = FALSE
        )
      }
    }
    adapted[[s]] <- cur
  }
  list(gates = adapted, qc = tibble::as_tibble(do.call(rbind, qc)))
}

# Interior local maxima of a Gaussian-kernel density; the n highest, by position.
.sw_gq_peaks <- function(x, bw, from, to, n_peaks) {
  # Gaussian kernel density on a grid, by linear binning and a direct kernel
  # sum. stats::density() uses an FFT whose rounding noise makes false ripples
  # in the empty tails, which would read as peaks.
  k <- 1024L
  gx <- seq(from, to, length.out = k)
  h <- gx[2] - gx[1]
  pos0 <- (x - from) / h + 1
  lo <- pmin(floor(pos0), k - 1L)
  fr <- pos0 - lo
  b <- rowsum(c(1 - fr, fr), c(lo, lo + 1L))
  w <- numeric(k)
  w[as.integer(rownames(b))] <- b[, 1]
  y <- as.vector(stats::dnorm(outer(gx, gx, "-"), sd = bw) %*% w) / length(x)
  d <- list(x = gx, y = y)
  i <- which(y[2:(k - 1)] > y[1:(k - 2)] & y[2:(k - 1)] >= y[3:k]) + 1L
  if (length(i) == 0) return(list(pos = numeric(), n_all = 0L))
  # Parabola through three grid points for a position between grid steps.
  h <- d$x[2] - d$x[1]
  denom <- y[i - 1L] - 2 * y[i] + y[i + 1L]
  off <- ifelse(denom < 0, 0.5 * (y[i - 1L] - y[i + 1L]) / denom, 0)
  pos <- d$x[i] + off * h
  top <- order(y[i], decreasing = TRUE)[seq_len(min(n_peaks, length(i)))]
  list(pos = sort(pos[top]), n_all = length(i))
}

# Landmark registration of one axis: how the sample moved against the reference.
.sw_gq_register <- function(xr, xs, bw, n_landmarks, allow_scale, max_shift) {
  xr <- xr[is.finite(xr)]
  xs <- xs[is.finite(xs)]
  none <- list(shift = NA_real_, scale = NA_real_, anchor = NA_real_, found = 0L, status = "refused_no_landmark")
  if (length(xr) < 2 || length(xs) < 2) return(none)
  from <- min(c(xr, xs)) - 3 * bw
  to <- max(c(xr, xs)) + 3 * bw
  pr <- .sw_gq_peaks(xr, bw, from, to, n_landmarks)
  ps <- .sw_gq_peaks(xs, bw, from, to, n_landmarks)
  found <- as.integer(min(length(pr$pos), length(ps$pos), n_landmarks))
  if (found < n_landmarks) {
    none$found <- found
    return(none)
  }
  r <- pr$pos
  s <- ps$pos
  if (n_landmarks == 2 && allow_scale) {
    scale <- (s[2] - s[1]) / (r[2] - r[1])
    shift <- s[1] - r[1]
    anchor <- r[1]
  } else {
    scale <- 1
    shift <- mean(s - r)
    anchor <- 0
  }
  status <- if (max(abs(s - r)) > max_shift) "refused_max_shift" else "adapted"
  list(shift = shift, scale = scale, anchor = anchor, found = found, status = status)
}

# new = anchor + shift + scale * (x - anchor), per gate axis.
.sw_gq_move_gate <- function(g, fits) {
  v <- g$vertices
  for (k in seq_along(g$channels)) {
    f <- fits[[g$channels[k]]]
    v[, k] <- f$anchor + f$shift + f$scale * (v[, k] - f$anchor)
  }
  if (g$type == "range") {
    sw_gate_range(g$id, g$channels, v["min", 1], v["max", 1], g$parent)
  } else {
    sw_gate_polygon(g$id, g$channels, v, g$parent)
  }
}
