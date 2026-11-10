# R/gate_qc.R
# Numbers that describe how a set of gates sits on each sample: how many parent
# events, how many are kept, how close they lie to the gate edge, and how far
# the in-gate centre moved against a reference sample. No pass or fail.

#' Gate QC numbers per sample and gate
#'
#' @description
#' Describe each gate on each sample with plain numbers. Gates are applied
#' hierarchically: the events of a gate are those inside its parent chain
#' (matched by gate id). The result is numbers only; no gate or sample is
#' marked bad.
#'
#' * `n_parent`: events inside the parent chain (all events for a root gate).
#' * `kept_share`: events inside the gate divided by `n_parent`.
#' * `edge_share`: share of parent events whose distance to the gate boundary
#'   is at most `edge_width`. For a range gate the distance is to the nearest
#'   finite bound. For a polygon it is the distance to the nearest edge
#'   segment, in raw units of the two channels. A high value means the gate
#'   cuts through dense events, so a small shift moves many events in or out.
#' * `kept_share_z`: robust z of `kept_share` across samples,
#'   `(v - median) / (1.4826 * MAD)`. `NA` with fewer than 3 samples or when
#'   the MAD is 0.
#' * `kept_share_log_ratio`: `log(kept_share / kept_share of the reference)`.
#'   Only with `reference`. Defined even when every sample keeps the same
#'   share, where `kept_share_z` is `NA`. Meaningful when the parent population
#'   has the same make-up as in the reference (same instrument and parent
#'   gate). Across instruments it measures make-up differences too: on the
#'   Glehr BAD data a failed shift-only Aurora-to-Fortessa transfer scored
#'   +0.45 and the good scaled one -0.46 to -1.46. There the failed gate was
#'   visible as `kept_share` = 1 in every sample (a gate that keeps all).
#' * `centroid_shift`: in-gate median of a gate channel minus the same median
#'   in the `reference` sample (signed, one row per channel).
#' * `centroid_dist`: Euclidean length of the `centroid_shift` values of a
#'   gate (`channel` is `NA`).
#'
#' @details
#' **Scale is yours.** `data` must be on the scale you gate on. This function
#' never transforms.
#'
#' **No defaults.** `edge_width` and `reference` change the numbers, so both
#' must be passed. The reference sample is gated with its own gate list (the
#' shared list, or its entry in a per-sample list).
#'
#' @param data A named list. Each element is a numeric matrix (or data frame)
#'   with one column per channel. Names are the sample ids.
#' @param gates A list of `sw_gate` objects shared by all samples, or a named
#'   list with one such list per sample (for example adapted gates from
#'   [sw_gate_transfer()]).
#' @param reference `NULL`, or the name of one sample in `data`. When given,
#'   the `kept_share_log_ratio`, `centroid_shift` and `centroid_dist` checks
#'   are added.
#' @param edge_width A positive number in channel units, or a named numeric
#'   vector with one value per gate channel. A polygon needs one width, so its
#'   two channels must carry the same value.
#'
#' @return A tibble with columns `sample`, `gate`, `check`, `channel` (`NA`
#'   when the check is not per channel) and `value`.
#'
#' @seealso [sw_gate_polygon()], [sw_gate_range()], [sw_gate_transfer()]
#'
#' @examples
#' set.seed(1)
#' mk <- function(shift) {
#'   matrix(stats::rnorm(400, mean = shift), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' }
#' d <- list(s1 = mk(0), s2 = mk(0.2), s3 = mk(0.1))
#' g <- sw_gate_polygon(
#'   "pop", c("A", "B"),
#'   rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1)),
#'   parent = NA
#' )
#' sw_gate_qc(d, list(g), reference = "s1", edge_width = 0.1) # illustrative
#'
#' @export
sw_gate_qc <- function(data, gates, reference, edge_width) {
  data <- .sw_gq_check_data(data)
  samples <- names(data)
  per_sample <- .sw_gq_per_sample(gates, samples)
  if (!is.null(reference) &&
    (!is.character(reference) || length(reference) != 1 || !reference %in% samples)) {
    stop("`reference` must be NULL or one sample name present in `data`.", call. = FALSE)
  }
  chans_all <- unique(unlist(lapply(per_sample, function(gl) {
    unlist(lapply(gl, `[[`, "channels"))
  })))
  edge_w <- .sw_gq_chan_arg(edge_width, chans_all, "edge_width")

  rows <- list()
  add <- function(s, g, check, channel, value) {
    rows[[length(rows) + 1L]] <<- data.frame(
      sample = s, gate = g, check = check, channel = channel,
      value = as.numeric(value), stringsAsFactors = FALSE
    )
  }
  medians <- list()
  for (s in samples) {
    m <- data[[s]]
    gl <- per_sample[[s]]
    masks <- .sw_gq_masks(gl, m)
    for (id in names(gl)) {
      g <- gl[[id]]
      pm <- masks[[id]]$parent
      ins <- masks[[id]]$inside
      n_par <- sum(pm)
      add(s, id, "n_parent", NA_character_, n_par)
      add(s, id, "kept_share", NA_character_, if (n_par > 0) sum(ins) / n_par else NA_real_)
      w <- .sw_gq_width(edge_w, g)
      dist <- .sw_gq_edge_dist(g, m[pm, , drop = FALSE])
      add(s, id, "edge_share", NA_character_, if (n_par > 0) sum(dist <= w, na.rm = TRUE) / n_par else NA_real_)
      medians[[paste(s, id, sep = "\r")]] <- vapply(g$channels, function(ch) {
        if (any(ins)) stats::median(m[ins, ch], na.rm = TRUE) else NA_real_
      }, numeric(1))
    }
  }
  out <- do.call(rbind, rows)

  # Robust z of kept_share across samples, per gate.
  ks <- out[out$check == "kept_share", ]
  z_rows <- lapply(split(ks, ks$gate), function(d) {
    v <- d$value
    ok <- is.finite(v)
    med <- stats::median(v[ok])
    sc <- 1.4826 * stats::median(abs(v[ok] - med))
    z <- if (sum(ok) >= 3 && is.finite(sc) && sc > 0) (v - med) / sc else rep(NA_real_, length(v))
    data.frame(
      sample = d$sample, gate = d$gate, check = "kept_share_z",
      channel = NA_character_, value = z, stringsAsFactors = FALSE
    )
  })
  out <- rbind(out, do.call(rbind, z_rows))

  if (!is.null(reference)) {
    # Kept share against the reference sample, as a log ratio. Unlike
    # kept_share_z it is defined when every sample keeps the same share: a
    # transferred gate that keeps all events in every sample (MAD 0) shows here.
    ref_ks <- ks[ks$sample == reference, c("gate", "value")]
    ratio <- merge(ks[c("sample", "gate", "value")], ref_ks, by = "gate", suffixes = c("", "_ref"))
    if (nrow(ratio) > 0) {
      out <- rbind(out, data.frame(
        sample = ratio$sample, gate = ratio$gate, check = "kept_share_log_ratio",
        channel = NA_character_,
        value = ifelse(ratio$value > 0 & ratio$value_ref > 0, log(ratio$value / ratio$value_ref), NA_real_),
        stringsAsFactors = FALSE
      ))
    }
    cen <- list()
    for (s in samples) {
      for (id in names(per_sample[[s]])) {
        mine <- medians[[paste(s, id, sep = "\r")]]
        ref_gl <- per_sample[[reference]]
        if (!id %in% names(ref_gl)) next
        ref <- medians[[paste(reference, id, sep = "\r")]]
        if (!identical(names(mine), names(ref))) next
        diffs <- mine - ref
        for (ch in names(diffs)) {
          cen[[length(cen) + 1L]] <- data.frame(
            sample = s, gate = id, check = "centroid_shift", channel = ch,
            value = unname(diffs[ch]), stringsAsFactors = FALSE
          )
        }
        cen[[length(cen) + 1L]] <- data.frame(
          sample = s, gate = id, check = "centroid_dist", channel = NA_character_,
          value = sqrt(sum(diffs^2)), stringsAsFactors = FALSE
        )
      }
    }
    if (length(cen) > 0) out <- rbind(out, do.call(rbind, cen))
  }

  lev <- c("n_parent", "kept_share", "edge_share", "kept_share_z", "kept_share_log_ratio", "centroid_shift", "centroid_dist")
  out <- out[order(match(out$sample, samples), match(out$gate, unique(out$gate)), match(out$check, lev)), ]
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

# ---- shared helpers for the gate modules ----------------------------------

# Data: a non-empty named list of numeric matrices with column names.
.sw_gq_check_data <- function(data) {
  if (!is.list(data) || length(data) == 0 || is.null(names(data)) ||
    anyNA(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    stop("`data` must be a named list with unique, non-empty names.", call. = FALSE)
  }
  lapply(data, .sw_gq_matrix)
}

.sw_gq_matrix <- function(x) {
  m <- as.matrix(x)
  if (!is.numeric(m) || is.null(colnames(m))) {
    stop("Each sample must be a numeric matrix with column names.", call. = FALSE)
  }
  m
}

.sw_gq_is_gate_list <- function(x) {
  is.list(x) && !inherits(x, "sw_gate") && length(x) > 0 &&
    all(vapply(x, inherits, logical(1), "sw_gate"))
}

# One gate list per sample, each ordered parents first and named by gate id.
.sw_gq_per_sample <- function(gates, samples) {
  if (.sw_gq_is_gate_list(gates)) {
    out <- rep(list(gates), length(samples))
    names(out) <- samples
  } else if (is.list(gates) && !is.null(names(gates)) &&
    all(vapply(gates, .sw_gq_is_gate_list, logical(1)))) {
    absent <- setdiff(samples, names(gates))
    if (length(absent) > 0) {
      stop("`gates` lacks entries for samples: ", paste(absent, collapse = ", "), call. = FALSE)
    }
    out <- gates[samples]
  } else {
    stop(
      "`gates` must be a list of sw_gate, or a named list of such lists (one per sample).",
      call. = FALSE
    )
  }
  lapply(out, .sw_gq_order)
}

# Order gates so each parent comes before its children; name the list by id.
.sw_gq_order <- function(gl) {
  if (!.sw_gq_is_gate_list(gl)) {
    stop("`gates` must be a non-empty list of sw_gate objects.", call. = FALSE)
  }
  ids <- vapply(gl, function(g) g$id, character(1))
  if (anyDuplicated(ids)) stop("Gate ids must be unique.", call. = FALSE)
  names(gl) <- ids
  parents <- vapply(gl, function(g) g$parent, character(1))
  bad <- parents[!is.na(parents) & !parents %in% ids]
  if (length(bad) > 0) {
    stop("Parent gate not found in `gates`: ", paste(unique(bad), collapse = ", "), call. = FALSE)
  }
  done <- character()
  ordered <- character()
  while (length(ordered) < length(ids)) {
    ready <- ids[!ids %in% ordered & (is.na(parents) | parents %in% ordered)]
    if (length(ready) == 0) stop("Gate parents form a cycle.", call. = FALSE)
    ordered <- c(ordered, ready)
  }
  gl[ordered]
}

# Parent and inside masks for each gate (ordered list) on one event matrix.
.sw_gq_masks <- function(gl, m) {
  out <- list()
  for (id in names(gl)) {
    g <- gl[[id]]
    absent <- setdiff(g$channels, colnames(m))
    if (length(absent) > 0) {
      stop("Data lack channels of gate '", id, "': ", paste(absent, collapse = ", "), call. = FALSE)
    }
    pm <- if (is.na(g$parent)) rep(TRUE, nrow(m)) else out[[g$parent]]$inside
    inside <- pm & sw_gate_apply(g, m)
    out[[id]] <- list(parent = pm, inside = inside)
  }
  out
}

# A single positive number, or a named vector with a value per channel.
.sw_gq_chan_arg <- function(x, channels, name) {
  if (!is.numeric(x) || length(x) == 0 || anyNA(x) || any(x <= 0)) {
    stop("`", name, "` must be positive numbers.", call. = FALSE)
  }
  if (length(x) == 1 && is.null(names(x))) {
    return(stats::setNames(rep(unname(x), length(channels)), channels))
  }
  if (is.null(names(x)) || !all(channels %in% names(x))) {
    stop(
      "`", name, "` must be one number, or a named vector covering: ",
      paste(channels, collapse = ", "), call. = FALSE
    )
  }
  x[channels]
}

.sw_gq_width <- function(w, g) {
  v <- unname(w[g$channels])
  if (length(unique(v)) > 1) {
    stop("A polygon gate needs one edge width; give both channels the same value.", call. = FALSE)
  }
  v[1]
}

# Distance of each event to the gate boundary.
.sw_gq_edge_dist <- function(g, m) {
  if (g$type == "range") {
    x <- m[, g$channels]
    return(pmin(abs(x - g$vertices["min", 1]), abs(x - g$vertices["max", 1])))
  }
  px <- m[, g$channels[1]]
  py <- m[, g$channels[2]]
  vx <- g$vertices[, 1]
  vy <- g$vertices[, 2]
  n <- length(vx)
  d <- rep(Inf, length(px))
  for (i in seq_len(n)) {
    j <- if (i == n) 1L else i + 1L
    dx <- vx[j] - vx[i]
    dy <- vy[j] - vy[i]
    len2 <- dx^2 + dy^2
    t <- if (len2 > 0) ((px - vx[i]) * dx + (py - vy[i]) * dy) / len2 else 0
    t <- pmin(pmax(t, 0), 1)
    d <- pmin(d, sqrt((px - vx[i] - t * dx)^2 + (py - vy[i] - t * dy)^2))
  }
  d
}
