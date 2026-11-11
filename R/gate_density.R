# R/gate_density.R
# Per-sample density pre-gates. Each sample gets its own cut from its own
# density, the idea of flowDensity (Malek et al. 2015, Bioinformatics,
# Artistic-2.0). The code is our own. Nothing here transforms data: pass values
# on the scale the gate should use. Every threshold rule parameter is a
# required argument.

#' One cut point from the density of a numeric vector
#'
#' @description
#' Returns one number, or `NA` with a reason when the rule cannot be applied.
#' It never guesses. Methods:
#'
#' * `"valley"`: Gaussian kernel density (`stats::density`, bandwidth
#'   `bandwidth` in data units, 1024 grid points from `min(x) - 3 * bandwidth`
#'   to `max(x) + 3 * bandwidth`). A peak is a grid point higher than its left
#'   neighbour and not lower than its right neighbour, with density at least
#'   `min_peak_height` times the highest density. Peaks closer than
#'   `min_peak_distance` to a higher accepted peak are dropped (highest first).
#'   The cut is the grid point of lowest density between the two highest
#'   accepted peaks. Fewer than two peaks gives `NA`.
#' * `"percentile"`: `stats::quantile(x, prob, type = 7)`.
#' * `"mad"`: `median(x) + k * mad(x)` for `mad_side = "upper"`, or
#'   `median(x) - k * mad(x)` for `"lower"`, with `stats::mad` at its
#'   constant 1.4826. A zero MAD gives `NA`.
#'
#' Non-finite values are dropped before the rule runs.
#'
#' @param x Numeric vector on the scale to gate on.
#' @param method `"valley"`, `"percentile"` or `"mad"`.
#' @param bandwidth Positive number, kernel bandwidth in data units.
#'   Required for `"valley"`.
#' @param min_peak_height Number in (0, 1]. Required for `"valley"`.
#' @param min_peak_distance Number >= 0 in data units. Required for
#'   `"valley"`.
#' @param prob Number between 0 and 1. Required for `"percentile"`.
#' @param k Positive number of MADs. Required for `"mad"`.
#' @param mad_side `"upper"` or `"lower"`. Required for `"mad"`.
#'
#' @return A number. Attribute `reason` is `NA_character_` when the rule
#'   worked and a short text when the result is `NA`. For `"valley"`,
#'   attribute `peaks` holds the x positions of the accepted peaks.
#' @seealso [sw_gate_density_range()], [sw_gate_density_apply()]
#' @export
sw_gate_density_cut <- function(
  x,
  method,
  bandwidth,
  min_peak_height,
  min_peak_distance,
  prob,
  k,
  mad_side
) {
  .sw_stop_missing(c(x = missing(x), method = missing(method)))
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(method) || length(method) != 1 ||
    !method %in% c("valley", "percentile", "mad")) {
    stop("`method` must be \"valley\", \"percentile\" or \"mad\".", call. = FALSE)
  }
  x <- x[is.finite(x)]

  if (method == "valley") {
    .sw_stop_missing(c(
      bandwidth = missing(bandwidth),
      min_peak_height = missing(min_peak_height),
      min_peak_distance = missing(min_peak_distance)
    ))
    .sw_check_num(bandwidth, "bandwidth", lo = 0, open_lo = TRUE)
    .sw_check_num(min_peak_height, "min_peak_height", lo = 0, hi = 1, open_lo = TRUE)
    .sw_check_num(min_peak_distance, "min_peak_distance", lo = 0)
    return(.sw_density_valley(x, bandwidth, min_peak_height, min_peak_distance))
  }
  if (method == "percentile") {
    .sw_stop_missing(c(prob = missing(prob)))
    .sw_check_num(prob, "prob", lo = 0, hi = 1)
    if (length(x) < 1) return(.sw_cut_na("no finite values"))
    return(.sw_cut_ok(unname(stats::quantile(x, prob, type = 7))))
  }
  .sw_stop_missing(c(k = missing(k), mad_side = missing(mad_side)))
  .sw_check_num(k, "k", lo = 0, open_lo = TRUE)
  if (!is.character(mad_side) || length(mad_side) != 1 ||
    !mad_side %in% c("upper", "lower")) {
    stop("`mad_side` must be \"upper\" or \"lower\".", call. = FALSE)
  }
  if (length(x) < 1) return(.sw_cut_na("no finite values"))
  s <- stats::mad(x)
  if (!(s > 0)) return(.sw_cut_na("MAD is zero"))
  sgn <- if (mad_side == "upper") 1 else -1
  .sw_cut_ok(stats::median(x) + sgn * k * s)
}

#' Range gate from a per-sample density cut
#'
#' @description
#' Computes the cut of `data[, channel]` with [sw_gate_density_cut()] and
#' returns a range gate on one side of it. `"above"` keeps events with
#' `x >= cut`; `"below"` keeps events with `x < cut` (the bounds of
#' [sw_gate_range()]). Pass events already inside the parent gate.
#'
#' @param data Numeric matrix or data frame with `channel` as a column.
#' @param channel One channel name.
#' @param id,parent See [sw_gate_range()].
#' @param side `"above"` or `"below"` the cut. Required.
#' @param method,bandwidth,... Passed to [sw_gate_density_cut()]; each is
#'   required by the method it names (`method` always).
#'
#' @return An `sw_gate` with attribute `cut` (the number), or, when the cut is
#'   `NA`, an `sw_gate_failed` list with `id`, `parent` and `reason`.
#' @export
sw_gate_density_range <- function(data, channel, id, parent, side, method, bandwidth, ...) {
  .sw_stop_missing(c(
    data = missing(data), channel = missing(channel), id = missing(id),
    parent = missing(parent), side = missing(side), method = missing(method)
  ))
  if (!is.character(side) || length(side) != 1 || !side %in% c("above", "below")) {
    stop("`side` must be \"above\" or \"below\".", call. = FALSE)
  }
  x <- .sw_density_col(data, channel)
  cut <- sw_gate_density_cut(x, method = method, bandwidth = bandwidth, ...)
  if (is.na(cut)) {
    return(.sw_gate_failed(id, parent, attr(cut, "reason")))
  }
  g <- if (side == "above") {
    sw_gate_range(id, channel, min = cut[[1]], max = Inf, parent = parent)
  } else {
    sw_gate_range(id, channel, min = -Inf, max = cut[[1]], parent = parent)
  }
  attr(g, "cut") <- cut[[1]]
  g
}

#' Polygon gate around the main population on two channels
#'
#' @description
#' Fits a 2D Gaussian kernel density (`MASS::kde2d`, bandwidths `bandwidth`
#' in data units, 128 x 128 grid over the data range widened by 3 bandwidths).
#' The contour level is the density `t` at which the grid cells with
#' density >= `t` hold a share `level` of the total grid mass. The gate is the
#' closed contour line at `t` that contains the highest-density grid point
#' (the main population). Other regions above `t`, such as a far cluster, are
#' left out. Use it for lymphocytes on FSC/SSC or for debris exclusion.
#'
#' @param data Numeric matrix or data frame with both `channels` as columns.
#' @param channels Two channel names (x, y).
#' @param id,parent See [sw_gate_polygon()].
#' @param level Number in (0, 1), the share of density mass inside the
#'   contour. Required.
#' @param bandwidth Numeric of length 2, positive kernel bandwidths for x and
#'   y in data units. Required.
#'
#' @return An `sw_gate` (polygon), or an `sw_gate_failed` list when fewer than
#'   3 finite events exist or no closed contour holds the main peak.
#' @export
sw_gate_density_2d <- function(data, channels, id, parent, level, bandwidth) {
  .sw_stop_missing(c(
    data = missing(data), channels = missing(channels), id = missing(id),
    parent = missing(parent), level = missing(level), bandwidth = missing(bandwidth)
  ))
  .sw_check_num(level, "level", lo = 0, hi = 1, open_lo = TRUE, open_hi = TRUE)
  if (!is.numeric(bandwidth) || length(bandwidth) != 2 || anyNA(bandwidth) ||
    any(bandwidth <= 0)) {
    stop("`bandwidth` must be two positive numbers (x, y).", call. = FALSE)
  }
  if (!is.character(channels) || length(channels) != 2) {
    stop("`channels` must be two channel names (x, y).", call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is needed for sw_gate_density_2d().", call. = FALSE)
  }
  x <- .sw_density_col(data, channels[1])
  y <- .sw_density_col(data, channels[2])
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) {
    return(.sw_gate_failed(id, parent, "fewer than 3 finite events"))
  }
  ng <- 128L
  lims <- c(
    min(x) - 3 * bandwidth[1], max(x) + 3 * bandwidth[1],
    min(y) - 3 * bandwidth[2], max(y) + 3 * bandwidth[2]
  )
  kd <- MASS::kde2d(x, y, h = bandwidth, n = ng, lims = lims)
  z <- kd$z
  zs <- sort(as.vector(z), decreasing = TRUE)
  cum <- cumsum(zs) / sum(zs)
  t_level <- zs[which(cum >= level)[1]]
  top <- which(z == max(z), arr.ind = TRUE)[1, ]
  mode_xy <- c(kd$x[top[1]], kd$y[top[2]])
  lines <- grDevices::contourLines(kd$x, kd$y, z, levels = t_level)
  best <- NULL
  best_area <- -Inf
  for (ln in lines) {
    closed <- length(ln$x) > 3 && ln$x[1] == ln$x[length(ln$x)] &&
      ln$y[1] == ln$y[length(ln$y)]
    if (!closed) next
    vx <- ln$x[-length(ln$x)]
    vy <- ln$y[-length(ln$y)]
    if (!.sw_point_in_polygon(mode_xy[1], mode_xy[2], vx, vy)) next
    a <- abs(.sw_poly_area(vx, vy))
    if (a > best_area) {
      best <- cbind(vx, vy)
      best_area <- a
    }
  }
  if (is.null(best)) {
    return(.sw_gate_failed(id, parent, "no closed contour around the main peak"))
  }
  g <- sw_gate_polygon(id, channels, best, parent)
  attr(g, "contour_density") <- t_level
  g
}

#' Describe one density gate for [sw_gate_density_apply()]
#'
#' @param fun A gate function taking `data` first, such as
#'   [sw_gate_density_range()] or [sw_gate_density_2d()].
#' @param ... Named arguments for `fun`, including `id` and `parent`.
#'
#' @return A list with `fun` and `args`.
#' @export
sw_gate_density_spec <- function(fun, ...) {
  if (!is.function(fun)) stop("`fun` must be a function.", call. = FALSE)
  args <- list(...)
  if (length(args) > 0 && (is.null(names(args)) || any(!nzchar(names(args))))) {
    stop("All arguments after `fun` must be named.", call. = FALSE)
  }
  list(fun = fun, args = args)
}

#' Apply density gates to every sample, each with its own cuts
#'
#' @description
#' Runs the gate specs in order on each sample. A gate whose `parent` is not
#' `NA` is fitted on the events inside that parent gate (the parent's own
#' parent chain included), so a cut comes from the parent population's
#' density. A gate whose parent failed or is unknown fails too, with a reason.
#' Nothing is guessed: a failed gate is left out of the sample's gate list and
#' listed in the table.
#'
#' @param data Named list of numeric matrices or data frames, one per sample,
#'   on the scale to gate on.
#' @param specs Named list of specs from [sw_gate_density_spec()], in the order
#'   to run. Names are unique. Each spec's `id` argument must equal its name,
#'   and a parent must come before its children.
#'
#' @return A named list (one entry per sample) of named lists of `sw_gate`.
#'   Attribute `cuts` is a tibble with columns `sample`, `gate`, `parent`,
#'   `type`, `cut` (range gates, else `NA`), `n_parent`, `n_inside`, `status`
#'   (`"ok"` or `"failed"`) and `reason`.
#' @export
sw_gate_density_apply <- function(data, specs) {
  .sw_stop_missing(c(data = missing(data), specs = missing(specs)))
  if (!is.list(data) || length(data) < 1 || is.null(names(data)) ||
    anyDuplicated(names(data)) || any(!nzchar(names(data)))) {
    stop("`data` must be a named list of samples with unique names.", call. = FALSE)
  }
  if (!is.list(specs) || length(specs) < 1 || is.null(names(specs)) ||
    anyDuplicated(names(specs)) || any(!nzchar(names(specs)))) {
    stop("`specs` must be a named list of specs with unique names.", call. = FALSE)
  }
  for (nm in names(specs)) {
    s <- specs[[nm]]
    if (!is.list(s) || !is.function(s$fun) || !is.list(s$args)) {
      stop("Spec '", nm, "' must come from sw_gate_density_spec().", call. = FALSE)
    }
    if (!identical(s$args$id, nm)) {
      stop("Spec '", nm, "' needs `id` equal to its name.", call. = FALSE)
    }
  }
  rows <- list()
  out <- lapply(names(data), function(sn) {
    d <- data[[sn]]
    if (is.data.frame(d)) d <- as.matrix(d)
    if (!is.matrix(d) || !is.numeric(d)) {
      stop("Sample '", sn, "' must be a numeric matrix or data frame.", call. = FALSE)
    }
    masks <- list()
    gates <- list()
    for (gn in names(specs)) {
      s <- specs[[gn]]
      par <- s$args$parent
      par_na <- length(par) == 1 && is.na(par)
      row <- list(
        sample = sn, gate = gn, parent = if (par_na) NA_character_ else par,
        type = NA_character_, cut = NA_real_, n_parent = NA_integer_,
        n_inside = NA_integer_, status = "failed", reason = NA_character_
      )
      if (par_na) {
        pm <- rep(TRUE, nrow(d))
      } else if (!is.null(masks[[par]])) {
        pm <- masks[[par]]
      } else {
        row$reason <- paste0("parent gate '", par, "' not available")
        rows[[length(rows) + 1]] <<- row
        next
      }
      row$n_parent <- sum(pm)
      g <- do.call(s$fun, c(list(data = d[pm, , drop = FALSE]), s$args))
      if (inherits(g, "sw_gate_failed")) {
        row$reason <- g$reason
      } else {
        row$status <- "ok"
        row$type <- g$type
        if (!is.null(attr(g, "cut"))) row$cut <- attr(g, "cut")
        m <- pm
        m[pm] <- sw_gate_apply(g, d[pm, , drop = FALSE])
        masks[[gn]] <- m
        row$n_inside <- sum(m)
        gates[[gn]] <- g
      }
      rows[[length(rows) + 1]] <<- row
    }
    gates
  })
  names(out) <- names(data)
  tab <- do.call(rbind, lapply(rows, tibble::as_tibble))
  attr(out, "cuts") <- tab
  out
}

# ---- helpers ---------------------------------------------------------------

.sw_cut_ok <- function(v) structure(v, reason = NA_character_)
.sw_cut_na <- function(reason) structure(NA_real_, reason = reason)

.sw_gate_failed <- function(id, parent, reason) {
  structure(list(id = id, parent = parent, reason = reason), class = "sw_gate_failed")
}

.sw_check_num <- function(v, name, lo = -Inf, hi = Inf, open_lo = FALSE, open_hi = FALSE) {
  bad <- !is.numeric(v) || length(v) != 1 || is.na(v) ||
    (if (open_lo) v <= lo else v < lo) || (if (open_hi) v >= hi else v > hi)
  if (bad) {
    stop(
      "`", name, "` must be one number in ", if (open_lo) "(" else "[", lo, ", ", hi,
      if (open_hi) ")" else "]", ".",
      call. = FALSE
    )
  }
}

.sw_density_col <- function(data, channel) {
  if (!is.character(channel) || length(channel) != 1) {
    stop("`channel` must be one channel name.", call. = FALSE)
  }
  if (!channel %in% colnames(data)) {
    stop("`data` lacks channel: ", channel, call. = FALSE)
  }
  as.numeric(data[, channel])
}

.sw_poly_area <- function(x, y) {
  j <- c(length(x), seq_len(length(x) - 1))
  sum(x[j] * y - x * y[j]) / 2
}

.sw_density_valley <- function(x, bw, min_h, min_d) {
  if (length(x) < 2 || length(unique(x)) < 2) {
    return(.sw_cut_na("fewer than two distinct finite values"))
  }
  d <- stats::density(
    x, bw = bw, n = 1024,
    from = min(x) - 3 * bw, to = max(x) + 3 * bw
  )
  gx <- d$x
  gy <- d$y
  n <- length(gy)
  i <- 2:(n - 1)
  is_peak <- rep(FALSE, n)
  is_peak[i] <- gy[i] > gy[i - 1] & gy[i] >= gy[i + 1]
  cand <- which(is_peak & gy >= min_h * max(gy))
  cand <- cand[order(gy[cand], decreasing = TRUE)]
  kept <- integer()
  for (i in cand) {
    if (length(kept) == 0 || all(abs(gx[i] - gx[kept]) >= min_d)) kept <- c(kept, i)
  }
  if (length(kept) < 2) {
    return(.sw_cut_na(sprintf("%d peak(s) found; valley needs two", length(kept))))
  }
  two <- sort(kept[1:2])
  idx <- two[1]:two[2]
  cut <- gx[idx[which.min(gy[idx])]]
  structure(cut, reason = NA_character_, peaks = gx[sort(kept)])
}
