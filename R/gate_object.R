# R/gate_object.R
# A gate stored as plain data: id, parent, one or two channels, and the vertices
# of a polygon (two channels) or the bounds of a range (one channel). Shared by
# gate QC, gate transfer, density pre-gates, plots and the import/export
# wrappers, so every gate can be saved, compared and moved the same way.

#' Make a polygon gate on two channels
#'
#' @param id Character, the gate (population) name.
#' @param channels Character of length 2, the x and y channel.
#' @param vertices Numeric matrix or data frame with two columns (x, y), one row
#'   per vertex, at least 3 rows. The polygon closes from the last vertex back
#'   to the first.
#' @param parent Character, the parent gate's id, or `NA` for a gate on all
#'   events.
#'
#' @return An object of class `sw_gate`.
#' @seealso [sw_gate_range()], [sw_gate_apply()]
#' @export
sw_gate_polygon <- function(id, channels, vertices, parent) {
  if (!is.character(channels) || length(channels) != 2) {
    stop("`channels` must be two channel names (x, y).", call. = FALSE)
  }
  v <- as.matrix(vertices)
  if (!is.numeric(v) || ncol(v) != 2 || nrow(v) < 3 || anyNA(v)) {
    stop("`vertices` must be a numeric two-column matrix with at least 3 rows.", call. = FALSE)
  }
  colnames(v) <- channels
  .sw_gate_new(id, parent, channels, "polygon", v)
}

#' Make a range gate on one channel
#'
#' @param id Character, the gate (population) name.
#' @param channel Character, one channel name.
#' @param min,max Numbers, the lower and upper bound (`-Inf` or `Inf` for an
#'   open side). Events with `min <= x < max` are inside.
#' @inheritParams sw_gate_polygon
#'
#' @return An object of class `sw_gate`.
#' @seealso [sw_gate_polygon()], [sw_gate_apply()]
#' @export
sw_gate_range <- function(id, channel, min, max, parent) {
  if (!is.character(channel) || length(channel) != 1) {
    stop("`channel` must be one channel name.", call. = FALSE)
  }
  if (!is.numeric(min) || !is.numeric(max) || length(min) != 1 || length(max) != 1 ||
    is.na(min) || is.na(max) || !(max > min)) {
    stop("`min` and `max` must be numbers with `max` above `min`.", call. = FALSE)
  }
  v <- matrix(c(min, max), ncol = 1, dimnames = list(c("min", "max"), channel))
  .sw_gate_new(id, parent, channel, "range", v)
}

#' Which events fall inside a gate
#'
#' Point-in-polygon by ray casting (polygon gates) or a bound check (range
#' gates). The parent gate is not applied here; pass events already inside the
#' parent, or use the logical result to combine gates.
#'
#' @param gate An `sw_gate`.
#' @param data A numeric matrix or data frame of events with the gate's
#'   channels as columns.
#'
#' @return A logical vector, one value per row of `data`. `NA` values in a gate
#'   channel give `FALSE`.
#' @export
sw_gate_apply <- function(gate, data) {
  if (!inherits(gate, "sw_gate")) {
    stop("`gate` must be an sw_gate.", call. = FALSE)
  }
  absent <- setdiff(gate$channels, colnames(data))
  if (length(absent) > 0) {
    stop("`data` lacks gate channels: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  if (gate$type == "range") {
    x <- data[, gate$channels]
    inside <- x >= gate$vertices["min", 1] & x < gate$vertices["max", 1]
  } else {
    inside <- .sw_point_in_polygon(
      data[, gate$channels[1]],
      data[, gate$channels[2]],
      gate$vertices[, 1],
      gate$vertices[, 2]
    )
  }
  inside & !is.na(inside)
}

#' @export
print.sw_gate <- function(x, ...) {
  cat(
    sprintf(
      "<sw_gate %s> %s on %s, parent %s, %d vertices\n",
      x$id,
      x$type,
      paste(x$channels, collapse = " x "),
      ifelse(is.na(x$parent), "none", x$parent),
      nrow(x$vertices)
    )
  )
  invisible(x)
}

.sw_gate_new <- function(id, parent, channels, type, vertices) {
  if (!is.character(id) || length(id) != 1 || is.na(id) || !nzchar(id)) {
    stop("`id` must be one non-empty name.", call. = FALSE)
  }
  if (length(parent) != 1 || !(is.na(parent) || is.character(parent))) {
    stop("`parent` must be one gate id or NA.", call. = FALSE)
  }
  structure(
    list(
      id = id,
      parent = as.character(parent),
      channels = channels,
      type = type,
      vertices = vertices
    ),
    class = "sw_gate"
  )
}

# Even-odd ray casting; vectorised over points, looped over polygon edges.
.sw_point_in_polygon <- function(px, py, vx, vy) {
  n <- length(vx)
  inside <- rep(FALSE, length(px))
  j <- n
  for (i in seq_len(n)) {
    crosses <- ((vy[i] > py) != (vy[j] > py)) &
      (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i]) + vx[i])
    inside <- xor(inside, crosses %in% TRUE)
    j <- i
  }
  inside
}
