#' Plant a known shift in transformed cytometry data
#'
#' @description
#' Applies a known, per-target shift to cells that are already on the asinh
#' scale. It builds a synthetic reference with a planted batch effect, so a
#' diagnostic can be scored against the truth. Four families are available.
#'
#' * `offset_asinh`: `x + d` on the transformed scale.
#' * `gain_linear`: back to the linear scale with `sinh(x) * cofactor`,
#'   multiplied by `g`, then `asinh(. / cofactor)`.
#' * `stretch_median`: `median + s * (x - median)`; the median of the group is
#'   fixed and the spread scales by `s`.
#' * `two_mode`: half of the cells get `+d` and half get `-d`, chosen at random.
#'
#' The caller transforms. Nothing is transformed inside this function (DR37).
#' The shifted data and the truth table come back as separate list elements, so
#' the truth never travels in the data.
#'
#' @details
#' `x` is one row per cell. It needs a `batch` column, a `population` column and
#' one numeric column per marker. A target row `(batch, marker, population, size)`
#' shifts every cell with that batch and population in that marker column. A
#' target that matches no cell, a marker that is not a column and a repeated
#' target are errors, so no shift is dropped silently.
#'
#' `size` is `d` for `offset_asinh` and `two_mode`, `g` (above 0) for
#' `gain_linear` and `s` (above 0) for `stretch_median`.
#'
#' The truth table holds, for each target, the planted `size` and the realised
#' `mean_change` and `median_change` (shifted minus unshifted cells), plus
#' `n_cells`. Only `two_mode` uses random numbers; `seed` is still required so a
#' run is always reproducible. The caller's RNG state is restored.
#'
#' @param x A data frame, one row per cell, with columns `batch`, `population`
#'   and numeric marker columns already on the asinh scale.
#' @param targets A data frame with columns `batch`, `marker`, `population` and
#'   `size`. The shift size is an uncalibrated choice that sets how large the
#'   planted effect is.
#' @param family One of `"offset_asinh"`, `"gain_linear"`, `"stretch_median"`
#'   or `"two_mode"`. The family is an uncalibrated choice: it decides what kind
#'   of batch effect is planted.
#' @param cofactor Positive number. The asinh cofactor that made `x`, recorded
#'   in the output. It is uncalibrated here and must match the caller's
#'   transform.
#' @param seed Whole number for the `two_mode` split. The seed is an
#'   uncalibrated choice; the caller's RNG state is restored.
#'
#' @return A list with `data` (the shifted copy of `x`), `truth` (a tibble with
#'   one row per target) and `call` (the family, cofactor, seed and
#'   `transform = "asinh"`). For `two_mode`, `truth` carries the attribute
#'   `cell_sign` with the sign given to each shifted cell.
#' @seealso [sw_batch_scores()]
#' @examples
#' set.seed(1)
#' x <- data.frame(
#'   batch = rep(c("b1", "b2"), each = 50),
#'   population = "CD4T",
#'   CD3 = rnorm(100)
#' )
#' targets <- data.frame(batch = "b2", marker = "CD3", population = "CD4T", size = 0.5)
#' out <- sw_synref_shift(
#'   x,
#'   targets,
#'   family = "offset_asinh",
#'   cofactor = 6000,
#'   seed = 1
#' )
#' out[["truth"]]
#' @export
sw_synref_shift <- function(x, targets, family, cofactor, seed) {
  .sw_stop_missing(c(
    x = missing(x),
    targets = missing(targets),
    family = missing(family),
    cofactor = missing(cofactor),
    seed = missing(seed)
  ))
  families <- c("offset_asinh", "gain_linear", "stretch_median", "two_mode")
  if (!is.character(family) || length(family) != 1 || !family %in% families) {
    stop(
      "`family` must be one of: ",
      paste(families, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (
    !is.numeric(cofactor) ||
      length(cofactor) != 1 ||
      is.na(cofactor) ||
      cofactor <= 0
  ) {
    stop("`cofactor` must be one positive number.", call. = FALSE)
  }
  if (
    !is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed != round(seed)
  ) {
    stop("`seed` must be one whole number.", call. = FALSE)
  }
  if (!is.data.frame(x) || !all(c("batch", "population") %in% names(x))) {
    stop(
      "`x` must be a data frame with columns `batch` and `population`.",
      call. = FALSE
    )
  }
  need <- c("batch", "marker", "population", "size")
  if (
    !is.data.frame(targets) ||
      !all(need %in% names(targets)) ||
      nrow(targets) == 0
  ) {
    stop(
      "`targets` must be a non-empty data frame with columns ",
      paste(need, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  targets <- tibble::as_tibble(targets)
  if (!is.numeric(targets[["size"]]) || anyNA(targets[["size"]])) {
    stop("`targets$size` must be numeric with no NA.", call. = FALSE)
  }
  if (
    family %in%
      c("gain_linear", "stretch_median") &&
      any(targets[["size"]] <= 0)
  ) {
    stop("`size` must be above 0 for family '", family, "'.", call. = FALSE)
  }
  key <- paste(
    targets[["batch"]],
    targets[["marker"]],
    targets[["population"]],
    sep = "\r"
  )
  if (anyDuplicated(key)) {
    stop(
      "`targets` repeats a batch, marker and population combination.",
      call. = FALSE
    )
  }
  absent <- setdiff(unique(as.character(targets[["marker"]])), names(x))
  if (length(absent) > 0) {
    stop(
      "Marker(s) not in `x`: ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  bad_type <- !vapply(
    unique(as.character(targets[["marker"]])),
    \(m) is.numeric(x[[m]]),
    logical(1)
  )
  if (any(bad_type)) {
    stop(
      "Marker column(s) not numeric: ",
      paste(names(bad_type)[bad_type], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (anyNA(x[["batch"]]) || anyNA(x[["population"]])) {
    stop(
      "`x` has NA in `batch` or `population`; those cells cannot be matched.",
      call. = FALSE
    )
  }
  bad_value <- vapply(
    unique(as.character(targets[["marker"]])),
    \(m) any(!is.finite(x[[m]])),
    logical(1)
  )
  if (any(bad_value)) {
    stop(
      "Marker column(s) with NA or non-finite values: ",
      paste(names(bad_value)[bad_value], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  out <- x
  rows <- vector("list", nrow(targets))
  signs <- vector("list", nrow(targets))
  .sw_with_seed(seed, {
    for (i in seq_len(nrow(targets))) {
      b <- as.character(targets[["batch"]][[i]])
      m <- as.character(targets[["marker"]][[i]])
      p <- as.character(targets[["population"]][[i]])
      size <- targets[["size"]][[i]]
      idx <- which(
        as.character(x[["batch"]]) == b & as.character(x[["population"]]) == p
      )
      if (length(idx) == 0) {
        stop(
          sprintf("Target (batch %s, population %s) matches no cell.", b, p),
          call. = FALSE
        )
      }
      old <- x[[m]][idx]
      sign_i <- NULL
      new <- switch(
        family,
        offset_asinh = old + size,
        gain_linear = asinh(sinh(old) * cofactor * size / cofactor),
        stretch_median = stats::median(old) + size * (old - stats::median(old)),
        two_mode = {
          sign_i <- rep(-1, length(idx))
          sign_i[sample.int(length(idx), ceiling(length(idx) / 2))] <- 1
          old + sign_i * size
        }
      )
      out[[m]][idx] <- new
      rows[[i]] <- tibble::tibble(
        batch = b,
        marker = m,
        population = p,
        family = family,
        size = size,
        n_cells = length(idx),
        mean_change = mean(new) - mean(old),
        median_change = stats::median(new) - stats::median(old)
      )
      if (!is.null(sign_i)) {
        signs[[i]] <- tibble::tibble(
          batch = b,
          marker = m,
          population = p,
          row = idx,
          sign = sign_i
        )
      }
    }
  })
  truth <- dplyr::bind_rows(rows)
  if (family == "two_mode") {
    attr(truth, "cell_sign") <- dplyr::bind_rows(signs)
  }
  list(
    data = out,
    truth = truth,
    call = list(
      family = family,
      cofactor = cofactor,
      seed = seed,
      transform = "asinh"
    )
  )
}

# Undo a shift from the shifted data and the truth table (offset, gain, stretch,
# and two_mode when the cell_sign attribute is present). Used to test invertibility.
.sw_synref_invert <- function(data, truth, cofactor) {
  out <- data
  signs <- attr(truth, "cell_sign")
  for (i in seq_len(nrow(truth))) {
    idx <- which(
      as.character(data[["batch"]]) == truth[["batch"]][[i]] &
        as.character(data[["population"]]) == truth[["population"]][[i]]
    )
    m <- truth[["marker"]][[i]]
    size <- truth[["size"]][[i]]
    y <- data[[m]][idx]
    out[[m]][idx] <- switch(
      truth[["family"]][[i]],
      offset_asinh = y - size,
      gain_linear = asinh(sinh(y) * cofactor / size / cofactor),
      stretch_median = stats::median(y) + (y - stats::median(y)) / size,
      two_mode = {
        s <- signs[
          signs[["batch"]] == truth[["batch"]][[i]] &
            signs[["marker"]] == m &
            signs[["population"]] == truth[["population"]][[i]],
        ]
        y - s[["sign"]][match(idx, s[["row"]])] * size
      }
    )
  }
  out
}
