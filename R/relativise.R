# R/relativise.R
# Reference-population (CD3) relativisation: per-channel offset and gain from a
# reference control, applied to linear intensities. Written from the formula in
# Glehr et al. 2026, EBioMedicine 131:106427 (Box 1), doi 10.1016/j.ebiom.2026.106427.
# No code from the authors' packages was used. Spec: plans/2026-11-09_reference-relativisation.md

#' Split reference controls into negative and positive events by k-means
#'
#' @description
#' Step 1 of the relativisation workflow. For each control in `control_map`,
#' take the channel that the control calibrates, apply `split_transform`, and
#' split the events into two groups with one-dimensional k-means (k = 2). The
#' group with the larger centre is "positive". Glehr et al. (2026) used k-means
#' with k = 2 on lymphocyte-gated CD3 single-stain controls.
#'
#' @details
#' **Gate first.** Pass controls that are already gated (for example on
#' lymphocytes). This function never gates.
#'
#' **No defaults.** The scale k-means runs on (`split_transform`), the `seed`
#' and `nstart` all change which events are called positive, so each must be
#' passed. The paper does not say which scale its k-means used. The caller's
#' RNG state is restored.
#'
#' A manual gate is an equal alternative: build the same named list of logical
#' vectors yourself and pass it to [sw_relativise_fit()].
#'
#' @param controls A named list of numeric matrices, one per control file, with
#'   one column per channel. Values are compensated or unmixed and linear.
#' @param control_map A data frame with columns `unit`, `channel` and `control`:
#'   which control calibrates which channel for which unit (a sample or a batch).
#' @param split_transform A function applied to the channel values before
#'   k-means, for example `function(v) asinh(v / 6000)` or `identity`.
#' @param seed Whole number for k-means starts.
#' @param nstart Whole number of k-means random starts.
#'
#' @return A named list of logical vectors, one per control used in
#'   `control_map`, `TRUE` for reference-positive events.
#'
#' @seealso [sw_relativise_fit()]
#'
#' @examples
#' set.seed(1)
#' ctrl <- list(c1 = cbind(B1 = c(stats::rnorm(200, 100, 20), stats::rnorm(200, 5000, 500))))
#' map <- data.frame(unit = "s1", channel = "B1", control = "c1")
#' pos <- sw_relativise_split_kmeans(ctrl, map, identity, seed = 1, nstart = 5)
#' table(pos$c1)
#'
#' @export
sw_relativise_split_kmeans <- function(controls, control_map, split_transform, seed, nstart) {
  .sw_relativise_check_map(controls, control_map)
  if (!is.function(split_transform)) {
    stop("`split_transform` must be a function, for example `identity`.", call. = FALSE)
  }
  .sw_check_count(seed, "seed")
  .sw_check_count(nstart, "nstart")

  # One channel per control: the channel it calibrates. A control used for two
  # channels would need two splits, which this step does not choose for you.
  per_control <- split(control_map$channel, control_map$control)
  multi <- names(per_control)[lengths(lapply(per_control, unique)) > 1]
  if (length(multi) > 0) {
    stop(
      "These controls calibrate more than one channel, so one k-means split is ",
      "ambiguous; gate them yourself: ",
      paste(multi, collapse = ", "),
      call. = FALSE
    )
  }

  .sw_with_seed(seed, {
    out <- lapply(names(per_control), function(ctrl) {
      channel <- per_control[[ctrl]][[1]]
      v <- split_transform(controls[[ctrl]][, channel])
      if (length(v) != nrow(controls[[ctrl]]) || anyNA(v)) {
        stop(
          "`split_transform` must return one non-missing value per event (control ",
          ctrl,
          ").",
          call. = FALSE
        )
      }
      km <- stats::kmeans(matrix(v), centers = 2, nstart = nstart)
      km$cluster == which.max(km$centers[, 1])
    })
  })
  stats::setNames(out, names(per_control))
}

#' Fit per-channel offset and gain from reference controls
#'
#' @description
#' Step 2 of the relativisation workflow. For each row of `control_map` (one
#' unit and channel), take the calibrating control's values in that channel and
#' compute
#'
#' * `m`: median over reference-negative events (the offset);
#' * `M`: median over reference-positive events;
#' * `gain`: `M - m`.
#'
#' These are the per-channel numbers of Glehr et al. (2026), Box 1.
#'
#' @details
#' Stops when a unit-channel has fewer than `min_events` events on either side,
#' or when `M <= m` (no positive population above the negative one), because the
#' rescaling would then divide by zero or flip the axis.
#'
#' @inheritParams sw_relativise_split_kmeans
#' @param positive A named list of logical vectors, one per control, `TRUE` for
#'   reference-positive events. From [sw_relativise_split_kmeans()] or a gate.
#' @param min_events Whole number. Fewest events allowed on each side.
#'
#' @return A tibble of class `sw_relativise_model` with columns `unit`,
#'   `channel`, `control`, `n_negative`, `n_positive`, `m`, `M` and `gain`.
#'
#' @seealso [sw_relativise_diagnose()], [sw_relativise_apply()]
#'
#' @examples
#' set.seed(1)
#' ctrl <- list(c1 = cbind(B1 = c(stats::rnorm(200, 100, 20), stats::rnorm(200, 5000, 500))))
#' map <- data.frame(unit = "s1", channel = "B1", control = "c1")
#' pos <- list(c1 = rep(c(FALSE, TRUE), each = 200))
#' sw_relativise_fit(ctrl, map, pos, min_events = 50)
#'
#' @export
sw_relativise_fit <- function(controls, control_map, positive, min_events) {
  .sw_relativise_check_map(controls, control_map)
  .sw_check_count(min_events, "min_events")
  missing_pos <- setdiff(unique(control_map$control), names(positive))
  if (length(missing_pos) > 0) {
    stop("`positive` has no entry for: ", paste(missing_pos, collapse = ", "), call. = FALSE)
  }

  rows <- lapply(seq_len(nrow(control_map)), function(i) {
    ctrl <- control_map$control[[i]]
    channel <- control_map$channel[[i]]
    v <- controls[[ctrl]][, channel]
    pos <- positive[[ctrl]]
    if (!is.logical(pos) || length(pos) != length(v) || anyNA(pos)) {
      stop(
        "`positive$",
        ctrl,
        "` must be TRUE/FALSE with one value per event.",
        call. = FALSE
      )
    }
    n_pos <- sum(pos)
    n_neg <- sum(!pos)
    if (n_pos < min_events || n_neg < min_events) {
      stop(
        sprintf(
          "Unit %s, channel %s: %d negative and %d positive events; `min_events` is %d.",
          control_map$unit[[i]],
          channel,
          n_neg,
          n_pos,
          min_events
        ),
        call. = FALSE
      )
    }
    m <- stats::median(v[!pos])
    big_m <- stats::median(v[pos])
    if (!(big_m > m)) {
      stop(
        sprintf(
          "Unit %s, channel %s: positive median (%g) is not above negative median (%g).",
          control_map$unit[[i]],
          channel,
          big_m,
          m
        ),
        call. = FALSE
      )
    }
    tibble::tibble(
      unit = as.character(control_map$unit[[i]]),
      channel = channel,
      control = ctrl,
      n_negative = n_neg,
      n_positive = n_pos,
      m = m,
      M = big_m,
      gain = big_m - m
    )
  })
  model <- dplyr::bind_rows(rows)
  dup <- duplicated(model[c("unit", "channel")])
  if (any(dup)) {
    stop(
      "`control_map` has more than one control for a unit and channel: ",
      paste(model$unit[dup], model$channel[dup], sep = "/", collapse = ", "),
      call. = FALSE
    )
  }
  class(model) <- c("sw_relativise_model", class(model))
  model
}

#' Diagnose a relativisation model: offset and gain across units
#'
#' @description
#' Step 3 of the relativisation workflow. The fitted gain per unit and channel
#' is itself a batch-effect measure: under the method's assumption that the
#' reference marker is stable across donors, the gain changes only with the
#' instrument, settings, reagent or control. Numbers only; nothing is marked bad.
#'
#' * `offset`: `m` per unit and channel.
#' * `gain`: `M - m` per unit and channel.
#' * `gain_z`: robust z of `log(gain)` across units within a channel,
#'   `(v - median) / (1.4826 * MAD)`. `NA` when the MAD is 0.
#' * `gain_batch_ratio`: a batch's median gain divided by the median gain over
#'   all units, per channel. Only when `batch` is given.
#' * `gain_order_rho`: Spearman correlation between `log(gain)` and acquisition
#'   order across units, per channel. Only when `order` is given and at least 4
#'   units have an order.
#'
#' @param model A `sw_relativise_model` from [sw_relativise_fit()].
#' @param batch `NULL`, or a named vector of batch labels, one per unit.
#' @param order `NULL`, or a named numeric vector of acquisition order or time,
#'   one per unit.
#'
#' @return A tibble with columns `check`, `channel`, `unit` (the unit, the batch
#'   for `gain_batch_ratio`, or `NA` for `gain_order_rho`) and `value`.
#'
#' @seealso [sw_relativise_fit()], [sw_batch_scores()]
#'
#' @export
sw_relativise_diagnose <- function(model, batch, order) {
  if (!inherits(model, "sw_relativise_model")) {
    stop("`model` must come from sw_relativise_fit().", call. = FALSE)
  }
  units <- unique(model$unit)
  for (arg in c("batch", "order")) {
    v <- get(arg)
    if (!is.null(v) && (is.null(names(v)) || !all(units %in% names(v)))) {
      stop("`", arg, "` must be NULL or named with every unit.", call. = FALSE)
    }
    if (!is.null(v) && anyDuplicated(names(v))) {
      stop("`", arg, "` has duplicated names.", call. = FALSE)
    }
  }

  long <- function(check, unit, value) {
    tibble::tibble(check = check, channel = model$channel, unit = unit, value = value)
  }
  out <- list(long("offset", model$unit, model$m), long("gain", model$unit, model$gain))

  log_gain <- log(model$gain)
  gain_z <- stats::ave(log_gain, model$channel, FUN = function(v) {
    s <- 1.4826 * stats::mad(v, constant = 1)
    if (s == 0) rep(NA_real_, length(v)) else (v - stats::median(v)) / s
  })
  out <- c(out, list(long("gain_z", model$unit, gain_z)))

  by_channel <- split(model, model$channel)
  if (!is.null(batch) && anyNA(batch[units])) {
    stop("`batch` is NA for: ", paste(units[is.na(batch[units])], collapse = ", "), call. = FALSE)
  }
  if (!is.null(batch)) {
    ratio <- lapply(by_channel, function(d) {
      b <- as.character(batch[d$unit])
      overall <- stats::median(d$gain)
      r <- tapply(d$gain, b, stats::median) / overall
      tibble::tibble(
        check = "gain_batch_ratio",
        channel = d$channel[[1]],
        unit = names(r),
        value = as.numeric(r)
      )
    })
    out <- c(out, ratio)
  }
  if (!is.null(order)) {
    rho <- lapply(by_channel, function(d) {
      o <- order[d$unit]
      ok <- !is.na(o)
      value <- if (sum(ok) >= 4) {
        suppressWarnings(stats::cor(log(d$gain[ok]), o[ok], method = "spearman"))
      } else {
        NA_real_
      }
      tibble::tibble(
        check = "gain_order_rho",
        channel = d$channel[[1]],
        unit = NA_character_,
        value = value
      )
    })
    out <- c(out, rho)
  }
  dplyr::bind_rows(out)
}

#' Apply relativisation to linear sample data
#'
#' @description
#' Step 4 of the relativisation workflow. For each file, look up its unit, then
#' rescale each channel in `channels`:
#'
#' `X~ = low + (X - m) / (M - m) * (high - low)`
#'
#' With `low = 0` and `high = 1e4` this is the formula of Glehr et al. (2026),
#' Box 1. The reference-negative median maps to `low` and the
#' reference-positive median maps to `high`.
#'
#' @details
#' **Linear values only.** The formula removes an additive offset and a gain,
#' which only holds on linear (compensated or unmixed, untransformed) values.
#' `input_scale` must be `"linear"`; it is a statement by the caller, because
#' the scale cannot be read reliably from the numbers. Transform afterwards.
#'
#' Channels not in `channels` are returned unchanged. A file whose unit has no
#' model row for a requested channel is an error, not a silent skip.
#'
#' @param data A named list of numeric matrices, one per sample file, with one
#'   column per channel. Linear values.
#' @param model A `sw_relativise_model` from [sw_relativise_fit()].
#' @param unit A named vector giving the unit of each file in `data`.
#' @param channels Character vector of channels to relativise.
#' @param low,high Numbers. Targets for the negative and positive medians
#'   (the paper used 0 and 1e4).
#' @param input_scale Must be the string `"linear"`.
#'
#' @return A named list of matrices like `data`, with attribute
#'   `relativisation` holding `low`, `high` and the model rows used.
#'
#' @seealso [sw_relativise_fit()]
#'
#' @export
sw_relativise_apply <- function(data, model, unit, channels, low, high, input_scale) {
  if (!inherits(model, "sw_relativise_model")) {
    stop("`model` must come from sw_relativise_fit().", call. = FALSE)
  }
  if (!identical(input_scale, "linear")) {
    stop(
      "`input_scale` must be \"linear\": relativisation is defined on linear ",
      "values, before any transform (Glehr et al. 2026, Box 1).",
      call. = FALSE
    )
  }
  for (arg in c("low", "high")) {
    v <- get(arg)
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v)) {
      stop("`", arg, "` must be one finite number.", call. = FALSE)
    }
  }
  if (!(high > low)) {
    stop("`high` must be above `low`.", call. = FALSE)
  }
  if (is.null(names(data)) || is.null(names(unit)) || !all(names(data) %in% names(unit))) {
    stop("`unit` must be named with every file in `data`.", call. = FALSE)
  }

  used <- list()
  out <- lapply(names(data), function(f) {
    x <- data[[f]]
    absent <- setdiff(channels, colnames(x))
    if (length(absent) > 0) {
      stop("File ", f, " lacks channels: ", paste(absent, collapse = ", "), call. = FALSE)
    }
    u <- as.character(unit[[f]])
    rows <- model[model$unit == u & model$channel %in% channels, ]
    no_ref <- setdiff(channels, rows$channel)
    if (length(no_ref) > 0) {
      stop(
        "Unit ",
        u,
        " (file ",
        f,
        ") has no reference for: ",
        paste(no_ref, collapse = ", "),
        call. = FALSE
      )
    }
    for (i in seq_len(nrow(rows))) {
      ch <- rows$channel[[i]]
      x[, ch] <- low + (x[, ch] - rows$m[[i]]) / rows$gain[[i]] * (high - low)
    }
    used[[f]] <<- rows
    x
  })
  names(out) <- names(data)
  attr(out, "relativisation") <- list(
    low = low,
    high = high,
    rows = dplyr::bind_rows(used, .id = "file")
  )
  out
}

#' Build a control map from per-batch reference controls
#'
#' @description
#' Extension E1 of the relativisation workflow, for panels too large to carry a
#' CD3 control per channel and sample (for example 40-colour spectral panels).
#' Each batch has one control per channel: a single stain of one stable
#' reference marker in that channel's fluorochrome, stained and acquired with
#' the batch. The spectral unmixing single stains qualify when they all use the
#' same marker per channel in every batch. The returned map has `unit = batch`,
#' so [sw_relativise_apply()] then takes the batch of each file as its unit.
#'
#' @details
#' **What this changes.** In Glehr et al. (2026) the reference is per sample.
#' Here it is per batch: it removes offset and gain between batches, not
#' between samples inside a batch. Different channels may use different
#' reference markers, but each channel must use the same marker in every batch,
#' otherwise the batches would be scaled to different biology. This function
#' stops when that fails, when a batch lacks a channel another batch has, or
#' when a batch has two controls for one channel.
#'
#' @param control_info A data frame with columns `batch`, `channel`, `control`
#'   (the name in the `controls` list) and `reference_marker`.
#'
#' @return A data frame with columns `unit`, `channel`, `control`, ready for
#'   [sw_relativise_split_kmeans()] and [sw_relativise_fit()].
#'
#' @seealso [sw_relativise_internal_map()], [sw_relativise_fit()]
#'
#' @examples
#' info <- data.frame(
#'   batch = rep(c("b1", "b2"), each = 2),
#'   channel = rep(c("BV421-A", "PE-A"), 2),
#'   control = c("b1_BV421", "b1_PE", "b2_BV421", "b2_PE"),
#'   reference_marker = rep(c("CD4", "CD8"), 2)
#' )
#' sw_relativise_batch_map(info)
#'
#' @export
sw_relativise_batch_map <- function(control_info) {
  need <- c("batch", "channel", "control", "reference_marker")
  if (!is.data.frame(control_info) || !all(need %in% names(control_info))) {
    stop(
      "`control_info` must be a data frame with columns ",
      paste(need, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  info <- control_info[need]
  key <- paste(info$batch, info$channel, sep = " / ")
  if (anyDuplicated(key)) {
    stop(
      "More than one control for: ",
      paste(unique(key[duplicated(key)]), collapse = ", "),
      call. = FALSE
    )
  }
  markers <- tapply(info$reference_marker, info$channel, \(v) length(unique(v)))
  mixed <- names(markers)[markers > 1]
  if (length(mixed) > 0) {
    stop(
      "These channels use different reference markers across batches: ",
      paste(mixed, collapse = ", "),
      call. = FALSE
    )
  }
  per_batch <- split(info$channel, info$batch)
  all_channels <- unique(info$channel)
  gaps <- vapply(per_batch, \(ch) paste(setdiff(all_channels, ch), collapse = ", "), "")
  gaps <- gaps[nzchar(gaps)]
  if (length(gaps) > 0) {
    stop(
      "Batches without a control for every channel: ",
      paste(names(gaps), gaps, sep = " lacks ", collapse = "; "),
      call. = FALSE
    )
  }
  data.frame(
    unit = as.character(info$batch),
    channel = info$channel,
    control = info$control,
    stringsAsFactors = FALSE
  )
}

#' Build a control map that uses each sample as its own reference
#'
#' @description
#' Extension E2 of the relativisation workflow. The fully stained sample is
#' its own control for the reference channel only (for example the CD3
#' channel): reference-negative and reference-positive events are split in
#' that channel, and only that channel is rescaled. No extra tube is needed.
#'
#' @details
#' **Only the reference channel is calibrated.** In the full stain, other
#' channels carry other markers, so their medians on reference-positive cells
#' are biology, not an instrument scale. Pass `channels = channel` to
#' [sw_relativise_apply()]. The negative median also depends on which cells are
#' reference-negative in that sample (for example more monocytes), which the
#' single-stain control in the paper avoids.
#'
#' @param data A named list of per-file matrices (the full stains).
#' @param channel One channel name: the reference channel.
#'
#' @return A data frame with columns `unit`, `channel`, `control`, where unit
#'   and control are the file names. Pass `data` itself as `controls`.
#'
#' @seealso [sw_relativise_batch_map()]
#'
#' @export
sw_relativise_internal_map <- function(data, channel) {
  if (!is.list(data) || is.null(names(data))) {
    stop("`data` must be a named list of matrices.", call. = FALSE)
  }
  if (!is.character(channel) || length(channel) != 1) {
    stop("`channel` must be one channel name.", call. = FALSE)
  }
  absent <- names(data)[!vapply(data, \(m) channel %in% colnames(m), logical(1))]
  if (length(absent) > 0) {
    stop("Files without channel ", channel, ": ", paste(absent, collapse = ", "), call. = FALSE)
  }
  data.frame(unit = names(data), channel = channel, control = names(data), stringsAsFactors = FALSE)
}

# Shared input check for the control list and control map.
.sw_relativise_check_map <- function(controls, control_map) {
  if (!is.list(controls) || is.null(names(controls))) {
    stop("`controls` must be a named list of matrices.", call. = FALSE)
  }
  need <- c("unit", "channel", "control")
  if (!is.data.frame(control_map) || !all(need %in% names(control_map))) {
    stop("`control_map` must be a data frame with columns unit, channel, control.", call. = FALSE)
  }
  absent <- setdiff(control_map$control, names(controls))
  if (length(absent) > 0) {
    stop("`controls` has no element: ", paste(unique(absent), collapse = ", "), call. = FALSE)
  }
  for (i in seq_len(nrow(control_map))) {
    if (!control_map$channel[[i]] %in% colnames(controls[[control_map$control[[i]]]])) {
      stop(
        "Control ",
        control_map$control[[i]],
        " has no column ",
        control_map$channel[[i]],
        ".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}
