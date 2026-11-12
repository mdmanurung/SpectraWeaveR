# R/batch_otd.R
# Joint (multivariate) batch distance by debiased Sinkhorn divergence, between
# batches in pairs or each batch against all others. Built on
# sw_sinkhorn_divergence() (R/otd.R); definitions as used by Glehr et al. 2026.

#' Optimal transport distance between batches
#'
#' @description
#' Where [sw_batch_scores()] looks at one marker at a time, this scores all
#' `markers` jointly. It draws events from each batch and computes the
#' debiased Sinkhorn divergence ([sw_sinkhorn_divergence()]) between batches.
#' With `mode = "pairwise"` every batch pair gets one row; with
#' `mode = "rest"` every batch is compared with the other batches pooled.
#'
#' @details
#' **Draws.** Each batch contributes `n_events` events, drawn in equal shares
#' from its files (`floor(n_events / n_files)` per file), so a file with many
#' events does not outweigh the others. Files are the replicates, not cells.
#' In `"rest"` mode the other batches each give `floor(n_events / (n_batches
#' - 1))` events, again in equal shares per file
#' (`floor(floor(n_events / (n_batches - 1)) / n_files)` from each file's
#' draw), so no file is left out of the pooled set. A file with fewer events
#' than its share is an error.
#'
#' **Scale is yours.** `data` must already be on the analysis scale (for
#' example arcsinh); `blur` is in those units. There is no pass or fail and no
#' roll-up: one row per batch or pair, ranked only within this call.
#'
#' **No defaults.** `n_events`, `blur`, `seed`, `max_iter`, `tol` and `mode`
#' change the numbers and are uncalibrated. The caller's RNG state is
#' restored.
#'
#' @param data A named list of numeric matrices, one per file, columns named.
#' @param batch A named vector of batch labels, one per file in `data`.
#' @param markers Character vector of columns to use jointly.
#' @param n_events Whole number, events drawn per batch (uncalibrated).
#' @param seed Whole number for the draws (uncalibrated).
#' @param mode `"pairwise"` or `"rest"` (uncalibrated).
#' @inheritParams sw_sinkhorn_divergence
#'
#' @return A tibble of class `sw_batch_otd` with columns `mode`, `batch` (the
#'   batch, or `"b1 ~ b2"` for a pair), `otd`, `converged`, `n_files` (files
#'   behind the row; for a pair, the smaller batch) and `n_events` (events in
#'   the first and second set, as `"n1 ~ n2"`).
#'
#' @seealso [sw_sinkhorn_divergence()], [sw_otd_scores()], [sw_batch_scores()]
#'
#' @examples
#' set.seed(1)
#' files <- lapply(1:6, function(i) {
#'   matrix(stats::rnorm(400) + (i > 4) * 2, ncol = 2, dimnames = list(NULL, c("A", "B")))
#' })
#' names(files) <- paste0("f", 1:6)
#' sw_batch_otd(
#'   files,
#'   batch = stats::setNames(c("b1", "b1", "b2", "b2", "b3", "b3"), names(files)),
#'   markers = c("A", "B"),
#'   n_events = 100, # illustrative
#'   blur = 1, # illustrative
#'   seed = 1,
#'   max_iter = 200, # illustrative
#'   tol = 1e-4, # illustrative
#'   mode = "rest"
#' )
#'
#' @export
sw_batch_otd <- function(
  data,
  batch,
  markers,
  n_events,
  blur,
  seed,
  max_iter,
  tol,
  mode
) {
  .sw_stop_missing(c(
    data = missing(data),
    batch = missing(batch),
    markers = missing(markers),
    n_events = missing(n_events),
    blur = missing(blur),
    seed = missing(seed),
    max_iter = missing(max_iter),
    tol = missing(tol),
    mode = missing(mode)
  ))
  if (
    !is.character(mode) || length(mode) != 1 || !mode %in% c("pairwise", "rest")
  ) {
    stop("`mode` must be \"pairwise\" or \"rest\".", call. = FALSE)
  }
  .sw_check_count(n_events, "n_events")
  .sw_check_count(seed, "seed")
  .sw_check_count(max_iter, "max_iter")
  if (!is.numeric(blur) || length(blur) != 1 || !is.finite(blur) || blur <= 0) {
    stop("`blur` must be one positive number.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be one positive number.", call. = FALSE)
  }
  files <- .sw_check_batch_files(data, batch, markers)
  batch <- as.character(batch[files])
  names(batch) <- files
  batches <- sort(unique(batch))
  if (length(batches) < 2) {
    stop("`batch` needs at least two batches.", call. = FALSE)
  }
  n_files <- table(batch)[batches]
  if (any(n_events < n_files)) {
    stop(
      "`n_events` must be at least the number of files in the largest batch.",
      call. = FALSE
    )
  }
  if (
    mode == "rest" &&
      floor(n_events / (length(batches) - 1)) < max(n_files)
  ) {
    stop(
      "`n_events` / (batches - 1) must be at least the number of files in the largest batch.",
      call. = FALSE
    )
  }

  # Draw each batch in equal shares per file, in sorted batch order.
  blocks <- .sw_with_seed(seed, {
    lapply(batches, \(b) {
      fs <- files[batch == b]
      per_file <- floor(n_events / length(fs))
      parts <- lapply(fs, \(f) {
        m <- as.matrix(data[[f]][, markers, drop = FALSE])
        if (nrow(m) < per_file) {
          stop(
            "File ",
            f,
            " has fewer events than its share (",
            per_file,
            ").",
            call. = FALSE
          )
        }
        m[sample.int(nrow(m), per_file), , drop = FALSE]
      })
      parts
    })
  })
  names(blocks) <- batches
  draws <- lapply(blocks, \(parts) do.call(rbind, parts))

  solve_pair <- \(x, y) {
    s <- sw_sinkhorn_divergence(
      x,
      y,
      blur = blur,
      max_iter = max_iter,
      tol = tol
    )
    list(otd = as.numeric(s), converged = isTRUE(attr(s, "converged")))
  }
  if (mode == "pairwise") {
    pairs <- utils::combn(batches, 2, simplify = FALSE)
    rows <- lapply(pairs, \(p) {
      r <- solve_pair(draws[[p[1]]], draws[[p[2]]])
      tibble::tibble(
        mode = mode,
        batch = paste(p, collapse = " ~ "),
        otd = r[["otd"]],
        converged = r[["converged"]],
        n_files = as.integer(min(n_files[[p[1]]], n_files[[p[2]]])),
        n_events = paste(nrow(draws[[p[1]]]), nrow(draws[[p[2]]]), sep = " ~ ")
      )
    })
  } else {
    per_other <- floor(n_events / (length(batches) - 1))
    # Rest set of b: an equal random share from every file block of each other
    # batch, so no file is left out because of its position in the draw.
    rests <- .sw_with_seed(seed, {
      lapply(batches, \(b) {
        shares <- lapply(blocks[setdiff(batches, b)], \(parts) {
          share <- floor(per_other / length(parts))
          lapply(parts, \(m) m[sample.int(nrow(m), share), , drop = FALSE])
        })
        do.call(rbind, unlist(shares, recursive = FALSE))
      })
    })
    names(rests) <- batches
    rows <- lapply(batches, \(b) {
      rest <- rests[[b]]
      r <- solve_pair(draws[[b]], rest)
      tibble::tibble(
        mode = mode,
        batch = b,
        otd = r[["otd"]],
        converged = r[["converged"]],
        n_files = as.integer(n_files[[b]]),
        n_events = paste(nrow(draws[[b]]), nrow(rest), sep = " ~ ")
      )
    })
  }
  out <- dplyr::bind_rows(rows)
  class(out) <- c("sw_batch_otd", class(out))
  out
}

# Check data, batch and markers; return the sorted file names.
.sw_check_batch_files <- function(data, batch, markers) {
  if (!is.list(data) || is.null(names(data)) || anyDuplicated(names(data))) {
    stop("`data` must be a named list of files.", call. = FALSE)
  }
  if (!is.character(markers) || length(markers) == 0) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  if (is.null(names(batch)) || !all(names(data) %in% names(batch))) {
    stop(
      "`batch` must be a named vector covering every file in `data`.",
      call. = FALSE
    )
  }
  files <- sort(names(data))
  if (anyNA(batch[files])) {
    stop("`batch` has missing labels.", call. = FALSE)
  }
  for (f in files) {
    absent <- setdiff(markers, colnames(data[[f]]))
    if (length(absent) > 0) {
      stop(
        "File ",
        f,
        " lacks markers: ",
        paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
    if (anyNA(data[[f]][, markers])) {
      stop("File ", f, " has missing values in `markers`.", call. = FALSE)
    }
  }
  files
}
