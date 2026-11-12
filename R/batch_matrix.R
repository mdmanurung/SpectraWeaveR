# R/batch_matrix.R
# Batch-by-marker view of any batch score or flag table. Ranks batches inside
# each check x marker column only; there is no roll-up over markers (Q14).

.sw_matrix_directions <- c("diag_key", "high", "low", "abs")

# Ranking key, larger = more unusual. `diag_key` follows the package's
# conventions for named checks (.sw_diag_key()).
.sw_matrix_key <- function(check, value, direction) {
  dir_row <- unname(direction[check])
  key <- rep(NA_real_, length(value))
  is_diag <- dir_row == "diag_key"
  key[is_diag] <- suppressWarnings(.sw_diag_key(check[is_diag], value[is_diag]))
  key[dir_row == "high"] <- value[dir_row == "high"]
  key[dir_row == "low"] <- -value[dir_row == "low"]
  key[dir_row == "abs"] <- abs(value[dir_row == "abs"])
  key[!is.finite(key)] <- NA_real_
  key
}

# Rank 1 = largest key inside one column; NA when the column has < 2 usable batches.
.sw_matrix_rank <- function(key) {
  if (sum(!is.na(key)) < 2) {
    return(rep(NA_integer_, length(key)))
  }
  as.integer(rank(-key, ties.method = "min", na.last = "keep"))
}

#' Batch-by-marker table of a batch score
#'
#' @description
#' Put any batch score or flag table into one long batch-by-marker table and
#' rank the batches inside each check x marker column. Nothing is rolled up
#' over markers, and no batch is marked bad: rank 1 is the batch with the most
#' unusual value of that marker for that check.
#'
#' @details
#' **Input.** `scores` is a data frame with columns `check`, `batch`, `marker`
#' and `value`, for example the result of [sw_batch_scores()] (also by
#' population, with a `population` column), a reference-consistency table, or
#' the result of [sw_flag_iqr()]. A table with `metric` and no `check` column
#' (the [sw_flag_iqr()] form) has `metric` used as `check`. Rows with no batch
#' (`median_range`, `eta2`) cannot sit in a batch column; they are left out and
#' listed in attribute `dropped`. A pair label such as `"b1 ~ b2"` stays as a
#' column of its own. Each check x marker x batch (x population) must occur
#' once. A [sw_flag_iqr()] table with several controls needs one control chosen
#' first, because averaging controls is a choice this function does not make.
#'
#' **Direction.** `direction` tells which end of a check is unusual. It is a
#' required named character vector with one entry per element of `check`:
#' * `"diag_key"`: the package convention of [sw_diagnosis_top()] for that
#'   check name (signed checks rank by size, `gain_batch_ratio` by
#'   `abs(log(value))`, `separation` from the lowest);
#' * `"high"`: the largest value is the most unusual;
#' * `"low"`: the smallest value is the most unusual;
#' * `"abs"`: the largest `abs(value)` is the most unusual.
#'
#' Ties share the lowest rank. A column with fewer than 2 batches, or an `NA`
#' value, gets rank `NA`.
#'
#' @param scores A data frame, see Details.
#' @param check Character vector of checks to keep (uncalibrated: it decides
#'   which columns of the matrix exist). Required.
#' @param direction Named character vector, one entry per `check`
#'   (uncalibrated: it decides which end of a check ranks first). Required.
#'
#' @return A tibble of class `sw_batch_matrix` with columns `check`, `batch`,
#'   `marker`, `population` (only when `scores` has it), `value` and `rank`.
#'   Attribute `wide` is the same values with one column per batch. Attribute
#'   `dropped` holds the rows left out for lack of a batch. Attribute `params`
#'   echoes `check` and `direction`.
#'
#' @seealso [sw_batch_scores()], [sw_flag_iqr()], [sw_plot_batch_matrix()]
#'
#' @examples
#' set.seed(1)
#' files <- lapply(1:6, function(i) {
#'   matrix(stats::rnorm(600) + (i > 3), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' })
#' names(files) <- paste0("f", 1:6)
#' sc <- sw_batch_scores(
#'   files,
#'   batch = stats::setNames(rep(c("b1", "b2"), each = 3), names(files)),
#'   markers = c("A", "B"),
#'   n_events = 100, # illustrative
#'   seed = 1,
#'   order = NULL,
#'   reference_group = NULL
#' )
#' sw_batch_matrix(
#'   sc,
#'   check = "spread_ratio",
#'   direction = c(spread_ratio = "diag_key")
#' )
#'
#' @export
sw_batch_matrix <- function(scores, check, direction) {
  .sw_stop_missing(c(
    scores = missing(scores),
    check = missing(check),
    direction = missing(direction)
  ))
  if (!is.data.frame(scores)) {
    stop("`scores` must be a data frame.", call. = FALSE)
  }
  scores <- tibble::as_tibble(scores)
  if (!"check" %in% names(scores) && "metric" %in% names(scores)) {
    scores[["check"]] <- scores[["metric"]]
  }
  need <- c("check", "batch", "marker", "value")
  if (!all(need %in% names(scores))) {
    stop(
      "`scores` needs columns: ",
      paste(setdiff(need, names(scores)), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.character(check) || length(check) == 0 || anyNA(check)) {
    stop("`check` must be a non-empty character vector.", call. = FALSE)
  }
  if (
    !is.character(direction) ||
      is.null(names(direction)) ||
      !setequal(names(direction), check) ||
      !all(direction %in% .sw_matrix_directions)
  ) {
    stop(
      "`direction` must be a named character vector with one entry per ",
      "`check`, each one of: ",
      paste(.sw_matrix_directions, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  absent <- setdiff(check, scores[["check"]])
  if (length(absent) > 0) {
    stop(
      "`check` not found in `scores`: ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  key_cols <- intersect(
    c("check", "marker", "population", "batch"),
    names(scores)
  )
  group_cols <- setdiff(key_cols, "batch")
  kept <- scores[scores[["check"]] %in% check, , drop = FALSE]
  dropped <- kept[is.na(kept[["batch"]]), , drop = FALSE]
  kept <- kept[!is.na(kept[["batch"]]), , drop = FALSE]
  kept[["batch"]] <- as.character(kept[["batch"]])
  if (anyDuplicated(kept[key_cols]) > 0) {
    repeated <- duplicated(kept[key_cols]) |
      duplicated(kept[key_cols], fromLast = TRUE)
    stop(
      "Each ",
      paste(key_cols, collapse = " x "),
      " must occur once; repeated in check: ",
      paste(unique(kept[["check"]][repeated]), collapse = ", "),
      ". Choose one control or reference first.",
      call. = FALSE
    )
  }

  key <- .sw_matrix_key(kept[["check"]], kept[["value"]], direction)
  kept[["rank"]] <- stats::ave(
    key,
    kept[group_cols],
    FUN = .sw_matrix_rank
  )
  kept[["rank"]] <- as.integer(kept[["rank"]])
  out <- kept[c(
    "check",
    "batch",
    "marker",
    intersect("population", names(kept)),
    "value",
    "rank"
  )]
  out <- out[
    order(
      factor(out[["check"]], levels = check),
      out[["marker"]],
      out[["batch"]]
    ),
  ]

  # Wide view: one row per check x marker (x population), one column per batch.
  row_id <- do.call(paste, c(out[group_cols], sep = "\r"))
  batches <- sort(unique(out[["batch"]]))
  wide <- out[!duplicated(row_id), group_cols]
  for (b in batches) {
    hit <- out[["batch"]] == b
    wide[[b]] <- out[["value"]][hit][match(
      do.call(paste, c(wide[group_cols], sep = "\r")),
      row_id[hit]
    )]
  }
  attr(out, "wide") <- tibble::as_tibble(wide)
  attr(out, "dropped") <- dropped
  attr(out, "params") <- list(check = check, direction = direction)
  class(out) <- c("sw_batch_matrix", class(out))
  out
}
