# R/import_integrity.R
# FCS integrity audit: is each file readable, does $TOT match the rows read, how
# long did the run last, which channels does it have.
# Clean rewrite from the plan description (plans/2026-11-13_import-step.md, I2).
# No upstream code is used. The source idea for the time-span part was not
# opened in the survey, so that part follows the plan text only.

#' Audit the integrity of FCS files
#'
#' @description
#' For each file, report four items: whether it can be read, the difference
#' between the header `$TOT` and the number of rows actually read, the time span
#' of the run, and the channel set. The result is a long table with one row per
#' file and item. Nothing is dropped and no file gets a pass or fail mark or an
#' overall score.
#'
#' The items are:
#'
#' * `readable`: `value` is 1 when the file could be read, 0 when it could not.
#'   When it could not, `reason` holds the error text. With `header_only = TRUE`
#'   no event is read, so the row is `not_assessed` (`value` `NA`).
#' * `tot_vs_rows`: `value` is header `$TOT` minus the rows read. `status` is
#'   `match` when the absolute difference is at most `tot_tolerance`, and
#'   `differs` otherwise. `$TOT` is never changed, in memory or on disk.
#' * `time_span`: `value` is the length of the run, from `time_keyword`.
#' * `channel_set`: `value` is the number of channels. `detail` lists each
#'   `$PnN` name with its `$PnS` marker.
#'
#' @details
#' **No defaults.** `tot_tolerance`, `time_keyword` and `header_only` have no
#' default, because each changes which files look different.
#'
#' **Time span.** `time_keyword` takes one of two forms. A single string is a
#' channel name (a `$PnN` value, such as `"Time"`). The span is the maximum minus
#' the minimum of its finite values, in the units of that channel. The count of
#' non-finite values is stated in `detail`, and `$TIMESTEP` is quoted there
#' when the header has it; the value is not rescaled. Two strings are a
#' begin and an end header keyword, such as `c("$BTIM", "$ETIM")`. The span is
#' the difference in seconds. Times are read as `HH:MM:SS`, with an optional
#' `.cc` (hundredths) or `:tt` (sixtieths) tail. An end earlier than the begin
#' is `not_assessed`, because no day is added silently. The source idea for the
#' time-span part was not opened in the survey, so this part follows the plan
#' text only.
#'
#' **Events.** Unless `header_only = TRUE`, every event of every file is read
#' (one file at a time) through [sw_io_read_fcs()], which sets
#' `truncate_max_range = FALSE`. With `header_only = TRUE` no event is read.
#' The `$TOT` check then has nothing to compare, so its row is `not_assessed`
#' with a reason that says `header_only = TRUE`. A channel-name `time_keyword`
#' also needs events and is `not_assessed` in the same way. The result carries
#' attribute `header_only`.
#'
#' **A wrong `$TOT`.** flowCore refuses to return events from a file whose
#' `$TOT` disagrees with its data section, and its error states the rows it
#' counted. The audit takes that count as the rows read, so `tot_vs_rows` is
#' still filled and `readable` says why no events came back. A channel
#' `time_keyword` is then `not_assessed` with that reason. The count comes from
#' the wording of the flowCore error; if that wording changes, the row becomes
#' `not_assessed` with the error text.
#'
#' **Unreadable files.** A file that cannot be read gives `not_assessed` rows
#' for every item, with the error text in `reason`. A file whose header reads
#' but whose events do not gets `not_assessed` for `readable`, `tot_vs_rows`
#' and a channel `time_keyword`; its header items are still filled.
#'
#' **Provenance.** Attribute `provenance` holds the package version (from
#' `.sw_pkg_info()`, the helper behind the unmix provenance block), the flowCore
#' version, R version, platform, time, `seed` (`NA`: the audit draws no random
#' numbers), a summary of `sessionInfo()` and the arguments.
#'
#' @param paths Character vector of FCS file paths. File names (without the
#'   folder) must be unique, because they label the rows. Required.
#' @param tot_tolerance Single number, 0 or more. The largest absolute
#'   difference between `$TOT` and the rows read that still counts as `match`.
#'   Required.
#' @param time_keyword One channel name, or two header keywords (begin and end)
#'   that give the run time. See Details. Required.
#' @param header_only `TRUE` or `FALSE`. `TRUE` reads headers only and skips the
#'   `$TOT` check; the output says so. Required.
#'
#' @return A tibble of class `sw_fcs_integrity_audit` with columns `file`,
#'   `item` (`readable`, `tot_vs_rows`, `time_span`, `channel_set`), `value`
#'   (numeric, `NA` when not assessed), `detail` (text on how the value arose),
#'   `status` (`assessed`, `match`, `differs` or `not_assessed`), `reason`
#'   (why a row is `not_assessed`, else `NA`) and `path`. Attributes
#'   `header_only` and `provenance`.
#'
#' @seealso [sw_fcs_header_audit()], [sw_fcs_signal_audit()]
#'
#' @examples
#' if (requireNamespace("flowCore", quietly = TRUE)) {
#'   d <- tempfile("integrity")
#'   dir.create(d)
#'   m <- cbind(Time = 1:20, "B1-A" = rnorm(20))
#'   p <- file.path(d, "a.fcs")
#'   flowCore::write.FCS(flowCore::flowFrame(m), p)
#'   # Values below are illustrative; choose your own.
#'   sw_fcs_integrity_audit(
#'     p,
#'     tot_tolerance = 0,
#'     time_keyword = "Time",
#'     header_only = FALSE
#'   )
#' }
#'
#' @export
sw_fcs_integrity_audit <- function(
  paths,
  tot_tolerance,
  time_keyword,
  header_only
) {
  .sw_stop_missing(c(
    paths = missing(paths),
    tot_tolerance = missing(tot_tolerance),
    time_keyword = missing(time_keyword),
    header_only = missing(header_only)
  ))
  if (!is.character(paths) || length(paths) == 0 || anyNA(paths)) {
    stop("`paths` must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(basename(paths))) {
    stop("File names in `paths` must be unique.", call. = FALSE)
  }
  if (
    !is.numeric(tot_tolerance) ||
      length(tot_tolerance) != 1 ||
      is.na(tot_tolerance) ||
      tot_tolerance < 0
  ) {
    stop("`tot_tolerance` must be one number, 0 or more.", call. = FALSE)
  }
  if (
    !is.character(time_keyword) ||
      !length(time_keyword) %in% 1:2 ||
      anyNA(time_keyword) ||
      !all(nzchar(time_keyword))
  ) {
    stop(
      "`time_keyword` must be one channel name or two header keywords ",
      "(begin, end).",
      call. = FALSE
    )
  }
  if (
    !is.logical(header_only) || length(header_only) != 1 || is.na(header_only)
  ) {
    stop("`header_only` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  out <- lapply(
    paths,
    .sw_int_file,
    tot_tolerance,
    time_keyword,
    header_only
  ) |>
    dplyr::bind_rows()
  out <- tibble::as_tibble(out)
  attr(out, "header_only") <- header_only
  attr(out, "provenance") <- .sw_int_provenance(
    list(
      tot_tolerance = tot_tolerance,
      time_keyword = time_keyword,
      header_only = header_only
    )
  )
  class(out) <- c("sw_fcs_integrity_audit", class(out))
  out
}

# ---- Internal helpers --------------------------------------------------------

.sw_int_items <- c("readable", "tot_vs_rows", "time_span", "channel_set")

# One row for one file and item.
.sw_int_row <- function(
  file,
  path,
  item,
  value = NA_real_,
  detail = NA_character_,
  status,
  reason = NA_character_
) {
  tibble::tibble(
    file = file,
    item = item,
    value = as.numeric(value),
    detail = detail,
    status = status,
    reason = reason,
    path = path
  )
}

# The four rows of one file. Every failure becomes a not_assessed row.
.sw_int_file <- function(path, tot_tolerance, time_keyword, header_only) {
  file <- basename(path)
  row <- function(item, ...) .sw_int_row(file, path, item, ...)
  skip <- function(item, why) row(item, status = "not_assessed", reason = why)

  hdr <- tryCatch(
    suppressWarnings(unlist(.sw_read_fcs_header(path)[[1]])),
    error = \(e) e
  )
  if (inherits(hdr, "error")) {
    why <- .sw_fcs_read_reason(path, hdr)
    return(dplyr::bind_rows(
      row("readable", value = 0, status = "not_assessed", reason = why),
      skip("tot_vs_rows", why),
      skip("time_span", why),
      skip("channel_set", why)
    ))
  }

  # Events: read once, only what the items need.
  events <- NULL
  read_err <- NULL
  n_rows <- NA_integer_
  # flowCore refuses to return events when $TOT disagrees with the data
  # section, but its error states the rows it counted. Use that count.
  tot_refusal <- FALSE
  if (!header_only) {
    ff <- suppressWarnings(tryCatch(sw_io_read_fcs(path), error = \(e) e))
    if (inherits(ff, "error")) {
      read_err <- .sw_int_msg(ff)
      counted <- regmatches(
        read_err,
        regexec(
          "data section \\(([0-9]+)\\) is not consistent with keyword '\\$TOT'",
          read_err
        )
      )[[1]]
      if (length(counted) == 2) {
        n_rows <- as.integer(counted[[2]])
        tot_refusal <- TRUE
      }
    } else {
      events <- flowCore::exprs(ff)
      n_rows <- nrow(events)
      rm(ff)
    }
  }
  no_events <- if (header_only) {
    "header_only = TRUE: no events were read"
  } else if (tot_refusal) {
    paste0(
      "flowCore returns no events because $TOT disagrees with the data ",
      "section (",
      n_rows,
      " rows counted)"
    )
  } else {
    paste("events unreadable:", read_err)
  }

  readable <- if (header_only) {
    # Only the header was parsed; nothing says the data section can be read.
    skip("readable", "header_only = TRUE: events not read, header parsed only")
  } else if (tot_refusal) {
    row(
      "readable",
      value = 1,
      detail = paste0(
        "data section holds ",
        n_rows,
        " rows; flowCore returns no events because $TOT disagrees"
      ),
      status = "assessed"
    )
  } else if (is.null(events)) {
    row("readable", value = 0, status = "not_assessed", reason = read_err)
  } else {
    row(
      "readable",
      value = 1,
      detail = sprintf("%d rows read", n_rows),
      status = "assessed"
    )
  }

  tot_row <- if (is.na(n_rows)) {
    skip("tot_vs_rows", no_events)
  } else if (
    !"$TOT" %in% names(hdr) ||
      is.na(suppressWarnings(as.numeric(hdr[["$TOT"]])))
  ) {
    skip("tot_vs_rows", "header has no readable $TOT")
  } else {
    tot <- as.numeric(hdr[["$TOT"]])
    diff <- tot - n_rows
    row(
      "tot_vs_rows",
      value = diff,
      detail = sprintf(
        "$TOT = %s; rows read = %d",
        hdr[["$TOT"]],
        n_rows
      ),
      status = if (abs(diff) <= tot_tolerance) "match" else "differs"
    )
  }

  span_row <- if (length(time_keyword) == 2) {
    .sw_int_span_keywords(hdr, time_keyword, row, skip)
  } else {
    .sw_int_span_channel(hdr, events, time_keyword, no_events, row, skip)
  }

  dplyr::bind_rows(
    readable,
    tot_row,
    span_row,
    .sw_int_channels(hdr, row, skip)
  )
}

# Time span from a begin and an end header keyword, in seconds.
.sw_int_span_keywords <- function(hdr, keys, row, skip) {
  missing_key <- keys[!keys %in% names(hdr)]
  if (length(missing_key) > 0) {
    return(skip(
      "time_span",
      paste("header lacks keyword(s):", paste(missing_key, collapse = ", "))
    ))
  }
  secs <- vapply(keys, \(k) .sw_int_clock(hdr[[k]]), numeric(1))
  if (anyNA(secs)) {
    bad <- keys[is.na(secs)]
    return(skip(
      "time_span",
      sprintf(
        "cannot read time from %s (%s)",
        paste(bad, collapse = ", "),
        paste(hdr[bad], collapse = ", ")
      )
    ))
  }
  if (secs[[2]] < secs[[1]]) {
    return(skip(
      "time_span",
      sprintf(
        "%s is earlier than %s; the run may cross midnight and no day is added",
        keys[[2]],
        keys[[1]]
      )
    ))
  }
  row(
    "time_span",
    value = secs[[2]] - secs[[1]],
    detail = sprintf("seconds, %s to %s", keys[[1]], keys[[2]]),
    status = "assessed"
  )
}

# Seconds since midnight from HH:MM:SS[.cc | :tt]; NA when the text does not fit.
.sw_int_clock <- function(x) {
  m <- regmatches(
    trimws(x),
    regexec(
      "^([0-9]{1,2}):([0-9]{2}):([0-9]{2})(?:([.:])([0-9]+))?$",
      trimws(x),
      perl = TRUE
    )
  )[[1]]
  if (length(m) == 0) {
    return(NA_real_)
  }
  frac <- if (nzchar(m[[6]])) {
    as.numeric(m[[6]]) / (if (m[[5]] == ".") 10^nchar(m[[6]]) else 60)
  } else {
    0
  }
  as.numeric(m[[2]]) *
    3600 +
    as.numeric(m[[3]]) * 60 +
    as.numeric(m[[4]]) +
    frac
}

# Time span from one channel of the events, in that channel's units.
.sw_int_span_channel <- function(hdr, events, chan, no_events, row, skip) {
  names_n <- unname(hdr[grep("^\\$P[0-9]+N$", names(hdr))])
  if (!chan %in% names_n) {
    return(skip(
      "time_span",
      sprintf("channel '%s' is not in the header", chan)
    ))
  }
  if (is.null(events)) {
    return(skip("time_span", no_events))
  }
  if (!chan %in% colnames(events)) {
    return(skip(
      "time_span",
      sprintf("channel '%s' is not in the events", chan)
    ))
  }
  v <- events[, chan]
  ok <- is.finite(v)
  if (!any(ok)) {
    return(skip("time_span", sprintf("channel '%s' has no finite value", chan)))
  }
  step <- if ("$TIMESTEP" %in% names(hdr)) {
    sprintf("; $TIMESTEP = %s", hdr[["$TIMESTEP"]])
  } else {
    ""
  }
  row(
    "time_span",
    value = max(v[ok]) - min(v[ok]),
    detail = sprintf(
      "units of channel '%s'%s; %d non-finite value(s)",
      chan,
      step,
      sum(!ok)
    ),
    status = "assessed"
  )
}

# Channel set from the header: $PnN with $PnS.
.sw_int_channels <- function(hdr, row, skip) {
  pn <- grep("^\\$P[0-9]+N$", names(hdr), value = TRUE)
  if (length(pn) == 0) {
    return(skip("channel_set", "header has no $PnN keyword"))
  }
  idx <- sub("^\\$P([0-9]+)N$", "\\1", pn)
  idx <- idx[order(as.integer(idx))]
  chan <- unname(hdr[paste0("$P", idx, "N")])
  marker <- unname(hdr[paste0("$P", idx, "S")])
  marker <- ifelse(is.na(marker) | !nzchar(trimws(marker)), "", trimws(marker))
  row(
    "channel_set",
    value = length(chan),
    detail = paste0(
      chan,
      ifelse(nzchar(marker), paste0(" (", marker, ")"), ""),
      collapse = "; "
    ),
    status = "assessed"
  )
}

# Error text on one line, without the call.
.sw_int_msg <- function(e) {
  gsub("\\s+", " ", trimws(conditionMessage(e)))
}

# Short one-line reason why a header read failed, shared by the import step.
# Decided from the file first: 'file not found', 'not an FCS file'. Otherwise
# 'header unreadable: ' plus the first line of the read error `e`, without the
# call prefix and without the path.
.sw_fcs_read_reason <- function(path, e = NULL) {
  if (!file.exists(path) || dir.exists(path)) {
    return("file not found")
  }
  magic <- tryCatch(
    readBin(path, what = "raw", n = 3L),
    error = function(e) raw(0),
    warning = function(w) raw(0)
  )
  if (!identical(magic, charToRaw("FCS"))) {
    return("not an FCS file")
  }
  msg <- if (is.null(e)) {
    ""
  } else {
    gsub(path, "", conditionMessage(e), fixed = TRUE)
  }
  lines <- trimws(sub("^Error in [^:]*:", "", trimws(strsplit(msg, "\n")[[1]])))
  lines <- lines[nzchar(lines)]
  paste0(
    "header unreadable: ",
    if (length(lines)) lines[[1]] else "unknown error"
  )
}

# Provenance block: the shared import block (R/import_util.R). The audit draws
# no random numbers, so the seed slot is NA.
.sw_int_provenance <- function(args) {
  .sw_import_provenance(args = args, seed = NA_integer_)
}
