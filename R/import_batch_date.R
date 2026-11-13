# R/import_batch_date.R
# Batch and acquisition-date check: declared batch (and run) from a manifest
# crossed with the date, begin time and serial number in each FCS header.
# Clean rewrite from a written method description (idea: batch labels that
# disagree with when and on which instrument the files were acquired); no
# upstream code is used.

#' Cross declared batches with acquisition date, time and instrument serial
#'
#' @description
#' Read only the header of each FCS file. Parse the acquisition date and begin
#' time, and read the instrument serial number. Put them next to the batch (and
#' run) the manifest declares for the file. The result is a set of tables of
#' numbers: you decide whether a batch label and the acquisition dates agree.
#' No file is marked pass or fail and no file gets a summary score.
#'
#' The result is a list of class `sw_fcs_batch_date_check` with these parts:
#'
#' * `files`: one row per file in `paths`, in the order of `paths`.
#' * `crosstab`: file counts for every batch and date, for a tile map.
#' * `groups`: one row per declared batch, and one per declared run when the
#'   manifest has a `run` column, with the number of dates, serials and
#'   acquisition sessions and the gaps between files.
#' * `unmatched`: files with no manifest row, and manifest rows with no file.
#' * `provenance`: the shared import keys (versions, R version, platform, time,
#'   seed (none), session, arguments). Also set as attribute `provenance`.
#'
#' @details
#' **No defaults.** Every argument is required, because each one changes the
#' dates, the sessions or the files that count. Header keys, date layouts and
#' time zones differ between instruments and export tools, so the package does
#' not guess them.
#'
#' **Matching files to the manifest.** A file matches the manifest row whose
#' `sample_id` equals the name of the element of `paths`. If `paths` has no
#' names, the file name without its extension is used. A file with no manifest
#' row stays in `files` with `batch` and `run` set to `NA`, is counted in
#' `crosstab` under `NA`, and is listed in `unmatched`.
#'
#' **Keywords.** `date_keywords` and `time_keywords` are tried in the order you
#' give them. The first keyword that holds a non-empty value is used, and its
#' name is written to `date_keyword` or `time_keyword`. Names are matched
#' without regard to case. The raw text of the value is always kept.
#'
#' **Date and time.** The date text is read with `date_format`, a
#' [strptime()] format; a year below 1000 (for example a two-digit year read
#' with `%Y`) is treated as unparseable; text after the part the format uses is ignored, and
#' month names follow the R locale. The time text must start with
#' `hh:mm:ss`; anything after the seconds (fractions or ticks) is ignored. The
#' date and time are joined in `time_zone`.
#'
#' **Sessions.** Within one declared batch (and, separately, within one
#' declared run), files with a date and time are sorted by time. A gap of more
#' than `run_gap_hours` hours between neighbours starts a new session.
#' `n_gaps_over` counts those gaps and `n_sessions` counts the sessions.
#' Files without a full date and time take no part in the gaps. They are
#' counted in `n_not_assessed`.
#'
#' **Not assessed.** A file that cannot be read, or whose date is missing or
#' does not parse, gets `status` `not_assessed` and a `reason`. A file with a
#' valid date but no valid time gets `date_only`. All other files get
#' `assessed`. Nothing is dropped: these rows stay in `files`, in `crosstab`
#' (under an `NA` date) and in the counts of `groups`.
#'
#' @param paths Character vector of FCS file paths. Named elements give the
#'   `sample_id` that is matched to the manifest; without names the file name
#'   without extension is used. Required.
#' @param manifest A [sw_sample_manifest()] object. Its `batch` column, and its
#'   `run` column when present, are the declared labels. Required.
#' @param date_keywords Character vector of header keywords that hold the date,
#'   for example `"$DATE"`. Tried in order. Required.
#' @param time_keywords Character vector of header keywords that hold the begin
#'   time, for example `"$BTIM"`. Tried in order. Required.
#' @param serial_keyword One header keyword that holds the instrument serial
#'   number, for example `"$CYTSN"`. Required.
#' @param date_format A [strptime()] format for the date text, for example
#'   `"%d-%b-%Y"`. Required.
#' @param time_zone Time zone name from [OlsonNames()] used to build the
#'   date-time, for example `"UTC"`. Required.
#' @param run_gap_hours One positive number: the gap in hours between two
#'   neighbouring files in a batch (or run) that starts a new session.
#'   Required.
#'
#' @return A list of class `sw_fcs_batch_date_check`.
#'
#'   `files` is a tibble with `file` (name with extension, as in the other
#'   import tables), `path`, `sample_id` (the manifest join key), `in_manifest`, `batch`,
#'   `run`, `date_keyword`, `date_raw`, `date`, `date_status`, `time_keyword`,
#'   `time_raw`, `time_status`, `datetime`, `serial`, `gap_before_hours` (hours
#'   since the previous file of the same declared batch), `session` (number of
#'   the session inside the declared batch), `status` and `reason`. The
#'   `*_status` columns are `parsed`, `missing_keyword`, `unparseable` or
#'   `unreadable`.
#'
#'   `crosstab` is a tibble with `batch`, `date` and `n_files`, with a row for
#'   every combination, zero counts included.
#'
#'   `groups` is a tibble with `level` (`batch` or `run`), `batch`, `run`,
#'   `n_files`, `n_assessed` (files with a date and time), `n_not_assessed`,
#'   `n_dates`, `dates`, `n_serials`, `serials`, `first_datetime`,
#'   `last_datetime`, `span_hours`, `max_gap_hours`, `n_gaps_over`,
#'   `n_sessions` and `run_gap_hours`.
#'
#'   `unmatched` is a tibble with `side` (`file_without_manifest_row` or
#'   `manifest_row_without_file`), `sample_id`, `path`, `batch` and `run`.
#'
#'   `provenance` is a list.
#'
#' @seealso [sw_sample_manifest()], [sw_fcs_header_audit()]
#'
#' @examples
#' \dontrun{
#' man <- sw_sample_manifest(
#'   data.frame(id = c("a", "b"), b = c("b1", "b2")),
#'   sample_id = "id",
#'   batch = "b"
#' )
#' res <- sw_fcs_batch_date_check(
#'   paths = c(a = "fcs/a.fcs", b = "fcs/b.fcs"),
#'   manifest = man,
#'   date_keywords = "$DATE",
#'   time_keywords = "$BTIM",
#'   serial_keyword = "$CYTSN",
#'   date_format = "%d-%b-%Y", # illustrative, choose your own
#'   time_zone = "UTC",
#'   run_gap_hours = 12 # illustrative, choose your own
#' )
#' res$crosstab
#' res$groups
#' }
#'
#' @export
sw_fcs_batch_date_check <- function(
  paths,
  manifest,
  date_keywords,
  time_keywords,
  serial_keyword,
  date_format,
  time_zone,
  run_gap_hours
) {
  .sw_stop_missing(c(
    paths = missing(paths),
    manifest = missing(manifest),
    date_keywords = missing(date_keywords),
    time_keywords = missing(time_keywords),
    serial_keyword = missing(serial_keyword),
    date_format = missing(date_format),
    time_zone = missing(time_zone),
    run_gap_hours = missing(run_gap_hours)
  ))
  ids <- .sw_bd_check_args(
    paths,
    manifest,
    date_keywords,
    time_keywords,
    serial_keyword,
    date_format,
    time_zone,
    run_gap_hours
  )
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  # One header read per file; an unreadable file keeps its row with a reason.
  heads <- lapply(unname(paths), .sw_bd_read_header)
  date_hit <- lapply(heads, .sw_bd_first_value, keywords = date_keywords)
  time_hit <- lapply(heads, .sw_bd_first_value, keywords = time_keywords)
  serial_hit <- lapply(heads, .sw_bd_first_value, keywords = serial_keyword)
  pick <- function(hits, field) {
    vapply(hits, function(x) x[[field]], character(1))
  }
  unreadable <- vapply(heads, function(x) is.na(x$reason) == FALSE, logical(1))

  date_raw <- pick(date_hit, "value")
  time_raw <- pick(time_hit, "value")
  date <- .sw_bd_parse_date(date_raw, date_format)
  time_ok <- .sw_bd_parse_time(time_raw)
  stamp <- ifelse(
    !is.na(date) & !is.na(time_ok),
    paste(format(date, "%Y-%m-%d"), time_ok),
    NA_character_
  )
  datetime <- as.POSIXct(stamp, format = "%Y-%m-%d %H:%M:%S", tz = time_zone)

  date_status <- .sw_bd_status(date_raw, date, unreadable)
  time_status <- .sw_bd_status(
    time_raw,
    datetime,
    unreadable,
    parsed_if = !is.na(time_ok)
  )
  serial <- pick(serial_hit, "value")

  # Row status and the reason, in words a reader can act on.
  reason <- vapply(
    seq_along(ids),
    function(i) {
      .sw_bd_reason(
        heads[[i]]$reason,
        date_status[i],
        time_status[i],
        date_keywords,
        time_keywords,
        is.na(serial[i])
      )
    },
    character(1)
  )
  status <- ifelse(
    unreadable | date_status != "parsed",
    "not_assessed",
    ifelse(time_status != "parsed", "date_only", "assessed")
  )

  # Declared labels from the manifest; files without a row get NA.
  man <- tibble::as_tibble(manifest)
  hit <- match(ids, as.character(man$sample_id))
  has_run <- "run" %in% names(man)
  batch <- as.character(man$batch)[hit]
  run <- if (has_run) {
    as.character(man$run)[hit]
  } else {
    rep(NA_character_, length(ids))
  }

  files <- tibble::tibble(
    file = basename(unname(paths)),
    path = unname(paths),
    sample_id = ids,
    in_manifest = !is.na(hit),
    batch = batch,
    run = run,
    date_keyword = pick(date_hit, "key"),
    date_raw = date_raw,
    date = date,
    date_status = date_status,
    time_keyword = pick(time_hit, "key"),
    time_raw = time_raw,
    time_status = time_status,
    datetime = datetime,
    serial = serial
  )
  sess <- .sw_bd_sessions(files$datetime, files$batch, run_gap_hours)
  files$gap_before_hours <- sess$gap_before_hours
  files$session <- sess$session
  files$status <- status
  files$reason <- reason

  groups <- .sw_bd_groups(files, has_run, run_gap_hours)
  out <- list(
    files = files,
    crosstab = .sw_bd_crosstab(files),
    groups = groups,
    unmatched = .sw_bd_unmatched(files, man, has_run),
    provenance = .sw_bd_provenance(
      list(
        n_paths = length(paths),
        date_keywords = date_keywords,
        time_keywords = time_keywords,
        serial_keyword = serial_keyword,
        date_format = date_format,
        time_zone = time_zone,
        run_gap_hours = run_gap_hours
      )
    )
  )
  # The attribute is the convention; the element stays for one release.
  attr(out, "provenance") <- out$provenance
  class(out) <- c("sw_fcs_batch_date_check", "list")
  out
}

# Argument checks. Returns the sample ids (one per path). No package is loaded.
.sw_bd_check_args <- function(
  paths,
  manifest,
  date_keywords,
  time_keywords,
  serial_keyword,
  date_format,
  time_zone,
  run_gap_hours
) {
  is_str <- function(x) {
    is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
  }
  is_strs <- function(x) {
    is.character(x) && length(x) > 0 && !anyNA(x) && all(nzchar(x))
  }
  if (!is.character(paths) || length(paths) == 0 || anyNA(paths)) {
    stop(
      "`paths` must be a non-empty character vector of FCS paths.",
      call. = FALSE
    )
  }
  if (!inherits(manifest, "sw_manifest")) {
    stop("`manifest` must come from sw_sample_manifest().", call. = FALSE)
  }
  if (!is_strs(date_keywords)) {
    stop("`date_keywords` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is_strs(time_keywords)) {
    stop("`time_keywords` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is_str(serial_keyword)) {
    stop("`serial_keyword` must be a single keyword.", call. = FALSE)
  }
  if (!is_str(date_format)) {
    stop("`date_format` must be a single strptime format.", call. = FALSE)
  }
  if (!is_str(time_zone) || !time_zone %in% OlsonNames()) {
    stop("`time_zone` must be one name from OlsonNames().", call. = FALSE)
  }
  if (
    !is.numeric(run_gap_hours) ||
      length(run_gap_hours) != 1 ||
      !is.finite(run_gap_hours) ||
      run_gap_hours <= 0
  ) {
    stop("`run_gap_hours` must be one positive number.", call. = FALSE)
  }
  nm <- names(paths)
  ids <- tools::file_path_sans_ext(basename(paths))
  if (!is.null(nm)) {
    ids <- ifelse(is.na(nm) | !nzchar(nm), ids, nm)
  }
  if (anyDuplicated(ids)) {
    stop(
      "Sample ids from `paths` must be unique; repeated: ",
      paste(unique(ids[duplicated(ids)]), collapse = ", "),
      ". Name the elements of `paths`.",
      call. = FALSE
    )
  }
  ids
}

# Header of the first data set of one file as a named character vector.
.sw_bd_read_header <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    return(list(h = NULL, reason = .sw_fcs_read_reason(path)))
  }
  tryCatch(
    {
      h <- unlist(.sw_read_fcs_header(path)[[1]])
      if (length(h) == 0) {
        stop("no header keywords found")
      }
      list(h = h, reason = NA_character_)
    },
    error = function(e) {
      list(h = NULL, reason = .sw_fcs_read_reason(path, e))
    },
    warning = function(w) {
      list(h = NULL, reason = .sw_fcs_read_reason(path, w))
    }
  )
}

# First keyword (in the order given) with a non-empty value.
.sw_bd_first_value <- function(x, keywords) {
  none <- list(key = NA_character_, value = NA_character_)
  if (is.null(x$h)) {
    return(none)
  }
  pos <- match(toupper(keywords), toupper(names(x$h)))
  val <- trimws(unname(x$h[pos]))
  ok <- !is.na(pos) & !is.na(val) & nzchar(val)
  if (!any(ok)) {
    return(none)
  }
  i <- which(ok)[1]
  list(key = keywords[i], value = val[i])
}

# Date text to Date; NA when the text is missing, does not fit the format, or
# gives a year below 1000 (a two-digit year read with %Y, such as 26 for 2026).
.sw_bd_parse_date <- function(raw, fmt) {
  lt <- strptime(raw, fmt, tz = "UTC")
  date <- as.Date(lt)
  date[!is.na(lt) & lt$year + 1900 < 1000] <- NA
  date
}

# Time text to "HH:MM:SS"; NA when missing or out of range.
.sw_bd_parse_time <- function(raw) {
  m <- regmatches(raw, regexec("^([0-9]{1,2}):([0-9]{2}):([0-9]{2})", raw))
  vapply(
    m,
    function(p) {
      if (length(p) != 4) {
        return(NA_character_)
      }
      v <- as.integer(p[2:4])
      if (v[1] > 23 || v[2] > 59 || v[3] > 59) {
        return(NA_character_)
      }
      sprintf("%02d:%02d:%02d", v[1], v[2], v[3])
    },
    character(1)
  )
}

# parsed / missing_keyword / unparseable / unreadable for one parsed field.
.sw_bd_status <- function(raw, parsed, unreadable, parsed_if = !is.na(parsed)) {
  ifelse(
    unreadable,
    "unreadable",
    ifelse(
      is.na(raw),
      "missing_keyword",
      ifelse(parsed_if, "parsed", "unparseable")
    )
  )
}

# Plain-text reason for one file; NA when nothing needs saying.
.sw_bd_reason <- function(
  read_reason,
  date_status,
  time_status,
  date_kw,
  time_kw,
  no_serial
) {
  if (!is.na(read_reason)) {
    return(read_reason)
  }
  msg <- c(
    if (date_status == "missing_keyword") {
      paste0("date: no value in ", paste(date_kw, collapse = ", "))
    },
    if (date_status == "unparseable") {
      "date: text does not fit date_format, or its year is below 1000"
    },
    if (date_status == "parsed" && time_status == "missing_keyword") {
      paste0("time: no value in ", paste(time_kw, collapse = ", "))
    },
    if (date_status == "parsed" && time_status == "unparseable") {
      "time: text is not hh:mm:ss"
    },
    if (no_serial) "serial: no value"
  )
  if (length(msg) == 0) NA_character_ else paste(msg, collapse = "; ")
}

# Hours between neighbours and session numbers for one set of times.
# Returns gaps (length n_ok - 1, in time order) and, per input row, the session
# number and the gap to the previous file (NA for the first or an unusable one).
.sw_bd_order <- function(dt, limit) {
  n <- length(dt)
  ok <- which(!is.na(dt))
  o <- ok[order(dt[ok])]
  gap <- rep(NA_real_, n)
  session <- rep(NA_integer_, n)
  if (length(o) > 0) {
    g <- c(NA_real_, diff(as.numeric(dt[o])) / 3600)
    gap[o] <- g
    session[o] <- as.integer(cumsum(c(TRUE, g[-1] > limit)))
  }
  list(gap = gap, session = session)
}

# Gap before each file and session number, inside each declared batch.
.sw_bd_sessions <- function(datetime, batch, limit) {
  gap <- rep(NA_real_, length(datetime))
  session <- rep(NA_integer_, length(datetime))
  for (b in unique(batch[!is.na(batch)])) {
    i <- which(batch %in% b)
    s <- .sw_bd_order(datetime[i], limit)
    gap[i] <- s$gap
    session[i] <- s$session
  }
  list(gap_before_hours = gap, session = session)
}

# One summary row for the files in `i` (row indices of `files`).
.sw_bd_group_row <- function(files, i, level, batch, run, limit) {
  f <- files[i, ]
  s <- .sw_bd_order(f$datetime, limit)
  have <- !is.na(f$datetime)
  dates <- sort(unique(f$date[!is.na(f$date)]))
  serials <- sort(unique(f$serial[!is.na(f$serial)]))
  gaps <- s$gap[!is.na(s$gap)]
  tibble::tibble(
    level = level,
    batch = batch,
    run = run,
    n_files = nrow(f),
    n_assessed = sum(have),
    n_not_assessed = sum(!have),
    n_dates = length(dates),
    dates = paste(format(dates), collapse = ";"),
    n_serials = length(serials),
    serials = paste(serials, collapse = ";"),
    first_datetime = if (any(have)) {
      min(f$datetime[have])
    } else {
      as.POSIXct(NA, tz = "UTC")
    },
    last_datetime = if (any(have)) {
      max(f$datetime[have])
    } else {
      as.POSIXct(NA, tz = "UTC")
    },
    span_hours = if (any(have)) {
      as.numeric(difftime(
        max(f$datetime[have]),
        min(f$datetime[have]),
        units = "hours"
      ))
    } else {
      NA_real_
    },
    max_gap_hours = if (length(gaps) > 0) max(gaps) else NA_real_,
    n_gaps_over = sum(gaps > limit),
    n_sessions = if (any(have)) max(s$session, na.rm = TRUE) else 0L,
    run_gap_hours = limit
  )
}

# Group rows: one per declared batch, then one per declared batch and run.
.sw_bd_groups <- function(files, has_run, limit) {
  rows <- list()
  for (b in unique(files$batch[!is.na(files$batch)])) {
    rows[[length(rows) + 1]] <- .sw_bd_group_row(
      files,
      which(files$batch %in% b),
      "batch",
      b,
      NA_character_,
      limit
    )
  }
  if (has_run) {
    keys <- unique(files[
      !is.na(files$batch) & !is.na(files$run),
      c("batch", "run")
    ])
    for (k in seq_len(nrow(keys))) {
      i <- which(files$batch %in% keys$batch[k] & files$run %in% keys$run[k])
      rows[[length(rows) + 1]] <- .sw_bd_group_row(
        files,
        i,
        "run",
        keys$batch[k],
        keys$run[k],
        limit
      )
    }
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    out <- .sw_bd_group_row(
      files,
      integer(0),
      "batch",
      NA_character_,
      NA_character_,
      limit
    )[0, ]
  }
  out
}

# File counts for every batch and date, empty cells included. NA is its own level.
.sw_bd_crosstab <- function(files) {
  batches <- unique(files$batch)
  dates <- sort(unique(files$date[!is.na(files$date)]))
  grid <- expand.grid(
    batch = batches,
    date = c(dates, as.Date(NA)),
    stringsAsFactors = FALSE
  )
  same <- function(x, y) {
    (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
  }
  grid$n_files <- vapply(
    seq_len(nrow(grid)),
    function(k) {
      sum(same(files$batch, grid$batch[k]) & same(files$date, grid$date[k]))
    },
    numeric(1)
  )
  out <- tibble::as_tibble(grid)
  # Keep the empty (batch, NA date) cells only where a file needs them.
  out <- out[order(is.na(out$date), out$date, out$batch), ]
  out[c("batch", "date", "n_files")]
}

# Files without a manifest row and manifest rows without a file.
.sw_bd_unmatched <- function(files, man, has_run) {
  lone_file <- files[!files$in_manifest, ]
  lone_row <- man[!as.character(man$sample_id) %in% files$sample_id, ]
  tibble::tibble(
    side = c(
      rep("file_without_manifest_row", nrow(lone_file)),
      rep("manifest_row_without_file", nrow(lone_row))
    ),
    sample_id = c(lone_file$sample_id, as.character(lone_row$sample_id)),
    path = c(lone_file$path, rep(NA_character_, nrow(lone_row))),
    batch = c(
      rep(NA_character_, nrow(lone_file)),
      as.character(lone_row$batch)
    ),
    run = c(
      rep(NA_character_, nrow(lone_file)),
      if (has_run) {
        as.character(lone_row$run)
      } else {
        rep(NA_character_, nrow(lone_row))
      }
    )
  )
}

# Provenance block: the shared import block (R/import_util.R). Nothing here is
# random, so the seed slot is NA.
.sw_bd_provenance <- function(args) {
  .sw_import_provenance(args = args, seed = NA_integer_)
}
