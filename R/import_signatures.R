# R/import_signatures.R
# Import step I3 (with the folded-in I8 panel match): group FCS files by identical
# header keyword values and channel sets, and list markers missing from and extra
# to reference panels. Clean rewrite from the idea in plans/2026-11-13_import-step.md;
# no upstream code is used.

#' Group FCS files by header signature and match them to marker panels
#'
#' @description
#' Read only the header of each FCS file. Files with identical values for the
#' keywords you name and an identical channel set share one signature. The
#' result lists every file, counts the files in each signature, and, when you
#' give a `panel`, lists the markers each file lacks and the markers it has
#' beyond each panel entry. No reference file is needed.
#'
#' Nothing is dropped and no file is marked pass or fail. A file that cannot be
#' read gets a `not_assessed` row with a reason, and the other files are still
#' processed. `sw_fcs_header_audit()` differs here: it stops the whole call when
#' one file cannot be read.
#'
#' @details
#' **No defaults.** `keywords` and `panel` change which files group together and
#' which markers count as missing, so you choose them. `panel = NULL` is an
#' explicit skip: the panel table is then `NULL`.
#'
#' **Signature.** A signature is the set of values of `keywords` plus the set of
#' channel names (`$PnN`), compared as text after trimming white space. Channel
#' order does not matter. A keyword the file lacks counts as its own value, so
#' two files that both lack it match on it. Keywords of the form `$PnX` (for
#' example `"$PnV"`) expand per channel, as in [sw_fcs_header_audit()]. Marker
#' names (`$PnS`) are not part of the signature; they enter only through
#' `panel`. Signature ids (`sig01`, `sig02`, ...) run from the largest group to
#' the smallest; ties follow the order of `paths`.
#'
#' **Panel match.** For each file and each panel entry, the markers of the file
#' are the non-empty `$PnS` values. Matching is exact text, case sensitive, with
#' no alias map. A file marker that is not in the entry is `extra`; an entry
#' marker the file lacks is `missing`. Channels with no `$PnS` (scatter, time)
#' are counted in `n_no_marker` and never called extra.
#'
#' **Several data sets.** An FCS file can hold more than one data set
#' (`$NEXTDATA`). Only data set 1 enters the signature. `n_datasets` in `files`
#' is 2 when `$NEXTDATA` is non-zero (the exact count is not read), and `reason` says so.
#'
#' **Keywords found nowhere.** The `keywords` table counts, per keyword, the
#' files that have it and the files that lack it. A keyword that no file has
#' (a typo, wrong case, a missing `$`) adds nothing to the grouping, and its
#' `reason` says so.
#'
#' **Scope.** Only headers are read, never events. Files are never changed.
#'
#' @param paths Character vector of FCS file paths. File names (without the
#'   folder) must be unique, because they label the rows. Required.
#' @param keywords Character vector of header keywords whose values define the
#'   signature, for example `c("$CYT", "$PnV")`. No empty strings. Required.
#' @param panel A named list of reference marker sets, for example
#'   `list(T = c("CD3", "CD4"), B = c("CD3", "CD19"))`, or `NULL` to skip the
#'   panel match. Required.
#'
#' @return A list of class `sw_fcs_header_signatures` with
#'   * `files`: one row per path with `file`, `path`, `status` (`assessed` or
#'     `not_assessed`), `reason`, `n_datasets`, `signature`, `n_channels` and
#'     `definition`. `signature`, `n_channels`, `n_datasets` and `definition` are
#'     `NA` when not assessed. `reason` is empty for an `assessed` file with one data
#'     set.
#'   * `counts`: one row per signature with `signature`, `n_files`, `n_channels`,
#'     `files` (names joined by `"; "`) and `definition`. Not-assessed files are
#'     not in this table; see `n_not_assessed`.
#'   * `keywords`: one row per keyword with `keyword`, `n_files_present`,
#'     `n_files_absent`, `n_not_assessed` and `reason` (`NA` when the keyword is found). A `$PnX` keyword counts as
#'     present in a file when any channel has a value.
#'   * `panel_match`: `NULL` when `panel = NULL`, else one row per file and panel
#'     entry with `file`, `panel`, `status`, `reason`, `n_panel`, `n_matched`,
#'     `n_missing`, `n_extra`, `n_no_marker`, `missing` and `extra` (names
#'     joined by `"; "`).
#'   * `n_not_assessed`: number of files without a signature.
#'   * `provenance`: package and R versions, `seed` (`NULL`, nothing is random),
#'     session summary and the arguments.
#'
#' @seealso [sw_fcs_header_audit()]
#'
#' @examples
#' \dontrun{
#' res <- sw_fcs_header_signatures(
#'   list.files("fcs", full.names = TRUE),
#'   keywords = c("$CYT", "$CYTSN", "$PnV"),
#'   panel = list(T = c("CD3", "CD4", "CD8"))
#' )
#' res$counts
#' res$panel_match
#' }
#'
#' @export
sw_fcs_header_signatures <- function(paths, keywords, panel) {
  .sw_stop_missing(c(
    paths = missing(paths),
    keywords = missing(keywords),
    panel = missing(panel)
  ))
  if (!is.character(paths) || length(paths) == 0 || anyNA(paths)) {
    stop(
      "`paths` must be a non-empty character vector of file paths.",
      call. = FALSE
    )
  }
  if (anyDuplicated(basename(paths))) {
    stop("File names in `paths` must be unique.", call. = FALSE)
  }
  if (
    !is.character(keywords) ||
      length(keywords) == 0 ||
      anyNA(keywords) ||
      !all(nzchar(trimws(keywords)))
  ) {
    stop(
      "`keywords` must be a non-empty character vector without empty strings.",
      call. = FALSE
    )
  }
  if (!is.null(panel)) {
    nm <- names(panel)
    ok <- is.list(panel) &&
      length(panel) > 0 &&
      !is.null(nm) &&
      !anyNA(nm) &&
      all(nzchar(nm)) &&
      !anyDuplicated(nm) &&
      all(vapply(
        panel,
        function(x) is.character(x) && length(x) > 0 && !anyNA(x),
        logical(1)
      ))
    if (!ok) {
      stop(
        "`panel` must be NULL or a named list with unique names, each a ",
        "non-empty character vector of markers.",
        call. = FALSE
      )
    }
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  labels <- basename(paths)

  # One row per file: its signature parts, or the reason it could not be read.
  parsed <- lapply(seq_along(paths), function(i) {
    .sw_sig_parse(paths[i], keywords)
  })
  ok <- vapply(parsed, function(x) is.na(x$reason), logical(1))

  keys <- vapply(parsed, function(x) x$key, character(1))
  defs <- vapply(parsed, function(x) x$definition, character(1))
  nchan <- vapply(parsed, function(x) x$n_channels, integer(1))
  keys[!ok] <- NA_character_
  ndata <- vapply(parsed, function(x) x$n_datasets, integer(1))

  # Signature ids: largest group first, ties by first appearance in `paths`.
  tab <- table(factor(keys[ok], levels = unique(keys[ok])))
  levs <- names(tab)[order(-as.integer(tab), seq_along(tab))]
  width <- max(2L, nchar(length(levs)))
  ids <- sprintf(paste0("sig%0", width, "d"), seq_along(levs))
  sig <- ids[match(keys, levs)]

  files <- tibble::tibble(
    file = labels,
    path = paths,
    status = ifelse(ok, "assessed", "not_assessed"),
    reason = vapply(parsed, function(x) x$reason, character(1)),
    n_datasets = ndata,
    signature = sig,
    n_channels = ifelse(ok, nchan, NA_integer_),
    definition = ifelse(ok, defs, NA_character_)
  )
  files$reason[ok] <- ifelse(
    ndata[ok] > 1L,
    "file has more than one data set ($NEXTDATA is not 0); only data set 1 is used",
    ""
  )

  counts <- tibble::tibble(
    signature = ids,
    n_files = as.integer(tab[levs]),
    n_channels = nchan[ok][match(levs, keys[ok])],
    files = vapply(
      levs,
      function(k) paste(labels[ok & keys %in% k], collapse = "; "),
      character(1)
    ),
    definition = defs[ok][match(levs, keys[ok])]
  )

  kw_names <- unique(keywords)
  present <- vapply(
    kw_names,
    function(k) {
      sum(vapply(parsed[ok], function(x) isTRUE(x$kw_present[[k]]), logical(1)))
    },
    integer(1)
  )
  kw_tab <- tibble::tibble(
    keyword = kw_names,
    n_files_present = unname(present),
    n_files_absent = sum(ok) - unname(present),
    n_not_assessed = sum(!ok),
    reason = ifelse(
      sum(ok) > 0 & unname(present) == 0L,
      "no file has this keyword; it adds nothing to the grouping",
      NA_character_
    )
  )

  panel_match <- NULL
  if (!is.null(panel)) {
    panel_match <- dplyr::bind_rows(lapply(seq_along(paths), function(i) {
      .sw_sig_panel_rows(labels[i], parsed[[i]], panel)
    }))
  }

  prov <- .sw_import_provenance(
    args = list(paths = paths, keywords = keywords, panel = panel),
    seed = NULL
  )
  out <- list(
    files = files,
    counts = counts,
    keywords = kw_tab,
    panel_match = panel_match,
    n_not_assessed = sum(!ok),
    provenance = prov
  )
  # The attribute is the convention; the element stays for one release.
  attr(out, "provenance") <- prov
  class(out) <- "sw_fcs_header_signatures"
  out
}

# Read one header and build its signature key; never stops on a bad file.
.sw_sig_parse <- function(path, keywords) {
  fail <- function(reason) {
    list(
      reason = reason,
      key = NA_character_,
      definition = NA_character_,
      n_channels = NA_integer_,
      n_datasets = NA_integer_,
      kw_present = list(),
      markers = character(),
      n_no_marker = NA_integer_
    )
  }
  if (!file.exists(path) || dir.exists(path)) {
    return(fail(.sw_fcs_read_reason(path)))
  }
  # flowCore may leave a connection open on a bad file and warn at garbage collection.
  # read.FCSheader() returns data set 1 only, so extra data sets come from $NEXTDATA.
  h <- tryCatch(.sw_bind_read_header(path), error = function(e) e)
  if (inherits(h, "error")) {
    return(fail(.sw_fcs_read_reason(path, h)))
  }
  h <- unlist(h)
  next_data <- if ("$NEXTDATA" %in% names(h)) {
    suppressWarnings(as.numeric(trimws(h[["$NEXTDATA"]])))
  } else {
    0
  }
  n_datasets <- if (!is.na(next_data) && next_data != 0) 2L else 1L
  items <- .sw_hdr_items(h, keywords)
  chan <- sort(items$key[items$kind == "channel"])
  if (length(chan) == 0) {
    return(fail("header has no $PnN channels"))
  }
  marker <- items$value[items$kind == "channel"]
  has_marker <- !is.na(marker) & nzchar(marker)

  kw <- items[items$kind != "channel", ]
  kw <- kw[order(kw$key, kw$channel, na.last = FALSE), ]
  kw_txt <- paste0(
    kw$key,
    ifelse(is.na(kw$channel), "", paste0("[", kw$channel, "]")),
    "=",
    ifelse(is.na(kw$value), "<absent>", kw$value)
  )
  definition <- paste0(
    paste(kw_txt, collapse = " | "),
    " | channels=",
    paste(chan, collapse = ",")
  )
  # The key is length-prefixed so no name or value can imitate a separator.
  enc <- function(x) ifelse(is.na(x), "~", paste0(nchar(x, "bytes"), ":", x))
  key <- paste(
    c(
      enc(kw$key),
      enc(kw$channel),
      enc(kw$value),
      paste0("n", length(chan)),
      enc(chan)
    ),
    collapse = "|"
  )
  kws <- unique(kw$key)
  kw_present <- as.list(vapply(
    kws,
    function(k) any(!is.na(kw$value[kw$key == k])),
    logical(1)
  ))
  list(
    reason = NA_character_,
    key = key,
    definition = definition,
    n_channels = length(chan),
    n_datasets = as.integer(n_datasets),
    kw_present = kw_present,
    markers = unique(marker[has_marker]),
    n_no_marker = sum(!has_marker)
  )
}

# Panel rows for one file: one per panel entry.
.sw_sig_panel_rows <- function(file, parsed, panel) {
  ok <- is.na(parsed$reason)
  rows <- lapply(names(panel), function(nm) {
    ref <- unique(trimws(panel[[nm]]))
    if (!ok) {
      return(tibble::tibble(
        file = file,
        panel = nm,
        status = "not_assessed",
        reason = parsed$reason,
        n_panel = length(ref),
        n_matched = NA_integer_,
        n_missing = NA_integer_,
        n_extra = NA_integer_,
        n_no_marker = NA_integer_,
        missing = NA_character_,
        extra = NA_character_
      ))
    }
    miss <- setdiff(ref, parsed$markers)
    extra <- setdiff(parsed$markers, ref)
    tibble::tibble(
      file = file,
      panel = nm,
      status = "assessed",
      reason = "",
      n_panel = length(ref),
      n_matched = length(intersect(ref, parsed$markers)),
      n_missing = length(miss),
      n_extra = length(extra),
      n_no_marker = parsed$n_no_marker,
      missing = paste(miss, collapse = "; "),
      extra = paste(extra, collapse = "; ")
    )
  })
  dplyr::bind_rows(rows)
}
