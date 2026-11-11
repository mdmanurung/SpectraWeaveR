# R/fcs_header_audit.R
# FCS header audit: compare header keywords and channel sets with a reference file.
# Clean-room implementation from a written method description (CytoTools
# getFCSmetadata idea: files from different instrument settings differ in header
# metadata); no upstream code is used.

#' Audit FCS file headers against a reference file
#'
#' @description
#' Read only the header of each FCS file and compare it with a reference file.
#' The result is a long table: one row per file and compared item. Nothing is
#' dropped and no file is marked pass or fail. The `status` column says what
#' differs, and you decide what to do about it.
#'
#' Three kinds of items are compared:
#'
#' * `channel`: each channel name (`$PnN`) and its marker (`$PnS`). A channel in
#'   the reference but not in the file is a `missing_in_file` row; a channel in
#'   the file but not in the reference is `absent_in_reference`.
#' * `keyword`: each plain keyword you name in `keywords`, for example
#'   `"$CYT"` or `"$DATE"`.
#' * `parameter_keyword`: a keyword of the form `$PnX` (for example `"$PnV"`
#'   for detector voltage). It is read once per channel name, so a voltage
#'   change on one detector shows as one row.
#'
#' @details
#' **No defaults.** `keywords` and `reference` have no default because each one
#' changes which files look different. Header keys differ between vendors, so
#' the package does not guess them. A keyword you name that the reference file
#' does not have is an error, since nothing could be compared with it.
#'
#' **Matching.** Values are compared as text after trimming white space. Key
#' names must match the header exactly as `flowCore::read.FCSheader()` returns
#' them (upper-case, leading `$` for standard keywords).
#'
#' **Scope.** Only headers are read, never events. The audit says nothing about
#' event data quality; use [sw_anomaly_features()] for that.
#'
#' @param paths Character vector of FCS file paths. File names (without the
#'   folder) must be unique, because they label the rows.
#' @param keywords Character vector of header keywords to compare. Plain names
#'   such as `"$CYT"`, or `"$PnX"` patterns that expand per channel. Required.
#' @param reference One path from `paths` to compare every file with. Required.
#'
#' @return A tibble of class `sw_fcs_header_audit` with columns `file`, `kind`
#'   (`channel`, `keyword` or `parameter_keyword`), `key`, `channel` (`NA` for
#'   plain keywords and for `channel` rows), `value`, `reference_value` and
#'   `status` (`match`, `differs`, `missing_in_file`, `absent_in_reference`).
#'   The reference file appears too, with every row `match`. Attribute
#'   `reference` holds the reference file label.
#'
#' @seealso [sw_fcs_header_summary()]
#'
#' @examples
#' \dontrun{
#' res <- sw_fcs_header_audit(
#'   list.files("fcs", full.names = TRUE),
#'   keywords = c("$CYT", "$CYTSN", "$PnV"),
#'   reference = "fcs/ref_sample.fcs"
#' )
#' sw_fcs_header_summary(res)
#' }
#'
#' @export
sw_fcs_header_audit <- function(paths, keywords, reference) {
  if (missing(paths)) {
    stop("`paths` is required.", call. = FALSE)
  }
  if (missing(keywords)) {
    stop("`keywords` is required and has no default.", call. = FALSE)
  }
  if (missing(reference)) {
    stop("`reference` is required and has no default.", call. = FALSE)
  }
  if (!is.character(paths) || length(paths) < 2 || anyNA(paths)) {
    stop("`paths` must be at least two FCS file paths.", call. = FALSE)
  }
  if (!is.character(keywords) || length(keywords) == 0 || anyNA(keywords)) {
    stop("`keywords` must be a non-empty character vector.", call. = FALSE)
  }
  if (
    !is.character(reference) ||
      length(reference) != 1 ||
      !reference %in% paths
  ) {
    stop("`reference` must be exactly one of `paths`.", call. = FALSE)
  }
  labels <- basename(paths)
  if (anyDuplicated(labels)) {
    stop("File names in `paths` must be unique.", call. = FALSE)
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  # One header per file, as a named character vector.
  headers <- flowCore::read.FCSheader(paths)
  headers <- lapply(headers, function(h) unlist(h))
  names(headers) <- labels
  ref_label <- basename(reference)

  items <- lapply(headers, .sw_hdr_items, keywords = keywords)
  ref_items <- items[[ref_label]]

  # The reference must hold every keyword, or nothing can be compared with it.
  ref_has <- tapply(
    !is.na(ref_items$value[ref_items$kind != "channel"]),
    ref_items$key[ref_items$kind != "channel"],
    any
  )
  lacking <- keywords[!ref_has[keywords] %in% TRUE]
  if (length(lacking) > 0) {
    stop(
      "Reference file lacks keyword(s): ",
      paste(lacking, collapse = ", "),
      call. = FALSE
    )
  }

  out <- lapply(labels, function(f) {
    .sw_hdr_compare(f, items[[f]], ref_items)
  }) |>
    dplyr::bind_rows()
  out <- tibble::as_tibble(out)
  attr(out, "reference") <- ref_label
  class(out) <- c("sw_fcs_header_audit", class(out))
  out
}

#' Count header differences per file
#'
#' @param audit An `sw_fcs_header_audit` result.
#'
#' @return A tibble with one row per file and one column per status, counting
#'   rows. Files with any status other than `match` come first.
#'
#' @examples
#' \dontrun{
#' sw_fcs_header_summary(res)
#' }
#'
#' @export
sw_fcs_header_summary <- function(audit) {
  if (!inherits(audit, "sw_fcs_header_audit")) {
    stop("`audit` must come from sw_fcs_header_audit().", call. = FALSE)
  }
  statuses <- c("match", "differs", "missing_in_file", "absent_in_reference")
  tab <- table(
    file = factor(audit$file, levels = unique(audit$file)),
    status = factor(audit$status, levels = statuses)
  )
  counts <- tibble::as_tibble(unclass(tab))
  counts <- tibble::add_column(counts, file = rownames(tab), .before = 1)
  n_diff <- counts$differs + counts$missing_in_file + counts$absent_in_reference
  counts[order(-n_diff, counts$file), ]
}

# Items of one header: channels, plain keywords and per-channel keywords.
.sw_hdr_items <- function(h, keywords) {
  pn <- grep("^\\$P[0-9]+N$", names(h), value = TRUE)
  idx <- sub("^\\$P([0-9]+)N$", "\\1", pn)
  chan <- unname(h[pn])
  marker <- unname(h[paste0("$P", idx, "S")])
  rows <- list(tibble::tibble(
    kind = "channel",
    key = chan,
    channel = NA_character_,
    value = trimws(marker)
  ))
  for (k in keywords) {
    if (grepl("^\\$P[nN][A-Za-z]$", k)) {
      suffix <- toupper(substr(k, 4, 4))
      v <- unname(h[paste0("$P", idx, suffix)])
      rows[[length(rows) + 1]] <- tibble::tibble(
        kind = "parameter_keyword",
        key = k,
        channel = chan,
        value = trimws(v)
      )
    } else {
      v <- if (k %in% names(h)) trimws(unname(h[[k]])) else NA_character_
      rows[[length(rows) + 1]] <- tibble::tibble(
        kind = "keyword",
        key = k,
        channel = NA_character_,
        value = v
      )
    }
  }
  dplyr::bind_rows(rows)
}

# Compare one file's items with the reference items; keep every row of both.
.sw_hdr_compare <- function(file, x, ref) {
  by <- c("kind", "key", "channel")
  # NA channel must match NA channel, so join on a text key.
  x$.id <- paste(x$kind, x$key, ifelse(is.na(x$channel), "", x$channel))
  ref$.id <- paste(
    ref$kind,
    ref$key,
    ifelse(is.na(ref$channel), "", ref$channel)
  )
  both <- dplyr::full_join(
    x[c(".id", by, "value")],
    ref[c(".id", "value")],
    by = ".id",
    suffix = c("", "_ref"),
    relationship = "one-to-one"
  )
  # Rows only in the reference lost their kind/key/channel; restore from ref.
  gap <- is.na(both$kind)
  both$kind[gap] <- ref$kind[match(both$.id[gap], ref$.id)]
  both$key[gap] <- ref$key[match(both$.id[gap], ref$.id)]
  both$channel[gap] <- ref$channel[match(both$.id[gap], ref$.id)]
  in_file <- both$.id %in% x$.id
  in_ref <- both$.id %in% ref$.id
  both$status <- dplyr::case_when(
    !in_ref ~ "absent_in_reference",
    !in_file ~ "missing_in_file",
    is.na(both$value) & is.na(both$value_ref) ~ "match",
    is.na(both$value) ~ "missing_in_file",
    .sw_hdr_equal(both$value, both$value_ref) ~ "match",
    .default = "differs"
  )
  tibble::tibble(
    file = file,
    kind = both$kind,
    key = both$key,
    channel = both$channel,
    value = both$value,
    reference_value = both$value_ref,
    status = both$status
  )
}

# Element-wise text equality that treats NA as unequal.
.sw_hdr_equal <- function(a, b) {
  !is.na(a) & !is.na(b) & a == b
}
