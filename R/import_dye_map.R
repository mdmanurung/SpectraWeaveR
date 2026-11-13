# R/import_dye_map.R
# Marker x batch map of dyes, read from FCS headers only.
# Clean-room implementation from a written idea (marker and dye are often
# packed into one channel label, and the dye behind a marker can change
# between batches); no upstream code is used.

#' Map markers to dyes from FCS channel labels, per batch
#'
#' @description
#' Read only the header of each FCS file, split each channel label into a
#' marker and a dye, and tabulate the dyes each marker carries in each batch.
#' A label that does not split is listed with its file and channel. It is
#' never guessed. The result has no pass or fail and no per-file score.
#'
#' @details
#' **What is split.** Only `$PnS` (the marker label) is split. A channel with
#' an empty `$PnS` keeps `$PnN` as its `name` (`source` is `$PnN`), is never
#' split and is listed as `not_split`, because `$PnN` is usually a detector
#' name such as `FSC-A` or `BV421-A`, not a marker label. Files that pack
#' marker and dye into `$PnN` only, or that hold the dye in `$PnN` and the
#' marker in `$PnS`, are therefore not mapped; which layout to support is
#' open and needs the user's decision. White space around the label is
#' trimmed. `separator` is a regular expression (perl-style). The
#' label must split into exactly two non-empty parts. The **first part is the
#' marker and the second is the dye**, as in `"CD3 BV421"` with
#' `separator = " "`. A label such as `"BV421-A"` (dye first, then a
#' detector suffix) has no marker, so `separator = "-"` would call `BV421` the
#' marker and `A` the dye. Choose a `separator` that fits your labels and
#' check the `unsplit` table.
#'
#' **No defaults.** `paths`, `batch` and `separator` change every row of the
#' result, so each is required.
#'
#' **Nothing dropped.** Every channel of the first dataset of every readable
#' file is a row of `map`. A file that cannot be read gives one `not_assessed`
#' row with the reason, and so does each further dataset of a multi-dataset
#' file (`channel` is `NA`). Rows that are not `split` are repeated in `unsplit`.
#'
#' **Dye counts.** `n_dyes_in_batch` is the number of different dyes a marker
#' has inside one batch. `n_dyes_marker` is the number of different dyes the
#' same marker has across all batches. A value above 1 means the marker
#' carries different dyes; the function does not decide what that means.
#'
#' @param paths Character vector of FCS file paths. File names (without the
#'   folder) must be unique, because they label the rows. Required.
#' @param batch Named character vector of batch labels. Names are the file
#'   names with extension (`basename(paths)`, the `file` column of the other
#'   import tables; not the stem used as `sample_id` elsewhere); every file needs one label, and a name that
#'   matches no file is an error. Required.
#' @param separator Single regular expression that splits a channel label
#'   into marker (first part) and dye (second part). Required.
#'
#' @return A list of class `sw_fcs_marker_dye_map` with tibbles
#'   \describe{
#'   \item{`map`}{one row per file and channel: `file`, `batch`,
#'     `channel_index`, `channel` (`$PnN`), `name` (the split label),
#'     `source` (`$PnS` or `$PnN`), `marker`, `dye`, `status` (`split`,
#'     `not_split` or `not_assessed`) and `reason`.}
#'   \item{`marker_batch`}{one row per marker and batch: `marker`, `batch`,
#'     `dyes` (sorted, joined by `; `), `n_files`, `n_dyes_in_batch` and
#'     `n_dyes_marker`.}
#'   \item{`unsplit`}{the rows of `map` whose status is not `split`.}
#'   }
#'   Attribute `provenance` starts with the shared import keys
#'   (`package_version`, `flowcore_version`, `r_version`, `platform`,
#'   `datetime`, `seed` (`NULL`, nothing is random), `session`, `args` with
#'   `paths`, `batch` and `separator`). It then adds the unmix-block keys
#'   `packages` (versions of SpectraWeaveR, AutoSpectral, AutoSpectralRcpp and
#'   flowCore), `files` (MD5 of each input) and the backend and baseline
#'   fields, plus `separator` and `n_files`.
#'
#' @seealso [sw_fcs_header_audit()]
#'
#' @examples
#' \dontrun{
#' paths <- list.files("fcs", pattern = "[.]fcs$", full.names = TRUE)
#' batch <- setNames(c("B1", "B1", "B2"), basename(paths))
#' res <- sw_fcs_marker_dye_map(paths, batch, separator = " ")
#' res$marker_batch[res$marker_batch$n_dyes_marker > 1, ]
#' res$unsplit
#' }
#'
#' @export
sw_fcs_marker_dye_map <- function(paths, batch, separator) {
  .sw_stop_missing(c(
    paths = missing(paths),
    batch = missing(batch),
    separator = missing(separator)
  ))
  if (!is.character(paths) || length(paths) == 0 || anyNA(paths)) {
    stop("`paths` must be a non-empty character vector.", call. = FALSE)
  }
  files <- basename(paths)
  if (anyDuplicated(files)) {
    stop("File names in `paths` must be unique.", call. = FALSE)
  }
  if (!is.character(batch) || is.null(names(batch)) || anyNA(names(batch))) {
    stop(
      "`batch` must be a character vector named by file name ",
      "(`basename(paths)`), e.g. `setNames(batch_labels, basename(paths))`.",
      call. = FALSE
    )
  }
  if (anyNA(batch)) {
    stop("`batch` has NA labels.", call. = FALSE)
  }
  no_label <- setdiff(files, names(batch))
  if (length(no_label) > 0) {
    stop(
      "`batch` has no label for: ",
      paste(no_label, collapse = ", "),
      call. = FALSE
    )
  }
  extra <- setdiff(names(batch), files)
  if (length(extra) > 0) {
    stop(
      "`batch` names match no file in `paths`: ",
      paste(extra, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(names(batch))) {
    stop("`batch` names must be unique.", call. = FALSE)
  }
  if (
    !is.character(separator) ||
      length(separator) != 1 ||
      is.na(separator) ||
      !nzchar(separator)
  ) {
    stop("`separator` must be one non-empty regular expression.", call. = FALSE)
  }
  bad_regex <- tryCatch(
    {
      grepl(separator, "", perl = TRUE)
      FALSE
    },
    error = \(e) TRUE,
    warning = \(w) TRUE
  )
  if (bad_regex) {
    stop("`separator` is not a valid regular expression.", call. = FALSE)
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  rows <- lapply(seq_along(paths), \(i) {
    .sw_dye_file_rows(paths[[i]], files[[i]], batch[[files[[i]]]], separator)
  })
  map <- dplyr::bind_rows(rows)
  unsplit <- map[map$status != "split", ]
  marker_batch <- .sw_dye_marker_batch(map[map$status == "split", ])

  out <- list(map = map, marker_batch = marker_batch, unsplit = unsplit)
  # Shared import block first; then the keys of the unmix block that are not in
  # it (package table, file MD5s, backend, baseline). Nothing is random.
  args <- list(paths = paths, batch = batch, separator = separator)
  unmix_prov <- .sw_unmix_provenance(inputs = paths, args = args, seed = NULL)
  prov <- .sw_import_provenance(
    args = args,
    seed = NULL,
    extra = c(
      unmix_prov[setdiff(
        names(unmix_prov),
        c("r_version", "platform", "datetime", "args", "seed")
      )],
      list(separator = separator, n_files = length(paths))
    )
  )
  attr(out, "provenance") <- prov
  class(out) <- c("sw_fcs_marker_dye_map", class(out))
  out
}

# One not_assessed row: a file, or a dataset inside a file, that was not mapped.
.sw_dye_na_row <- function(file, batch, reason) {
  tibble::tibble(
    file = file,
    batch = batch,
    channel_index = NA_integer_,
    channel = NA_character_,
    name = NA_character_,
    source = NA_character_,
    marker = NA_character_,
    dye = NA_character_,
    status = "not_assessed",
    reason = reason
  )
}

# Rows of one file: one per channel of the first dataset, one not_assessed row
# per further dataset, or one not_assessed row if the file is unreadable.
.sw_dye_file_rows <- function(path, file, batch, separator) {
  # flowCore warns about files it cannot open; the not_assessed row says why.
  hl <- if (file.exists(path)) {
    tryCatch(
      suppressWarnings(.sw_read_fcs_header(path)),
      error = \(e) e
    )
  } else {
    NULL
  }
  h <- if (is.list(hl) && length(hl) > 0) unlist(hl[[1]]) else NULL
  if (is.null(names(h)) || !any(grepl("^\\$P[0-9]+N$", names(h)))) {
    reason <- if (is.null(names(h))) {
      .sw_fcs_read_reason(path, if (inherits(hl, "condition")) hl)
    } else {
      "header has no channels"
    }
    return(.sw_dye_na_row(file, batch, reason))
  }

  # Channels come back in key order; sort them by their real $Pn number.
  ch <- .sw_hdr_items(h, keywords = character())
  pn_num <- as.integer(sub(
    "^\\$P([0-9]+)N$",
    "\\1",
    grep("^\\$P[0-9]+N$", names(h), value = TRUE)
  ))
  ch <- ch[order(pn_num), ]
  pn_num <- sort(pn_num)
  label <- ch$value
  # An empty $PnS has no marker label. $PnN is then a detector name, so it is
  # listed and never split.
  use_pnn <- is.na(label) | !nzchar(label)
  name <- if_else_chr(use_pnn, trimws(ch$key), label)
  parts <- strsplit(name, separator, perl = TRUE)
  ok <- !use_pnn &
    lengths(parts) == 2 &
    vapply(parts, \(p) all(!is.na(p) & nzchar(trimws(p))), logical(1))
  pick <- \(k) {
    vapply(
      parts,
      \(p) if (length(p) >= k) trimws(p[[k]]) else NA_character_,
      ""
    )
  }
  rows <- tibble::tibble(
    file = file,
    batch = batch,
    channel_index = pn_num,
    channel = ch$key,
    name = name,
    source = ifelse(use_pnn, "$PnN", "$PnS"),
    marker = ifelse(ok, pick(1), NA_character_),
    dye = ifelse(ok, pick(2), NA_character_),
    status = ifelse(ok, "split", "not_split"),
    reason = ifelse(
      use_pnn,
      "$PnS is empty; $PnN is a detector name and is not split",
      ifelse(
        ok,
        NA_character_,
        sprintf(
          "label has %d part(s) for separator '%s'",
          lengths(parts),
          separator
        )
      )
    )
  )
  if (length(hl) == 1) {
    return(rows)
  }
  skipped <- vapply(
    seq_along(hl)[-1],
    \(k) {
      sprintf(
        "dataset %d of %d not read; only the first is mapped",
        k,
        length(hl)
      )
    },
    ""
  )
  dplyr::bind_rows(rows, .sw_dye_na_row(file, batch, skipped))
}

# Plain vectorised if-else for character vectors.
if_else_chr <- function(test, yes, no) {
  out <- no
  out[test] <- yes[test]
  out
}

# One row per marker and batch, with dye counts.
.sw_dye_marker_batch <- function(split_rows) {
  if (nrow(split_rows) == 0) {
    return(tibble::tibble(
      marker = character(),
      batch = character(),
      dyes = character(),
      n_files = integer(),
      n_dyes_in_batch = integer(),
      n_dyes_marker = integer()
    ))
  }
  per <- split_rows |>
    dplyr::summarise(
      dyes = paste(sort(unique(dye)), collapse = "; "),
      n_files = dplyr::n_distinct(file),
      n_dyes_in_batch = dplyr::n_distinct(dye),
      .by = c(marker, batch)
    )
  all_dyes <- split_rows |>
    dplyr::summarise(n_dyes_marker = dplyr::n_distinct(dye), .by = marker)
  per |>
    dplyr::left_join(all_dyes, by = "marker", relationship = "many-to-one") |>
    dplyr::arrange(marker, batch)
}
