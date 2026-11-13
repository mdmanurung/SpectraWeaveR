# R/import_bind.R
# Import step I1: bind a sample manifest to FCS files by an exact key.
# Clean-room rewrite from the plan (plans/2026-11-13_import-step.md, I1). The
# matching idea follows .sw_ref_manifest() (R/ref_vs_study.R) and the join in
# sw_correct_prepare() (R/batch_correct.R); neither is changed, and unlike
# .sw_ref_manifest() nothing is dropped without a row.

#' Bind a sample manifest to FCS files
#'
#' @description
#' Join the rows of a [sw_sample_manifest()] to FCS files by an exact key. Only
#' file names and FCS headers are read, never events, and no file is changed.
#' Nothing is dropped without a row: every manifest row is in `bound` or
#' `manifest_no_file`, and every file is in `bound` or `file_no_manifest`.
#'
#' @details
#' **No defaults.** The match key, the case rule, recursion and the three
#' policies each change which files bind, so each is a required argument.
#'
#' **Key.** `key = "file_stem"` matches `sample_id` against the file name
#' without its extension (`a.fcs` has key `a`). `key = "file_name"` matches it
#' against the whole file name (`a.fcs`). The match is exact text after the case
#' rule: no pattern, no partial match, no fuzzy match.
#'
#' **Duplicates.** Two files with the same key (in different folders, or
#' differing only in case when `case = "ignore"`) cannot both be the file for
#' one sample. Two manifest `sample_id` values that fold to the same key (when
#' `case = "ignore"`) always stop, because the manifest itself is then
#' ambiguous. The same path given twice counts once.
#'
#' **Links.** The key comes from the name a file has in its folder, not from
#' the target of a symbolic link. Each link name is its own file, so links to one
#' real file each bind to their own row. Because those rows read the same
#' events, they are listed in `provenance$same_target` (`path`, `target`). A
#' link whose target is missing goes to `not_fcs` with reason `broken_link`.
#'
#' **`$NEXTDATA`.** A non-zero `$NEXTDATA` means the file holds more than one
#' data set. Only the first is read by the rest of the package. The value is
#' read from the header of every listed FCS file. A header that cannot be read,
#' or a `$NEXTDATA` that is missing or not a number, gives a `not_assessed`
#' row in `multi_dataset` with a reason. That file keeps its match, and the
#' `nextdata` policy does not act on it because its value is unknown.
#'
#' @param manifest A [sw_sample_manifest()] object. Required. It must not have
#'   a column named `path`, `file`, `nextdata`, `header_status` or `reason`,
#'   because the result adds those; the call stops and names them.
#' @param paths Character vector of FCS file paths and/or folders. A folder is
#'   listed (see `recursive`); every entry ending in `.fcs` (any case) is an FCS
#'   file and every other entry is returned in `not_fcs`. A path that does not
#'   exist stops the call. A file given directly is used whatever its extension.
#'   Required.
#' @param key How a file maps to `manifest$sample_id`: `"file_stem"` (name
#'   without extension) or `"file_name"` (name with extension). Required.
#' @param case `"exact"` (upper and lower case differ) or `"ignore"` (compare
#'   after lower-casing both sides; `bound` keeps the manifest spelling).
#'   Required.
#' @param recursive `TRUE` lists the sub-folders of each folder in `paths`,
#'   `FALSE` lists only its top level; a sub-folder is then not entered and
#'   appears in `not_fcs` with reason `directory_not_listed`. Hidden entries
#'   (names starting with a dot) are listed like any other. Required.
#' @param on_duplicate_file `"stop"` stops and names each duplicated key with
#'   its paths. `"unbind"` binds none of the copies: they appear in
#'   `file_no_manifest` and the manifest row appears in `manifest_no_file`,
#'   both with reason `duplicate_file_key`. Required.
#' @param on_extra_file `"list"` returns files with no manifest row in
#'   `file_no_manifest`. `"stop"` stops and names them. Required.
#' @param nextdata What to do with a matched file whose `$NEXTDATA` is not
#'   zero. `"bind"` binds it and lists it in `multi_dataset`. `"unbind"` does
#'   not bind it: it appears in `multi_dataset` (`bound = FALSE`), in
#'   `file_no_manifest` and in `manifest_no_file`, with reason
#'   `multi_dataset_unbound`. `"stop"` stops and names the files. Required.
#'
#' @return A list with six elements.
#'   * `bound`: a `sw_manifest` with one row per bound manifest row, in
#'     manifest order: every manifest column, then `path`, `file`, `nextdata` (the header text, `NA`
#'     when unread), `header_status` (`assessed` or `not_assessed`) and `reason`.
#'   * `manifest_no_file`: manifest rows with no bound file, plus `reason`
#'     (`manifest_row_without_file`, `duplicate_file_key` or `multi_dataset_unbound`).
#'   * `file_no_manifest`: files not bound, with `path`, `file`, `sample_id`
#'     (`NA` unless the file matched a row but was unbound) and `reason`
#'     (`file_without_manifest_row`, `duplicate_file_key` or `multi_dataset_unbound`).
#'   * `multi_dataset`: files with a non-zero `$NEXTDATA` or an unreadable one:
#'     `sample_id` (`NA` when unmatched), `path`, `file`, `nextdata`, `status`
#'     (`multi_dataset` or `not_assessed`), `bound` and `reason`.
#'   * `not_fcs`: entries of listed folders that are not used as FCS files, with
#'     `path` and `reason`: `not_fcs_extension`, `broken_link`,
#'     `directory_not_listed`, `directory_named_fcs` or `directory_unreadable`.
#'   * `provenance`: a list with the arguments, `n_manifest_rows`, `n_files`,
#'     `n_not_fcs`, `same_target` (files whose links share one real file), R version, platform, time, package versions
#'     (`SpectraWeaveR`, `flowCore`), the attached packages and `seed` (`NULL`,
#'     the function is deterministic).
#'
#' @seealso [sw_sample_manifest()], [sw_fcs_header_audit()]
#'
#' @examples
#' \dontrun{
#' man <- sw_sample_manifest(
#'   data.frame(id = c("s1", "s2"), run = c("b1", "b2")),
#'   sample_id = "id",
#'   batch = "run"
#' )
#' res <- sw_manifest_bind_files(
#'   man,
#'   paths = "fcs",
#'   key = "file_stem",
#'   case = "exact",
#'   recursive = FALSE,
#'   on_duplicate_file = "stop",
#'   on_extra_file = "list",
#'   nextdata = "bind"
#' )
#' res$manifest_no_file
#' }
#'
#' @export
sw_manifest_bind_files <- function(
  manifest,
  paths,
  key,
  case,
  recursive,
  on_duplicate_file,
  on_extra_file,
  nextdata
) {
  # unlist(list()) because c(recursive = ) would set c()'s own argument.
  # Required arguments first, then values, then dependencies.
  .sw_stop_missing(unlist(list(
    manifest = missing(manifest),
    paths = missing(paths),
    key = missing(key),
    case = missing(case),
    recursive = missing(recursive),
    on_duplicate_file = missing(on_duplicate_file),
    on_extra_file = missing(on_extra_file),
    nextdata = missing(nextdata)
  )))
  if (!inherits(manifest, "sw_manifest")) {
    stop("`manifest` must come from sw_sample_manifest().", call. = FALSE)
  }
  if (!is.character(paths) || length(paths) == 0 || anyNA(paths)) {
    stop("`paths` must be a non-empty character vector.", call. = FALSE)
  }
  .sw_bind_choice(key, "key", c("file_stem", "file_name"))
  .sw_bind_choice(case, "case", c("exact", "ignore"))
  .sw_bind_choice(on_duplicate_file, "on_duplicate_file", c("stop", "unbind"))
  .sw_bind_choice(on_extra_file, "on_extra_file", c("list", "stop"))
  .sw_bind_choice(nextdata, "nextdata", c("bind", "unbind", "stop"))
  if (!is.logical(recursive) || length(recursive) != 1 || is.na(recursive)) {
    stop("`recursive` must be TRUE or FALSE.", call. = FALSE)
  }
  clash <- intersect(
    c("path", "file", "nextdata", "header_status", "reason"),
    names(manifest)
  )
  if (length(clash) > 0) {
    stop(
      "`manifest` has column(s) that the result adds: ",
      paste(clash, collapse = ", "),
      ". Rename them first.",
      call. = FALSE
    )
  }
  absent <- paths[!file.exists(paths)]
  if (length(absent) > 0) {
    stop("Path(s) not found: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  # Manifest keys: the manifest must be unambiguous under the case rule.
  man_id <- as.character(manifest$sample_id)
  man_cmp <- .sw_bind_fold(man_id, case)
  if (anyDuplicated(man_cmp)) {
    stop(
      "`manifest$sample_id` values collide under case = \"",
      case,
      "\": ",
      paste(
        unique(man_id[man_cmp %in% man_cmp[duplicated(man_cmp)]]),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  # List files: folders are listed, files are taken as given.
  listed <- .sw_bind_list(paths, recursive)
  files <- listed$fcs
  file_cmp <- .sw_bind_fold(.sw_bind_key(files, key), case)

  # Duplicate file keys.
  dup <- file_cmp %in% file_cmp[duplicated(file_cmp)]
  if (any(dup) && on_duplicate_file == "stop") {
    stop(
      "Duplicate file key(s):\n",
      paste(
        vapply(
          unique(file_cmp[dup]),
          function(k) {
            paste0("  ", k, ": ", paste(files[file_cmp == k], collapse = ", "))
          },
          character(1)
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  # $NEXTDATA from the header of every listed file.
  hdr <- .sw_bind_nextdata(files)
  multi <- hdr$status == "multi_dataset"

  # Which file each manifest row would take.
  hit <- match(man_cmp, file_cmp)
  hit[!is.na(hit) & dup[hit]] <- NA_integer_
  man_reason <- rep(NA_character_, length(man_id))
  man_reason[man_cmp %in% file_cmp[dup]] <- "duplicate_file_key"
  man_reason[is.na(man_reason) & is.na(hit)] <- "manifest_row_without_file"

  unbound_multi <- !is.na(hit) & multi[hit]
  if (nextdata == "stop" && any(unbound_multi)) {
    stop(
      "File(s) with non-zero $NEXTDATA (more than one data set): ",
      paste(files[hit[unbound_multi]], collapse = ", "),
      call. = FALSE
    )
  }
  if (nextdata == "unbind") {
    man_reason[unbound_multi] <- "multi_dataset_unbound"
    hit[unbound_multi] <- NA_integer_
  }

  # Files not bound: no manifest row, duplicate key, or unbound multi-dataset.
  matched_file <- !is.na(match(file_cmp, man_cmp))
  file_reason <- rep(NA_character_, length(files))
  file_reason[!matched_file] <- "file_without_manifest_row"
  file_reason[dup] <- "duplicate_file_key"
  file_reason[
    is.na(file_reason) & !seq_along(files) %in% hit
  ] <- "multi_dataset_unbound"
  if (
    on_extra_file == "stop" && any(file_reason %in% "file_without_manifest_row")
  ) {
    stop(
      "File(s) with no manifest row: ",
      paste(
        files[file_reason %in% "file_without_manifest_row"],
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  file_sample <- man_id[match(file_cmp, man_cmp)]
  file_sample[file_reason %in% c("file_without_manifest_row")] <- NA_character_
  file_sample[file_reason %in% "duplicate_file_key"] <- NA_character_

  # Output tables.
  man_tbl <- tibble::as_tibble(unclass(manifest))
  bound_rows <- !is.na(hit)
  bound <- man_tbl[bound_rows, ]
  bound$path <- files[hit[bound_rows]]
  bound$file <- basename(bound$path)
  bound$nextdata <- hdr$nextdata[hit[bound_rows]]
  bound$header_status <- ifelse(
    hdr$status[hit[bound_rows]] == "not_assessed",
    "not_assessed",
    "assessed"
  )
  bound$reason <- hdr$reason[hit[bound_rows]]
  # Still a manifest, so later steps (e.g. sw_fcs_batch_date_check) accept it.
  class(bound) <- c("sw_manifest", class(bound))

  manifest_no_file <- man_tbl[!bound_rows, ]
  manifest_no_file$reason <- man_reason[!bound_rows]

  not_bound <- !is.na(file_reason)
  file_no_manifest <- tibble::tibble(
    path = files[not_bound],
    file = basename(files[not_bound]),
    sample_id = file_sample[not_bound],
    reason = file_reason[not_bound]
  )

  keep <- hdr$status != "ok"
  multi_dataset <- tibble::tibble(
    sample_id = file_sample[keep],
    path = files[keep],
    file = basename(files[keep]),
    nextdata = hdr$nextdata[keep],
    status = hdr$status[keep],
    bound = seq_along(files)[keep] %in% hit,
    reason = hdr$reason[keep]
  )
  multi_dataset$sample_id[is.na(multi_dataset$sample_id)] <-
    man_id[match(
      .sw_bind_fold(.sw_bind_key(multi_dataset$path, key), case),
      man_cmp
    )][
      is.na(multi_dataset$sample_id)
    ]

  not_fcs <- tibble::tibble(
    path = listed$other,
    reason = listed$other_reason
  )

  prov <- .sw_bind_provenance(
    args = list(
      key = key,
      case = case,
      recursive = recursive,
      on_duplicate_file = on_duplicate_file,
      on_extra_file = on_extra_file,
      nextdata = nextdata,
      paths = paths
    ),
    n_manifest_rows = nrow(man_tbl),
    n_files = length(files),
    n_not_fcs = length(listed$other),
    same_target = listed$same_target
  )

  # The attribute is the convention; the element stays for one release.
  out <- list(
    bound = bound,
    manifest_no_file = manifest_no_file,
    file_no_manifest = file_no_manifest,
    multi_dataset = multi_dataset,
    not_fcs = not_fcs,
    provenance = prov
  )
  attr(out, "provenance") <- prov
  class(out) <- c("sw_manifest_bind_files", "list")
  out
}

# One allowed value, checked by name.
.sw_bind_choice <- function(x, name, allowed) {
  if (!is.character(x) || length(x) != 1 || !x %in% allowed) {
    stop(
      "`",
      name,
      "` must be one of: ",
      paste0("\"", allowed, "\"", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(x)
}

.sw_bind_fold <- function(x, case) {
  if (case == "ignore") tolower(x) else x
}

# Key text of each path under the key rule.
.sw_bind_key <- function(paths, key) {
  nm <- basename(paths)
  if (key == "file_name") nm else tools::file_path_sans_ext(nm)
}

# Folders are listed (hidden entries too); files are used as given. Every
# entry that is not used as an FCS file is returned with a reason. The key
# comes from the name in the folder; each link name is its own file. Only a
# repeated identical path counts once. Names sharing one real file are returned
# in `same_target`.
.sw_bind_list <- function(paths, recursive) {
  fcs <- character()
  other <- character()
  why <- character()
  for (p in paths) {
    if (!dir.exists(p)) {
      fcs <- c(fcs, p)
      next
    }
    all <- list.files(
      p,
      full.names = TRUE,
      recursive = recursive,
      all.files = TRUE,
      include.dirs = TRUE,
      no.. = TRUE
    )
    is_dir <- dir.exists(all)
    reason <- rep(NA_character_, length(all))
    is_fcs <- grepl("\\.fcs$", all, ignore.case = TRUE)
    reason[is_dir & !recursive] <- "directory_not_listed"
    reason[
      is_dir & recursive & file.access(all, 4) != 0
    ] <- "directory_unreadable"
    reason[is_dir & recursive & is.na(reason) & is_fcs] <- "directory_named_fcs"
    reason[!is_dir & !file.exists(all)] <- "broken_link"
    reason[!is_dir & is.na(reason) & !is_fcs] <- "not_fcs_extension"
    # Readable sub-folders under recursive = TRUE are entered, not skipped.
    skip <- is_dir & recursive & is.na(reason)
    fcs <- c(fcs, all[!skip & is.na(reason)])
    other <- c(other, all[!skip & !is.na(reason)])
    why <- c(why, reason[!skip & !is.na(reason)])
  }
  fcs_abs <- .sw_bind_abs(fcs)
  ord <- order(fcs_abs)
  keep <- !duplicated(fcs_abs[ord])
  fcs_use <- fcs[ord][keep]
  fcs_real <- normalizePath(fcs_use, mustWork = TRUE)
  shared <- fcs_real %in% fcs_real[duplicated(fcs_real)]
  other_abs <- .sw_bind_abs(other)
  ord_o <- order(other_abs)
  keep_o <- !duplicated(other_abs[ord_o])
  list(
    fcs = fcs_abs[ord][keep],
    same_target = tibble::tibble(
      path = fcs_abs[ord][keep][shared],
      target = fcs_real[shared]
    ),
    other = other_abs[ord_o][keep_o],
    other_reason = why[ord_o][keep_o]
  )
}

# Absolute path that keeps the file's own name (a link is not resolved).
.sw_bind_abs <- function(x) {
  file.path(
    normalizePath(dirname(x), mustWork = FALSE),
    basename(x)
  )
}

# First header of a file. flowCore warns about extra data sets; that is the
# fact this function reports, so the warning is muffled and nothing else is.
.sw_bind_read_header <- function(f) {
  withCallingHandlers(
    .sw_read_fcs_header(f)[[1]],
    warning = function(w) {
      if (grepl("additional data segment", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# $NEXTDATA of each file's header; an unreadable header is not_assessed.
.sw_bind_nextdata <- function(files) {
  one <- function(f) {
    h <- tryCatch(
      .sw_bind_read_header(f),
      error = function(e) e
    )
    # The header reader can return the error instead of throwing it.
    if (!inherits(h, "condition")) {
      h <- unlist(h)
    }
    if (inherits(h, "condition") || is.null(names(h))) {
      return(list(
        nextdata = NA_character_,
        status = "not_assessed",
        reason = .sw_fcs_read_reason(f, if (inherits(h, "condition")) h)
      ))
    }
    val <- if ("$NEXTDATA" %in% names(h)) {
      trimws(h[["$NEXTDATA"]])
    } else {
      NA_character_
    }
    num <- suppressWarnings(as.numeric(val))
    if (is.na(num)) {
      return(list(
        nextdata = val,
        status = "not_assessed",
        reason = "$NEXTDATA missing or not a number"
      ))
    }
    list(
      nextdata = val,
      status = if (num != 0) "multi_dataset" else "ok",
      reason = NA_character_
    )
  }
  res <- lapply(files, one)
  tibble::tibble(
    nextdata = vapply(res, `[[`, character(1), "nextdata"),
    status = vapply(res, `[[`, character(1), "status"),
    reason = vapply(res, `[[`, character(1), "reason")
  )
}

# Provenance block: reuses .sw_pkg_info() from R/unmix_report.R.
.sw_bind_provenance <- function(args, ...) {
  .sw_import_provenance(
    args = args,
    seed = NULL,
    extra = c(list(...), list(packages = .sw_import_packages()))
  )
}
