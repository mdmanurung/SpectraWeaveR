# R/import_compensation.R
# Compensation state of FCS files, read from the header only. Clean-room
# implementation from the written idea (check whether a spillover matrix is
# stored, whether it is the identity, and what the "applied" flag says); no
# upstream code is used. The FCS spillover string layout follows the FCS 3.1
# standard: n, then n channel names, then n * n numbers in row order. The older
# "$COMP" layout (n, then n * n numbers, no names) is read too.

#' Report the compensation state of FCS files from their headers
#'
#' @description
#' Read only the header of each FCS file and report three separate facts about
#' compensation: is a spillover matrix stored, is that matrix the identity
#' matrix, and what does the "compensation applied" flag say. The facts stay in
#' separate columns. The function applies no precedence rule between them and
#' gives no combined verdict, no pass or fail and no per-file score. You read
#' the columns and decide.
#'
#' @details
#' **No defaults.** Every argument is required, because each one changes what
#' the table says. Header keys differ between vendors and software, so you name
#' the candidate keys. A key counts as present when it matches a header key
#' exactly, as `flowCore::read.FCSheader()` returns it (upper-case, leading `$`
#' for standard keywords). Values are trimmed of white space first.
#'
#' **Several keys.** If more than one of `matrix_keywords` is in a file, the
#' function does not pick one. `matrix_keys_matched` lists them all, and the
#' matrix columns for that file are `NA` with `matrix_status = "not_assessed"`
#' and a reason. The same holds for `apply_keywords`. Narrow the keyword vector
#' to see one reading.
#'
#' **Matrix string.** A stored matrix is parsed as n, then n channel names,
#' then n * n numbers (row order), separated by commas. The older layout without
#' names (n, then n * n numbers) is accepted. Anything else is malformed: the
#' matrix is not read, `matrix_status` is `"not_assessed"` and `matrix_reason`
#' says why (wrong count, text where a number belongs, empty or repeated
#' channel name, non-finite number, empty value).
#'
#' **Identity.** `matrix_max_dev` is the largest absolute difference between
#' the matrix and the identity matrix. `matrix_identity` is `TRUE` when it is at
#' most `zero_tol`.
#'
#' **Apply flag.** The raw text is kept in `apply_raw`. `apply_parsed` is `TRUE`
#' when the trimmed text equals one word of `apply_true` and `FALSE` when it
#' equals one word of `apply_false` (case ignored). You choose both word lists,
#' because vendors use flags like `1` and `0` with different meanings. Any other
#' text gives `NA` and `apply_status = "not_assessed"`; the raw text stays.
#'
#' **Nothing is dropped.** The result has one row per element of `paths`, in the
#' same order. A file that cannot be read gives a row with
#' `header_status = "not_assessed"` and a reason. A file named in
#' `expected_state` but not in `paths` also gets a row, with a reason.
#'
#' **Scope.** Only headers are read, never events, and no file is changed. The
#' function says nothing about whether the stored matrix is correct.
#'
#' @param paths Character vector of FCS file paths. File names (without the
#'   folder) must be unique, because they label the rows. Required.
#' @param matrix_keywords Character vector of header keys that may hold a
#'   spillover matrix, for example `c("$SPILLOVER", "$SPILL", "$COMP")`.
#'   Required.
#' @param apply_keywords Character vector of header keys that say compensation
#'   was applied. Required.
#' @param apply_true Character vector of flag words that mean "compensation
#'   applied", for example `c("true", "yes")`. Required.
#' @param apply_false Character vector of flag words that mean "not applied".
#'   No word may be in both lists. Required.
#' @param zero_tol Largest absolute distance of the stored matrix from the
#'   identity matrix that still counts as identity. One non-negative number.
#'   Required.
#' @param expected_state `NULL` to skip, or a data frame with one row per file.
#'   It needs a `file` column (file name as in `paths`, each once) and at least
#'   one of the logical columns `matrix_present`, `matrix_identity` and
#'   `apply_parsed`. An `NA` means no expectation for that item. Required.
#'
#' @return A tibble of class `sw_fcs_compensation_state` with one row per file
#'   and columns `file`, `path`, `header_status` (`read` or `not_assessed`),
#'   `header_reason`; for the matrix `matrix_status` (`assessed`, `not_present`
#'   or `not_assessed`), `matrix_reason`, `matrix_present`,
#'   `matrix_keys_matched` (all matching keys, joined by `"; "`), `matrix_key`
#'   (the single key read, else `NA`), `matrix_n` (matrix size), `matrix_max_dev`
#'   and `matrix_identity`; for the flag `apply_status`, `apply_reason`,
#'   `apply_present`, `apply_keys_matched`, `apply_key`, `apply_raw` and
#'   `apply_parsed`, then `status` (`"assessed"`, or `"not_assessed"` when the
#'   header, the matrix or the flag is) and `reason` (the reasons of the parts
#'   not assessed, joined by `"; "`, else `NA`). With `expected_state` it adds `expected_<item>` and
#'   `<item>_expected_vs_observed` for each item given, with values `same`,
#'   `differs`, `not_expected` (no expectation) or `not_assessed` (observed
#'   value is `NA`). Attribute `provenance` holds package versions, R version,
#'   date, `seed` (`NULL`, nothing is random here) and the arguments used.
#'
#' @seealso [sw_fcs_header_audit()]
#'
#' @examples
#' \dontrun{
#' res <- sw_fcs_compensation_state(
#'   list.files("fcs", pattern = "\\.fcs$", full.names = TRUE),
#'   matrix_keywords = c("$SPILLOVER", "$SPILL", "$COMP"),
#'   apply_keywords = c("$APPLY_COMP"), # illustrative; use your vendor's key
#'   apply_true = c("true"), # illustrative; use your vendor's words
#'   apply_false = c("false"),
#'   zero_tol = 1e-6, # illustrative, choose your own
#'   expected_state = NULL
#' )
#' res[, c("file", "matrix_present", "matrix_identity", "apply_parsed")]
#' }
#'
#' @export
sw_fcs_compensation_state <- function(
  paths,
  matrix_keywords,
  apply_keywords,
  apply_true,
  apply_false,
  zero_tol,
  expected_state
) {
  .sw_stop_missing(c(
    paths = missing(paths),
    matrix_keywords = missing(matrix_keywords),
    apply_keywords = missing(apply_keywords),
    apply_true = missing(apply_true),
    apply_false = missing(apply_false),
    zero_tol = missing(zero_tol),
    expected_state = missing(expected_state)
  ))
  if (!is.character(paths) || length(paths) == 0 || anyNA(paths)) {
    stop("`paths` must be a non-empty character vector.", call. = FALSE)
  }
  labels <- basename(paths)
  if (anyDuplicated(labels)) {
    stop("File names in `paths` must be unique.", call. = FALSE)
  }
  for (arg in c("matrix_keywords", "apply_keywords")) {
    v <- get(arg)
    if (!is.character(v) || length(v) == 0 || anyNA(v) || !all(nzchar(v))) {
      stop(
        sprintf("`%s` must be a non-empty character vector.", arg),
        call. = FALSE
      )
    }
    if (anyDuplicated(v)) {
      stop(sprintf("`%s` must hold unique keys.", arg), call. = FALSE)
    }
  }
  for (arg in c("apply_true", "apply_false")) {
    v <- get(arg)
    if (!is.character(v) || length(v) == 0 || anyNA(v) || !all(nzchar(v))) {
      stop(
        sprintf("`%s` must be a non-empty character vector.", arg),
        call. = FALSE
      )
    }
  }
  apply_true <- tolower(trimws(apply_true))
  apply_false <- tolower(trimws(apply_false))
  if (length(both <- intersect(apply_true, apply_false))) {
    stop(
      "`apply_true` and `apply_false` both hold: ",
      paste(both, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (
    !is.numeric(zero_tol) ||
      length(zero_tol) != 1 ||
      !is.finite(zero_tol) ||
      zero_tol < 0
  ) {
    stop("`zero_tol` must be one non-negative number.", call. = FALSE)
  }
  exp_items <- .sw_comp_check_expected(expected_state)
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required.", call. = FALSE)
  }

  rows <- lapply(seq_along(paths), function(i) {
    .sw_comp_one(
      paths[i],
      labels[i],
      matrix_keywords,
      apply_keywords,
      apply_true,
      apply_false,
      zero_tol
    )
  })
  # Files only named in expected_state stay in the table, with a reason.
  if (!is.null(exp_items)) {
    ghosts <- setdiff(expected_state$file, labels)
    rows <- c(
      rows,
      lapply(ghosts, function(f) {
        .sw_comp_blank(
          f,
          NA_character_,
          "named in `expected_state` but not in `paths`"
        )
      })
    )
  }
  out <- tibble::as_tibble(dplyr::bind_rows(rows))
  out <- .sw_comp_add_status(out)
  if (!is.null(exp_items)) {
    out <- .sw_comp_add_expected(out, expected_state, exp_items)
  }

  attr(out, "provenance") <- .sw_import_provenance(
    args = list(
      matrix_keywords = matrix_keywords,
      apply_keywords = apply_keywords,
      apply_true = apply_true,
      apply_false = apply_false,
      zero_tol = zero_tol,
      expected_state_given = !is.null(expected_state)
    ),
    seed = NULL,
    extra = list(packages = .sw_import_packages())
  )
  class(out) <- c("sw_fcs_compensation_state", class(out))
  out
}

# One status and reason per file, from the header, matrix and flag parts:
# not_assessed when any part is (not_present still counts as assessed).
.sw_comp_add_status <- function(out) {
  part_bad <- cbind(
    out$header_status == "not_assessed",
    out$matrix_status == "not_assessed",
    out$apply_status == "not_assessed"
  )
  reasons <- cbind(out$header_reason, out$matrix_reason, out$apply_reason)
  reasons[!part_bad] <- NA_character_
  out$status <- dplyr::if_else(
    rowSums(part_bad) > 0,
    "not_assessed",
    "assessed"
  )
  out$reason <- apply(reasons, 1, \(r) {
    r <- unique(r[!is.na(r)])
    if (length(r)) paste(r, collapse = "; ") else NA_character_
  })
  out
}

# Items an expected_state table may hold, and the output column each compares.
.sw_comp_items <- c(
  matrix_present = "matrix_present",
  matrix_identity = "matrix_identity",
  apply_parsed = "apply_parsed"
)

# Check expected_state; return the item columns it holds, or NULL when skipped.
.sw_comp_check_expected <- function(expected_state) {
  if (is.null(expected_state)) {
    return(NULL)
  }
  if (!is.data.frame(expected_state) || !"file" %in% names(expected_state)) {
    stop(
      "`expected_state` must be NULL or a data frame with a `file` column.",
      call. = FALSE
    )
  }
  items <- intersect(names(.sw_comp_items), names(expected_state))
  if (length(items) == 0) {
    stop(
      "`expected_state` needs at least one of: ",
      paste(names(.sw_comp_items), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.character(expected_state$file) || anyNA(expected_state$file)) {
    stop("`expected_state$file` must be file names, no NA.", call. = FALSE)
  }
  if (anyDuplicated(expected_state$file)) {
    stop("Each file may appear once in `expected_state`.", call. = FALSE)
  }
  for (it in items) {
    if (!is.logical(expected_state[[it]])) {
      stop(
        sprintf(
          "`expected_state$%s` must be logical (NA = no expectation).",
          it
        ),
        call. = FALSE
      )
    }
  }
  items
}

# Add expected_<item> and <item>_expected_vs_observed columns.
.sw_comp_add_expected <- function(out, expected_state, items) {
  at <- match(out$file, expected_state$file)
  for (it in items) {
    want <- expected_state[[it]][at]
    got <- out[[.sw_comp_items[[it]]]]
    out[[paste0("expected_", it)]] <- want
    out[[paste0(it, "_expected_vs_observed")]] <- dplyr::case_when(
      is.na(want) ~ "not_expected",
      is.na(got) ~ "not_assessed",
      want == got ~ "same",
      .default = "differs"
    )
  }
  out
}

# One output row for a file whose header could not be used.
.sw_comp_blank <- function(file, path, reason) {
  tibble::tibble(
    file = file,
    path = path,
    header_status = "not_assessed",
    header_reason = reason,
    matrix_status = "not_assessed",
    matrix_reason = "header not read",
    matrix_present = NA,
    matrix_keys_matched = NA_character_,
    matrix_key = NA_character_,
    matrix_n = NA_real_,
    matrix_max_dev = NA_real_,
    matrix_identity = NA,
    apply_status = "not_assessed",
    apply_reason = "header not read",
    apply_present = NA,
    apply_keys_matched = NA_character_,
    apply_key = NA_character_,
    apply_raw = NA_character_,
    apply_parsed = NA
  )
}

# Read one header and fill one row.
.sw_comp_one <- function(
  path,
  label,
  matrix_keywords,
  apply_keywords,
  apply_true,
  apply_false,
  zero_tol
) {
  # flowCore leaves a connection open on a non-FCS file, so check the magic bytes.
  magic <- tryCatch(
    readBin(path, what = "raw", n = 3L),
    error = function(e) raw(0),
    warning = function(w) raw(0)
  )
  if (!identical(rawToChar(magic), "FCS")) {
    return(.sw_comp_blank(label, path, .sw_fcs_read_reason(path)))
  }
  # A malformed file past the magic bytes can still leave a connection open;
  # close any connection that this read opened.
  conns_before <- getAllConnections()
  h <- tryCatch(
    suppressWarnings(unlist(flowCore::read.FCSheader(path)[[1]])),
    error = function(e) e
  )
  for (id in setdiff(getAllConnections(), conns_before)) {
    try(close(getConnection(id)), silent = TRUE)
  }
  if (inherits(h, "error")) {
    return(.sw_comp_blank(label, path, .sw_fcs_read_reason(path, h)))
  }
  if (!is.character(h) || !"$PAR" %in% names(h)) {
    return(.sw_comp_blank(
      label,
      path,
      "header unreadable: no FCS keywords found"
    ))
  }
  h <- stats::setNames(trimws(unname(h)), names(h))

  row <- .sw_comp_blank(label, path, NA_character_)
  row$header_status <- "read"

  # Matrix: present, which key, parsed, distance from identity.
  hits <- unique(matrix_keywords[matrix_keywords %in% names(h)])
  row$matrix_present <- length(hits) > 0
  row$matrix_keys_matched <- if (length(hits)) {
    paste(hits, collapse = "; ")
  } else {
    NA_character_
  }
  if (length(hits) == 0) {
    row$matrix_status <- "not_present"
    row$matrix_reason <- NA_character_
  } else if (length(hits) > 1) {
    row$matrix_reason <- paste0(
      "several matrix keys matched: ",
      row$matrix_keys_matched
    )
  } else {
    row$matrix_key <- hits
    parsed <- .sw_parse_spillover(h[[hits]])
    if (parsed$ok) {
      dev <- max(abs(parsed$matrix - diag(parsed$n)))
      row$matrix_status <- "assessed"
      row$matrix_reason <- NA_character_
      row$matrix_n <- parsed$n
      row$matrix_max_dev <- dev
      row$matrix_identity <- dev <= zero_tol
    } else {
      row$matrix_reason <- parsed$reason
    }
  }

  # Apply flag: raw text and parsed value.
  hits <- unique(apply_keywords[apply_keywords %in% names(h)])
  row$apply_present <- length(hits) > 0
  row$apply_keys_matched <- if (length(hits)) {
    paste(hits, collapse = "; ")
  } else {
    NA_character_
  }
  if (length(hits) == 0) {
    row$apply_status <- "not_present"
    row$apply_reason <- NA_character_
  } else if (length(hits) > 1) {
    row$apply_reason <- paste0(
      "several apply keys matched: ",
      row$apply_keys_matched
    )
  } else {
    row$apply_key <- hits
    row$apply_raw <- h[[hits]]
    flag <- .sw_parse_flag(row$apply_raw, apply_true, apply_false)
    row$apply_parsed <- flag
    if (is.na(flag)) {
      row$apply_reason <- sprintf(
        "flag value '%s' is not a recognised true or false",
        row$apply_raw
      )
    } else {
      row$apply_status <- "assessed"
      row$apply_reason <- NA_character_
    }
  }
  row
}

# Text to TRUE / FALSE using the caller's word lists, NA when not listed.
.sw_parse_flag <- function(x, true_words, false_words) {
  v <- tolower(trimws(x))
  if (v %in% tolower(trimws(true_words))) {
    return(TRUE)
  }
  if (v %in% tolower(trimws(false_words))) {
    return(FALSE)
  }
  NA
}

# Parse an FCS spillover string. Returns list(ok, reason, n, names, matrix);
# `matrix` and `names` are NULL when not ok (or no names in the old layout).
.sw_parse_spillover <- function(x) {
  fail <- function(reason) {
    list(ok = FALSE, reason = reason, n = NA_real_, names = NULL, matrix = NULL)
  }
  if (length(x) != 1 || is.na(x) || !nzchar(trimws(x))) {
    return(fail("matrix value is empty"))
  }
  tok <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  # strsplit() drops a final empty field; a trailing comma is a malformed field.
  if (endsWith(trimws(x), ",")) {
    tok <- c(tok, "")
  }
  n <- suppressWarnings(as.numeric(tok[1]))
  if (!is.finite(n) || n < 1 || n != round(n)) {
    return(fail(sprintf(
      "first field '%s' is not a whole number of channels",
      tok[1]
    )))
  }
  named_len <- 1 + n + n^2
  plain_len <- 1 + n^2
  if (length(tok) == named_len) {
    nm <- tok[2:(n + 1)]
    vals <- tok[(n + 2):length(tok)]
    if (any(!nzchar(nm))) {
      return(fail("a channel name is empty"))
    }
    if (anyDuplicated(nm)) {
      return(fail("a channel name is repeated"))
    }
  } else if (length(tok) == plain_len) {
    nm <- NULL
    vals <- tok[-1]
  } else {
    return(fail(sprintf(
      "expected %.0f values (n, %.0f names, %.0f numbers) or %.0f (no names), got %d",
      named_len,
      n,
      n^2,
      plain_len,
      length(tok)
    )))
  }
  num <- suppressWarnings(as.numeric(vals))
  if (any(!is.finite(num))) {
    bad <- which(!is.finite(num))[1]
    return(fail(sprintf(
      "matrix value %d ('%s') is not a finite number",
      bad,
      vals[bad]
    )))
  }
  list(
    ok = TRUE,
    reason = NA_character_,
    n = n,
    names = nm,
    matrix = matrix(
      num,
      nrow = n,
      ncol = n,
      byrow = TRUE,
      dimnames = if (!is.null(nm)) list(nm, nm)
    )
  )
}
