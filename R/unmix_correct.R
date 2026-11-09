#' @title Thin Wrappers for AutoSpectral Unmixing-Error Correction
#'
#' @description
#' Wrappers around two AutoSpectral functions that correct reference spectra
#' from stained samples. `sw_unmix_correct_signatures()` wraps
#' `AutoSpectral::correct.unmixing.signatures()`. `sw_unmix_correct_spillover()`
#' wraps `AutoSpectral::fix.my.unmix()`. `sw_unmix_correct()` runs them in
#' AutoSpectral's stated order. No AutoSpectral code is copied here.
#'
#' @details
#' Signatures and return values were read from AutoSpectral version 1.8.8.
#' AutoSpectral is AGPL licensed and sits in Suggests.
#'
#' The corrected spectra are a candidate. They are not a known improvement.
#' The AutoSpectral documentation calls error correction fairly speculative,
#' and it is not adaptive unmixing. Score the corrected spectra against a run
#' of the unmodified AutoSpectral baseline on the same events, controls and
#' panel before you use them (ROADMAP.md, "Gold standard").
#'
#' Orientation follows the package convention. Spectra sit in rows and
#' detectors sit in columns. Row names and column names must be present and
#' unique. Input is checked before the dependency check, so argument errors
#' show up even when AutoSpectral is missing.
#'
#' Thresholds are never set here. Every tuning value goes through `...` under
#' AutoSpectral's own name, and AutoSpectral's default applies when you leave
#' it out.
#'
#' @name unmix_correct
#' @keywords internal
NULL


# ---- Internal helpers --------------------------------------------------------

.sw_correct_require_autospectral <- function() {
  if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
    stop("Package 'AutoSpectral' is required.", call. = FALSE)
  }
  invisible(TRUE)
}

# Call an exported AutoSpectral function with a named argument list.
.sw_correct_call <- function(fn, args) {
  do.call(getExportedValue("AutoSpectral", fn), args)
}

.sw_correct_version <- function() {
  as.character(getNamespaceVersion("AutoSpectral"))
}

# Numeric matrix, fluorophores in rows, detectors in columns, unique names.
.sw_correct_check_spectra <- function(x, arg) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(
      sprintf(
        "`%s` must be a numeric matrix (fluorophores in rows, detectors in columns).",
        arg
      ),
      call. = FALSE
    )
  }
  rn <- rownames(x)
  cn <- colnames(x)
  if (is.null(rn) || is.null(cn) || any(!nzchar(rn)) || any(!nzchar(cn))) {
    stop(
      sprintf(
        "`%s` must have row names (fluorophores) and column names (detectors).",
        arg
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(rn) > 0L || anyDuplicated(cn) > 0L) {
    stop(
      sprintf("`%s` must have unique row names and unique column names.", arg),
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", arg), call. = FALSE)
  }
  invisible(x)
}

.sw_correct_check_files <- function(x, arg) {
  if (!is.character(x) || length(x) == 0L || anyNA(x) || any(!nzchar(x))) {
    stop(
      sprintf("`%s` must be a character vector of file paths.", arg),
      call. = FALSE
    )
  }
  missing_files <- x[!file.exists(x)]
  if (length(missing_files) > 0L) {
    stop(
      sprintf(
        "`%s` file(s) not found: %s",
        arg,
        paste(missing_files, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(x)
}

# Dots must be named and must not repeat an argument the wrapper sets.
.sw_correct_check_dots <- function(dots, reserved) {
  if (
    length(dots) > 0L && (is.null(names(dots)) || any(!nzchar(names(dots))))
  ) {
    stop(
      "Arguments in `...` must be named (AutoSpectral's own names).",
      call. = FALSE
    )
  }
  clash <- intersect(names(dots), reserved)
  if (length(clash) > 0L) {
    stop(
      sprintf(
        "Do not pass %s through `...`; use the wrapper's own argument.",
        paste0("`", clash, "`", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(dots)
}

# File-based entry needs an existing unstained sample and a flow.control list.
.sw_correct_check_file_args <- function(dots) {
  if (is.null(dots[["unstained.sample"]])) {
    stop("File input needs `unstained.sample` in `...`.", call. = FALSE)
  }
  .sw_correct_check_files(dots[["unstained.sample"]], "unstained.sample")
  if (length(dots[["unstained.sample"]]) != 1L) {
    stop("`unstained.sample` must be a single file path.", call. = FALSE)
  }
  if (!is.list(dots[["flow.control"]])) {
    stop("File input needs `flow.control` (a list) in `...`.", call. = FALSE)
  }
  invisible(dots)
}

# Named logical over fluorophores, aligned to `fluorophores` by name.
.sw_correct_align_accepted <- function(accepted, fluorophores) {
  out <- rep(NA, length(fluorophores))
  names(out) <- fluorophores
  if (!is.null(accepted) && !is.null(names(accepted))) {
    hit <- intersect(names(accepted), fluorophores)
    out[hit] <- as.logical(accepted[hit])
  }
  out
}

# Build the sw_unmix_correction object from one or more step results.
.sw_correct_package <- function(
  original,
  corrected,
  accepted,
  steps,
  raw,
  call_args
) {
  changes <- sw_unmix_compare_spectra(original, corrected)
  changes[["accepted"]] <- unname(accepted[changes[["fluorophore"]]])
  changes[["steps_run"]] <- paste(steps, collapse = "+")
  structure(
    list(
      spectra = corrected[rownames(original), colnames(original), drop = FALSE],
      original = original,
      changes = changes,
      fit_logs = list(autospectral = raw),
      note = "rerun spectral variants with corrected spectra",
      provenance = list(
        autospectral_version = .sw_correct_version(),
        steps = steps,
        call_args = call_args
      )
    ),
    class = "sw_unmix_correction"
  )
}

# Run correct.unmixing.signatures(); return the raw list and accepted.
.sw_correct_run_signatures <- function(raw_data_or_files, spectra, asp, dots) {
  args <- c(list(spectra = spectra, asp = asp), dots)
  if (is.matrix(raw_data_or_files)) {
    args[["raw.data"]] <- raw_data_or_files
  } else {
    args[["fully.stained.sample"]] <- raw_data_or_files
  }
  raw <- .sw_correct_call("correct.unmixing.signatures", args)
  list(
    raw = raw,
    spectra = raw[["spectra"]],
    accepted = .sw_correct_align_accepted(raw[["accepted"]], rownames(spectra)),
    arg_names = names(args)
  )
}

# Run fix.my.unmix(); accepted comes from signature.log by fluorophore name.
.sw_correct_run_spillover <- function(
  stained_files_or_data,
  spectra,
  asp,
  dots
) {
  args <- c(
    list(
      spectra = spectra,
      fully.stained.sample = stained_files_or_data,
      asp = asp
    ),
    dots
  )
  raw <- .sw_correct_call("fix.my.unmix", args)
  sig_log <- raw[["signature.log"]]
  accepted <- NULL
  if (
    !is.null(sig_log) && all(c("fluorophore", "accepted") %in% names(sig_log))
  ) {
    accepted <- stats::setNames(sig_log[["accepted"]], sig_log[["fluorophore"]])
  }
  list(
    raw = raw,
    spectra = raw[["spectra"]],
    accepted = .sw_correct_align_accepted(accepted, rownames(spectra)),
    arg_names = names(args)
  )
}

.sw_correct_validate_signatures <- function(
  raw_data_or_files,
  spectra,
  asp,
  dots
) {
  .sw_correct_check_spectra(spectra, "spectra")
  if (!is.null(asp) && !is.list(asp)) {
    stop(
      "`asp` must be a list (from get.autospectral.param()) or NULL.",
      call. = FALSE
    )
  }
  .sw_correct_check_dots(
    dots,
    c("spectra", "asp", "raw.data", "fully.stained.sample")
  )
  if (is.matrix(raw_data_or_files)) {
    if (
      !is.numeric(raw_data_or_files) || is.null(colnames(raw_data_or_files))
    ) {
      stop(
        "`raw_data_or_files` as a matrix must be numeric with detector column names.",
        call. = FALSE
      )
    }
    if (!identical(colnames(raw_data_or_files), colnames(spectra))) {
      stop(
        "`raw_data_or_files` columns must equal `colnames(spectra)`, in the same order.",
        call. = FALSE
      )
    }
  } else if (is.character(raw_data_or_files)) {
    .sw_correct_check_files(raw_data_or_files, "raw_data_or_files")
    if (length(raw_data_or_files) != 1L) {
      stop(
        "correct.unmixing.signatures() reads one fully stained file; got several.",
        call. = FALSE
      )
    }
    .sw_correct_check_file_args(dots)
  } else {
    stop(
      "`raw_data_or_files` must be a numeric matrix (events x detectors) or a file path.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.sw_correct_validate_spillover <- function(
  stained_files_or_data,
  spectra,
  asp,
  dots
) {
  .sw_correct_check_spectra(spectra, "spectra")
  if (!is.list(asp)) {
    stop("`asp` must be a list (from get.autospectral.param()).", call. = FALSE)
  }
  .sw_correct_check_dots(dots, c("spectra", "asp", "fully.stained.sample"))
  if (!is.character(stained_files_or_data)) {
    stop(
      "`stained_files_or_data` must be file paths: fix.my.unmix() reads files only.",
      call. = FALSE
    )
  }
  .sw_correct_check_files(stained_files_or_data, "stained_files_or_data")
  .sw_correct_check_file_args(dots)
  invisible(TRUE)
}


# ---- Compare -----------------------------------------------------------------

#' Per-fluorophore change between two spectra matrices
#'
#' Compares an original and a corrected spectra matrix row by row. Rows and
#' columns are aligned by name, never by position.
#'
#' @details
#' Both matrices hold fluorophores in rows and detectors in columns. They must
#' have the same row names and the same column names, in any order. A missing
#' or extra name is an error. `corrected` is reordered to match `original`.
#'
#' `angle_deg` is the angle between the two rows, from the cosine of the raw
#' rows (`acos`, clamped to `[-1, 1]`). Identical rows give exactly 0.
#' `max_abs_delta` is the largest absolute detector difference after each row
#' is divided by its own peak absolute value. `changed` is `angle_deg > 0`, so
#' a row that is only rescaled counts as unchanged.
#'
#' @param original Numeric matrix of starting spectra.
#' @param corrected Numeric matrix of corrected spectra.
#'
#' @return A tibble with one row per fluorophore of `original` and columns
#'   `fluorophore`, `angle_deg`, `max_abs_delta` and `changed`.
#'
#' @examples
#' orig <- matrix(c(1, 0, 0, 1), 2, dimnames = list(c("A", "B"), c("d1", "d2")))
#' corr <- orig[c("B", "A"), ]
#' corr["A", ] <- c(1, 0.2)
#' sw_unmix_compare_spectra(orig, corr)
#'
#' @export
sw_unmix_compare_spectra <- function(original, corrected) {
  .sw_correct_check_spectra(original, "original")
  .sw_correct_check_spectra(corrected, "corrected")
  row_diff <- union(
    setdiff(rownames(original), rownames(corrected)),
    setdiff(rownames(corrected), rownames(original))
  )
  col_diff <- union(
    setdiff(colnames(original), colnames(corrected)),
    setdiff(colnames(corrected), colnames(original))
  )
  if (length(row_diff) > 0L) {
    stop(
      sprintf(
        "`original` and `corrected` fluorophore names differ: %s",
        paste(row_diff, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(col_diff) > 0L) {
    stop(
      sprintf(
        "`original` and `corrected` detector names differ: %s",
        paste(col_diff, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  corrected <- corrected[rownames(original), colnames(original), drop = FALSE]

  angle <- vapply(
    rownames(original),
    \(f) {
      a <- original[f, ]
      b <- corrected[f, ]
      if (identical(a, b)) {
        return(0)
      }
      cos_ab <- sum(a * b) / sqrt(sum(a * a) * sum(b * b))
      acos(min(1, max(-1, cos_ab))) * 180 / pi
    },
    numeric(1)
  )
  peak_norm <- function(m) m / apply(abs(m), 1, max)
  delta <- apply(abs(peak_norm(original) - peak_norm(corrected)), 1, max)

  tibble::tibble(
    fluorophore = rownames(original),
    angle_deg = unname(angle),
    max_abs_delta = unname(delta),
    changed = unname(angle) > 0
  )
}


# ---- Signatures --------------------------------------------------------------

#' Correct spectrum shape errors with AutoSpectral
#'
#' Wraps `AutoSpectral::correct.unmixing.signatures()`. It estimates shape
#' errors in each fluorophore row from stained data and updates rows that pass
#' AutoSpectral's acceptance gates.
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8. AutoSpectral states
#' the order: run this first, then `sw_unmix_correct_spillover()` on the
#' spectra it returns.
#'
#' The corrected spectra are a candidate to test against the unmodified
#' AutoSpectral baseline (ROADMAP.md, "Gold standard"). This wrapper makes no
#' claim that they improve unmixing.
#'
#' Input. A numeric matrix (events x detectors) goes to `raw.data`. Its column
#' names must equal `colnames(spectra)`. A single file path goes to
#' `fully.stained.sample`; then `unstained.sample` and `flow.control` are
#' required in `...`. `asp` is optional in AutoSpectral (used only for
#' `asp$bird.seed`), so `NULL` is allowed.
#'
#' Thresholds keep AutoSpectral's defaults unless passed in `...`. In 1.8.8
#' these include `min.events = 200L`, `min.span = 5`, `min.explained = 0.5`,
#' `min.gain = 0.002`, `max.step = 0.15`, `max.span.drift = 1.10`,
#' `max.bg.alignment = -0.9`, `unstained.threshold = 0.99`,
#' `unstained.margin = 1.3`, `n.iter = 6L` and `bg.mode = "global.mean"`.
#' See the AutoSpectral help page for the full list.
#'
#' @param raw_data_or_files Numeric matrix (events x detectors) or one fully
#'   stained FCS file path.
#' @param spectra Numeric matrix of starting spectra (fluorophores x
#'   detectors), L-infinity normalised.
#' @param asp AutoSpectral parameter list from
#'   `AutoSpectral::get.autospectral.param()`, or `NULL`.
#' @param ... Passed to `AutoSpectral::correct.unmixing.signatures()` under
#'   its own argument names (for example `variants`, `af.name`,
#'   `unstained.sample`, `flow.control`, `min.gain`).
#'
#' @return An `sw_unmix_correction` object: a list with `spectra` (corrected),
#'   `original`, `changes` (tibble from [sw_unmix_compare_spectra()] plus
#'   `accepted`, from AutoSpectral's `accepted`, and `steps_run`), `fit_logs`
#'   (the full AutoSpectral return under `autospectral$signatures`), `note`
#'   and `provenance`.
#'
#' @examples
#' \dontrun{
#' res <- sw_unmix_correct_signatures(
#'   "stained.fcs",
#'   spectra,
#'   asp,
#'   unstained.sample = "unstained.fcs",
#'   flow.control = flow_control,
#'   variants = variants
#' )
#' res$changes
#' }
#'
#' @export
sw_unmix_correct_signatures <- function(raw_data_or_files, spectra, asp, ...) {
  dots <- list(...)
  .sw_correct_validate_signatures(raw_data_or_files, spectra, asp, dots)
  .sw_correct_require_autospectral()

  run <- .sw_correct_run_signatures(raw_data_or_files, spectra, asp, dots)
  .sw_correct_package(
    original = spectra,
    corrected = run[["spectra"]],
    accepted = run[["accepted"]],
    steps = "signatures",
    raw = list(signatures = run[["raw"]]),
    call_args = list(signatures = run[["arg_names"]])
  )
}


# ---- Spillover ---------------------------------------------------------------

#' Correct in-span spillover error with AutoSpectral
#'
#' Wraps `AutoSpectral::fix.my.unmix()`. It estimates the residual spillover
#' matrix from fully stained samples and re-measures rows that pass
#' AutoSpectral's gates.
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8. AutoSpectral states
#' the order: run `sw_unmix_correct_signatures()` first and pass its spectra
#' here. `fix.my.unmix()` does not recompute the spillover spread under the
#' corrected spectra. Re-extract spectral variants with the corrected spectra
#' before using them.
#'
#' The corrected spectra are a candidate to test against the unmodified
#' AutoSpectral baseline (ROADMAP.md, "Gold standard"). This wrapper makes no
#' claim that they improve unmixing.
#'
#' Input. `fix.my.unmix()` reads files only in 1.8.8; it has no in-memory
#' entry point. `stained_files_or_data` must be one or more fully stained FCS
#' paths. `unstained.sample` and `flow.control` are required in `...`, and
#' `asp` must be a list.
#'
#' Side effects. AutoSpectral's defaults `figures = TRUE` and `save = TRUE`
#' write a spillover heatmap and CSV files to `asp$fix.unmixing.dir`. Pass
#' `figures = FALSE` and `save = FALSE` in `...` to stop that.
#'
#' Thresholds keep AutoSpectral's defaults unless passed in `...`. In 1.8.8
#' these include `unstained.threshold = 0.95`, `unstained.margin = 1.3`,
#' `spread.kappa = 3`, `max.iter = 20L`, `max.coefficient = 0.5`,
#' `convergence.threshold = 0.003`, `max.angle = 30`, `max.resid = 0.03`,
#' `min.explained = 0.8`, `smooth.updates = TRUE`, `smooth.rounds = 5L` and
#' `bg.mode = "per.cell"` (which needs `af.spectra`). See the AutoSpectral
#' help page for the full list.
#'
#' @param stained_files_or_data Character vector of one or more fully stained
#'   FCS file paths.
#' @param spectra Numeric matrix of starting spectra (fluorophores x
#'   detectors), L-infinity normalised.
#' @param asp AutoSpectral parameter list from
#'   `AutoSpectral::get.autospectral.param()`.
#' @param ... Passed to `AutoSpectral::fix.my.unmix()` under its own argument
#'   names (for example `unstained.sample`, `flow.control`, `variants`,
#'   `af.spectra`, `figures`, `save`).
#'
#' @return An `sw_unmix_correction` object, as for
#'   [sw_unmix_correct_signatures()]. `accepted` comes from the
#'   `fluorophore` and `accepted` columns of AutoSpectral's `signature.log`.
#'   The full AutoSpectral return sits in `fit_logs$autospectral$spillover`.
#'
#' @examples
#' \dontrun{
#' res <- sw_unmix_correct_spillover(
#'   c("stained_1.fcs", "stained_2.fcs"),
#'   spectra,
#'   asp,
#'   unstained.sample = "unstained.fcs",
#'   flow.control = flow_control,
#'   variants = variants,
#'   af.spectra = af_spectra,
#'   figures = FALSE
#' )
#' }
#'
#' @export
sw_unmix_correct_spillover <- function(
  stained_files_or_data,
  spectra,
  asp,
  ...
) {
  dots <- list(...)
  .sw_correct_validate_spillover(stained_files_or_data, spectra, asp, dots)
  .sw_correct_require_autospectral()

  run <- .sw_correct_run_spillover(stained_files_or_data, spectra, asp, dots)
  .sw_correct_package(
    original = spectra,
    corrected = run[["spectra"]],
    accepted = run[["accepted"]],
    steps = "spillover",
    raw = list(spillover = run[["raw"]]),
    call_args = list(spillover = run[["arg_names"]])
  )
}


# ---- Two-step convenience ----------------------------------------------------

#' Run AutoSpectral's two correction steps in order
#'
#' Runs `sw_unmix_correct_signatures()` and then
#' `sw_unmix_correct_spillover()` on the spectra the first step returns.
#'
#' @details
#' `steps` has no default. Name the steps you want. Allowed values are
#' `"signatures"` and `"spillover"`. With both, the order must be
#' `c("signatures", "spillover")`, the order AutoSpectral states; the reverse
#' is an error.
#'
#' Each step takes its own argument list. The two AutoSpectral functions
#' share argument names with different meanings and defaults (for example
#' `bg.mode`, `unstained.threshold`, `min.explained`), so one shared `...`
#' would be wrong for one of them.
#'
#' With both steps, `stained` must be a single fully stained file, because
#' `correct.unmixing.signatures()` reads one file and `fix.my.unmix()` reads
#' files only. Run the two wrappers by hand to use a matrix for step one.
#'
#' In `changes`, `angle_deg` and `max_abs_delta` compare the input spectra with
#' the final spectra. `accepted` is `TRUE` when any step accepted the row,
#' `FALSE` when every reporting step rejected it, and `NA` when no step
#' reported it. Per-step AutoSpectral results stay in `fit_logs$autospectral`.
#'
#' `fix.my.unmix()` does not recompute the spillover spread. Re-extract
#' spectral variants with the corrected spectra. The result is a candidate to
#' test against the unmodified AutoSpectral baseline (ROADMAP.md, "Gold
#' standard"); no improvement is claimed.
#'
#' @param stained Input for the first step run: see
#'   [sw_unmix_correct_signatures()] and [sw_unmix_correct_spillover()].
#' @param spectra Numeric matrix of starting spectra (fluorophores x
#'   detectors).
#' @param asp AutoSpectral parameter list.
#' @param steps Character, required. One or both of `"signatures"` and
#'   `"spillover"`.
#' @param signatures_args Named list passed to
#'   `AutoSpectral::correct.unmixing.signatures()`.
#' @param spillover_args Named list passed to `AutoSpectral::fix.my.unmix()`.
#'
#' @return An `sw_unmix_correction` object.
#'
#' @examples
#' \dontrun{
#' res <- sw_unmix_correct(
#'   "stained.fcs",
#'   spectra,
#'   asp,
#'   steps = c("signatures", "spillover"),
#'   signatures_args = list(
#'     unstained.sample = "unstained.fcs",
#'     flow.control = flow_control,
#'     variants = variants
#'   ),
#'   spillover_args = list(
#'     unstained.sample = "unstained.fcs",
#'     flow.control = flow_control,
#'     variants = variants,
#'     af.spectra = af_spectra,
#'     figures = FALSE,
#'     save = FALSE
#'   )
#' )
#' print(res)
#' }
#'
#' @export
sw_unmix_correct <- function(
  stained,
  spectra,
  asp,
  steps,
  signatures_args = list(),
  spillover_args = list()
) {
  allowed <- c("signatures", "spillover")
  if (missing(steps)) {
    stop(
      "`steps` is required: name \"signatures\", \"spillover\" or both.",
      call. = FALSE
    )
  }
  if (
    !is.character(steps) ||
      length(steps) == 0L ||
      anyNA(steps) ||
      anyDuplicated(steps) > 0L ||
      !all(steps %in% allowed)
  ) {
    stop(
      "`steps` must be unique values from \"signatures\" and \"spillover\".",
      call. = FALSE
    )
  }
  if (identical(steps, rev(allowed))) {
    stop(
      "Run \"signatures\" before \"spillover\": use steps = c(\"signatures\", \"spillover\").",
      call. = FALSE
    )
  }
  if (!is.list(signatures_args) || !is.list(spillover_args)) {
    stop("`signatures_args` and `spillover_args` must be lists.", call. = FALSE)
  }
  if ("signatures" %in% steps) {
    .sw_correct_validate_signatures(stained, spectra, asp, signatures_args)
  }
  if ("spillover" %in% steps) {
    .sw_correct_validate_spillover(stained, spectra, asp, spillover_args)
  }
  .sw_correct_require_autospectral()

  current <- spectra
  accepted <- list()
  raw <- list()
  call_args <- list()
  for (step in steps) {
    run <- if (step == "signatures") {
      .sw_correct_run_signatures(stained, current, asp, signatures_args)
    } else {
      .sw_correct_run_spillover(stained, current, asp, spillover_args)
    }
    current <- run[["spectra"]]
    accepted[[step]] <- run[["accepted"]][rownames(spectra)]
    raw[[step]] <- run[["raw"]]
    call_args[[step]] <- run[["arg_names"]]
  }

  # Any step accepted -> TRUE; all reporting steps rejected -> FALSE; else NA.
  acc_mat <- do.call(cbind, accepted)
  acc <- apply(acc_mat, 1, \(x) if (all(is.na(x))) NA else any(x, na.rm = TRUE))
  names(acc) <- rownames(spectra)

  .sw_correct_package(
    original = spectra,
    corrected = current,
    accepted = acc,
    steps = steps,
    raw = raw,
    call_args = call_args
  )
}


# ---- Print -------------------------------------------------------------------

#' @export
print.sw_unmix_correction <- function(x, ...) {
  ch <- x[["changes"]]
  cat("<sw_unmix_correction>\n")
  cat(
    "Steps: ",
    paste(x[["provenance"]][["steps"]], collapse = " -> "),
    " (AutoSpectral ",
    x[["provenance"]][["autospectral_version"]],
    ")\n",
    sep = ""
  )
  cat(
    sprintf(
      "Fluorophores: %d; changed: %d; accepted: %d; accepted NA: %d\n",
      nrow(ch),
      sum(ch[["changed"]]),
      sum(ch[["accepted"]], na.rm = TRUE),
      sum(is.na(ch[["accepted"]]))
    )
  )
  moved <- ch[ch[["changed"]], , drop = FALSE]
  if (nrow(moved) > 0L) {
    moved <- moved[order(-moved[["angle_deg"]]), , drop = FALSE]
    cat("Largest changes (angle_deg, max_abs_delta):\n")
    shown <- moved[seq_len(min(10L, nrow(moved))), , drop = FALSE]
    for (i in seq_len(nrow(shown))) {
      cat(
        sprintf(
          "  %s: %.3f, %.4f\n",
          shown[["fluorophore"]][i],
          shown[["angle_deg"]][i],
          shown[["max_abs_delta"]][i]
        )
      )
    }
  }
  cat("Note:", x[["note"]], "\n")
  cat("Candidate only: compare with the unmodified AutoSpectral baseline.\n")
  invisible(x)
}
