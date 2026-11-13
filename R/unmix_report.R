#' @title Unmixing Results, Diagnostics and Reports
#'
#' @description
#' Collect the pieces of a spectral unmixing run into one object, compute
#' table-based diagnostics, and write a self-contained report.
#'
#' * [sw_unmix_result()] wraps the output of [sw_unmix_pipeline()] or pieces
#'   built by hand, checks names, and stores a provenance block.
#' * [sw_unmix_diagnose()] computes checks into one tidy table. A verdict
#'   appears only where you supply a limit.
#' * [sw_unmix_report()] writes an HTML page and CSV tables to a new folder.
#'
#' @details
#' The reference baseline is unmodified AutoSpectral 1.8.8 at commit
#' `3c504c4560cf1f0fd5dce298eaabaafaebfae9be` (see `ROADMAP.md`, "Gold
#' standard"). The provenance block records whether the installed copy
#' matches it.
#'
#' Orientation follows AutoSpectral: spectra in rows, detectors in columns.
#' Every matrix needs row and column names. AF and variant matrices are
#' reordered to the detector order of `spectra` by name.
#'
#' The report has tables only. It draws no figures.
#'
#' @name unmix_report
#' @keywords internal
NULL


# ---- Constants ---------------------------------------------------------------

.sw_ref_autospectral_version <- "1.8.8"
.sw_ref_autospectral_commit <- "3c504c4560cf1f0fd5dce298eaabaafaebfae9be"

# Limit name -> metric it judges. The prefix of the name ("max" or "min")
# says whether larger or smaller values are worse.
.sw_unmix_limit_metrics <- c(
  max_cosine = "max_cosine",
  max_sif = "max_sif",
  max_condition_number = "condition_number",
  max_cv = "cv",
  max_pct_negative = "pct_negative",
  max_control_errors = "n_control_errors",
  max_control_warnings = "n_control_warnings",
  min_events = "n_events",
  min_af_spectra = "n_af_spectra",
  min_variants = "n_variants",
  max_share_lost = "share_lost",
  min_median_recovery = "median_recovery",
  max_share_moved = "share_moved",
  max_delta_squeezed = "delta_squeezed",
  max_pct_below_cut = "pct_below_cut",
  max_pct_below_floor = "pct_below_floor",
  max_spread_q = "spread_q",
  max_hyper_share = "hyper_share",
  min_stain_index = "si_mad",
  min_stain_index_q = "si_q",
  max_abs_slope = "abs_slope",
  max_abs_rho = "abs_rho",
  max_heldout_spread = "heldout_spread_q",
  max_heldout_bias = "heldout_bias"
)
.sw_unmix_limit_spec <- data.frame(
  limit = names(.sw_unmix_limit_metrics),
  metric = unname(.sw_unmix_limit_metrics),
  direction = sub("_.*$", "", names(.sw_unmix_limit_metrics)),
  stringsAsFactors = FALSE
)

# Delimiter of the `item` column for spike-in and screen rows:
# "setup|marker|pop".
.sw_unmix_item_sep <- "|"

.sw_unmix_statuses <- c("pass", "warn", "fail", "not_judged", "not_run")

# Files sw_unmix_report() writes; overwrite = TRUE replaces only these.
.sw_unmix_report_files <- c(
  "report.html",
  "checks.csv",
  "spectra.csv",
  "provenance.csv",
  "files.csv",
  "report.json",
  "spike_summary.csv",
  "baseline_check.csv",
  "marker_screen.csv",
  "event_channels.csv",
  "event_pairs.csv",
  "heldout_scores.csv",
  "heldout_repro.csv"
)


# ---- Internal helpers: names and shapes --------------------------------------

# Numeric matrix with unique, non-empty row and column names.
.sw_check_named_matrix <- function(x, arg) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(
      sprintf(
        "`%s` must be a numeric matrix (spectra in rows, detectors in columns).",
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
        "`%s` must have row names (spectra) and column names (detectors).",
        arg
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(rn) > 0L || anyDuplicated(cn) > 0L) {
    stop(
      sprintf("`%s` has duplicated row or column names.", arg),
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` must not contain NA, NaN or Inf.", arg), call. = FALSE)
  }
  invisible(x)
}

# Same detectors as `spectra`, returned in the order of `spectra`.
.sw_align_detectors <- function(m, detectors, arg) {
  if (!setequal(colnames(m), detectors) || ncol(m) != length(detectors)) {
    stop(
      sprintf(
        "`%s` detectors do not match `colnames(spectra)`.",
        arg
      ),
      call. = FALSE
    )
  }
  m[, detectors, drop = FALSE]
}

# AF spectra as a named list of aligned matrices (one per AF set or tissue).
.sw_normalise_af <- function(af, detectors) {
  if (is.null(af)) {
    return(NULL)
  }
  if (is.matrix(af)) {
    af <- list(af = af)
  }
  if (!is.list(af) || length(af) == 0L) {
    stop(
      "`af_spectra` must be a matrix or a named list of matrices.",
      call. = FALSE
    )
  }
  if (
    is.null(names(af)) || any(!nzchar(names(af))) || anyDuplicated(names(af))
  ) {
    stop("`af_spectra` list must have unique, non-empty names.", call. = FALSE)
  }
  out <- lapply(names(af), \(nm) {
    arg <- sprintf("af_spectra[[\"%s\"]]", nm)
    .sw_check_named_matrix(af[[nm]], arg)
    .sw_align_detectors(af[[nm]], detectors, arg)
  })
  names(out) <- names(af)
  out
}

# Variants as a named list of aligned matrices. Accepts an sw_variants object,
# a list with a `variants` element, or a bare named list. A bare named vector
# (AutoSpectral's single-spectrum early return) becomes a one-row matrix.
.sw_normalise_variants <- function(v, fluors, detectors) {
  if (is.null(v)) {
    return(NULL)
  }
  if (inherits(v, "sw_variants") || (is.list(v) && "variants" %in% names(v))) {
    v <- v[["variants"]]
  }
  if (
    !is.list(v) ||
      length(v) == 0L ||
      is.null(names(v)) ||
      any(!nzchar(names(v))) ||
      anyDuplicated(names(v))
  ) {
    stop(
      "`spectra_variants` must be a named list of variant matrices, or an sw_variants object.",
      call. = FALSE
    )
  }
  unknown <- setdiff(names(v), fluors)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "`spectra_variants` names not in `rownames(spectra)`: %s.",
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  out <- lapply(names(v), \(nm) {
    m <- v[[nm]]
    arg <- sprintf("spectra_variants[[\"%s\"]]", nm)
    if (is.numeric(m) && is.null(dim(m)) && !is.null(names(m))) {
      m <- matrix(m, nrow = 1L, dimnames = list(paste0(nm, "_1"), names(m)))
    }
    if (!is.matrix(m) || !is.numeric(m) || is.null(colnames(m))) {
      stop(
        sprintf("`%s` must be a numeric matrix with detector names.", arg),
        call. = FALSE
      )
    }
    .sw_align_detectors(m, detectors, arg)
  })
  names(out) <- names(v)
  out
}

# Cosine similarity between the rows of a matrix; NA for zero-norm rows.
.sw_row_cosine <- function(m) {
  norms <- sqrt(rowSums(m^2))
  norms[norms == 0] <- NA_real_
  unit <- m / norms
  unit %*% t(unit)
}


# ---- Internal helpers: provenance --------------------------------------------

# Package version and RemoteSha, NA when the package is absent.
.sw_pkg_info <- function(pkg) {
  desc <- tryCatch(
    suppressWarnings(utils::packageDescription(pkg)),
    error = \(e) NA
  )
  if (!is.list(desc)) {
    return(tibble::tibble(
      package = pkg,
      version = NA_character_,
      remote_sha = NA_character_
    ))
  }
  sha <- desc[["RemoteSha"]]
  tibble::tibble(
    package = pkg,
    version = as.character(desc[["Version"]]),
    remote_sha = if (is.null(sha)) NA_character_ else as.character(sha)
  )
}

# Inventory of files: path, role, existence, size and MD5.
.sw_file_inventory <- function(paths, role) {
  if (length(paths) == 0L) {
    return(tibble::tibble(
      role = character(),
      path = character(),
      exists = logical(),
      size_bytes = numeric(),
      md5 = character()
    ))
  }
  exists <- file.exists(paths) & !dir.exists(paths)
  size <- rep(NA_real_, length(paths))
  md5 <- rep(NA_character_, length(paths))
  size[exists] <- file.size(paths[exists])
  md5[exists] <- unname(tools::md5sum(paths[exists]))
  tibble::tibble(
    role = role,
    path = paths,
    exists = exists,
    size_bytes = size,
    md5 = md5
  )
}

# Does the installed AutoSpectral match the reference pin?
.sw_baseline_check <- function(as_info) {
  ref <- sprintf(
    "AutoSpectral %s @ %s",
    .sw_ref_autospectral_version,
    .sw_ref_autospectral_commit
  )
  ver <- as_info$version
  sha <- as_info$remote_sha
  msg <- if (is.na(ver)) {
    "AutoSpectral is not installed, so the run cannot be tied to the reference."
  } else if (ver != .sw_ref_autospectral_version) {
    sprintf(
      "Installed AutoSpectral %s differs from the reference %s.",
      ver,
      .sw_ref_autospectral_version
    )
  } else if (is.na(sha)) {
    "AutoSpectral 1.8.8 is installed but records no RemoteSha, so the commit cannot be confirmed."
  } else if (sha != .sw_ref_autospectral_commit) {
    sprintf(
      "Installed AutoSpectral commit %s differs from the reference commit %s.",
      sha,
      .sw_ref_autospectral_commit
    )
  } else {
    "Installed AutoSpectral matches the reference version and commit."
  }
  list(
    reference = ref,
    version = .sw_ref_autospectral_version,
    commit = .sw_ref_autospectral_commit,
    matches_reference = !is.na(ver) &&
      ver == .sw_ref_autospectral_version &&
      !is.na(sha) &&
      sha == .sw_ref_autospectral_commit,
    message = msg
  )
}

# Provenance block stored in every sw_unmix_result.
.sw_unmix_provenance <- function(
  inputs = NULL,
  outputs = NULL,
  args = NULL,
  seed = NULL
) {
  pkgs <- do.call(
    dplyr::bind_rows,
    lapply(
      c("SpectraWeaveR", "AutoSpectral", "AutoSpectralRcpp", "flowCore"),
      .sw_pkg_info
    )
  )
  as_info <- pkgs[pkgs$package == "AutoSpectral", ]
  rcpp_installed <- !is.na(pkgs$version[pkgs$package == "AutoSpectralRcpp"])
  list(
    r_version = R.version.string,
    platform = R.version$platform,
    datetime = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    packages = pkgs,
    autospectral_commit = as_info$remote_sha,
    autospectralrcpp_installed = rcpp_installed,
    backend = if (rcpp_installed) {
      "AutoSpectralRcpp installed (compiled backend available)"
    } else {
      "AutoSpectralRcpp not installed (R backend)"
    },
    args = args,
    seed = seed,
    files = dplyr::bind_rows(
      .sw_file_inventory(unique(inputs), "input"),
      .sw_file_inventory(unique(outputs), "output")
    ),
    baseline = .sw_baseline_check(as_info)
  )
}


# ---- Constructor -------------------------------------------------------------

#' Collect an unmixing run into one object
#'
#' Wraps the list returned by [sw_unmix_pipeline()], or pieces you assembled
#' by hand after running the steps separately. It checks names and shapes and
#' attaches a provenance block.
#'
#' @details
#' Give either `x` or the separate pieces, not both for the same piece.
#'
#' Name checks follow AutoSpectral's orientation (spectra in rows, detectors
#' in columns). `spectra` needs unique row and column names and finite values.
#' Each AF and variant matrix must have the same detector set as `spectra`; it
#' is reordered to the column order of `spectra` by name. Variant list names
#' must be row names of `spectra`.
#'
#' The provenance block holds the R version, platform, time, the versions of
#' SpectraWeaveR, AutoSpectral, AutoSpectralRcpp and flowCore (`NA` if
#' absent), the AutoSpectral `RemoteSha` (`NA` if not recorded), whether
#' AutoSpectralRcpp is installed, `args`, `seed`, and an inventory (path,
#' size, MD5) of input and output files. Installed is not the same as used:
#' AutoSpectral picks its backend at run time.
#'
#' `baseline$matches_reference` is `TRUE` only when the installed AutoSpectral
#' is version 1.8.8 and its `RemoteSha` is
#' `3c504c4560cf1f0fd5dce298eaabaafaebfae9be`. Otherwise it is `FALSE` and
#' `baseline$message` says why. The check reads the package installed now,
#' which may differ from the one that ran.
#'
#' @param x Optional list from [sw_unmix_pipeline()] with elements `setup`,
#'   `flow_control`, `spectra`, `af_spectra`, `spectra_variants` and
#'   `unmixed` (from [sw_unmix_run()]).
#' @param spectra Numeric matrix, fluorophores x detectors.
#' @param af_spectra Optional AF matrix (AF spectra x detectors) or a named
#'   list of such matrices (one per tissue).
#' @param spectra_variants Optional `sw_variants` object from
#'   [sw_unmix_extract_variants()], or a named list of variant matrices (row 1
#'   is the reference spectrum).
#' @param unmixed_files Optional character vector of unmixed FCS paths.
#' @param setup Optional `sw_setup` object from [sw_unmix_setup()].
#' @param flow_control Optional AutoSpectral `flow.control` list.
#' @param method Optional unmixing method name.
#' @param args Optional named list of the run's arguments, stored as given.
#' @param seed Optional seed used for the run.
#' @param inputs Optional character vector of input file paths (controls,
#'   unstained, samples) for the file inventory. The control file in `setup`
#'   is added automatically.
#' @param run_id Optional run identifier. Without one, a time stamp plus the
#'   first 8 characters of the MD5 of `spectra` is used.
#'
#' @return An object of class `sw_unmix_result`: a list with `run_id`,
#'   `spectra`, `af_spectra` (named list or `NULL`), `variants` (named list
#'   or `NULL`), `spectra_variants` (as given), `unmixed_files`, `method`,
#'   `setup`, `flow_control` and `provenance`.
#'
#' @examples
#' spectra <- matrix(
#'   c(1, 0.5, 0.1, 0, 0.1, 1, 0.4, 0.1, 0, 0.2, 1, 0.6),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("BV421", "PE", "APC"), paste0("D", 1:4))
#' )
#' res <- sw_unmix_result(spectra = spectra, seed = 1L, run_id = "demo")
#' res
#'
#' @export
sw_unmix_result <- function(
  x = NULL,
  spectra = NULL,
  af_spectra = NULL,
  spectra_variants = NULL,
  unmixed_files = NULL,
  setup = NULL,
  flow_control = NULL,
  method = NULL,
  args = NULL,
  seed = NULL,
  inputs = NULL,
  run_id = NULL
) {
  if (!is.null(x)) {
    if (inherits(x, "sw_unmix_result")) {
      stop("`x` is already an sw_unmix_result.", call. = FALSE)
    }
    if (!is.list(x) || is.null(x[["spectra"]])) {
      stop(
        "`x` must be the list from sw_unmix_pipeline() with a `spectra` element.",
        call. = FALSE
      )
    }
    from_x <- list(
      spectra = x[["spectra"]],
      af_spectra = x[["af_spectra"]],
      spectra_variants = x[["spectra_variants"]],
      unmixed_files = x[["unmixed"]][["output_files"]],
      setup = x[["setup"]],
      flow_control = x[["flow_control"]],
      method = x[["unmixed"]][["method"]]
    )
    given <- list(
      spectra = spectra,
      af_spectra = af_spectra,
      spectra_variants = spectra_variants,
      unmixed_files = unmixed_files,
      setup = setup,
      flow_control = flow_control,
      method = method
    )
    both <- names(given)[
      !vapply(given, is.null, logical(1)) &
        !vapply(from_x, is.null, logical(1))
    ]
    if (length(both) > 0L) {
      stop(
        sprintf(
          "Given both in `x` and as an argument: %s.",
          paste(both, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    # Fill each piece from `x` only where the argument was not given.
    for (nm in names(from_x)) {
      if (is.null(given[[nm]])) given[[nm]] <- from_x[[nm]]
    }
    spectra <- given$spectra
    af_spectra <- given$af_spectra
    spectra_variants <- given$spectra_variants
    unmixed_files <- given$unmixed_files
    setup <- given$setup
    flow_control <- given$flow_control
    method <- given$method
  }

  if (is.null(spectra)) {
    stop("`spectra` is required.", call. = FALSE)
  }
  .sw_check_named_matrix(spectra, "spectra")
  detectors <- colnames(spectra)

  af <- .sw_normalise_af(af_spectra, detectors)
  variants <- .sw_normalise_variants(
    spectra_variants,
    rownames(spectra),
    detectors
  )

  if (
    !is.null(unmixed_files) &&
      (!is.character(unmixed_files) || anyNA(unmixed_files))
  ) {
    stop(
      "`unmixed_files` must be a character vector of paths.",
      call. = FALSE
    )
  }
  if (!is.null(inputs) && (!is.character(inputs) || anyNA(inputs))) {
    stop("`inputs` must be a character vector of paths.", call. = FALSE)
  }
  if (!is.null(setup) && !inherits(setup, "sw_setup")) {
    stop("`setup` must be an sw_setup object.", call. = FALSE)
  }
  if (
    !is.null(method) &&
      (!is.character(method) || length(method) != 1L || is.na(method))
  ) {
    stop("`method` must be a single string.", call. = FALSE)
  }
  if (
    !is.null(args) &&
      (!is.list(args) || is.null(names(args)) || any(!nzchar(names(args))))
  ) {
    stop("`args` must be a named list.", call. = FALSE)
  }
  if (
    !is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed))
  ) {
    stop("`seed` must be a single number.", call. = FALSE)
  }
  if (
    !is.null(run_id) &&
      (!is.character(run_id) ||
        length(run_id) != 1L ||
        is.na(run_id) ||
        !nzchar(run_id))
  ) {
    stop("`run_id` must be a single non-empty string.", call. = FALSE)
  }

  if (is.null(run_id)) {
    # Hash the spectra so two runs started in the same second differ.
    tmp <- tempfile()
    saveRDS(spectra, tmp)
    run_id <- paste0(
      "unmix-",
      format(Sys.time(), "%Y%m%d-%H%M%S"),
      "-",
      substr(unname(tools::md5sum(tmp)), 1L, 8L)
    )
    unlink(tmp)
  }

  control_file <- if (!is.null(setup)) setup[["control_file"]] else NULL

  structure(
    list(
      run_id = run_id,
      spectra = spectra,
      af_spectra = af,
      variants = variants,
      spectra_variants = spectra_variants,
      unmixed_files = if (is.null(unmixed_files)) {
        character()
      } else {
        unmixed_files
      },
      method = method,
      setup = setup,
      flow_control = flow_control,
      provenance = .sw_unmix_provenance(
        inputs = c(control_file, inputs),
        outputs = unmixed_files,
        args = args,
        seed = seed
      )
    ),
    class = "sw_unmix_result"
  )
}

#' @export
print.sw_unmix_result <- function(x, ...) {
  cat("<sw_unmix_result>", x$run_id, "\n")
  cat(
    "  spectra:",
    nrow(x$spectra),
    "fluorophores x",
    ncol(x$spectra),
    "detectors\n"
  )
  cat(
    "  AF sets:",
    if (is.null(x$af_spectra)) 0L else length(x$af_spectra),
    "| variant sets:",
    if (is.null(x$variants)) 0L else length(x$variants),
    "| unmixed files:",
    length(x$unmixed_files),
    "\n"
  )
  cat("  reference:", x$provenance$baseline$message, "\n")
  invisible(x)
}


# ---- Diagnostics: internal helpers -------------------------------------------

# One or more check rows; status is set later from the limits.
.sw_check_rows <- function(
  check_id,
  scope,
  item,
  metric,
  value,
  reason = "",
  run = TRUE
) {
  tibble::tibble(
    check_id = check_id,
    scope = scope,
    item = item,
    metric = metric,
    value = as.numeric(value),
    limit = NA_character_,
    status = if (run) "not_judged" else "not_run",
    reason = reason
  )
}

# Validate `limits` into a data frame of limit, direction, warn, fail.
.sw_parse_limits <- function(limits) {
  if (!is.list(limits)) {
    stop("`limits` must be a named list.", call. = FALSE)
  }
  if (length(limits) == 0L) {
    return(data.frame(
      limit = character(),
      direction = character(),
      warn = numeric(),
      fail = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  if (is.null(names(limits)) || any(!nzchar(names(limits)))) {
    stop("Every element of `limits` must be named.", call. = FALSE)
  }
  unknown <- setdiff(names(limits), .sw_unmix_limit_spec$limit)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown limit name(s): %s. Allowed: %s.",
        paste(unknown, collapse = ", "),
        paste(.sw_unmix_limit_spec$limit, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(names(limits))) {
    stop("`limits` has duplicated names.", call. = FALSE)
  }
  rows <- lapply(names(limits), \(nm) {
    v <- limits[[nm]]
    if (!is.numeric(v) || length(v) == 0L || length(v) > 2L || anyNA(v)) {
      stop(
        sprintf("`limits$%s` must be one or two non-missing numbers.", nm),
        call. = FALSE
      )
    }
    # An unnamed single number is the fail limit.
    if (is.null(names(v)) && length(v) == 1L) {
      names(v) <- "fail"
    }
    if (
      is.null(names(v)) ||
        !all(names(v) %in% c("warn", "fail")) ||
        anyDuplicated(names(v))
    ) {
      stop(
        sprintf(
          "`limits$%s` must be named `warn` and/or `fail`, e.g. c(warn = 1, fail = 2).",
          nm
        ),
        call. = FALSE
      )
    }
    direction <- .sw_unmix_limit_spec$direction[
      .sw_unmix_limit_spec$limit == nm
    ]
    warn <- if ("warn" %in% names(v)) v[["warn"]] else NA_real_
    fail <- if ("fail" %in% names(v)) v[["fail"]] else NA_real_
    bad_order <- !is.na(warn) &&
      !is.na(fail) &&
      ((direction == "max" && warn > fail) ||
        (direction == "min" && warn < fail))
    if (bad_order) {
      stop(
        sprintf(
          "`limits$%s`: warn must be %s fail.",
          nm,
          if (direction == "max") "<=" else ">="
        ),
        call. = FALSE
      )
    }
    data.frame(
      limit = nm,
      direction = direction,
      warn = warn,
      fail = fail,
      stringsAsFactors = FALSE
    )
  })
  do.call(dplyr::bind_rows, rows)
}

# Set limit text and status on rows whose metric has a supplied limit.
.sw_apply_limits <- function(checks, lims) {
  for (i in seq_len(nrow(lims))) {
    metric <- .sw_unmix_limit_spec$metric[
      .sw_unmix_limit_spec$limit == lims$limit[i]
    ]
    op <- if (lims$direction[i] == "max") ">" else "<"
    breach <- if (lims$direction[i] == "max") `>` else `<`
    parts <- c(
      if (!is.na(lims$warn[i])) paste("warn if", op, lims$warn[i]),
      if (!is.na(lims$fail[i])) paste("fail if", op, lims$fail[i])
    )
    rows <- checks$metric == metric & checks$status != "not_run"
    checks$limit[rows] <- paste(parts, collapse = "; ")
    for (r in which(rows)) {
      v <- checks$value[r]
      if (is.na(v)) {
        checks$reason[r] <- trimws(paste(
          checks$reason[r],
          "Value is NA, so the limit was not applied."
        ))
        next
      }
      checks$status[r] <- if (!is.na(lims$fail[i]) && breach(v, lims$fail[i])) {
        "fail"
      } else if (!is.na(lims$warn[i]) && breach(v, lims$warn[i])) {
        "warn"
      } else {
        "pass"
      }
    }
  }
  checks
}

# Mean, sd, cv and percent negative per channel of a numeric matrix (no flag).
.sw_channel_stats <- function(mat, channels) {
  do.call(
    dplyr::bind_rows,
    lapply(channels, \(ch) {
      vals <- mat[, ch]
      ch_mean <- mean(vals, na.rm = TRUE)
      ch_sd <- stats::sd(vals, na.rm = TRUE)
      tibble::tibble(
        channel = ch,
        mean = ch_mean,
        sd = ch_sd,
        cv = if (abs(ch_mean) > .Machine$double.eps) {
          ch_sd / abs(ch_mean)
        } else {
          NA_real_
        },
        pct_negative = 100 * sum(vals < 0, na.rm = TRUE) / length(vals)
      )
    })
  )
}

# Control-file checks through sw_unmix_check_controls().
.sw_diag_controls <- function(setup, dots) {
  have <- !is.null(setup) &&
    !is.null(setup[["control_dir"]]) &&
    !is.null(setup[["control_file"]]) &&
    is.list(setup[["asp"]])
  if (!have) {
    return(.sw_check_rows(
      "controls",
      "run",
      "control_file",
      c("n_control_errors", "n_control_warnings"),
      NA,
      "No setup with control_dir, control_file and asp was given.",
      run = FALSE
    ))
  }
  if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
    return(.sw_check_rows(
      "controls",
      "run",
      "control_file",
      c("n_control_errors", "n_control_warnings"),
      NA,
      "AutoSpectral is not installed.",
      run = FALSE
    ))
  }
  issues <- tryCatch(
    do.call(
      sw_unmix_check_controls,
      c(
        list(
          control_dir = setup[["control_dir"]],
          control_def_file = setup[["control_file"]],
          asp = setup[["asp"]]
        ),
        dots
      )
    ),
    error = \(e) e
  )
  if (inherits(issues, "error")) {
    return(.sw_check_rows(
      "controls",
      "run",
      "control_file",
      c("n_control_errors", "n_control_warnings"),
      NA,
      paste("Control check stopped:", conditionMessage(issues)),
      run = FALSE
    ))
  }
  sev <- if (nrow(issues) > 0L) issues$severity else character()
  run_rows <- .sw_check_rows(
    "controls",
    "run",
    "control_file",
    c("n_control_errors", "n_control_warnings"),
    c(sum(sev == "error"), sum(sev == "warning")),
    "Counts of issues from AutoSpectral::check.control.file()."
  )
  if (nrow(issues) == 0L) {
    return(run_rows)
  }
  # One row per issue; AutoSpectral's severity is reported, not used as a verdict.
  issue_rows <- .sw_check_rows(
    "controls",
    "file",
    ifelse(is.na(issues$filename), "control_file", issues$filename),
    paste0("control_issue:", issues$rule),
    NA,
    sprintf("AutoSpectral severity %s: %s", issues$severity, issues$message)
  )
  dplyr::bind_rows(run_rows, issue_rows)
}

# Per-fluorophore maximum cosine similarity to any other spectrum.
.sw_diag_similarity <- function(spectra) {
  if (nrow(spectra) < 2L) {
    return(.sw_check_rows(
      "spectra_similarity",
      "fluorophore",
      rownames(spectra),
      "max_cosine",
      NA,
      "Needs at least 2 spectra.",
      run = FALSE
    ))
  }
  cos <- .sw_row_cosine(spectra)
  diag(cos) <- NA_real_
  best <- apply(cos, 1L, \(r) if (all(is.na(r))) NA_integer_ else which.max(r))
  value <- vapply(
    seq_len(nrow(cos)),
    \(i) if (is.na(best[i])) NA_real_ else cos[i, best[i]],
    numeric(1)
  )
  .sw_check_rows(
    "spectra_similarity",
    "fluorophore",
    rownames(spectra),
    "max_cosine",
    value,
    ifelse(
      is.na(best),
      "Spectrum has zero norm.",
      ifelse(
        value <= 0,
        "No overlap with any other spectrum.",
        paste("Closest spectrum:", rownames(spectra)[best])
      )
    )
  )
}

# Spreading inflation factor (SIF) of Mage et al. 2026: the diagonal of the hotspot
# matrix sqrt(abs(pinv(S))), S being the cosine-similarity matrix. This is the
# matrix of AutoSpectral::calculate.hotspot.matrix() in 1.8.8. The pseudo-inverse is
# built from base svd() with the tolerance of MASS::ginv(). The reason names the
# largest off-diagonal entry in the row, which is the partner that drives the SIF.
.sw_diag_sif <- function(spectra) {
  skip <- function(why) {
    .sw_check_rows(
      "spectra_sif",
      "fluorophore",
      rownames(spectra),
      "max_sif",
      NA,
      why,
      run = FALSE
    )
  }
  if (nrow(spectra) < 2L) {
    return(skip("Needs at least 2 spectra."))
  }
  cos <- .sw_row_cosine(spectra)
  if (anyNA(cos)) {
    return(skip("A spectrum has zero norm."))
  }
  sv <- svd(cos)
  keep <- sv$d > max(sqrt(.Machine$double.eps) * sv$d[1L], 0)
  inv <- sv$v[, keep, drop = FALSE] %*%
    (t(sv$u[, keep, drop = FALSE]) / sv$d[keep])
  hot <- sqrt(abs(inv))
  off <- hot
  diag(off) <- NA_real_
  partner <- apply(off, 1L, which.max)
  .sw_check_rows(
    "spectra_sif",
    "fluorophore",
    rownames(spectra),
    "max_sif",
    diag(hot),
    paste0(
      "Diagonal of the hotspot matrix. Largest off-diagonal entry: ",
      rownames(spectra)[partner],
      " (",
      signif(off[cbind(seq_len(nrow(off)), partner)], 3),
      ")."
    )
  )
}

# Condition number of the spectra matrix by base svd().
.sw_diag_condition <- function(spectra) {
  if (nrow(spectra) < 2L) {
    return(.sw_check_rows(
      "condition_number",
      "run",
      "spectra",
      "condition_number",
      NA,
      "Needs at least 2 spectra.",
      run = FALSE
    ))
  }
  d <- svd(spectra)$d
  .sw_check_rows(
    "condition_number",
    "run",
    "spectra",
    "condition_number",
    max(d) / min(d),
    "max/min singular value from base svd(); not rounded."
  )
}

# AF library: count and minimum pairwise cosine per AF set.
.sw_diag_af <- function(af) {
  if (is.null(af)) {
    return(.sw_check_rows(
      "af_library",
      "run",
      "af",
      c("n_af_spectra", "min_af_cosine"),
      NA,
      "No AF spectra were given.",
      run = FALSE
    ))
  }
  do.call(
    dplyr::bind_rows,
    lapply(names(af), \(nm) {
      m <- af[[nm]]
      n_rows <- .sw_check_rows(
        "af_library",
        "run",
        nm,
        "n_af_spectra",
        nrow(m)
      )
      if (nrow(m) < 2L) {
        cos_row <- .sw_check_rows(
          "af_library",
          "run",
          nm,
          "min_af_cosine",
          NA,
          "Needs at least 2 AF spectra.",
          run = FALSE
        )
      } else {
        cos <- .sw_row_cosine(m)
        cos_row <- .sw_check_rows(
          "af_library",
          "run",
          nm,
          "min_af_cosine",
          min(cos[upper.tri(cos)], na.rm = FALSE),
          "Minimum cosine similarity between AF spectra; no limit option."
        )
      }
      dplyr::bind_rows(n_rows, cos_row)
    })
  )
}

# Variants: number per fluorophore (rows minus the reference row).
.sw_diag_variants <- function(variants, fluors) {
  if (is.null(variants)) {
    return(.sw_check_rows(
      "variants",
      "fluorophore",
      fluors,
      "n_variants",
      NA,
      "No spectral variants were given.",
      run = FALSE
    ))
  }
  have <- fluors %in% names(variants)
  value <- rep(NA_real_, length(fluors))
  value[have] <- vapply(
    variants[fluors[have]],
    \(m) nrow(m) - 1,
    numeric(1)
  )
  rows <- .sw_check_rows(
    "variants",
    "fluorophore",
    fluors,
    "n_variants",
    value,
    ifelse(have, "Rows after the reference spectrum (row 1).", "")
  )
  rows$status[!have] <- "not_run"
  rows$reason[!have] <- "No variant matrix for this fluorophore."
  rows
}

# Unmixed channel for each fluorophore: exact name, else name + "-A".
.sw_match_channels <- function(fluors, cols) {
  vapply(
    fluors,
    \(f) {
      if (f %in% cols) {
        f
      } else if (paste0(f, "-A") %in% cols) {
        paste0(f, "-A")
      } else {
        NA_character_
      }
    },
    character(1)
  )
}

# Per-file event count and per-channel statistics.
.sw_diag_one_file <- function(name, mat, stats_fn, fluors) {
  file_row <- .sw_check_rows(
    "unmixed_file",
    "file",
    name,
    "n_events",
    nrow(mat)
  )
  chan <- .sw_match_channels(fluors, colnames(mat))
  found <- !is.na(chan)
  miss_rows <- if (any(!found)) {
    .sw_check_rows(
      "unmixed_channel",
      "channel",
      paste0(name, " | ", fluors[!found]),
      "channel_stats",
      NA,
      sprintf("No channel named %s or %s-A.", fluors[!found], fluors[!found]),
      run = FALSE
    )
  } else {
    NULL
  }
  if (!any(found)) {
    return(dplyr::bind_rows(file_row, miss_rows))
  }
  st <- stats_fn(chan[found])
  item <- paste0(name, " | ", chan[found])
  # Long form: one row per channel and statistic.
  stat_rows <- do.call(
    dplyr::bind_rows,
    lapply(c("mean", "sd", "cv", "pct_negative"), \(s) {
      .sw_check_rows(
        "unmixed_channel",
        "channel",
        item,
        s,
        st[[s]][match(chan[found], st$channel)]
      )
    })
  )
  dplyr::bind_rows(file_row, stat_rows, miss_rows)
}

.sw_diag_files <- function(unmixed_files, unmixed_data, fluors) {
  if (!is.null(unmixed_data)) {
    return(do.call(
      dplyr::bind_rows,
      lapply(names(unmixed_data), \(nm) {
        m <- unmixed_data[[nm]]
        .sw_diag_one_file(nm, m, \(ch) .sw_channel_stats(m, ch), fluors)
      })
    ))
  }
  if (length(unmixed_files) == 0L) {
    return(.sw_check_rows(
      "unmixed_file",
      "file",
      "unmixed",
      "n_events",
      NA,
      "No unmixed files or data were given.",
      run = FALSE
    ))
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    return(.sw_check_rows(
      "unmixed_file",
      "file",
      unmixed_files,
      "n_events",
      NA,
      "flowCore is not installed; pass `unmixed_data` instead.",
      run = FALSE
    ))
  }
  do.call(
    dplyr::bind_rows,
    lapply(unmixed_files, \(f) {
      if (!file.exists(f)) {
        return(.sw_check_rows(
          "unmixed_file",
          "file",
          f,
          "n_events",
          NA,
          "File does not exist.",
          run = FALSE
        ))
      }
      ff <- flowCore::read.FCS(
        f,
        transformation = FALSE,
        truncate_max_range = FALSE
      )
      mat <- flowCore::exprs(ff)
      .sw_diag_one_file(basename(f), mat, \(ch) .sw_channel_stats(mat, ch), fluors)
    })
  )
}

# Donor range text for the reason column.
.sw_donor_reason <- function(n, lo, hi) {
  sprintf("%d donor(s), range %s to %s", n, signif(lo, 3), signif(hi, 3))
}

# Spike-in recovery rows from an sw_spike_score. One row per setup x set x
# pop for share_lost and median_recovery (spiked levels pooled per donor),
# and one per setup x set x group for share_moved (worst channel per donor).
# The value is the donor median.
.sw_diag_spike <- function(spike) {
  if (is.null(spike)) {
    return(NULL)
  }
  if (!inherits(spike, "sw_spike_score")) {
    stop("`spike` must be an sw_spike_score.", call. = FALSE)
  }
  sep <- .sw_unmix_item_sep
  coupled <- spike$file_status |>
    dplyr::summarise(
      n_coupled = sum(.data$file_status == "coupled"),
      .by = c("setup", "set")
    )
  per_donor <- spike$cells[spike$cells$level > 0, ] |>
    dplyr::summarise(
      share_lost = mean(.data$recovery < spike$loss_cut),
      median_recovery = stats::median(.data$recovery),
      levels = paste(sort(unique(.data$level)), collapse = ", "),
      .by = c("setup", "set", "pop", "sample")
    )
  recovery <- lapply(c("share_lost", "median_recovery"), \(metric) {
    per_donor |>
      dplyr::summarise(
        value = stats::median(.data[[metric]]),
        lo = min(.data[[metric]]),
        hi = max(.data[[metric]]),
        n = dplyr::n(),
        levels = dplyr::first(.data$levels),
        .by = c("setup", "set", "pop")
      ) |>
      dplyr::left_join(
        coupled,
        by = c("setup", "set"),
        relationship = "many-to-one"
      ) |>
      dplyr::mutate(metric = .env$metric)
  }) |>
    dplyr::bind_rows()
  rows_recovery <- .sw_check_rows(
    "spike_recovery",
    "setup",
    paste(recovery$setup, recovery$set, recovery$pop, sep = sep),
    recovery$metric,
    recovery$value,
    reason = sprintf(
      "%s; levels %s pooled; loss cut %s; %d coupled file(s).",
      .sw_donor_reason(recovery$n, recovery$lo, recovery$hi),
      recovery$levels,
      format(spike$loss_cut),
      dplyr::coalesce(recovery$n_coupled, 0L)
    )
  )
  moved <- spike$baseline_check |>
    dplyr::summarise(
      share_moved = max(.data$share_moved, na.rm = TRUE),
      .by = c("setup", "set", "group", "sample")
    ) |>
    dplyr::summarise(
      value = stats::median(.data$share_moved),
      lo = min(.data$share_moved),
      hi = max(.data$share_moved),
      n = dplyr::n(),
      .by = c("setup", "set", "group")
    )
  rows_moved <- .sw_check_rows(
    "spike_recovery",
    "setup",
    paste(moved$setup, moved$set, moved$group, sep = sep),
    "share_moved",
    moved$value,
    reason = sprintf(
      "%s; worst channel per donor; tolerance %s arcsinh.",
      .sw_donor_reason(moved$n, moved$lo, moved$hi),
      format(spike$baseline_tol)
    )
  )
  dplyr::bind_rows(rows_recovery, rows_moved)
}

# Marker screen rows from an sw_marker_screen: one row per setup x channel x
# pop; the value is the donor-median rise in squeezed share.
.sw_diag_screen <- function(screen) {
  if (is.null(screen)) {
    return(NULL)
  }
  if (!inherits(screen, "sw_marker_screen")) {
    stop("`screen` must be an sw_marker_screen.", call. = FALSE)
  }
  ranges <- screen$per_donor |>
    dplyr::filter(!is.na(.data$d_squeezed)) |>
    dplyr::summarise(
      lo = min(.data$d_squeezed),
      hi = max(.data$d_squeezed),
      .by = c("setup", "pop", "channel")
    )
  f <- screen$flags |>
    dplyr::left_join(
      ranges,
      by = c("setup", "pop", "channel"),
      relationship = "one-to-one"
    )
  .sw_check_rows(
    "marker_screen",
    "setup",
    paste(f$setup, f$channel, f$pop, sep = .sw_unmix_item_sep),
    "delta_squeezed",
    f$median_d_squeezed,
    reason = sprintf(
      "%s; screen flag %s (reference %s).",
      .sw_donor_reason(f$n_donors, f$lo, f$hi),
      ifelse(is.na(f$flag), "NA", f$flag),
      screen$reference
    )
  )
}


# ---- Diagnose ----------------------------------------------------------------

#' SIF limits from the AutoSpectral workflow article
#'
#' Returns the `limits` entry for [sw_unmix_diagnose()] that the user chose on
#' 2026-11-09. Values from 4 to 6 are "worth reviewing" (`warn`), and values
#' above 6 are "most likely to cause unmixing problems" (`fail`). Source:
#' AutoSpectral article `01_Automated_Spectral_Workflow`, which cites Mage et
#' al. (2026) and gives no formula and no source for 4 and 6. Mage et al. define
#' the spreading inflation factor (SIF) as the diagonal of the hotspot matrix
#' and give 3 and 8 as rough guidelines (below 3 minimal, 3 to 8 moderate,
#' above 8 severe). This function applies the AutoSpectral cutoffs to the SIF,
#' not the Mage ones. Mage et al. also say the SIF can overstate measured
#' spreading about two-fold.
#'
#' The function is a separate call so that the limit is never a hidden
#' default of [sw_unmix_diagnose()]. Pass it explicitly.
#'
#' @return A named list for the `limits` argument of [sw_unmix_diagnose()].
#'
#' @examples
#' sw_unmix_limits_autospectral()
#'
#' @export
sw_unmix_limits_autospectral <- function() {
  list(max_sif = c(warn = 4, fail = 6))
}

#' Diagnose an unmixing run
#'
#' Computes table-based checks on an [sw_unmix_result()] and returns them in
#' one tidy table. Each check runs only when its inputs exist; otherwise its
#' status is `"not_run"` with the reason.
#'
#' @details
#' **No cutoff is built in.** A status of `pass`, `warn` or `fail` appears
#' only for metrics whose limit you supply in `limits`. Without a limit the
#' status is `"not_judged"` and the value is still reported. The numbers are
#' not calibrated against ground truth; a limit you pass is your choice.
#'
#' Checks:
#'
#' * `controls`: [sw_unmix_check_controls()] on the setup's control folder
#'   and file. Rows count AutoSpectral errors and warnings, plus one row per
#'   issue with AutoSpectral's severity in `reason` (not used as a verdict).
#'   Needs `result$setup` and AutoSpectral. `...` goes to that function.
#' * `spectra_similarity`: per fluorophore, the largest cosine similarity to
#'   any other row of `spectra` (base R).
#' * `spectra_sif`: per fluorophore, the spreading inflation factor of Mage
#'   et al. (2026): the diagonal of the hotspot matrix `sqrt(abs(pinv(S)))`,
#'   where `S` is the cosine-similarity matrix. This is the matrix of
#'   `AutoSpectral::calculate.hotspot.matrix()` in 1.8.8, computed in base R.
#'   The reason names the largest off-diagonal entry of the row, the partner
#'   that drives the value. Cutoffs are in [sw_unmix_limits_autospectral()].
#' * `condition_number`: `max(d) / min(d)` of `svd(spectra)$d` in base R. This
#'   is the formula of `AutoSpectral::calculate.condition.number()` in 1.8.8,
#'   which also rounds to 2 decimals; this value is not rounded. Base R is
#'   used so the check runs without AutoSpectral.
#' * `af_library`: per AF set, the number of AF spectra and the minimum
#'   pairwise cosine between them. The minimum cosine has no limit option,
#'   because neither direction is clearly worse.
#' * `variants`: per fluorophore, the number of variants (rows after the
#'   reference row 1). A fluorophore without a matrix is `not_run`.
#' * `unmixed_file` and `unmixed_channel`: per file, the event count; per
#'   fluorophore channel, `mean`, `sd`, `cv` (`sd / |mean|`) and
#'   `pct_negative`. A channel is the fluorophore name, or the name plus
#'   `"-A"` (AutoSpectral's unmixed FCS naming). With `unmixed_data` the
#'   statistics come from those matrices. Otherwise files are read with
#'   flowCore (`truncate_max_range = FALSE`) and summarised by
#'   the formulas of the removed `sw_unmix_quality()`, without a flag.
#' * `spike_recovery` (only with `spike`): from an [sw_unmix_spike_score()]
#'   result. `item` is `"setup|marker|pop"`, where the marker slot holds the
#'   spike set (the marker for single spikes, `cospike_<dye2>_<level2>` for
#'   co-spikes). `share_lost` and `median_recovery` pool each donor's spiked
#'   cells (levels above 0); `share_moved` is the worst channel per donor of
#'   the independence check, with `"level0"` or `"untouched"` in the pop slot.
#' * `marker_screen` (only with `screen`): from [sw_unmix_marker_screen()].
#'   `item` is `"setup|channel|pop"`; `delta_squeezed` is the rise in squeezed
#'   share against the reference setup. The screen's own flag is in `reason`;
#'   the status still follows `limits` only.
#'
#' For both, the value is the donor median and `reason` gives the number of
#' donors and their range.
#' * `event_channel` (only with `event_metrics`): from
#'   [sw_unmix_event_metrics()], `item` is `"file|channel"`; metrics
#'   `pct_below_cut`, `pct_below_floor`, `neg_median`, `neg_rsd`, `spread_q`,
#'   `pileup`, `hyper_share`, `si_mad`, `si_q`.
#' * `event_pair` (only with `event_metrics` that has pairs): `item` is
#'   `"file|channel_i|channel_j"`; `abs_slope` and `abs_rho` of the line check.
#'   These rows go to checks.csv and event_pairs.csv, not to the HTML table.
#' * `heldout` (only with `heldout`): from [sw_unmix_heldout()], `item` is
#'   `"control|channel"`; `heldout_spread_q` and `heldout_bias` (linear).
#'
#' `limits` names and the metric each judges. `max_*` limits fail when the
#' value is larger than the limit; `min_*` limits fail when it is smaller. A
#' value equal to the limit passes.
#'
#' | limit | metric |
#' |---|---|
#' | `max_cosine` | `max_cosine` |
#' | `max_sif` | `max_sif` |
#' | `max_condition_number` | `condition_number` |
#' | `max_cv` | `cv` |
#' | `max_pct_negative` | `pct_negative` |
#' | `max_control_errors` | `n_control_errors` |
#' | `max_control_warnings` | `n_control_warnings` |
#' | `min_events` | `n_events` |
#' | `min_af_spectra` | `n_af_spectra` |
#' | `min_variants` | `n_variants` |
#' | `max_share_lost` | `share_lost` |
#' | `min_median_recovery` | `median_recovery` |
#' | `max_share_moved` | `share_moved` |
#' | `max_delta_squeezed` | `delta_squeezed` |
#' | `max_pct_below_cut` | `pct_below_cut` |
#' | `max_pct_below_floor` | `pct_below_floor` |
#' | `max_spread_q` | `spread_q` |
#' | `max_hyper_share` | `hyper_share` |
#' | `min_stain_index` | `si_mad` |
#' | `min_stain_index_q` | `si_q` |
#' | `max_abs_slope` | `abs_slope` |
#' | `max_abs_rho` | `abs_rho` |
#' | `max_heldout_spread` | `heldout_spread_q` |
#' | `max_heldout_bias` | `heldout_bias` |
#'
#' Each element is `c(warn = , fail = )`, either one alone, or a single
#' unnamed number, which is read as the `fail` limit. Unknown names error.
#'
#' @param result An `sw_unmix_result`.
#' @param limits Named list of limits (see Details). Empty by default, so
#'   nothing is judged.
#' @param unmixed_files Optional character vector of unmixed FCS paths. It
#'   replaces `result$unmixed_files` for this call.
#' @param unmixed_data Optional named list of numeric matrices (events x
#'   channels, with column names), one per file. It takes precedence over
#'   files and needs no flowCore.
#' @param spike Optional `sw_spike_score` from [sw_unmix_spike_score()].
#' @param screen Optional `sw_marker_screen` from [sw_unmix_marker_screen()].
#' @param event_metrics Optional `sw_event_metrics` from
#'   [sw_unmix_event_metrics()].
#' @param heldout Optional `sw_heldout` from [sw_unmix_heldout()].
#' @param ... Passed to [sw_unmix_check_controls()].
#'
#' @return An object of class `sw_unmix_diagnostics`: a list with `run_id`,
#'   `checks` (tibble: `check_id`, `scope`, `item`, `metric`, `value`,
#'   `limit`, `status`, `reason`), `summary` (tibble of counts for all five
#'   statuses), `limits` (as parsed), `notes`, and the `spike` and `screen`
#'   inputs (`NULL` when not given) for [sw_unmix_report()].
#'
#' @examples
#' spectra <- matrix(
#'   c(1, 0.5, 0.1, 0, 0.1, 1, 0.4, 0.1, 0, 0.2, 1, 0.6),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("BV421", "PE", "APC"), paste0("D", 1:4))
#' )
#' res <- sw_unmix_result(spectra = spectra, run_id = "demo")
#' set.seed(1)
#' dat <- list(s1 = matrix(rnorm(300, 100, 30), ncol = 3,
#'                         dimnames = list(NULL, paste0(rownames(spectra), "-A"))))
#' # Without limits every value is "not_judged".
#' sw_unmix_diagnose(res, unmixed_data = dat)
#' # With a limit you chose (uncalibrated example values).
#' dx <- sw_unmix_diagnose(res, limits = list(max_cosine = c(warn = 0.5, fail = 0.9)))
#' dx$checks[dx$checks$metric == "max_cosine", ]
#'
#' @export
sw_unmix_diagnose <- function(
  result,
  limits = list(),
  unmixed_files = NULL,
  unmixed_data = NULL,
  spike = NULL,
  screen = NULL,
  event_metrics = NULL,
  heldout = NULL,
  ...
) {
  if (!inherits(result, "sw_unmix_result")) {
    stop("`result` must be an sw_unmix_result.", call. = FALSE)
  }
  lims <- .sw_parse_limits(limits)
  if (
    !is.null(unmixed_files) &&
      (!is.character(unmixed_files) || anyNA(unmixed_files))
  ) {
    stop(
      "`unmixed_files` must be a character vector of paths.",
      call. = FALSE
    )
  }
  if (!is.null(unmixed_data)) {
    ok <- is.list(unmixed_data) &&
      length(unmixed_data) > 0L &&
      !is.null(names(unmixed_data)) &&
      all(nzchar(names(unmixed_data))) &&
      !anyDuplicated(names(unmixed_data)) &&
      all(vapply(
        unmixed_data,
        \(m) is.matrix(m) && is.numeric(m) && !is.null(colnames(m)),
        logical(1)
      ))
    if (!ok) {
      stop(
        "`unmixed_data` must be a uniquely named list of numeric matrices with column names.",
        call. = FALSE
      )
    }
  }
  files <- if (is.null(unmixed_files)) result$unmixed_files else unmixed_files
  fluors <- rownames(result$spectra)

  checks <- dplyr::bind_rows(
    .sw_diag_controls(result$setup, list(...)),
    .sw_diag_similarity(result$spectra),
    .sw_diag_sif(result$spectra),
    .sw_diag_condition(result$spectra),
    .sw_diag_af(result$af_spectra),
    .sw_diag_variants(result$variants, fluors),
    .sw_diag_files(files, unmixed_data, fluors),
    .sw_diag_spike(spike),
    .sw_diag_screen(screen),
    .sw_diag_event_metrics(event_metrics),
    .sw_diag_heldout(heldout)
  )
  checks <- .sw_apply_limits(checks, lims)

  counts <- table(factor(checks$status, levels = .sw_unmix_statuses))
  summary_tbl <- tibble::tibble(
    status = .sw_unmix_statuses,
    n = as.integer(counts)
  )

  notes <- c(
    "No cutoff is built in: pass, warn and fail appear only for metrics with a supplied limit.",
    "Values are not calibrated against ground truth.",
    "Condition number: max/min singular value from base svd(), unrounded (AutoSpectral 1.8.8 rounds to 2 decimals).",
    if (nrow(lims) == 0L) {
      "No limits were supplied, so every computed value is not_judged."
    }
  )

  structure(
    list(
      run_id = result$run_id,
      checks = checks,
      summary = summary_tbl,
      limits = lims,
      notes = c(notes, if (!is.null(spike)) spike$notes),
      spike = spike,
      screen = screen,
      event_metrics = event_metrics,
      heldout = heldout
    ),
    class = "sw_unmix_diagnostics"
  )
}

#' @export
print.sw_unmix_diagnostics <- function(x, ...) {
  cat("<sw_unmix_diagnostics>", x$run_id, "\n")
  cat("  checks:", nrow(x$checks), "rows\n")
  cat(
    "  ",
    paste(sprintf("%s %d", x$summary$status, x$summary$n), collapse = " | "),
    "\n",
    sep = ""
  )
  if (nrow(x$limits) == 0L) {
    cat("  no limits supplied; nothing judged\n")
  }
  invisible(x)
}


# ---- Report: internal helpers ------------------------------------------------

.sw_html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

# Wrapped so tests can mock a missing jsonlite.
.sw_has_jsonlite <- function() requireNamespace("jsonlite", quietly = TRUE)

# Readable cell text: numbers to 4 significant digits.
.sw_cell_text <- function(x) {
  if (is.numeric(x)) {
    # Per element, so 2e6 events print as 2000000, not 2e+06.
    out <- vapply(
      x,
      \(v) trimws(formatC(signif(v, 4), format = "fg", digits = 4)),
      character(1)
    )
    out[is.na(x)] <- "NA"
    return(out)
  }
  as.character(x)
}

.sw_html_table <- function(df) {
  if (nrow(df) == 0L) {
    return("<p><em>None.</em></p>")
  }
  head <- paste0(
    "<tr>",
    paste0("<th>", .sw_html_escape(names(df)), "</th>", collapse = ""),
    "</tr>"
  )
  cells <- lapply(df, \(col) .sw_html_escape(.sw_cell_text(col)))
  body <- vapply(
    seq_len(nrow(df)),
    \(i) {
      paste0(
        "<tr>",
        paste0(
          "<td>",
          vapply(cells, `[`, character(1), i),
          "</td>",
          collapse = ""
        ),
        "</tr>"
      )
    },
    character(1)
  )
  paste0("<table>", head, paste(body, collapse = ""), "</table>")
}

.sw_html_list <- function(items) {
  paste0(
    "<ul>",
    paste0("<li>", .sw_html_escape(items), "</li>", collapse = ""),
    "</ul>"
  )
}

# Is `child` equal to or inside `parent`? Both normalised, separator appended.
.sw_path_within <- function(child, parent) {
  # paste0() turns a zero-length input into "/", which would match everything.
  if (length(child) == 0L || length(parent) == 0L) {
    return(logical())
  }
  child <- paste0(sub("/+$", "", child), "/")
  parent <- paste0(sub("/+$", "", parent), "/")
  startsWith(child, parent)
}

# Provenance as long key-value rows for provenance.csv.
.sw_provenance_long <- function(prov, run_id) {
  pk <- prov$packages
  args_txt <- if (is.null(prov$args)) {
    character()
  } else {
    vapply(prov$args, \(a) paste(deparse(a), collapse = " "), character(1))
  }
  tibble::tibble(
    section = c(
      "run",
      rep("system", 3L),
      rep("package_version", nrow(pk)),
      rep("package_remote_sha", nrow(pk)),
      rep("backend", 2L),
      rep("baseline", 4L),
      "seed",
      rep("args", length(args_txt))
    ),
    field = c(
      "run_id",
      "r_version",
      "platform",
      "datetime",
      pk$package,
      pk$package,
      "autospectralrcpp_installed",
      "backend",
      "reference",
      "autospectral_commit",
      "matches_reference",
      "message",
      "seed",
      names(args_txt)
    ),
    value = c(
      run_id,
      prov$r_version,
      prov$platform,
      prov$datetime,
      pk$version,
      pk$remote_sha,
      as.character(prov$autospectralrcpp_installed),
      prov$backend,
      prov$baseline$reference,
      prov$autospectral_commit,
      as.character(prov$baseline$matches_reference),
      prov$baseline$message,
      if (is.null(prov$seed)) NA_character_ else as.character(prov$seed),
      unname(args_txt)
    )
  )
}

.sw_report_css <- paste(
  "body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;",
  "padding:0 16px;color:#1a1a1a;background:#fff;line-height:1.45}",
  "h1{font-size:1.6em}h2{font-size:1.25em;border-bottom:1px solid #ccc;",
  "padding-bottom:.2em;margin-top:1.8em}",
  "table{border-collapse:collapse;margin:.6em 0;font-size:.9em}",
  "th,td{border:1px solid #ccc;padding:.25em .5em;text-align:left;",
  "vertical-align:top}th{background:#f2f2f2}",
  ".ok{color:#1b5e20}.bad{color:#b71c1c}code{font-size:.95em}"
)


# ---- Report ------------------------------------------------------------------

#' Write an unmixing report
#'
#' Writes a self-contained HTML page and CSV tables for an unmixing run to a
#' new folder. The report has tables only; it draws no figures. Sections
#' "Spike-in recovery" and "Marker screen" appear when the diagnostics carry
#' a spike score or a marker screen.
#'
#' @details
#' Files written to `path`:
#'
#' * `report.html`: one page with inline CSS, no JavaScript and no external
#'   resources. All values are HTML-escaped; numbers show 4 significant
#'   digits.
#' * `checks.csv`: the diagnostics table at full precision.
#' * `spectra.csv`: the spectra matrix in long form (`fluorophore`,
#'   `detector`, `value`).
#' * `provenance.csv`: provenance as `section`, `field`, `value` rows.
#' * `files.csv`: the file inventory (role, path, size, MD5).
#' * `spike_summary.csv` and `baseline_check.csv`: only when `diagnostics`
#'   was made with `spike`. The per-donor recovery summary and the
#'   independence check of [sw_unmix_spike_score()].
#' * `marker_screen.csv`: only when `diagnostics` was made with `screen`. The
#'   flags table of [sw_unmix_marker_screen()].
#' * `report.json`: everything above, only when jsonlite is installed.
#'   Otherwise it is skipped and the HTML says so.
#'
#' Safety: `path` has no default. It must not exist unless
#' `overwrite = TRUE`; then only the files above are replaced. `path` is
#' refused when it equals, sits inside, or contains the folder of any file in
#' the inventory, or the setup's control folder, so the report cannot
#' overwrite an input.
#'
#' When `diagnostics` is `NULL`, [sw_unmix_diagnose()] runs with no limits, so
#' every value is `not_judged`.
#'
#' @param result An `sw_unmix_result`.
#' @param diagnostics Optional `sw_unmix_diagnostics` for the same run.
#' @param path Folder to create for the report. Required.
#' @param title Optional page title.
#' @param overwrite Write into an existing `path`. Default `FALSE`.
#' @param notes Optional character vector, shown under Limitations. Use it to
#'   state the real-data run status (for example, "Not yet run on real FCS
#'   data").
#'
#' @return The paths of the files written, invisibly, as a named character
#'   vector.
#'
#' @examples
#' spectra <- matrix(
#'   c(1, 0.5, 0.1, 0, 0.1, 1, 0.4, 0.1, 0, 0.2, 1, 0.6),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("BV421", "PE", "APC"), paste0("D", 1:4))
#' )
#' res <- sw_unmix_result(spectra = spectra, run_id = "demo")
#' out <- sw_unmix_report(res, path = file.path(tempdir(), "unmix_report_demo"),
#'                        notes = "Example on mock spectra only.")
#' basename(out)
#' unlink(dirname(out[[1]]), recursive = TRUE)
#'
#' @export
sw_unmix_report <- function(
  result,
  diagnostics = NULL,
  path,
  title = NULL,
  overwrite = FALSE,
  notes = NULL
) {
  if (!inherits(result, "sw_unmix_result")) {
    stop("`result` must be an sw_unmix_result.", call. = FALSE)
  }
  if (missing(path)) {
    stop(
      "`path` has no default. Choose a new folder for the report.",
      call. = FALSE
    )
  }
  if (
    !is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)
  ) {
    stop("`path` must be a single non-empty string.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(title) && (!is.character(title) || length(title) != 1L)) {
    stop("`title` must be NULL or a single string.", call. = FALSE)
  }
  if (!is.null(notes) && !is.character(notes)) {
    stop("`notes` must be NULL or a character vector.", call. = FALSE)
  }
  if (is.null(diagnostics)) {
    diagnostics <- sw_unmix_diagnose(result)
  }
  if (!inherits(diagnostics, "sw_unmix_diagnostics")) {
    stop("`diagnostics` must be an sw_unmix_diagnostics.", call. = FALSE)
  }
  if (!identical(diagnostics$run_id, result$run_id)) {
    stop("`diagnostics` belongs to a different run_id.", call. = FALSE)
  }

  # Refuse a report folder that overlaps any input or output folder.
  prov <- result$provenance
  target <- normalizePath(path, mustWork = FALSE)
  guard_dirs <- unique(c(
    dirname(prov$files$path),
    if (!is.null(result$setup)) result$setup[["control_dir"]]
  ))
  guard_dirs <- normalizePath(guard_dirs, mustWork = FALSE)
  clash <- guard_dirs[
    .sw_path_within(target, guard_dirs) | .sw_path_within(guard_dirs, target)
  ]
  if (length(clash) > 0L) {
    stop(
      sprintf(
        "`path` overlaps the folder of an input or output file: %s.",
        paste(clash, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (file.exists(path) && !overwrite) {
    stop(
      sprintf("`path` already exists: %s. Use overwrite = TRUE.", path),
      call. = FALSE
    )
  }
  if (file.exists(path) && !dir.exists(path)) {
    stop(sprintf("`path` is a file, not a folder: %s.", path), call. = FALSE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  # Remove only this function's own files, so stale ones cannot survive.
  unlink(file.path(path, .sw_unmix_report_files))

  spectra <- result$spectra
  out <- c(
    checks = file.path(path, "checks.csv"),
    spectra = file.path(path, "spectra.csv"),
    provenance = file.path(path, "provenance.csv"),
    files = file.path(path, "files.csv")
  )

  # Tables for the CSV files and the page.
  spectra_long <- tibble::tibble(
    fluorophore = rep(rownames(spectra), times = ncol(spectra)),
    detector = rep(colnames(spectra), each = nrow(spectra)),
    value = as.vector(spectra)
  )
  prov_long <- .sw_provenance_long(prov, result$run_id)
  checks <- diagnostics$checks

  utils::write.csv(checks, out[["checks"]], row.names = FALSE)
  utils::write.csv(spectra_long, out[["spectra"]], row.names = FALSE)
  utils::write.csv(prov_long, out[["provenance"]], row.names = FALSE)
  utils::write.csv(prov$files, out[["files"]], row.names = FALSE)
  spike <- diagnostics$spike
  screen <- diagnostics$screen
  if (!is.null(spike)) {
    out[["spike_summary"]] <- file.path(path, "spike_summary.csv")
    out[["baseline_check"]] <- file.path(path, "baseline_check.csv")
    utils::write.csv(spike$summary, out[["spike_summary"]], row.names = FALSE)
    utils::write.csv(
      spike$baseline_check,
      out[["baseline_check"]],
      row.names = FALSE
    )
  }
  if (!is.null(screen)) {
    out[["marker_screen"]] <- file.path(path, "marker_screen.csv")
    utils::write.csv(screen$flags, out[["marker_screen"]], row.names = FALSE)
  }
  em <- diagnostics$event_metrics
  ho <- diagnostics$heldout
  if (!is.null(em)) {
    out[["event_channels"]] <- file.path(path, "event_channels.csv")
    utils::write.csv(em$channels, out[["event_channels"]], row.names = FALSE)
    if (!is.null(em$pairs)) {
      out[["event_pairs"]] <- file.path(path, "event_pairs.csv")
      utils::write.csv(em$pairs, out[["event_pairs"]], row.names = FALSE)
    }
  }
  if (!is.null(ho)) {
    out[["heldout_scores"]] <- file.path(path, "heldout_scores.csv")
    utils::write.csv(ho$scores, out[["heldout_scores"]], row.names = FALSE)
    if (!is.null(ho$repro)) {
      out[["heldout_repro"]] <- file.path(path, "heldout_repro.csv")
      utils::write.csv(ho$repro, out[["heldout_repro"]], row.names = FALSE)
    }
  }

  has_json <- .sw_has_jsonlite()
  if (has_json) {
    out[["json"]] <- file.path(path, "report.json")
    jsonlite::write_json(
      list(
        run_id = result$run_id,
        provenance = prov_long,
        files = prov$files,
        spectra = spectra_long,
        checks = checks,
        summary = diagnostics$summary,
        notes = c(diagnostics$notes, notes),
        spike_summary = spike$summary,
        baseline_check = spike$baseline_check,
        marker_screen = screen$flags,
        event_channels = em$channels,
        heldout_scores = ho$scores,
        heldout_repro = ho$repro
      ),
      out[["json"]],
      digits = NA,
      auto_unbox = TRUE,
      pretty = TRUE,
      na = "null"
    )
  }
  out[["html"]] <- file.path(path, "report.html")

  # Spectra summary: peak detector and closest other spectrum.
  sim <- checks[checks$check_id == "spectra_similarity", ]
  spec_tbl <- tibble::tibble(
    fluorophore = rownames(spectra),
    peak_detector = colnames(spectra)[apply(spectra, 1L, which.max)],
    max_cosine_to_other = sim$value[match(rownames(spectra), sim$item)],
    note = sim$reason[match(rownames(spectra), sim$item)]
  )
  cond <- checks[checks$check_id == "condition_number", ]
  af_tbl <- checks[
    checks$check_id == "af_library",
    c("item", "metric", "value", "status", "reason")
  ]
  var_tbl <- checks[
    checks$check_id == "variants",
    c("item", "value", "status", "reason")
  ]
  names(var_tbl) <- c("fluorophore", "n_variants", "status", "reason")

  bl <- prov$baseline
  page_title <- if (is.null(title)) "Spectral unmixing report" else title
  limitations <- c(
    "Thresholds are uncalibrated. Any limit used here was supplied by the user.",
    "The diagnostics are not validated against ground truth.",
    paste0(
      "The baseline is unmodified ",
      bl$reference,
      "; extensions must be scored against it."
    ),
    if (is.null(notes)) {
      "Real-data run status: not stated by the caller."
    } else {
      paste("Real-data run status:", notes)
    }
  )

  html <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\"><head><meta charset=\"utf-8\">",
    sprintf("<title>%s</title>", .sw_html_escape(page_title)),
    sprintf("<style>%s</style>", .sw_report_css),
    "</head><body>",
    sprintf("<h1>%s</h1>", .sw_html_escape(page_title)),
    sprintf(
      "<p>Run id: <code>%s</code>. Method: %s.</p>",
      .sw_html_escape(result$run_id),
      .sw_html_escape(
        if (is.null(result$method)) "not recorded" else result$method
      )
    ),
    "<h2>Provenance and reference check</h2>",
    sprintf(
      "<p class=\"%s\">matches_reference: <strong>%s</strong>. %s</p>",
      if (bl$matches_reference) "ok" else "bad",
      bl$matches_reference,
      .sw_html_escape(bl$message)
    ),
    .sw_html_table(prov_long),
    "<h2>Inputs</h2>",
    .sw_html_table(prov$files),
    "<h2>Spectra summary</h2>",
    sprintf(
      "<p>%d fluorophores x %d detectors. Condition number: %s. %s</p>",
      nrow(spectra),
      ncol(spectra),
      .sw_html_escape(.sw_cell_text(cond$value)),
      .sw_html_escape(cond$reason)
    ),
    .sw_html_table(spec_tbl),
    "<h2>AF library</h2>",
    .sw_html_table(af_tbl),
    "<h2>Variants</h2>",
    .sw_html_table(var_tbl),
    if (!is.null(spike)) {
      c(
        "<h2>Spike-in recovery</h2>",
        sprintf(
          "<p>Design <code>%s</code> (%s). Loss cut %s; independence tolerance %s arcsinh. Item is setup|marker|pop; value is the donor median.</p>",
          .sw_html_escape(spike$design_id),
          .sw_html_escape(spike$type),
          .sw_html_escape(format(spike$loss_cut)),
          .sw_html_escape(format(spike$baseline_tol))
        ),
        .sw_html_table(checks[
          checks$check_id == "spike_recovery",
          c("item", "metric", "value", "limit", "status", "reason")
        ]),
        "<p>File status (independence check):</p>",
        .sw_html_table(spike$file_status),
        "<p>Per-donor summary by level: spike_summary.csv.</p>"
      )
    },
    if (!is.null(screen)) {
      c(
        "<h2>Marker screen</h2>",
        sprintf(
          "<p>Reference setup <code>%s</code>. Rule: %s.</p>",
          .sw_html_escape(screen$reference),
          .sw_html_escape(screen$rule$text)
        ),
        .sw_html_table(screen$flags)
      )
    },
    .sw_report_event_sections(em, ho),
    "<h2>Diagnostics</h2>",
    "<p><strong>Verdicts appear only where you supplied a limit; no cutoff is built in.</strong></p>",
    .sw_html_table(diagnostics$summary),
    .sw_html_list(diagnostics$notes),
    .sw_html_table(checks[checks$check_id != "event_pair", ]),
    if (any(checks$check_id == "event_pair")) {
      "<p>Pair rows (event_pair) are in checks.csv and event_pairs.csv only.</p>"
    },
    "<h2>Limitations</h2>",
    .sw_html_list(limitations),
    "<h2>Output files</h2>",
    .sw_html_list(basename(out)),
    if (!has_json) {
      "<p>report.json was not written because jsonlite is not installed.</p>"
    },
    "</body></html>"
  )
  writeLines(html, out[["html"]], useBytes = TRUE)

  invisible(out)
}
