#' @title Thin Wrappers for AutoSpectral Unmixing Diagnostics
#'
#' @description
#' Wrappers around the diagnostic functions that AutoSpectral exports. Each
#' wrapper checks its arguments, checks that AutoSpectral is installed, calls
#' the AutoSpectral function, and returns a tibble where AutoSpectral returns a
#' matrix or data frame. No AutoSpectral code is copied here.
#'
#' @details
#' Signatures and return values were read from AutoSpectral version 1.8.8.
#' AutoSpectral is AGPL licensed and sits in Suggests.
#'
#' Orientation follows the package convention. Spectra sit in rows and
#' detectors sit in columns. Every matrix needs row names and column names.
#' Input is checked before the dependency check, so argument errors show up
#' even when AutoSpectral is missing.
#'
#' Thresholds are never defaulted here. Where AutoSpectral defines a default,
#' the wrapper passes nothing and AutoSpectral applies its own default. The
#' default values are listed in each help page.
#'
#' @name unmix_assess
#' @keywords internal
NULL


# ---- Internal helpers --------------------------------------------------------

.sw_assess_require_autospectral <- function() {
  if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
    stop("Package 'AutoSpectral' is required.", call. = FALSE)
  }
  invisible(TRUE)
}

# Numeric matrix with row and column names (spectra x detectors).
.sw_assess_check_matrix <- function(x, arg, min_rows = 1L) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(
      sprintf(
        "`%s` must be a numeric matrix (spectra in rows, detectors in columns).",
        arg
      ),
      call. = FALSE
    )
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(
      sprintf(
        "`%s` must have row names (spectra) and column names (detectors).",
        arg
      ),
      call. = FALSE
    )
  }
  if (nrow(x) < min_rows) {
    stop(
      sprintf(
        "`%s` needs at least %d row(s); it has %d.",
        arg,
        min_rows,
        nrow(x)
      ),
      call. = FALSE
    )
  }
  invisible(x)
}

# Named list of numeric matrices with detector names (row 1 = reference).
.sw_assess_check_variants <- function(x, arg) {
  if (
    !is.list(x) || length(x) == 0L || is.null(names(x)) || any(names(x) == "")
  ) {
    stop(
      sprintf("`%s` must be a named list of spectral variant matrices.", arg),
      call. = FALSE
    )
  }
  ok <- vapply(
    x,
    \(m) is.matrix(m) && is.numeric(m) && !is.null(colnames(m)),
    logical(1)
  )
  if (!all(ok)) {
    stop(
      sprintf(
        "`%s` elements must be numeric matrices with detector column names; check: %s.",
        arg,
        paste(names(x)[!ok], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(x)
}

.sw_assess_check_scalar_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be a single non-empty string.", arg), call. = FALSE)
  }
  invisible(x)
}

.sw_assess_check_file <- function(x, arg) {
  .sw_assess_check_scalar_string(x, arg)
  if (!file.exists(x)) {
    stop(sprintf("`%s` does not exist: %s", arg, x), call. = FALSE)
  }
  invisible(x)
}

# A spectrum as a named numeric vector or a single-row named matrix.
.sw_assess_check_spectrum <- function(x, spectra, arg) {
  nms <- if (is.matrix(x)) {
    if (nrow(x) != 1L) {
      stop(
        sprintf("`%s` must be a named vector or a single-row matrix.", arg),
        call. = FALSE
      )
    }
    colnames(x)
  } else {
    names(x)
  }
  if (!is.numeric(x) || is.null(nms)) {
    stop(
      sprintf("`%s` must be numeric with detector names.", arg),
      call. = FALSE
    )
  }
  if (!all(nms %in% colnames(spectra))) {
    stop(
      sprintf(
        "`%s` has detector names that are not columns of `spectra`.",
        arg
      ),
      call. = FALSE
    )
  }
  invisible(x)
}

# Call an AutoSpectral function with named args plus pass-through dots.
.sw_assess_call <- function(fn, args, dots) {
  do.call(getExportedValue("AutoSpectral", fn), c(args, dots))
}


# ---- Mismatch ----------------------------------------------------------------

#' Spectral mismatch between reference and test variant lists
#'
#' Wraps `AutoSpectral::assess.mismatch()` and
#' `AutoSpectral::assess.mismatch.angle()`. For each fluorophore in both lists,
#' it compares row 1 of the reference matrix with row 1 of the test matrix.
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8.
#'
#' Orientation: each element of `reference_variants` and `test_variants` is a
#' matrix with variants in rows (row 1 is the reference spectrum) and detectors
#' in columns. Both lists need fluorophore names. AutoSpectral aligns the two
#' lists by detector name and warns if the detector sets differ. These lists
#' come from `AutoSpectral::get.spectral.variants()` (`variants` element).
#'
#' This wrapper sets no cutoff. Judging a "bad" cosine or angle is up to you.
#'
#' @param reference_variants Named list of variant matrices (reference).
#' @param test_variants Named list of variant matrices (to compare).
#' @param metric `"cosine"` for cosine similarity or `"angle"` for spectral
#'   angle in degrees. Required, no default.
#' @param ... Passed to the AutoSpectral function. Neither function takes extra
#'   arguments in 1.8.8.
#'
#' @return A tibble with columns `fluorophore` and `cosine_similarity` (or
#'   `angle_degrees` when `metric = "angle"`).
#'
#' @examples
#' \dontrun{
#' vars <- AutoSpectral::get.spectral.variants(...)$variants
#' sw_unmix_assess_mismatch(vars_cells, vars_beads, metric = "angle")
#' }
#'
#' @keywords internal
sw_unmix_assess_mismatch <- function(
  reference_variants,
  test_variants,
  metric,
  ...
) {
  .sw_assess_check_variants(reference_variants, "reference_variants")
  .sw_assess_check_variants(test_variants, "test_variants")
  if (
    missing(metric) ||
      !is.character(metric) ||
      length(metric) != 1L ||
      !metric %in% c("cosine", "angle")
  ) {
    stop("`metric` must be \"cosine\" or \"angle\".", call. = FALSE)
  }
  if (
    length(intersect(names(reference_variants), names(test_variants))) == 0L
  ) {
    stop(
      "`reference_variants` and `test_variants` share no fluorophore names.",
      call. = FALSE
    )
  }
  .sw_assess_require_autospectral()

  fn <- if (metric == "cosine") "assess.mismatch" else "assess.mismatch.angle"
  out <- .sw_assess_call(
    fn,
    list(
      reference.variants = reference_variants,
      test.variants = test_variants
    ),
    list(...)
  )

  value_name <- if (metric == "cosine") "cosine_similarity" else "angle_degrees"
  res <- tibble::tibble(fluorophore = rownames(out), value = out[, 1L])
  names(res)[2L] <- value_name
  res
}

#' Cluster-based permutation test for mismatch regions
#'
#' Wraps `AutoSpectral::assess.mismatch.clusters()`. It finds runs of adjacent
#' detectors where the bead-vs-cell mismatch is larger than chance.
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8.
#'
#' Orientation: `dist_mat` and `ref_spectra` are fluorophores in rows and
#' detectors in columns. They must share row and column names. Column order
#' must be true detector order, because contiguity is defined by that order.
#'
#' Tuning values are passed through `...` under AutoSpectral's own names.
#' This wrapper does not set them. AutoSpectral's defaults apply when you leave
#' them out: `peak.threshold = 0.1`, `min.n = 3L`, `cluster.threshold = 1.96`,
#' `n.perm = 2000L`, `bird.seed = NULL` (no seed, so results vary between
#' runs). Pass `bird.seed` for a reproducible p-value.
#'
#' @param dist_mat Numeric matrix of signed per-detector differences
#'   (fluorophore x detector).
#' @param ref_spectra Numeric matrix of reference intensities
#'   (fluorophore x detector).
#' @param ... Passed to `AutoSpectral::assess.mismatch.clusters()`.
#'
#' @return A tibble with one row per cluster and columns `start`, `end`,
#'   `n.detectors`, `mass` and `p.value`. It has zero rows if no cluster
#'   formed. The full AutoSpectral list (`stat`, `mask`, `n.detector`,
#'   `clusters`, `null.max.mass`) is in `attr(x, "autospectral")`.
#'
#' @examples
#' \dontrun{
#' sw_unmix_assess_mismatch_clusters(dist_mat, ref_spectra, bird.seed = 1L)
#' }
#'
#' @keywords internal
sw_unmix_assess_mismatch_clusters <- function(dist_mat, ref_spectra, ...) {
  .sw_assess_check_matrix(dist_mat, "dist_mat")
  .sw_assess_check_matrix(ref_spectra, "ref_spectra")
  if (
    length(intersect(rownames(dist_mat), rownames(ref_spectra))) == 0L ||
      length(intersect(colnames(dist_mat), colnames(ref_spectra))) == 0L
  ) {
    stop(
      "`dist_mat` and `ref_spectra` share no row or column names.",
      call. = FALSE
    )
  }
  .sw_assess_require_autospectral()

  raw <- .sw_assess_call(
    "assess.mismatch.clusters",
    list(dist.mat = dist_mat, ref.spectra = ref_spectra),
    list(...)
  )
  out <- tibble::as_tibble(raw$clusters)
  attr(out, "autospectral") <- raw
  out
}


# ---- Variability -------------------------------------------------------------

#' Within-fluorophore spectral variability
#'
#' Wraps `AutoSpectral::assess.variability()` (one summary row per
#' fluorophore) and `AutoSpectral::assess.variability.mad()` (a per-detector
#' MAD profile).
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8.
#'
#' Orientation: each element of `variant_list` is a matrix with variants in
#' rows (row 1 is the reference spectrum, rows 2 and up are variants) and
#' detectors in columns. The list is named by fluorophore.
#'
#' `peak_height` is a threshold. AutoSpectral defines a default of `0.1`. When
#' `peak_height` is `NULL`, nothing is passed and AutoSpectral applies `0.1`.
#' Any other value is passed unchanged.
#'
#' @param variant_list Named list of variant matrices, as in
#'   `get.spectral.variants()$variants`.
#' @param what `"summary"` for `assess.variability()` or `"mad"` for
#'   `assess.variability.mad()`. Required, no default.
#' @param peak_height `NULL` (use AutoSpectral's default) or a number in
#'   `[0, 1]`.
#' @param ... Passed to the AutoSpectral function.
#'
#' @return For `"summary"`, a tibble with `fluorophore`, `rSD`, `Scaled`,
#'   `Sim_rSD` and `Denoised_rSD`. For `"mad"`, a tibble with `fluorophore` and
#'   one column per detector. A zero-row tibble comes back if AutoSpectral
#'   returns `NULL` (no usable variability data).
#'
#' @examples
#' \dontrun{
#' sw_unmix_assess_variability(vars, what = "summary")
#' sw_unmix_assess_variability(vars, what = "mad", peak_height = 0.2)
#' }
#'
#' @keywords internal
sw_unmix_assess_variability <- function(
  variant_list,
  what,
  peak_height = NULL,
  ...
) {
  .sw_assess_check_variants(variant_list, "variant_list")
  if (
    missing(what) ||
      !is.character(what) ||
      length(what) != 1L ||
      !what %in% c("summary", "mad")
  ) {
    stop("`what` must be \"summary\" or \"mad\".", call. = FALSE)
  }
  if (
    !is.null(peak_height) &&
      (!is.numeric(peak_height) ||
        length(peak_height) != 1L ||
        is.na(peak_height) ||
        peak_height < 0 ||
        peak_height > 1)
  ) {
    stop(
      "`peak_height` must be NULL or a single number in [0, 1].",
      call. = FALSE
    )
  }
  .sw_assess_require_autospectral()

  args <- list(variant.list = variant_list)
  if (!is.null(peak_height)) {
    args$peak.height <- peak_height
  }
  fn <- if (what == "summary") {
    "assess.variability"
  } else {
    "assess.variability.mad"
  }
  out <- .sw_assess_call(fn, args, list(...))

  if (is.null(out)) {
    return(tibble::tibble(fluorophore = character()))
  }
  tibble::as_tibble(out, rownames = "fluorophore")
}

#' Alignment between variability and mismatch profiles
#'
#' Wraps `AutoSpectral::assess.variability.alignment()`. It gives the cosine
#' similarity between each fluorophore's per-detector MAD profile and the
#' absolute size of its bead-vs-cell mismatch.
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8.
#'
#' Orientation: both matrices have fluorophores in rows and detectors in
#' columns. AutoSpectral aligns them by name and warns if detector sets differ.
#' The result lies in `[0, 1]`. This wrapper sets no cutoff.
#'
#' @param variability_mad Matrix from [sw_unmix_assess_variability()] with
#'   `what = "mad"`, converted back with `as.matrix()` and row names, or from
#'   `AutoSpectral::assess.variability.mad()`.
#' @param mismatch_dist Numeric matrix of signed per-detector mismatch
#'   (fluorophore x detector).
#' @param ... Passed to `AutoSpectral::assess.variability.alignment()`.
#'
#' @return A tibble with `fluorophore` and `alignment`.
#'
#' @examples
#' \dontrun{
#' sw_unmix_assess_alignment(mad_matrix, dist_mat)
#' }
#'
#' @keywords internal
sw_unmix_assess_alignment <- function(variability_mad, mismatch_dist, ...) {
  .sw_assess_check_matrix(variability_mad, "variability_mad")
  .sw_assess_check_matrix(mismatch_dist, "mismatch_dist")
  if (
    length(intersect(rownames(variability_mad), rownames(mismatch_dist))) == 0L
  ) {
    stop(
      "`variability_mad` and `mismatch_dist` share no fluorophore names.",
      call. = FALSE
    )
  }
  .sw_assess_require_autospectral()

  out <- .sw_assess_call(
    "assess.variability.alignment",
    list(variability.mad = variability_mad, mismatch.dist = mismatch_dist),
    list(...)
  )
  tibble::tibble(fluorophore = rownames(out), alignment = out[, 1L])
}


# ---- Condition number --------------------------------------------------------

#' Condition number of a spectral panel
#'
#' Wraps `AutoSpectral::calculate.condition.number()`, the mixing matrix
#' condition number (Complexity Index).
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8.
#'
#' Orientation: `spectra` has fluorophores in rows and detectors in columns,
#' normalised between 0 and 1 (AutoSpectral's documented input). It needs more
#' than one row. This wrapper sets no alarm level for the result.
#'
#' @param spectra Numeric matrix, fluorophores x detectors.
#' @param ... Passed to `AutoSpectral::calculate.condition.number()`.
#'
#' @return A single number.
#'
#' @examples
#' \dontrun{
#' sw_unmix_condition_number(spectra)
#' }
#'
#' @export
sw_unmix_condition_number <- function(spectra, ...) {
  .sw_assess_check_matrix(spectra, "spectra", min_rows = 2L)
  if (anyNA(spectra) || any(!is.finite(spectra))) {
    stop("`spectra` must not contain NA, NaN or Inf.", call. = FALSE)
  }
  .sw_assess_require_autospectral()

  .sw_assess_call(
    "calculate.condition.number",
    list(spectra = spectra),
    list(...)
  )
}


# ---- Control checks ----------------------------------------------------------

#' Check single-stained control set-up
#'
#' Wraps four AutoSpectral checks. The control file check always runs. The
#' duplicate, gate and channel checks run only when you supply their inputs.
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8.
#'
#' * `AutoSpectral::check.control.file()` always runs, with `strict` left at
#'   AutoSpectral's default (`FALSE`), so issues come back as a table.
#'   Thresholds `min.event.warning` (default 5000) and `min.event.error`
#'   (default 1000) are AutoSpectral defaults. Change them through `...`.
#' * `AutoSpectral::check.spectra.duplicates()` runs when `spectra` is given.
#'   It stops on a duplicated fluorophore. Without a `"fluorophore"` attribute
#'   it only warns about rows with cosine similarity above 0.95 (AutoSpectral's
#'   default).
#' * `AutoSpectral::check.gates()` runs when `gate_list` and `control_table`
#'   are both given. It stops on a bad gate.
#' * `AutoSpectral::check.channels()` runs when `spectral_channels` is given.
#'   It reorders names into excitation and emission order.
#'
#' Orientation: `spectra` has fluorophores in rows and detectors in columns.
#'
#' `AutoSpectral::check.consistency()` is not wrapped. It is a one-line helper
#' for gate settings and has no diagnostic output of its own.
#'
#' @param control_dir Directory of single-stained control FCS files.
#' @param control_def_file Path to the control definition CSV.
#' @param asp AutoSpectral parameter list from
#'   `AutoSpectral::get.autospectral.param()`.
#' @param spectra Optional numeric matrix, fluorophores x detectors.
#' @param gate_list Optional named list of gates.
#' @param control_table Control table. Needed with `gate_list`.
#' @param spectral_channels Optional character vector of channel names.
#' @param ... Passed to `AutoSpectral::check.control.file()`.
#'
#' @return A tibble of issues from the control file check. It has zero rows
#'   when AutoSpectral finds none. Attributes: `"spectra_duplicates_ok"`,
#'   `"gates_ok"` (each `TRUE` when run and passed, otherwise `NA`) and
#'   `"channel_order"` (character vector, or `NULL`).
#'
#' @examples
#' \dontrun{
#' asp <- AutoSpectral::get.autospectral.param(cytometer = "aurora")
#' sw_unmix_check_controls("controls/", "control_def.csv", asp, spectra = spectra)
#' }
#'
#' @export
sw_unmix_check_controls <- function(
  control_dir,
  control_def_file,
  asp,
  spectra = NULL,
  gate_list = NULL,
  control_table = NULL,
  spectral_channels = NULL,
  ...
) {
  .sw_assess_check_scalar_string(control_dir, "control_dir")
  if (!dir.exists(control_dir)) {
    stop(
      sprintf("`control_dir` is not a directory: %s", control_dir),
      call. = FALSE
    )
  }
  .sw_assess_check_file(control_def_file, "control_def_file")
  if (!is.list(asp)) {
    stop("`asp` must be the AutoSpectral parameter list.", call. = FALSE)
  }
  if (!is.null(spectra)) {
    .sw_assess_check_matrix(spectra, "spectra", min_rows = 2L)
  }
  if (!is.null(gate_list) != !is.null(control_table)) {
    stop(
      "`gate_list` and `control_table` must be given together.",
      call. = FALSE
    )
  }
  if (!is.null(gate_list) && !is.list(gate_list)) {
    stop("`gate_list` must be a named list of gates.", call. = FALSE)
  }
  if (
    !is.null(spectral_channels) &&
      (!is.character(spectral_channels) || length(spectral_channels) == 0L)
  ) {
    stop(
      "`spectral_channels` must be a character vector of channel names.",
      call. = FALSE
    )
  }
  .sw_assess_require_autospectral()

  issues <- .sw_assess_call(
    "check.control.file",
    list(
      control.dir = control_dir,
      control.def.file = control_def_file,
      asp = asp
    ),
    list(...)
  )
  out <- if (is.data.frame(issues)) {
    tibble::as_tibble(issues)
  } else {
    tibble::tibble()
  }

  attr(out, "spectra_duplicates_ok") <- NA
  attr(out, "gates_ok") <- NA
  attr(out, "channel_order") <- NULL

  if (!is.null(spectra)) {
    attr(out, "spectra_duplicates_ok") <-
      isTRUE(AutoSpectral::check.spectra.duplicates(spectra))
  }
  if (!is.null(gate_list)) {
    attr(out, "gates_ok") <-
      isTRUE(AutoSpectral::check.gates(gate_list, control_table, asp))
  }
  if (!is.null(spectral_channels)) {
    attr(out, "channel_order") <- AutoSpectral::check.channels(
      spectral_channels,
      asp
    )
  }
  out
}


# ---- Comparisons (these draw or save plots) ----------------------------------

#' Compare unmixing with two spectral references
#'
#' Wraps `AutoSpectral::compare.unmix()`. It unmixes one single-stained file
#' with a reference and a test spectrum and tabulates SSI and spillover.
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8.
#'
#' **Side effect:** AutoSpectral draws a biplot to the active graphics device
#' and saves a JPEG to `plot_dir`. AutoSpectral cannot skip the file, so
#' `plot_dir` is required. This wrapper adds no plotting.
#'
#' Orientation: `spectra` is fluorophores in rows and detectors in columns,
#' normalised 0 to 1. `ref_spectrum` and `test_spectrum` are named vectors or
#' single-row matrices whose names are detectors in `spectra`. All other
#' AutoSpectral arguments (`ref.spectra`, `test.spectra`, `gate.bound`,
#' `ref.label`, `test.label`, axis settings) go through `...` with
#' AutoSpectral's defaults.
#'
#' @param single_stained_fcs Path to the single-stained FCS file.
#' @param unstained_fcs Path to the unstained FCS file.
#' @param fluorophore Name of the target fluorophore.
#' @param spectra Numeric matrix, fluorophores x detectors.
#' @param ref_spectrum Reference spectrum.
#' @param test_spectrum Test spectrum.
#' @param cytometer Cytometer name for `AutoSpectral::get.autospectral.param()`.
#' @param ... Passed to `AutoSpectral::compare.unmix()`.
#'
#' @return A tibble with `Fluorophore`, `Reference.SSI`, `Test.SSI`,
#'   `Reference.Spill` and `Test.Spill`.
#'
#' @examples
#' \dontrun{
#' sw_unmix_compare("BV421.fcs", "unstained.fcs", "BV421", spectra,
#'                  ref_spectrum = cell_bv421, test_spectrum = bead_bv421,
#'                  cytometer = "aurora")
#' }
#'
#' @keywords internal
sw_unmix_compare <- function(
  single_stained_fcs,
  unstained_fcs,
  fluorophore,
  spectra,
  ref_spectrum,
  test_spectrum,
  cytometer,
  plot_dir,
  ...
) {
  if (missing(plot_dir)) {
    stop(
      "`plot_dir` has no default: AutoSpectral always saves a JPEG there. ",
      "Choose where it goes.",
      call. = FALSE
    )
  }
  .sw_assess_check_scalar_string(plot_dir, "plot_dir")
  .sw_assess_check_file(single_stained_fcs, "single_stained_fcs")
  .sw_assess_check_file(unstained_fcs, "unstained_fcs")
  .sw_assess_check_scalar_string(fluorophore, "fluorophore")
  .sw_assess_check_matrix(spectra, "spectra")
  .sw_assess_check_spectrum(ref_spectrum, spectra, "ref_spectrum")
  .sw_assess_check_spectrum(test_spectrum, spectra, "test_spectrum")
  .sw_assess_check_scalar_string(cytometer, "cytometer")
  .sw_assess_require_autospectral()

  out <- .sw_assess_call(
    "compare.unmix",
    list(
      single.stained.fcs = single_stained_fcs,
      unstained.fcs = unstained_fcs,
      fluorophore = fluorophore,
      spectra = spectra,
      ref.spectrum = ref_spectrum,
      test.spectrum = test_spectrum,
      cytometer = cytometer,
      plot.dir = plot_dir
    ),
    list(...)
  )
  tibble::as_tibble(out)
}

#' Compare already-unmixed data sets as biplot grids
#'
#' Wraps `AutoSpectral::compare.unmixed.data()`. It is a plotting function.
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8.
#'
#' **Side effect:** AutoSpectral prints the grid. Unlike AutoSpectral, this
#' wrapper writes no file unless you set `output_dir`. This wrapper adds no
#' plotting.
#'
#' Orientation: each data set is events in rows and channels in columns, with
#' column names. Plot settings go through `...` with AutoSpectral's defaults.
#'
#' @param data_list List of two or more matrices or data frames.
#' @param channels Character vector of at least two channel names.
#' @param asp AutoSpectral parameter list.
#' @param output_dir Directory for the saved JPEG. `NULL` (default) saves nothing.
#' @param ... Passed to `AutoSpectral::compare.unmixed.data()`.
#'
#' @return The combined `cowplot` object, invisibly. It is not a tibble.
#'
#' @examples
#' \dontrun{
#' asp <- AutoSpectral::get.autospectral.param(cytometer = "aurora")
#' sw_unmix_compare_data(list(cells = cells, beads = beads),
#'                       channels = c("CD3", "CD4"), asp = asp, save = FALSE)
#' }
#'
#' @keywords internal
sw_unmix_compare_data <- function(
  data_list,
  channels,
  asp,
  output_dir = NULL,
  ...
) {
  if (!is.list(data_list) || length(data_list) < 2L) {
    stop("`data_list` must be a list of two or more data sets.", call. = FALSE)
  }
  has_names <- vapply(
    data_list,
    \(d) (is.matrix(d) || is.data.frame(d)) && !is.null(colnames(d)),
    logical(1)
  )
  if (!all(has_names)) {
    stop(
      "Every element of `data_list` must be a matrix or data frame with column names.",
      call. = FALSE
    )
  }
  if (!is.character(channels) || length(channels) < 2L) {
    stop(
      "`channels` must be a character vector of at least two channel names.",
      call. = FALSE
    )
  }
  shared <- Reduce(intersect, lapply(data_list, colnames))
  if (sum(channels %in% shared) < 2L) {
    stop(
      "Fewer than two `channels` are present in every element of `data_list`.",
      call. = FALSE
    )
  }
  if (!is.list(asp)) {
    stop("`asp` must be the AutoSpectral parameter list.", call. = FALSE)
  }
  .sw_assess_require_autospectral()

  .sw_assess_call(
    "compare.unmixed.data",
    c(
      list(data.list = data_list, channels = channels, asp = asp),
      if (is.null(output_dir)) {
        list(save = FALSE)
      } else {
        list(save = TRUE, output.dir = output_dir)
      }
    ),
    list(...)
  )
}

#' Compare autofluorescence spectra sets
#'
#' Wraps `AutoSpectral::compare.af()`. It scores candidate AF spectra against
#' one unstained FCS file.
#'
#' @details
#' Signatures were read from AutoSpectral version 1.8.8.
#'
#' **Side effect:** AutoSpectral saves a bar chart PDF only if you set
#' `plot_dir`. The default `NULL` writes nothing (AutoSpectral's own default
#' would write to `figure_af_accuracy`). This wrapper adds no plotting.
#'
#' Orientation: `spectra` and every matrix in `af_spectra_list` have spectra
#' in rows and detectors in columns. Detector names of each AF matrix must
#' match `colnames(spectra)`. Row 1 of each AF matrix is its mean AF spectrum.
#' `n.downsample` (default 1000), `assign.fn` and `title` go through `...` with
#' AutoSpectral's defaults. The random subsample is not seeded here.
#'
#' @param unstained_fcs Path to the unstained FCS file.
#' @param spectra Numeric matrix, fluorophores x detectors.
#' @param af_spectra_list List of AF spectra matrices.
#' @param plot_dir Folder for the summary PDF. `NULL` (default) saves nothing.
#' @param ... Passed to `AutoSpectral::compare.af()`.
#'
#' @return A tibble with one row per candidate (plus `"baseline"`) and columns
#'   `candidate`, `mean_sim`, `rsd_sim` and `n_variants`. The full list with
#'   per-cell `Assignments` and `Similarity` is in `attr(x, "autospectral")`.
#'
#' @examples
#' \dontrun{
#' sw_unmix_compare_af("unstained.fcs", spectra, list(k5 = af_k5, k10 = af_k10),
#'                     plot_dir = NULL)
#' }
#'
#' @keywords internal
sw_unmix_compare_af <- function(
  unstained_fcs,
  spectra,
  af_spectra_list,
  plot_dir = NULL,
  ...
) {
  .sw_assess_check_file(unstained_fcs, "unstained_fcs")
  .sw_assess_check_matrix(spectra, "spectra")
  if (!is.list(af_spectra_list) || length(af_spectra_list) == 0L) {
    stop(
      "`af_spectra_list` must be a non-empty list of AF spectra matrices.",
      call. = FALSE
    )
  }
  af_ok <- vapply(
    af_spectra_list,
    \(m) {
      is.matrix(m) && is.numeric(m) && setequal(colnames(m), colnames(spectra))
    },
    logical(1)
  )
  if (!all(af_ok)) {
    stop(
      "Each `af_spectra_list` element must be a numeric matrix whose columns match ",
      "`colnames(spectra)`.",
      call. = FALSE
    )
  }
  .sw_assess_require_autospectral()

  raw <- .sw_assess_call(
    "compare.af",
    list(
      unstained.fcs = unstained_fcs,
      spectra = spectra,
      af.spectra.list = af_spectra_list,
      plot.dir = plot_dir
    ),
    list(...)
  )
  out <- tibble::tibble(
    candidate = names(raw),
    mean_sim = vapply(raw, \(x) x$Mean_Sim, numeric(1)),
    rsd_sim = vapply(raw, \(x) x$rSD_Sim, numeric(1)),
    n_variants = vapply(raw, \(x) as.integer(x$n.variants), integer(1))
  )
  attr(out, "autospectral") <- raw
  out
}
