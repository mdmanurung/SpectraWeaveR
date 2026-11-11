# R/cytoscan_wrap.R
# Wrappers for the GPL (>= 2) packages CytoScan and CytoTools (AUMC-HEMA).
# Nothing is copied: each function calls the installed package and reshapes the
# result into a tibble. Tests replace the lookup helper .sw_cytoscan_fn().

# Fetch an exported function from CytoScan or CytoTools, or stop with install help.
.sw_cytoscan_fn <- function(name, pkg = "CytoScan") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '",
      pkg,
      "' is required. Install it with ",
      "remotes::install_github('AUMC-HEMA/",
      pkg,
      "').",
      call. = FALSE
    )
  }
  getExportedValue(pkg, name)
}

#' Preprocess function that leaves a flowFrame unchanged
#'
#' CytoScan applies compensation and an arcsinh transform when it reads a file,
#' unless it is given another preprocess function. Pass this one when your files
#' are already on the analysis scale, or give your own function.
#'
#' @param ff A `flowFrame`.
#'
#' @return `ff`, unchanged.
#'
#' @examples
#' \dontrun{
#' sw_cytoscan_run(..., preprocess = sw_cytoscan_preprocess_none)
#' }
#'
#' @export
sw_cytoscan_preprocess_none <- function(ff) {
  ff
}

#' Run CytoScan anomaly detection on FCS files
#'
#' @description
#' Wrap the CytoScan package: read files, build features, and flag outliers or
#' novelties. The result is a tibble with one row per test file. This calls
#' CytoScan itself; use [sw_anomaly_features()] and [sw_anomaly_flag_outlier()]
#' for the package's own light versions that need no CytoScan install.
#'
#' @details
#' **No defaults.** Every argument below changes which files are flagged, so
#' none has one. CytoScan's own defaults are conventions (for example, a
#' quantile distance of 0.1 and an outlier cutoff near 0.5), not validated
#' values.
#'
#' **Scale.** CytoScan reads files with a preprocess function. Its built-in
#' one compensates and applies arcsinh with cofactor 150, which suits
#' conventional flow. This wrapper always replaces it with `preprocess`, so
#' you choose the scale. Pass [sw_cytoscan_preprocess_none()] for data that is
#' already transformed.
#'
#' **Aurora files.** `truncate_max_range` must be `FALSE` for Aurora FCS files
#' (intensities near 4e6 exceed flowCore's default range).
#'
#' @param test_paths Character vector of FCS paths to check.
#' @param reference_paths Character vector of reference FCS paths. Needed for
#'   `flag_method = "novelty"` and for aggregate features built from the
#'   reference; use `character(0)` when none.
#' @param channels Character vector of channel names to use.
#' @param feat_method One of `"quantiles"` or `"EMD"` (CytoScan also offers
#'   `"summary"`, `"binning"` and `"fingerprint"`; these are passed through
#'   unchanged if you name them).
#' @param flag_method `"outlier"` or `"novelty"`.
#' @param n Whole number, events read per file.
#' @param quantile_dist Number, quantile spacing for quantile features.
#' @param agg_size Whole number, events in the pooled sample (EMD, binning).
#' @param outlier_threshold Isolation-score cutoff in (0, 1), for outliers.
#' @param novelty_threshold GMM log-density quantile in (0, 1), for novelties.
#' @param ntrees Whole number, trees for the outlier forest.
#' @param preprocess Function taking and returning a `flowFrame`.
#' @param truncate_max_range `TRUE` or `FALSE`, passed to
#'   `flowCore::read.FCS()` through CytoScan.
#' @param cores Whole number of cores for CytoScan's parallel reading.
#'
#' @return A tibble of class `sw_cytoscan` with `file`, `flag` (`TRUE` for a
#'   flagged file) and `kind` (`outlier` or `novelty`). Attribute `params`
#'   echoes the arguments and `features` holds CytoScan's feature table.
#'
#' @seealso [sw_anomaly_features()], [sw_anomaly_flag_outlier()]
#'
#' @examples
#' \dontrun{
#' sw_cytoscan_run(
#'   test_paths = files, reference_paths = character(0),
#'   channels = chans, feat_method = "quantiles", flag_method = "outlier",
#'   n = 1000, quantile_dist = 0.1, agg_size = 10000,
#'   outlier_threshold = 0.6, novelty_threshold = 0.05, ntrees = 1000,
#'   preprocess = sw_cytoscan_preprocess_none,
#'   truncate_max_range = FALSE, cores = 1
#' )
#' }
#'
#' @export
sw_cytoscan_run <- function(
  test_paths,
  reference_paths,
  channels,
  feat_method,
  flag_method,
  n,
  quantile_dist,
  agg_size,
  outlier_threshold,
  novelty_threshold,
  ntrees,
  preprocess,
  truncate_max_range,
  cores
) {
  .sw_stop_missing(c(
    test_paths = missing(test_paths),
    reference_paths = missing(reference_paths),
    channels = missing(channels),
    feat_method = missing(feat_method),
    flag_method = missing(flag_method),
    n = missing(n),
    quantile_dist = missing(quantile_dist),
    agg_size = missing(agg_size),
    outlier_threshold = missing(outlier_threshold),
    novelty_threshold = missing(novelty_threshold),
    ntrees = missing(ntrees),
    preprocess = missing(preprocess),
    truncate_max_range = missing(truncate_max_range),
    cores = missing(cores)
  ))
  if (!is.character(test_paths) || length(test_paths) < 2) {
    stop("`test_paths` must be at least two file paths.", call. = FALSE)
  }
  if (!is.character(reference_paths)) {
    stop("`reference_paths` must be a character vector.", call. = FALSE)
  }
  if (!flag_method %in% c("outlier", "novelty")) {
    stop("`flag_method` must be \"outlier\" or \"novelty\".", call. = FALSE)
  }
  if (flag_method == "novelty" && length(reference_paths) == 0) {
    stop("`flag_method = \"novelty\"` needs `reference_paths`.", call. = FALSE)
  }
  if (!is.function(preprocess)) {
    stop("`preprocess` must be a function.", call. = FALSE)
  }
  if (!is.logical(truncate_max_range) || length(truncate_max_range) != 1) {
    stop("`truncate_max_range` must be TRUE or FALSE.", call. = FALSE)
  }
  .sw_check_count(n, "n")
  .sw_check_count(agg_size, "agg_size")
  .sw_check_count(ntrees, "ntrees")
  .sw_check_count(cores, "cores")
  for (nm in c("outlier_threshold", "novelty_threshold")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1 || v <= 0 || v >= 1) {
      stop("`", nm, "` must be one number in (0, 1).", call. = FALSE)
    }
  }

  cs <- .sw_cytoscan_fn("CytoScan")(truncate_max_range = truncate_max_range)
  # Replace CytoScan's built-in compensation + arcsinh 150 with the caller's step.
  cs[["preprocessFunction"]] <- preprocess
  cs <- .sw_cytoscan_fn("addTestdata")(cs, test_paths)
  if (length(reference_paths) > 0) {
    cs <- .sw_cytoscan_fn("addReferencedata")(cs, reference_paths)
  }
  cs <- .sw_cytoscan_fn("generateFeatures")(
    cs,
    channels = channels,
    featMethod = feat_method,
    n = n,
    quantileDist = quantile_dist,
    aggSize = agg_size,
    cores = cores
  )
  cs <- .sw_cytoscan_fn("Flag")(
    cs,
    featMethod = feat_method,
    flagMethod = flag_method,
    outlier_threshold = outlier_threshold,
    novelty_threshold = novelty_threshold,
    ntrees = ntrees
  )

  flags <- if (flag_method == "outlier") {
    as.logical(as.character(cs$outliers[[feat_method]]))
  } else {
    as.logical(cs$novelties[[feat_method]])
  }
  out <- tibble::tibble(
    file = test_paths,
    flag = flags,
    kind = flag_method
  )
  attr(out, "params") <- list(
    channels = channels,
    feat_method = feat_method,
    flag_method = flag_method,
    n = n,
    quantile_dist = quantile_dist,
    agg_size = agg_size,
    outlier_threshold = outlier_threshold,
    novelty_threshold = novelty_threshold,
    ntrees = ntrees,
    truncate_max_range = truncate_max_range
  )
  attr(out, "features") <- cs$features$test[[feat_method]]
  class(out) <- c("sw_cytoscan", class(out))
  out
}

#' FCS header table from CytoTools
#'
#' Wrap `CytoTools::getFCSmetadata()`, which reads cytometer id, export time and
#' the detector voltage string of each file. For a comparison with a reference
#' file and a long result, use [sw_fcs_header_audit()], which needs no extra
#' package.
#'
#' @param paths Character vector of FCS paths.
#'
#' @return A tibble, one row per file, as returned by CytoTools.
#'
#' @seealso [sw_fcs_header_audit()]
#'
#' @examples
#' \dontrun{
#' sw_cytotools_header(list.files("fcs", full.names = TRUE))
#' }
#'
#' @export
sw_cytotools_header <- function(paths) {
  if (missing(paths)) {
    stop("`paths` is required.", call. = FALSE)
  }
  if (!is.character(paths) || length(paths) == 0) {
    stop("`paths` must be a character vector of file paths.", call. = FALSE)
  }
  tibble::as_tibble(.sw_cytoscan_fn("getFCSmetadata", "CytoTools")(paths))
}
