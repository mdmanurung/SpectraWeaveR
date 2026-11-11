# R/anomaly_flag.R
# Anomaly flags on per-file features: outliers (isolation forest) and novelties
# (Gaussian mixture fitted on reference files). Clean-room implementation from the
# CytoScan paper's written description; no upstream code is used. Both flags call
# installed packages (isotree, mclust) and add only the bookkeeping around them.

#' Flag outlier files with an isolation forest
#'
#' @description
#' Fit an isolation forest on the feature rows of all files and flag files
#' whose anomaly score reaches `cutoff`. An outlier is a file that looks
#' unlike the other files in the same set. No reference set is needed.
#'
#' @details
#' **No defaults.** `ntrees`, `sample_size`, `ndim`, `cutoff` and `seed` each
#' change which files are flagged, so each is required. The isolation forest
#' paper suggests 0.5 as the score above which points are likely anomalous;
#' that is a convention, not a calibrated cutoff. CytoScan reports that
#' raising sensitivity raises the false positive rate. Treat the flag as a
#' prompt to look at a file, not as a verdict.
#'
#' **Scores are the result.** `score` is returned for every file, so you can
#' rank files without any cutoff.
#'
#' @param features An `sw_anomaly_features` tibble, or a data frame with a
#'   `file` column and numeric feature columns.
#' @param ntrees Whole number, trees in the forest.
#' @param sample_size Whole number, rows drawn per tree. Not more than the
#'   number of files.
#' @param ndim Whole number, columns combined per split. 1 is the standard
#'   isolation forest.
#' @param cutoff Number in (0, 1). Files with `score >= cutoff` are flagged.
#' @param seed Whole number. The caller's RNG state is not touched.
#'
#' @return A tibble of class `sw_anomaly_flag` with `file`, `score`, `flag`
#'   (`TRUE` at or above `cutoff`). Attribute `params` echoes the arguments.
#'
#' @seealso [sw_anomaly_features()], [sw_anomaly_flag_novelty()]
#'
#' @examples
#' \dontrun{
#' sw_anomaly_flag_outlier(
#'   feats,
#'   ntrees = 500, sample_size = 32, ndim = 1, cutoff = 0.6, seed = 1
#' )
#' }
#'
#' @export
sw_anomaly_flag_outlier <- function(
  features,
  ntrees,
  sample_size,
  ndim,
  cutoff,
  seed
) {
  .sw_stop_missing(c(
    features = missing(features),
    ntrees = missing(ntrees),
    sample_size = missing(sample_size),
    ndim = missing(ndim),
    cutoff = missing(cutoff),
    seed = missing(seed)
  ))
  if (!requireNamespace("isotree", quietly = TRUE)) {
    stop("Package 'isotree' is required for outlier flags.", call. = FALSE)
  }
  x <- .sw_flag_matrix(features)
  .sw_check_count(ntrees, "ntrees")
  .sw_check_count(sample_size, "sample_size")
  .sw_check_count(ndim, "ndim")
  if (sample_size > nrow(x)) {
    stop("`sample_size` is larger than the number of files.", call. = FALSE)
  }
  if (ndim > ncol(x)) {
    stop("`ndim` is larger than the number of features.", call. = FALSE)
  }
  if (
    !is.numeric(cutoff) || length(cutoff) != 1 || cutoff <= 0 || cutoff >= 1
  ) {
    stop("`cutoff` must be one number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1 || seed != round(seed)) {
    stop("`seed` must be one whole number.", call. = FALSE)
  }

  forest <- isotree::isolation.forest(
    x,
    ntrees = ntrees,
    sample_size = sample_size,
    ndim = ndim,
    seed = seed,
    nthreads = 1
  )
  score <- as.numeric(stats::predict(forest, x))
  out <- tibble::tibble(
    file = rownames(x),
    score = score,
    flag = score >= cutoff
  )
  attr(out, "params") <- list(
    kind = "outlier",
    ntrees = ntrees,
    sample_size = sample_size,
    ndim = ndim,
    cutoff = cutoff,
    seed = seed
  )
  class(out) <- c("sw_anomaly_flag", class(out))
  out
}

#' Flag novel files with a Gaussian mixture fitted on reference files
#'
#' @description
#' Fit a Gaussian mixture on the features of the reference files, find the
#' `quantile` of the reference log-density, and flag every other file whose
#' log-density is at or below it. A novelty is a file unlike the reference
#' set, even when many such files exist.
#'
#' @details
#' **Fewer files than features stops the run.** With as many features as
#' reference files, a mixture fit is not stable: on synthetic data with 60 files
#' and 252 features `mclust` returned 9 components. The function stops and
#' names both counts. It never reduces the features for you, because a
#' reduction is a model-input transform that you approve.
#'
#' **Model choice.** `mclust::Mclust()` picks the number of components and the
#' covariance form by BIC. The chosen model is stored in attribute `model`.
#'
#' **No defaults.** `reference` and `quantile` are required. CytoScan uses a
#' 5% quantile; a higher value flags more files and raises false positives.
#'
#' @param features An `sw_anomaly_features` tibble, or a data frame with a
#'   `file` column and numeric feature columns.
#' @param reference Character vector of file ids that define "normal".
#' @param quantile Number in (0, 1). Test files at or below this quantile of
#'   the reference log-density are flagged.
#'
#' @return A tibble of class `sw_anomaly_flag` with `file`, `role`
#'   (`reference` or `test`), `log_density` and `flag`. `flag` is `NA` for
#'   reference files, since their threshold comes from themselves.
#'   Attributes `params` and `model` hold the arguments and the chosen
#'   mixture.
#'
#' @seealso [sw_anomaly_features()], [sw_anomaly_flag_outlier()]
#'
#' @examples
#' \dontrun{
#' sw_anomaly_flag_novelty(feats, reference = ref_files, quantile = 0.05)
#' }
#'
#' @export
sw_anomaly_flag_novelty <- function(features, reference, quantile) {
  .sw_stop_missing(c(
    features = missing(features),
    reference = missing(reference),
    quantile = missing(quantile)
  ))
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for novelty flags.", call. = FALSE)
  }
  x <- .sw_flag_matrix(features)
  if (!is.character(reference) || !all(reference %in% rownames(x))) {
    stop("`reference` must name files present in `features`.", call. = FALSE)
  }
  if (anyDuplicated(reference)) {
    stop("`reference` must not repeat a file.", call. = FALSE)
  }
  if (
    !is.numeric(quantile) ||
      length(quantile) != 1 ||
      quantile <= 0 ||
      quantile >= 1
  ) {
    stop("`quantile` must be one number in (0, 1).", call. = FALSE)
  }
  if (length(reference) <= ncol(x)) {
    stop(
      "Need more reference files than features: ",
      length(reference),
      " reference files, ",
      ncol(x),
      " features.",
      call. = FALSE
    )
  }

  ref_x <- x[reference, , drop = FALSE]
  # Mclust() evaluates mclustBIC() in the caller's frame, so bind it here
  # instead of attaching the package.
  mclustBIC <- mclust::mclustBIC
  fit <- mclust::Mclust(ref_x, verbose = FALSE)
  if (is.null(fit$parameters)) {
    stop("The mixture fit failed; check for constant features.", call. = FALSE)
  }
  logdens <- function(m) {
    log(mclust::dens(
      data = m,
      modelName = fit$modelName,
      parameters = fit$parameters
    ))
  }
  ref_ld <- logdens(ref_x)
  cutoff <- stats::quantile(ref_ld, quantile, names = FALSE)
  all_ld <- logdens(x)
  is_ref <- rownames(x) %in% reference
  out <- tibble::tibble(
    file = rownames(x),
    role = ifelse(is_ref, "reference", "test"),
    log_density = all_ld,
    flag = ifelse(is_ref, NA, all_ld <= cutoff)
  )
  attr(out, "params") <- list(
    kind = "novelty",
    reference = reference,
    quantile = quantile,
    log_density_cutoff = cutoff
  )
  attr(out, "model") <- list(model = fit$modelName, G = fit$G)
  class(out) <- c("sw_anomaly_flag", class(out))
  out
}

# Numeric feature matrix with file ids as row names.
.sw_flag_matrix <- function(features) {
  if (!is.data.frame(features) || !"file" %in% names(features)) {
    stop("`features` must be a data frame with a `file` column.", call. = FALSE)
  }
  if (anyDuplicated(features$file)) {
    stop("`features$file` must be unique.", call. = FALSE)
  }
  num <- features[setdiff(names(features), "file")]
  if (ncol(num) == 0 || !all(vapply(num, is.numeric, logical(1)))) {
    stop("All feature columns must be numeric.", call. = FALSE)
  }
  x <- as.matrix(as.data.frame(num))
  if (anyNA(x)) {
    stop("Features must not contain NA.", call. = FALSE)
  }
  rownames(x) <- as.character(features$file)
  x
}
