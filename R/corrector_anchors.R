# Internal anchors for the heb engine (plan A1, A2, A7): per sample and marker, a two-mode split.

#' Internal anchor modes per sample and marker
#'
#' Splits each sample's cells into two modes per marker with one-dimensional
#' k-means (k = 2) and returns the median of each mode. A marker is eligible
#' when its median separation over samples, (high centre - low centre) /
#' pooled within-mode SD, is at least \code{min_separation}. The low mode of
#' every eligible marker is used as an anchor; the high mode is used only for
#' \code{positive_markers} (a constitutively bright, donor-stable population,
#' such as CD3 on T cells; Glehr et al. 2026, eBioMedicine). Pass the result
#' to \code{\link{sw_corrector_fit}} as \code{anchors} with \code{engine = "heb"}.
#'
#' Only give cells that may calibrate: never held-out reference samples, and
#' only samples present at fit time, so the frozen model stays reproducible.
#'
#' @param x Matrix or data frame of transformed values, cells by markers.
#' @param sample_id,batch Per-cell sample and batch labels.
#' @param markers Marker columns to split.
#' @param positive_markers Markers whose high mode is also an anchor.
#' @param min_separation Required; minimum median separation in pooled SDs.
#' @param seed Required; seed for k-means.
#' @return A list: \code{modes} (sample_id, batch, marker, mode "low"/"high",
#'   value, var_mean, n_cells, use) and \code{markers} (marker, median
#'   separation, eligible).
#' @export
sw_corrector_anchors <- function(x, sample_id, batch, markers, positive_markers, min_separation, seed) {
  if (missing(min_separation) || !is.numeric(min_separation) || length(min_separation) != 1) {
    stop("'min_separation' is required (no default).", call. = FALSE)
  }
  if (missing(seed) || !is.numeric(seed) || length(seed) != 1) {
    stop("'seed' is required (no default).", call. = FALSE)
  }
  if (missing(positive_markers) || !all(positive_markers %in% markers)) {
    stop("'positive_markers' must be a subset of 'markers' (no default).", call. = FALSE)
  }
  xm <- as.matrix(x[, markers, drop = FALSE])
  stopifnot("x, sample_id and batch have the same rows" = nrow(xm) == length(sample_id) && nrow(xm) == length(batch))
  samples <- unique(sample_id)
  rows <- lapply(samples, function(s) {
    ii <- which(sample_id == s)
    lapply(markers, function(m) {
      v <- xm[ii, m]
      if (length(unique(v)) < 2) {
        return(NULL)
      }
      set.seed(seed)
      km <- stats::kmeans(v, centers = range(v))
      hi <- which.max(km$centers)
      grp <- ifelse(km$cluster == hi, "high", "low")
      pooled_sd <- sqrt(sum(tapply(v, grp, function(g) sum((g - mean(g))^2))) / (length(v) - 2))
      tibble::tibble(
        sample_id = s, batch = batch[ii[1]], marker = m, mode = c("low", "high"),
        value = c(stats::median(v[grp == "low"]), stats::median(v[grp == "high"])),
        var_mean = c(stats::var(v[grp == "low"]) / sum(grp == "low"), stats::var(v[grp == "high"]) / sum(grp == "high")),
        n_cells = c(sum(grp == "low"), sum(grp == "high")),
        separation = diff(sort(as.vector(km$centers))) / pooled_sd
      )
    })
  })
  modes <- dplyr::bind_rows(unlist(rows, recursive = FALSE))
  sep <- dplyr::summarise(modes, separation = stats::median(separation), .by = marker)
  sep$eligible <- sep$separation >= min_separation
  modes$use <- modes$marker %in% sep$marker[sep$eligible] &
    (modes$mode == "low" | modes$marker %in% positive_markers) &
    modes$n_cells >= 2
  list(modes = dplyr::select(modes, -separation), markers = sep)
}
