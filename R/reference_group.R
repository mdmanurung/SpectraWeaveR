# R/reference_group.R
# Shared helpers for the `reference_group` argument: robust z-scores taken within
# groups of files (for example the same donor) instead of over all files.

# NULL, or a named character vector with one non-missing label per file.
# Returns NULL, or the labels in the order of `files`.
.sw_check_reference_group <- function(reference_group, files) {
  if (is.null(reference_group)) {
    return(NULL)
  }
  if (
    !is.character(reference_group) ||
      is.null(names(reference_group)) ||
      anyDuplicated(names(reference_group)) ||
      !all(files %in% names(reference_group))
  ) {
    stop(
      "`reference_group` must be NULL or a named character vector with one ",
      "label for every file in `data` (names = file ids).",
      call. = FALSE
    )
  }
  g <- unname(reference_group[files])
  if (anyNA(g)) {
    stop("`reference_group` has missing labels.", call. = FALSE)
  }
  g
}

# Robust z (x - median) / (1.4826 * raw MAD). With `group` (a vector the length
# of x), the median and MAD come from the members of each group only. A group
# with fewer than 3 non-missing values, or with MAD 0, gives NA.
.sw_robust_z <- function(x, group = NULL) {
  one <- function(v) {
    s <- 1.4826 * stats::mad(v, constant = 1, na.rm = TRUE)
    if (sum(!is.na(v)) < 3 || is.na(s) || s == 0) {
      return(rep(NA_real_, length(v)))
    }
    (v - stats::median(v, na.rm = TRUE)) / s
  }
  if (is.null(group)) {
    return(one(x))
  }
  out <- rep(NA_real_, length(x))
  for (g in unique(group)) {
    i <- which(group == g)
    out[i] <- one(x[i])
  }
  out
}
