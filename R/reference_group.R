# R/reference_group.R
# Shared helpers for the `reference_group` argument: robust z-scores taken within
# groups of files (for example the same donor) instead of over all files, or on
# the residuals of a two-way median polish (for example donor and batch).

# NULL; a named character vector with one label per file (one-way, returns the
# labels in the order of `files`); or a data frame with a `file` column and two
# more character columns (two-way, returns list(a, b) in the order of `files`
# with class "sw_two_way").
.sw_check_reference_group <- function(reference_group, files) {
  if (is.null(reference_group)) {
    return(NULL)
  }
  if (is.data.frame(reference_group)) {
    rg <- reference_group
    bad <- !"file" %in% names(rg) ||
      ncol(rg) != 3 ||
      !all(vapply(rg[setdiff(names(rg), "file")], is.character, logical(1))) ||
      !is.character(rg$file) ||
      anyDuplicated(rg$file) ||
      !all(files %in% rg$file)
    if (bad) {
      stop(
        "A data-frame `reference_group` needs a character `file` column that ",
        "covers every file in `data` once, and exactly two more character ",
        "columns (for example donor and batch).",
        call. = FALSE
      )
    }
    fac <- setdiff(names(rg), "file")
    i <- match(files, rg$file)
    a <- rg[[fac[1]]][i]
    b <- rg[[fac[2]]][i]
    if (anyNA(a) || anyNA(b)) {
      stop("`reference_group` has missing labels.", call. = FALSE)
    }
    return(structure(list(a = a, b = b), class = "sw_two_way"))
  }
  if (
    !is.character(reference_group) ||
      is.null(names(reference_group)) ||
      anyDuplicated(names(reference_group)) ||
      !all(files %in% names(reference_group))
  ) {
    stop(
      "`reference_group` must be NULL, a named character vector with one ",
      "label for every file in `data` (names = file ids), or a data frame ",
      "with a `file` column and two character columns.",
      call. = FALSE
    )
  }
  g <- unname(reference_group[files])
  if (anyNA(g)) {
    stop("`reference_group` has missing labels.", call. = FALSE)
  }
  g
}

# One-way labels only: the vector, or NULL when the design is NULL or two-way.
.sw_one_way <- function(grp) if (is.character(grp)) grp else NULL

# Residual of each value after a two-way median polish on factors a and b (one
# value per file; files sharing a cell enter through the cell median). Rows or
# columns with no finite value are left out and their files get NA.
.sw_medpolish_resid <- function(x, a, b) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) < 3) {
    return(out)
  }
  ra <- sort(unique(a[ok]))
  cb <- sort(unique(b[ok]))
  tab <- matrix(NA_real_, length(ra), length(cb), dimnames = list(ra, cb))
  ia <- match(a[ok], ra)
  ib <- match(b[ok], cb)
  cell_id <- paste(ia, ib)
  for (k in unique(cell_id)) {
    sel <- cell_id == k
    tab[ia[sel][1], ib[sel][1]] <- stats::median(x[ok][sel])
  }
  mp <- stats::medpolish(tab, maxiter = 10, trace.iter = FALSE, na.rm = TRUE)
  out[ok] <- x[ok] - mp$overall - mp$row[ra[ia]] - mp$col[cb[ib]]
  out
}

# Robust z (x - median) / (1.4826 * raw MAD). With one-way `group` (a vector the
# length of x), the median and MAD come from the members of each group only. With
# a two-way design, the z is taken over all files on the residuals of a median
# polish. Fewer than 3 finite values, or MAD 0, gives NA.
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
  if (inherits(group, "sw_two_way")) {
    return(one(.sw_medpolish_resid(x, group$a, group$b)))
  }
  out <- rep(NA_real_, length(x))
  for (g in unique(group)) {
    i <- which(group == g)
    out[i] <- one(x[i])
  }
  out
}
