# R/signal_reference.R
# Frozen reference for the per-file signal audit: per-channel centre and scale
# of the file-level median and SD, fitted once on clean files and reused to
# score later files. Pattern and digest follow sw_frozen_model (R/corrector.R).

.sw_signal_reference_schema <- "sw_signal_reference/1"

#' Fit a frozen signal reference from clean files
#'
#' @description
#' Summarise each reference file by its median and SD per channel. For each
#' channel, the reference keeps the centre (median over files) and the scale
#' (`1.4826 * MAD` over files) of those file-level medians, and of those
#' file-level SDs. Pass the result to the `reference` argument of
#' [sw_fcs_signal_audit()]: its `median_z` and `sd_z` then measure a file
#' against these clean files instead of against the files scored with it. A
#' fault that hits every scored file then still shows.
#'
#' The reference is only as clean as the files you give it. Which files count
#' as clean is your choice.
#'
#' @param data A named list. Each element is a numeric matrix or data frame
#'   with one column per channel, as for [sw_fcs_signal_audit()]. Names are the
#'   file ids. At least 3 files.
#' @param channels Character vector of channel names to keep in the reference.
#'
#' @return An object of class `sw_signal_reference`: a list with
#'   `schema_version`, `channels` (a tibble with `channel`, `median_centre`,
#'   `median_scale`, `sd_centre`, `sd_scale`), `files` (reference file ids),
#'   `n_files`, `created` (text) and `digest` (MD5 of every other field).
#'
#' @seealso [sw_signal_reference_save()], [sw_fcs_signal_audit()]
#'
#' @examples
#' set.seed(1)
#' files <- lapply(1:4, function(i) {
#'   matrix(stats::rnorm(400), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' })
#' names(files) <- paste0("f", 1:4)
#' ref <- sw_signal_reference_fit(files, channels = c("A", "B"))
#' ref$channels
#'
#' @export
sw_signal_reference_fit <- function(data, channels) {
  .sw_stop_missing(c(data = missing(data), channels = missing(channels)))
  if (!is.list(data) || is.null(names(data)) || anyDuplicated(names(data))) {
    stop(
      "`data` must be a named list of files with unique names.",
      call. = FALSE
    )
  }
  if (length(data) < 3) {
    stop("A reference needs at least 3 files.", call. = FALSE)
  }
  if (
    !is.character(channels) || length(channels) < 1 || anyDuplicated(channels)
  ) {
    stop("`channels` must be unique channel names.", call. = FALSE)
  }
  files <- sort(names(data))
  med_mat <- sd_mat <- matrix(
    NA_real_,
    length(files),
    length(channels),
    dimnames = list(files, channels)
  )
  for (f in files) {
    m <- as.matrix(data[[f]])
    if (!all(channels %in% colnames(m))) {
      stop(
        "File '",
        f,
        "' lacks column(s): ",
        paste(setdiff(channels, colnames(m)), collapse = ", "),
        call. = FALSE
      )
    }
    for (ch in channels) {
      v <- m[, ch]
      v <- v[is.finite(v)]
      if (length(v) > 1) {
        med_mat[f, ch] <- stats::median(v)
        sd_mat[f, ch] <- stats::sd(v)
      }
    }
  }
  centre <- function(mat) apply(mat, 2, stats::median, na.rm = TRUE)
  scale <- function(mat) {
    1.4826 * apply(mat, 2, stats::mad, constant = 1, na.rm = TRUE)
  }
  ref <- structure(
    list(
      schema_version = .sw_signal_reference_schema,
      channels = tibble::tibble(
        channel = channels,
        median_centre = unname(centre(med_mat)),
        median_scale = unname(scale(med_mat)),
        sd_centre = unname(centre(sd_mat)),
        sd_scale = unname(scale(sd_mat))
      ),
      files = files,
      n_files = length(files),
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      digest = NA_character_
    ),
    class = "sw_signal_reference"
  )
  ref$digest <- .sw_model_digest(ref)
  ref
}

#' @export
print.sw_signal_reference <- function(x, ...) {
  cat(
    "<sw_signal_reference>",
    x$n_files,
    "files,",
    nrow(x$channels),
    "channels | created",
    x$created,
    "| digest:",
    x$digest,
    "\n"
  )
  invisible(x)
}

#' Save and load a frozen signal reference
#'
#' Stores the reference as an RDS file and checks its class, schema version
#' and digest on load. The digest is an MD5 of a canonical text form of every
#' field (doubles written in hex). Load fails when the file was changed after
#' saving. Only load files you trust: `readRDS()` can run code stored in a
#' crafted file.
#'
#' @param ref An `sw_signal_reference` from [sw_signal_reference_fit()].
#' @param path File path.
#'
#' @return `sw_signal_reference_save()` returns `path` invisibly;
#'   `sw_signal_reference_load()` returns the `sw_signal_reference`.
#'
#' @examples
#' set.seed(1)
#' files <- lapply(1:4, function(i) {
#'   matrix(stats::rnorm(400), ncol = 2, dimnames = list(NULL, c("A", "B")))
#' })
#' names(files) <- paste0("f", 1:4)
#' ref <- sw_signal_reference_fit(files, channels = c("A", "B"))
#' path <- tempfile(fileext = ".rds")
#' sw_signal_reference_save(ref, path)
#' identical(sw_signal_reference_load(path)$digest, ref$digest)
#'
#' @export
sw_signal_reference_save <- function(ref, path) {
  .sw_stop_missing(c(ref = missing(ref), path = missing(path)))
  .sw_validate_signal_reference(ref)
  saveRDS(ref, path)
  invisible(path)
}

#' @rdname sw_signal_reference_save
#' @export
sw_signal_reference_load <- function(path) {
  .sw_stop_missing(c(path = missing(path)))
  ref <- readRDS(path)
  .sw_validate_signal_reference(ref)
  ref
}

# Stop unless `ref` has the current schema and an intact digest.
.sw_validate_signal_reference <- function(ref) {
  if (!inherits(ref, "sw_signal_reference")) {
    stop("'ref' must be an sw_signal_reference.", call. = FALSE)
  }
  need <- c(
    "schema_version",
    "channels",
    "files",
    "n_files",
    "created",
    "digest"
  )
  if (!all(need %in% names(ref))) {
    stop(
      "Signal reference is missing field(s): ",
      paste(setdiff(need, names(ref)), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!identical(ref$schema_version, .sw_signal_reference_schema)) {
    stop(
      sprintf(
        "Unsupported schema_version '%s'; this package reads '%s'.",
        as.character(ref$schema_version),
        .sw_signal_reference_schema
      ),
      call. = FALSE
    )
  }
  if (!identical(ref$digest, .sw_model_digest(ref))) {
    stop(
      "Digest mismatch: the reference changed after it was fitted. ",
      "Refit, or load the original file.",
      call. = FALSE
    )
  }
  invisible(ref)
}
