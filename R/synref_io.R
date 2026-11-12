# Synthetic references: manifest rows, save and load -------------------------------

# Stop unless `x` is an sw_synref with the current schema and an intact digest.
.sw_validate_synref <- function(x) {
  if (!inherits(x, "sw_synref")) {
    stop("'synref' must be an sw_synref.", call. = FALSE)
  }
  need <- c(
    "schema_version",
    "events",
    "provenance",
    "recipe",
    "realised",
    "scope",
    "floor",
    "certainty",
    "digest"
  )
  if (!all(need %in% names(x))) {
    stop(
      "Synthetic reference is missing field(s): ",
      paste(setdiff(need, names(x)), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!identical(x[["schema_version"]], .sw_synref_schema)) {
    stop(
      sprintf(
        "Unsupported schema_version '%s'; this package reads '%s'.",
        as.character(x[["schema_version"]]),
        .sw_synref_schema
      ),
      call. = FALSE
    )
  }
  if (!identical(x[["digest"]], .sw_model_digest(x))) {
    stop(
      "Digest mismatch: the synthetic reference changed after it was built. ",
      "Rebuild it, or load the original.",
      call. = FALSE
    )
  }
  invisible(x)
}

#' Manifest rows for a synthetic reference
#'
#' Returns one row per synthetic reference with `reference_role =
#' "synthetic_reference"` and a `source_batch` note. A synthetic reference is
#' never a measured reference and never bridges batches in the design audit.
#' The rows are a plain tibble. Bind them to the study manifest table, then
#' call [sw_sample_manifest()] with `reference_role`, `reference_id` and
#' `source_batch` mapped. [sw_design_audit()] lists them under
#' `synthetic_references` and does not let them bridge batches.
#'
#' @param synref An `sw_synref`.
#' @param manifest Data frame with `sample_id` and `batch`. Its `study` and
#'   `run` columns, when present, are copied from the source sample.
#'
#' @return A tibble with `sample_id` (the synref id), `batch`, `reference_role`,
#'   `reference_id`, `source_batch` and any `study` or `run` of the source.
#'
#' @examples
#' set.seed(1)
#' d <- list(a = cbind(M1 = rnorm(100)))
#' man <- data.frame(sample_id = "a", batch = "b1")
#' ref <- sw_synref_downsample(
#'   d, man, groups = list(r = "a"), strata = "file", labels = NULL,
#'   n_events = 20, n_draws = 1, seed = 1
#' )
#' sw_synref_as_manifest(ref, man)
#'
#' @export
sw_synref_as_manifest <- function(synref, manifest) {
  .sw_stop_missing(c(synref = missing(synref), manifest = missing(manifest)))
  .sw_validate_synref(synref)
  if (!identical(synref[["scope"]], "within_batch")) {
    stop(
      "Only within-batch synthetic references become manifest rows.",
      call. = FALSE
    )
  }
  prov <- synref[["provenance"]]
  .sw_synref_check_manifest(manifest, unique(prov[["source_sample"]]))
  rows <- lapply(unique(prov[["synref_id"]]), function(id) {
    p <- prov[prov[["synref_id"]] == id, , drop = FALSE]
    batch <- unique(as.character(p[["batch"]]))
    if (length(batch) != 1) {
      stop(
        sprintf("Synthetic reference '%s' spans batches.", id),
        call. = FALSE
      )
    }
    m <- manifest[
      match(p[["source_sample"]][1], as.character(manifest[["sample_id"]])),
      ,
      drop = FALSE
    ]
    out <- tibble::tibble(
      sample_id = id,
      batch = batch,
      reference_role = "synthetic_reference",
      reference_id = id,
      source_batch = batch
    )
    for (col in intersect(c("study", "run"), names(manifest))) {
      out[[col]] <- m[[col]]
    }
    out
  })
  dplyr::bind_rows(rows)
}

#' Save and load a synthetic reference
#'
#' `sw_synref_save()` writes the object, a readable provenance table and the
#' digest into a new directory. It never overwrites: a directory that exists
#' stops the call. `sw_synref_load()` stops when the digest does not match.
#'
#' @param synref An `sw_synref`.
#' @param dir Directory to create. It must not exist. For `sw_synref_load()`,
#'   the directory written by `sw_synref_save()`.
#'
#' @return `sw_synref_save()` returns `dir` invisibly. `sw_synref_load()`
#'   returns the `sw_synref`.
#'
#' @examples
#' set.seed(1)
#' d <- list(a = cbind(M1 = rnorm(100)))
#' man <- data.frame(sample_id = "a", batch = "b1")
#' ref <- sw_synref_downsample(
#'   d, man, groups = list(r = "a"), strata = "file", labels = NULL,
#'   n_events = 20, n_draws = 1, seed = 1
#' )
#' dir <- file.path(tempdir(), "synref_example")
#' sw_synref_save(ref, dir)
#' identical(sw_synref_load(dir)$digest, ref$digest)
#'
#' @export
sw_synref_save <- function(synref, dir) {
  .sw_stop_missing(c(synref = missing(synref), dir = missing(dir)))
  .sw_validate_synref(synref)
  if (!is.character(dir) || length(dir) != 1 || is.na(dir)) {
    stop("'dir' must be one path.", call. = FALSE)
  }
  if (file.exists(dir)) {
    stop(
      sprintf(
        "'%s' already exists; a synthetic reference is never overwritten.",
        dir
      ),
      call. = FALSE
    )
  }
  dir.create(dir, recursive = TRUE)
  saveRDS(synref, file.path(dir, "synref.rds"))
  utils::write.csv(
    synref[["provenance"]],
    file.path(dir, "provenance.csv"),
    row.names = FALSE
  )
  writeLines(synref[["digest"]], file.path(dir, "digest.txt"))
  invisible(dir)
}

#' @rdname sw_synref_save
#' @export
sw_synref_load <- function(dir) {
  .sw_stop_missing(c(dir = missing(dir)))
  path <- file.path(dir, "synref.rds")
  if (!file.exists(path)) {
    stop(sprintf("No synref.rds in '%s'.", dir), call. = FALSE)
  }
  x <- readRDS(path)
  .sw_validate_synref(x)
  stamp <- file.path(dir, "digest.txt")
  if (
    file.exists(stamp) &&
      !identical(trimws(readLines(stamp, warn = FALSE)[1]), x[["digest"]])
  ) {
    stop(
      "Digest mismatch: digest.txt does not match synref.rds.",
      call. = FALSE
    )
  }
  x
}
