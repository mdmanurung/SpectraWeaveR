#' Build a validated sample manifest
#'
#' One row per sample, with standard column names. The manifest is the single
#' sample table shared by the design audit (\code{\link{sw_design_audit}}),
#' the correctors (\code{sw_corrector_*}) and study intake. It replaces the
#' \code{sample}/\code{batch}/\code{condition} table that
#' \code{\link{sw_correct_prepare}} reads.
#'
#' Standard columns in the result: \code{sample_id}, \code{batch}, and any of
#' \code{study}, \code{run}, \code{instrument}, \code{donor}, \code{aliquot},
#' \code{condition}, \code{reference_role}, \code{reference_id},
#' \code{source_batch} that you map.
#' Extra columns of \code{df} are kept after them.
#'
#' Checks, all of which stop with a message:
#' \itemize{
#'   \item \code{sample_id} is unique and not missing.
#'   \item \code{batch} is not missing.
#'   \item Each batch belongs to one \code{study} (batches are nested in
#'     studies), and each \code{run} belongs to one batch.
#'   \item \code{reference_role} and \code{reference_id} come together, and
#'     the role is one of \code{none}, \code{technical_replicate},
#'     \code{anchor}, \code{biological_pseudo_control},
#'     \code{synthetic_reference}.
#'   \item A sample with a measured reference role has a \code{reference_id}.
#'   \item A sample with role \code{synthetic_reference} has a
#'     \code{source_batch}: the batch its cells were drawn from. It must be
#'     one of the batches in \code{batch}; it may equal the sample's own
#'     batch. Like
#'     \code{biological_pseudo_control}, it is never a measured reference and
#'     never bridges batches in \code{\link{sw_design_audit}}.
#' }
#'
#' @param df A data frame with one row per sample.
#' @param sample_id,batch Names of the columns in \code{df} that hold the
#'   sample id and the batch. Required, no default, because the source
#'   column names differ between projects.
#' @param study,run,instrument,donor,aliquot,condition,reference_role,reference_id,source_batch
#'   Optional names of columns in \code{df}. \code{NULL} means the column is
#'   not available.
#'
#' @return A tibble of class \code{sw_manifest}.
#'
#' @seealso \code{\link{sw_manifest_from_sample_meta}}
#'
#' @examples
#' df <- data.frame(
#'   id = c("s1", "s2", "s3", "s4"),
#'   b = c("b1", "b1", "b2", "b2"),
#'   cohort = c("A", "A", "B", "B"),
#'   donor_id = c("d1", "d2", "d1", "d3")
#' )
#' sw_sample_manifest(df, sample_id = "id", batch = "b", study = "cohort",
#'                    donor = "donor_id")
#'
#' @export
sw_sample_manifest <- function(df,
                               sample_id,
                               batch,
                               study = NULL,
                               run = NULL,
                               instrument = NULL,
                               donor = NULL,
                               aliquot = NULL,
                               condition = NULL,
                               reference_role = NULL,
                               reference_id = NULL,
                               source_batch = NULL) {
  if (missing(sample_id) || missing(batch)) {
    stop("'sample_id' and 'batch' are required: name the columns of 'df'.",
         call. = FALSE)
  }
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame.", call. = FALSE)
  }
  map <- list(
    sample_id = sample_id, batch = batch, study = study, run = run,
    instrument = instrument, donor = donor, aliquot = aliquot,
    condition = condition, reference_role = reference_role,
    reference_id = reference_id, source_batch = source_batch
  )
  map <- Filter(Negate(is.null), map)
  bad_arg <- names(map)[!vapply(
    map, function(v) is.character(v) && length(v) == 1 && !is.na(v), logical(1)
  )]
  if (length(bad_arg) > 0) {
    stop("Column name for ", paste(bad_arg, collapse = ", "),
         " must be a single string.", call. = FALSE)
  }
  absent <- unlist(map)[!unlist(map) %in% names(df)]
  if (length(absent) > 0) {
    stop("Column(s) not found in 'df': ", paste(absent, collapse = ", "),
         call. = FALSE)
  }

  # Standard columns first, then the rest of df under their own names.
  std <- tibble::as_tibble(lapply(map, function(col) {
    v <- df[[col]]
    if (is.factor(v)) as.character(v) else v
  }))
  rest <- setdiff(names(df), unlist(map))
  clash <- intersect(rest, c(
    "sample_id", "batch", "study", "run", "instrument", "donor", "aliquot",
    "condition", "reference_role", "reference_id", "source_batch"
  ))
  if (length(clash) > 0) {
    stop("Unmapped column(s) of 'df' share a standard name: ",
         paste(clash, collapse = ", "), ". Map or rename them.",
         call. = FALSE)
  }
  out <- tibble::as_tibble(c(as.list(std), as.list(df[rest])))

  if (anyNA(out$sample_id)) {
    stop("'sample_id' has missing values.", call. = FALSE)
  }
  if (anyDuplicated(out$sample_id) > 0) {
    stop("'sample_id' must be unique; repeated: ",
         paste(unique(out$sample_id[duplicated(out$sample_id)]), collapse = ", "),
         call. = FALSE)
  }
  if (anyNA(out$batch)) {
    stop("'batch' has missing values.", call. = FALSE)
  }
  .sw_manifest_check_nesting(out, child = "batch", parent = "study")
  .sw_manifest_check_nesting(out, child = "run", parent = "batch")

  has_ref <- c("reference_role", "reference_id") %in% names(out)
  if (any(has_ref) && !all(has_ref)) {
    stop("'reference_role' and 'reference_id' must be given together.",
         call. = FALSE)
  }
  if (all(has_ref)) {
    roles <- c("none", "technical_replicate", "anchor",
               "biological_pseudo_control", "synthetic_reference")
    bad <- setdiff(unique(stats::na.omit(as.character(out$reference_role))),
                   roles)
    if (length(bad) > 0) {
      stop("Unknown reference_role value(s): ", paste(bad, collapse = ", "),
           ". Use ", paste(roles, collapse = ", "), ".", call. = FALSE)
    }
    measured <- out$reference_role %in% c("technical_replicate", "anchor")
    if (any(measured & is.na(out$reference_id))) {
      stop("Samples with a measured reference role need a 'reference_id'.",
           call. = FALSE)
    }
    synthetic <- out$reference_role %in% "synthetic_reference"
    if (any(synthetic) && (!"source_batch" %in% names(out) ||
                           any(is.na(out$source_batch[synthetic])))) {
      stop("Samples with role 'synthetic_reference' need a 'source_batch'.",
           call. = FALSE)
    }
    unknown_src <- setdiff(
      unique(as.character(out$source_batch[synthetic])),
      unique(as.character(out$batch))
    )
    if (length(unknown_src) > 0) {
      stop("'source_batch' value(s) not found in 'batch': ",
           paste(unknown_src, collapse = ", "), ".", call. = FALSE)
    }
  }

  structure(out, class = c("sw_manifest", class(out)))
}

# Stop when a child level maps to more than one parent level.
.sw_manifest_check_nesting <- function(df, child, parent) {
  if (!all(c(child, parent) %in% names(df))) {
    return(invisible(NULL))
  }
  d <- df[!is.na(df[[child]]) & !is.na(df[[parent]]), ]
  n_parent <- tapply(as.character(d[[parent]]), as.character(d[[child]]),
                     function(v) length(unique(v)))
  multi <- names(n_parent)[n_parent > 1]
  if (length(multi) > 0) {
    stop(child, " level(s) ", paste(multi, collapse = ", "),
         " belong to more than one ", parent, ".", call. = FALSE)
  }
  invisible(NULL)
}

#' Convert the older sample table to a manifest
#'
#' \code{\link{sw_correct_prepare}} reads a table with \code{sample},
#' \code{batch} and optionally \code{condition}. This turns it into a
#' manifest with \code{sample_id}, so the same table can feed the design
#' audit and the correctors.
#'
#' @param sample_meta A data frame with columns \code{sample} and
#'   \code{batch}, and optionally \code{condition}.
#' @param ... Further column mappings passed to
#'   \code{\link{sw_sample_manifest}} (for example \code{donor = "donor"}).
#'
#' @return A tibble of class \code{sw_manifest}.
#'
#' @examples
#' meta <- data.frame(sample = c("a", "b"), batch = c("b1", "b2"),
#'                    condition = c("ctrl", "case"))
#' sw_manifest_from_sample_meta(meta, condition = "condition")
#'
#' @export
sw_manifest_from_sample_meta <- function(sample_meta, ...) {
  if (!is.data.frame(sample_meta) ||
      !all(c("sample", "batch") %in% names(sample_meta))) {
    stop("'sample_meta' must be a data frame with 'sample' and 'batch'.",
         call. = FALSE)
  }
  sw_sample_manifest(sample_meta, sample_id = "sample", batch = "batch", ...)
}

#' @export
print.sw_manifest <- function(x, ...) {
  cat("<sw_manifest>", nrow(x), "samples,",
      length(unique(x$batch)), "batches")
  if ("study" %in% names(x)) cat(",", length(unique(x$study)), "studies")
  cat("\n")
  NextMethod()
}
