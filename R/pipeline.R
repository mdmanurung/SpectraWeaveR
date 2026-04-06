#' @title End-to-End Pipeline Orchestrator
#'
#' @description
#' Runs the complete SpectraWeaveR pipeline: unmixed FCS loading, margin
#' removal, gating, signal QC, batch correction, and clustering.
#'
#' @name pipeline
#' @keywords internal
NULL

#' Run the Full SpectraWeaveR Pipeline
#'
#' Executes all five steps of the spectral flow cytometry analysis pipeline
#' in validated order:
#' \enumerate{
#'   \item Load unmixed FCS files
#'   \item Remove margin events
#'   \item Apply automated gating (optional)
#'   \item Signal quality control (PeacoQC)
#'   \item Batch correction (cyCombine)
#'   \item Clustering (FlowSOM)
#' }
#'
#' @param fcs_dir Character path to directory containing unmixed FCS files.
#' @param sample_meta A \code{data.frame} or \code{tibble} with columns
#'   \code{file} (matching FCS file names), \code{sample}, \code{batch},
#'   and optionally \code{condition}.
#' @param markers Character vector of marker/channel names to use throughout
#'   the pipeline.
#' @param lineage_markers Character vector of lineage markers for clustering
#'   (subset of \code{markers}).
#' @param gating_template Character path to a CSV gating template file, or
#'   \code{NULL} to skip gating (default: \code{NULL}).
#' @param gate_node Character string specifying which gate node to extract
#'   (default: \code{"lymphocytes"}).
#' @param output_dir Character path for output files (plots, reports).
#'   Default: \code{"SpectraWeaveR_output"}.
#' @param cofactor Numeric; arcsinh cofactor (default: 6000).
#' @param n_metaclusters Integer; number of metaclusters for FlowSOM
#'   (default: 20).
#' @param seed Integer; random seed for reproducibility (default: 42).
#' @param ... Additional arguments.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{flowset}}{The loaded \code{flowSet}}
#'     \item{\code{gating_set}}{The \code{GatingSet} (if gating was applied)}
#'     \item{\code{qc_results}}{QC summary from PeacoQC}
#'     \item{\code{uncorrected}}{Pre-correction tibble}
#'     \item{\code{corrected}}{Post-correction tibble}
#'     \item{\code{correction_eval}}{Batch correction evaluation metrics}
#'     \item{\code{fsom}}{FlowSOM result object}
#'     \item{\code{cluster_assignments}}{Per-cell metacluster assignments}
#'     \item{\code{cluster_mfis}}{MFI table per metacluster}
#'   }
#'
#' @export
run_pipeline <- function(fcs_dir,
                         sample_meta,
                         markers,
                         lineage_markers,
                         gating_template = NULL,
                         gate_node = "lymphocytes",
                         output_dir = "SpectraWeaveR_output",
                         cofactor = 6000,
                         n_metaclusters = 20,
                         seed = 42,
                         ...) {
  # --- Input validation ---
  if (!is.character(fcs_dir) || length(fcs_dir) != 1) {
    stop("'fcs_dir' must be a single directory path.", call. = FALSE)
  }

  if (!dir.exists(fcs_dir)) {
    stop("FCS directory does not exist: ", fcs_dir, call. = FALSE)
  }

  if (!is.data.frame(sample_meta)) {
    stop("'sample_meta' must be a data.frame or tibble.", call. = FALSE)
  }

  required_meta_cols <- c("file", "sample", "batch")
  missing_meta <- setdiff(required_meta_cols, names(sample_meta))
  if (length(missing_meta) > 0) {
    stop("'sample_meta' is missing required column(s): ",
         paste(missing_meta, collapse = ", "), call. = FALSE)
  }

  if (!is.character(markers) || length(markers) == 0) {
    stop("'markers' must be a non-empty character vector.", call. = FALSE)
  }

  if (!is.character(lineage_markers) || length(lineage_markers) == 0) {
    stop("'lineage_markers' must be a non-empty character vector.",
         call. = FALSE)
  }

  non_lineage <- setdiff(lineage_markers, markers)
  if (length(non_lineage) > 0) {
    stop("'lineage_markers' must be a subset of 'markers'. Not found: ",
         paste(non_lineage, collapse = ", "), call. = FALSE)
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  results <- list()

  # --- Step 1: Load unmixed FCS files ---
  message("\n=== Step 1/6: Loading unmixed FCS files ===")
  fs <- sw_load_unmixed(fcs_dir)
  results$flowset <- fs

  # --- Step 2: Remove margin events ---
  message("\n=== Step 2/6: Removing margin events ===")
  ff_list <- list()
  for (i in seq_along(fs)) {
    sn <- flowCore::sampleNames(fs)[i]
    message("  Processing: ", sn)
    ff_list[[sn]] <- sw_remove_margins(fs[[i]])
  }

  # --- Step 3: Gating (optional) ---
  if (!is.null(gating_template)) {
    message("\n=== Step 3/6: Automated gating ===")
    fcs_files <- file.path(fcs_dir,
                           list.files(fcs_dir, pattern = "\\.fcs$",
                                      ignore.case = TRUE))
    gs <- sw_gate(fcs_files, gating_template)
    results$gating_set <- gs

    # Extract gated populations
    ff_list <- sw_extract_gated(gs, gate_node)
    message("Extracted gated populations from node: ", gate_node)
  } else {
    message("\n=== Step 3/6: Gating skipped (no template provided) ===")
    results$gating_set <- NULL
  }

  # --- Step 4: Signal QC ---
  message("\n=== Step 4/6: Signal quality control (PeacoQC) ===")
  qc_output_dir <- file.path(output_dir, "QC")
  qc_batch_result <- sw_signal_qc_batch(ff_list,
                                         output_dir = qc_output_dir)
  results$qc_results <- qc_batch_result

  # Check for high-removal samples
  qc_summary <- sw_qc_summary(qc_batch_result)
  message("QC summary: ", sum(qc_summary$flagged), " samples flagged")

  ff_list_clean <- qc_batch_result$cleaned

  # --- Step 5: Batch correction ---
  message("\n=== Step 5/6: Batch correction (cyCombine) ===")

  # Prepare the sample_meta for batch correction (use 'sample' column)
  correction_meta <- sample_meta[, intersect(names(sample_meta),
                                              c("sample", "batch",
                                                "condition")),
                                  drop = FALSE]

  uncorrected <- sw_prepare_for_correction(
    ff_list_clean,
    sample_meta = correction_meta,
    markers = markers,
    cofactor = cofactor
  )
  results$uncorrected <- uncorrected

  covar <- if ("condition" %in% names(sample_meta)) "condition" else NULL

  corrected <- sw_batch_correct(
    uncorrected,
    markers = markers,
    covar = covar,
    seed = seed
  )
  results$corrected <- corrected

  # Evaluate correction
  eval_result <- sw_evaluate_correction(uncorrected, corrected,
                                         markers = markers)
  results$correction_eval <- eval_result

  # --- Step 6: Clustering ---
  message("\n=== Step 6/6: Clustering (FlowSOM) ===")

  fsom <- sw_cluster(
    corrected,
    lineage_markers = lineage_markers,
    n_metaclusters = n_metaclusters,
    seed = seed
  )
  results$fsom <- fsom

  results$cluster_assignments <- sw_get_cluster_assignments(fsom)
  results$cluster_mfis <- sw_cluster_mfis(fsom)

  # Save cluster plots
  plot_file <- file.path(output_dir, "FlowSOM_clusters.pdf")
  tryCatch({
    sw_plot_clusters(fsom, plot_file)
  }, error = function(e) {
    warning("Could not generate cluster plots: ", e$message, call. = FALSE)
  })

  message("\n=== Pipeline complete ===")
  message("Results contain: ", paste(names(results), collapse = ", "))

  results
}
