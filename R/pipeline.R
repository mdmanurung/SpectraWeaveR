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
#'   \item Clustering (kohonen SOM / FastPG)
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
#' @param n_metaclusters Integer; number of metaclusters for clustering
#'   (default: 20).
#' @param seed Integer; random seed for reproducibility (default: 42).
#' @param unmix_from_raw Logical; if \code{TRUE}, run the AutoSpectral
#'   unmixing pipeline on raw spectral FCS files before loading. Requires
#'   \code{control_dir} and \code{unstained_fcs}. Default: \code{FALSE}.
#' @param control_dir Character path to single-stain control FCS files.
#'   Required when \code{unmix_from_raw = TRUE}.
#' @param unstained_fcs Path(s) to unstained FCS file(s) for AF extraction.
#'   Required when \code{unmix_from_raw = TRUE}.
#' @param control_file Character path to a control CSV, or \code{NULL} to
#'   auto-generate. Used when \code{unmix_from_raw = TRUE}.
#' @param cytometer Character; cytometer type (default: \code{"aurora"}).
#'   Used when \code{unmix_from_raw = TRUE}.
#' @param unmix_method Character; unmixing method (default:
#'   \code{"AutoSpectral"}). Used when \code{unmix_from_raw = TRUE}.
#' @param ... Additional arguments.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{unmix_result}}{Unmixing results (if
#'       \code{unmix_from_raw = TRUE})}
#'     \item{\code{flowset}}{The loaded \code{flowSet}}
#'     \item{\code{gating_set}}{The \code{GatingSet} (if gating was applied)}
#'     \item{\code{qc_results}}{QC summary from PeacoQC}
#'     \item{\code{uncorrected}}{Pre-correction tibble}
#'     \item{\code{corrected}}{Post-correction tibble}
#'     \item{\code{correction_eval}}{Batch correction evaluation metrics}
#'     \item{\code{cluster_result}}{Clustering result object (sw_cluster_result)}
#'     \item{\code{cluster_assignments}}{Per-cell cluster assignments}
#'     \item{\code{cluster_mfis}}{MFI table per cluster}
#'   }
#'
#' @examples
#' \dontrun{
#' meta <- data.frame(
#'   file = c("sample1.fcs", "sample2.fcs"),
#'   sample = c("S1", "S2"),
#'   batch = c("B1", "B2")
#' )
#' result <- run_pipeline(
#'   fcs_dir = "path/to/fcs",
#'   sample_meta = meta,
#'   markers = c("CD3", "CD4", "CD8"),
#'   lineage_markers = c("CD3", "CD4", "CD8")
#' )
#' names(result)
#' }
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
                         unmix_from_raw = FALSE,
                         control_dir = NULL,
                         unstained_fcs = NULL,
                         control_file = NULL,
                         cytometer = "aurora",
                         unmix_method = "AutoSpectral",
                         ...) {
  # --- Input validation ---
  if (!is.character(fcs_dir) || length(fcs_dir) != 1) {
    stop("'fcs_dir' must be a single directory path.", call. = FALSE)
  }

  if (!dir.exists(fcs_dir)) {
    stop("FCS directory does not exist: ", fcs_dir, call. = FALSE)
  }

  .validate_df(sample_meta, "sample_meta")

  required_meta_cols <- c("file", "sample", "batch")
  missing_meta <- setdiff(required_meta_cols, names(sample_meta))
  if (length(missing_meta) > 0) {
    stop("'sample_meta' is missing required column(s): ",
         paste(missing_meta, collapse = ", "), call. = FALSE)
  }

  dup_samples <- sample_meta$sample[duplicated(sample_meta$sample)]
  if (length(dup_samples) > 0) {
    stop("'sample_meta' contains duplicate sample name(s): ",
         paste(unique(dup_samples), collapse = ", "), call. = FALSE)
  }

  .validate_markers(markers)

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

  # --- Step 0 (optional): Unmix from raw spectral data ---
  if (unmix_from_raw) {
    message("\n=== Step 0: Unmixing raw spectral data ===")
    if (is.null(control_dir)) {
      stop("'control_dir' is required when unmix_from_raw = TRUE.",
           call. = FALSE)
    }
    unmix_output_dir <- file.path(output_dir, "unmixing")
    unmix_result <- sw_unmix_pipeline(
      control_dir = control_dir,
      sample_input = fcs_dir,
      unstained_fcs = unstained_fcs,
      cytometer = cytometer,
      control_file = control_file,
      method = unmix_method,
      output_dir = unmix_output_dir
    )
    results$unmix_result <- unmix_result

    # Update fcs_dir to point to unmixed output
    fcs_dir <- unmix_result$unmixed$output_dir
    message("Using unmixed files from: ", fcs_dir)
  }

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

    # Build a flowSet from margin-removed flowFrames so gating operates
    # on cleaned data rather than re-reading raw FCS files from disk.
    fs_clean <- flowCore::flowSet(ff_list)
    gs <- flowWorkspace::GatingSet(fs_clean)

    # Apply gating template
    if (!file.exists(gating_template)) {
      stop("Gating template file not found: ", gating_template, call. = FALSE)
    }
    gt <- openCyto::gatingTemplate(gating_template)
    openCyto::gt_gating(gt, gs)
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
  message("\n=== Step 6/6: Clustering ===")

  cluster_result <- sw_cluster(
    corrected,
    lineage_markers = lineage_markers,
    n_metaclusters = n_metaclusters,
    seed = seed
  )
  results$cluster_result <- cluster_result

  results$cluster_assignments <- sw_get_cluster_assignments(cluster_result)
  results$cluster_mfis <- sw_cluster_mfis(cluster_result)

  # Save cluster plots
  plot_file <- file.path(output_dir, "cluster_heatmap.pdf")
  tryCatch({
    sw_plot_clusters(cluster_result, plot_file)
  }, error = function(e) {
    warning("Could not generate cluster plots: ", e$message, call. = FALSE)
  })

  message("\n=== Pipeline complete ===")
  message("Results contain: ", paste(names(results), collapse = ", "))

  results
}
