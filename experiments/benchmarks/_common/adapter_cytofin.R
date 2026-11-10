# Purpose: wrap the published CytofIn package (cytofin) behind the same inputs and return shape as
#          05_cytether-synthetic/adapters.R (run_cycombine, run_cytonorm).
# Run:     sourced by a benchmark run script; nothing is run here. Self-test at the bottom runs only with
#          Rscript adapter_cytofin.R --selftest
# Inputs:  x (cells by markers, arcsinh cofactor 6000), meta (per cell: batch, sample_id, donor_id,
#          is_reference), markers, anchor, mode
# Outputs: list(x, supported, abstained, reason, fits_on_heldout, mode, hist_delta)
# Pinned:  cytofin 0.0.0.9000, bennyyclo/Cytofin cb3e24b (MIT), wrapped through exported functions only.
#
# How cytofin works (read from the installed function bodies and the README at cb3e24b):
# - It reads and writes FCS files on disk; it has no in-memory entry point. We write temp FCS files.
# - It needs a metadata table (filename, cohort, plate_number, patient_id, condition, is_anchor,
#   validation) and a panel table (metal_name, antigen_name, antigen_pattern, lineage, functional, general).
# - cytofin_prep_anchors() pools every anchor file and stores the pooled per-marker mean and variance
#   (the "universal" statistics) on the asinh(a + b * x) scale.
# - cytofin_normalize() finds, for each file, the anchor with the same cohort and plate_number, moves that
#   anchor's per-marker mean (and variance, per mode) to the universal one, and applies the same map to
#   every cell of the plate. One anchor file per plate: we pool each batch's reference cells into one file.
# - Transform: it applies asinh(shift_factor + scale_factor * x) to raw counts internally (defaults 0 and
#   0.2, the CyTOF cofactor 5), then undoes it with sinh before writing. Our x is already arcsinh(raw / 6000).
#   So we write raw = 6000 * sinh(x) and pass scale_factor = 1 / 6000: the internal scale then equals ours
#   and the round trip returns the corrected value on our scale.
# - Modes: meanshift, meanshift_bulk, variance, z_score, beadlike. README gives no default and says the user
#   should compare them; its worked example uses "meanshift". The function signature lists meanshift first
#   only because match.arg takes the first entry. So mode is a required argument here.

# Parameters ----------------------------------------------------------------------
CYTOFIN_COFACTOR <- 6000 # spectral arcsinh cofactor (CLAUDE.md); cytofin's own default is 5 (scale 0.2)
CYTOFIN_COHORT <- "bench" # cytofin matches anchors on cohort and plate_number; one cohort here

# CytofIn: reference mode, per-batch anchors pooled into one universal target ---------
run_cytofin <- function(x, meta, markers, anchor, mode, ...) {
  stopifnot(
    "mode is required: one of meanshift, meanshift_bulk, variance, z_score, beadlike" =
      !missing(mode) && length(mode) == 1 &&
        mode %in% c("meanshift", "meanshift_bulk", "variance", "z_score", "beadlike"),
    "meta rows match x rows" = nrow(meta) == nrow(x)
  )
  # Only reference cells that are not held out may calibrate. Held-out columns are honoured when present.
  calibrates <- meta$is_reference
  if ("is_holdout_reference" %in% names(meta)) calibrates <- calibrates & !meta$is_holdout_reference
  if ("arrival" %in% names(meta)) calibrates <- calibrates & meta$arrival == "initial"
  known <- meta$batch %in% unique(meta$batch[calibrates]) # batches that own an anchor
  if (!any(calibrates)) {
    return(list(
      x = NULL, supported = FALSE, abstained = NA_integer_, reason = "no reference samples",
      fits_on_heldout = FALSE, mode = mode, hist_delta = NA_real_
    ))
  }

  scale_factor <- 1 / CYTOFIN_COFACTOR
  raw <- sinh(x[, markers, drop = FALSE]) * CYTOFIN_COFACTOR
  colnames(raw) <- markers

  # One FCS per sample; one anchor FCS per batch, built from its calibrating cells only.
  sample_ids <- unique(meta$sample_id[known])
  sample_info <- tibble::tibble(
    sample_id = sample_ids,
    batch = meta$batch[match(sample_ids, meta$sample_id)],
    filename = paste0("s", seq_along(sample_ids), ".fcs"),
    is_anchor = 0L
  )
  anchor_batches <- unique(meta$batch[calibrates])
  anchor_info <- tibble::tibble(
    sample_id = NA_character_,
    batch = anchor_batches,
    filename = paste0("anchor", seq_along(anchor_batches), ".fcs"),
    is_anchor = 1L
  )

  # Fit on calibrating cells, normalize all cells of known batches. `rows_in` limits which samples run.
  cytofin_once <- function(sample_info, tmp) {
    in_dir <- file.path(tmp, "in")
    out_dir <- file.path(tmp, "out")
    dir.create(in_dir, recursive = TRUE)
    write_fcs <- function(m, filename) {
      # flowFrame() re-reads its own temp file and warns about $PnR truncation of that copy only;
      # the file we write keeps every value (checked: residual 2e-7 on a +5 shift).
      ff <- suppressWarnings(flowCore::flowFrame(m))
      flowCore::write.FCS(ff, file.path(in_dir, filename))
    }
    for (i in seq_len(nrow(sample_info))) {
      write_fcs(raw[meta$sample_id == sample_info$sample_id[[i]], , drop = FALSE], sample_info$filename[[i]])
    }
    anchors <- anchor_info[anchor_info$batch %in% sample_info$batch, ]
    for (i in seq_len(nrow(anchors))) {
      write_fcs(raw[calibrates & meta$batch == anchors$batch[[i]], , drop = FALSE], anchors$filename[[i]])
    }
    files <- dplyr::bind_rows(sample_info, anchors) |>
      dplyr::transmute(
        filename,
        cohort = CYTOFIN_COHORT,
        plate_number = batch,
        patient_id = dplyr::coalesce(sample_id, "anchor"),
        condition = "none",
        is_anchor,
        validation = filename
      )
    metadata_path <- file.path(tmp, "metadata.csv")
    panel_path <- file.path(tmp, "panel.csv")
    utils::write.csv(files, metadata_path, row.names = FALSE)
    utils::write.csv(
      data.frame(
        metal_name = markers, antigen_name = markers, antigen_pattern = markers,
        lineage = 1L, functional = 0L, general = 0L
      ),
      panel_path, row.names = FALSE
    )
    # input_prefix = "" because files are written with their final names.
    stats <- cytofin::cytofin_prep_anchors(
      metadata_path, panel_path, in_dir, input_prefix = "", output_path = "none",
      shift_factor = 0, scale_factor = scale_factor
    )
    cytofin::cytofin_normalize(
      metadata_path, panel_path, stats, in_dir, out_dir, mode = mode, input_prefix = "",
      output_prefix = "normalized_", shift_factor = 0, scale_factor = scale_factor
    ) |>
      invisible()
    out <- x
    for (i in seq_len(nrow(sample_info))) {
      rows <- which(meta$sample_id == sample_info$sample_id[[i]])
      ff <- flowCore::read.FCS(
        file.path(out_dir, paste0("normalized_", sample_info$filename[[i]])),
        transformation = FALSE, truncate_max_range = FALSE
      )
      back <- asinh(flowCore::exprs(ff)[, markers, drop = FALSE] / CYTOFIN_COFACTOR)
      stopifnot("cytofin changed the number of cells in a sample" = nrow(back) == length(rows))
      out[rows, markers] <- back
    }
    out
  }

  out <- cytofin_once(sample_info, withr::local_tempdir())

  # Historical invariance, measured: refit on the initial samples and anchors alone, compare their rows.
  hist_delta <- NA_real_
  if ("arrival" %in% names(meta)) {
    old_samples <- sample_info[sample_info$sample_id %in% unique(meta$sample_id[meta$arrival == "initial"]), ]
    out_old <- cytofin_once(old_samples, withr::local_tempdir())
    old_rows <- meta$sample_id %in% old_samples$sample_id
    hist_delta <- max(abs(out_old[old_rows, markers] - out[old_rows, markers]))
  }

  list(
    x = out, supported = TRUE, abstained = sum(!known),
    reason = paste0(
      "mode=", mode, "; per-batch anchor = pooled calibrating reference cells; ",
      "universal target = pool of all batch anchors (the anchor batch ", anchor, " is not special)"
    ),
    fits_on_heldout = FALSE, mode = mode, hist_delta = hist_delta
  )
}

# Self-test: 3 batches x 2 donors, 300 cells each; donor 2 is held out, donor 1 is the reference ----
if (sys.nframe() == 0 && "--selftest" %in% commandArgs(TRUE)) {
  suppressPackageStartupMessages({
    library(cytofin)
    library(dplyr)
  })
  set.seed(1)
  markers <- c("CD3", "CD4", "CD8", "CD19")
  meta <- expand.grid(batch = c("b1", "b2", "b3"), donor_id = c("d1", "d2"), stringsAsFactors = FALSE) |>
    mutate(sample_id = paste(batch, donor_id, sep = "_"), is_reference = donor_id == "d1") |>
    slice(rep(row_number(), each = 300))
  batch_shift <- c(b1 = 0, b2 = 1, b3 = -0.8)
  x <- matrix(rnorm(nrow(meta) * 4, mean = 3), ncol = 4, dimnames = list(NULL, markers)) +
    batch_shift[meta$batch]

  for (m in c("meanshift", "meanshift_bulk", "variance", "z_score", "beadlike")) {
    res <- run_cytofin(x, meta, markers, anchor = "b1", mode = m)
    mean_by_batch <- \(z) sapply(split(seq_len(nrow(z)), meta$batch), \(i) mean(z[i, ]))
    cat(m, "| spread of batch means before:", round(diff(range(mean_by_batch(x))), 3),
        "after:", round(diff(range(mean_by_batch(res$x))), 3), "| abstained:", res$abstained, "\n")
  }

  # Held-out cells must not fit: changing donor d2 cells leaves donor d1 output unchanged.
  x2 <- x
  x2[!meta$is_reference, ] <- x2[!meta$is_reference, ] + 5
  a <- run_cytofin(x, meta, markers, "b1", mode = "meanshift")$x
  b <- run_cytofin(x2, meta, markers, "b1", mode = "meanshift")$x
  cat("max change in reference-cell output when held-out cells shift by 5:",
      signif(max(abs(a[meta$is_reference, ] - b[meta$is_reference, ])), 3), "\n")

  # Meanshift check: output - input is one constant per marker and batch.
  d <- (a - x)[meta$batch == "b2", ]
  cat("meanshift constant offset per marker (sd across cells, should be ~0):",
      signif(max(apply(d, 2, sd)), 3), "\n")

  # A batch without reference cells keeps its values and counts as abstained.
  meta_nb <- meta |> mutate(is_reference = is_reference & batch != "b3")
  r <- run_cytofin(x, meta_nb, markers, "b1", mode = "meanshift")
  cat("batch b3 without reference: abstained =", r$abstained, "unchanged =",
      identical(r$x[meta$batch == "b3", ], x[meta$batch == "b3", ]), "\n")

  # Missing mode must fail.
  cat("missing mode errors:", inherits(try(run_cytofin(x, meta, markers, "b1"), silent = TRUE), "try-error"), "\n")
}
