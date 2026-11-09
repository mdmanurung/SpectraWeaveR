# Purpose: wrap published competitors (cyCombine, CytoNorm) behind the same inputs as run_method.R.
# Run:     sourced by run_method.R; nothing is run here.
# Inputs:  x (cells by markers), meta (cell_metadata with arrival, is_reference, is_holdout_reference, batch, sample_id), anchor
# Outputs: each function returns list(x, supported, abstained, reason, hist_delta)
# Settings approved by the user: H1 cyCombine package defaults with ref.batch = anchor, no covar;
# H2 CytoNorm 2.0.12 defaults (15x15 SOM, 10 clusters, 99 quantiles) with seed 431.
# Pinned: cyCombine fcd42db (0.3.0), CytoNorm 278f352 (2.0.12), both wrapped, never copied.

# Parameters ----------------------------------------------------------------------
CYCOMBINE_SEED <- 473 # cyCombine package default
CYTONORM_SEED <- 431 # bundle simulation seed; the package default is NULL
CYTONORM_SOM <- 15L # package default grid (CytoNorm 2.0 Table 1)
CYTONORM_CLUSTERS <- 10L # package default
CYTONORM_QUANTILES <- 99L # package default

# cyCombine: refits on all cells it is given, so it cannot freeze --------------------
.cycombine_once <- function(x, batch, markers, ref_batch) {
  df <- tibble::as_tibble(x)
  df$batch <- batch
  df$id <- seq_len(nrow(df))
  res <- cyCombine::batch_correct(
    df, markers = markers, norm_method = "scale", ref.batch = ref_batch,
    seed = CYCOMBINE_SEED, covar = NULL
  )
  stopifnot("cyCombine changed the number of cells" = nrow(res) == nrow(df))
  out <- as.matrix(res[order(res$id), markers])
  # cyCombine leaves NA for cells in clusters too small to correct: keep the original value (abstain).
  attr(out, "n_na") <- sum(is.na(out))
  out[is.na(out)] <- x[is.na(out)]
  out
}

# H1b: reference mode first (ref.batch = anchor); if ComBat errors because a cluster lacks the anchor
# batch, run with no reference batch and record the error in `reason`.
run_cycombine <- function(x, meta, markers, anchor) {
  ref_batch <- anchor
  reason <- paste0("ref.batch=", anchor)
  out_all <- tryCatch(
    .cycombine_once(x, meta$batch, markers, ref_batch),
    error = function(e) {
      ref_batch <<- NULL
      reason <<- paste0("ref.batch=NULL; reference mode failed: ", gsub("[\r\n,]+", " ", conditionMessage(e)))
      .cycombine_once(x, meta$batch, markers, NULL)
    }
  )
  # Historical invariance: correct the initial samples alone and compare with their rows in the full run.
  old <- meta$arrival == "initial"
  out_old <- tryCatch(
    .cycombine_once(x[old, , drop = FALSE], meta$batch[old], markers, ref_batch),
    error = function(e) NULL
  )
  list(
    x = out_all, supported = TRUE, abstained = attr(out_all, "n_na"), reason = reason,
    hist_delta = if (is.null(out_old)) NA_real_ else max(abs(c(out_old) - c(out_all[old, , drop = FALSE])))
  )
}

# CytoNorm: train on non-held-out initial reference cells, apply to every sample whose batch it saw --
run_cytonorm <- function(x, meta, markers) {
  calibrates <- meta$is_reference & !meta$is_holdout_reference & meta$arrival == "initial"
  if (!any(calibrates)) {
    return(list(x = NULL, supported = FALSE, abstained = NA_integer_, reason = "no reference samples"))
  }
  as_set <- function(rows) {
    frames <- lapply(split(which(rows), meta$sample_id[rows]), function(i) {
      flowCore::flowFrame(x[i, , drop = FALSE])
    })
    methods::as(frames, "flowSet")
  }
  train_ids <- unique(meta$sample_id[calibrates])
  train_batch <- meta$batch[match(train_ids, meta$sample_id)]
  train_set <- as_set(calibrates)[train_ids]
  tmp <- tempfile("cytonorm_")
  dir.create(tmp)
  model <- CytoNorm::CytoNorm.train(
    files = train_set, labels = train_batch, channels = markers, transformList = NULL,
    outputDir = tmp, FlowSOM.params = list(nCells = 1e6, xdim = CYTONORM_SOM, ydim = CYTONORM_SOM,
                                           nClus = CYTONORM_CLUSTERS, scale = FALSE),
    normParams = list(nQ = CYTONORM_QUANTILES), seed = CYTONORM_SEED, clean = TRUE, plot = FALSE, verbose = FALSE
  )
  # Samples from a batch with no training reference keep their values (abstain).
  known <- meta$batch %in% unique(train_batch)
  all_ids <- unique(meta$sample_id[known])
  all_set <- as_set(known)[all_ids]
  all_batch <- meta$batch[match(all_ids, meta$sample_id)]
  normed <- CytoNorm::CytoNorm.normalize(
    model, files = all_set, labels = all_batch, transformList = NULL, transformList.reverse = NULL,
    outputDir = tmp, write = FALSE, clean = TRUE, verbose = FALSE
  )
  out <- x
  stopifnot("CytoNorm returned a different number of samples" = length(normed) == length(all_ids))
  for (i in seq_along(all_ids)) {
    rows <- which(meta$sample_id == all_ids[[i]])
    ff <- normed[[i]] # same order as all_ids
    stopifnot("CytoNorm changed the number of cells in a sample" = nrow(ff) == length(rows))
    out[rows, ] <- flowCore::exprs(ff)[, markers, drop = FALSE]
  }
  list(x = out, supported = TRUE, abstained = sum(!known), reason = "train on initial non-holdout references",
       hist_delta = 0) # applying a frozen model to old samples cannot depend on new samples
}
