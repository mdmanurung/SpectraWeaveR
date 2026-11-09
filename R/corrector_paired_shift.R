# Paired reference mean shift engine for the corrector contract.
#
# The idea (shift each batch onto an anchor batch by the mean paired
# difference of repeated reference donors) comes from the cyTether
# implementation bundle v1 (MIT), r-package/cyTether/R/tether.R and
# docs/03_statistical_spec.md. The code is re-derived. It does not inherit
# the bundle's aliquot collapse, cell-count variance, hidden constants or
# default anchor.

#' Paired reference mean shift engine
#'
#' \code{engine = "paired_shift"} in \code{\link{sw_corrector_fit}}. It
#' moves every batch onto the anchor batch by an additive shift per
#' (batch, population, marker), estimated only from reference samples that
#' were measured in both that batch and the anchor.
#'
#' @details
#' Estimation, on data already on the model's transformed scale:
#' \enumerate{
#'   \item Aliquot summary: the mean of the cells of one (reference id,
#'     aliquot, batch, population), per marker. Summaries with fewer than
#'     \code{min_cells} cells are excluded and listed in \code{excluded}
#'     with the reason. Aliquots are never merged: two aliquots of one
#'     reference in one batch give two summaries. When the manifest has no
#'     \code{aliquot} column the aliquot key falls back to the reference id,
#'     and the model's \code{aliquot_source} says so.
#'   \item Donor summary: the mean of the aliquot summaries of one reference
#'     in one batch, each aliquot weighted equally. Cells are measurements,
#'     not replicates.
#'   \item Effect of batch \eqn{b}, population \eqn{k}, marker \eqn{m}: the
#'     mean over references \eqn{r} with a donor summary in both \eqn{b} and
#'     the anchor of \eqn{z_{rbkm} - z_{r,anchor,km}}. The anchor effect is
#'     zero. \code{se_donor} is the standard deviation of those differences
#'     over references divided by \eqn{\sqrt{n_{donors}}}, and \code{NA} for
#'     one reference. No cell-count variance and no variance floor are used.
#'     The optional column \code{df} (\code{n_donors - 1}) gives the degrees
#'     of freedom of \code{se_donor}, for t-based intervals.
#'   \item Correction: \code{corrected = x - effect} for every cell of that
#'     batch and population, reference or participant. Participant cells are
#'     never used as calibration evidence.
#' }
#' A new batch at \code{\link{sw_corrector_integrate}} is estimated the same
#' way from its own incoming reference cells against the frozen anchor
#' donor summaries (\code{model$targets}). Effects of batches already in the
#' model are never re-estimated from new rows.
#'
#' Limits. The model is a pure location shift. On asinh-transformed data an
#' additive shift equals a multiplicative gain change only where the
#' population is bright (values well above the cofactor, where asinh is
#' close to a log). For dim or negative populations, where asinh is close
#' to linear, a gain change is not additive, and one shift per population
#' will under- or over-correct parts of the distribution. Scale changes,
#' nonlinear drift and preparation effects are not modelled. A population
#' with no reference cells in the anchor gets no effect and is abstained.
#' \code{hierarchy} must be \code{NULL}: this engine has no study or run
#' terms.
#'
#' @name sw_corrector_paired_shift
#' @seealso \code{\link{sw_corrector_fit}}, \code{\link{sw_corrector_integrate}}
NULL

# Aliquot summaries in long form, plus the summaries below min_cells.
.sw_ps_summaries <- function(xm, cells, min_cells) {
  markers <- colnames(xm)
  ref <- cells[cells$is_reference, , drop = FALSE]
  empty <- tibble::tibble(
    reference_id = character(),
    aliquot = character(),
    batch = character(),
    population = character(),
    marker = character(),
    value = double(),
    n_cells = integer()
  )
  if (nrow(ref) == 0) {
    return(list(summ = empty, excluded = .sw_empty_excluded()))
  }
  key <- paste(
    ref$reference_id,
    ref$aliquot,
    ref$batch,
    ref$population,
    sep = "\r"
  )
  groups <- split(ref$row, key)
  info <- ref[
    match(names(groups), key),
    c(
      "reference_id",
      "aliquot",
      "batch",
      "population"
    )
  ]
  n <- lengths(groups)
  keep <- n >= min_cells

  excluded <- tibble::tibble(
    info[!keep, , drop = FALSE],
    n_cells = as.integer(n[!keep]),
    reason = sprintf(
      "%d cell(s), below min_cells = %d; excluded from the effect estimate.",
      n[!keep],
      min_cells
    )
  )
  if (!any(keep)) {
    return(list(summ = empty, excluded = excluded))
  }

  # Mean of sorted values, so a summary does not depend on row order.
  means <- do.call(
    rbind,
    lapply(groups[keep], function(ii) {
      apply(xm[ii, , drop = FALSE], 2, function(v) mean(sort(v)))
    })
  )
  k <- which(keep)
  n_m <- length(markers)
  summ <- tibble::tibble(
    reference_id = rep(info$reference_id[k], each = n_m),
    aliquot = rep(info$aliquot[k], each = n_m),
    batch = rep(info$batch[k], each = n_m),
    population = rep(info$population[k], each = n_m),
    marker = rep(markers, length(k)),
    value = as.vector(t(means)),
    n_cells = rep(as.integer(n[k]), each = n_m)
  )
  list(summ = summ, excluded = excluded)
}

# Donor summaries: equal-weight mean of aliquot summaries per reference and batch.
.sw_ps_donor_level <- function(summ) {
  key <- paste(
    summ$reference_id,
    summ$batch,
    summ$population,
    summ$marker,
    sep = "\r"
  )
  g <- unname(split(seq_len(nrow(summ)), key))
  first <- vapply(g, `[`, integer(1), 1)
  tibble::tibble(
    summ[first, c("reference_id", "batch", "population", "marker")],
    value = vapply(g, function(ii) mean(sort(summ$value[ii])), numeric(1)),
    n_aliquots = as.integer(lengths(g))
  )
}

# Effects: mean over shared references of (batch summary - anchor summary).
.sw_ps_effects <- function(donor, targets, source) {
  paired <- dplyr::inner_join(
    donor,
    targets,
    by = c("reference_id", "population", "marker"),
    relationship = "many-to-one"
  )
  if (nrow(paired) == 0) {
    return(.sw_empty_effects())
  }
  diff <- paired$value - paired$target
  key <- paste(paired$batch, paired$population, paired$marker, sep = "\r")
  g <- unname(split(seq_len(nrow(paired)), key))
  first <- vapply(g, `[`, integer(1), 1)
  tibble::tibble(
    batch = paired$batch[first],
    population = paired$population[first],
    marker = paired$marker[first],
    estimate = vapply(g, function(ii) mean(sort(diff[ii])), numeric(1)),
    se_donor = vapply(
      g,
      function(ii) {
        if (length(ii) < 2) {
          NA_real_
        } else {
          stats::sd(sort(diff[ii])) / sqrt(length(ii))
        }
      },
      numeric(1)
    ),
    # Degrees of freedom of se_donor: the SD over n references has n - 1.
    df = ifelse(lengths(g) >= 2, lengths(g) - 1, NA_real_),
    n_donors = as.integer(lengths(g)),
    n_aliquots = vapply(g, function(ii) sum(paired$n_aliquots[ii]), integer(1)),
    n_aliquots_anchor = vapply(
      g,
      function(ii) sum(paired$n_aliquots_anchor[ii]),
      integer(1)
    ),
    support = "paired_reference",
    source = source
  )
}

# Engine fit: anchor targets from reference donors, then one shift per batch.
.sw_ps_fit <- function(xm, cells, anchor, min_cells, hierarchy, seed, ...) {
  if (!is.null(hierarchy)) {
    stop("paired_shift has no hierarchy; set hierarchy = NULL.", call. = FALSE)
  }
  s <- .sw_ps_summaries(xm, cells, min_cells)
  if (nrow(s$summ) == 0) {
    stop(
      "No reference summary has at least min_cells cells. paired_shift needs ",
      "samples with reference_role 'technical_replicate' or 'anchor' and a ",
      "reference_id.",
      call. = FALSE
    )
  }
  donor <- .sw_ps_donor_level(s$summ)
  anc <- donor[donor$batch == anchor, , drop = FALSE]
  if (nrow(anc) == 0) {
    stop(
      sprintf(
        "Anchor batch '%s' has no reference summary with at least min_cells cells.",
        anchor
      ),
      call. = FALSE
    )
  }
  targets <- tibble::tibble(
    reference_id = anc$reference_id,
    population = anc$population,
    marker = anc$marker,
    target = anc$value,
    n_aliquots_anchor = anc$n_aliquots
  )
  effects <- .sw_ps_effects(donor, targets, source = "fit")
  effects$support[effects$batch == anchor] <- "anchor"
  list(effects = effects, targets = targets, excluded = s$excluded)
}

# Engine apply: estimate only unseen batches from their own reference cells.
.sw_ps_apply <- function(model, xm, cells, new_batches) {
  bridging <- cells$batch %in%
    new_batches &
    cells$is_reference &
    cells$reference_id %in% model$targets$reference_id
  s <- .sw_ps_summaries(xm, cells[bridging, , drop = FALSE], model$min_cells)
  if (nrow(s$summ) == 0) {
    return(list(effects = .sw_empty_effects(), excluded = s$excluded))
  }
  donor <- .sw_ps_donor_level(s$summ)
  list(
    effects = .sw_ps_effects(donor, model$targets, source = "integrate"),
    excluded = s$excluded
  )
}

sw_corrector_register_engine(
  "paired_shift",
  fit = .sw_ps_fit,
  apply = .sw_ps_apply,
  can_freeze = TRUE,
  description = "Paired reference mean shift onto an anchor batch (additive, per batch, population and marker)."
)
