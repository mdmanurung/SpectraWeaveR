# Event-level unmixing metrics: negatives, spread, stain index, pairwise line
# check (sw_unmix_event_metrics), the edges they rest on (sw_unmix_edges), and
# held-out control error with extraction reproducibility (sw_unmix_heldout).
# Ported from benchmark 05 scripts 101f, 103c, 103d, 103f and 103g; plan
# plans/2026-11-15_unmix-event-metrics.md. AutoSpectral (AGPL) is called, never
# copied: the spread-aware edge is AutoSpectral::get.spread.thresholds().


# ---- Internal helpers ---------------------------------------------------------

# Arguments of `call_args` the caller did not pass, for one stop message.
.sw_em_missing <- function(fn, call_env, names) {
  miss <- names[vapply(names, \(n) eval(call("missing", as.name(n)), call_env), logical(1))]
  if (length(miss) > 0L) {
    stop(
      sprintf(
        "%s() needs every argument; missing: %s. No default is built in (see ?%s).",
        fn,
        paste(miss, collapse = ", "),
        fn
      ),
      call. = FALSE
    )
  }
}

.sw_em_scalar <- function(x, arg, lower = -Inf, integer = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x > lower
  if (ok && integer) ok <- x == round(x)
  if (!ok) {
    stop(
      sprintf(
        "`%s` must be one finite %s%s.",
        arg,
        if (integer) "whole number" else "number",
        if (is.finite(lower)) paste(" above", lower) else ""
      ),
      call. = FALSE
    )
  }
  invisible(x)
}

.sw_em_probs <- function(p, arg) {
  if (!is.numeric(p) || length(p) != 2L || anyNA(p) || any(p <= 0 | p >= 1) || p[[1]] >= p[[2]]) {
    stop(sprintf("`%s` must be two increasing probabilities in (0, 1).", arg), call. = FALSE)
  }
  invisible(p)
}

# Named list of numeric event matrices from a list, flowFrame or flowSet.
.sw_em_as_list <- function(data) {
  if (methods::is(data, "flowFrame")) {
    return(stats::setNames(list(flowCore::exprs(data)), flowCore::identifier(data)))
  }
  if (methods::is(data, "flowSet")) {
    return(stats::setNames(
      lapply(seq_along(data), \(i) flowCore::exprs(data[[i]])),
      flowCore::sampleNames(data)
    ))
  }
  ok <- is.list(data) &&
    length(data) > 0L &&
    !is.null(names(data)) &&
    all(nzchar(names(data))) &&
    !anyDuplicated(names(data)) &&
    all(vapply(data, \(m) is.matrix(m) && is.numeric(m) && !is.null(colnames(m)), logical(1)))
  if (!ok) {
    stop(
      "`data` must be a flowFrame, a flowSet, or a uniquely named list of numeric matrices with column names.",
      call. = FALSE
    )
  }
  data
}

# Robust slope: median of pairwise slopes (Theil-Sen).
.sw_theil_sen <- function(x, y) {
  if (length(x) < 2L) {
    return(NA_real_)
  }
  ij <- utils::combn(length(x), 2L)
  dx <- x[ij[2L, ]] - x[ij[1L, ]]
  ok <- dx != 0
  if (!any(ok)) {
    return(NA_real_)
  }
  stats::median((y[ij[2L, ]] - y[ij[1L, ]])[ok] / dx[ok])
}

# Spread-aware edge matrix from AutoSpectral (one place, so tests can mock it).
.sw_em_spread_edges <- function(x, thresholds, spread, kappa, side) {
  if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
    stop(
      "Spread-aware edges need AutoSpectral (>= 1.8.8): remotes::install_github('DrCytometer/AutoSpectral'). Pass spread = NULL for flat edges.",
      call. = FALSE
    )
  }
  AutoSpectral::get.spread.thresholds(
    unmixed = x,
    thresholds = thresholds,
    spillover.spread = spread,
    spread.kappa = kappa,
    side = side,
    verbose = FALSE
  )
}

# One AutoSpectral call by name (one place, so tests can mock it).
.sw_ho_call <- function(fn, args) {
  do.call(getExportedValue("AutoSpectral", fn), args)
}

# Width of the central interval between two standard-normal quantiles.
.sw_z_width <- function(p) stats::qnorm(p[[2]]) - stats::qnorm(p[[1]])


# ---- Edges --------------------------------------------------------------------

#' Per-channel floor and positive edge for unmix event metrics
#'
#' Gives the linear-scale negative floor and positive edge of each channel, as
#' AutoSpectral 1.8.8 defines them: quantiles of the unmixed unstained control
#' (`get.spectral.variants()` stores q0.005 as `neg.thresholds` and q0.995 as
#' `thresholds`). Read them from a variants object, or compute them from an
#' unmixed unstained matrix at quantiles you choose.
#'
#' @param variants An AutoSpectral spectral-variants list with named numeric
#'   `thresholds` and `neg.thresholds`, or `NULL`.
#' @param unstained A numeric matrix (events x channels, column names, linear
#'   scale) of an unmixed unstained control, or `NULL`. Give exactly one of
#'   `variants` and `unstained`.
#' @param probs With `unstained`: named numeric `c(floor = , edge = )`
#'   quantiles (AutoSpectral uses 0.005 and 0.995). Must be `NULL` with
#'   `variants`.
#'
#' @return A tibble with `channel`, `floor`, `edge` (linear) and `source`.
#' @seealso [sw_unmix_event_metrics()]
#' @examples
#' set.seed(1)
#' un <- matrix(rnorm(2000, 0, 300), ncol = 2, dimnames = list(NULL, c("BV421", "PE")))
#' sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
#' @export
sw_unmix_edges <- function(variants = NULL, unstained = NULL, probs = NULL) {
  if (is.null(variants) == is.null(unstained)) {
    stop("Give exactly one of `variants` or `unstained`.", call. = FALSE)
  }
  if (!is.null(variants)) {
    if (!is.null(probs)) {
      stop("`probs` is used only with `unstained`; `variants` already holds AutoSpectral's quantiles.", call. = FALSE)
    }
    th <- variants[["thresholds"]]
    nt <- variants[["neg.thresholds"]]
    ok <- is.numeric(th) && is.numeric(nt) && !is.null(names(th)) && setequal(names(th), names(nt))
    if (!ok) {
      stop("`variants` needs named numeric `thresholds` and `neg.thresholds` over the same channels.", call. = FALSE)
    }
    out <- tibble::tibble(channel = names(th), floor = unname(nt[names(th)]), edge = unname(th), source = "variants")
  } else {
    if (!is.matrix(unstained) || !is.numeric(unstained) || is.null(colnames(unstained))) {
      stop("`unstained` must be a numeric matrix with column names.", call. = FALSE)
    }
    ok <- is.numeric(probs) && identical(sort(names(probs)), c("edge", "floor")) &&
      !anyNA(probs) && all(probs > 0 & probs < 1) && probs[["floor"]] < probs[["edge"]]
    if (!ok) {
      stop("`probs` must be named c(floor = , edge = ) with 0 < floor < edge < 1.", call. = FALSE)
    }
    q <- apply(unstained, 2L, stats::quantile, probs = unname(probs[c("floor", "edge")]), names = FALSE)
    out <- tibble::tibble(
      channel = colnames(unstained),
      floor = q[1L, ],
      edge = q[2L, ],
      source = sprintf("unstained q%s / q%s", format(probs[["floor"]]), format(probs[["edge"]]))
    )
  }
  bad <- out$channel[is.na(out$floor) | is.na(out$edge)]
  if (length(bad) > 0L) {
    stop(sprintf("Floor or edge is NA for: %s.", paste(bad, collapse = ", ")), call. = FALSE)
  }
  out
}


# ---- Event metrics --------------------------------------------------------------

#' Event-level unmixing metrics: negatives, spread, stain index, line check
#'
#' Scores unmixed events per file and channel, and optionally per channel pair.
#' Positives are events at or above the channel edge, negatives the rest; with
#' `spread` and `kappa` the edge (and floor) widen per event as in
#' `AutoSpectral::get.spread.thresholds()`. Masks use the linear scale; every
#' other metric uses `asinh(x / cofactor)`.
#'
#' Per file x channel:
#' * `pct_below_cut`: % of events with asinh value below `negative_cut`.
#' * `pct_below_floor`: % of events below the floor (AutoSpectral q0.005 of the
#'   unstained when the edges come from variants).
#' * `neg_median`, `neg_rsd` (MAD) of negatives; `spread_q` = width of the
#'   `spread_probs` interval of negatives divided by its standard-normal width
#'   (3.92 for 0.025 / 0.975); `pileup` = `neg_rsd / spread_q`.
#' * `hyper_share`: share of negatives below `neg_median - hyper_k * neg_rsd`.
#' * `si_mad` = (median pos - median neg) / (2 MAD neg), the convention of
#'   `AutoSpectral::calculate.ssi()` with the sample's own negatives;
#'   `si_q` = the same over `2 * spread_q`. Both `NA` below `min_pos` positives.
#'
#' Per file x pair (`pairs`): in events positive for `channel_i` and negative
#' for `channel_j`, the Theil-Sen slope of the `channel_j` medians on the
#' `channel_i` medians of `n_bins` quantile bins, and the Spearman rho of the
#' events. A flat, uncorrelated line (slope and rho near 0) is what adequate
#' unmixing gives; `NA` below `min_cells` events.
#'
#' No cutoff is built in; [sw_unmix_diagnose()] judges these only against
#' limits you pass.
#'
#' @param data Unmixed events, linear scale: a flowFrame, a flowSet, or a
#'   uniquely named list of numeric matrices with column names.
#' @param edges Output of [sw_unmix_edges()]. `channel` matches a data column
#'   as is or with `-A` appended.
#' @param cofactor arcsinh cofactor (6000 for Aurora data in this project).
#' @param negative_cut Fixed cut on the asinh scale (e.g. -1.5).
#' @param spread_probs Two probabilities for `spread_q` (e.g. 0.025, 0.975).
#' @param min_pos Fewest positives for a stain index.
#' @param hyper_k Robust SDs below the negative median that count as
#'   hypernegative.
#' @param spread `NULL` for flat edges, or AutoSpectral's spillover-spread
#'   matrix (`variants$spillover.spread`) for spread-aware edges.
#' @param kappa `NULL` with `spread = NULL`; otherwise the `spread.kappa` of
#'   `get.spread.thresholds()` (AutoSpectral's default is 2).
#' @param pairs `NULL` (skip the line check), `"all"` (every ordered pair of
#'   matched channels) or a data frame with columns `i` and `j` (edge names).
#' @param n_bins Quantile bins for the slope (`NULL` when `pairs` is `NULL`).
#' @param min_cells Fewest i+ j- events for a slope (`NULL` when `pairs` is
#'   `NULL`).
#'
#' @return An `sw_event_metrics` list: `channels` (tibble, one row per file x
#'   channel), `pairs` (tibble or `NULL`), `missing` (file x channel with no
#'   matching column) and `settings`.
#' @seealso [sw_unmix_edges()], [sw_unmix_diagnose()], [sw_unmix_heldout()]
#' @examples
#' set.seed(1)
#' ch <- c("BV421", "PE")
#' un <- matrix(rnorm(4000, 0, 300), ncol = 2, dimnames = list(NULL, ch))
#' edges <- sw_unmix_edges(unstained = un, probs = c(floor = 0.005, edge = 0.995))
#' s1 <- rbind(un, cbind(rnorm(500, 4e4, 8e3), rnorm(500, 0, 300)))
#' em <- sw_unmix_event_metrics(
#'   list(s1 = s1), edges, cofactor = 6000, negative_cut = -1.5,
#'   spread_probs = c(0.025, 0.975), min_pos = 100, hyper_k = 3,
#'   spread = NULL, kappa = NULL, pairs = "all", n_bins = 10, min_cells = 100
#' )
#' em$channels
#' em$pairs
#' @export
sw_unmix_event_metrics <- function(
  data,
  edges,
  cofactor,
  negative_cut,
  spread_probs,
  min_pos,
  hyper_k,
  spread,
  kappa,
  pairs,
  n_bins,
  min_cells
) {
  .sw_em_missing(
    "sw_unmix_event_metrics",
    environment(),
    c("data", "edges", "cofactor", "negative_cut", "spread_probs", "min_pos", "hyper_k", "spread", "kappa", "pairs", "n_bins", "min_cells")
  )
  data <- .sw_em_as_list(data)
  if (!is.data.frame(edges) || !all(c("channel", "floor", "edge") %in% names(edges)) || anyDuplicated(edges$channel)) {
    stop("`edges` must come from sw_unmix_edges(): unique `channel`, `floor`, `edge`.", call. = FALSE)
  }
  .sw_em_scalar(cofactor, "cofactor", lower = 0)
  .sw_em_scalar(negative_cut, "negative_cut")
  .sw_em_probs(spread_probs, "spread_probs")
  .sw_em_scalar(min_pos, "min_pos", lower = 0, integer = TRUE)
  .sw_em_scalar(hyper_k, "hyper_k", lower = 0)
  if (is.null(spread) != is.null(kappa)) {
    stop("Give `spread` and `kappa` together (spread-aware edges) or both as NULL (flat edges).", call. = FALSE)
  }
  if (!is.null(spread)) {
    if (!is.matrix(spread) || !is.numeric(spread) || is.null(colnames(spread))) {
      stop("`spread` must be AutoSpectral's numeric spillover-spread matrix with dimnames.", call. = FALSE)
    }
    .sw_em_scalar(kappa, "kappa", lower = 0)
    absent <- setdiff(edges$channel, colnames(spread))
    if (length(absent) > 0L) {
      stop(sprintf("`spread` has no column for: %s.", paste(absent, collapse = ", ")), call. = FALSE)
    }
  }
  if (!is.null(pairs)) {
    .sw_em_scalar(n_bins, "n_bins", lower = 1, integer = TRUE)
    .sw_em_scalar(min_cells, "min_cells", lower = 0, integer = TRUE)
    if (min_cells < n_bins) {
      stop("`min_cells` must be at least `n_bins`.", call. = FALSE)
    }
    if (is.character(pairs)) {
      if (!identical(pairs, "all")) stop("`pairs` must be NULL, \"all\" or a data frame with `i` and `j`.", call. = FALSE)
    } else if (!is.data.frame(pairs) || !all(c("i", "j") %in% names(pairs))) {
      stop("`pairs` must be NULL, \"all\" or a data frame with `i` and `j`.", call. = FALSE)
    } else {
      unknown <- setdiff(c(pairs$i, pairs$j), edges$channel)
      if (length(unknown) > 0L) {
        stop(sprintf("`pairs` names channels not in `edges`: %s.", paste(unique(unknown), collapse = ", ")), call. = FALSE)
      }
      if (any(pairs$i == pairs$j)) stop("`pairs` has rows with i == j.", call. = FALSE)
    }
  } else if (!is.null(n_bins) || !is.null(min_cells)) {
    stop("`n_bins` and `min_cells` must be NULL when `pairs` is NULL.", call. = FALSE)
  }
  zw <- .sw_z_width(spread_probs)

  per_file <- lapply(names(data), \(nm) {
    m <- data[[nm]]
    col <- .sw_match_channels(edges$channel, colnames(m))
    found <- !is.na(col)
    if (!any(found)) {
      stop(sprintf("File %s: no column matches any edge channel.", nm), call. = FALSE)
    }
    e <- edges[found, ]
    x <- m[, col[found], drop = FALSE]
    colnames(x) <- e$channel
    up <- stats::setNames(e$edge, e$channel)
    lo <- stats::setNames(e$floor, e$channel)
    if (is.null(spread)) {
      pos <- sweep(x, 2L, up, `>=`)
      below <- sweep(x, 2L, lo, `<`)
    } else {
      pos <- x >= .sw_em_spread_edges(x, up, spread, kappa, "upper")
      below <- x < .sw_em_spread_edges(x, lo, spread, kappa, "lower")
    }
    a <- asinh(x / cofactor)
    ch_rows <- lapply(e$channel, \(k) {
      xa <- a[, k]
      neg <- xa[!pos[, k]]
      pv <- xa[pos[, k]]
      n_neg <- length(neg)
      med <- if (n_neg > 0L) stats::median(neg) else NA_real_
      rsd <- if (n_neg > 1L) stats::mad(neg) else NA_real_
      sq <- if (n_neg > 1L) unname(diff(stats::quantile(neg, spread_probs))) / zw else NA_real_
      dmed <- if (length(pv) >= min_pos) stats::median(pv) - med else NA_real_
      tibble::tibble(
        file = nm,
        channel = k,
        column = unname(col[found][e$channel == k]),
        n = length(xa),
        pct_below_cut = 100 * mean(xa < negative_cut),
        pct_below_floor = 100 * mean(below[, k]),
        n_neg = n_neg,
        neg_median = med,
        neg_rsd = rsd,
        spread_q = sq,
        pileup = rsd / sq,
        hyper_share = if (n_neg > 1L) mean(neg < med - hyper_k * rsd) else NA_real_,
        n_pos = length(pv),
        si_mad = dmed / (2 * rsd),
        si_q = dmed / (2 * sq)
      )
    })
    pair_rows <- NULL
    if (!is.null(pairs)) {
      pr <- if (is.character(pairs)) {
        g <- expand.grid(i = e$channel, j = e$channel, stringsAsFactors = FALSE)
        g[g$i != g$j, ]
      } else {
        as.data.frame(pairs)[pairs$i %in% e$channel & pairs$j %in% e$channel, c("i", "j")]
      }
      pair_rows <- lapply(seq_len(nrow(pr)), \(r) {
        i <- pr$i[[r]]
        j <- pr$j[[r]]
        idx <- pos[, i] & !pos[, j]
        n <- sum(idx)
        slope <- NA_real_
        rho <- NA_real_
        if (n >= min_cells) {
          xi <- a[idx, i]
          xj <- a[idx, j]
          br <- unique(stats::quantile(xi, seq(0, 1, length.out = n_bins + 1L), names = FALSE))
          if (length(br) >= 3L) {
            b <- cut(xi, br, include.lowest = TRUE)
            mx <- tapply(xi, b, stats::median)
            my <- tapply(xj, b, stats::median)
            keep <- !is.na(mx) & !is.na(my)
            slope <- .sw_theil_sen(unname(mx[keep]), unname(my[keep]))
          }
          rho <- suppressWarnings(stats::cor(xi, xj, method = "spearman"))
        }
        tibble::tibble(file = nm, channel_i = i, channel_j = j, n = n, slope = slope, rho = rho)
      })
    }
    list(
      channels = dplyr::bind_rows(ch_rows),
      pairs = if (!is.null(pair_rows)) dplyr::bind_rows(pair_rows),
      missing = tibble::tibble(file = rep(nm, sum(!found)), channel = edges$channel[!found])
    )
  })

  structure(
    list(
      channels = dplyr::bind_rows(lapply(per_file, `[[`, "channels")),
      pairs = if (!is.null(pairs)) dplyr::bind_rows(lapply(per_file, `[[`, "pairs")),
      missing = dplyr::bind_rows(lapply(per_file, `[[`, "missing")),
      settings = list(
        cofactor = cofactor,
        negative_cut = negative_cut,
        spread_probs = spread_probs,
        min_pos = min_pos,
        hyper_k = hyper_k,
        edge_type = if (is.null(spread)) "flat" else "spread-aware",
        kappa = kappa,
        edge_source = unique(edges$source),
        n_bins = n_bins,
        min_cells = min_cells
      )
    ),
    class = "sw_event_metrics"
  )
}

#' @export
print.sw_event_metrics <- function(x, ...) {
  cat("<sw_event_metrics>", length(unique(x$channels$file)), "file(s),", length(unique(x$channels$channel)), "channel(s)\n")
  cat("  edges:", x$settings$edge_type, "| negative cut", x$settings$negative_cut, "asinh, cofactor", x$settings$cofactor, "\n")
  cat("  pairs:", if (is.null(x$pairs)) "not run" else paste(nrow(x$pairs), "rows"), "\n")
  if (nrow(x$missing) > 0L) cat("  missing channel x file:", nrow(x$missing), "\n")
  invisible(x)
}


# ---- Held-out control error --------------------------------------------------------

#' Held-out control error and extraction reproducibility (AutoSpectral)
#'
#' The G1 check of benchmark 05: split every single-stain control 50/50 (or at
#' `train_fraction`), extract spectra, AF spectra and variants from the train
#' halves with AutoSpectral, unmix the test halves (plus the unstained AF file,
#' optionally gated), and score each off-target channel on the target's
#' positive events. A control is positive above the `pos_quantile` of its
#' universal negative's test half in the target channel, as named in the
#' AutoSpectral control file (`fluorophore`, `universal.negative`); negatives,
#' the `AF` row and the unstained file are scored on all events. Scores are on
#' the linear scale: `spread_q` (interval of `spread_probs` over its normal
#' width) and `bias` = |median|. With `repro = TRUE` the full controls are
#' extracted twice with the same seed and the largest absolute differences are
#' reported.
#'
#' Every AutoSpectral setting is passed by you; nothing is defaulted.
#'
#' @param control_dir Folder of single-stain control FCS files.
#' @param control_file AutoSpectral control file (csv) for `control_dir`.
#' @param unstained_file Unstained FCS file used for AF extraction.
#' @param out_dir New folder for the split, unmixed files and work files; must
#'   not exist.
#' @param cytometer Passed to `AutoSpectral::get.autospectral.param()`.
#' @param seed Seed for the split and for each extraction.
#' @param train_fraction Share of each control in the train half.
#' @param pos_quantile Quantile of the universal negative that marks
#'   positives (benchmark 05 used 0.995).
#' @param spread_probs Two probabilities for `spread_q`.
#' @param cell_gate `NULL` or a named list of two-number ranges (e.g.
#'   `` list(`FSC-A` = c(6e5, 4e6), `SSC-A` = c(1e5, 1.6e6)) ``) applied to the
#'   unstained file before it is unmixed and scored.
#' @param extract_args,af_args,variant_args,unmix_args Named lists of further
#'   arguments for `get.spectra.automated()`, `get.af.spectra()`,
#'   `get.spectral.variants()` and `unmix.folder()` (for example `refine`,
#'   `method`, `pipeline`, `speed`, `n.passes`). Use `list()` for none.
#' @param repro `TRUE` to run the same-seed reproducibility check.
#'
#' @return An `sw_heldout` list: `scores` (tibble: `control`, `type`,
#'   `target`, `channel`, `n`, `spread_q`, `bias`), `repro` (tibble or `NULL`),
#'   `split` (file, events, train, test) and `settings`.
#' @seealso [sw_unmix_event_metrics()], [sw_unmix_diagnose()]
#' @export
sw_unmix_heldout <- function(
  control_dir,
  control_file,
  unstained_file,
  out_dir,
  cytometer,
  seed,
  train_fraction,
  pos_quantile,
  spread_probs,
  cell_gate,
  extract_args,
  af_args,
  variant_args,
  unmix_args,
  repro
) {
  .sw_em_missing(
    "sw_unmix_heldout",
    environment(),
    c("control_dir", "control_file", "unstained_file", "out_dir", "cytometer", "seed", "train_fraction", "pos_quantile",
      "spread_probs", "cell_gate", "extract_args", "af_args", "variant_args", "unmix_args", "repro")
  )
  for (a in c("control_dir", "control_file", "unstained_file", "out_dir", "cytometer")) {
    v <- get(a)
    if (!is.character(v) || length(v) != 1L || is.na(v) || !nzchar(v)) {
      stop(sprintf("`%s` must be one non-empty string.", a), call. = FALSE)
    }
  }
  if (!dir.exists(control_dir)) stop("`control_dir` does not exist: ", control_dir, call. = FALSE)
  if (!file.exists(control_file)) stop("`control_file` does not exist: ", control_file, call. = FALSE)
  if (!file.exists(unstained_file)) stop("`unstained_file` does not exist: ", unstained_file, call. = FALSE)
  if (file.exists(out_dir)) stop("`out_dir` already exists: ", out_dir, ". Use a new folder.", call. = FALSE)
  .sw_em_scalar(seed, "seed", integer = TRUE)
  .sw_em_scalar(train_fraction, "train_fraction", lower = 0)
  if (train_fraction >= 1) stop("`train_fraction` must be below 1.", call. = FALSE)
  .sw_em_scalar(pos_quantile, "pos_quantile", lower = 0)
  if (pos_quantile >= 1) stop("`pos_quantile` must be below 1.", call. = FALSE)
  .sw_em_probs(spread_probs, "spread_probs")
  if (!is.null(cell_gate)) {
    ok <- is.list(cell_gate) && length(cell_gate) > 0L && !is.null(names(cell_gate)) &&
      all(vapply(cell_gate, \(r) is.numeric(r) && length(r) == 2L && !anyNA(r) && r[[1]] < r[[2]], logical(1)))
    if (!ok) stop("`cell_gate` must be NULL or a named list of increasing two-number ranges.", call. = FALSE)
  }
  for (a in c("extract_args", "af_args", "variant_args", "unmix_args")) {
    v <- get(a)
    if (!is.list(v) || (length(v) > 0L && (is.null(names(v)) || any(!nzchar(names(v)))))) {
      stop(sprintf("`%s` must be a named list (or list()).", a), call. = FALSE)
    }
  }
  if (!is.logical(repro) || length(repro) != 1L || is.na(repro)) stop("`repro` must be TRUE or FALSE.", call. = FALSE)
  if (!requireNamespace("AutoSpectral", quietly = TRUE) || !requireNamespace("flowCore", quietly = TRUE)) {
    stop("sw_unmix_heldout() needs AutoSpectral (>= 1.8.8) and flowCore.", call. = FALSE)
  }

  ctl <- utils::read.csv(control_file, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("filename", "fluorophore", "control.type", "universal.negative")
  if (!all(need %in% names(ctl))) {
    stop(sprintf("`control_file` needs columns %s.", paste(need, collapse = ", ")), call. = FALSE)
  }
  files <- sort(list.files(control_dir, pattern = "\\.fcs$"))
  absent <- setdiff(ctl$filename, files)
  if (length(absent) > 0L) stop(sprintf("Control files missing from `control_dir`: %s.", paste(absent, collapse = ", ")), call. = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dirs <- file.path(out_dir, c("train", "test", "unmixed"))
  for (d in dirs) dir.create(d, recursive = TRUE)
  read_raw <- \(f) flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)

  # Split every control; the train and test folders get a copy of the control file.
  set.seed(seed)
  split <- dplyr::bind_rows(lapply(ctl$filename, \(f) {
    ff <- read_raw(file.path(control_dir, f))
    n <- nrow(ff)
    tr <- sort(sample(n, floor(n * train_fraction)))
    flowCore::write.FCS(ff[tr, ], file.path(dirs[[1]], f))
    flowCore::write.FCS(ff[-tr, ], file.path(dirs[[2]], f))
    tibble::tibble(file = f, events = n, train = length(tr), test = n - length(tr))
  }))
  for (d in dirs[1:2]) {
    if (!file.copy(control_file, file.path(d, basename(control_file)))) stop("Could not copy the control file to ", d, call. = FALSE)
  }
  uf <- read_raw(unstained_file)
  if (!is.null(cell_gate)) uf <- flowCore::Subset(uf, do.call(flowCore::rectangleGate, cell_gate))
  un_name <- paste0("unstained__", basename(unstained_file))
  flowCore::write.FCS(uf, file.path(dirs[[2]], un_name))

  extract <- \(cdir, work) {
    dir.create(work, recursive = TRUE)
    asp <- .sw_ho_call("get.autospectral.param", list(cytometer = cytometer, figures = FALSE))
    df <- names(asp)[grepl("dir", names(asp)) & vapply(asp, is.character, logical(1))]
    asp[df] <- lapply(asp[df], \(d) file.path(work, d))
    cf <- file.path(cdir, basename(control_file))
    set.seed(seed)
    sp <- .sw_ho_call("get.spectra.automated", c(list(control.dir = cdir, control.def.file = cf, asp = asp), extract_args))
    af <- .sw_ho_call("get.af.spectra", c(list(unstained.sample = unstained_file, asp = asp, spectra = sp), af_args))
    va <- .sw_ho_call(
      "get.spectral.variants",
      c(list(control.dir = cdir, control.def.file = cf, asp = asp, spectra = sp, output.dir = work), variant_args)
    )
    list(spectra = sp, af = af, variants = va, asp = asp, cf = cf)
  }
  ex <- extract(dirs[[1]], file.path(out_dir, "work_train"))
  .sw_ho_call("unmix.folder", c(
    list(
      fcs.dir = dirs[[2]],
      spectra = ex$spectra,
      asp = ex$asp,
      flow.control = .sw_ho_call("reload.flow.control", list(dirs[[1]], ex$cf, ex$asp)),
      af.spectra = ex$af,
      spectra.variants = ex$variants,
      output.dir = dirs[[3]]
    ),
    unmix_args
  ))

  # Score off-target channels of every unmixed test file.
  um <- list.files(dirs[[3]], pattern = "\\.fcs$", full.names = TRUE)
  find_um <- \(f) {
    hit <- um[startsWith(basename(um), paste0(tools::file_path_sans_ext(f), " "))]
    if (length(hit) != 1L) stop(sprintf("Expected one unmixed file for %s, found %d.", f, length(hit)), call. = FALSE)
    flowCore::exprs(read_raw(hit))
  }
  fluors <- rownames(ex$spectra) # unmixed fluorophore channels
  negatives <- unique(ctl$universal.negative)
  neg_tbl <- stats::setNames(lapply(negatives, find_um), negatives)
  rows <- rbind(
    data.frame(control = ctl$filename, type = ctl$control.type, fluor = ctl$fluorophore, negative = ctl$universal.negative),
    data.frame(control = un_name, type = "unstained", fluor = "AF", negative = NA_character_)
  )
  zw <- .sw_z_width(spread_probs)
  scores <- dplyr::bind_rows(lapply(seq_len(nrow(rows)), \(r) {
    e <- find_um(rows$control[[r]])
    cols <- .sw_match_channels(fluors, colnames(e))
    cols <- cols[!is.na(cols)]
    is_neg <- rows$fluor[[r]] %in% c("AF", "Negative") || rows$control[[r]] %in% negatives
    target <- if (is_neg) NA_character_ else unname(.sw_match_channels(rows$fluor[[r]], colnames(e)))
    if (!is_neg && is.na(target)) stop(sprintf("No unmixed column for %s.", rows$fluor[[r]]), call. = FALSE)
    keep <- if (is_neg) {
      rep(TRUE, nrow(e))
    } else {
      e[, target] > stats::quantile(neg_tbl[[rows$negative[[r]]]][, target], pos_quantile, names = FALSE)
    }
    dplyr::bind_rows(lapply(setdiff(cols, target), \(j) {
      v <- e[keep, j]
      tibble::tibble(
        control = rows$control[[r]],
        type = rows$type[[r]],
        target = target,
        channel = j,
        n = sum(keep),
        spread_q = if (length(v) > 1L) unname(diff(stats::quantile(v, spread_probs))) / zw else NA_real_,
        bias = if (length(v) > 0L) abs(stats::median(v)) else NA_real_
      )
    }))
  }))

  rep_tbl <- NULL
  if (repro) {
    r1 <- extract(control_dir, file.path(out_dir, "work_rep1"))
    r2 <- extract(control_dir, file.path(out_dir, "work_rep2"))
    max_diff <- \(a, b) if (is.numeric(a) && identical(dim(a), dim(b))) max(abs(a - b)) else NA_real_
    rep_tbl <- tibble::tibble(
      object = c("spectra", "af_spectra", "variants"),
      identical = c(identical(r1$spectra, r2$spectra), identical(r1$af, r2$af), identical(r1$variants, r2$variants)),
      max_abs_diff = c(max_diff(r1$spectra, r2$spectra), max_diff(r1$af, r2$af), NA_real_)
    )
  }

  structure(
    list(
      scores = scores,
      repro = rep_tbl,
      split = split,
      settings = list(
        seed = seed,
        train_fraction = train_fraction,
        pos_quantile = pos_quantile,
        spread_probs = spread_probs,
        cell_gate = cell_gate,
        extract_args = extract_args,
        af_args = af_args,
        variant_args = variant_args,
        unmix_args = unmix_args,
        out_dir = out_dir,
        autospectral = as.character(utils::packageVersion("AutoSpectral"))
      )
    ),
    class = "sw_heldout"
  )
}

#' @export
print.sw_heldout <- function(x, ...) {
  cat("<sw_heldout>", length(unique(x$scores$control)), "control(s),", nrow(x$scores), "score rows\n")
  cat("  repro:", if (is.null(x$repro)) "not run" else paste(x$repro$object, x$repro$identical, collapse = "; "), "\n")
  invisible(x)
}


# ---- Diagnose and report hooks ----------------------------------------------------

# Check rows from an sw_event_metrics: event_channel and event_pair.
.sw_diag_event_metrics <- function(em) {
  if (is.null(em)) {
    return(NULL)
  }
  if (!inherits(em, "sw_event_metrics")) {
    stop("`event_metrics` must be an sw_event_metrics.", call. = FALSE)
  }
  sep <- .sw_unmix_item_sep
  st <- em$settings
  ch <- em$channels
  item <- paste(ch$file, ch$channel, sep = sep)
  metrics <- c("pct_below_cut", "pct_below_floor", "neg_median", "neg_rsd", "spread_q", "pileup", "hyper_share", "si_mad", "si_q")
  ch_rows <- dplyr::bind_rows(lapply(metrics, \(m) {
    reason <- sprintf("%d events, %d negative, %d positive; %s edges.", ch$n, ch$n_neg, ch$n_pos, st$edge_type)
    if (m %in% c("si_mad", "si_q")) {
      reason <- ifelse(ch$n_pos < st$min_pos, sprintf("Fewer than %d positives (%d).", st$min_pos, ch$n_pos), reason)
    }
    if (m == "pct_below_cut") {
      reason <- sprintf("asinh(x / %s) < %s.", format(st$cofactor), format(st$negative_cut))
    }
    .sw_check_rows("event_channel", "channel", item, m, ch[[m]], reason = reason)
  }))
  miss_rows <- if (nrow(em$missing) > 0L) {
    .sw_check_rows(
      "event_channel",
      "channel",
      paste(em$missing$file, em$missing$channel, sep = sep),
      "pct_below_cut",
      NA,
      "No column matches this edge channel.",
      run = FALSE
    )
  }
  pair_rows <- if (!is.null(em$pairs) && nrow(em$pairs) > 0L) {
    pr <- em$pairs
    pitem <- paste(pr$file, pr$channel_i, pr$channel_j, sep = sep)
    reason <- ifelse(
      pr$n < st$min_cells,
      sprintf("%d i+ j- events, below min_cells %d.", pr$n, st$min_cells),
      sprintf("%d i+ j- events; %d bins.", pr$n, st$n_bins)
    )
    dplyr::bind_rows(
      .sw_check_rows("event_pair", "pair", pitem, "abs_slope", abs(pr$slope), reason = reason),
      .sw_check_rows("event_pair", "pair", pitem, "abs_rho", abs(pr$rho), reason = reason)
    )
  }
  dplyr::bind_rows(ch_rows, miss_rows, pair_rows)
}

# Check rows from an sw_heldout.
.sw_diag_heldout <- function(ho) {
  if (is.null(ho)) {
    return(NULL)
  }
  if (!inherits(ho, "sw_heldout")) {
    stop("`heldout` must be an sw_heldout.", call. = FALSE)
  }
  sc <- ho$scores
  item <- paste(sc$control, sc$channel, sep = .sw_unmix_item_sep)
  reason <- sprintf(
    "%s control; %d scored events; target %s.",
    sc$type,
    sc$n,
    dplyr::coalesce(sc$target, "none (all events)")
  )
  dplyr::bind_rows(
    .sw_check_rows("heldout", "control", item, "heldout_spread_q", sc$spread_q, reason = reason),
    .sw_check_rows("heldout", "control", item, "heldout_bias", sc$bias, reason = reason)
  )
}

# HTML sections of sw_unmix_report() for event metrics and held-out scores.
.sw_report_event_sections <- function(em, ho) {
  med <- \(x) stats::median(x, na.rm = TRUE)
  c(
    if (!is.null(em)) {
      st <- em$settings
      by_ch <- em$channels |>
        dplyr::summarise(
          n_files = dplyr::n(),
          pct_below_cut = med(.data$pct_below_cut),
          pct_below_floor = med(.data$pct_below_floor),
          spread_q = med(.data$spread_q),
          pileup = med(.data$pileup),
          hyper_share = med(.data$hyper_share),
          si_mad = med(.data$si_mad),
          si_q = med(.data$si_q),
          .by = "channel"
        ) |>
        dplyr::arrange(dplyr::desc(.data$pct_below_cut))
      c(
        "<h2>Event metrics</h2>",
        sprintf(
          "<p>Median over files per channel, sorted by %% below the cut. Cut asinh(x / %s) &lt; %s; %s edges (%s); spread_q over the %s / %s interval; hypernegative below median - %s rSD.</p>",
          .sw_html_escape(format(st$cofactor)),
          .sw_html_escape(format(st$negative_cut)),
          .sw_html_escape(st$edge_type),
          .sw_html_escape(paste(st$edge_source, collapse = ", ")),
          .sw_html_escape(format(st$spread_probs[[1]])),
          .sw_html_escape(format(st$spread_probs[[2]])),
          .sw_html_escape(format(st$hyper_k))
        ),
        .sw_html_table(by_ch),
        if (!is.null(em$pairs)) {
          top <- em$pairs[!is.na(em$pairs$slope), ]
          top <- utils::head(top[order(-abs(top$slope)), ], 25L)
          c(
            sprintf(
              "<p>Line check: the 25 file x pair rows with the largest |slope| (of %d; all in event_pairs.csv). Rows below %d i+ j- events have no slope.</p>",
              nrow(em$pairs),
              st$min_cells
            ),
            .sw_html_table(top)
          )
        }
      )
    },
    if (!is.null(ho)) {
      by_type <- ho$scores |>
        dplyr::summarise(
          n_controls = dplyr::n_distinct(.data$control),
          spread_q = med(.data$spread_q),
          bias = med(.data$bias),
          .by = "type"
        )
      c(
        "<h2>Held-out control error</h2>",
        sprintf(
          "<p>Median over off-target channels by control type (linear). Train fraction %s, seed %s, positives above q%s of the universal negative.</p>",
          .sw_html_escape(format(ho$settings$train_fraction)),
          .sw_html_escape(format(ho$settings$seed)),
          .sw_html_escape(format(ho$settings$pos_quantile))
        ),
        .sw_html_table(by_type),
        if (!is.null(ho$repro)) c("<p>Same-seed reproducibility:</p>", .sw_html_table(ho$repro))
      )
    }
  )
}
