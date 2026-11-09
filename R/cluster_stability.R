# R/cluster_stability.R
# How much a k-means clustering of the pooled events changes with the seed and
# with the number of events. The check comes from O. Burton's Colibri Cytometry
# post "The peculiarities of Phenograph" (cluster labels moved for 40% of cells
# between two seeds); the code is our own and uses stats::kmeans, as
# sw_cluster_scores() does.

#' Check how stable a clustering is across seeds and event counts
#'
#' @description
#' Pool `n_events` events from each file, then:
#'
#' * cluster the pool with every seed in `seeds` and compare each pair of runs
#'   (`ari_seed_pair`, `moved_share_seed_pair`);
#' * for each size in `sizes`, cluster a subsample of that size with the first
#'   seed, assign every pooled event to the nearest centre, and compare with
#'   the full-pool run of the first seed (`ari_size`, `moved_share_size`).
#'
#' `ari` is the adjusted Rand index (1 = same partition, about 0 = chance).
#' `moved_share` is the share of events whose cluster differs after the best
#' one-to-one matching of cluster labels (greedy on the contingency table).
#' Numbers only; nothing is marked unstable.
#'
#' @details
#' **No defaults.** `k`, `seeds`, `sizes`, `n_events`, `draw_seed` and `nstart`
#' change the numbers, so each must be passed. The caller's RNG state is
#' restored.
#'
#' **Scale is yours.** Cluster on the scale you analyse; this function never
#' transforms or scales.
#'
#' @param data A named list of numeric matrices, one per file.
#' @param markers Character vector of columns to cluster on.
#' @param k Whole number of clusters.
#' @param seeds Whole numbers, at least two, one k-means run each.
#' @param sizes Whole numbers, subsample sizes (each below the pool size).
#' @param n_events Whole number of events drawn per file into the pool.
#' @param draw_seed Whole number for the event draws.
#' @param nstart Whole number of k-means random starts per run.
#'
#' @return A tibble with columns `check`, `item` (the seed pair `"s1 ~ s2"` or
#'   the size), and `value`. Attribute `pool_size` holds the pool's event count.
#'
#' @seealso [sw_cluster_scores()]
#'
#' @export
sw_cluster_stability <- function(data, markers, k, seeds, sizes, n_events, draw_seed, nstart) {
  if (!is.list(data) || is.null(names(data))) {
    stop("`data` must be a named list of matrices.", call. = FALSE)
  }
  .sw_check_count(k, "k")
  .sw_check_count(n_events, "n_events")
  .sw_check_count(draw_seed, "draw_seed")
  .sw_check_count(nstart, "nstart")
  if (length(seeds) < 2) {
    stop("`seeds` must hold at least two seeds.", call. = FALSE)
  }
  for (s in seeds) .sw_check_count(s, "seeds")
  for (s in sizes) .sw_check_count(s, "sizes")
  for (f in names(data)) {
    absent <- setdiff(markers, colnames(data[[f]]))
    if (length(absent) > 0) {
      stop("File ", f, " lacks markers: ", paste(absent, collapse = ", "), call. = FALSE)
    }
    if (nrow(data[[f]]) < n_events) {
      stop("File ", f, " has fewer than `n_events` events.", call. = FALSE)
    }
  }
  pool <- .sw_with_seed(draw_seed, {
    do.call(
      rbind,
      lapply(data, \(m) as.matrix(m[sample.int(nrow(m), n_events), markers, drop = FALSE]))
    )
  })
  if (any(sizes >= nrow(pool))) {
    stop("Every value in `sizes` must be below the pool size (", nrow(pool), ").", call. = FALSE)
  }

  run <- function(x, seed) {
    .sw_with_seed(seed, suppressWarnings(stats::kmeans(x, centers = k, nstart = nstart)))
  }
  full <- lapply(seeds, \(s) run(pool, s)$cluster)
  pairs <- utils::combn(seq_along(seeds), 2)
  seed_rows <- lapply(seq_len(ncol(pairs)), function(j) {
    a <- full[[pairs[1, j]]]
    b <- full[[pairs[2, j]]]
    item <- paste(seeds[pairs[1, j]], seeds[pairs[2, j]], sep = " ~ ")
    tibble::tibble(
      check = c("ari_seed_pair", "moved_share_seed_pair"),
      item = item,
      value = c(.sw_ari(a, b), .sw_moved_share(a, b))
    )
  })
  size_rows <- lapply(sizes, function(n) {
    idx <- .sw_with_seed(seeds[[1]], sample.int(nrow(pool), n))
    fit <- run(pool[idx, , drop = FALSE], seeds[[1]])
    # Nearest centre for every pooled event.
    d2 <- outer(rowSums(pool^2), rowSums(fit$centers^2), "+") - 2 * tcrossprod(pool, fit$centers)
    assigned <- max.col(-d2, ties.method = "first")
    tibble::tibble(
      check = c("ari_size", "moved_share_size"),
      item = as.character(n),
      value = c(.sw_ari(full[[1]], assigned), .sw_moved_share(full[[1]], assigned))
    )
  })
  out <- dplyr::bind_rows(c(seed_rows, size_rows))
  attr(out, "pool_size") <- nrow(pool)
  out
}

# Adjusted Rand index from the contingency table (Hubert and Arabie 1985).
.sw_ari <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(n) n * (n - 1) / 2
  sum_ij <- sum(comb2(tab))
  sum_a <- sum(comb2(rowSums(tab)))
  sum_b <- sum(comb2(colSums(tab)))
  expected <- sum_a * sum_b / comb2(sum(tab))
  max_index <- (sum_a + sum_b) / 2
  if (max_index == expected) {
    return(1)
  }
  (sum_ij - expected) / (max_index - expected)
}

# Share of events whose label differs after greedy one-to-one label matching.
.sw_moved_share <- function(a, b) {
  tab <- unclass(table(a, b))
  matched <- 0
  while (length(tab) > 0 && max(tab) > 0) {
    hit <- which(tab == max(tab), arr.ind = TRUE)[1, ]
    matched <- matched + tab[hit[1], hit[2]]
    tab <- tab[-hit[1], -hit[2], drop = FALSE]
  }
  1 - matched / length(a)
}
