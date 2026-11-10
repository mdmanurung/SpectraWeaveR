# tests/testthat/test-cluster_scores.R
# Unit tests for R/cluster_scores.R on three simulated populations.

# Populations at -4, 0 and 4 (marker A); shares 0.5, 0.3, 0.2. 14 files.
make_pops <- function(n = 600, seed = 1) {
  set.seed(seed)
  f <- lapply(1:14, function(i) {
    pop <- sample(1:3, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
    cbind(
      A = c(-4, 0, 4)[pop] + stats::rnorm(n),
      B = stats::rnorm(n)
    )
  })
  names(f) <- paste0("f", 1:14)
  f
}

cs <- function(d, ...) {
  args <- list(
    data = d,
    markers = c("A", "B"),
    k = 3,
    n_events = 500,
    seed = 1,
    nstart = 5,
    min_cluster_events = 30,
    reference = NULL,
    reference_group = NULL
  )
  args[names(list(...))] <- list(...)
  do.call(sw_cluster_scores, args)
}

test_that("required arguments have no default", {
  expect_error(sw_cluster_scores(), "no default")
  d <- make_pops()
  expect_error(
    sw_cluster_scores(
      d,
      c("A", "B"),
      k = 3,
      n_events = 100,
      seed = 1,
      nstart = 2, reference_group = NULL
    ),
    "min_cluster_events"
  )
})

test_that("a shift confined to the rarest population shows in that cluster only", {
  d <- make_pops()
  # Shift marker B by 1 in population 3 (A near 4) for file f5 only.
  hit <- d$f5[, "A"] > 2
  d$f5[hit, "B"] <- d$f5[hit, "B"] + 1
  res <- cs(d)
  expect_s3_class(res, "sw_cluster_scores")
  m <- res[res$check == "median_z" & res$marker == "B", ]
  # The cluster whose centre has the highest A is population 3.
  pop3 <- which.max(attr(res, "centers")[, "A"])
  z5 <- m$value[m$file == "f5" & m$cluster == pop3]
  expect_gt(z5, 4)
  # Whole-file median of B barely moves: shift 1 on 20% of events.
  expect_lt(abs(stats::median(d$f5[, "B"])), 0.3)
  # Other clusters of the same file are not flagged.
  others <- m$value[m$file == "f5" & m$cluster != pop3]
  expect_true(all(abs(others) < 3, na.rm = TRUE))
})

test_that("a composition change shows in freq_z", {
  d <- make_pops()
  # Replace half of the population-1 events of f9 by population-3 events.
  p1 <- which(d$f9[, "A"] < -2)
  swap <- p1[seq_len(floor(length(p1) / 2))]
  d$f9[swap, "A"] <- d$f9[swap, "A"] + 8
  res <- cs(d)
  fz <- res[res$check == "freq_z", ]
  big <- fz[fz$file == "f9", ]
  expect_gt(max(abs(big$value)), 4)
  expect_true(is.na(unique(big$marker)))
})

test_that("clean files give moderate z-scores and share sums to one", {
  res <- cs(make_pops())
  expect_lt(max(abs(res$value), na.rm = TRUE), 6)
  s <- res[res$check == "freq_z", ]
  expect_equal(sum(s$share[s$file == "f1"]), 1)
})

test_that("small clusters get NA medians and input errors are clear", {
  d <- make_pops()
  res <- cs(d, min_cluster_events = 500, reference = NULL)
  expect_true(all(is.na(res$value[res$check == "median_z"])))
  expect_error(cs(d, k = 1), "at least 2")
  expect_error(cs(d, n_events = 5000), "fewer than `n_events`")
  expect_error(cs(d, markers = "Z"), "lacks marker")
  expect_error(cs(d[1:2]), "at least three files")
})

test_that("results do not depend on file order and keep the RNG state", {
  d <- make_pops()
  set.seed(8)
  a <- stats::runif(1)
  set.seed(8)
  r1 <- cs(d)
  expect_equal(stats::runif(1), a)
  r2 <- cs(d[rev(names(d))])
  expect_equal(r1$value, r2$value)
})

test_that("kmeans warnings are kept in an attribute and re-issued once", {
  res <- cs(make_pops())
  expect_type(attr(res, "kmeans_warnings"), "character")
  # A large pool of overlapping clusters can trigger the quick-transfer warning;
  # whether it fires depends on the data, so only the contract is tested: if
  # any warning was kept, it was also signalled as one clear message.
  set.seed(3)
  big <- lapply(1:12, function(i) {
    cbind(A = stats::rnorm(20000), B = stats::rnorm(20000))
  })
  names(big) <- paste0("g", 1:12)
  got <- NULL
  out <- withCallingHandlers(
    sw_cluster_scores(
      big,
      c("A", "B"),
      k = 30,
      n_events = 20000,
      seed = 1,
      nstart = 1,
      min_cluster_events = 50, reference = NULL, reference_group = NULL
    ),
    warning = function(w) {
      got <<- c(got, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  kept <- attr(out, "kmeans_warnings")
  expect_equal(length(got) > 0, length(kept) > 0)
  if (length(kept) > 0) {
    expect_length(got, 1)
    expect_match(got, "k-means reported")
  }
})

test_that("reference fits clusters on the given files and assigns the rest", {
  d <- make_pops()
  base <- cs(d)
  # All files as reference gives the same clusters as NULL (converged k-means
  # labels are the nearest centres).
  all_ref <- cs(d, reference = names(d))
  expect_equal(all_ref$value, base$value)
  # A reference subset still scores every file.
  sub <- cs(d, reference = paste0("f", 1:7))
  expect_setequal(unique(sub$file), names(d))
  expect_equal(nrow(sub), nrow(base))
  # A file whose rare population moved is still found with clusters fitted on
  # clean files only.
  d$f12[d$f12[, "A"] > 2, "B"] <- d$f12[d$f12[, "A"] > 2, "B"] + 1.5
  r <- cs(d, reference = paste0("f", 1:7))
  mz <- r[r$check == "median_z" & r$marker == "B", ]
  expect_equal(mz$file[which.max(abs(mz$value))], "f12")
  expect_error(cs(d, reference = "nope"), "reference")
  expect_error(cs(d, reference = c("f1", "f1")), "reference")
})

# reference_group ---------------------------------------------------------------

# Three donors of four files, two populations. The share of the right-hand
# population is 0.2, 0.5 and 0.8 by donor; file d2_4 has 0.35 instead of 0.5.
make_donor_shares <- function(seed = 4, n = 500) {
  set.seed(seed)
  ids <- paste0("d", rep(1:3, each = 4), "_", rep(1:4, 3))
  d <- lapply(ids, function(i) {
    donor <- as.integer(substr(i, 2, 2))
    p <- c(0.2, 0.5, 0.8)[donor]
    if (i == "d2_4") p <- 0.35
    right <- stats::runif(n) < p
    cbind(A = stats::rnorm(n, ifelse(right, 4, -4)), B = stats::rnorm(n))
  })
  names(d) <- ids
  list(data = d, donor = stats::setNames(substr(ids, 1, 2), ids))
}

test_that("reference_group finds a composition change that is unusual only within its donor", {
  x <- make_donor_shares()
  top_z <- function(res) {
    r <- res[res$check == "freq_z" & res$file == "d2_4", ]
    max(abs(r$value), na.rm = TRUE)
  }
  all_files <- cs(x$data, k = 2)
  by_donor <- cs(x$data, k = 2, reference_group = x$donor)
  expect_lt(top_z(all_files), 2)
  expect_gt(top_z(by_donor), 3)
  f <- by_donor[by_donor$check == "freq_z", ]
  expect_equal(f$file[which.max(abs(f$value))], "d2_4")
  # A group of 2 files gives NA for every row of those files.
  g <- x$donor
  g[names(g) %in% c("d3_1", "d3_2")] <- "other"
  small <- cs(x$data, k = 2, reference_group = g)
  expect_true(all(is.na(small$value[small$file %in% c("d3_1", "d3_2")])))
})

test_that("cluster reference_group is required and must cover every file", {
  x <- make_donor_shares()
  expect_error(
    sw_cluster_scores(x$data, c("A", "B"), 2, 100, 1, 3, 10, reference = NULL),
    "no default"
  )
  expect_error(cs(x$data, k = 2, reference_group = x$donor[-1]), "label for every file")
  expect_error(cs(x$data, k = 2, reference_group = unname(x$donor)), "label for every file")
})
