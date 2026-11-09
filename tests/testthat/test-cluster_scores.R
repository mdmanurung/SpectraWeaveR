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
    min_cluster_events = 30
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
      nstart = 2
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
  res <- cs(d, min_cluster_events = 500)
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
      min_cluster_events = 50
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
