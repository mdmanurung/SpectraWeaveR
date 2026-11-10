set.seed(11)
bimodal <- function(shift = 0, n = 2000) {
  c(rnorm(n, 1 + shift, 0.4), rnorm(n, 5 + shift, 0.5))
}

test_that("valley cut sits between the two modes", {
  x <- bimodal()
  cut <- sw_gate_density_cut(x, "valley", bandwidth = 0.3,
    min_peak_height = 0.1, min_peak_distance = 1)
  expect_false(is.na(cut))
  expect_true(is.na(attr(cut, "reason")))
  expect_gt(cut, 2)
  expect_lt(cut, 4)
  expect_length(attr(cut, "peaks"), 2)
})

test_that("valley returns NA with a reason for one peak", {
  x <- rnorm(3000)
  cut <- sw_gate_density_cut(x, "valley", bandwidth = 0.3,
    min_peak_height = 0.1, min_peak_distance = 1)
  expect_true(is.na(cut))
  expect_match(attr(cut, "reason"), "valley needs two")
})

test_that("min_peak_height drops a tiny second peak", {
  x <- c(rnorm(5000, 0, 0.5), rnorm(60, 6, 0.3))
  low <- sw_gate_density_cut(x, "valley", bandwidth = 0.3,
    min_peak_height = 0.01, min_peak_distance = 1)
  high <- sw_gate_density_cut(x, "valley", bandwidth = 0.3,
    min_peak_height = 0.5, min_peak_distance = 1)
  expect_false(is.na(low))
  expect_true(is.na(high))
})

test_that("percentile and mad cuts are exact", {
  x <- c(1:9, 100)
  p <- sw_gate_density_cut(x, "percentile", prob = 0.5)
  expect_equal(as.numeric(p), 5.5)
  expect_true(is.na(attr(p, "reason")))
  m <- sw_gate_density_cut(x, "mad", k = 2, mad_side = "upper")
  expect_equal(as.numeric(m), median(x) + 2 * mad(x))
  lo <- sw_gate_density_cut(x, "mad", k = 2, mad_side = "lower")
  expect_equal(as.numeric(lo), median(x) - 2 * mad(x))
  z <- sw_gate_density_cut(rep(3, 10), "mad", k = 2, mad_side = "upper")
  expect_true(is.na(z))
  expect_match(attr(z, "reason"), "MAD")
})

test_that("required arguments have no default", {
  x <- rnorm(100)
  expect_error(sw_gate_density_cut(x, "valley", bandwidth = 0.3), "min_peak_height")
  expect_error(sw_gate_density_cut(x, "valley", min_peak_height = 0.1,
    min_peak_distance = 1), "bandwidth")
  expect_error(sw_gate_density_cut(x, "percentile"), "prob")
  expect_error(sw_gate_density_cut(x, "mad", k = 2), "mad_side")
  expect_error(sw_gate_density_cut(x), "method")
  expect_error(sw_gate_density_range(data.frame(a = x), "a", "g", NA,
    method = "percentile", prob = 0.5), "side")
})

test_that("range gate keeps the side of the cut", {
  x <- bimodal()
  d <- cbind(CD3 = x)
  g <- sw_gate_density_range(d, "CD3", "cd3", NA, side = "above", method = "valley",
    bandwidth = 0.3, min_peak_height = 0.1, min_peak_distance = 1)
  expect_s3_class(g, "sw_gate")
  inside <- sw_gate_apply(g, d)
  expect_gt(mean(inside), 0.45)
  expect_lt(mean(inside), 0.55)
  b <- sw_gate_density_range(d, "CD3", "neg", NA, side = "below",
    method = "percentile", prob = 0.25)
  expect_equal(mean(sw_gate_apply(b, d)), 0.25, tolerance = 0.01)
  f <- sw_gate_density_range(cbind(CD3 = rnorm(500)), "CD3", "cd3", NA, side = "above",
    method = "valley", bandwidth = 0.3, min_peak_height = 0.1, min_peak_distance = 1)
  expect_s3_class(f, "sw_gate_failed")
  expect_true(nzchar(f$reason))
})

test_that("2D contour gate keeps the blob and excludes a far cluster", {
  blob <- cbind(FSC = rnorm(3000, 100, 10), SSC = rnorm(3000, 50, 8))
  far <- cbind(FSC = rnorm(300, 200, 5), SSC = rnorm(300, 150, 5))
  d <- rbind(blob, far)
  g <- sw_gate_density_2d(d, c("FSC", "SSC"), "lymph", NA, level = 0.9,
    bandwidth = c(4, 3))
  expect_s3_class(g, "sw_gate")
  inside <- sw_gate_apply(g, d)
  expect_gt(mean(inside[1:3000]), 0.8)
  expect_equal(sum(inside[3001:3300]), 0)
  expect_error(sw_gate_density_2d(d, c("FSC", "SSC"), "l", NA, bandwidth = c(4, 3)),
    "level")
  tiny <- sw_gate_density_2d(d[1:2, ], c("FSC", "SSC"), "l", NA, level = 0.9,
    bandwidth = c(4, 3))
  expect_s3_class(tiny, "sw_gate_failed")
})

test_that("per-sample application gives own cuts and reports failures", {
  mk <- function(shift) {
    cbind(CD3 = bimodal(shift), FSC = rnorm(4000, 100, 10), SSC = rnorm(4000, 50, 8))
  }
  samples <- list(
    a = mk(0), b = mk(2),
    c = cbind(CD3 = rnorm(4000), FSC = rnorm(4000, 100, 10), SSC = rnorm(4000, 50, 8))
  )
  vopts <- list(side = "above", method = "valley", bandwidth = 0.3,
    min_peak_height = 0.1, min_peak_distance = 1)
  specs <- list(
    lymph = sw_gate_density_spec(sw_gate_density_2d, channels = c("FSC", "SSC"),
      id = "lymph", parent = NA, level = 0.99, bandwidth = c(4, 3)),
    cd3 = do.call(sw_gate_density_spec, c(list(fun = sw_gate_density_range,
      channel = "CD3", id = "cd3", parent = "lymph"), vopts)),
    orphan = do.call(sw_gate_density_spec, c(list(fun = sw_gate_density_range,
      channel = "CD3", id = "orphan", parent = "missing"), vopts))
  )
  res <- sw_gate_density_apply(samples, specs)
  tab <- attr(res, "cuts")
  expect_named(res, c("a", "b", "c"))
  cuts <- tab[tab$gate == "cd3", ]
  ca <- cuts$cut[cuts$sample == "a"]
  cb <- cuts$cut[cuts$sample == "b"]
  expect_gt(cb - ca, 1.5)
  expect_equal(cuts$status[cuts$sample == "c"], "failed")
  expect_match(cuts$reason[cuts$sample == "c"], "valley needs two")
  expect_false("cd3" %in% names(res$c))
  expect_true(all(tab$status[tab$gate == "orphan"] == "failed"))
  expect_match(tab$reason[tab$gate == "orphan"][1], "not available")
  expect_error(sw_gate_density_apply(samples, list(x = specs$cd3)), "equal to its name")
})
