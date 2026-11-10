make_sample <- function(shift = 0, width = 1, n = 3000, seed = 1) {
  set.seed(seed)
  x <- c(stats::rnorm(round(n * 0.6), 0, 1), stats::rnorm(round(n * 0.4), 6, 1))
  # scale about 0, then shift
  cbind(A = x * width + shift, B = stats::rnorm(length(x)))
}

ref_gate <- function() sw_gate_range("pop", "A", -2, 2, parent = NA)

run <- function(data, gates = list(ref_gate()), ref = make_sample(seed = 10),
                bandwidth = 0.3, n_landmarks = 1, allow_scale = FALSE, max_shift = 4) {
  sw_gate_transfer(gates, ref, data, bandwidth, n_landmarks, allow_scale, max_shift)
}

test_that("a known shift is recovered and kept_share matches the reference", {
  ref <- make_sample(seed = 10)
  smp <- make_sample(shift = 1.5, seed = 2)
  res <- run(list(s = smp), ref = ref)
  q <- res$qc
  expect_equal(q$status, "adapted")
  expect_equal(q$shift, 1.5, tolerance = 0.15)
  expect_equal(q$scale, 1)
  g <- res$gates$s$pop
  expect_equal(unname(g$vertices["min", 1]), -2 + q$shift)
  ks <- function(d, gate) mean(sw_gate_apply(gate, d))
  after <- ks(smp, g)
  before <- ks(smp, ref_gate())
  target <- ks(ref, ref_gate())
  expect_lt(abs(after - target), 0.03)
  expect_gt(abs(before - target), 0.1)
})

test_that("max_shift refuses and returns the reference gate unchanged", {
  res <- run(list(s = make_sample(shift = 1.5, seed = 2)), max_shift = 0.5)
  expect_equal(res$qc$status, "refused_max_shift")
  expect_identical(res$gates$s$pop, ref_gate())
})

test_that("a missing landmark is refused", {
  set.seed(3)
  uni <- cbind(A = stats::rnorm(2000), B = stats::rnorm(2000))
  res <- run(list(s = uni), n_landmarks = 2, bandwidth = 1.2)
  expect_equal(res$qc$status, "refused_no_landmark")
  expect_equal(res$qc$n_landmarks_found, 1L)
  expect_true(is.na(res$qc$shift))
  expect_identical(res$gates$s$pop, ref_gate())
  few <- run(list(s = cbind(A = 1, B = 1)))
  expect_equal(few$qc$status, "refused_no_landmark")
})

test_that("allow_scale recovers a known width change", {
  ref <- make_sample(seed = 10)
  smp <- make_sample(shift = 0.5, width = 1.25, seed = 4)
  res <- run(list(s = smp), ref = ref, n_landmarks = 2, allow_scale = TRUE, bandwidth = 0.4)
  q <- res$qc
  expect_equal(q$status, "adapted")
  expect_equal(q$scale, 1.25, tolerance = 0.05)
  # peaks at 0 and 6 move to 0.5 and 8, so s1 - r1 = 0.5
  expect_lt(abs(q$shift - 0.5), 0.25)
  g <- res$gates$s$pop
  after <- mean(sw_gate_apply(g, smp))
  target <- mean(sw_gate_apply(ref_gate(), ref))
  expect_lt(abs(after - target), 0.03)
  # a shift-only fit on the same data is worse
  fixed <- run(list(s = smp), ref = ref, n_landmarks = 2, allow_scale = FALSE, bandwidth = 0.4)
  worse <- mean(sw_gate_apply(fixed$gates$s$pop, smp))
  expect_gt(abs(worse - target), abs(after - target))
})

test_that("hierarchical gates: child follows the adapted parent", {
  ref <- make_sample(seed = 10)
  smp <- make_sample(shift = 1.5, seed = 2)
  parent <- sw_gate_range("lo", "A", -3, 3, parent = NA)
  child <- sw_gate_polygon("kid", c("A", "B"),
    rbind(c(-2, -2), c(2, -2), c(2, 2), c(-2, 2)), parent = "lo")
  res <- run(list(s = smp), gates = list(child, parent), ref = ref)
  expect_equal(res$qc$gate, c("lo", "kid", "kid"))
  expect_equal(res$qc$status, rep("adapted", 3))
  expect_equal(res$qc$shift[res$qc$gate == "lo"], 1.5, tolerance = 0.15)
  kid <- res$gates$s$kid
  expect_equal(
    unname(kid$vertices[, 1]),
    c(-2, 2, 2, -2) + res$qc$shift[res$qc$gate == "kid" & res$qc$channel == "A"]
  )
  expect_lt(abs(res$qc$shift[res$qc$channel == "B"]), 0.3)
})

test_that("per-channel bandwidth and max_shift are accepted", {
  smp <- make_sample(shift = 1.5, seed = 2)
  g <- sw_gate_polygon("kid", c("A", "B"), rbind(c(-2, -2), c(2, -2), c(2, 2), c(-2, 2)), parent = NA)
  res <- run(list(s = smp), gates = list(g),
    bandwidth = c(A = 0.3, B = 0.3), max_shift = c(A = 4, B = 0.001))
  expect_equal(res$qc$status, rep("refused_max_shift", 2))
})

test_that("sw_gate_transfer checks arguments", {
  smp <- make_sample(seed = 2)
  expect_error(run(list(s = smp), n_landmarks = 3), "n_landmarks")
  expect_error(run(list(s = smp), allow_scale = NA), "allow_scale")
  expect_error(run(list(s = smp), bandwidth = -1), "positive")
  expect_error(sw_gate_transfer(list(ref_gate()), make_sample(), list(s = smp)), "bandwidth")
})
