# Square gate on a 10 x 10 grid of points: hand-checkable numbers.
grid_events <- function() {
  g <- expand.grid(A = seq(0.5, 9.5, by = 1), B = seq(0.5, 9.5, by = 1))
  as.matrix(g)
}

test_that("sw_gate_qc gives hand-checked numbers for a square gate", {
  m <- grid_events()
  gate <- sw_gate_polygon("sq", c("A", "B"), rbind(c(2, 2), c(8, 2), c(8, 8), c(2, 8)), parent = NA)
  out <- sw_gate_qc(list(s1 = m), list(gate), reference = NULL, edge_width = 0.6)
  val <- function(chk) out$value[out$check == chk]
  expect_named(out, c("sample", "gate", "check", "channel", "value"))
  expect_equal(val("n_parent"), 100)
  # inside: coordinates 2.5..7.5 on both axes -> 6 x 6 points
  expect_equal(val("kept_share"), 36 / 100)
  # Distance to the outline of [2,8]^2 is <= 0.6 for a point when one coordinate
  # is within 0.6 of 2 or 8 (values 1.5, 2.5, 7.5, 8.5) and the other coordinate
  # lies in 1.5..8.5. Count: 16 (x in band, y in 3.5..6.5) + 16 (swapped) + 12 (both in band, minus the 4 outer corners at 0.71).
  expect_equal(val("edge_share"), 44 / 100)
  expect_true(is.na(val("kept_share_z")))
})

test_that("range gate edge distance uses the finite bounds", {
  m <- cbind(A = c(0.5, 1.5, 2.5, 3.5, 4.5))
  g <- sw_gate_range("r", "A", 1, Inf, parent = NA)
  out <- sw_gate_qc(list(s = m), list(g), reference = NULL, edge_width = 0.6)
  expect_equal(out$value[out$check == "kept_share"], 4 / 5)
  # 0.5 and 1.5 are within 0.6 of the bound 1; the open side has no edge.
  expect_equal(out$value[out$check == "edge_share"], 2 / 5)
})

test_that("centroid shift against a reference", {
  mk <- function(k) {
    m <- grid_events()
    m[, "A"] <- m[, "A"] + k
    m
  }
  data <- list(r = mk(0), a = mk(0), b = mk(0), c = mk(1))
  gate <- sw_gate_polygon("sq", c("A", "B"), rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10)), parent = NA)
  out <- sw_gate_qc(data, list(gate), reference = "r", edge_width = 0.1)
  ks <- out[out$check == "kept_share", ]
  expect_equal(ks$value[ks$sample == "c"], 0.9)
  z <- out[out$check == "kept_share_z", ]
  # kept shares 1, 1, 1, 0.9 have MAD 0, so z is NA by rule
  expect_true(all(is.na(z$value)))
  cs <- out[out$check == "centroid_shift" & out$sample == "c", ]
  expect_equal(cs$channel, c("A", "B"))
  expect_equal(cs$value[cs$channel == "B"], 0)
  expect_gt(cs$value[cs$channel == "A"], 0)
  cd <- out[out$check == "centroid_dist" & out$sample == "c", ]
  expect_true(is.na(cd$channel))
  expect_equal(cd$value, abs(cs$value[cs$channel == "A"]))
  expect_equal(out$value[out$check == "centroid_shift" & out$sample == "r"], c(0, 0))
})

test_that("robust z uses the MAD formula", {
  mk <- function(n_out) {
    m <- cbind(A = seq(0.5, 9.5, by = 1))
    m[seq_len(n_out), 1] <- 20
    m
  }
  data <- list(a = mk(0), b = mk(1), c = mk(2), d = mk(3), e = mk(6))
  g <- sw_gate_range("r", "A", 0, 10, parent = NA)
  out <- sw_gate_qc(data, list(g), reference = NULL, edge_width = 0.1)
  v <- c(1, 0.9, 0.8, 0.7, 0.4)
  expected <- (v - median(v)) / (1.4826 * median(abs(v - median(v))))
  expect_equal(out$value[out$check == "kept_share_z"], expected)
})

test_that("gates apply hierarchically through the parent chain", {
  m <- grid_events()
  parent <- sw_gate_range("p", "A", 0, 5, parent = NA)
  child <- sw_gate_range("c", "B", 0, 5, parent = "p")
  # child listed first on purpose: order must not matter
  out <- sw_gate_qc(list(s = m), list(child, parent), reference = NULL, edge_width = 0.1)
  expect_equal(out$value[out$gate == "c" & out$check == "n_parent"], 50)
  expect_equal(out$value[out$gate == "c" & out$check == "kept_share"], 0.5)
  expect_equal(out$value[out$gate == "p" & out$check == "n_parent"], 100)
})

test_that("per-sample gate lists are used for their own sample", {
  m <- grid_events()
  g1 <- sw_gate_range("g", "A", 0, 5, parent = NA)
  g2 <- sw_gate_range("g", "A", 0, 3, parent = NA)
  out <- sw_gate_qc(list(a = m, b = m), list(a = list(g1), b = list(g2)),
    reference = NULL, edge_width = 0.1)
  ks <- out[out$check == "kept_share", ]
  expect_equal(ks$value[ks$sample == "a"], 0.5)
  expect_equal(ks$value[ks$sample == "b"], 0.3)
})

test_that("sw_gate_qc stops on bad input", {
  m <- grid_events()
  g <- sw_gate_range("g", "A", 0, 5, parent = NA)
  expect_error(sw_gate_qc(m, list(g), reference = NULL, edge_width = 1), "named list")
  expect_error(sw_gate_qc(list(s = m), list(g), reference = "zz", edge_width = 1), "reference")
  expect_error(sw_gate_qc(list(s = m), list(g), reference = NULL, edge_width = -1), "positive")
  expect_error(sw_gate_qc(list(s = m), list(g), reference = NULL), "edge_width")
  orphan <- sw_gate_range("o", "A", 0, 5, parent = "nope")
  expect_error(sw_gate_qc(list(s = m), list(orphan), reference = NULL, edge_width = 1), "Parent gate")
})

test_that("kept_share_log_ratio flags a gate that keeps everything when kept_share_z cannot", {
  set.seed(8)
  # Reference: two populations, the gate keeps the upper one (about 30%).
  ref <- cbind(A = c(stats::rnorm(700, 0, 0.3), stats::rnorm(300, 3, 0.3)))
  # Three samples where a mis-placed gate keeps every event.
  others <- lapply(1:3, function(i) cbind(A = stats::rnorm(500, 5, 0.3)))
  data <- c(list(r = ref), stats::setNames(others, c("a", "b", "c")))
  g <- sw_gate_range("pos", "A", 1.5, Inf, parent = NA)
  out <- sw_gate_qc(data, list(g), reference = "r", edge_width = 0.2)
  ks <- out[out$check == "kept_share", ]
  expect_equal(ks$value[ks$sample != "r"], c(1, 1, 1))
  z <- out[out$check == "kept_share_z", ]
  expect_true(all(is.na(z$value)))
  lr <- out[out$check == "kept_share_log_ratio", ]
  expect_equal(lr$value[lr$sample == "r"], 0)
  ref_share <- mean(ref[, "A"] >= 1.5)
  expect_equal(lr$value[lr$sample == "a"], log(1 / ref_share))
  expect_gt(min(lr$value[lr$sample != "r"]), 1)
})
