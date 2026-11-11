# tests/testthat/test-gating_utils.R
# Unit tests for R/gating_utils.R — Gating and event filtering utilities

# =========================================================================
# sw_gate_singlets
# =========================================================================

test_that("sw_gate_singlets rejects non-flowFrame", {
  skip_if_not_installed("flowCore")
  expect_error(sw_gate_singlets("not_ff", nmad = 4), "flowFrame")
})

test_that("sw_gate_singlets rejects invalid channel1", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_gate_singlets(ff, channel1 = 42, nmad = 4), "single channel name")
})

test_that("sw_gate_singlets rejects invalid channel2", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_gate_singlets(ff, channel2 = c("a", "b"), nmad = 4),
               "single channel name")
})

test_that("sw_gate_singlets rejects invalid nmad", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_gate_singlets(ff, nmad = -1), "positive number")
})

test_that("sw_gate_singlets rejects missing channel", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_gate_singlets(ff, channel1 = "MISSING", nmad = 4),
               "not found")
})

test_that("sw_gate_singlets returns polygonGate", {
  skip_if_not_installed("flowCore")
  set.seed(42)
  ch1 <- abs(rnorm(500, mean = 100000, sd = 20000))
  ch2 <- ch1 * 0.7 + rnorm(500, sd = 5000)
  mat <- matrix(c(ch1, ch2), ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  gate <- sw_gate_singlets(ff, nmad = 4)
  expect_true(methods::is(gate, "polygonGate"))
  expect_equal(gate@filterId, "Singlets")
})

test_that("sw_gate_singlets accepts custom filter_id", {
  skip_if_not_installed("flowCore")
  set.seed(42)
  ch1 <- abs(rnorm(500, mean = 100000, sd = 20000))
  ch2 <- ch1 * 0.7 + rnorm(500, sd = 5000)
  mat <- matrix(c(ch1, ch2), ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  gate <- sw_gate_singlets(ff, filter_id = "MySinglets", nmad = 4)
  expect_equal(gate@filterId, "MySinglets")
})

# =========================================================================
# sw_filter_doublets
# =========================================================================

test_that("sw_filter_doublets rejects non-flowFrame", {
  skip_if_not_installed("flowCore")
  expect_error(sw_filter_doublets("not_ff", nmads = 4), "flowFrame")
})

test_that("sw_filter_doublets rejects mismatched channel lengths", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(400)) + 1, ncol = 4,
                dimnames = list(NULL, c("FSC-A", "FSC-H", "SSC-A", "SSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_doublets(ff,
                       area_channels = c("FSC-A", "SSC-A"),
                       height_channels = c("FSC-H"),
                       nmads = c(4, 4)),
    "same length"
  )
})

test_that("sw_filter_doublets rejects mismatched nmads length", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(400)) + 1, ncol = 4,
                dimnames = list(NULL, c("FSC-A", "FSC-H", "SSC-A", "SSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_doublets(ff,
                       area_channels = c("FSC-A"),
                       height_channels = c("FSC-H"),
                       nmads = c(3, 5)),
    "same length"
  )
})

test_that("sw_filter_doublets rejects >2 channel pairs", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(600)) + 1, ncol = 6,
                dimnames = list(NULL, c("FSC-A", "FSC-H", "SSC-A",
                                        "SSC-H", "X-A", "X-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_doublets(ff,
                       area_channels = c("FSC-A", "SSC-A", "X-A"),
                       height_channels = c("FSC-H", "SSC-H", "X-H"),
                       nmads = c(4, 4, 4)),
    "length 1 or 2"
  )
})

test_that("sw_filter_doublets returns flowFrame with fewer events", {
  skip_if_not_installed("flowCore")
  set.seed(42)
  n <- 1000
  # Normal cells: linear FSC-A/FSC-H ratio
  fsc_a <- abs(rnorm(n, mean = 100000, sd = 20000))
  fsc_h <- fsc_a * 0.65 + rnorm(n, sd = 3000)
  # Add some doublets
  n_doublets <- 50
  fsc_a_d <- abs(rnorm(n_doublets, mean = 200000, sd = 20000))
  fsc_h_d <- fsc_a_d * 0.3 + rnorm(n_doublets, sd = 3000)

  mat <- matrix(
    c(c(fsc_a, fsc_a_d), c(fsc_h, fsc_h_d)),
    ncol = 2,
    dimnames = list(NULL, c("FSC-A", "FSC-H"))
  )
  ff <- flowCore::flowFrame(mat)
  ff_clean <- sw_filter_doublets(ff,
                                  area_channels = "FSC-A",
                                  height_channels = "FSC-H",
                                  nmads = 4)
  expect_true(methods::is(ff_clean, "flowFrame"))
  expect_true(nrow(ff_clean) <= nrow(ff))
})

# =========================================================================
# sw_filter_doublets_peacoqc
# =========================================================================

test_that("sw_filter_doublets_peacoqc rejects non-flowFrame", {
  skip_if_not_installed("PeacoQC")
  expect_error(sw_filter_doublets_peacoqc("not_ff", nmads = 4), "flowFrame")
})

test_that("sw_filter_doublets_peacoqc rejects mismatched channels", {
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(400)) + 1, ncol = 4,
                dimnames = list(NULL, c("FSC-A", "FSC-H", "SSC-A", "SSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_doublets_peacoqc(ff,
                                area_channels = c("FSC-A", "SSC-A"),
                                height_channels = "FSC-H",
                                nmads = c(4, 4)),
    "same length"
  )
})

test_that("sw_filter_doublets_peacoqc handles empty flowFrame", {
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")
  # flowCore's constructor cannot build a 0-event frame (it round-trips an FCS
  # file and $TOT = 0 reads as corrupt), so subset a 1-event frame instead.
  mat <- matrix(1, nrow = 1, ncol = 4,
                dimnames = list(NULL, c("FSC-A", "FSC-H", "SSC-A", "SSC-H")))
  ff <- flowCore::flowFrame(mat)[integer(0), ]
  result <- sw_filter_doublets_peacoqc(ff, nmads = c(4, 4))
  expect_equal(unname(nrow(result)), 0) # flowCore's nrow() is named "events"
  expect_equal(attr(result, "events_removed"), 0L)
})

# =========================================================================
# sw_filter_debris
# =========================================================================

test_that("sw_filter_debris rejects non-flowFrame", {
  skip_if_not_installed("flowCore")
  expect_error(sw_filter_debris("not_ff", gate_data = 1:6), "flowFrame")
})

test_that("sw_filter_debris rejects invalid fsc_channel", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "SSC-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_debris(ff, fsc_channel = 42, gate_data = 1:6),
    "single channel name"
  )
})

test_that("sw_filter_debris rejects odd-length gate_data", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "SSC-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_debris(ff, gate_data = 1:5),
    "even length"
  )
})

test_that("sw_filter_debris rejects too-short gate_data", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "SSC-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_debris(ff, gate_data = 1:4),
    "even length >= 6"
  )
})

test_that("sw_filter_debris rejects missing channel", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "SSC-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_debris(ff, fsc_channel = "MISSING",
                           gate_data = c(0, 1, 2, 0, 1, 2)),
    "not found"
  )
})

test_that("sw_filter_debris returns flowFrame", {
  skip_if_not_installed("flowCore")
  set.seed(42)
  n <- 500
  fsc <- abs(rnorm(n, mean = 100000, sd = 30000))
  ssc <- abs(rnorm(n, mean = 80000, sd = 25000))
  mat <- matrix(c(fsc, ssc), ncol = 2,
                dimnames = list(NULL, c("FSC-A", "SSC-A")))
  ff <- flowCore::flowFrame(mat)

  # Define a generous gate
  gate_coords <- c(50000, 50000, 200000, 200000,  # FSC coords
                    30000, 200000, 200000, 30000)   # SSC coords
  ff_clean <- sw_filter_debris(ff, gate_data = gate_coords)
  expect_true(methods::is(ff_clean, "flowFrame"))
  expect_true(nrow(ff_clean) <= nrow(ff))
})

# =========================================================================
# required arguments and events_removed
# =========================================================================

test_that("nmad and nmads have no default", {
  skip_if_not_installed("flowCore")
  expect_true(identical(formals(sw_gate_singlets)$nmad, quote(expr = )))
  expect_true(identical(formals(sw_filter_doublets)$nmads, quote(expr = )))
  expect_true(identical(formals(sw_filter_doublets_peacoqc)$nmads, quote(expr = )))
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_gate_singlets(ff), "'nmad' is required")
  expect_error(sw_filter_doublets(ff, area_channels = "FSC-A",
                                  height_channels = "FSC-H"),
               "'nmads' is required")
  expect_error(sw_filter_doublets_peacoqc(ff, area_channels = "FSC-A",
                                          height_channels = "FSC-H"),
               "'nmads' is required")
})

test_that("sw_filter_debris reports events_removed on a synthetic frame", {
  skip_if_not_installed("flowCore")
  # 10 events: 7 inside the square gate, 3 outside (FSC below 100)
  fsc <- c(rep(500, 7), rep(10, 3))
  ssc <- rep(500, 10)
  mat <- matrix(c(fsc, ssc), ncol = 2,
                dimnames = list(NULL, c("FSC-A", "SSC-A")))
  ff <- flowCore::flowFrame(mat)
  gate_coords <- c(100, 100, 1000, 1000, 100, 1000, 1000, 100)
  expect_silent(out <- sw_filter_debris(ff, gate_data = gate_coords))
  expect_equal(nrow(out), 7, ignore_attr = TRUE)
  expect_equal(attr(out, "events_in"), 10L)
  expect_equal(attr(out, "events_removed"), 3L)
  expect_message(
    sw_filter_debris(ff, gate_data = gate_coords, verbose = TRUE),
    "removed 3 / 10 events \\(30%\\)"
  )
})

test_that("sw_filter_doublets reports events_removed on a synthetic frame", {
  skip_if_not_installed("flowCore")
  set.seed(1)
  n <- 200
  a <- runif(n, 90000, 110000)
  h <- a * 0.7 + rnorm(n, sd = 100)
  # 5 clear doublets: height is far below 0.7 * area
  a_d <- rep(200000, 5)
  h_d <- rep(40000, 5)
  mat <- matrix(c(c(a, a_d), c(h, h_d)), ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  out <- sw_filter_doublets(ff, area_channels = "FSC-A",
                            height_channels = "FSC-H", nmads = 4)
  expect_equal(attr(out, "events_in"), nrow(ff), ignore_attr = TRUE)
  expect_equal(attr(out, "events_removed"), nrow(ff) - nrow(out),
               ignore_attr = TRUE)
  expect_gte(attr(out, "events_removed"), 5L)
  expect_message(
    sw_filter_doublets(ff, area_channels = "FSC-A",
                       height_channels = "FSC-H", nmads = 4,
                       verbose = TRUE),
    "sw_filter_doublets: removed"
  )
})

# =========================================================================
# sw_filter_margins_peacoqc
# =========================================================================

test_that("sw_filter_margins_peacoqc rejects non-flowFrame/flowSet", {
  skip_if_not_installed("PeacoQC")
  expect_error(sw_filter_margins_peacoqc("not_ff"), "flowFrame or flowSet")
})

test_that("sw_filter_margins_peacoqc rejects non-list channel_specifications", {
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(300)) + 1, ncol = 3,
                dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_margins_peacoqc(ff, channel_specifications = "not_a_list"),
    "list of lists"
  )
})

test_that("sw_filter_margins_peacoqc rejects wrong-length specifications", {
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(300)) + 1, ncol = 3,
                dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_margins_peacoqc(ff,
      channel_specifications = list("FSC-A" = c(0, 100, 200))),
    "exactly 2 values"
  )
})
