# tests/testthat/test-gating_utils.R
# Unit tests for R/gating_utils.R — Gating and event filtering utilities

# =========================================================================
# sw_gate_singlets
# =========================================================================

test_that("sw_gate_singlets rejects non-flowFrame", {
  skip_if_not_installed("flowCore")
  expect_error(sw_gate_singlets("not_ff"), "flowFrame")
})

test_that("sw_gate_singlets rejects invalid channel1", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_gate_singlets(ff, channel1 = 42), "single channel name")
})

test_that("sw_gate_singlets rejects invalid channel2", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(200)) + 1, ncol = 2,
                dimnames = list(NULL, c("FSC-A", "FSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_gate_singlets(ff, channel2 = c("a", "b")),
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
  expect_error(sw_gate_singlets(ff, channel1 = "MISSING"),
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
  gate <- sw_gate_singlets(ff)
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
  gate <- sw_gate_singlets(ff, filter_id = "MySinglets")
  expect_equal(gate@filterId, "MySinglets")
})

# =========================================================================
# sw_filter_doublets
# =========================================================================

test_that("sw_filter_doublets rejects non-flowFrame", {
  skip_if_not_installed("flowCore")
  expect_error(sw_filter_doublets("not_ff"), "flowFrame")
})

test_that("sw_filter_doublets rejects mismatched channel lengths", {
  skip_if_not_installed("flowCore")
  mat <- matrix(abs(rnorm(400)) + 1, ncol = 4,
                dimnames = list(NULL, c("FSC-A", "FSC-H", "SSC-A", "SSC-H")))
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_filter_doublets(ff,
                       area_channels = c("FSC-A", "SSC-A"),
                       height_channels = c("FSC-H")),
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
                       height_channels = c("FSC-H", "SSC-H", "X-H")),
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
  expect_error(sw_filter_doublets_peacoqc("not_ff"), "flowFrame")
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
                                height_channels = "FSC-H"),
    "same length"
  )
})

test_that("sw_filter_doublets_peacoqc handles empty flowFrame", {
  skip_if_not_installed("PeacoQC")
  skip_if_not_installed("flowCore")
  mat <- matrix(numeric(0), ncol = 4,
                dimnames = list(NULL, c("FSC-A", "FSC-H", "SSC-A", "SSC-H")))
  ff <- flowCore::flowFrame(mat)
  result <- sw_filter_doublets_peacoqc(ff)
  expect_equal(nrow(result), 0)
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
