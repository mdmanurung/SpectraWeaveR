# tests/testthat/test-utils.R
# Unit tests for R/utils.R — Format conversion utilities

# ---- Helper: create a mock flowFrame for testing ----
# Since flowCore may not be available, tests are guarded with skip_if_not_installed

test_that("sw_io_read_fcs rejects empty file vector", {
  expect_error(sw_io_read_fcs(character(0)), "non-empty character vector")
})

test_that("sw_io_read_fcs rejects non-character input", {
  expect_error(sw_io_read_fcs(42), "non-empty character vector")
})

test_that("sw_io_read_fcs rejects missing files", {
  expect_error(sw_io_read_fcs("/nonexistent/file.fcs"), "not found")
})

test_that("sw_io_ff_to_tibble converts flowFrame to tibble", {
  skip_if_not_installed("flowCore")

  mat <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A"))
  )
  ff <- flowCore::flowFrame(mat)
  result <- sw_io_ff_to_tibble(
    ff,
    sample = "S1",
    batch = "B1",
    condition = "Stim"
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 100)
  expect_true(all(
    c("FSC-A", "SSC-A", "BV421-A", "sample", "batch", "condition") %in%
      names(result)
  ))
  expect_equal(result$sample[1], "S1")
  expect_equal(result$batch[1], "B1")
  expect_equal(result$condition[1], "Stim")
})

test_that("sw_io_ff_to_tibble works without metadata", {
  skip_if_not_installed("flowCore")

  mat <- matrix(
    rnorm(200),
    ncol = 2,
    dimnames = list(NULL, c("FSC-A", "BV421-A"))
  )
  ff <- flowCore::flowFrame(mat)
  result <- sw_io_ff_to_tibble(ff)

  expect_s3_class(result, "tbl_df")
  expect_false("sample" %in% names(result))
  expect_false("batch" %in% names(result))
})

test_that("sw_io_ff_to_tibble rejects non-flowFrame", {
  skip_if_not_installed("flowCore")

  expect_error(sw_io_ff_to_tibble(data.frame(x = 1)), "flowFrame")
})

test_that("sw_io_tibble_to_ff converts tibble to flowFrame", {
  skip_if_not_installed("flowCore")

  df <- tibble::tibble(CD3 = rnorm(50), CD4 = rnorm(50), sample = "S1")
  ff <- sw_io_tibble_to_ff(df, markers = c("CD3", "CD4"))

  expect_true(methods::is(ff, "flowFrame"))
  expect_equal(nrow(ff), 50, ignore_attr = TRUE)
  expect_equal(ncol(ff), 2, ignore_attr = TRUE)
  expect_true(all(c("CD3", "CD4") %in% flowCore::colnames(ff)))
})

test_that("sw_io_tibble_to_ff uses all numeric cols when markers=NULL", {
  skip_if_not_installed("flowCore")

  df <- tibble::tibble(CD3 = rnorm(30), CD4 = rnorm(30))
  ff <- sw_io_tibble_to_ff(df)

  expect_true(methods::is(ff, "flowFrame"))
  expect_equal(ncol(ff), 2, ignore_attr = TRUE)
})

test_that("sw_io_tibble_to_ff rejects missing markers", {
  skip_if_not_installed("flowCore")

  df <- tibble::tibble(CD3 = rnorm(10))
  expect_error(sw_io_tibble_to_ff(df, markers = c("CD3", "CD99")), "not found")
})

test_that("sw_io_tibble_to_ff rejects non-data.frame", {
  skip_if_not_installed("flowCore")

  expect_error(sw_io_tibble_to_ff("not_a_df"), "data.frame")
})

test_that("sw_io_exprs_to_tibble converts matrix to tibble", {
  mat <- matrix(1:12, ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
  result <- sw_io_exprs_to_tibble(mat, sample = "S1", batch = 1)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 4)
  expect_true(all(c("A", "B", "C", "sample", "batch") %in% names(result)))
  expect_equal(result$sample[1], "S1")
})

test_that("sw_io_exprs_to_tibble renames columns", {
  mat <- matrix(1:6, ncol = 2)
  result <- sw_io_exprs_to_tibble(mat, colnames = c("X", "Y"))

  expect_true(all(c("X", "Y") %in% names(result)))
})

test_that("sw_io_exprs_to_tibble rejects mismatched colnames", {
  mat <- matrix(1:6, ncol = 2)
  expect_error(sw_io_exprs_to_tibble(mat, colnames = c("X", "Y", "Z")), "ncol")
})

test_that("sw_io_exprs_to_tibble rejects non-matrix input", {
  expect_error(sw_io_exprs_to_tibble("not_a_matrix"), "matrix or data.frame")
})

test_that("sw_channel_get_fluor excludes scatter and time", {
  skip_if_not_installed("flowCore")

  mat <- matrix(
    rnorm(500),
    ncol = 5,
    dimnames = list(NULL, c("FSC-A", "SSC-A", "Time", "BV421-A", "PE-A"))
  )
  ff <- flowCore::flowFrame(mat)
  fluor <- sw_channel_get_fluor(ff)

  expect_equal(sort(fluor), sort(c("BV421-A", "PE-A")))
  expect_false("FSC-A" %in% fluor)
  expect_false("SSC-A" %in% fluor)
  expect_false("Time" %in% fluor)
})

test_that("sw_channel_get_fluor rejects invalid input", {
  skip_if_not_installed("flowCore")

  expect_error(sw_channel_get_fluor("not_a_ff"), "flowFrame or flowSet")
})

test_that("sw_channel_set_markers renames channels", {
  skip_if_not_installed("flowCore")

  mat <- matrix(
    rnorm(200),
    ncol = 2,
    dimnames = list(NULL, c("BV421-A", "PE-A"))
  )
  ff <- flowCore::flowFrame(mat)
  marker_map <- c("BV421-A" = "CD3", "PE-A" = "CD4")
  ff2 <- sw_channel_set_markers(ff, marker_map)

  pdata <- flowCore::pData(flowCore::parameters(ff2))
  expect_equal(pdata$desc[pdata$name == "BV421-A"], "CD3", ignore_attr = TRUE)
  expect_equal(pdata$desc[pdata$name == "PE-A"], "CD4", ignore_attr = TRUE)
})

test_that("sw_channel_set_markers warns on missing channel", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(100), ncol = 1, dimnames = list(NULL, c("BV421-A")))
  ff <- flowCore::flowFrame(mat)
  marker_map <- c("NONEXISTENT" = "CD3")

  expect_warning(sw_channel_set_markers(ff, marker_map), "not found")
})

test_that("sw_channel_set_markers rejects unnamed vector", {
  skip_if_not_installed("flowCore")

  mat <- matrix(rnorm(100), ncol = 1, dimnames = list(NULL, c("BV421-A")))
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_channel_set_markers(ff, c("CD3")), "named character vector")
})

# Round-trip test: flowFrame → tibble → flowFrame
test_that("round-trip flowFrame -> tibble -> flowFrame preserves data", {
  skip_if_not_installed("flowCore")

  mat <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  ff_orig <- flowCore::flowFrame(mat)

  tbl <- sw_io_ff_to_tibble(ff_orig)
  ff_back <- sw_io_tibble_to_ff(tbl, markers = c("CD3", "CD4", "CD8"))

  expect_equal(
    flowCore::exprs(ff_back),
    flowCore::exprs(ff_orig),
    tolerance = 1e-10
  )
})

# =========================================================================
# sw_channel_is_signal / sw_channel_is_fluor
# =========================================================================

test_that("sw_channel_is_signal excludes metadata channels", {
  skip_if_not_installed("flowCore")
  mat <- matrix(
    rnorm(500),
    ncol = 5,
    dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A", "Time", "Original_ID"))
  )
  ff <- flowCore::flowFrame(mat)
  result <- sw_channel_is_signal(ff)
  expect_true(result["FSC-A"])
  expect_true(result["SSC-A"])
  expect_true(result["BV421-A"])
  expect_false(result["Time"])
  expect_false(result["Original_ID"])
})

test_that("sw_channel_is_fluor excludes scatter and metadata", {
  skip_if_not_installed("flowCore")
  mat <- matrix(
    rnorm(500),
    ncol = 5,
    dimnames = list(NULL, c("FSC-A", "SSC-A", "BV421-A", "Time", "APC-A"))
  )
  ff <- flowCore::flowFrame(mat)
  result <- sw_channel_is_fluor(ff)
  expect_false(result["FSC-A"])
  expect_false(result["SSC-A"])
  expect_true(result["BV421-A"])
  expect_false(result["Time"])
  expect_true(result["APC-A"])
})

test_that("AF Index (AutoSpectral AF label) is neither signal nor fluor", {
  skip_if_not_installed("flowCore")
  mat <- matrix(
    rnorm(400),
    ncol = 4,
    dimnames = list(NULL, c("FSC-A", "AF-A", "AF Index", "BV421-A"))
  )
  ff <- flowCore::flowFrame(mat)
  expect_false(sw_channel_is_signal(ff)[["AF Index"]])
  expect_false(sw_channel_is_fluor(ff)[["AF Index"]])
  expect_true(sw_channel_is_fluor(ff)[["AF-A"]])
})

test_that("sw_channel_is_signal rejects non-flowFrame/flowSet", {
  skip_if_not_installed("flowCore")
  expect_error(sw_channel_is_signal("not_ff"), "flowFrame or flowSet")
})

test_that("sw_channel_is_signal handles custom patterns", {
  skip_if_not_installed("flowCore")
  mat <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("FSC-A", "MyCustom", "BV421-A"))
  )
  ff <- flowCore::flowFrame(mat)
  result <- sw_channel_is_signal(ff, exclude_patterns = c("MyCustom"))
  expect_true(result["FSC-A"])
  expect_false(result["MyCustom"])
  expect_true(result["BV421-A"])
})

# =========================================================================
# sw_io_subsample
# =========================================================================

test_that("sw_io_subsample rejects non-flowSet", {
  skip_if_not_installed("flowCore")
  expect_error(
    sw_io_subsample("not_fs", 100, setup = "forceNEvent", seed = 1),
    "flowSet or flowFrame"
  )
})

test_that("sw_io_subsample rejects invalid n_total_events", {
  skip_if_not_installed("flowCore")
  mat <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  ff <- flowCore::flowFrame(mat)
  expect_error(
    sw_io_subsample(ff, -1, setup = "forceNEvent", seed = 1),
    "positive integer"
  )
})

test_that("sw_io_subsample works with flowFrame input", {
  skip_if_not_installed("flowCore")
  mat <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  ff <- flowCore::flowFrame(mat)
  result <- sw_io_subsample(
    ff,
    n_total_events = 50,
    setup = "forceNEvent",
    seed = 42
  )
  expect_true(methods::is(result, "flowFrame"))
  expect_true(nrow(result) <= 50)
  expect_true("File" %in% flowCore::colnames(result))
  expect_true("Original_ID" %in% flowCore::colnames(result))
})

test_that("sw_io_subsample works with flowSet input", {
  skip_if_not_installed("flowCore")
  mat1 <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  mat2 <- matrix(
    rnorm(450),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  ff1 <- flowCore::flowFrame(mat1)
  ff2 <- flowCore::flowFrame(mat2)
  fs <- flowCore::flowSet(ff1, ff2)
  result <- sw_io_subsample(
    fs,
    n_total_events = 100,
    setup = "forceNEvent",
    seed = 42
  )
  expect_true(methods::is(result, "flowFrame"))
  expect_true(nrow(result) <= 100) # requested total
})

test_that("sw_io_subsample forceBalance strategy", {
  skip_if_not_installed("flowCore")
  mat1 <- matrix(
    rnorm(60),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  mat2 <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  ff1 <- flowCore::flowFrame(mat1)
  ff2 <- flowCore::flowFrame(mat2)
  fs <- flowCore::flowSet(ff1, ff2)
  result <- sw_io_subsample(
    fs,
    n_total_events = 100,
    setup = "forceBalance",
    seed = 42
  )
  expect_true(methods::is(result, "flowFrame"))
  # forceBalance: max 20*2=40 (min of 20 and 100 per frame)
  expect_true(nrow(result) <= 40)
})

test_that("sw_io_subsample errors when setup or seed is missing", {
  skip_if_not_installed("flowCore")
  mat <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  ff <- flowCore::flowFrame(mat)
  expect_error(sw_io_subsample(ff, 50, seed = 1), "'setup'")
  expect_error(sw_io_subsample(ff, 50, setup = "forceNEvent"), "'seed'")
  expect_error(sw_io_subsample(ff, 50), "'setup', 'seed'")
  expect_error(sw_io_subsample(ff, 50, setup = "bad", seed = 1), "forceNEvent")
})

test_that("sw_io_subsample is seed reproducible and restores the RNG", {
  skip_if_not_installed("flowCore")
  cols <- list(NULL, c("CD3", "CD4", "CD8"))
  fs <- flowCore::flowSet(
    flowCore::flowFrame(matrix(rnorm(600), ncol = 3, dimnames = cols)),
    flowCore::flowFrame(matrix(rnorm(900), ncol = 3, dimnames = cols))
  )
  set.seed(7)
  before <- get(".Random.seed", envir = .GlobalEnv)
  a <- sw_io_subsample(fs, 100, setup = "forceNEvent", seed = 42)
  after <- get(".Random.seed", envir = .GlobalEnv)
  b <- sw_io_subsample(fs, 100, setup = "forceNEvent", seed = 42)
  d <- sw_io_subsample(fs, 100, setup = "forceNEvent", seed = 43)
  ids <- function(x) flowCore::exprs(x)[, "Original_ID"]
  expect_identical(before, after)
  expect_identical(ids(a), ids(b))
  expect_false(identical(ids(a), ids(d)))
})

test_that("sw_io_subsample seed = NULL uses and advances the caller's RNG stream", {
  skip_if_not_installed("flowCore")
  mat <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  ff <- flowCore::flowFrame(mat)
  ids <- function(x) flowCore::exprs(x)[, "Original_ID"]
  set.seed(1)
  before <- get(".Random.seed", envir = .GlobalEnv)
  a <- sw_io_subsample(ff, 50, setup = "forceNEvent", seed = NULL)
  after <- get(".Random.seed", envir = .GlobalEnv)
  # same caller stream gives the same draw; a fixed seed = 1 draw is the same
  # only if the code seeded itself, so it must differ from the stream draw
  set.seed(1)
  b <- sw_io_subsample(ff, 50, setup = "forceNEvent", seed = NULL)
  fixed <- sw_io_subsample(ff, 50, setup = "forceNEvent", seed = 99)
  expect_false(identical(before, after))
  expect_identical(ids(a), ids(b))
  expect_false(identical(ids(a), ids(fixed)))
  # a second call continues the stream, so it differs from the first
  c2 <- sw_io_subsample(ff, 50, setup = "forceNEvent", seed = NULL)
  expect_false(identical(ids(b), ids(c2)))
})

test_that("sw_io_subsample removes .Random.seed again when none existed", {
  skip_if_not_installed("flowCore")
  mat <- matrix(
    rnorm(300),
    ncol = 3,
    dimnames = list(NULL, c("CD3", "CD4", "CD8"))
  )
  ff <- flowCore::flowFrame(mat)
  set.seed(5)
  saved <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
  rm(".Random.seed", envir = .GlobalEnv)
  sw_io_subsample(ff, 50, setup = "forceNEvent", seed = 3)
  expect_false(exists(".Random.seed", envir = .GlobalEnv))
})

# =========================================================================
# sw_io_event_audit
# =========================================================================

test_that("sw_io_event_audit rejects non-list", {
  expect_error(sw_io_event_audit("not_a_list"), "non-empty named list")
})

test_that("sw_io_event_audit rejects empty list", {
  expect_error(sw_io_event_audit(list()), "non-empty named list")
})

test_that("sw_io_event_audit rejects unnamed list", {
  expect_error(
    sw_io_event_audit(list(100, 90)),
    "named list"
  )
})

test_that("sw_io_event_audit works with numeric inputs", {
  result <- sw_io_event_audit(list(
    load = 1000,
    qc = 900,
    gate = 800,
    cluster = 750
  ))

  expect_equal(nrow(result), 4, ignore_attr = TRUE)
  expect_equal(
    result$step,
    c("load", "qc", "gate", "cluster"),
    ignore_attr = TRUE
  )
  expect_equal(result$n_events, c(1000L, 900L, 800L, 750L), ignore_attr = TRUE)
  expect_equal(result$pct_of_initial[1], 100, ignore_attr = TRUE)
  expect_equal(result$pct_of_initial[4], 75, ignore_attr = TRUE)
  expect_equal(result$pct_of_previous[1], 100, ignore_attr = TRUE)
})

test_that("sw_io_event_audit works with data.frames", {
  result <- sw_io_event_audit(list(
    step1 = data.frame(x = 1:100),
    step2 = data.frame(x = 1:80)
  ))
  expect_equal(result$n_events, c(100L, 80L), ignore_attr = TRUE)
  expect_equal(result$pct_of_initial[2], 80, ignore_attr = TRUE)
})

test_that("sw_io_event_audit works with flowFrame input", {
  skip_if_not_installed("flowCore")
  mat1 <- matrix(rnorm(300), ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
  mat2 <- matrix(rnorm(240), ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
  ff1 <- flowCore::flowFrame(mat1)
  ff2 <- flowCore::flowFrame(mat2)

  result <- sw_io_event_audit(list(
    load = ff1,
    filter = ff2
  ))
  expect_equal(result$n_events, c(100L, 80L), ignore_attr = TRUE)
})

test_that("sw_io_event_audit handles list of flowFrames", {
  skip_if_not_installed("flowCore")
  mat1 <- matrix(rnorm(150), ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
  mat2 <- matrix(rnorm(120), ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
  ff1 <- flowCore::flowFrame(mat1)
  ff2 <- flowCore::flowFrame(mat2)

  result <- sw_io_event_audit(list(
    step1 = list(sample1 = ff1, sample2 = ff2)
  ))
  expect_equal(result$n_events, 90L, ignore_attr = TRUE) # 50 + 40
})
