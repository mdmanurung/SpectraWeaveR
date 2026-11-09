# tests/testthat/test-cytoscan_wrap.R
# Wrapper tests. CytoScan and CytoTools are replaced by fakes through
# .sw_cytoscan_fn(), so the tests need neither package.

fake_cytoscan <- function(seen) {
  list(
    CytoScan = function(truncate_max_range) {
      seen$truncate <- truncate_max_range
      structure(list(), class = "CytoScan")
    },
    addTestdata = function(CS, input) {
      CS$paths$test <- input
      CS
    },
    addReferencedata = function(CS, input) {
      CS$paths$reference <- input
      CS
    },
    generateFeatures = function(
      CS,
      channels,
      featMethod,
      n,
      quantileDist,
      aggSize,
      cores
    ) {
      seen$gen <- list(
        featMethod = featMethod,
        n = n,
        q = quantileDist,
        agg = aggSize
      )
      CS$features$test[[featMethod]] <- data.frame(x = seq_along(CS$paths$test))
      CS
    },
    Flag = function(
      CS,
      featMethod,
      flagMethod,
      outlier_threshold,
      novelty_threshold,
      ntrees
    ) {
      seen$flag <- list(
        m = flagMethod,
        o = outlier_threshold,
        n = novelty_threshold,
        t = ntrees
      )
      if (flagMethod == "outlier") {
        CS$outliers[[featMethod]] <- factor(c(FALSE, TRUE, FALSE)[seq_along(
          CS$paths$test
        )])
      } else {
        v <- c(FALSE, TRUE, FALSE)[seq_along(CS$paths$test)]
        names(v) <- CS$paths$test
        CS$novelties[[featMethod]] <- v
      }
      CS
    }
  )
}

run_wrap <- function(...) {
  args <- list(
    test_paths = c("a.fcs", "b.fcs", "c.fcs"),
    reference_paths = character(0),
    channels = c("A", "B"),
    feat_method = "quantiles",
    flag_method = "outlier",
    n = 100,
    quantile_dist = 0.1,
    agg_size = 1000,
    outlier_threshold = 0.6,
    novelty_threshold = 0.05,
    ntrees = 50,
    preprocess = sw_cytoscan_preprocess_none,
    truncate_max_range = FALSE,
    cores = 1
  )
  args[names(list(...))] <- list(...)
  do.call(sw_cytoscan_run, args)
}

test_that("required arguments have no default", {
  expect_error(sw_cytoscan_run(), "no default")
  expect_error(
    sw_cytoscan_run(test_paths = c("a", "b")),
    "reference_paths"
  )
})

test_that("it stops with install help when CytoScan is absent", {
  skip_if(requireNamespace("CytoScan", quietly = TRUE))
  expect_error(run_wrap(), "remotes::install_github")
})

test_that("outlier run returns one flag per file and swaps in the caller's preprocess", {
  seen <- new.env()
  fake <- fake_cytoscan(seen)
  local_mocked_bindings(.sw_cytoscan_fn = function(name, pkg = "CytoScan") {
    fake[[name]]
  })
  mine <- function(ff) ff
  res <- run_wrap(preprocess = mine)
  expect_s3_class(res, "sw_cytoscan")
  expect_equal(res$file, c("a.fcs", "b.fcs", "c.fcs"))
  expect_equal(res$flag, c(FALSE, TRUE, FALSE))
  expect_equal(unique(res$kind), "outlier")
  expect_false(seen$truncate)
  expect_equal(seen$gen$agg, 1000)
  expect_equal(seen$flag$o, 0.6)
  expect_equal(seen$flag$t, 50)
  expect_equal(nrow(attr(res, "features")), 3)
})

test_that("novelty needs reference files and returns logical flags", {
  seen <- new.env()
  fake <- fake_cytoscan(seen)
  local_mocked_bindings(.sw_cytoscan_fn = function(name, pkg = "CytoScan") {
    fake[[name]]
  })
  expect_error(run_wrap(flag_method = "novelty"), "needs `reference_paths`")
  res <- run_wrap(
    flag_method = "novelty",
    reference_paths = c("r1.fcs", "r2.fcs")
  )
  expect_equal(res$flag, c(FALSE, TRUE, FALSE))
  expect_equal(unique(res$kind), "novelty")
})

test_that("argument checks fire before CytoScan is touched", {
  expect_error(run_wrap(flag_method = "other"), "outlier")
  expect_error(run_wrap(outlier_threshold = 1.2), "outlier_threshold")
  expect_error(run_wrap(preprocess = "arcsinh"), "must be a function")
  expect_error(run_wrap(n = 0), "whole number")
})

test_that("preprocess_none returns the same object", {
  x <- list(a = 1)
  expect_identical(sw_cytoscan_preprocess_none(x), x)
})

test_that("sw_cytotools_header returns a tibble from CytoTools", {
  local_mocked_bindings(
    .sw_cytoscan_fn = function(name, pkg = "CytoScan") {
      expect_equal(pkg, "CytoTools")
      function(paths) data.frame(file = paths, voltage_group = "VOLTAGEGROUP1")
    }
  )
  res <- sw_cytotools_header(c("a.fcs", "b.fcs"))
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 2)
  expect_error(sw_cytotools_header(), "required")
})
