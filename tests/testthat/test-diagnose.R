# tests/testthat/test-diagnose.R
# Tests for R/diagnose.R: the stacked table and the per-check ranking.

make_dx <- function(seed = 1, n = 300) {
  set.seed(seed)
  ids <- paste0("f", 1:15)
  batch <- stats::setNames(rep(c("b1", "b2", "b3"), each = 5), ids)
  d <- lapply(ids, function(i) {
    cbind(
      Time = sort(stats::runif(n, 0, 100)),
      A = stats::rnorm(n) + if (batch[[i]] == "b2") 1 else 0,
      B = stats::rnorm(n),
      C = stats::rnorm(n)
    )
  })
  names(d) <- ids
  d$f3[1:90, "A"] <- 6 # saturated file
  list(data = d, batch = batch)
}

signal_set <- list(
  channels = c("A", "B", "C"),
  range_max = 6,
  range_min = NULL,
  time_channel = "Time",
  time_bin = 5
)

test_that("every argument must be given", {
  expect_error(sw_diagnose(), "no default")
  x <- make_dx()
  expect_error(
    sw_diagnose(x$data, NULL, NULL, NULL, NULL),
    "at least one module"
  )
  expect_error(
    sw_diagnose(NULL, NULL, signal_set, NULL, NULL),
    "`data` is needed"
  )
  expect_error(
    sw_diagnose(x$data, "paths", NULL, NULL, NULL),
    "must be NULL or a list"
  )
})

test_that("signal and batch modules stack into one table", {
  x <- make_dx()
  dx <- sw_diagnose(
    x$data,
    header = NULL,
    signal = signal_set,
    anomaly = NULL,
    batch = list(
      batch = x$batch,
      markers = c("A", "B", "C"),
      n_events = 200,
      seed = 1,
      order = NULL
    )
  )
  expect_s3_class(dx, "sw_diagnosis")
  expect_equal(attr(dx, "modules"), c("signal", "batch"))
  expect_setequal(
    names(dx),
    c("module", "subject", "check", "item", "value", "text")
  )
  expect_setequal(unique(dx$module), c("signal", "batch"))
  sat <- dx[dx$check == "saturation_share" & dx$item == "A", ]
  expect_equal(sat$subject[which.max(sat$value)], "f3")
  eta <- dx[dx$check == "eta2", ]
  expect_equal(eta$item[which.max(eta$value)], "A")
})

test_that("the anomaly module keeps forest scores, and novelty when asked", {
  skip_if_not_installed("isotree")
  skip_if_not_installed("mclust")
  x <- make_dx()
  an <- list(
    channels = c("A", "B", "C"),
    features = "quantile",
    quantile_step = 0.25,
    n_events = 200,
    seed = 1,
    scale = "test",
    ntrees = 100,
    sample_size = 10,
    ndim = 1
  )
  dx <- sw_diagnose(x$data, NULL, NULL, an, NULL)
  expect_equal(unique(dx$check), "outlier_score")
  expect_equal(nrow(dx), 15)
  an$novelty_reference <- paste0("f", c(1:2, 4:8, 10:15)) # 13 files > 9 features
  an$novelty_quantile <- 0.05
  dx2 <- sw_diagnose(x$data, NULL, NULL, an, NULL)
  expect_setequal(unique(dx2$check), c("outlier_score", "novelty_log_density"))
})

test_that("the header module reads FCS headers", {
  skip_if_not_installed("flowCore")
  d <- withr::local_tempdir()
  mk <- function(name, cyt) {
    ff <- flowCore::flowFrame(matrix(
      stats::rnorm(40),
      ncol = 2,
      dimnames = list(NULL, c("A", "B"))
    ))
    flowCore::keyword(ff) <- list("$CYT" = cyt)
    p <- file.path(d, name)
    flowCore::write.FCS(ff, p)
    p
  }
  paths <- c(mk("r.fcs", "X"), mk("a.fcs", "X"), mk("b.fcs", "Y"))
  dx <- sw_diagnose(
    NULL,
    header = list(paths = paths, keywords = "$CYT", reference = paths[1]),
    signal = NULL,
    anomaly = NULL,
    batch = NULL
  )
  hit <- dx[dx$subject == "b.fcs" & dx$check == "keyword $CYT", ]
  expect_equal(hit$value, 1)
  expect_match(hit$text, "differs: Y vs X")
  expect_equal(dx$value[dx$subject == "a.fcs" & dx$check == "keyword $CYT"], 0)
})

test_that("sw_diagnosis_top ranks per check and needs n", {
  x <- make_dx()
  dx <- sw_diagnose(x$data, NULL, signal_set, NULL, NULL)
  expect_error(sw_diagnosis_top(dx), "n")
  expect_error(sw_diagnosis_top(data.frame(), n = 2), "must come from")
  top <- sw_diagnosis_top(dx, n = 2)
  expect_true(all(top$rank <= 2))
  sat <- top[top$check == "saturation_share" & top$item == "A", ]
  expect_equal(nrow(sat), 2)
  expect_equal(sat$subject[sat$rank == 1], "f3")
})
