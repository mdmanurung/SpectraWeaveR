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
  time_bin = 5,
  area_height = NULL,
  area_height_quantile = NULL,
  reference_group = NULL,
  reference = NULL
)

test_that("every argument must be given", {
  expect_error(sw_diagnose(), "no default")
  x <- make_dx()
  expect_error(
    sw_diagnose(x$data, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL),
    "at least one module"
  )
  expect_error(
    sw_diagnose(
      NULL,
      NULL,
      signal_set,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
    "`data` is needed"
  )
  expect_error(
    sw_diagnose(
      x$data,
      "paths",
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
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
      order = NULL,
      reference_group = NULL
    ),
    cluster = NULL,
    relativise = NULL,
    otd = NULL,
    negative = NULL,
    gate = NULL
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
    reference_group = NULL,
    ntrees = 100,
    sample_size = 10,
    ndim = 1
  )
  dx <- sw_diagnose(x$data, NULL, NULL, an, NULL, NULL, NULL, NULL, NULL, NULL)
  expect_equal(unique(dx$check), "outlier_score")
  expect_equal(nrow(dx), 15)
  an$novelty_reference <- paste0("f", c(1:2, 4:8, 10:15)) # 13 files > 9 features
  an$novelty_quantile <- 0.05
  dx2 <- sw_diagnose(x$data, NULL, NULL, an, NULL, NULL, NULL, NULL, NULL, NULL)
  expect_setequal(unique(dx2$check), c("outlier_score", "novelty_log_density"))
  # With EMD features, one emd_z row per file and channel; a split marker tops it.
  d3 <- x$data
  d3$f7[, "B"] <- d3$f7[, "B"] + sample(c(-2, 2), nrow(d3$f7), replace = TRUE)
  an_emd <- an[c(
    "channels",
    "n_events",
    "seed",
    "scale",
    "ntrees",
    "sample_size",
    "ndim"
  )]
  an_emd <- c(
    an_emd,
    list(
      features = "emd",
      quantile_step = NULL,
      aggregate_from = "all",
      agg_size = 1500,
      reference = NULL,
      reference_group = NULL
    )
  )
  dx3 <- sw_diagnose(d3, NULL, NULL, an_emd, NULL, NULL, NULL, NULL, NULL, NULL)
  ez <- dx3[dx3$check == "emd_z", ]
  expect_equal(nrow(ez), 15 * 3)
  expect_equal(ez$subject[which.max(ez$value)], "f7")
  expect_equal(ez$item[which.max(ez$value)], "B")
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
    batch = NULL,
    cluster = NULL,
    relativise = NULL,
    otd = NULL,
    negative = NULL,
    gate = NULL
  )
  hit <- dx[dx$subject == "b.fcs" & dx$check == "keyword $CYT", ]
  expect_equal(hit$value, 1)
  expect_match(hit$text, "differs: Y vs X")
  expect_equal(dx$value[dx$subject == "a.fcs" & dx$check == "keyword $CYT"], 0)
})

test_that("sw_diagnosis_top ranks per check and needs n", {
  x <- make_dx()
  dx <- sw_diagnose(
    x$data,
    NULL,
    signal_set,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
  )
  expect_error(sw_diagnosis_top(dx), "n")
  expect_error(sw_diagnosis_top(data.frame(), n = 2), "must come from")
  top <- sw_diagnosis_top(dx, n = 2)
  expect_true(all(top$rank <= 2))
  sat <- top[top$check == "saturation_share" & top$item == "A", ]
  expect_equal(nrow(sat), 2)
  expect_equal(sat$subject[sat$rank == 1], "f3")
})

test_that("signed checks rank by size, so a large negative z comes first", {
  x <- make_dx()
  # Lower marker B of one file by 3 SD: its median_z is large and negative.
  x$data$f7[, "B"] <- x$data$f7[, "B"] - 3
  dx <- sw_diagnose(
    x$data,
    NULL,
    signal_set,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
  )
  z <- dx[dx$check == "median_z" & dx$item == "B", ]
  expect_lt(z$value[z$subject == "f7"], -3)
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(top$subject[top$check == "median_z" & top$item == "B"], "f7")
})

test_that("the cluster module ranks each marker and cluster as its own item", {
  x <- make_dx()
  dx <- sw_diagnose(
    x$data,
    header = NULL,
    signal = NULL,
    anomaly = NULL,
    batch = NULL,
    cluster = list(
      markers = c("A", "B", "C"),
      k = 3,
      n_events = 200,
      seed = 1,
      nstart = 3,
      min_cluster_events = 10,
      reference = NULL,
      reference_group = NULL
    ),
    relativise = NULL,
    otd = NULL,
    negative = NULL,
    gate = NULL
  )
  expect_equal(attr(dx, "modules"), "cluster")
  expect_setequal(unique(dx$check), c("freq_z", "median_z"))
  expect_true(all(grepl("^n_events=", dx$text)))
  expect_setequal(unique(dx$item[dx$check == "freq_z"]), paste("cluster", 1:3))
  expect_setequal(
    unique(dx$item[dx$check == "median_z"]),
    paste(rep(c("A", "B", "C"), 3), "/ cluster", rep(1:3, each = 3))
  )
  # Each marker and cluster is ranked on its own: n = 1 gives one row per item.
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(sum(top$check == "median_z"), 9)
})

test_that("ranking keys: signed, low-first and ratio checks, without spurious warnings", {
  key <- .sw_diag_key(
    c(
      "neg_corr",
      "separation",
      "gain_batch_ratio",
      "offset",
      "novelty_log_density"
    ),
    c(-0.9, 0.5, 0.5, -3, -10)
  )
  expect_equal(key, c(0.9, -0.5, abs(log(0.5)), -3, 10))
  expect_no_warning(.sw_diag_key(c("offset", "gain_batch_ratio"), c(-5, 2)))
})

test_that("bimodality_z ranks by size; spread_z and rest_emd_z stay one-sided (G38)", {
  # A split blurs two populations, so the coefficient drops: a large negative z is the fault.
  key <- .sw_diag_key(
    c("bimodality_z", "spread_z", "rest_emd_z"),
    c(-4, -4, -4)
  )
  expect_equal(key, c(4, -4, -4))
})

test_that("the gate module stacks sw_gate_qc rows with gate/channel items", {
  set.seed(3)
  files <- lapply(c(0, 0, 0, 0, 1.5), function(shift) {
    cbind(A = stats::rnorm(500, 5 + shift), B = stats::rnorm(500, 5))
  })
  names(files) <- paste0("f", 1:5)
  g <- sw_gate_range("pos", "A", 4, Inf, parent = NA)
  dx <- sw_diagnose(
    files,
    header = NULL,
    signal = NULL,
    anomaly = NULL,
    batch = NULL,
    cluster = NULL,
    relativise = NULL,
    otd = NULL,
    negative = NULL,
    gate = list(gates = list(g), reference = "f1", edge_width = 0.5)
  )
  expect_equal(attr(dx, "modules"), "gate")
  expect_true(all(dx$module == "gate"))
  expect_true("pos/A" %in% dx$item[dx$check == "centroid_shift"])
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(top$subject[top$check == "kept_share_z"], "f5")
})
