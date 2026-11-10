# Tests for R/relativise.R (reference-population relativisation).

make_control <- function(neg_med, pos_med, n = 300, sd = 0.05) {
  # Two log-normal groups with known medians on the linear scale.
  v <- c(neg_med * exp(stats::rnorm(n, 0, sd)), pos_med * exp(stats::rnorm(n, 0, sd)))
  cbind(B1 = v, B2 = rev(v))
}

test_that("fit recovers the negative and positive medians", {
  set.seed(1)
  ctrl <- list(c1 = make_control(100, 5000))
  map <- data.frame(unit = "s1", channel = "B1", control = "c1")
  pos <- list(c1 = rep(c(FALSE, TRUE), each = 300))
  model <- sw_relativise_fit(ctrl, map, pos, min_events = 50)
  expect_s3_class(model, "sw_relativise_model")
  expect_equal(model$m, stats::median(ctrl$c1[1:300, "B1"]))
  expect_equal(model$M, stats::median(ctrl$c1[301:600, "B1"]))
  expect_equal(model$gain, model$M - model$m)
  expect_equal(c(model$n_negative, model$n_positive), c(300L, 300L))
})

test_that("k-means split finds the bright group and restores the RNG state", {
  set.seed(2)
  ctrl <- list(c1 = make_control(100, 5000))
  map <- data.frame(unit = "s1", channel = "B1", control = "c1")
  before <- .Random.seed
  pos <- sw_relativise_split_kmeans(ctrl, map, function(v) asinh(v / 150), seed = 1, nstart = 3)
  expect_identical(.Random.seed, before)
  expect_equal(pos$c1, rep(c(FALSE, TRUE), each = 300))
})

test_that("k-means split refuses a control that calibrates two channels", {
  set.seed(3)
  ctrl <- list(c1 = make_control(100, 5000))
  map <- data.frame(unit = "s1", channel = c("B1", "B2"), control = "c1")
  expect_error(sw_relativise_split_kmeans(ctrl, map, identity, 1, 3), "more than one channel")
})

test_that("apply maps the negative median to low and the positive median to high", {
  set.seed(4)
  ctrl <- list(c1 = make_control(100, 5000), c2 = make_control(300, 9000))
  map <- data.frame(unit = c("s1", "s2"), channel = "B1", control = c("c1", "c2"))
  pos <- list(c1 = rep(c(FALSE, TRUE), each = 300), c2 = rep(c(FALSE, TRUE), each = 300))
  model <- sw_relativise_fit(ctrl, map, pos, min_events = 50)
  data <- list(
    f1 = cbind(B1 = c(model$m[1], model$M[1]), B2 = c(1, 2)),
    f2 = cbind(B1 = c(model$m[2], model$M[2]), B2 = c(1, 2))
  )
  out <- sw_relativise_apply(
    data,
    model,
    unit = c(f1 = "s1", f2 = "s2"),
    channels = "B1",
    low = 0,
    high = 1e4,
    input_scale = "linear"
  )
  expect_equal(out$f1[, "B1"], c(0, 1e4))
  expect_equal(out$f2[, "B1"], c(0, 1e4))
  expect_equal(out$f1[, "B2"], c(1, 2)) # not requested, unchanged
  expect_equal(nrow(attr(out, "relativisation")$rows), 2)
})

test_that("apply removes a per-unit offset and gain from the same biology", {
  set.seed(5)
  truth <- stats::rexp(500, 1 / 2000)
  # Two instruments: offset a, gain g. The control sees the same instrument.
  ctrl <- list(
    c1 = make_control(50 + 1 * 100, 50 + 1 * 5000),
    c2 = make_control(400 + 2.5 * 100, 400 + 2.5 * 5000)
  )
  map <- data.frame(unit = c("i1", "i2"), channel = "B1", control = c("c1", "c2"))
  pos <- list(c1 = rep(c(FALSE, TRUE), each = 300), c2 = rep(c(FALSE, TRUE), each = 300))
  model <- sw_relativise_fit(ctrl, map, pos, min_events = 50)
  data <- list(
    f1 = cbind(B1 = 50 + 1 * truth, B2 = 0),
    f2 = cbind(B1 = 400 + 2.5 * truth, B2 = 0)
  )
  out <- sw_relativise_apply(data, model, c(f1 = "i1", f2 = "i2"), "B1", 0, 1e4, "linear")
  # After relativisation the two files agree up to control noise.
  expect_lt(max(abs(out$f1[, "B1"] - out$f2[, "B1"])) / 1e4, 0.05)
  expect_gt(max(abs(data$f1[, "B1"] - data$f2[, "B1"])) / 1e4, 1)
})

test_that("fit and apply stop on the documented conditions", {
  set.seed(6)
  ctrl <- list(c1 = make_control(100, 5000))
  map <- data.frame(unit = "s1", channel = "B1", control = "c1")
  pos <- list(c1 = rep(c(FALSE, TRUE), each = 300))
  expect_error(sw_relativise_fit(ctrl, map, pos, min_events = 301), "min_events")
  expect_error(
    sw_relativise_fit(ctrl, map, list(c1 = rep(c(TRUE, FALSE), each = 300)), 50),
    "not above"
  )
  expect_error(sw_relativise_fit(ctrl, map, list(c1 = TRUE), 50), "one value per event")
  expect_error(sw_relativise_fit(ctrl, map, list(), 50), "no entry")
  bad_map <- data.frame(unit = "s1", channel = "B9", control = "c1")
  expect_error(sw_relativise_fit(ctrl, bad_map, pos, 50), "no column")
  dup_map <- rbind(map, map)
  expect_error(sw_relativise_fit(ctrl, dup_map, pos, 50), "more than one control")

  model <- sw_relativise_fit(ctrl, map, pos, min_events = 50)
  data <- list(f1 = cbind(B1 = 1:3, B2 = 1:3))
  expect_error(
    sw_relativise_apply(data, model, c(f1 = "s1"), "B1", 0, 1e4, "arcsinh"),
    "linear"
  )
  expect_error(sw_relativise_apply(data, model, c(f1 = "s2"), "B1", 0, 1e4, "linear"), "no reference")
  expect_error(sw_relativise_apply(data, model, c(f1 = "s1"), "B1", 1, 1, "linear"), "above")
  expect_error(sw_relativise_apply(data, model, c(f9 = "s1"), "B1", 0, 1, "linear"), "named")
})

test_that("diagnose flags a unit and a batch with a doubled gain", {
  set.seed(7)
  units <- paste0("s", 1:8)
  gains <- c(rep(1, 6), 2, 2)
  ctrl <- lapply(gains, \(g) make_control(100 * g, 5000 * g * exp(stats::rnorm(1, 0, 0.02))))
  names(ctrl) <- units
  map <- data.frame(unit = units, channel = "B1", control = units)
  pos <- lapply(ctrl, \(x) rep(c(FALSE, TRUE), each = 300))
  model <- sw_relativise_fit(ctrl, map, pos, min_events = 50)
  dx <- sw_relativise_diagnose(
    model,
    batch = stats::setNames(rep(c("b1", "b2"), c(6, 2)), units),
    order = stats::setNames(seq_along(units), units)
  )
  z <- dx[dx$check == "gain_z", ]
  expect_true(all(z$unit[order(-abs(z$value))][1:2] %in% c("s7", "s8")))
  ratio <- dx[dx$check == "gain_batch_ratio", ]
  expect_gt(ratio$value[ratio$unit == "b2"], 1.8)
  expect_equal(nrow(dx[dx$check == "gain_order_rho", ]), 1)
  no_extra <- sw_relativise_diagnose(model, batch = NULL, order = NULL)
  expect_setequal(unique(no_extra$check), c("offset", "gain", "separation", "gain_z"))
  expect_error(sw_relativise_diagnose(model, batch = c(s1 = "b1"), order = NULL), "every unit")
})

test_that("batch map builds unit = batch and enforces one marker per channel", {
  info <- data.frame(
    batch = rep(c("b1", "b2"), each = 2),
    channel = rep(c("B1", "B2"), 2),
    control = c("b1_B1", "b1_B2", "b2_B1", "b2_B2"),
    reference_marker = rep(c("CD4", "CD8"), 2)
  )
  map <- sw_relativise_batch_map(info)
  expect_equal(map$unit, c("b1", "b1", "b2", "b2"))
  expect_equal(names(map), c("unit", "channel", "control"))

  mixed <- info
  mixed$reference_marker[4] <- "CD3"
  expect_error(sw_relativise_batch_map(mixed), "different reference markers")
  expect_error(sw_relativise_batch_map(info[-4, ]), "b2 lacks B2")
  expect_error(sw_relativise_batch_map(rbind(info, info[1, ])), "More than one control")
  expect_error(sw_relativise_batch_map(info[1:3]), "columns")
})

test_that("internal map uses each file as its own control for one channel", {
  set.seed(8)
  data <- list(s1 = make_control(100, 5000), s2 = make_control(200, 9000))
  map <- sw_relativise_internal_map(data, "B1")
  expect_equal(map$unit, c("s1", "s2"))
  expect_equal(map$control, c("s1", "s2"))
  pos <- sw_relativise_split_kmeans(data, map, \(v) asinh(v / 150), seed = 1, nstart = 3)
  model <- sw_relativise_fit(data, map, pos, min_events = 50)
  out <- sw_relativise_apply(data, model, c(s1 = "s1", s2 = "s2"), "B1", 0, 1e4, "linear")
  expect_equal(stats::median(out$s1[301:600, "B1"]), 1e4)
  expect_equal(out$s1[, "B2"], data$s1[, "B2"])
  expect_error(sw_relativise_internal_map(data, "B9"), "without channel")
})

test_that("sw_diagnose runs the relativise module without data", {
  set.seed(9)
  units <- paste0("s", 1:5)
  ctrl <- lapply(c(1, 1, 1, 1, 2), \(g) make_control(100 * g, 5000 * g))
  names(ctrl) <- units
  map <- data.frame(unit = units, channel = "B1", control = units)
  pos <- lapply(ctrl, \(x) rep(c(FALSE, TRUE), each = 300))
  model <- sw_relativise_fit(ctrl, map, pos, min_events = 50)
  dx <- sw_diagnose(
    NULL,
    header = NULL,
    signal = NULL,
    anomaly = NULL,
    batch = NULL,
    cluster = NULL,
    relativise = list(model = model, batch = NULL, order = stats::setNames(1:5, units)),
    otd = NULL,
    negative = NULL
  )
  expect_equal(attr(dx, "modules"), "relativise")
  expect_setequal(unique(dx$check), c("offset", "gain", "separation", "gain_z", "gain_order_rho"))
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(top$subject[top$check == "gain_z"], "s5")
})

test_that("separation is gain over the summed robust SDs and is low for overlapping sides", {
  set.seed(5)
  wide <- c(stats::rnorm(400, 100, 10), stats::rnorm(400, 5000, 10))
  close <- c(stats::rnorm(400, 100, 300), stats::rnorm(400, 600, 300))
  ctrl <- list(a = cbind(B1 = wide), b = cbind(B1 = close))
  map <- data.frame(unit = c("a", "b"), channel = "B1", control = c("a", "b"))
  pos <- list(a = rep(c(FALSE, TRUE), each = 400), b = rep(c(FALSE, TRUE), each = 400))
  model <- sw_relativise_fit(ctrl, map, pos, min_events = 50)
  exp_sep <- (model$M[1] - model$m[1]) / (stats::mad(wide[1:400]) + stats::mad(wide[401:800]))
  expect_equal(model$separation[1], exp_sep)
  expect_gt(model$separation[1], 100)
  expect_lt(model$separation[2], 2)
  dx <- sw_relativise_diagnose(model, batch = NULL, order = NULL)
  expect_equal(dx$value[dx$check == "separation"], model$separation)
})
