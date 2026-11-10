# tests/testthat/test-fcs_signal_audit.R
# Unit tests for R/fcs_signal_audit.R on small simulated files.

make_sig <- function(n_files = 5, n = 400, seed = 1) {
  set.seed(seed)
  f <- lapply(seq_len(n_files), function(i) {
    cbind(
      Time = sort(stats::runif(n, 0, 100)),
      A = stats::rnorm(n),
      B = stats::rnorm(n),
      C = stats::rnorm(n)
    )
  })
  names(f) <- paste0("f", seq_len(n_files))
  f
}

audit <- function(d, ...) {
  args <- list(
    data = d,
    channels = c("A", "B", "C"),
    range_max = NULL,
    range_min = NULL,
    time_channel = NULL,
    time_bin = NULL,
    area_height = NULL,
    area_height_quantile = NULL,
    reference_group = NULL,
    reference = NULL
  )
  args[names(list(...))] <- list(...)
  do.call(sw_fcs_signal_audit, args)
}

val <- function(res, file, check, channel = NA) {
  r <- res[res$file == file & res$check == check, ]
  if (!is.na(channel)) {
    r <- r[r$channel == channel, ]
  }
  r$value
}

test_that("required arguments have no default", {
  expect_error(sw_fcs_signal_audit(), "no default")
  d <- make_sig()
  expect_error(
    sw_fcs_signal_audit(
      d,
      channels = c("A", "B"),
      area_height = NULL,
      area_height_quantile = NULL,
      reference_group = NULL,
      reference = NULL
    ),
    "range_max"
  )
})

test_that("clean files give the base checks only", {
  res <- audit(make_sig())
  expect_s3_class(res, "sw_fcs_signal_audit")
  expect_setequal(
    unique(res$check),
    c(
      "n_events",
      "nonfinite_share",
      "zero_variance",
      "negative_share",
      "median_z",
      "sd_z",
      "bimodality",
      "bimodality_z",
      "corr_shift",
      "corr_pair_z",
      "global_sd_z"
    )
  )
  expect_equal(val(res, "f1", "n_events"), 400)
  expect_equal(val(res, "f1", "zero_variance", "A"), 0)
})

test_that("median_z and sd_z flag a shifted and a stretched marker", {
  d <- make_sig(n_files = 12, n = 800)
  d$f2[, "A"] <- d$f2[, "A"] + 0.5
  d$f3[, "B"] <- mean(d$f3[, "B"]) + (d$f3[, "B"] - mean(d$f3[, "B"])) * 1.5
  # A split into two modes: +-0.5 on unit noise raises the SD by about 12%.
  set.seed(4)
  d$f4[, "C"] <- d$f4[, "C"] + sample(c(-0.5, 0.5), 800, replace = TRUE)
  res <- audit(d)
  expect_gt(val(res, "f2", "median_z", "A"), 3)
  expect_gt(val(res, "f3", "sd_z", "B"), 3)
  expect_gt(val(res, "f4", "sd_z", "C"), 2)
  # The same files are not flagged on the other channels.
  expect_lt(abs(val(res, "f2", "median_z", "B")), 3)
  expect_lt(abs(val(res, "f1", "sd_z", "B")), 3)
  # Fewer than 3 files: no z-score.
  two <- audit(make_sig(n_files = 2))
  expect_true(all(is.na(two$value[two$check %in% c("median_z", "sd_z")])))
})

test_that("a dead channel, non-finite values and saturation are found", {
  d <- make_sig()
  d$f2[, "B"] <- 0
  d$f3[1:40, "C"] <- NA
  d$f4[1:100, "A"] <- 1000
  res <- audit(d, range_max = 1000, range_min = -1000)
  expect_equal(val(res, "f2", "zero_variance", "B"), 1)
  expect_equal(val(res, "f1", "zero_variance", "B"), 0)
  expect_equal(val(res, "f3", "nonfinite_share", "C"), 0.1)
  expect_equal(val(res, "f4", "saturation_share", "A"), 0.25)
  expect_equal(val(res, "f1", "saturation_share", "A"), 0)
  expect_equal(val(res, "f1", "floor_share", "A"), 0)
})

test_that("rate_cv counts events per bin, so thinning a window raises it", {
  # Regression: rate_cv once counted occupied bins, not events, so a window
  # with 80% fewer events (no empty bin) looked like a clean file.
  d <- make_sig(n = 2000)
  win <- d$f5[, "Time"] > 40 & d$f5[, "Time"] < 70
  set.seed(9)
  drop <- win & stats::runif(nrow(d$f5)) < 0.8
  d$f5 <- d$f5[!drop, ]
  res <- audit(d, time_channel = "Time", time_bin = 5)
  clean <- c(val(res, "f1", "rate_cv"), val(res, "f2", "rate_cv"))
  expect_true(val(res, "f5", "rate_cv") > 2 * max(clean))
  # Exact value on a hand-made file: bins of 1 unit hold 3, 1, 0, 2 events.
  m <- cbind(Time = c(0.1, 0.2, 0.3, 1.5, 3.2, 3.9), A = 1:6, B = 6:1)
  d2 <- list(x = m, y = m)
  r2 <- sw_fcs_signal_audit(
    d2,
    c("A", "B"),
    NULL,
    NULL,
    "Time",
    1,
    area_height = NULL,
    area_height_quantile = NULL,
    reference_group = NULL,
    reference = NULL
  )
  counts <- c(3, 1, 0, 2)
  expect_equal(
    r2$value[r2$check == "rate_cv"][1],
    stats::sd(counts) / mean(counts)
  )
})

test_that("a marker drifting with time shows in time_trend for that channel only", {
  d <- make_sig(n = 600)
  d$f2[, "A"] <- d$f2[, "A"] + 0.03 * d$f2[, "Time"] # about 3 SD over the run
  res <- audit(d, time_channel = "Time", time_bin = 5)
  # A ramp 0 to 3 on unit noise correlates at about sqrt(0.75 / 1.75) = 0.65.
  expect_equal(val(res, "f2", "time_trend", "A"), 0.65, tolerance = 0.1)
  expect_lt(abs(val(res, "f2", "time_trend", "B")), 0.2)
  expect_lt(abs(val(res, "f1", "time_trend", "A")), 0.2)
  # Without a time channel the check is absent.
  expect_false("time_trend" %in% audit(d)$check)
})

test_that("a gap in event times raises rate_cv and max_gap", {
  d <- make_sig()
  keep <- !(d$f5[, "Time"] > 40 & d$f5[, "Time"] < 70)
  d$f5 <- d$f5[keep, ]
  res <- audit(d, time_channel = "Time", time_bin = 5)
  expect_true(
    val(res, "f5", "rate_cv") >
      2 * max(val(res, "f1", "rate_cv"), val(res, "f2", "rate_cv"))
  )
  expect_true(val(res, "f5", "max_gap") > 10 * val(res, "f1", "max_gap"))
})

test_that("a changed channel correlation is named by pair", {
  d <- make_sig(n = 600)
  d$f3[, "B"] <- d$f3[, "A"] * 0.9 + stats::rnorm(600, sd = 0.3)
  res <- audit(d)
  r <- res[res$check == "corr_shift", ]
  expect_equal(r$file[which.max(r$value)], "f3")
  expect_equal(r$detail[r$file == "f3"], "A ~ B")
  expect_true(r$value[r$file == "f3"] > 0.5)
})

test_that("limits can be per channel and must cover every channel", {
  d <- make_sig()
  expect_error(audit(d, range_max = c(A = 1, B = 1)), "name every channel")
  res <- audit(d, range_max = c(A = 1000, B = 1000, C = 1000))
  expect_true("saturation_share" %in% res$check)
  expect_error(audit(d, time_channel = "Time"), "time_bin")
  expect_error(audit(d, time_channel = "T2", time_bin = 1), "lacks column")
})

test_that("area_height_z finds the file with more doublets", {
  set.seed(41)
  files <- lapply(1:10, function(i) {
    n <- 2000
    h <- stats::rlnorm(n, log(1e5), 0.2)
    a <- h * stats::rnorm(n, 1.2, 0.03)
    if (i == 4) {
      # 15% doublets: two cells pass together, area adds up, height barely rises.
      k <- sample.int(n, 0.15 * n)
      a[k] <- a[k] + a[sample.int(n, length(k))]
    }
    cbind(`FSC-A` = a, `FSC-H` = h, X = stats::rnorm(n), Y = stats::rnorm(n))
  })
  names(files) <- paste0("f", 1:10)
  au <- sw_fcs_signal_audit(
    files,
    c("X", "Y"),
    NULL,
    NULL,
    time_channel = NULL,
    time_bin = NULL,
    area_height = c(`FSC-A` = "FSC-H"),
    area_height_quantile = 0.95,
    reference_group = NULL,
    reference = NULL
  )
  z <- au[au$check == "area_height_z", ]
  expect_equal(nrow(z), 10)
  expect_equal(z$file[which.max(z$value)], "f4")
  expect_gt(max(z$value), 5)
  expect_equal(unique(z$detail), "FSC-A / FSC-H")
  expect_error(
    sw_fcs_signal_audit(
      files,
      c("X", "Y"),
      NULL,
      NULL,
      NULL,
      NULL,
      c(`FSC-A` = "FSC-H"),
      NULL,
      reference_group = NULL,
      reference = NULL
    ),
    "area_height_quantile"
  )
})

test_that("global_sd_z finds a file widened a little in every channel", {
  set.seed(31)
  files <- lapply(1:12, function(i) {
    w <- if (i == 7) 1.15 else 1
    m <- matrix(stats::rnorm(800 * 6, sd = w), ncol = 6)
    colnames(m) <- paste0("c", 1:6)
    m
  })
  names(files) <- paste0("f", 1:12)
  au <- sw_fcs_signal_audit(
    files,
    paste0("c", 1:6),
    NULL,
    NULL,
    time_channel = NULL,
    time_bin = NULL,
    area_height = NULL,
    area_height_quantile = NULL,
    reference_group = NULL,
    reference = NULL
  )
  g <- au[au$check == "global_sd_z", ]
  expect_equal(nrow(g), 12)
  expect_true(all(is.na(g$channel)))
  expect_equal(g$file[which.max(g$value)], "f7")
  expect_gt(max(g$value), 3)
})

test_that("median_z is (x - median) / (1.4826 * MAD) exactly; corr_pair_z has one row per pair", {
  set.seed(21)
  meds <- c(0, 1, 2, 3, 10)
  files <- lapply(meds, \(m) {
    cbind(
      A = m + stats::rnorm(1001, 0, 1e-9),
      B = stats::rnorm(1001),
      C = stats::rnorm(1001)
    )
  })
  names(files) <- paste0("f", 1:5)
  au <- sw_fcs_signal_audit(
    files,
    c("A", "B", "C"),
    NULL,
    NULL,
    time_channel = NULL,
    time_bin = NULL,
    area_height = NULL,
    area_height_quantile = NULL,
    reference_group = NULL,
    reference = NULL
  )
  z <- au$value[au$check == "median_z" & au$channel == "A"]
  # Median of the file medians is 2 and their raw MAD is 1, so z = (m - 2) / 1.4826.
  expect_equal(z, (meds - 2) / 1.4826, tolerance = 1e-6)
  pz <- au[au$check == "corr_pair_z", ]
  expect_equal(nrow(pz), 5 * 3)
  expect_setequal(unique(pz$channel), c("A ~ B", "A ~ C", "B ~ C"))
})

# reference_group ---------------------------------------------------------------

# Three donors of four files. Donors differ in the median of A by 0, 4 and 8.
# File d2_4 is shifted by 0.5 on top of its donor.
make_donors <- function(seed = 3) {
  set.seed(seed)
  ids <- paste0("d", rep(1:3, each = 4), "_", rep(1:4, 3))
  d <- lapply(ids, function(i) {
    donor <- as.integer(substr(i, 2, 2))
    cbind(
      A = stats::rnorm(400) + 4 * (donor - 1) + if (i == "d2_4") 0.5 else 0,
      B = stats::rnorm(400),
      C = stats::rnorm(400)
    )
  })
  names(d) <- ids
  donor <- stats::setNames(substr(ids, 1, 2), ids)
  list(data = d, donor = donor)
}

test_that("reference_group finds a file that is unusual only within its donor", {
  x <- make_donors()
  all_files <- audit(x$data)
  by_donor <- audit(x$data, reference_group = x$donor)
  z_all <- val(all_files, "d2_4", "median_z", "A")
  z_grp <- val(by_donor, "d2_4", "median_z", "A")
  expect_lt(abs(z_all), 1)
  expect_gt(abs(z_grp), 3)
  # The faulty file ranks first on A only with groups.
  rank_of <- function(res) {
    r <- res[res$check == "median_z" & res$channel == "A", ]
    r$file[which.max(abs(r$value))]
  }
  expect_equal(rank_of(by_donor), "d2_4")
  # The result does not depend on the order of the labels.
  shuffled <- audit(x$data, reference_group = x$donor[rev(names(x$donor))])
  expect_equal(shuffled$value, by_donor$value)
  # sd_z, global_sd_z and corr_pair_z keep their rows.
  expect_setequal(unique(by_donor$check), unique(all_files$check))
})

test_that("a group with fewer than 3 files gives NA z-scores", {
  x <- make_donors()
  g <- x$donor
  g[names(g) %in% c("d3_1", "d3_2")] <- "other"
  res <- audit(x$data, reference_group = g)
  zrows <- res[
    res$check %in% c("median_z", "sd_z", "global_sd_z", "corr_pair_z"),
  ]
  small <- zrows$file %in% c("d3_1", "d3_2")
  expect_true(all(is.na(zrows$value[small])))
  expect_false(all(is.na(zrows$value[!small])))
})

test_that("reference_group is required and must cover every file", {
  x <- make_donors()
  expect_error(
    sw_fcs_signal_audit(
      x$data,
      c("A", "B"),
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
    "no default"
  )
  expect_error(
    sw_fcs_signal_audit(
      x$data,
      c("A", "B"),
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
    "reference_group"
  )
  expect_error(
    audit(x$data, reference_group = x$donor[-1]),
    "label for every file"
  )
  expect_error(
    audit(x$data, reference_group = unname(x$donor)),
    "label for every file"
  )
  expect_error(
    audit(
      x$data,
      reference_group = stats::setNames(rep(1, 12), names(x$donor))
    ),
    "label for every file"
  )
})

# reference_group, two-way (data frame) --------------------------------------------------

# Three donors x four batches, one file per cell. The median of A carries a donor
# effect (0, 4, 8) and a batch effect (0, 3, 6, 9). File d2_b3 is shifted by 1 more.
make_two_way <- function(seed = 8) {
  set.seed(seed)
  g <- expand.grid(
    donor = paste0("d", 1:3),
    batch = paste0("b", 1:4),
    stringsAsFactors = FALSE
  )
  g$file <- paste(g$donor, g$batch, sep = "_")
  d <- lapply(seq_len(nrow(g)), function(i) {
    cbind(
      A = stats::rnorm(400) +
        c(0, 4, 8)[as.integer(substr(g$donor[i], 2, 2))] +
        c(0, 3, 6, 9)[as.integer(substr(g$batch[i], 2, 2))] +
        if (g$file[i] == "d2_b3") 1 else 0,
      B = stats::rnorm(400),
      C = stats::rnorm(400)
    )
  })
  names(d) <- g$file
  list(
    data = d,
    design = g[c("file", "donor", "batch")],
    donor = stats::setNames(g$donor, g$file)
  )
}

test_that("a two-way reference_group finds a file unusual only after donor and batch effects", {
  x <- make_two_way()
  top <- function(res) {
    r <- res[res$check == "median_z" & res$channel == "A", ]
    c(file = r$file[which.max(abs(r$value))], z = max(abs(r$value)))
  }
  none <- top(audit(x$data))
  one_way <- top(audit(x$data, reference_group = x$donor))
  two_way <- top(audit(x$data, reference_group = x$design))
  expect_false(none[["file"]] == "d2_b3")
  expect_false(one_way[["file"]] == "d2_b3")
  expect_equal(two_way[["file"]], "d2_b3")
  expect_gt(as.numeric(two_way[["z"]]), 5)
  # The row order of the design does not matter.
  shuffled <- audit(
    x$data,
    reference_group = x$design[rev(seq_len(nrow(x$design))), ]
  )
  expect_equal(shuffled$value, audit(x$data, reference_group = x$design)$value)
})

test_that("a two-way reference_group needs a file column, two character columns and every file", {
  x <- make_two_way()
  expect_error(
    audit(x$data, reference_group = x$design[c("donor", "batch")]),
    "file. column"
  )
  expect_error(
    audit(x$data, reference_group = x$design["file"]),
    "exactly two more"
  )
  expect_error(
    audit(x$data, reference_group = cbind(x$design, extra = "e")),
    "exactly two more"
  )
  expect_error(
    audit(x$data, reference_group = x$design[-1, ]),
    "covers every file"
  )
  expect_error(
    audit(x$data, reference_group = rbind(x$design, x$design[1, ])),
    "once"
  )
  num <- x$design
  num$batch <- 1
  expect_error(audit(x$data, reference_group = num), "character")
})

test_that("the bimodality coefficient is 1/3 for a normal sample and above 5/9 for a two-mode mix", {
  set.seed(11)
  # Normal: skewness and excess kurtosis are 0, so BC = 1 / 3. At n = 1e5 the
  # SD of BC is about 0.002; 0.01 is five SD.
  expect_equal(.sw_bimodality(stats::rnorm(1e5)), 1 / 3, tolerance = 0.01)
  # Equal mix of N(-2, 1) and N(+2, 1): population BC is about 0.60. At
  # n = 2000 the SD of BC is about 0.01, so 5 / 9 is cleared by several SD.
  mix <- c(stats::rnorm(1000, -2), stats::rnorm(1000, 2))
  expect_gt(.sw_bimodality(mix), 5 / 9)
  # Too few values, or no variance: NA.
  expect_true(is.na(.sw_bimodality(c(1, 2, 3))))
  expect_true(is.na(.sw_bimodality(rep(1, 50))))
})

test_that("bimodality_z is largest for the one file with a bimodal channel", {
  d <- make_sig(n_files = 8, n = 800)
  set.seed(5)
  d$f3[, "B"] <- c(stats::rnorm(400, -2), stats::rnorm(400, 2))
  res <- audit(d)
  b <- res[res$check == "bimodality_z" & res$channel == "B", ]
  expect_equal(nrow(b), 8)
  expect_identical(b$file[which.max(b$value)], "f3")
  expect_gt(max(b$value), 3)
  # One row per file and channel for the raw coefficient too.
  expect_equal(sum(res$check == "bimodality"), 8 * 3)
  expect_gt(val(res, "f3", "bimodality", "B"), 5 / 9)
  # Fewer than 3 files: no z-score.
  two <- audit(make_sig(n_files = 2))
  expect_true(all(is.na(two$value[two$check == "bimodality_z"])))
})

# reference (frozen) -------------------------------------------------------------

test_that("reference = NULL keeps the cohort z and is required", {
  d <- make_sig(6, seed = 4)
  res <- audit(d)
  med <- vapply(d, \(m) stats::median(m[, "A"]), numeric(1))
  hand <- (med - stats::median(med)) / (1.4826 * stats::mad(med, constant = 1))
  r <- res[res$check == "median_z" & res$channel == "A", ]
  expect_equal(stats::setNames(r$value, r$file)[names(d)], hand)
  expect_error(
    sw_fcs_signal_audit(
      d,
      c("A", "B"),
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      reference_group = NULL
    ),
    "'reference'"
  )
})
