# Tests for R/negative_audit.R and R/cluster_stability.R.

skip_if_not_installed("mclust")

make_neg_file <- function(
  n = 1500,
  neg_sd = 0.5,
  neg_shift = 0,
  hyper = 0,
  angle = 0
) {
  # Channels A, B, C: negative N(0, neg_sd) with a donor-like share, positive N(3, 0.8).
  share <- stats::runif(1, 0.5, 0.7)
  x <- sapply(c("A", "B", "C"), function(ch) {
    neg <- stats::runif(n) < share
    ifelse(neg, stats::rnorm(n, 0, 0.5), stats::rnorm(n, 3, 0.8))
  })
  neg_b <- x[, "B"] < 1.5
  x[neg_b, "B"] <- x[neg_b, "B"] * neg_sd / 0.5 + neg_shift
  neg_a <- x[, "A"] < 1.5
  both <- neg_a & neg_b
  x[both, "B"] <- x[both, "B"] + angle * x[both, "A"]
  if (hyper > 0) {
    idx <- which(neg_b)[seq_len(round(hyper * sum(neg_b)))]
    x[idx, "B"] <- -4
  }
  x
}

make_cohort <- function(...) {
  files <- lapply(1:10, \(i) make_neg_file())
  files[[10]] <- make_neg_file(...)
  stats::setNames(files, paste0("f", 1:10))
}

test_that("a widened negative in one file has the largest neg_sd_z", {
  set.seed(1)
  files <- make_cohort(neg_sd = 1)
  na <- sw_negative_audit(
    files,
    c("A", "B", "C"),
    1000,
    seed = 1,
    hyper_k = 4,
    pairs = NULL,
    spread_quantile = NULL,
    min_separation = 0
  )
  expect_s3_class(na, "sw_negative_audit")
  z <- na[na$check == "neg_sd_z" & na$channel == "B", ]
  expect_equal(z$file[which.max(z$value)], "f10")
  expect_gt(max(z$value), 3)
  sep <- na[na$check == "separation", ]
  expect_true(all(sep$value > 2))
})

test_that("a shifted negative and hyper-negatives are found", {
  set.seed(2)
  files <- make_cohort(neg_shift = 0.4)
  na <- sw_negative_audit(
    files,
    c("A", "B"),
    1000,
    seed = 1,
    hyper_k = 4,
    pairs = NULL,
    spread_quantile = NULL,
    min_separation = 0
  )
  z <- na[na$check == "neg_median_z" & na$channel == "B", ]
  expect_equal(z$file[which.max(abs(z$value))], "f10")

  set.seed(3)
  files <- make_cohort(hyper = 0.05)
  na <- sw_negative_audit(
    files,
    c("A", "B"),
    1000,
    seed = 1,
    hyper_k = 4,
    pairs = NULL,
    spread_quantile = NULL,
    min_separation = 0
  )
  h <- na[na$check == "hyper_negative_share" & na$channel == "B", ]
  expect_equal(h$file[which.max(h$value)], "f10")
  expect_gt(max(h$value), 0.01)
})

test_that("an angled negative shows in neg_corr_z for that pair", {
  set.seed(4)
  files <- make_cohort(angle = 0.4)
  na <- sw_negative_audit(
    files,
    c("A", "B", "C"),
    1000,
    seed = 1,
    hyper_k = 4,
    pairs = "all",
    spread_quantile = NULL,
    min_separation = 0
  )
  z <- na[na$check == "neg_corr_z" & na$channel == "A ~ B", ]
  expect_equal(z$file[which.max(z$value)], "f10")
  expect_setequal(
    unique(na$channel[na$check == "neg_corr"]),
    c("A ~ B", "A ~ C", "B ~ C")
  )
})

test_that("negative audit restores the RNG state and checks its input", {
  set.seed(5)
  files <- make_cohort()
  before <- .Random.seed
  sw_negative_audit(
    files,
    "A",
    500,
    seed = 1,
    hyper_k = 4,
    pairs = NULL,
    spread_quantile = NULL,
    min_separation = 0
  )
  expect_identical(.Random.seed, before)
  expect_error(
    sw_negative_audit(
      files,
      "A",
      500,
      1,
      0,
      NULL,
      spread_quantile = NULL,
      min_separation = 0
    ),
    "hyper_k"
  )
  expect_error(
    sw_negative_audit(
      files,
      "Z",
      500,
      1,
      4,
      NULL,
      spread_quantile = NULL,
      min_separation = 0
    ),
    "lacks"
  )
  expect_error(
    sw_negative_audit(
      files,
      "A",
      5000,
      1,
      4,
      NULL,
      spread_quantile = NULL,
      min_separation = 0
    ),
    "fewer"
  )
  expect_error(
    sw_negative_audit(
      files,
      "A",
      500,
      1,
      4,
      "all",
      spread_quantile = NULL,
      min_separation = 0
    ),
    "two channels"
  )
  expect_error(
    sw_negative_audit(
      files,
      c("A", "B"),
      500,
      1,
      4,
      matrix(c("A", "Q"), 1),
      spread_quantile = NULL,
      min_separation = 0
    ),
    "not in"
  )
})

test_that("a constant channel gives NA with a reason, not an error", {
  set.seed(6)
  files <- make_cohort()
  files$f1[, "C"] <- 0
  na <- sw_negative_audit(
    files,
    c("A", "C"),
    500,
    seed = 1,
    hyper_k = 4,
    pairs = NULL,
    spread_quantile = NULL,
    min_separation = 0
  )
  row <- na[na$file == "f1" & na$channel == "C" & na$check == "neg_median", ]
  expect_true(is.na(row$value))
  expect_equal(row$detail, "constant channel")
  files$f2[1, "A"] <- NA
  na <- sw_negative_audit(
    files,
    c("A", "C"),
    1500,
    seed = 1,
    hyper_k = 4,
    pairs = NULL,
    spread_quantile = NULL,
    min_separation = 0
  )
  row <- na[na$file == "f2" & na$channel == "A" & na$check == "neg_median", ]
  expect_true(is.na(row$value))
  expect_match(row$detail, "non-finite")
})

test_that("sw_diagnose runs the negative module", {
  set.seed(7)
  files <- make_cohort()
  dx <- sw_diagnose(
    files,
    header = NULL,
    signal = NULL,
    anomaly = NULL,
    batch = NULL,
    cluster = NULL,
    relativise = NULL,
    otd = NULL,
    negative = list(
      channels = c("A", "B"),
      n_events = 500,
      seed = 1,
      hyper_k = 4,
      pairs = "all",
      spread_quantile = 0.9,
      min_separation = 3
    ),
    gate = NULL
  )
  expect_equal(attr(dx, "modules"), "negative")
  expect_true(all(
    c("neg_sd_z", "neg_corr_z", "hyper_negative_share", "spread_z") %in%
      dx$check
  ))
})

test_that("cluster stability: identical runs score 1, and well-separated clusters are stable", {
  set.seed(8)
  files <- lapply(1:4, function(i) {
    centre <- sample(c(-4, 0, 4), 600, replace = TRUE)
    cbind(A = centre + stats::rnorm(600, 0, 0.3), B = stats::rnorm(600, 0, 0.3))
  })
  names(files) <- paste0("f", 1:4)
  before <- .Random.seed
  st <- sw_cluster_stability(
    files,
    c("A", "B"),
    k = 3,
    seeds = 1:3,
    sizes = c(300, 1000),
    n_events = 500,
    draw_seed = 1,
    nstart = 3,
    iter_max = 50
  )
  expect_identical(.Random.seed, before)
  expect_equal(attr(st, "pool_size"), 2000)
  expect_true(all(st$value[st$check == "ari_seed_pair"] > 0.99))
  expect_true(all(st$value[st$check == "moved_share_seed_pair"] < 0.01))
  expect_equal(sort(unique(st$item[st$check == "ari_size"])), c("1000", "300"))
  expect_error(
    sw_cluster_stability(files, "A", 3, 1, 300, 500, 1, 3, 50),
    "two different"
  )
  expect_error(
    sw_cluster_stability(files, "A", 3, c(1, 1), 300, 500, 1, 3, 50),
    "two different"
  )
  expect_error(
    sw_cluster_stability(files, "A", 3, 1:2, 2, 500, 1, 3, 50),
    "at least `k`"
  )
  expect_error(
    sw_cluster_stability(files, "A", 3, 1:2, 5000, 500, 1, 3, 50),
    "below the pool"
  )
})

test_that("ARI and moved share agree with hand-worked cases", {
  a <- c(1, 1, 2, 2, 3, 3)
  expect_equal(.sw_ari(a, c(2, 2, 3, 3, 1, 1)), 1)
  expect_equal(.sw_moved_share(a, c(2, 2, 3, 3, 1, 1)), 0)
  expect_equal(.sw_moved_share(a, c(1, 1, 2, 2, 3, 1)), 1 / 6)
  expect_lt(.sw_ari(a, c(1, 2, 3, 1, 2, 3)), 0)
})

test_that("cluster stability reports k-means non-convergence instead of hiding it", {
  set.seed(10)
  files <- list(
    a = matrix(stats::rnorm(4000), ncol = 2),
    b = matrix(stats::rnorm(4000), ncol = 2)
  )
  files <- lapply(files, \(m) {
    colnames(m) <- c("A", "B")
    m
  })
  expect_warning(
    st <- sw_cluster_stability(
      files,
      c("A", "B"),
      8,
      1:2,
      500,
      2000,
      1,
      1,
      iter_max = 1
    ),
    "k-means warning"
  )
  expect_gt(length(attr(st, "kmeans_warnings")), 0)
})

test_that("spill_corr finds a leak from a bright channel into one that should be negative", {
  set.seed(11)
  files <- lapply(1:10, \(i) make_neg_file())
  # Breakdown-like leak in f10: A-positive events push signal into B where B is negative.
  x <- files[[10]]
  leak <- x[, "A"] > 1.5 & x[, "B"] < 1.5
  x[leak, "B"] <- x[leak, "B"] + 0.15 * x[leak, "A"]
  files[[10]] <- x
  names(files) <- paste0("f", 1:10)
  na <- sw_negative_audit(
    files,
    c("A", "B"),
    1500,
    seed = 1,
    hyper_k = 4,
    pairs = "all",
    spread_quantile = NULL,
    min_separation = 0
  )
  z <- na[na$check == "spill_corr_z" & na$channel == "A > B", ]
  expect_equal(z$file[which.max(z$value)], "f10")
  expect_setequal(
    unique(na$channel[na$check == "spill_corr"]),
    c("A > B", "B > A")
  )
})

test_that("spread_z finds a file whose negative widens with the source signal", {
  set.seed(12)
  files <- lapply(1:10, \(i) make_neg_file())
  # Spread in f10: B-negative events get zero-mean noise growing with A, so the
  # median and the correlation stay put while the width rises with A.
  x <- files[[10]]
  neg_b <- x[, "B"] < 1.5
  x[neg_b, "B"] <- x[neg_b, "B"] +
    stats::rnorm(sum(neg_b)) * 0.3 * pmax(x[neg_b, "A"], 0)
  files[[10]] <- x
  names(files) <- paste0("f", 1:10)
  na <- sw_negative_audit(
    files,
    c("A", "B"),
    1500,
    seed = 1,
    hyper_k = 4,
    pairs = "all",
    spread_quantile = 0.9,
    min_separation = 0
  )
  z <- na[na$check == "spread_z" & na$channel == "A > B", ]
  expect_equal(nrow(z), 10)
  expect_equal(z$file[which.max(z$value)], "f10")
  expect_gt(max(z$value), 3)
  r <- na[na$check == "spread_log_ratio" & na$channel == "A > B", ]
  expect_gt(r$value[r$file == "f10"], log(1.5))
  expect_error(
    sw_negative_audit(
      files,
      c("A", "B"),
      1500,
      1,
      4,
      NULL,
      spread_quantile = 0.9,
      min_separation = 0
    ),
    "needs `pairs`"
  )
  expect_error(
    sw_negative_audit(
      files,
      c("A", "B"),
      1500,
      1,
      4,
      "all",
      spread_quantile = 0.4,
      min_separation = 0
    ),
    "spread_quantile"
  )
  expect_error(
    sw_negative_audit(files, c("A", "B"), 1500, 1, 4, "all"),
    "no default"
  )
})

# Cohort: A and B bimodal, U unimodal (one population, no negative).
make_gate_cohort <- function(widen_b = FALSE) {
  files <- lapply(1:10, function(i) {
    x <- make_neg_file()
    cbind(x[, c("A", "B")], U = stats::rnorm(nrow(x), 0, 1))
  })
  if (widen_b) {
    x <- files[[10]]
    neg_b <- x[, "B"] < 1.5
    x[neg_b, "B"] <- x[neg_b, "B"] * 3
    files[[10]] <- x
  }
  stats::setNames(files, paste0("f", 1:10))
}

gate_args <- list(
  n_events = 1000,
  seed = 1,
  hyper_k = 4,
  pairs = "all",
  spread_quantile = 0.9
)

test_that("a unimodal channel is gated to NA with a reason, and split_ok is 0", {
  set.seed(20)
  files <- make_gate_cohort()
  na <- do.call(
    sw_negative_audit,
    c(list(files, c("A", "B", "U"), min_separation = 3), gate_args)
  )
  gated_checks <- c(
    "neg_median",
    "neg_sd",
    "neg_share",
    "hyper_negative_share",
    "neg_median_z",
    "neg_sd_z"
  )
  u <- na[na$channel == "U" & na$check %in% gated_checks, ]
  expect_equal(nrow(u), 6 * 10)
  expect_true(all(is.na(u$value)))
  expect_true(all(grepl(
    "^gated: median separation .* < min_separation 3",
    u$detail
  )))
  # Separation keeps its value; split_ok is 0 for U and 1 for A and B.
  sep <- na[na$check == "separation" & na$channel == "U", ]
  expect_false(anyNA(sep$value))
  ok <- na[na$check == "split_ok", ]
  expect_equal(nrow(ok), 3 * 10)
  expect_true(all(ok$value[ok$channel == "U"] == 0))
  expect_true(all(ok$value[ok$channel %in% c("A", "B")] == 1))
  # Pair checks with U are gated, including the directed and z rows.
  pair_checks <- c(
    "neg_corr",
    "neg_corr_z",
    "spill_corr",
    "spill_corr_z",
    "spread_log_ratio",
    "spread_z"
  )
  withu <- na[na$check %in% pair_checks & grepl("U", na$channel), ]
  expect_gt(nrow(withu), 0)
  expect_true(all(is.na(withu$value)))
  expect_true(all(grepl("^gated:", withu$detail)))
  # The A ~ B pair is scored.
  ab <- na[na$check == "neg_corr" & na$channel == "A ~ B", ]
  expect_false(anyNA(ab$value))
})

test_that("a bimodal channel is unchanged by the gate", {
  set.seed(21)
  files <- make_gate_cohort()
  gated <- do.call(
    sw_negative_audit,
    c(list(files, c("A", "B", "U"), min_separation = 3), gate_args)
  )
  open <- do.call(
    sw_negative_audit,
    c(list(files, c("A", "B", "U"), min_separation = 0), gate_args)
  )
  keep <- \(d) d[!grepl("U", d$channel) & d$check != "split_ok", ]
  expect_equal(as.data.frame(keep(gated)), as.data.frame(keep(open)))
  expect_true(all(open$value[open$check == "split_ok"] == 1))
})

test_that("one widened file in a bimodal channel stays scored", {
  set.seed(22)
  files <- make_gate_cohort(widen_b = TRUE)
  na <- do.call(
    sw_negative_audit,
    c(list(files, c("A", "B"), min_separation = 3), gate_args)
  )
  z <- na[na$check == "neg_sd_z" & na$channel == "B", ]
  expect_false(anyNA(z$value))
  expect_equal(z$file[which.max(z$value)], "f10")
})

test_that("min_separation is required and checked", {
  set.seed(23)
  files <- make_gate_cohort()
  expect_error(
    sw_negative_audit(files, "A", 500, 1, 4, NULL, spread_quantile = NULL),
    "no default"
  )
  expect_error(
    sw_negative_audit(
      files,
      "A",
      500,
      1,
      4,
      NULL,
      spread_quantile = NULL,
      min_separation = -1
    ),
    "min_separation"
  )
})
