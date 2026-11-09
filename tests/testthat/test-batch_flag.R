# tests/testthat/test-batch_flag.R
# Unit tests for R/batch_flag.R. Pure data frames; no flow packages needed.

# Events for one control with fixed negative and positive values per batch.
make_fixed <- function(
  batches,
  neg,
  pos,
  n_neg = 5,
  n_pos = 5,
  control = "c1"
) {
  dplyr::bind_rows(lapply(seq_along(batches), function(i) {
    tibble::tibble(
      batch = batches[i],
      control = control,
      M1 = c(
        rep(neg[i], n_neg[(i - 1) %% length(n_neg) + 1]),
        rep(pos[i], n_pos[(i - 1) %% length(n_pos) + 1])
      )
    )
  }))
}

flag_fixed <- function(cells, quantile_type = 7, min_cells = 1, ...) {
  sw_flag_iqr(
    cells,
    markers = "M1",
    thresholds = c(M1 = 1),
    iqr_multiplier = 1.5,
    min_cells_per_side = min_cells,
    quantile_type = quantile_type,
    ...
  )
}

# Seeded simulation: two markers, one control, n batches.
# Pure rnorm() noise trips a Tukey fence on 8 batches by chance (seed 42 gave
# 5 noise flags), so each population is a shuffled normal-quantile grid. Its
# median is then exactly `centre + offset`, and the per-batch offsets are an
# evenly spaced, shuffled grid that never crosses the fence.
simulate_cells <- function(n_batches, n_per_batch = 2000) {
  grid_offsets <- function() {
    sample(seq(-0.05, 0.05, length.out = n_batches))
  }
  population <- function(k, centre, offset) {
    centre + offset + sample(stats::qnorm(stats::ppoints(k))) * 0.3
  }
  off <- list(grid_offsets(), grid_offsets(), grid_offsets(), grid_offsets())
  n_a <- 0.4 * n_per_batch
  n_b <- 0.6 * n_per_batch
  dplyr::bind_rows(lapply(seq_len(n_batches), function(b) {
    tibble::tibble(
      batch = sprintf("b%02d", b),
      control = "ref",
      # Fixed positive fractions, so percent positive does not vary by batch.
      CD3 = sample(c(
        population(n_a, 3, off[[1]][b]),
        population(n_per_batch - n_a, 0, off[[2]][b])
      )),
      CD19 = sample(c(
        population(n_b, 4, off[[3]][b]),
        population(n_per_batch - n_b, 0, off[[4]][b])
      ))
    )
  }))
}

flag_sim <- function(cells, min_cells = 100) {
  sw_flag_iqr(
    cells,
    markers = c("CD3", "CD19"),
    thresholds = c(CD3 = 1.5, CD19 = 2),
    iqr_multiplier = 1.5,
    min_cells_per_side = min_cells,
    quantile_type = 7
  )
}

test_that("each required argument without a default errors by name", {
  cells <- make_fixed(paste0("b", 1:4), rep(0, 4), rep(2, 4))
  thr <- c(M1 = 1)
  expect_error(sw_flag_iqr(), "`cells`")
  expect_error(
    sw_flag_iqr(
      cells,
      thresholds = thr,
      iqr_multiplier = 1.5,
      min_cells_per_side = 1,
      quantile_type = 7
    ),
    "`markers`"
  )
  expect_error(
    sw_flag_iqr(
      cells,
      "M1",
      iqr_multiplier = 1.5,
      min_cells_per_side = 1,
      quantile_type = 7
    ),
    "`thresholds`"
  )
  expect_error(
    sw_flag_iqr(cells, "M1", thr, min_cells_per_side = 1, quantile_type = 7),
    "`iqr_multiplier`"
  )
  expect_error(
    sw_flag_iqr(cells, "M1", thr, iqr_multiplier = 1.5, quantile_type = 7),
    "`min_cells_per_side`"
  )
  expect_error(
    sw_flag_iqr(cells, "M1", thr, iqr_multiplier = 1.5, min_cells_per_side = 1),
    "`quantile_type`"
  )
})

test_that("bad inputs fail loudly", {
  cells <- make_fixed(paste0("b", 1:4), rep(0, 4), rep(2, 4))
  expect_error(
    sw_flag_iqr(
      cells,
      "M1",
      c(M2 = 1),
      iqr_multiplier = 1.5,
      min_cells_per_side = 1,
      quantile_type = 7
    ),
    "No threshold"
  )
  expect_error(
    sw_flag_iqr(
      cells,
      "M9",
      c(M9 = 1),
      iqr_multiplier = 1.5,
      min_cells_per_side = 1,
      quantile_type = 7
    ),
    "missing"
  )
  cells$M1[1] <- NA
  expect_error(flag_fixed(cells), "NA values")
  expect_error(
    flag_fixed(
      make_fixed(paste0("b", 1:4), rep(0, 4), rep(2, 4)),
      quantile_type = 10
    ),
    "1 to 9"
  )
})

test_that("negative and positive split at the threshold, with ties negative", {
  cells <- tibble::tibble(
    batch = rep(paste0("b", 1:4), each = 6),
    control = "c1",
    # Threshold 1: 0, 1 and 1 are negative; 1.01, 3 and 5 are positive.
    M1 = rep(c(0, 1, 1, 1.01, 3, 5), 4)
  )
  res <- flag_fixed(cells)
  b1 <- res[res$batch == "b1", ]
  expect_equal(b1$value[b1$metric == "neg_median"], 1)
  expect_equal(b1$n_cells[b1$metric == "neg_median"], 3L)
  expect_equal(b1$value[b1$metric == "pos_median"], 3)
  expect_equal(b1$n_cells[b1$metric == "pos_median"], 3L)
  expect_equal(b1$value[b1$metric == "pct_positive"], 50)
  expect_equal(b1$n_cells[b1$metric == "pct_positive"], 6L)
})

test_that("percent positive follows the share of events above threshold", {
  cells <- make_fixed(
    paste0("b", 1:4),
    rep(0, 4),
    rep(2, 4),
    n_neg = c(9, 3, 1, 5),
    n_pos = c(1, 1, 3, 15)
  )
  res <- flag_fixed(cells)
  pct <- res[res$metric == "pct_positive", ]
  expect_equal(pct$value, c(10, 25, 75, 75))
  expect_equal(pct$n_cells, c(10L, 4L, 4L, 20L))
})

test_that("a shifted batch is flagged on the right marker and metric only", {
  set.seed(42)
  cells <- simulate_cells(8)
  # Shift the CD3-positive cells of batch b05 upward by 1.5 units.
  shift <- cells$batch == "b05" & cells$CD3 > 1.5
  cells$CD3[shift] <- cells$CD3[shift] + 1.5
  res <- flag_sim(cells)
  expect_true(all(res$status == "ok"))
  hits <- res[res$flag, ]
  expect_equal(nrow(hits), 1L)
  expect_equal(hits$batch, "b05")
  expect_equal(hits$marker, "CD3")
  expect_equal(hits$metric, "pos_median")
  expect_equal(attr(res, "params")$iqr_multiplier, 1.5)
})

test_that("three batches are ineligible, not 'no flags'", {
  set.seed(7)
  res <- flag_sim(simulate_cells(3))
  expect_true(all(res$status == "ineligible_n_batches"))
  expect_true(all(is.na(res$flag)))
  expect_true(all(is.na(res$upper)))
  expect_false(anyNA(res$value))
})

test_that("a control missing from one batch makes its panels ineligible", {
  set.seed(3)
  cells <- simulate_cells(5)
  other <- simulate_cells(5)
  other$control <- "ref2"
  cells <- dplyr::bind_rows(cells, other[other$batch != "b03", ])
  res <- flag_sim(cells)
  ref2 <- res[res$control == "ref2", ]
  expect_true(all(ref2$status == "ineligible_missing_control"))
  expect_true(all(is.na(ref2$flag)))
  absent <- ref2[ref2$batch == "b03", ]
  expect_true(all(is.na(absent$value)))
  expect_true(all(absent$n_cells == 0L))
  expect_match(ref2$reason[1], "b03")
  expect_true(all(res$status[res$control == "ref"] == "ok"))
})

test_that("a min-cells shortfall hits only the panels on that side", {
  set.seed(11)
  cells <- simulate_cells(6)
  # Batch b02 keeps only 20 CD3-positive events.
  keep <- !(cells$batch == "b02" & cells$CD3 > 1.5) |
    seq_len(nrow(cells)) %in%
      which(cells$batch == "b02" & cells$CD3 > 1.5)[1:20]
  cells <- cells[keep, ]
  res <- flag_sim(cells)
  cd3 <- res[res$marker == "CD3", ]
  expect_true(all(
    cd3$status[cd3$metric == "pos_median"] == "ineligible_min_cells"
  ))
  expect_true(all(is.na(cd3$flag[cd3$metric == "pos_median"])))
  expect_match(cd3$reason[cd3$metric == "pos_median"][1], "b02")
  expect_true(all(cd3$status[cd3$metric != "pos_median"] == "ok"))
  expect_true(all(res$status[res$marker == "CD19"] == "ok"))
})

test_that("an empty side with min_cells_per_side = 0 is undefined, not 0", {
  cells <- make_fixed(
    paste0("b", 1:4),
    rep(0, 4),
    rep(2, 4),
    n_pos = c(0, 5, 5, 5)
  )
  res <- flag_fixed(cells, min_cells = 0)
  pos <- res[res$metric == "pos_median", ]
  expect_true(all(pos$status == "undefined"))
  expect_true(is.na(pos$value[pos$batch == "b1"]))
  expect_true(all(is.na(pos$flag)))
})

test_that("quantile_type changes the fence on a tiny example", {
  # Per-batch positive medians 1.5, 2.5, 3.5 and 10.5 (threshold 1).
  cells <- make_fixed(paste0("b", 1:4), rep(0, 4), c(1.5, 2.5, 3.5, 10.5))
  pos7 <- flag_fixed(cells, quantile_type = 7)
  pos7 <- pos7[pos7$metric == "pos_median", ]
  pos6 <- flag_fixed(cells, quantile_type = 6)
  pos6 <- pos6[pos6$metric == "pos_median", ]
  expect_equal(pos7$upper[1], 9.75)
  expect_equal(pos6$upper[1], 19.25)
  expect_equal(pos7$flag, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(pos6$flag, c(FALSE, FALSE, FALSE, FALSE))
})

test_that("a repeated control in one batch errors when sample_col is given", {
  cells <- make_fixed(paste0("b", 1:4), rep(0, 4), rep(2, 4))
  cells$sample <- paste0("s_", cells$batch)
  cells$sample[cells$batch == "b2"][1:3] <- "s_b2_rerun"
  expect_error(flag_fixed(cells, sample_col = "sample"), "c1 in b2")
  cells$sample <- paste0("s_", cells$batch)
  expect_s3_class(flag_fixed(cells, sample_col = "sample"), "sw_flag_iqr")
})

test_that("summaries count flags per batch and per marker", {
  set.seed(42)
  cells <- simulate_cells(8)
  shift <- cells$batch == "b05" & cells$CD3 > 1.5
  cells$CD3[shift] <- cells$CD3[shift] + 1.5
  res <- flag_sim(cells)
  by_batch <- sw_flag_iqr_summary(res, by = "batch")
  expect_equal(by_batch$n_flagged[by_batch$batch == "b05"], 1L)
  expect_equal(sum(by_batch$n_flagged), 1L)
  expect_equal(unique(by_batch$n_assessed), 6L)
  by_marker <- sw_flag_iqr_summary(res, by = "marker")
  expect_equal(by_marker$n_flagged[by_marker$marker == "CD3"], 1L)
  s <- summary(res)
  expect_named(s, c("by_batch", "by_marker", "status"))
  expect_equal(s$status$n, 48L)
})
