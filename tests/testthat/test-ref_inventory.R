# tests/testthat/test-ref_inventory.R
# sw_ref_inventory() on mock data.

make_inv <- function() {
  set.seed(1)
  d <- list(
    r1 = cbind(A = rnorm(100), B = rnorm(100)),
    r2 = cbind(A = rnorm(40), B = rnorm(40)),
    r3 = cbind(A = rnorm(100), B = NA_real_),
    s1 = cbind(A = rnorm(100), B = rnorm(100))
  )
  man <- data.frame(
    sample_id = names(d),
    batch = c("b1", "b2", "b3", "b1"),
    reference_role = c(rep("technical_replicate", 3), "study"),
    reference_id = c("ref1", "ref1", "ref1", NA)
  )
  lab <- list(
    r1 = rep(c("T", "B"), c(70, 30)),
    r2 = rep("T", 40),
    r3 = rep(c("T", "B"), c(50, 50))
  )
  list(d = d, man = man, lab = lab)
}

inv <- function(x, ...) {
  args <- list(
    data = x$d,
    manifest = x$man,
    markers = c("A", "B"),
    labels = x$lab,
    min_events = 50,
    roles = "technical_replicate"
  )
  args[names(list(...))] <- list(...)
  do.call(sw_ref_inventory, args)
}

test_that("required arguments have no default", {
  x <- make_inv()
  expect_error(sw_ref_inventory(), "no default")
  expect_error(
    sw_ref_inventory(
      x$d,
      x$man,
      markers = "A",
      min_events = 5,
      roles = "anchor"
    ),
    "labels"
  )
  expect_error(inv(x, min_events = 0), "min_events")
})

test_that("counts, markers and under_min match the planted data", {
  x <- make_inv()
  res <- inv(x)
  expect_s3_class(res, "sw_ref_inventory")
  # Only measured reference samples are listed.
  expect_setequal(res$sample_id, c("r1", "r2", "r3"))
  mk <- res[res$row_type == "marker", ]
  expect_equal(unique(mk$n_events[mk$sample_id == "r2"]), 40L)
  expect_true(all(mk$under_min[mk$sample_id == "r2"]))
  expect_false(any(mk$under_min[mk$sample_id == "r1"]))
  # The all-NA marker is not present.
  expect_false(mk$present[mk$sample_id == "r3" & mk$marker == "B"])
  expect_true(mk$present[mk$sample_id == "r3" & mk$marker == "A"])
  # A population a sample lacks has a count of 0.
  pp <- res[res$row_type == "population", ]
  expect_equal(
    pp$n_population_events[pp$sample_id == "r2" & pp$population == "B"],
    0L
  )
  expect_equal(
    pp$n_population_events[pp$sample_id == "r1" & pp$population == "T"],
    70L
  )
  expect_true(pp$under_min[pp$sample_id == "r1" & pp$population == "B"])
  expect_false(any(c("pass", "ok", "verdict") %in% names(res)))
})

test_that("labels = NULL skips population rows; bad inputs stop", {
  x <- make_inv()
  res <- inv(x, labels = NULL)
  expect_false("population" %in% res$row_type)
  bad <- x
  bad$lab$r1 <- bad$lab$r1[-1]
  expect_error(inv(bad), "labels")
  expect_error(inv(x, roles = "nothing"), "No manifest row")
  gone <- x
  gone$d$r3 <- NULL
  expect_error(inv(gone), "not in `data`")
})
