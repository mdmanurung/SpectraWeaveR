# tests/testthat/test-ref_composition.R

make_rc <- function(p_b3 = 0.3, n = 2000, seed = 1) {
  set.seed(seed)
  man <- data.frame(
    id = c("r1", "r2", "r3"),
    b = c("b1", "b2", "b3"),
    role = "anchor",
    ref = "R"
  )
  man <- sw_sample_manifest(
    man,
    sample_id = "id",
    batch = "b",
    reference_role = "role",
    reference_id = "ref"
  )
  p <- c(b1 = 0.2, b2 = 0.2, b3 = p_b3)
  labs <- lapply(p, function(pb) {
    sample(c("B", "T"), n, replace = TRUE, prob = c(pb, 1 - pb))
  })
  names(labs) <- c("r1", "r2", "r3")
  list(labels = labs, manifest = man)
}

rcomp <- function(x, ...) {
  args <- list(
    labels = x$labels,
    manifest = x$manifest,
    reference_roles = "anchor",
    min_population_events = 50,
    label_source = "gate"
  )
  args[names(list(...))] <- list(...)
  do.call(sw_ref_composition, args)
}

test_that("a planted share shift is found and the floor matches the multinomial formula", {
  res <- rcomp(make_rc())
  b <- res[res$population == "B", ]
  expect_equal(b$batch, c("b1", "b2", "b3"))
  expect_gt(b$share[b$batch == "b3"], 0.27)
  expect_lt(b$share[b$batch == "b1"], 0.23)
  expect_equal(b$floor, sqrt(b$share * (1 - b$share) / b$n_events))
  expect_equal(sum(res$share[res$batch == "b1"]), 1)
  # the shift is larger than the floor
  expect_gt(abs(b$share[3] - b$share[1]), 3 * max(b$floor))
  flat <- rcomp(make_rc(p_b3 = 0.2))
  expect_lt(flat$cv_across_batches[1], b$cv_across_batches[1])
  expect_equal(attr(res, "label_source"), "gate")
})

test_that("a rare population is not_assessed and a missing one counts as zero", {
  x <- make_rc()
  x$labels$r1[1:3] <- "rare"
  res <- rcomp(x)
  rare <- res[res$population == "rare", ]
  expect_equal(nrow(rare), 3)
  expect_true(all(rare$status == "not_assessed"))
  expect_true(all(is.na(rare$cv_across_batches)))
  expect_equal(rare$n_population_events, c(3L, 0L, 0L))
})

test_that("it is deterministic and needs no RNG", {
  x <- make_rc()
  set.seed(3)
  before <- .Random.seed
  r1 <- rcomp(x)
  expect_identical(.Random.seed, before)
  expect_identical(r1, rcomp(x))
})

test_that("required arguments and bad inputs stop", {
  expect_error(sw_ref_composition(), "Required argument")
  x <- make_rc()
  expect_error(rcomp(x, min_population_events = 0), "whole number")
  expect_error(rcomp(x, reference_roles = "none"), "No labelled sample")
  x$labels$r2[1] <- NA
  expect_error(rcomp(x), "must not contain NA")
})

test_that("a batch with a zero-event file stays in the grid as not_assessed", {
  x <- make_rc()
  x$labels[["r3"]] <- character()
  res <- rcomp(x)
  empty <- res[res$batch == "b3", ]
  expect_true(nrow(empty) > 0)
  expect_true(all(empty$n_events == 0))
  expect_true(all(is.na(empty$share)))
  expect_true(all(empty$status == "not_assessed"))
  b <- res[res$population == "B" & res$batch != "b3", ]
  expect_equal(b$cv_across_batches[[1]], stats::sd(b$share) / mean(b$share))
})
