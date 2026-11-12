# tests/testthat/test-ref_consistency.R
# sw_ref_consistency() on mock data. ref1 sits in 3 batches and is shifted in
# marker A in b2; ref2 sits in 2 batches only.

make_ref <- function(shift = 1, n = 400, seed = 1) {
  set.seed(seed)
  rows <- data.frame(
    id = c("r1_b1", "r1_b2", "r1_b3", "r2_b1", "r2_b2"),
    batch = c("b1", "b2", "b3", "b1", "b2"),
    ref = c("ref1", "ref1", "ref1", "ref2", "ref2")
  )
  d <- lapply(seq_len(nrow(rows)), function(i) {
    cbind(
      A = rnorm(n) + if (rows$id[i] == "r1_b2") shift else 0,
      B = rnorm(n)
    )
  })
  names(d) <- rows$id
  man <- data.frame(
    sample_id = rows$id,
    batch = rows$batch,
    reference_role = "technical_replicate",
    reference_id = rows$ref
  )
  list(d = d, man = man)
}

cons <- function(x, ...) {
  args <- list(
    data = x$d,
    manifest = x$man,
    markers = c("A", "B"),
    n_events = 300,
    seed = 1,
    order = NULL,
    min_batches_per_reference = 3,
    min_events = 100,
    roles = "technical_replicate",
    transform = "asinh",
    cofactor = 6000
  )
  args[names(list(...))] <- list(...)
  do.call(sw_ref_consistency, args)
}

pick <- function(res, ref, check, marker = NA, batch = NA) {
  r <- res[res$reference_id == ref & res$check == check, ]
  if (!is.na(marker)) {
    r <- r[r$marker == marker, ]
  }
  if (!is.na(batch)) {
    r <- r[r$batch == batch, ]
  }
  r
}

test_that("required arguments have no default", {
  x <- make_ref()
  expect_error(sw_ref_consistency(), "no default")
  expect_error(
    sw_ref_consistency(
      x$d,
      x$man,
      c("A", "B"),
      n_events = 100,
      seed = 1,
      order = NULL,
      min_batches_per_reference = 3,
      min_events = 50,
      roles = "technical_replicate",
      transform = "asinh"
    ),
    "cofactor"
  )
  expect_error(cons(x, transform = 1), "transform")
  expect_error(cons(x, cofactor = -1), "cofactor")
})

test_that("a planted shift ranks first for its batch and marker only", {
  res <- cons(make_ref(shift = 1))
  expect_s3_class(res, "sw_ref_consistency")
  rest <- pick(res, "ref1", "rest_emd", "A")
  expect_equal(rest$batch[which.max(rest$value)], "b2")
  # With 3 batches the robust z is capped near 2; the shifted batch is the top.
  for (chk in c("rest_emd_z", "median_shift_z")) {
    z <- pick(res, "ref1", chk, "A")
    expect_equal(z$batch[which.max(z$value)], "b2")
    expect_gt(max(z$value), 1.5)
  }
  expect_equal(
    pick(res, "ref1", "pair_emd", "A", "b1 ~ b2")$value,
    1,
    tolerance = 0.2
  )
  expect_lt(pick(res, "ref1", "pair_emd", "A", "b1 ~ b3")$value, 0.3)
  expect_lt(max(pick(res, "ref1", "rest_emd", "B")$value), 0.3)
  # No per-reference verdict column.
  expect_false(any(c("pass", "consistent", "verdict", "flag") %in% names(res)))
})

test_that("the seed reproduces the table and the caller's RNG state is restored", {
  x <- make_ref()
  set.seed(99)
  before <- .Random.seed
  a <- cons(x)
  expect_identical(.Random.seed, before)
  b <- cons(x)
  expect_equal(a, b)
  expect_false(isTRUE(all.equal(a$value, cons(x, seed = 2)$value)))
})

test_that("thin rows are not_assessed with a reason", {
  res <- cons(make_ref())
  # ref2 has 2 batches: pair and rest EMD are filled, z rows are not.
  expect_equal(pick(res, "ref2", "pair_emd", "A")$status, "assessed")
  expect_equal(pick(res, "ref2", "rest_emd", "A", "b1")$status, "assessed")
  z <- pick(res, "ref2", "rest_emd_z", "A")
  expect_true(all(z$status == "not_assessed" & is.na(z$value)))
  expect_match(z$reason[1], "min_batches_per_reference")
  # A file under min_events leaves its batch unscored but visible.
  x <- make_ref()
  x$d$r1_b3 <- x$d$r1_b3[1:60, ]
  thin <- cons(x)
  row <- pick(thin, "ref1", "rest_emd", "A", "b3")
  expect_equal(row$status, "not_assessed")
  expect_match(row$reason, "min_events")
  expect_equal(pick(thin, "ref1", "count_ratio", batch = "b3")$value, 60 / 400)
  expect_true(all(is.na(pick(thin, "ref1", "count_ratio")$marker)))
  # Two of three batches thin: nothing is left to compare.
  x$d$r1_b2 <- x$d$r1_b2[1:60, ]
  one <- cons(x)
  expect_true(all(pick(one, "ref1", "rest_emd")$status == "not_assessed"))
  expect_equal(nrow(pick(one, "ref1", "pair_emd")), 0)
})

test_that("order_rho appears with order and settings are recorded", {
  set.seed(3)
  x <- make_ref()
  # Four files of ref1 in b1 with a drift of marker A along the order.
  extra <- paste0("r1_b1_", 1:3)
  for (i in seq_along(extra)) {
    x$d[[extra[i]]] <- cbind(A = rnorm(400) + i, B = rnorm(400))
  }
  x$man <- rbind(
    x$man,
    data.frame(
      sample_id = extra,
      batch = "b1",
      reference_role = "technical_replicate",
      reference_id = "ref1"
    )
  )
  ord <- stats::setNames(seq_along(x$d), names(x$d))
  ord[c("r1_b1", extra)] <- 1:4
  res <- cons(x, order = ord)
  expect_equal(pick(res, "ref1", "order_rho", "A", "b1")$value, 1)
  expect_equal(nrow(pick(cons(x), "ref1", "order_rho")), 0)
  st <- attr(res, "settings")
  expect_equal(st$transform, "asinh")
  expect_equal(st$cofactor, 6000)
  expect_equal(st$min_batches_per_reference, 3)
})
