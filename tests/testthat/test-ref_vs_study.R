# tests/testthat/test-ref_vs_study.R
# Mock data: 3 batches, 5 donors each, one reference id. The reference is the
# same in b1 and b2 and shifted by `shift` in b3.

make_rvs <- function(shift = 3, n = 300, seed = 1, n_donors = 5) {
  set.seed(seed)
  batches <- c("b1", "b2", "b3")
  rows <- list()
  dat <- list()
  for (b in batches) {
    rid <- paste0("ref_", b)
    dat[[rid]] <- cbind(
      A = stats::rnorm(n) + if (b == "b3") shift else 0,
      B = stats::rnorm(n)
    )
    rows[[rid]] <- data.frame(
      id = rid,
      b = b,
      donor = NA,
      role = "anchor",
      ref = "R"
    )
    for (d in seq_len(n_donors)) {
      id <- paste0(b, "_d", d)
      dat[[id]] <- cbind(
        A = stats::rnorm(n, mean = d * 0.2),
        B = stats::rnorm(n, mean = d * 0.2)
      )
      rows[[id]] <- data.frame(
        id = id,
        b = b,
        donor = paste0("p", d),
        role = "none",
        ref = NA
      )
    }
  }
  man <- sw_sample_manifest(
    do.call(rbind, rows),
    sample_id = "id",
    batch = "b",
    donor = "donor",
    reference_role = "role",
    reference_id = "ref"
  )
  list(data = dat, manifest = man)
}

rvs <- function(x, ...) {
  args <- list(
    data = x$data,
    manifest = x$manifest,
    markers = c("A", "B"),
    n_events = 200,
    seed = 1,
    reference_roles = "anchor",
    min_donors_per_batch = 4,
    noise_floor = NULL,
    transform = "asinh",
    cofactor = 6000
  )
  args[names(list(...))] <- list(...)
  do.call(sw_ref_vs_study, args)
}

test_that("a planted reference shift gives a high ratio only for the shifted batch", {
  res <- rvs(make_rvs())
  a <- res[res$marker == "A", ]
  r12 <- a$ratio[a$batch_a == "b1" & a$batch_b == "b2"]
  r13 <- a$ratio[a$batch_a == "b1" & a$batch_b == "b3"]
  expect_lt(r12, 0.5)
  expect_gt(r13, 1)
  # marker B was not shifted
  b <- res[res$marker == "B" & res$batch_a == "b1" & res$batch_b == "b3", ]
  expect_lt(b$ratio, 0.5)
  expect_equal(nrow(res), 6)
  expect_true(all(res$status == "ok"))
  expect_equal(attr(res, "cofactor"), 6000)
})

test_that("thin batches are not_assessed", {
  res <- rvs(make_rvs(n_donors = 3))
  expect_true(all(res$status == "not_assessed"))
  expect_true(all(is.na(res$ratio)))
  expect_match(res$reason[[1]], "fewer than min_donors_per_batch")
})

test_that("a seed reproduces and the RNG state is restored", {
  x <- make_rvs()
  set.seed(99)
  before <- .Random.seed
  r1 <- rvs(x, seed = 5)
  expect_identical(.Random.seed, before)
  r2 <- rvs(x, seed = 5)
  r3 <- rvs(x, seed = 6)
  expect_identical(r1, r2)
  expect_false(isTRUE(all.equal(r1$ref_emd, r3$ref_emd)))
})

test_that("a noise floor is joined by reference and marker", {
  fl <- data.frame(
    reference_id = "R",
    marker = c("A", "B"),
    floor_emd = c(0.05, 0.06)
  )
  res <- rvs(make_rvs(), noise_floor = fl)
  expect_setequal(unique(res$floor_emd), c(0.05, 0.06))
  dup <- rbind(fl, fl)
  expect_error(rvs(make_rvs(), noise_floor = dup), "noise_floor")
})

test_that("required arguments and bad inputs stop", {
  expect_error(sw_ref_vs_study(), "Required argument")
  x <- make_rvs()
  expect_error(rvs(x, min_donors_per_batch = 1), "at least 2")
  expect_error(rvs(x, n_events = 1000), "fewer than `n_events`")
  expect_error(
    rvs(x, reference_roles = "technical_replicate"),
    "No sample has a role"
  )
  expect_error(rvs(x, cofactor = "a"), "cofactor")
  one <- x$manifest[x$manifest$batch == "b1", ]
  expect_error(
    rvs(list(data = x$data[one$sample_id], manifest = one)),
    "two or more batches"
  )
})

test_that("a noise floor that misses a result key stops and names it", {
  fl <- data.frame(reference_id = "R", marker = "A", floor_emd = 0.05)
  expect_error(rvs(make_rvs(), noise_floor = fl), "R/B")
  typo <- data.frame(
    reference_id = "Rx",
    marker = c("A", "B"),
    floor_emd = c(0.05, 0.06)
  )
  expect_error(rvs(make_rvs(), noise_floor = typo), "no row for")
})
