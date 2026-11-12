# tests/testthat/test-synref_io.R
# sw_synref_save(), sw_synref_load() and sw_synref_as_manifest().

make_ref <- function() {
  set.seed(1)
  d <- list(
    f1 = cbind(M1 = stats::rnorm(100), M2 = stats::rnorm(100)),
    f2 = cbind(M1 = stats::rnorm(100), M2 = stats::rnorm(100))
  )
  man <- data.frame(
    sample_id = c("f1", "f2"),
    batch = c("b1", "b2"),
    study = "s1",
    donor = c("d1", "d2")
  )
  ref <- sw_synref_downsample(
    d,
    man,
    groups = list(r1 = "f1", r2 = "f2"),
    strata = "file",
    labels = NULL,
    n_events = 20,
    n_draws = 2,
    seed = 4
  )
  list(ref = ref, man = man)
}

test_that("save then load returns an identical object", {
  x <- make_ref()
  dir <- file.path(withr::local_tempdir(), "ref")
  sw_synref_save(x$ref, dir)
  expect_true(all(
    c("synref.rds", "provenance.csv", "digest.txt") %in% list.files(dir)
  ))
  back <- sw_synref_load(dir)
  expect_identical(back, x$ref)
})

test_that("save never overwrites an existing directory", {
  x <- make_ref()
  dir <- withr::local_tempdir()
  expect_error(sw_synref_save(x$ref, dir), "never overwritten")
  new <- file.path(dir, "ref")
  sw_synref_save(x$ref, new)
  expect_error(sw_synref_save(x$ref, new), "never overwritten")
})

test_that("load stops when the object was changed after the build", {
  x <- make_ref()
  dir <- file.path(withr::local_tempdir(), "ref")
  sw_synref_save(x$ref, dir)
  tampered <- readRDS(file.path(dir, "synref.rds"))
  tampered$events[[1]][1, 1] <- tampered$events[[1]][1, 1] + 1
  saveRDS(tampered, file.path(dir, "synref.rds"))
  expect_error(sw_synref_load(dir), "Digest mismatch")
  expect_error(sw_synref_load(file.path(dir, "none")), "No synref.rds")
})

test_that("manifest rows carry the role, the source batch and the study", {
  x <- make_ref()
  rows <- sw_synref_as_manifest(x$ref, x$man)
  expect_equal(nrow(rows), 4)
  expect_true(all(rows$reference_role == "synthetic_reference"))
  expect_equal(rows$source_batch, rows$batch)
  expect_equal(sort(unique(rows$batch)), c("b1", "b2"))
  expect_equal(rows$sample_id, rows$reference_id)
  expect_true(all(rows$study == "s1"))
  expect_false(any(rows$sample_id %in% x$man$sample_id))
})

test_that("manifest rows pass sw_sample_manifest and the design audit", {
  x <- make_ref()
  rows <- sw_synref_as_manifest(x$ref, x$man)
  study <- data.frame(
    sample_id = x$man$sample_id,
    batch = x$man$batch,
    study = "s1",
    donor = x$man$donor,
    reference_role = NA_character_,
    reference_id = NA_character_,
    source_batch = NA_character_
  )
  both <- dplyr::bind_rows(study, as.data.frame(rows))
  both$donor[is.na(both$donor)] <- "synthetic"
  m <- sw_sample_manifest(
    both,
    sample_id = "sample_id",
    batch = "batch",
    study = "study",
    donor = "donor",
    reference_role = "reference_role",
    reference_id = "reference_id",
    source_batch = "source_batch"
  )
  expect_s3_class(m, "data.frame")
  audit <- sw_design_audit(
    m,
    minimum_batches = 2,
    minimum_batches_origin = "test",
    require_shared_reference = FALSE
  )
  expect_equal(sort(audit$synthetic_references$sample_id), sort(rows$sample_id))
  nobatch <- both
  nobatch$source_batch <- NA_character_
  expect_error(
    sw_sample_manifest(
      nobatch,
      sample_id = "sample_id",
      batch = "batch",
      study = "study",
      donor = "donor",
      reference_role = "reference_role",
      reference_id = "reference_id",
      source_batch = "source_batch"
    ),
    "source_batch"
  )
})

test_that("missing arguments stop", {
  x <- make_ref()
  expect_error(sw_synref_as_manifest(x$ref), "'manifest'")
  expect_error(sw_synref_save(x$ref), "'dir'")
  expect_error(sw_synref_load(), "'dir'")
  expect_error(sw_synref_as_manifest(list(), x$man), "sw_synref")
})
