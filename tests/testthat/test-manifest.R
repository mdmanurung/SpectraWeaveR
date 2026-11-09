# tests/testthat/test-manifest.R
# Unit tests for R/manifest.R. Pure data-frame logic.

df <- data.frame(
  id = c("s1", "s2", "s3", "s4"),
  b = c("b1", "b1", "b2", "b2"),
  cohort = c("A", "A", "B", "B"),
  donor_id = c("d1", "d2", "d1", "d3"),
  extra = 1:4
)

test_that("columns are standardised and extras kept", {
  m <- sw_sample_manifest(df, sample_id = "id", batch = "b", study = "cohort",
                          donor = "donor_id")
  expect_s3_class(m, "sw_manifest")
  expect_equal(names(m), c("sample_id", "batch", "study", "donor", "extra"))
  expect_equal(m$sample_id, df$id)
})

test_that("sample_id and batch have no default", {
  expect_error(sw_sample_manifest(df), "required")
  expect_error(sw_sample_manifest(df, sample_id = "id"), "required")
})

test_that("bad column names and inputs are rejected", {
  expect_error(sw_sample_manifest(df, sample_id = "nope", batch = "b"), "not found")
  expect_error(sw_sample_manifest(df, sample_id = c("id", "b"), batch = "b"),
               "single string")
  expect_error(sw_sample_manifest(1, sample_id = "id", batch = "b"), "data frame")
})

test_that("duplicate or missing ids and missing batch stop", {
  d <- df
  d$id[2] <- "s1"
  expect_error(sw_sample_manifest(d, "id", "b"), "unique")
  d <- df
  d$id[1] <- NA
  expect_error(sw_sample_manifest(d, "id", "b"), "missing")
  d <- df
  d$b[1] <- NA
  expect_error(sw_sample_manifest(d, "id", "b"), "'batch' has missing")
})

test_that("batches must nest in studies and runs in batches", {
  d <- df
  d$cohort[1] <- "B"
  expect_error(sw_sample_manifest(d, "id", "b", study = "cohort"), "more than one study")
  d <- df
  d$r <- c("r1", "r2", "r1", "r3")
  expect_error(sw_sample_manifest(d, "id", "b", run = "r"), "more than one batch")
})

test_that("reference columns come as a pair, with valid roles and ids", {
  d <- df
  d$role <- c("anchor", "none", "anchor", "none")
  d$ref <- c("r1", NA, "r1", NA)
  expect_error(sw_sample_manifest(d, "id", "b", reference_role = "role"), "together")
  m <- sw_sample_manifest(d, "id", "b", reference_role = "role", reference_id = "ref")
  expect_true(all(c("reference_role", "reference_id") %in% names(m)))
  d$role[2] <- "pseudo"
  expect_error(
    sw_sample_manifest(d, "id", "b", reference_role = "role", reference_id = "ref"),
    "Unknown reference_role"
  )
  d$role[2] <- "none"
  d$ref[1] <- NA
  expect_error(
    sw_sample_manifest(d, "id", "b", reference_role = "role", reference_id = "ref"),
    "need a 'reference_id'"
  )
})

test_that("an unmapped column that shares a standard name is refused", {
  d <- df
  d$study <- "x"
  expect_error(sw_sample_manifest(d, "id", "b", donor = "donor_id"), "standard name")
  d2 <- df
  d2$donor <- "x"
  expect_error(sw_sample_manifest(d2, "id", "b"), "standard name")
})

test_that("old sample_meta converts and round-trips", {
  meta <- data.frame(sample = c("a", "b"), batch = c("b1", "b2"),
                     condition = c("ctrl", "case"))
  m <- sw_manifest_from_sample_meta(meta, condition = "condition")
  expect_equal(m$sample_id, meta$sample)
  expect_equal(m$condition, meta$condition)
  expect_error(sw_manifest_from_sample_meta(data.frame(x = 1)), "'sample' and 'batch'")
})

test_that("a manifest feeds the design audit directly", {
  m <- sw_sample_manifest(df, "id", "b", study = "cohort", donor = "donor_id")
  res <- sw_design_audit(
    m, minimum_batches = 2, minimum_batches_origin = "test",
    require_shared_reference = FALSE
  )
  expect_equal(res$n_batches, 2)
})

test_that("print method runs", {
  m <- sw_sample_manifest(df, "id", "b", study = "cohort")
  expect_output(print(m), "sw_manifest")
})
