# tests/testthat/test-design_audit.R
# Unit tests for R/design_audit.R. Pure data-frame logic; no flow packages needed.

make_samples <- function(batch, condition, role = "none", ref = NA_character_) {
  tibble::tibble(
    sample_id = paste0("s", seq_along(batch)),
    batch = batch,
    condition = condition,
    reference_role = role,
    reference_id = ref
  )
}

audit <- function(s, ...) {
  sw_design_audit(
    s,
    minimum_batches = 3,
    minimum_batches_origin = "test",
    require_shared_reference = TRUE,
    protected = "condition",
    ...
  )
}

test_that("the floor arguments have no default", {
  s <- make_samples(c("a", "b"), c("x", "y"))
  expect_error(sw_design_audit(s), "no default")
  expect_error(sw_design_audit(s, minimum_batches = 2), "no default")
  expect_error(
    sw_design_audit(
      s,
      minimum_batches = "3",
      minimum_batches_origin = "t",
      require_shared_reference = TRUE
    ),
    "single number"
  )
  expect_error(
    sw_design_audit(
      s,
      minimum_batches = 3,
      minimum_batches_origin = "",
      require_shared_reference = TRUE
    ),
    "non-empty"
  )
})

test_that("missing columns and duplicate ids are rejected", {
  s <- make_samples(c("a", "b"), c("x", "y"))
  expect_error(audit(s[, -3]), "missing")
  s$sample_id <- "dup"
  expect_error(audit(s), "unique")
})

test_that("a shared anchor with overlapping biology is reference_supported", {
  s <- make_samples(
    batch = rep(c("b1", "b2"), each = 3),
    condition = rep(c("ctrl", "case", "ctrl"), 2),
    role = rep(c("anchor", "none", "none"), 2),
    ref = rep(c("r1", NA, NA), 2)
  )
  res <- audit(s)
  expect_equal(res$comparisons$support, "reference_supported")
  expect_equal(res$comparisons$shared_reference_ids[[1]], "r1")
})

test_that("a biological pseudo-control never counts as a measured reference", {
  s <- make_samples(
    batch = rep(c("b1", "b2"), each = 2),
    condition = rep(c("ctrl", "case"), 2),
    role = rep(c("biological_pseudo_control", "none"), 2),
    ref = rep(c("r1", NA), 2)
  )
  expect_equal(
    audit(s)$comparisons$support,
    "biologically_stratified_conditional"
  )
})

test_that("perfect confounding is non_identifiable even with a shared reference", {
  s <- make_samples(
    batch = rep(c("b1", "b2"), each = 2),
    condition = rep(c("ctrl", "case"), each = 2),
    role = rep(c("anchor", "none"), 2),
    ref = rep(c("r1", NA), 2)
  )
  res <- audit(s)
  expect_equal(res$comparisons$support, "non_identifiable")
  expect_match(res$comparisons$reason, "does not change that")
  expect_equal(res$confounded_factors, "condition")
  expect_false(res$eligibility$eligible[res$eligibility$name == "cycombine"])
  expect_false(res$eligibility$override_allowed[
    res$eligibility$name == "cycombine"
  ])
})

test_that("no reference and no protected factor is descriptive_only", {
  s <- make_samples(c("b1", "b2"), c("x", "x"))
  res <- sw_design_audit(
    s,
    minimum_batches = 2,
    minimum_batches_origin = "test",
    require_shared_reference = FALSE
  )
  expect_equal(res$comparisons$support, "descriptive_only")
})

test_that("the design floor reports every reason", {
  s <- make_samples(c("b1", "b2"), c("x", "y"))
  res <- audit(s)
  expect_true(res$below_design_floor)
  expect_match(res$design_floor_reason, "declared minimum is 3")
  expect_match(res$design_floor_reason, "No measured reference")
})

test_that("the IQR flag needs 4 batches and a reference in every batch", {
  s <- make_samples(
    batch = rep(c("b1", "b2", "b3"), each = 2),
    condition = rep(c("ctrl", "case"), 3),
    role = c("anchor", "none", "anchor", "none", "none", "none"),
    ref = c("r1", NA, "r1", NA, NA, NA)
  )
  el <- audit(s)$eligibility
  row <- el[el$name == "iqr_flag", ]
  expect_false(row$eligible)
  expect_match(row$reason_code, "batch_count_below_minimum")
  expect_match(row$reason_code, "missing_reference_in_some_batch")
  expect_true(el$eligible[el$name == "none"])
})

test_that("disconnected reference graph is reported", {
  s <- make_samples(
    batch = rep(c("b1", "b2", "b3", "b4"), each = 2),
    condition = rep(c("ctrl", "case"), 4),
    role = rep(c("anchor", "none"), 4),
    ref = rep(c("r1", NA, "r2", NA), each = 2)[c(1, 2, 3, 4, 5, 6, 7, 8)]
  )
  s$reference_id <- c("r1", NA, "r1", NA, "r2", NA, "r2", NA)
  res <- audit(s)
  expect_length(res$reference_components, 2)
  expect_true(any(grepl("disconnected", res$warnings)))
})

test_that("population-level pairs are classified separately", {
  s <- make_samples(
    batch = rep(c("b1", "b2"), each = 2),
    condition = rep(c("ctrl", "case"), 2),
    role = rep(c("anchor", "none"), 2),
    ref = rep(c("r1", NA), 2)
  )
  counts <- tibble::tibble(
    sample_id = rep(s$sample_id, 2),
    population = rep(c("T", "pDC"), each = 4),
    n = c(100, 100, 100, 100, 0, 0, 5, 5)
  )
  res <- audit(s, populations = c("T", "pDC"), population_counts = counts)
  expect_equal(nrow(res$comparisons), 3)
  pdc <- res$comparisons[res$comparisons$population %in% "pDC", ]
  expect_equal(pdc$support, "descriptive_only")
})

test_that("print method runs", {
  s <- make_samples(c("b1", "b2"), c("x", "y"))
  expect_output(print(audit(s)), "sw_design_audit")
})

test_that("reference columns come as a pair and roles are checked", {
  s <- make_samples(c("b1", "b2"), c("x", "y"))
  expect_error(audit(s[, names(s) != "reference_id"]), "together")
  s$reference_role <- "pseudo"
  expect_error(audit(s), "Unknown reference_role")
})

test_that("one batch and missing batch values are handled", {
  s <- make_samples(c("b1", "b1", NA), c("x", "y", "x"))
  res <- audit(s)
  expect_equal(res$n_batches, 1)
  expect_equal(nrow(res$comparisons), 0)
  expect_true(any(grepl("no value for 'batch'", res$warnings)))
  expect_true(res$below_design_floor)
})

test_that("a donor column that is all NA warns", {
  s <- make_samples(c("b1", "b2"), c("x", "y"))
  s$donor <- NA_character_
  expect_true(any(grepl("donor", audit(s, donor_col = "donor")$warnings)))
})

test_that("nested batches do not trigger a rank warning; non-nested ones warn", {
  s <- make_samples(
    batch = rep(c("b1", "b2", "b3"), each = 2),
    condition = rep(c("ctrl", "case"), 3)
  )
  s$study <- rep(c("s1", "s1", "s2"), each = 2)
  res <- audit(s, nested_within = "study")
  expect_true(res$design_full_rank)
  expect_false(any(grepl("not nested", res$warnings)))
  s$study[1] <- "s2"
  expect_true(any(grepl(
    "not nested",
    audit(s, nested_within = "study")$warnings
  )))
})

test_that("no disconnected-graph warning when no reference exists", {
  s <- make_samples(c("b1", "b2", "b3"), c("x", "y", "x"))
  expect_false(any(grepl("disconnected", audit(s)$warnings)))
})

test_that("the floor origin is stored", {
  s <- make_samples(c("b1", "b2"), c("x", "y"))
  expect_equal(audit(s)$minimum_batches_origin, "test")
})

# Calibration graph, anchor and paired-shift eligibility (C1) ----------------

# Three batches; b1 is the anchor. r1 links b1-b2, r2 links b2-b3 only.
make_chain <- function() {
  s <- make_samples(
    batch = rep(c("b1", "b2", "b3"), each = 2),
    condition = rep(c("ctrl", "case"), 3),
    role = rep(c("anchor", "none"), 3),
    ref = c("r1", NA, "r1", NA, "r2", NA)
  )
  s <- dplyr::bind_rows(
    s,
    tibble::tibble(
      sample_id = "s7",
      batch = "b2",
      condition = "ctrl",
      reference_role = "anchor",
      reference_id = "r2"
    )
  )
  s
}

test_that("without an anchor, paired_shift is ineligible and bridge is empty", {
  res <- audit(make_chain())
  row <- res$eligibility[res$eligibility$name == "paired_shift", ]
  expect_false(row$eligible)
  expect_equal(row$reason_code, "anchor_not_declared")
  expect_false(row$override_allowed)
  expect_equal(row$sw_function, "sw_corrector_fit(engine = 'paired_shift')")
  expect_equal(nrow(res$anchor_bridge), 0)
  expect_null(res$anchor)
  expect_output(print(res), "Anchor: none declared")
})

test_that("an anchor missing from the data is an error", {
  expect_error(audit(make_chain(), anchor = "b9"), "not present")
  expect_error(audit(make_chain(), anchor = c("b1", "b2")), "one batch")
})

test_that("a batch bridged only through another batch is indirect_only", {
  res <- audit(make_chain(), anchor = "b1")
  br <- res$anchor_bridge
  expect_equal(br$batch, c("b2", "b3"))
  expect_equal(br$n_shared_with_anchor, c(1L, 0L))
  expect_equal(br$shared_ids[[1]], "r1")
  expect_equal(br$bridged, c(TRUE, FALSE))
  expect_equal(br$indirect_only, c(FALSE, TRUE))
  expect_length(res$reference_components, 1)
  row <- res$eligibility[res$eligibility$name == "paired_shift", ]
  expect_false(row$eligible)
  expect_equal(row$reason_code, "batch_without_anchor_bridge")
  expect_match(row$reason, "b3")
  expect_match(row$reason, "mean of \\(batch - anchor\\)")
  expect_true(any(grepl(
    "b3 share no measured reference with anchor",
    res$warnings
  )))
  expect_output(print(res), "Anchor: b1; 1 of 2")
})

test_that("every batch sharing a reference with the anchor is eligible", {
  s <- make_chain()
  s$reference_id[s$sample_id == "s5"] <- "r1"
  res <- audit(s, anchor = "b1")
  expect_true(all(res$anchor_bridge$bridged))
  expect_false(any(res$anchor_bridge$indirect_only))
  row <- res$eligibility[res$eligibility$name == "paired_shift", ]
  expect_true(row$eligible)
  expect_equal(row$reason_code, "eligible")
})

test_that("an anchor with no references bridges nothing", {
  s <- make_chain()
  s$reference_role[s$batch == "b1"] <- "none"
  s$reference_id[s$batch == "b1"] <- NA
  res <- audit(s, anchor = "b1")
  expect_equal(res$anchor_bridge$bridged, c(FALSE, FALSE))
  expect_equal(res$anchor_bridge$indirect_only, c(FALSE, FALSE))
  expect_true(any(grepl("holds no measured reference", res$warnings)))
  row <- res$eligibility[res$eligibility$name == "paired_shift", ]
  expect_match(row$reason, "b2, b3")
})

test_that("paired_shift needs at least two batches", {
  s <- make_samples(c("b1", "b1"), c("x", "y"), role = "anchor", ref = "r1")
  res <- audit(s, anchor = "b1")
  row <- res$eligibility[res$eligibility$name == "paired_shift", ]
  expect_equal(row$reason_code, "batch_count_below_minimum")
  expect_true(row$override_allowed)
})

test_that("the reference graph keeps aliquots of one reference id apart", {
  s <- make_chain()
  s$donor <- c("dR1", "d1", "dR1", "d2", "dR2", "d3", "dR2")
  s$aliquot <- c("a1", NA, "a1", NA, "a3", NA, "a4")
  s <- dplyr::bind_rows(
    s,
    tibble::tibble(
      sample_id = "s8",
      batch = "b2",
      condition = "case",
      reference_role = "anchor",
      reference_id = "r1",
      donor = "dR1",
      aliquot = "a2"
    )
  )
  res <- audit(s, anchor = "b1", donor_col = "donor", aliquot_col = "aliquot")
  g <- res$reference_graph
  expect_named(g, c("batch", "reference_id", "donor", "aliquot", "n_samples"))
  b2r1 <- g[g$batch == "b2" & g$reference_id == "r1", ]
  expect_equal(b2r1$aliquot, c("a1", "a2"))
  expect_equal(b2r1$n_samples, c(1L, 1L))
  expect_equal(sum(g$n_samples), 5L)
  expect_equal(res$anchor_bridge$n_shared_with_anchor[1], 1L)
  apd <- res$aliquots_per_donor
  expect_equal(apd$n_aliquots[apd$donor == "dR1"], 2L)
  dpb <- res$donors_per_batch
  expect_equal(dpb$n_donors, c(2L, 3L, 2L))
  expect_equal(dpb$n_reference_ids, c(1L, 2L, 1L))
})

test_that("without donor or aliquot columns the graph holds NA, not a level", {
  res <- audit(make_chain())
  expect_true(all(is.na(res$reference_graph$donor)))
  expect_true(all(is.na(res$reference_graph$aliquot)))
  expect_true(all(is.na(res$donors_per_batch$n_donors)))
  expect_equal(nrow(res$aliquots_per_donor), 0)
  expect_true(is.na(res$n_studies))
  expect_true(is.na(res$n_runs))
})

test_that("two studies with runs are counted and checked", {
  s <- make_samples(
    batch = rep(c("b1", "b2", "b3", "b4"), each = 2),
    condition = rep(c("ctrl", "case"), 4)
  )
  s$study <- rep(c("S1", "S2"), each = 4)
  s$run <- paste0("R", rep(1:4, each = 2))
  res <- audit(s, study_col = "study", run_col = "run")
  expect_equal(res$n_studies, 2L)
  expect_equal(res$n_runs, 4L)
  expect_true(res$design_full_rank)
  expect_false(any(grepl("not nested", res$warnings)))
  s$study[1] <- "S2"
  s$run[3] <- "R1"
  w <- audit(s, study_col = "study", run_col = "run")$warnings
  expect_true(any(grepl("Batch\\(es\\) b1 span more than one 'study'", w)))
  expect_true(any(grepl("Run\\(s\\) R1 span more than one 'batch'", w)))
})

test_that("study_col and nested_within on the same column warn once", {
  s <- make_samples(rep(c("b1", "b2"), each = 2), rep(c("x", "y"), 2))
  s$study <- c("S1", "S2", "S2", "S2")
  w <- audit(s, study_col = "study", nested_within = "study")$warnings
  expect_equal(sum(grepl("not nested", w)), 1)
})

test_that("crossed factors are tabulated and confounding is flagged apart", {
  s <- make_samples(
    batch = rep(c("b1", "b2", "b3"), each = 2),
    condition = rep(c("ctrl", "case"), 3)
  )
  s$instrument <- rep(c("A", "B", "C"), each = 2)
  s$stain_lot <- c("L1", "L2", "L1", "L2", "L1", "L2")
  res <- audit(s, crossed = c("instrument", "stain_lot"))
  expect_equal(res$confounded_crossed, "instrument")
  expect_length(res$confounded_factors, 0)
  expect_true(res$eligibility$eligible[res$eligibility$name == "cycombine"])
  expect_named(res$crossed_contingency, c("factor", "batch", "level", "n"))
  expect_equal(sum(res$crossed_contingency$n), 12L)
  expect_true(any(grepl("Crossed factor\\(s\\) instrument", res$warnings)))
  expect_error(audit(s, crossed = "lot"), "missing")
})

test_that("a manifest is accepted and fills its standard columns", {
  df <- tibble::tibble(
    id = paste0("s", 1:4),
    b = c("b1", "b1", "b2", "b2"),
    st = c("S1", "S1", "S1", "S1"),
    d = c("dR", "d1", "dR", "d2"),
    al = c("a1", NA, "a2", NA),
    cond = c("ctrl", "case", "ctrl", "case"),
    role = c("anchor", "none", "anchor", "none"),
    ref = c("r1", NA, "r1", NA)
  )
  m <- sw_sample_manifest(
    df,
    sample_id = "id",
    batch = "b",
    study = "st",
    donor = "d",
    aliquot = "al",
    condition = "cond",
    reference_role = "role",
    reference_id = "ref"
  )
  res <- audit(m, anchor = "b1")
  expect_equal(res$n_studies, 1L)
  expect_equal(res$reference_graph$donor, c("dR", "dR"))
  expect_equal(res$reference_graph$aliquot, c("a1", "a2"))
  expect_true(res$eligibility$eligible[res$eligibility$name == "paired_shift"])
})

test_that("synthetic_reference never bridges and eligibility is unchanged", {
  base <- make_samples(
    c("a", "a", "b", "b", "c", "c"),
    rep(c("x", "y"), 3),
    role = c("technical_replicate", "none", "technical_replicate", "none", "none", "none"),
    ref = c("R1", NA, "R1", NA, NA, NA)
  )
  base$source_batch <- NA_character_
  with_syn <- base
  with_syn$reference_role[c(2, 6)] <- "synthetic_reference"
  with_syn$reference_id[c(2, 6)] <- "S1"
  with_syn$source_batch[c(2, 6)] <- "a"

  a0 <- audit(base, anchor = "a")
  a1 <- audit(with_syn, anchor = "a")
  # Synthetic S1 is shared by batches a and c but must not connect them.
  expect_false(a1$anchor_bridge$bridged[a1$anchor_bridge$batch == "c"])
  expect_equal(a1$reference_components, a0$reference_components)
  expect_equal(a1$reference_coverage, a0$reference_coverage)
  expect_equal(a1$reference_graph, a0$reference_graph)
  expect_equal(a1$eligibility, a0$eligibility)
  expect_equal(nrow(a0$synthetic_references), 0)
  expect_equal(a1$synthetic_references$sample_id, c("s2", "s6"))
  expect_equal(a1$synthetic_references$source_batch, c("a", "a"))
})

test_that("synthetic_reference needs a source_batch", {
  s <- make_samples(c("a", "b", "c"), c("x", "y", "x"),
                    role = c("synthetic_reference", "none", "none"),
                    ref = c("S1", NA, NA))
  expect_error(audit(s), "source_batch")
  s$source_batch <- NA_character_
  expect_error(audit(s), "source_batch")
})
