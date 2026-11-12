# tests/testthat/test-batch_decide.R
# Unit tests for R/batch_decide.R. Pure data-frame logic; no flow packages needed.

thresholds <- c(
  "min_coverage",
  "min_emd_reduction",
  "max_mad_compression_fraction",
  "max_joint_structure_change",
  "no_effect_emd"
)

make_policy <- function(
  calibrated = FALSE,
  require_ref = TRUE,
  contrast = NULL,
  min_donors = NULL,
  spread_ratio = if (is.null(contrast)) NULL else 2
) {
  origin <- if (calibrated) {
    "calibrated on spike-in study"
  } else {
    "uncalibrated engineering default"
  }
  all_thresholds <- c(thresholds, if (!is.null(contrast)) "max_contrast_change")
  artifacts <- if (calibrated) {
    stats::setNames(paste0("calib/", all_thresholds, ".json"), all_thresholds)
  } else {
    character(0)
  }
  sw_batch_policy(
    min_coverage = 0.6,
    min_emd_reduction = 0.2,
    max_mad_compression_fraction = 0.25,
    max_joint_structure_change = 0.1,
    no_effect_emd = 0.05,
    require_reference_supported = require_ref,
    max_contrast_change = contrast,
    min_donors_per_condition = min_donors,
    min_contrast_spread_ratio = spread_ratio,
    threshold_origins = stats::setNames(
      rep(origin, length(all_thresholds)),
      all_thresholds
    ),
    calibration_artifacts = artifacts
  )
}

# Hand-built audit so this file does not depend on design_audit.R.
make_design <- function(
  support = "reference_supported",
  below_floor = FALSE,
  confounded = character(0)
) {
  structure(
    list(
      batch_col = "batch",
      confounded_factors = confounded,
      below_design_floor = below_floor,
      design_floor_reason = if (below_floor) {
        "2 batch(es) present; the declared minimum is 3."
      } else {
        NA_character_
      },
      comparisons = tibble::tibble(
        comparison_id = paste0("c", seq_along(support)),
        population = NA_character_,
        support = support
      )
    ),
    class = "sw_design_audit"
  )
}

# One good candidate "cc" next to the identity "none", both on comparison c1.
make_candidates <- function(...) {
  base <- tibble::tibble(
    candidate_id = c("none", "cc"),
    status = "ok",
    emd_median = c(0.8, 0.3),
    emd_worst_marker_median = c(1.2, 0.4),
    emd_relative_reduction = c(0, 0.6),
    coverage = c(1, 1),
    mad_compression_fraction = c(0, 0.05),
    joint_structure_distance_median = c(0, 0.02),
    evidence_comparison_ids = list("c1", "c1")
  )
  changes <- list(...)
  for (nm in names(changes)) {
    base[[nm]][2] <- changes[[nm]]
  }
  base
}

decide <- function(
  candidates = make_candidates(),
  design = make_design(),
  policy = make_policy(calibrated = TRUE),
  label_source = "manual_gating",
  contrast_records = NULL
) {
  sw_batch_decide(
    design,
    candidates,
    policy,
    identity_id = "none",
    label_source = label_source,
    contrast_records = contrast_records
  )
}

row_of <- function(res, id = "cc") {
  res$assessments[res$assessments$candidate_id == id, ]
}

state_of <- function(res, id = "cc") {
  res$assessments$state[res$assessments$candidate_id == id]
}

test_that("every policy argument has no default", {
  expect_error(sw_batch_policy(), "min_coverage")
  expect_error(
    sw_batch_policy(
      min_coverage = 0.6,
      min_emd_reduction = 0.2,
      max_mad_compression_fraction = 0.25,
      max_joint_structure_change = 0.1,
      no_effect_emd = 0.05,
      require_reference_supported = TRUE,
      max_contrast_change = NULL,
      min_donors_per_condition = NULL,
      min_contrast_spread_ratio = NULL,
      calibration_artifacts = character(0)
    ),
    "threshold_origins"
  )
})

test_that("threshold_origins must name all five thresholds", {
  expect_error(
    sw_batch_policy(
      min_coverage = 0.6,
      min_emd_reduction = 0.2,
      max_mad_compression_fraction = 0.25,
      max_joint_structure_change = 0.1,
      no_effect_emd = 0.05,
      require_reference_supported = TRUE,
      max_contrast_change = NULL,
      min_donors_per_condition = NULL,
      min_contrast_spread_ratio = NULL,
      threshold_origins = c(min_coverage = "x"),
      calibration_artifacts = character(0)
    ),
    "one non-empty origin"
  )
})

test_that("calibration needs a non-uncalibrated origin and an artifact for every threshold", {
  expect_true(make_policy(calibrated = FALSE)$any_threshold_uncalibrated)
  expect_false(make_policy(calibrated = TRUE)$any_threshold_uncalibrated)

  # A clean origin without its artifact still counts as uncalibrated.
  p <- sw_batch_policy(
    min_coverage = 0.6,
    min_emd_reduction = 0.2,
    max_mad_compression_fraction = 0.25,
    max_joint_structure_change = 0.1,
    no_effect_emd = 0.05,
    require_reference_supported = TRUE,
    max_contrast_change = NULL,
    min_donors_per_condition = NULL,
    min_contrast_spread_ratio = NULL,
    threshold_origins = stats::setNames(rep("from study X", 5), thresholds),
    calibration_artifacts = stats::setNames(rep("a.json", 4), thresholds[1:4])
  )
  expect_true(p$any_threshold_uncalibrated)
  expect_false(p$settings$calibrated[p$settings$threshold == "no_effect_emd"])
})

test_that("decide requires identity_id and label_source, and checks label_source", {
  expect_error(
    sw_batch_decide(make_design(), make_candidates(), make_policy()),
    "identity_id"
  )
  expect_error(
    sw_batch_decide(
      make_design(),
      make_candidates(),
      make_policy(),
      identity_id = "none"
    ),
    "label_source"
  )
  expect_error(decide(label_source = "guess"), "must be one of")
})

test_that("the identity candidate must be present", {
  cands <- make_candidates()[2, ]
  expect_error(decide(cands), "always be assessed")
})

test_that("a calibrated policy with every gate passed recommends the candidate", {
  res <- decide()
  expect_equal(state_of(res), "recommended_under_policy")
  expect_equal(state_of(res, "none"), "eligible_for_review")
  expect_equal(res$recommended_candidate_id, "cc")
  expect_true(is.na(res$recommendation_withheld_reason))
})

test_that("an uncalibrated policy never reaches recommended_under_policy", {
  res <- decide(policy = make_policy(calibrated = FALSE))
  expect_equal(state_of(res), "eligible_for_review")
  expect_match(row_of(res)$reason, "uncalibrated")
  expect_true(is.na(res$recommended_candidate_id))
  expect_true(any(grepl("uncalibrated", res$warnings)))
})

test_that("clustering labels cap a passing candidate at eligible_for_review", {
  res <- decide(label_source = "clustering")
  expect_equal(state_of(res), "eligible_for_review")
  expect_match(row_of(res)$reason, "clustering")
})

test_that("identity is not_needed when the baseline shows no batch effect", {
  cands <- make_candidates()
  cands$emd_worst_marker_median[1] <- 0.01
  res <- decide(cands)
  expect_equal(state_of(res, "none"), "not_needed")
  expect_match(
    res$assessments$reason[res$assessments$candidate_id == "none"],
    "0.01"
  )
})

test_that("a failed backend is execution_failed with its own reason", {
  cands <- make_candidates(status = "failed")
  cands$status_reason <- c(NA, "ComBat did not converge")
  res <- decide(cands)
  expect_equal(state_of(res), "execution_failed")
  expect_equal(row_of(res)$reason, "ComBat did not converge")
})

test_that("ineligible and not_available give insufficient_evidence", {
  expect_equal(
    state_of(decide(make_candidates(status = "ineligible"))),
    "insufficient_evidence"
  )
  expect_equal(
    state_of(decide(make_candidates(status = "not_available"))),
    "insufficient_evidence"
  )
})

test_that("a design below its floor sets below_design_floor for every candidate", {
  res <- decide(design = make_design(below_floor = TRUE))
  expect_equal(res$assessments$state, rep("below_design_floor", 2))
  expect_match(res$assessments$reason[1], "declared minimum is 3")
})

test_that("all non-identifiable contrasts set non_identifiable, identity included", {
  res <- decide(
    design = make_design(support = c("non_identifiable", "non_identifiable"))
  )
  expect_equal(res$assessments$state, rep("non_identifiable", 2))
})

test_that("a confounded factor makes corrections non_identifiable but keeps the identity", {
  design <- make_design(
    support = c("reference_supported", "non_identifiable"),
    confounded = "condition"
  )
  res <- decide(
    make_candidates(mad_compression_fraction = 0.5),
    design = design
  )
  expect_equal(state_of(res), "non_identifiable")
  expect_equal(state_of(res, "none"), "eligible_for_review")
  # The guardrail is still recorded under the design gate.
  expect_length(row_of(res)$guardrails_failed[[1]], 1)
})

test_that("compression or joint change above the limit is rejected_distortion", {
  res <- decide(make_candidates(
    mad_compression_fraction = 0.3,
    joint_structure_distance_median = 0.2
  ))
  expect_equal(state_of(res), "rejected_distortion")
  expect_length(row_of(res)$guardrails_failed[[1]], 2)
})

test_that("missing reduction gives insufficient_evidence", {
  res <- decide(make_candidates(emd_relative_reduction = NA_real_))
  expect_equal(state_of(res), "insufficient_evidence")
  expect_match(row_of(res)$reason, "relative reduction")
})

test_that("missing and low coverage give insufficient_evidence with different reasons", {
  na_cov <- decide(make_candidates(coverage = NA_real_))
  low_cov <- decide(make_candidates(coverage = 0.5))
  expect_equal(state_of(na_cov), "insufficient_evidence")
  expect_equal(state_of(low_cov), "insufficient_evidence")
  expect_match(row_of(na_cov)$reason, "unavailable")
  expect_match(row_of(low_cov)$reason, "0.50 is below")
})

test_that("a missing guardrail metric is insufficient_evidence, not a pass", {
  res <- decide(make_candidates(joint_structure_distance_median = NA_real_))
  expect_equal(state_of(res), "insufficient_evidence")
  expect_match(row_of(res)$reason, "joint_structure_distance")
})

test_that("a small reduction is eligible_for_review", {
  res <- decide(make_candidates(emd_relative_reduction = 0.1))
  expect_equal(state_of(res), "eligible_for_review")
  expect_match(row_of(res)$reason, "below the declared minimum")
})

test_that("the reference gate uses the candidate's own evidence comparisons", {
  design <- make_design(
    support = c("reference_supported", "biologically_stratified_conditional")
  )
  cands <- make_candidates()
  cands$evidence_comparison_ids[[2]] <- c("c1", "c2")
  strict <- decide(cands, design = design)
  relaxed <- decide(
    cands,
    design = design,
    policy = make_policy(TRUE, require_ref = FALSE)
  )
  expect_equal(state_of(strict), "eligible_for_review")
  expect_match(row_of(strict)$reason, "measured reference")
  expect_equal(state_of(relaxed), "recommended_under_policy")
})

test_that("a tie names no candidate and states why", {
  cands <- dplyr::bind_rows(make_candidates(), make_candidates()[2, ])
  cands$candidate_id[3] <- "cc2"
  res <- decide(cands)
  expect_equal(sum(res$assessments$state == "recommended_under_policy"), 2)
  expect_true(is.na(res$recommended_candidate_id))
  expect_match(
    res$recommendation_withheld_reason,
    "2 candidates cleared every gate under this policy \\(cc, cc2\\)"
  )
})

test_that("no recommendation says retaining the uncorrected data is valid", {
  res <- decide(policy = make_policy(calibrated = FALSE))
  expect_match(
    res$recommendation_withheld_reason,
    "Retaining the uncorrected data"
  )
})

test_that("the Pareto front drops dominated candidates", {
  cands <- dplyr::bind_rows(make_candidates(), make_candidates()[2, ])
  cands$candidate_id[3] <- "worse"
  cands$emd_relative_reduction[3] <- 0.5
  cands$mad_compression_fraction[3] <- 0.1
  res <- decide(cands)
  expect_true(res$pareto_comparable)
  expect_equal(res$pareto_front, c("cc", "none"))
})

test_that("candidates on different comparisons give no Pareto front", {
  cands <- make_candidates()
  cands$evidence_comparison_ids[[2]] <- "c2"
  res <- decide(
    cands,
    design = make_design(support = rep("reference_supported", 2))
  )
  expect_false(res$pareto_comparable)
  expect_true(any(grepl("comparable evidence", res$warnings)))
})

test_that("policy settings are echoed with their origins", {
  res <- decide(policy = make_policy(calibrated = FALSE))
  expect_equal(res$policy_settings$threshold, thresholds)
  expect_true(all(
    res$policy_settings$origin == "uncalibrated engineering default"
  ))
})

test_that("candidates_from_metrics summarises records and computes reductions", {
  records <- tibble::tibble(
    candidate_id = rep(c("none", "cc"), each = 5),
    metric = rep(
      c(
        "emd_within_population",
        "emd_within_population",
        "emd_within_population",
        "mad_compression_fraction",
        "joint_structure_distance"
      ),
      2
    ),
    marker = rep(c("CD3", "CD4", "CD8", NA, NA), 2),
    comparison_id = "c1",
    value = c(1, 2, 1, 0, 0, 0.5, 1, NA, 0.1, 0.02),
    status = c(rep("ok", 7), "failed", "ok", "ok"),
    weight = 1
  )
  out <- sw_batch_candidates_from_metrics(records, identity_id = "none")
  cc <- out[out$candidate_id == "cc", ]
  # Shared keys: CD3 and CD4; identity sum 3, candidate sum 1.5.
  expect_equal(cc$emd_relative_reduction, 0.5)
  expect_equal(cc$coverage, 2 / 3)
  expect_equal(cc$emd_worst_marker_median, 1)
  expect_equal(out$emd_worst_marker_median[out$candidate_id == "none"], 2)
  expect_equal(cc$mad_compression_fraction, 0.1)
  expect_equal(cc$evidence_comparison_ids[[1]], "c1")

  res <- decide(out)
  expect_equal(state_of(res), "recommended_under_policy")
})

test_that("candidates_from_metrics marks status and rejects duplicate keys", {
  records <- tibble::tibble(
    candidate_id = "none",
    metric = "emd_within_population",
    marker = c("CD3", "CD3"),
    comparison_id = "c1",
    value = 1,
    status = "ok",
    weight = 1
  )
  expect_error(sw_batch_candidates_from_metrics(records, "none"), "duplicate")
  expect_error(sw_batch_candidates_from_metrics(records), "identity_id")
  out <- sw_batch_candidates_from_metrics(
    records[1, ],
    "none",
    candidate_status = c(none = "failed")
  )
  expect_equal(out$status, "failed")
})

# Protected contrast ---------------------------------------------------------

# Planted case - control gap of 1 on CD3 in the identity data, 6 donors per
# condition, 20 events per donor. "keep" shifts every value by 0.3 (the gap is
# kept); "erase" sets every value to one level (the gap is gone).
make_contrast_records <- function(n_donors = 6, with_condition = TRUE) {
  set.seed(1)
  grid <- expand.grid(
    donor = paste0("d", seq_len(n_donors)),
    condition = c("ctrl", "case"),
    event = 1:20,
    stringsAsFactors = FALSE
  )
  grid$donor <- paste(grid$donor, grid$condition, sep = "_")
  level <- ifelse(grid$condition == "case", 3, 2)
  base <- level + stats::rnorm(nrow(grid), sd = 0.05)
  one <- function(id, value) {
    tibble::tibble(
      candidate_id = id,
      donor = grid$donor,
      condition = if (with_condition) grid$condition else NA_character_,
      marker = "CD3",
      population = "T",
      value = value
    )
  }
  dplyr::bind_rows(
    one("none", base),
    one("keep", base + 0.3),
    one("erase", 2.5 + stats::rnorm(nrow(grid), sd = 0.05))
  )
}

contrast_candidates <- function() {
  tibble::tibble(
    candidate_id = c("erase", "keep", "none"),
    status = "ok",
    emd_median = c(0.3, 0.3, 0.8),
    emd_worst_marker_median = c(0.4, 0.4, 1.2),
    emd_relative_reduction = c(0.6, 0.6, 0),
    coverage = 1,
    mad_compression_fraction = c(0.05, 0.05, 0),
    joint_structure_distance_median = c(0.02, 0.02, 0),
    evidence_comparison_ids = list("c1", "c1", "c1")
  )
}

test_that("contrast gate arguments are required and validated", {
  base <- list(
    min_coverage = 0.6,
    min_emd_reduction = 0.2,
    max_mad_compression_fraction = 0.25,
    max_joint_structure_change = 0.1,
    no_effect_emd = 0.05,
    require_reference_supported = TRUE,
    threshold_origins = stats::setNames(rep("uncalibrated", 5), thresholds),
    calibration_artifacts = character(0)
  )
  expect_error(do.call(sw_batch_policy, base), "max_contrast_change")
  expect_error(
    do.call(
      sw_batch_policy,
      c(
        base,
        list(
          max_contrast_change = 0.2,
          min_donors_per_condition = NULL,
          min_contrast_spread_ratio = 2
        )
      )
    ),
    "min_donors_per_condition"
  )
  expect_error(
    do.call(
      sw_batch_policy,
      c(
        base,
        list(
          max_contrast_change = 0.2,
          min_donors_per_condition = 4,
          min_contrast_spread_ratio = 2
        )
      )
    ),
    "one non-empty origin"
  )
  expect_error(
    do.call(
      sw_batch_policy,
      c(
        base,
        list(
          max_contrast_change = NULL,
          min_donors_per_condition = 4,
          min_contrast_spread_ratio = NULL
        )
      )
    ),
    "must be NULL"
  )
  expect_null(make_policy()$max_contrast_change)
  expect_equal(
    make_policy(contrast = 0.2, min_donors = 4)$max_contrast_change,
    0.2
  )
})

test_that("the metric finds a planted erased contrast and a kept one", {
  res <- sw_protected_contrast_change(
    make_contrast_records(),
    identity_id = "none",
    min_donors_per_condition = 4,
    min_contrast_spread_ratio = 2
  )
  worst <- stats::setNames(
    res$summary$worst_abs_change,
    res$summary$candidate_id
  )
  # keep: the gap stays 1, so the change is about 0; erase: the gap is near 0.
  expect_lt(worst[["keep"]], 0.05)
  expect_gt(worst[["erase"]], 0.9)
  expect_equal(res$detail$min_donors[1], 6L)
})

test_that("medians are taken per donor first, not per cell", {
  rec <- make_contrast_records()
  # 2000 extra cells in one ctrl donor of the candidate; pooled cells would
  # move the condition median, a per-donor median barely moves.
  extra <- rec[rec$candidate_id == "keep" & rec$donor == "d1_ctrl", ][1, ]
  extra <- extra[rep(1, 2000), ]
  extra$value <- 10
  res <- sw_protected_contrast_change(
    dplyr::bind_rows(rec, extra),
    identity_id = "none",
    min_donors_per_condition = 4,
    min_contrast_spread_ratio = 2
  )
  expect_lt(
    res$summary$worst_abs_change[res$summary$candidate_id == "keep"],
    0.5
  )
})

test_that("thin conditions are not assessed and the note says why", {
  res <- sw_protected_contrast_change(
    make_contrast_records(n_donors = 3),
    identity_id = "none",
    min_donors_per_condition = 4,
    min_contrast_spread_ratio = 2
  )
  expect_true(all(is.na(res$summary$worst_abs_change)))
  expect_true(all(grepl("not_assessed", res$detail$status)))
  expect_match(res$summary$note[2], "at least 4 donors")
})

test_that("metric errors without required arguments and is deterministic", {
  expect_error(sw_protected_contrast_change(), "no default")
  expect_error(
    sw_protected_contrast_change(
      make_contrast_records(),
      identity_id = "none",
      min_donors_per_condition = 1,
      min_contrast_spread_ratio = 2
    ),
    ">= 2"
  )
  a <- sw_protected_contrast_change(make_contrast_records(), "none", 4, 2)
  b <- sw_protected_contrast_change(make_contrast_records(), "none", 4, 2)
  expect_identical(a, b)
})

test_that("the gate rejects an erased contrast and passes a kept one", {
  res <- decide(
    contrast_candidates(),
    policy = make_policy(calibrated = TRUE, contrast = 0.2, min_donors = 4),
    contrast_records = make_contrast_records()
  )
  expect_equal(state_of(res, "erase"), "rejected_distortion")
  expect_match(row_of(res, "erase")$reason, "protected_contrast_change")
  expect_equal(state_of(res, "keep"), "recommended_under_policy")
  expect_gt(row_of(res, "erase")$protected_contrast_change, 0.9)
  expect_s3_class(res$contrast_detail, "tbl_df")
})

test_that("the gate abstains with a reason when donor or condition is missing", {
  pol <- make_policy(calibrated = TRUE, contrast = 0.2, min_donors = 4)
  for (rec in list(
    NULL,
    make_contrast_records(with_condition = FALSE),
    dplyr::select(make_contrast_records(), -"donor")
  )) {
    res <- decide(contrast_candidates(), policy = pol, contrast_records = rec)
    expect_equal(state_of(res, "keep"), "insufficient_evidence")
    expect_match(row_of(res, "keep")$reason, "protected-contrast gate abstains")
  }
})

test_that("the identity candidate is kept and a gate-off policy ignores records", {
  pol <- make_policy(calibrated = TRUE, contrast = 0.2, min_donors = 4)
  res <- decide(contrast_candidates(), policy = pol, contrast_records = NULL)
  expect_equal(state_of(res, "none"), "eligible_for_review")
  off <- decide(
    contrast_candidates(),
    contrast_records = make_contrast_records()
  )
  expect_equal(state_of(off, "erase"), "recommended_under_policy")
  expect_null(off$contrast_detail)
})

# Two markers: M1 has a real contrast and is kept; M2 has none (noise only).
make_two_marker_records <- function(drop_m2_candidate = FALSE) {
  set.seed(2)
  grid <- expand.grid(
    donor = paste0("d", 1:6),
    condition = c("ctrl", "case"),
    event = 1:10,
    stringsAsFactors = FALSE
  )
  grid$donor <- paste(grid$donor, grid$condition, sep = "_")
  one <- function(id, marker, value) {
    tibble::tibble(
      candidate_id = id,
      donor = grid$donor,
      condition = grid$condition,
      marker = marker,
      population = "T",
      value = value
    )
  }
  m1 <- ifelse(grid$condition == "case", 3, 2) +
    stats::rnorm(nrow(grid), sd = 0.05)
  donor_shift <- stats::rnorm(length(unique(grid$donor)), sd = 0.05)
  m2 <- 1 + donor_shift[match(grid$donor, unique(grid$donor))]
  out <- dplyr::bind_rows(
    one("none", "M1", m1),
    one("none", "M2", m2),
    one("shift", "M1", m1),
    one("shift", "M2", m2 + 0.5 * (grid$condition == "case"))
  )
  if (drop_m2_candidate) {
    out <- out[!(out$candidate_id == "shift" & out$marker == "M2"), ]
  }
  out
}

test_that("a marker with no baseline contrast cannot set the worst change", {
  res <- sw_protected_contrast_change(
    make_two_marker_records(),
    identity_id = "none",
    min_donors_per_condition = 4,
    min_contrast_spread_ratio = 3
  )
  det <- res$detail
  expect_match(det$status[det$marker == "M2"], "below donor spread")
  expect_equal(det$status[det$marker == "M1"], "assessed")
  expect_lt(
    res$summary$worst_abs_change[res$summary$candidate_id == "shift"],
    0.05
  )
  expect_equal(
    res$summary$n_not_assessed[res$summary$candidate_id == "shift"],
    1L
  )
  expect_match(res$summary$note[2], "1 not assessed")
  # With the ratio at zero the noise marker is judged and drives the worst value.
  open_gate <- sw_protected_contrast_change(
    make_two_marker_records(),
    identity_id = "none",
    min_donors_per_condition = 4,
    min_contrast_spread_ratio = 0
  )
  expect_gt(open_gate$summary$worst_abs_change[2], 1)
})

test_that("the decision records rows the gate could not judge", {
  pol <- make_policy(calibrated = TRUE, contrast = 0.2, min_donors = 4)
  cand <- contrast_candidates()[0, ]
  cand <- tibble::tibble(
    candidate_id = c("shift", "none"),
    status = "ok",
    emd_median = c(0.3, 0.8),
    emd_worst_marker_median = c(0.4, 1.2),
    emd_relative_reduction = c(0.6, 0),
    coverage = 1,
    mad_compression_fraction = c(0.05, 0),
    joint_structure_distance_median = c(0.02, 0),
    evidence_comparison_ids = list("c1", "c1")
  )
  res <- decide(
    cand,
    policy = pol,
    contrast_records = make_two_marker_records()
  )
  row <- row_of(res, "shift")
  expect_equal(row$protected_contrast_n_not_assessed, 1L)
  expect_match(row$protected_contrast_note, "1 not assessed")
  expect_match(row$reason, "1 protected-contrast row\\(s\\) were not assessed")
})

test_that("a candidate absent from identity-only records gets the missing-record reason", {
  pol <- make_policy(calibrated = TRUE, contrast = 0.2, min_donors = 4)
  only_identity <- make_contrast_records()
  only_identity <- only_identity[only_identity$candidate_id == "none", ]
  res <- decide(
    contrast_candidates(),
    policy = pol,
    contrast_records = only_identity
  )
  expect_match(row_of(res, "keep")$reason, "no donor and condition records")
  expect_false(grepl("is the baseline", row_of(res, "keep")$reason))
  null_res <- decide(
    contrast_candidates(),
    policy = pol,
    contrast_records = NULL
  )
  expect_match(row_of(null_res, "keep")$reason, "no donor or condition data")
})

test_that("the spread ratio is required and checked", {
  expect_error(
    sw_protected_contrast_change(make_contrast_records(), "none", 4),
    "min_contrast_spread_ratio"
  )
  expect_error(
    sw_protected_contrast_change(make_contrast_records(), "none", 4, -1),
    ">= 0"
  )
})
