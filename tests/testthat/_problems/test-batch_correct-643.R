# Extracted from test-batch_correct.R:643

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SpectraWeaveR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
make_gate_audit <- function(batches) {
  samples <- data.frame(
    sample_id = paste0("s", seq_along(batches)),
    batch = batches,
    condition = rep(c("ctrl", "case"), length.out = length(batches)),
    reference_role = "none",
    reference_id = NA_character_
  )
  sw_design_audit(
    samples,
    minimum_batches = 1,
    minimum_batches_origin = "test fixture",
    require_shared_reference = FALSE,
    protected = "condition"
  )
}
gate_df <- tibble::tibble(
  CD3 = c(1, 2, 3, 4),
  batch = c("B1", "B1", "B2", "B2"),
  sample = c("S1", "S1", "S2", "S2")
)

# test -------------------------------------------------------------------------
skip_if_not_installed("cyCombine")
audit <- make_gate_audit(c("b1", "b1", "b1", "b1"))
row <- audit$eligibility[audit$eligibility$name == "cycombine", ]
expect_false(row$eligible)
expect_error(
    sw_correct_run(gate_df, markers = "CD3", design_audit = audit),
    "Design audit row 'cycombine' \\(sw_correct_run\\) is ineligible"
  )
