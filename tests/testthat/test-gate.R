# tests/testthat/test-gate.R
# Unit tests for R/gate.R — Automated gating (openCyto)

tmpl_args <- list(debris_min = 50000, singlet_wider_gate = TRUE,
                  flowclust_k = 2, flowclust_target = 1, pool_samples = TRUE)
make_template <- function(file, ...) {
  args <- utils::modifyList(tmpl_args, list(...))
  do.call(sw_gate_template, c(list(file), args))
}

test_that("sw_gate_template creates valid CSV", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  result <- make_template(tmp_file)

  expect_true(file.exists(tmp_file))
  expect_equal(result, tmp_file)

  # Read and validate structure
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 3)

  # Required columns for openCyto gating template
  required_cols <- c("alias", "pop", "parent", "dims", "gating_method")
  expect_true(all(required_cols %in% names(gt)))

  # Check aliases
  expect_true("nonDebris" %in% gt$alias)
  expect_true("singlets" %in% gt$alias)
  expect_true("lymphocytes" %in% gt$alias)
})

test_that("sw_gate_template creates parent directory", {
  tmp_dir <- file.path(tempdir(), "nested", "gating_test")
  tmp_file <- file.path(tmp_dir, "template.csv")
  on.exit(unlink(file.path(tempdir(), "nested"), recursive = TRUE))

  make_template(tmp_file)
  expect_true(file.exists(tmp_file))
})

test_that("sw_gate_template rejects invalid template_type", {
  expect_error(sw_gate_template(tempfile(), template_type = "invalid", debris_min = 1,
                                   singlet_wider_gate = TRUE, pool_samples = TRUE),
               "Unknown template_type")
})

test_that("sw_gate_template rejects non-character output_file", {
  expect_error(sw_gate_template(42),
               "single file path")
})

# =========================================================================
# New template types
# =========================================================================

test_that("sw_gate_template creates myeloid template", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  make_template(tmp_file, template_type = "myeloid")
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 4)
  expect_true("nonDebris" %in% gt$alias)
  expect_true("singlets" %in% gt$alias)
  expect_true("CD3neg" %in% gt$alias)
  expect_true("HLADRpos" %in% gt$alias)
  # Verify parent hierarchy: every parent must be root or an alias
  valid_parents <- c("root", gt$alias)
  expect_true(all(gt$parent %in% valid_parents))
})

test_that("sw_gate_template creates nk template", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  make_template(tmp_file, template_type = "nk")
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 4)
  expect_true("CD56pos" %in% gt$alias)
  valid_parents <- c("root", gt$alias)
  expect_true(all(gt$parent %in% valid_parents))
})

test_that("sw_gate_template creates treg template", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  make_template(tmp_file, template_type = "treg")
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 5)
  expect_true("CD3pos" %in% gt$alias)
  expect_true("CD4pos" %in% gt$alias)
  expect_true("Treg" %in% gt$alias)
  valid_parents <- c("root", gt$alias)
  expect_true(all(gt$parent %in% valid_parents))
})

test_that("sw_gate_template creates full_pbmc template", {
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))

  make_template(tmp_file, template_type = "full_pbmc")
  gt <- utils::read.csv(tmp_file, stringsAsFactors = FALSE)
  expect_true(nrow(gt) >= 8)
  # Should have branching: CD3pos and CD3neg both under singlets
  expect_true("CD3pos" %in% gt$alias)
  expect_true("CD3neg" %in% gt$alias)
  expect_true("CD19pos" %in% gt$alias)
  expect_true("CD14pos" %in% gt$alias)
  valid_parents <- c("root", gt$alias)
  expect_true(all(gt$parent %in% valid_parents))
})

test_that("sw_gate_run rejects empty fcs_files", {
  expect_error(sw_gate_run(character(0), "template.csv"),
               "non-empty character vector")
})

test_that("sw_gate_run rejects missing template file", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")

  # Create a dummy FCS file path (won't actually be read due to template check)
  expect_error(sw_gate_run("dummy.fcs", "/nonexistent/template.csv"),
               "not found")
})

test_that("sw_gate_extract rejects non-GatingSet input", {
  skip_if_not_installed("flowWorkspace")

  expect_error(sw_gate_extract("not_a_gs", "/singlets"),
               "GatingSet")
})

test_that("sw_gate_extract rejects non-character node", {
  skip_if_not_installed("flowWorkspace")

  expect_error(sw_gate_extract("not_a_gs", 42),
               "single character string")
})

# =========================================================================
# Integration tests — actual gating execution
# =========================================================================

test_that("sw_gate_run runs end-to-end with synthetic FCS", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("openCyto")

  set.seed(42)
  n <- 500
  # Bimodal FSC-A: debris (low) + cells (high)
  fsc_a <- c(abs(rnorm(100, 20000, 5000)),
             abs(rnorm(400, 150000, 30000)))
  fsc_h <- fsc_a * 0.7 + rnorm(n, sd = 5000)
  ssc_a <- abs(rnorm(n, 50000, 15000))

  mat <- cbind("FSC-A" = fsc_a, "FSC-H" = fsc_h, "SSC-A" = ssc_a)
  storage.mode(mat) <- "double"
  ff <- flowCore::flowFrame(mat)

  # Write to temp FCS
  tmp_fcs <- tempfile(fileext = ".fcs")
  flowCore::write.FCS(ff, tmp_fcs)
  on.exit(unlink(tmp_fcs))

  # Build template
  tmp_template <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_template), add = TRUE)
  make_template(tmp_template)

  gs <- sw_gate_run(tmp_fcs, tmp_template, transform_channels = NULL)
  expect_true(methods::is(gs, "GatingSet"))

  # Verify gate nodes exist
  nodes <- flowWorkspace::gs_get_pop_paths(gs)
  expect_true(length(nodes) >= 2)
})

test_that("sw_gate_extract returns named list of flowFrames", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("openCyto")

  set.seed(42)
  n <- 500
  fsc_a <- c(abs(rnorm(100, 20000, 5000)),
             abs(rnorm(400, 150000, 30000)))
  fsc_h <- fsc_a * 0.7 + rnorm(n, sd = 5000)
  ssc_a <- abs(rnorm(n, 50000, 15000))

  mat <- cbind("FSC-A" = fsc_a, "FSC-H" = fsc_h, "SSC-A" = ssc_a)
  storage.mode(mat) <- "double"
  ff <- flowCore::flowFrame(mat)

  tmp_fcs <- tempfile(fileext = ".fcs")
  flowCore::write.FCS(ff, tmp_fcs)
  on.exit(unlink(tmp_fcs))

  tmp_template <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_template), add = TRUE)
  make_template(tmp_template)

  gs <- sw_gate_run(tmp_fcs, tmp_template, transform_channels = NULL)
  nodes <- flowWorkspace::gs_get_pop_paths(gs)

  # Extract from first non-root node
  leaf <- nodes[length(nodes)]
  ff_list <- sw_gate_extract(gs, leaf)
  expect_true(is.list(ff_list))
  expect_true(length(ff_list) >= 1)
  expect_true(methods::is(ff_list[[1]], "flowFrame"))
})

test_that("sw_gate_template requires every cutoff argument", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  expect_error(sw_gate_template(tmp), "required")
  expect_error(sw_gate_template(tmp, debris_min = 1, singlet_wider_gate = TRUE),
               "pool_samples")
  expect_error(
    sw_gate_template(tmp, debris_min = 1, singlet_wider_gate = TRUE,
                     pool_samples = TRUE),
    "flowclust_k"
  )
  # templates without flowClust do not need K or target
  expect_no_error(
    sw_gate_template(tmp, "nk", debris_min = 1, singlet_wider_gate = TRUE,
                     pool_samples = TRUE)
  )
})

test_that("sw_gate_template writes the supplied values into the CSV", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  make_template(tmp, debris_min = 12345, singlet_wider_gate = FALSE,
                flowclust_k = 3, flowclust_target = 2, pool_samples = FALSE)
  gt <- utils::read.csv(tmp, stringsAsFactors = FALSE)
  expect_equal(gt$gating_args[gt$alias == "nonDebris"], "min = 12345")
  expect_equal(gt$gating_args[gt$alias == "singlets"], "wider_gate = FALSE")
  expect_equal(gt$gating_args[gt$alias == "lymphocytes"],
               "K = 3, target = 2")
  expect_false(any(gt$collapseDataForGating))

  make_template(tmp, pool_samples = TRUE)
  gt <- utils::read.csv(tmp, stringsAsFactors = FALSE)
  expect_true(all(gt$collapseDataForGating))
})

test_that("sw_gate_run needs cofactor when it transforms with arcsinh", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("openCyto")
  tmpl <- tempfile(fileext = ".csv")
  on.exit(unlink(tmpl))
  make_template(tmpl)
  expect_error(
    sw_gate_run("x.fcs", tmpl, transform_channels = "FSC-A"),
    "'cofactor' is required"
  )
  expect_false("..." %in% names(formals(sw_gate_run)))
})

# =========================================================================
# Pooling: sw_gate_run sets groupBy when collapseDataForGating is TRUE
# =========================================================================

# Three samples whose FSC-A peaks sit at different places, so the valley of
# each sample differs. Debris (low FSC-A) and doublets (low FSC-H) give the
# nonDebris and singlets gates something to cut.
pooling_fs <- function() {
  one <- function(seed, scale) {
    set.seed(seed)
    n_cells <- 3000
    n_debris <- 200
    n_doublet <- 200
    a <- abs(c(
      rnorm(n_cells, 6e5 * scale, 3e4),
      rnorm(n_debris, 5e4 * scale, 1.5e4),
      rnorm(n_doublet, 6e5 * scale, 3e4)
    ))
    h <- a * 0.8 + rnorm(length(a), 0, 8000)
    h[length(a) - n_doublet + seq_len(n_doublet)] <-
      a[length(a) - n_doublet + seq_len(n_doublet)] * 0.45
    flowCore::flowFrame(cbind(
      "FSC-A" = a,
      "FSC-H" = h,
      "SSC-A" = abs(rnorm(length(a), 1e5, 2e4))
    ))
  }
  fs <- flowCore::flowSet(list(one(1, 0.9), one(2, 1), one(3, 1.1)))
  flowCore::sampleNames(fs) <- c("a", "b", "c")
  fs
}

# Template with the first two rows of sw_gate_template() (nonDebris,
# singlets), after edit() has changed the data frame.
pooling_template <- function(pool_samples, edit = identity) {
  file <- tempfile(fileext = ".csv")
  suppressMessages(sw_gate_template(
    file,
    template_type = "myeloid",
    debris_min = 1e4,
    singlet_wider_gate = TRUE,
    pool_samples = pool_samples
  ))
  tmpl <- utils::read.csv(file, stringsAsFactors = FALSE)
  tmpl <- edit(tmpl[tmpl$alias %in% c("nonDebris", "singlets"), ])
  utils::write.csv(tmpl, file, row.names = FALSE)
  file
}

# One entry per sample: the nonDebris cut and the singlets polygon.
pooling_gates <- function(file) {
  gs <- suppressMessages(sw_gate_run(pooling_fs(), file))
  sn <- flowWorkspace::sampleNames(gs)
  gate_of <- function(node) {
    lapply(sn, function(s) {
      g <- flowWorkspace::gh_pop_get_gate(gs[[s]], node)
      if (methods::is(g, "rectangleGate")) unname(g@min) else g@boundaries
    })
  }
  list(
    cut = unlist(gate_of("nonDebris")),
    singlets = gate_of("singlets")
  )
}

skip_if_no_gate_stack <- function() {
  for (p in c("flowCore", "flowWorkspace", "openCyto")) {
    skip_if_not_installed(p)
  }
}

test_that("sw_gate_run pools when pool_samples is TRUE in the template", {
  skip_if_no_gate_stack()
  file <- pooling_template(TRUE)
  on.exit(unlink(file))
  # The template from sw_gate_template() leaves groupBy empty.
  expect_true(all(is.na(utils::read.csv(file)$groupBy)))
  gates <- pooling_gates(file)
  expect_length(unique(gates$cut), 1)
  expect_length(unique(gates$singlets), 1)
})

test_that("sw_gate_run fits one gate per sample when pool_samples is FALSE", {
  skip_if_no_gate_stack()
  file <- pooling_template(FALSE)
  on.exit(unlink(file))
  gates <- pooling_gates(file)
  expect_length(unique(gates$cut), 3)
  expect_length(unique(gates$singlets), 3)
})

test_that("sw_gate_run leaves an explicit groupBy alone", {
  skip_if_no_gate_stack()
  # groupBy = 1 puts each sample in its own group, so the cuts differ although
  # collapseDataForGating is TRUE. Setting it to 3 would pool them.
  file <- pooling_template(TRUE, function(tmpl) {
    tmpl$groupBy <- 1L
    tmpl
  })
  on.exit(unlink(file))
  expect_identical(.sw_gate_pool_template(file, 3L), file)
  gates <- pooling_gates(file)
  expect_length(unique(gates$cut), 3)
})

test_that(".sw_gate_pool_template changes only rows that pool with an empty groupBy", {
  tmpl <- data.frame(
    alias = c("one", "two", "three", "four"),
    collapseDataForGating = c(TRUE, FALSE, TRUE, TRUE),
    groupBy = c(NA, NA, "2", ""),
    stringsAsFactors = FALSE
  )
  file <- tempfile(fileext = ".csv")
  on.exit(unlink(c(file, out)))
  utils::write.csv(tmpl, file, row.names = FALSE)
  out <- suppressMessages(.sw_gate_pool_template(file, 5L))
  expect_false(identical(out, file))
  got <- utils::read.csv(out, stringsAsFactors = FALSE, colClasses = "character")
  expect_equal(got$groupBy, c("5", NA, "2", "5"))
  expect_equal(got$collapseDataForGating, c("TRUE", "FALSE", "TRUE", "TRUE"))
  # The input file is not edited.
  kept <- utils::read.csv(file, stringsAsFactors = FALSE, colClasses = "character")
  expect_equal(kept$groupBy, c(NA, NA, "2", ""))
})
