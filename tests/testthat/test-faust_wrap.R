# tests/testthat/test-faust_wrap.R
# faust is replaced by a fake through .sw_faust_fn(), so the tests need only
# flowCore and flowWorkspace. The fake writes what the wrapper reads from disk:
# faustData/faustCountMatrix.rds and sampleData/<sample>/faustAnnotation.csv.

skip_if_not_installed("flowCore")
skip_if_not_installed("flowWorkspace")

faust_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  mk <- function() {
    m <- matrix(
      rnorm(n * 3),
      n,
      3,
      dimnames = list(NULL, c("ch1", "ch2", "ch3"))
    )
    m
  }
  list(f1 = mk(), f2 = mk(), f3 = mk(), f4 = mk())
}

faust_table <- function() {
  data.frame(
    sample = c("f1", "f2", "f3", "f4"),
    batch = c("B1", "B1", "B2", "B2"),
    donor = c("d1", "d2", "d3", "d4"),
    role = c("ref", "donor", "ref", "donor")
  )
}

faust_markers <- c(ch1 = "M1", ch2 = "M2")

# Fake faust(): labels repeat a fixed pattern. `labels` is the pattern, `args`
# collects what the wrapper passed.
fake_faust <- function(labels, seen = new.env()) {
  function(name, pkg = "faust") {
    stopifnot(name == "faust")
    function(gatingSet, projectPath, annotationsApproved, ...) {
      seen$args <- c(list(projectPath = projectPath), list(...))
      kept <- attr(gatingSet, "kept_rows")
      root <- file.path(projectPath, "faustData")
      counts <- list()
      for (s in names(kept)) {
        lab <- rep(labels, length.out = length(kept[[s]]))
        dir.create(file.path(root, "sampleData", s), recursive = TRUE)
        utils::write.table(
          data.frame(lab),
          file.path(root, "sampleData", s, "faustAnnotation.csv"),
          sep = ",",
          row.names = FALSE,
          col.names = FALSE,
          quote = FALSE
        )
        counts[[s]] <- table(lab)
      }
      if (annotationsApproved) {
        all_labs <- sort(unique(unlist(lapply(counts, names))))
        cm <- t(vapply(
          counts,
          \(t) as.numeric(t[all_labs]),
          numeric(length(all_labs))
        ))
        cm[is.na(cm)] <- 0
        dimnames(cm) <- list(names(kept), all_labs)
        saveRDS(cm, file.path(root, "faustCountMatrix.rds"))
      }
      NULL
    }
  }
}

run_args <- function(gs, path, ...) {
  args <- list(
    gs = gs,
    project_path = path,
    starting_pop = "root",
    markers = c("M1", "M2"),
    depth_score_threshold = 0.01,
    selection_quantile = 0.5,
    experimental_unit = "",
    imputation_hierarchy = "",
    channel_bounds = "",
    name_occurrence_num = 0,
    supervised_list = NA,
    seed = 1,
    threads = 1,
    density_subsample_threshold = 1e6,
    density_subsample_size = 1e6,
    density_subsample_iterations = 1,
    annotation_forest_depth = 3,
    annotations_approved = TRUE,
    draw_histograms = FALSE
  )
  utils::modifyList(args, list(...))
}

make_gs <- function(data = faust_data(), n_events = 50) {
  sw_faust_gatingset(
    data,
    faust_markers,
    faust_table(),
    n_events = n_events,
    subsample_seed = 1
  )
}

make_run <- function(gs, labels, tag) {
  path <- file.path(withr::local_tempdir(.local_envir = parent.frame()), tag)
  testthat::local_mocked_bindings(
    .sw_faust_fn = fake_faust(labels),
    .env = parent.frame()
  )
  do.call(sw_faust_run, run_args(gs, path))
}

test_that("sw_faust_gatingset stops on a missing argument", {
  expect_error(
    sw_faust_gatingset(
      faust_data(),
      faust_markers,
      faust_table(),
      n_events = 5
    ),
    "subsample_seed"
  )
})

test_that("sw_faust_gatingset stops when files and sample table differ", {
  tbl <- faust_table()
  tbl$sample[2] <- "other"
  expect_error(
    sw_faust_gatingset(
      faust_data(),
      faust_markers,
      tbl,
      n_events = 5,
      subsample_seed = 1
    ),
    "names\\(data\\) equals sample_table\\$sample"
  )
})

test_that("sw_faust_gatingset keeps markers by name and rows in order", {
  gs <- make_gs()
  expect_s4_class(gs, "GatingSet")
  expect_equal(unname(flowWorkspace::markernames(gs)), c("M1", "M2"))
  kept <- attr(gs, "kept_rows")
  expect_length(kept[["f1"]], 50)
  expect_false(is.unsorted(kept[["f1"]]))
  # A file at or under n_events keeps every row.
  expect_equal(attr(make_gs(n_events = 60), "kept_rows")[["f1"]], 1:60)
  # Same data shape and seed give the same rows.
  expect_equal(attr(make_gs(faust_data(seed = 2)), "kept_rows"), kept)
})

test_that("sw_faust_run stops on a missing argument", {
  gs <- make_gs()
  args <- run_args(gs, withr::local_tempdir())
  args$seed <- NULL
  expect_error(do.call(sw_faust_run, args), "'seed'")
})

test_that("sw_faust_run stops when faustData exists", {
  gs <- make_gs()
  path <- withr::local_tempdir()
  dir.create(file.path(path, "faustData"))
  expect_error(do.call(sw_faust_run, run_args(gs, path)), "faustData")
})

test_that("sw_faust_run stops when markers do not match desc", {
  gs <- make_gs()
  expect_error(
    do.call(
      sw_faust_run,
      run_args(
        gs,
        file.path(withr::local_tempdir(), "p"),
        markers = c("M1", "CD3")
      )
    ),
    "CD3"
  )
})

test_that("sw_faust_run passes the arguments through and echoes them", {
  gs <- make_gs()
  seen <- new.env()
  local_mocked_bindings(
    .sw_faust_fn = fake_faust(c("M1~1~2~", "M1~2~2~"), seen)
  )
  run <- do.call(
    sw_faust_run,
    run_args(gs, file.path(withr::local_tempdir(), "p"), seed = 7, threads = 3)
  )
  expect_s3_class(run, "sw_faust_run")
  expect_equal(attr(run, "status"), "complete")
  expect_equal(seen$args$seedValue, 7)
  expect_equal(seen$args$threadNum, 3)
  expect_equal(seen$args$plottingDevice, "png")
  expect_equal(seen$args$archDescriptionList, list(targetArch = "singleCPU"))
  expect_equal(attr(run, "params")$seed, 7)
})

test_that("a thresholds_only run stops in sw_faust_labels", {
  gs <- make_gs()
  local_mocked_bindings(.sw_faust_fn = fake_faust(c("M1~1~2~", "M1~2~2~")))
  run <- do.call(
    sw_faust_run,
    run_args(
      gs,
      file.path(withr::local_tempdir(), "p"),
      annotations_approved = FALSE
    )
  )
  expect_equal(attr(run, "status"), "thresholds_only")
  expect_error(sw_faust_labels(run, faust_data()), "scoreLines")
})

test_that("sw_faust_labels returns cells and proportions that add up", {
  gs <- make_gs()
  run <- make_run(gs, c("M1~1~2~", "M1~2~2~", "0_0"), "a")
  lab <- sw_faust_labels(run, faust_data())
  expect_s3_class(lab, "sw_faust_labels")
  expect_named(lab$cells, c("sample", "cell", "phenotype"))
  expect_named(lab$proportions, c("sample", "phenotype", "n", "proportion"))
  expect_equal(nrow(lab$cells), 4 * 50)
  expect_equal(
    unname(tapply(lab$proportions$proportion, lab$proportions$sample, sum)),
    rep(1, 4),
    ignore_attr = TRUE
  )
  expect_true("0_0" %in% lab$proportions$phenotype)
})

test_that("sw_faust_labels stops when label counts disagree with the count matrix", {
  gs <- make_gs()
  run <- make_run(gs, c("M1~1~2~", "M1~2~2~"), "a")
  cm <- readRDS(file.path(unclass(run), "faustData", "faustCountMatrix.rds"))
  cm[1, 1] <- cm[1, 1] + 1
  saveRDS(cm, file.path(unclass(run), "faustData", "faustCountMatrix.rds"))
  expect_error(sw_faust_labels(run, faust_data()), "label counts equal")
})

compare_args <- function(lab_b, lab_a, ...) {
  data <- faust_data()
  tbl <- faust_table()
  args <- list(
    before = lab_b,
    after = lab_a,
    data_before = data,
    data_after = data,
    batch = stats::setNames(tbl$batch, tbl$sample),
    role = stats::setNames(tbl$role, tbl$sample),
    markers = c("ch1", "ch2"),
    match_by = "exact",
    ari_scope = "file",
    median_scope = "file"
  )
  utils::modifyList(args, list(...))
}

test_that("identical labels give ARI 1 and a long table in the diagnose shape", {
  gs <- make_gs()
  lab_b <- sw_faust_labels(
    make_run(gs, c("M1~1~2~", "M1~2~2~", "0_0"), "b"),
    faust_data()
  )
  lab_a <- sw_faust_labels(
    make_run(gs, c("M1~1~2~", "M1~2~2~", "0_0"), "a"),
    faust_data()
  )
  out <- do.call(sw_faust_compare, compare_args(lab_b, lab_a))
  expect_s3_class(out, "sw_faust_compare")
  expect_named(
    out,
    c("module", "subject", "check", "item", "phenotype", "value", "text")
  )
  expect_true(all(out$module == "faust"))
  expect_equal(out$value[out$check == "ari"], rep(1, 4), ignore_attr = TRUE)
  expect_setequal(
    unique(out$check),
    c(
      "ari",
      "prop_before",
      "prop_after",
      "prop_diff",
      "median_before",
      "median_after",
      "median_diff",
      "prop_by_batch_before",
      "prop_by_batch_after"
    )
  )
  expect_true(all(out$value[out$check %in% c("prop_diff", "median_diff")] == 0))
  # One row per file, phenotype and marker for the medians.
  expect_equal(sum(out$check == "median_before"), 4 * 3 * 2)
  # Roles are kept apart in the batch spread.
  expect_setequal(
    out$subject[out$check == "prop_by_batch_after"],
    c("ref", "donor")
  )
  expect_false(any(c("pass", "fail", "flag") %in% names(out)))
})

test_that("a phenotype on one side only gets NA on the other and the zero bin stays", {
  gs <- make_gs()
  lab_b <- sw_faust_labels(
    make_run(gs, c("M1~1~2~", "M1~2~2~", "0_0_0"), "b"),
    faust_data()
  )
  lab_a <- sw_faust_labels(
    make_run(gs, c("M1~1~2~", "M1~2~2~", "M2~2~2~"), "a"),
    faust_data()
  )
  out <- do.call(sw_faust_compare, compare_args(lab_b, lab_a))
  only_b <- out[out$check == "prop_after" & out$phenotype == "0_0_0", ]
  only_a <- out[out$check == "prop_before" & out$phenotype == "M2~2~2~", ]
  expect_equal(nrow(only_b), 4)
  expect_true(all(is.na(only_b$value)))
  expect_true(all(is.na(only_a$value)))
  expect_true(all(
    out$value[out$check == "prop_before" & out$phenotype == "0_0_0"] > 0
  ))
  expect_true(all(is.na(out$value[
    out$check == "prop_diff" & out$phenotype == "0_0_0"
  ])))
  expect_true(all(.sw_faust_is_unclassified(c("0_0_0", "0_0_0_0_0", "0"))))
  expect_false(.sw_faust_is_unclassified("M1~1~2~"))
})

test_that("changed labels lower the ARI and pooled scope gives one row", {
  gs <- make_gs()
  lab_b <- sw_faust_labels(
    make_run(gs, c("M1~1~2~", "M1~2~2~"), "b"),
    faust_data()
  )
  lab_a <- sw_faust_labels(
    make_run(gs, c("M1~1~2~", "M1~1~2~", "M1~2~2~"), "a"),
    faust_data()
  )
  out <- do.call(
    sw_faust_compare,
    compare_args(lab_b, lab_a, ari_scope = "pooled", median_scope = "pooled")
  )
  ari <- out$value[out$check == "ari"]
  expect_length(ari, 1)
  expect_lt(ari, 1)
  expect_equal(unique(out$subject[out$check == "median_before"]), "pooled")
})

test_that("sw_faust_compare stops when the kept rows differ", {
  lab_b <- sw_faust_labels(
    make_run(make_gs(n_events = 50), c("M1~1~2~", "M1~2~2~"), "b"),
    faust_data()
  )
  lab_a <- sw_faust_labels(
    make_run(make_gs(n_events = 40), c("M1~1~2~", "M1~2~2~"), "a"),
    faust_data()
  )
  expect_error(
    do.call(sw_faust_compare, compare_args(lab_b, lab_a)),
    "same kept rows"
  )
})

test_that("sw_faust_compare stops on a missing argument", {
  gs <- make_gs()
  lab <- sw_faust_labels(
    make_run(gs, c("M1~1~2~", "M1~2~2~"), "b"),
    faust_data()
  )
  args <- compare_args(lab, lab)
  args$median_scope <- NULL
  expect_error(do.call(sw_faust_compare, args), "'median_scope'")
  expect_error(
    do.call(sw_faust_compare, compare_args(lab, lab, match_by = "fuzzy")),
    "exact"
  )
})
