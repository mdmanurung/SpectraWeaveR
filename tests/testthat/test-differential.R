# tests/testthat/test-differential.R
# Unit tests for R/differential.R — Differential Abundance & Expression

# ---------------------------------------------------------------------------
# Helper: synthetic cell-level tibble
# ---------------------------------------------------------------------------
.make_corrected <- function(n_per_group = 50, seed = 42) {
  set.seed(seed)
  # 2 groups × 3 samples × 2 clusters
  groups  <- rep(c("ctrl", "trt"), each = 3L)
  samples <- paste0("S", seq_len(6L))
  clusters <- c(1L, 2L)

  rows <- do.call(rbind, lapply(seq_along(samples), function(si) {
    sn  <- samples[si]
    grp <- groups[si]
    do.call(rbind, lapply(clusters, function(cl) {
      n    <- n_per_group
      # cluster 2 has elevated CD4 in trt group
      cd4_mu <- if (cl == 2L && grp == "trt") 2.5 else 0.5
      data.frame(
        sample    = sn,
        batch     = if (si <= 3L) "B1" else "B2",
        condition = grp,
        cluster   = cl,
        CD3       = stats::rnorm(n, mean = 1.5, sd = 0.3),
        CD4       = stats::rnorm(n, mean = cd4_mu, sd = 0.3),
        CD8a      = stats::rnorm(n, mean = 0.4, sd = 0.2),
        stringsAsFactors = FALSE
      )
    }))
  }))
  tibble::as_tibble(rows)
}

.make_meta <- function() {
  data.frame(
    sample    = paste0("S", seq_len(6L)),
    condition = rep(c("ctrl", "trt"), each = 3L),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# .build_count_matrix (tested indirectly through sw_differential_abundance)
# ---------------------------------------------------------------------------

test_that(".compute_cluster_medians computes correct medians", {
  df <- data.frame(
    sample  = c("S1", "S1", "S2", "S2"),
    cluster = c(1L,    1L,   2L,    2L),
    CD3     = c(1,     3,    2,     4),
    stringsAsFactors = FALSE
  )
  # S1 / cluster 1: median of c(1, 3) = 2
  # S2 / cluster 2: median of c(2, 4) = 3
  med <- SpectraWeaveR:::.compute_cluster_medians(
    df, cluster_col = "cluster", markers = "CD3",
    sample_col = "sample", min_cells = 1L
  )
  expect_equal(med[["1"]]["S1", "CD3"], 2)
  expect_equal(med[["2"]]["S2", "CD3"], 3)
})

test_that(".compute_cluster_medians returns NA for samples below min_cells", {
  df <- data.frame(
    sample  = c("S1", "S2"),
    cluster = c(1L,    1L),
    CD3     = c(5,     5),
    stringsAsFactors = FALSE
  )
  # min_cells = 2, each sample has only 1 event → should be NA
  med <- SpectraWeaveR:::.compute_cluster_medians(
    df, cluster_col = "cluster", markers = "CD3",
    sample_col = "sample", min_cells = 2L
  )
  expect_true(all(is.na(med[["1"]])))
})

# ---------------------------------------------------------------------------
# sw_differential_abundance — validation
# ---------------------------------------------------------------------------

test_that("sw_differential_abundance rejects non-data.frame x", {
  expect_error(
    sw_differential_abundance("bad", .make_meta(), group_col = "condition"),
    "data.frame or tibble"
  )
})

test_that("sw_differential_abundance rejects x without cluster column", {
  x <- tibble::tibble(sample = "S1", CD3 = 1)
  expect_error(
    sw_differential_abundance(x, .make_meta(), group_col = "condition"),
    "'cluster' column"
  )
})

test_that("sw_differential_abundance rejects x without sample column", {
  x <- tibble::tibble(cluster = 1L, CD3 = 1)
  expect_error(
    sw_differential_abundance(x, .make_meta(), group_col = "condition"),
    "'sample' column"
  )
})

test_that("sw_differential_abundance rejects meta missing group_col", {
  x    <- .make_corrected()
  meta <- data.frame(sample = paste0("S", 1:6), stringsAsFactors = FALSE)
  expect_error(
    sw_differential_abundance(x, meta, group_col = "condition"),
    "missing required column"
  )
})

test_that("sw_differential_abundance rejects invalid min_prop", {
  x <- .make_corrected()
  expect_error(
    sw_differential_abundance(x, .make_meta(), group_col = "condition",
                               min_prop = 1.5),
    "min_prop"
  )
})

test_that("sw_differential_abundance rejects invalid method", {
  x <- .make_corrected()
  expect_error(
    sw_differential_abundance(x, .make_meta(), group_col = "condition",
                               method = "bad"),
    "should be one of"
  )
})

# ---------------------------------------------------------------------------
# sw_differential_abundance — integration (requires edgeR)
# ---------------------------------------------------------------------------

test_that("sw_differential_abundance edgeR returns correct output structure", {
  skip_if_not_installed("edgeR")

  x    <- .make_corrected()
  meta <- .make_meta()
  result <- sw_differential_abundance(x, meta, group_col = "condition",
                                       method = "edgeR", min_prop = 0)

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("cluster", "logFC", "pvalue", "adj_pvalue",
                     "mean_prop_A", "mean_prop_B") %in% names(result)))
  expect_equal(nrow(result), 2L)  # 2 clusters
  expect_true(all(result$adj_pvalue >= 0 & result$adj_pvalue <= 1))
  expect_equal(attr(result, "groups"), c("ctrl", "trt"))
})

test_that("sw_differential_abundance voom returns correct output structure", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")

  x    <- .make_corrected()
  meta <- .make_meta()
  result <- sw_differential_abundance(x, meta, group_col = "condition",
                                       method = "voom", min_prop = 0)

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("cluster", "logFC", "pvalue", "adj_pvalue") %in%
                    names(result)))
  expect_equal(nrow(result), 2L)
})

test_that("sw_differential_abundance min_prop filters clusters", {
  skip_if_not_installed("edgeR")

  x    <- .make_corrected()
  meta <- .make_meta()
  # Set min_prop so high that all clusters are excluded
  expect_error(
    sw_differential_abundance(x, meta, group_col = "condition",
                               min_prop = 0.99),
    "No clusters pass"
  )
})

test_that("sw_differential_abundance adj_pvalue is BH corrected", {
  skip_if_not_installed("edgeR")

  x    <- .make_corrected()
  meta <- .make_meta()
  result <- sw_differential_abundance(x, meta, group_col = "condition",
                                       method = "edgeR", min_prop = 0)
  # BH is monotone: adj_pvalue >= pvalue
  expect_true(all(result$adj_pvalue >= result$pvalue - 1e-10))
})

# ---------------------------------------------------------------------------
# sw_differential_expression — validation
# ---------------------------------------------------------------------------

test_that("sw_differential_expression rejects non-data.frame corrected", {
  expect_error(
    sw_differential_expression("bad", .make_meta(), group_col = "condition"),
    "data.frame or tibble"
  )
})

test_that("sw_differential_expression rejects corrected without cluster col", {
  x <- tibble::tibble(sample = "S1", CD3 = 1)
  expect_error(
    sw_differential_expression(x, .make_meta(), group_col = "condition"),
    "'cluster' column"
  )
})

test_that("sw_differential_expression rejects corrected without sample col", {
  x <- tibble::tibble(cluster = 1L, CD3 = 1)
  expect_error(
    sw_differential_expression(x, .make_meta(), group_col = "condition"),
    "'sample' column"
  )
})

test_that("sw_differential_expression rejects meta missing group_col", {
  x    <- .make_corrected()
  meta <- data.frame(sample = paste0("S", 1:6), stringsAsFactors = FALSE)
  expect_error(
    sw_differential_expression(x, meta, group_col = "condition"),
    "missing required column"
  )
})

test_that("sw_differential_expression rejects markers not in corrected", {
  x    <- .make_corrected()
  meta <- .make_meta()
  expect_error(
    sw_differential_expression(x, meta, group_col = "condition",
                                markers = c("CD3", "NOT_REAL")),
    "Markers not found"
  )
})

test_that("sw_differential_expression rejects invalid min_cells", {
  x    <- .make_corrected()
  meta <- .make_meta()
  expect_error(
    sw_differential_expression(x, meta, group_col = "condition",
                                min_cells = 0L),
    "min_cells"
  )
})

# ---------------------------------------------------------------------------
# sw_differential_expression — integration (requires limma)
# ---------------------------------------------------------------------------

test_that("sw_differential_expression returns long tibble with correct columns", {
  skip_if_not_installed("limma")

  x    <- .make_corrected(n_per_group = 20)
  meta <- .make_meta()
  result <- sw_differential_expression(
    x, meta, group_col = "condition",
    markers   = c("CD3", "CD4", "CD8a"),
    min_cells = 5L
  )

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("cluster", "marker", "logFC", "pvalue", "adj_pvalue",
                     "mean_A", "mean_B") %in% names(result)))
  # 2 clusters × 3 markers = 6 rows (if both clusters have enough cells)
  expect_true(nrow(result) > 0L)
  expect_true(all(result$adj_pvalue >= 0 & result$adj_pvalue <= 1))
})

test_that("sw_differential_expression adj_pvalue exists as BH-corrected col", {
  skip_if_not_installed("limma")

  x    <- .make_corrected(n_per_group = 20)
  meta <- .make_meta()
  result <- sw_differential_expression(
    x, meta, group_col = "condition",
    markers = c("CD3", "CD4", "CD8a"), min_cells = 5L
  )
  expect_true("adj_pvalue" %in% names(result))
})

test_that("sw_differential_expression detects known group difference in CD4", {
  skip_if_not_installed("limma")

  # Use a larger effect size to ensure significance
  set.seed(99)
  n  <- 100L
  df <- do.call(rbind, list(
    # ctrl samples — cluster 1
    data.frame(sample = "S1", condition = "ctrl", cluster = 1L,
               CD4 = stats::rnorm(n, 0.3, 0.1), stringsAsFactors = FALSE),
    data.frame(sample = "S2", condition = "ctrl", cluster = 1L,
               CD4 = stats::rnorm(n, 0.3, 0.1), stringsAsFactors = FALSE),
    data.frame(sample = "S3", condition = "ctrl", cluster = 1L,
               CD4 = stats::rnorm(n, 0.3, 0.1), stringsAsFactors = FALSE),
    # trt samples — cluster 1, elevated CD4
    data.frame(sample = "S4", condition = "trt", cluster = 1L,
               CD4 = stats::rnorm(n, 3.0, 0.1), stringsAsFactors = FALSE),
    data.frame(sample = "S5", condition = "trt", cluster = 1L,
               CD4 = stats::rnorm(n, 3.0, 0.1), stringsAsFactors = FALSE),
    data.frame(sample = "S6", condition = "trt", cluster = 1L,
               CD4 = stats::rnorm(n, 3.0, 0.1), stringsAsFactors = FALSE)
  ))
  meta   <- data.frame(sample = paste0("S", 1:6),
                        condition = rep(c("ctrl", "trt"), each = 3),
                        stringsAsFactors = FALSE)
  result <- sw_differential_expression(
    tibble::as_tibble(df), meta,
    group_col = "condition",
    markers   = "CD4",
    min_cells = 5L
  )
  expect_true(nrow(result) >= 1L)
  expect_true(result$logFC[result$marker == "CD4"] > 0)
})

test_that("sw_differential_expression min_cells filter removes sparse pairs", {
  skip_if_not_installed("limma")

  x    <- .make_corrected(n_per_group = 5)
  meta <- .make_meta()
  # With min_cells = 100, no sample will contribute enough cells
  expect_error(
    sw_differential_expression(x, meta, group_col = "condition",
                                markers = c("CD3", "CD4"), min_cells = 100L),
    "enough samples"
  )
})

test_that("sw_differential_expression groups attribute is correct", {
  skip_if_not_installed("limma")

  x    <- .make_corrected(n_per_group = 20)
  meta <- .make_meta()
  result <- sw_differential_expression(
    x, meta, group_col = "condition",
    markers = c("CD3", "CD4"), min_cells = 5L
  )
  expect_equal(attr(result, "groups"), c("ctrl", "trt"))
})

# ---------------------------------------------------------------------------
# sw_plot_volcano — validation + basic output
# ---------------------------------------------------------------------------

test_that("sw_plot_volcano rejects non-data.frame de_result", {
  expect_error(sw_plot_volcano("bad"), "sw_differential_expression")
})

test_that("sw_plot_volcano rejects de_result missing required columns", {
  df <- data.frame(cluster = 1L, marker = "CD4", logFC = 1)
  expect_error(sw_plot_volcano(df), "columns")
})

test_that("sw_plot_volcano rejects unknown cluster", {
  df <- data.frame(
    cluster    = c(1L, 1L),
    marker     = c("CD3", "CD4"),
    logFC      = c(1, -1),
    pvalue     = c(0.01, 0.5),
    adj_pvalue = c(0.05, 0.5),
    stringsAsFactors = FALSE
  )
  expect_error(sw_plot_volcano(df, cluster = 99L), "No rows found")
})

test_that("sw_plot_volcano returns ggplot", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(
    cluster    = c(1L, 1L, 1L),
    marker     = c("CD3", "CD4", "CD8a"),
    logFC      = c(1.5, -0.3, 0.8),
    pvalue     = c(0.001, 0.4, 0.01),
    adj_pvalue = c(0.01, 0.5, 0.05),
    stringsAsFactors = FALSE
  )
  p <- sw_plot_volcano(df)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# sw_plot_boxplots — validation + basic output
# ---------------------------------------------------------------------------

test_that("sw_plot_boxplots rejects non-data.frame corrected", {
  expect_error(sw_plot_boxplots("bad", "CD3", "condition"), "data.frame")
})

test_that("sw_plot_boxplots rejects missing markers", {
  df <- data.frame(condition = "ctrl", cluster = 1L, stringsAsFactors = FALSE)
  expect_error(sw_plot_boxplots(df, "CD3", "condition"), "not found")
})

test_that("sw_plot_boxplots rejects missing group_col", {
  df <- data.frame(cluster = 1L, CD3 = 1, stringsAsFactors = FALSE)
  expect_error(sw_plot_boxplots(df, "CD3", "condition"), "not found")
})

test_that("sw_plot_boxplots rejects unknown cluster", {
  df <- data.frame(cluster = 1L, CD3 = 1, condition = "ctrl",
                   stringsAsFactors = FALSE)
  expect_error(sw_plot_boxplots(df, "CD3", "condition", cluster = 99L),
               "No rows found")
})

test_that("sw_plot_boxplots returns ggplot", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(
    cluster    = c(1L, 1L, 1L, 1L),
    condition  = c("ctrl", "ctrl", "trt", "trt"),
    CD3        = c(1, 1.2, 2, 2.3),
    CD4        = c(0.5, 0.7, 1.5, 1.8),
    stringsAsFactors = FALSE
  )
  p <- sw_plot_boxplots(df, c("CD3", "CD4"), group_col = "condition")
  expect_s3_class(p, "ggplot")
})

test_that("sw_plot_boxplots cluster filter works", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(
    cluster   = c(1L, 1L, 2L, 2L),
    condition = c("ctrl", "trt", "ctrl", "trt"),
    CD3       = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )
  # Should not error — only cluster 1 rows used
  p <- sw_plot_boxplots(df, "CD3", group_col = "condition", cluster = 1L)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# sw_plot_abundance — validation + basic output
# ---------------------------------------------------------------------------

test_that("sw_plot_abundance rejects non-data.frame da_result", {
  expect_error(sw_plot_abundance("bad"), "sw_differential_abundance")
})

test_that("sw_plot_abundance rejects da_result missing required columns", {
  df <- data.frame(cluster = 1L, logFC = 0.5, stringsAsFactors = FALSE)
  expect_error(sw_plot_abundance(df), "columns")
})

test_that("sw_plot_abundance(type='bar') returns ggplot", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(
    cluster     = c(1L, 2L),
    logFC       = c(0.5, -0.3),
    pvalue      = c(0.01, 0.2),
    adj_pvalue  = c(0.04, 0.3),
    mean_prop_A = c(0.3, 0.4),
    mean_prop_B = c(0.5, 0.3),
    stringsAsFactors = FALSE
  )
  attr(df, "groups")        <- c("ctrl", "trt")
  attr(df, "fdr_threshold") <- 0.05
  p <- sw_plot_abundance(df, type = "bar")
  expect_s3_class(p, "ggplot")
})

test_that("sw_plot_abundance(type='box') returns ggplot", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(
    cluster     = c(1L, 2L),
    logFC       = c(0.5, -0.3),
    pvalue      = c(0.01, 0.2),
    adj_pvalue  = c(0.04, 0.3),
    mean_prop_A = c(0.3, 0.4),
    mean_prop_B = c(0.5, 0.3),
    stringsAsFactors = FALSE
  )
  attr(df, "groups")        <- c("ctrl", "trt")
  attr(df, "fdr_threshold") <- 0.05
  p <- sw_plot_abundance(df, type = "box")
  expect_s3_class(p, "ggplot")
})

test_that("sw_plot_abundance uses fallback group names when attr missing", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(
    cluster     = 1L,
    logFC       = 0.5,
    pvalue      = 0.01,
    adj_pvalue  = 0.04,
    mean_prop_A = 0.3,
    mean_prop_B = 0.5,
    stringsAsFactors = FALSE
  )
  # No groups attribute — should not error
  p <- sw_plot_abundance(df)
  expect_s3_class(p, "ggplot")
})
