# tests/testthat/test-annotate.R
# Unit tests for R/annotate.R — Cell Type Annotation

# ---------------------------------------------------------------------------
# Helper: minimal synthetic MFI tibble
# ---------------------------------------------------------------------------
.make_mfi <- function(n_clusters = 5) {
  # Each row = one cluster, columns = markers
  # Marker names match the built-in reference
  tibble::tibble(
    cluster    = seq_len(n_clusters),
    CD3        = c(3, 3, 0, 0, 3),
    CD4        = c(3, 0, 0, 0, 3),
    CD8a       = c(0, 3, 0, 0, 0),
    CD19       = c(0, 0, 3, 0, 0),
    CD20       = c(0, 0, 3, 0, 0),
    CD56       = c(0, 0, 0, 2, 0),
    CD16       = c(0, 0, 0, 2, 0),
    CD14       = c(0, 0, 0, 0, 0),
    HLA_DR     = c(0, 0, 2, 0, 0),
    TCRab      = c(3, 3, 0, 0, 3),
    CD25       = c(1, 0, 0, 0, 3),
    CD127      = c(2, 2, 0, 0, 1),
    FoxP3      = c(0, 0, 0, 0, 3)
  )
}

# ---------------------------------------------------------------------------
# sw_load_reference
# ---------------------------------------------------------------------------

test_that("sw_load_reference returns a list with matrix and mask", {
  ref <- sw_load_reference("pbmc")
  expect_type(ref, "list")
  expect_named(ref, c("matrix", "mask"))
  expect_true(is.matrix(ref$matrix))
  expect_true(is.matrix(ref$mask))
})

test_that("sw_load_reference matrix and mask have identical dimensions", {
  ref <- sw_load_reference("pbmc")
  expect_identical(dim(ref$matrix), dim(ref$mask))
  expect_identical(dimnames(ref$matrix), dimnames(ref$mask))
})

test_that("sw_load_reference reference has expected populations", {
  ref <- sw_load_reference("pbmc")
  pops <- rownames(ref$matrix)
  expect_true("CD4_T"     %in% pops)
  expect_true("CD8_T"     %in% pops)
  expect_true("B_cell"    %in% pops)
  expect_true("NK"        %in% pops)
  expect_true("CD14_Mono" %in% pops)
  expect_true("Treg"      %in% pops)
  expect_gte(nrow(ref$matrix), 30L)
})

test_that("sw_load_reference matrix values are in 0-3 range", {
  ref <- sw_load_reference("pbmc")
  expect_true(all(ref$matrix >= 0))
  expect_true(all(ref$matrix <= 3))
})

test_that("sw_load_reference mask values are binary 0/1", {
  ref <- sw_load_reference("pbmc")
  expect_true(all(ref$mask %in% c(0L, 1L)))
})

test_that("sw_load_reference rejects unsupported name", {
  expect_error(sw_load_reference("t_cell"), "should be one of")
})

test_that("sw_load_reference rejects non-character name", {
  expect_error(sw_load_reference(1L), "single non-empty character")
})

# ---------------------------------------------------------------------------
# .check_annotation_input (via sw_annotate_clusters)
# ---------------------------------------------------------------------------

test_that("sw_annotate_clusters rejects non-data.frame non-cluster_result", {
  expect_error(sw_annotate_clusters("not_valid"), "sw_cluster_result")
})

test_that("sw_annotate_clusters rejects data.frame missing cluster column", {
  df <- data.frame(CD3 = 1, CD4 = 2)
  expect_error(sw_annotate_clusters(df), "'cluster' column")
})

test_that("sw_annotate_clusters rejects data.frame with no numeric markers", {
  df <- data.frame(cluster = 1L, label = "A", stringsAsFactors = FALSE)
  expect_error(sw_annotate_clusters(df), "numeric marker")
})

test_that("sw_annotate_clusters rejects invalid reference structure", {
  mfi <- .make_mfi(2)
  expect_error(
    sw_annotate_clusters(mfi, reference = list(not_matrix = 1)),
    "'matrix' and 'mask'"
  )
})

test_that("sw_annotate_clusters rejects non-matrix reference$matrix", {
  mfi <- .make_mfi(2)
  ref <- list(matrix = data.frame(CD3 = 1), mask = matrix(1L, 1, 1))
  expect_error(sw_annotate_clusters(mfi, reference = ref), "must be matrices")
})

test_that("sw_annotate_clusters rejects min_score outside [0,1]", {
  mfi <- .make_mfi(2)
  expect_error(sw_annotate_clusters(mfi, min_score = -0.1), "\\[0, 1\\]")
  expect_error(sw_annotate_clusters(mfi, min_score = 1.5),  "\\[0, 1\\]")
})

test_that("sw_annotate_clusters rejects markers not in MFI data", {
  mfi <- .make_mfi(2)
  expect_error(
    sw_annotate_clusters(mfi, markers = c("CD3", "NOT_A_MARKER")),
    "not found"
  )
})

test_that("sw_annotate_clusters rejects too few overlapping markers", {
  mfi <- tibble::tibble(cluster = 1L, XYZ_fake = 1)
  mat  <- matrix(1, 1, 1, dimnames = list("Pop", "ABC_fake"))
  mask <- matrix(1L, 1, 1, dimnames = list("Pop", "ABC_fake"))
  ref  <- list(matrix = mat, mask = mask)
  expect_error(sw_annotate_clusters(mfi, reference = ref),
               "Fewer than 2 markers")
})

# ---------------------------------------------------------------------------
# sw_annotate_clusters — correctness
# ---------------------------------------------------------------------------

test_that("sw_annotate_clusters returns tibble with expected columns", {
  mfi    <- .make_mfi(5)
  result <- sw_annotate_clusters(mfi)
  expect_s3_class(result, "tbl_df")
  expect_named(result,
    c("cluster", "cell_type", "score", "second_type", "second_score",
      "n_markers_used"),
    ignore.order = FALSE
  )
  expect_equal(nrow(result), 5L)
})

test_that("sw_annotate_clusters assigns CD4-related type to a CD4+ cluster", {
  mfi    <- .make_mfi(1)  # first row: CD3+ CD4+ TCRab+ pattern
  result <- sw_annotate_clusters(mfi, min_score = 0)
  expect_false(is.na(result$cell_type[1L]))
  # Should match a T cell or CD4 subtype
  expect_true(
    grepl("CD4|Treg|T", result$cell_type[1L], ignore.case = TRUE)
  )
})

test_that("sw_annotate_clusters assigns B-cell type to a B-cell cluster", {
  mfi <- tibble::tibble(
    cluster = 1L,
    CD3 = 0, CD4 = 0, CD8a = 0,
    CD19 = 3, CD20 = 3, CD56 = 0, CD16 = 0, CD14 = 0,
    HLA_DR = 2, TCRab = 0, CD25 = 0, CD127 = 0, FoxP3 = 0
  )
  result <- sw_annotate_clusters(mfi, min_score = 0)
  expect_false(is.na(result$cell_type[1L]))
  expect_true(grepl("B", result$cell_type[1L], ignore.case = TRUE))
})

test_that("sw_annotate_clusters sets NA when score below min_score", {
  mfi <- tibble::tibble(
    cluster = 1L,
    CD3 = 1.5, CD4 = 1.5, CD19 = 1.5, CD56 = 1.5, CD14 = 1.5
  )
  # Using a very high min_score should force NA
  result <- sw_annotate_clusters(mfi, min_score = 0.9999)
  expect_true(is.na(result$cell_type[1L]))
})

test_that("sw_annotate_clusters respects user-supplied markers subset", {
  mfi    <- .make_mfi(3)
  result <- sw_annotate_clusters(mfi, markers = c("CD3", "CD4", "CD8a",
                                                    "CD19", "CD20"))
  expect_equal(result$n_markers_used[1L], 5L)
})

test_that("sw_annotate_clusters works with a custom reference", {
  mfi <- tibble::tibble(
    cluster = 1L,
    MarkerA = 3,
    MarkerB = 0,
    MarkerC = 3
  )
  mat  <- matrix(c(3, 0, 3, 0, 3, 0), nrow = 2,
                 dimnames = list(c("Pop1", "Pop2"),
                                 c("MarkerA", "MarkerB", "MarkerC")))
  mask <- matrix(c(1L, 1L, 1L, 1L, 1L, 1L), nrow = 2,
                 dimnames = list(c("Pop1", "Pop2"),
                                 c("MarkerA", "MarkerB", "MarkerC")))
  ref    <- list(matrix = mat, mask = mask)
  result <- sw_annotate_clusters(mfi, reference = ref, min_score = 0)
  expect_equal(result$cell_type[1L], "Pop1")
})

test_that("sw_annotate_clusters score is in [0, 1]", {
  mfi    <- .make_mfi(5)
  result <- sw_annotate_clusters(mfi, min_score = 0)
  expect_true(all(result$score >= 0, na.rm = TRUE))
  expect_true(all(result$score <= 1, na.rm = TRUE))
})

test_that("sw_annotate_clusters score >= second_score", {
  mfi    <- .make_mfi(5)
  result <- sw_annotate_clusters(mfi, min_score = 0)
  expect_true(all(result$score >= result$second_score, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# sw_annotate_manual
# ---------------------------------------------------------------------------

test_that("sw_annotate_manual rejects non-named annotation_map", {
  mfi <- .make_mfi(2)
  expect_error(sw_annotate_manual(mfi, c("CD4+ T", "B cell")),
               "named character vector")
})

test_that("sw_annotate_manual rejects non-character annotation_map", {
  mfi <- .make_mfi(2)
  expect_error(sw_annotate_manual(mfi, c("1" = 1L, "2" = 2L)),
               "named character vector")
})

test_that("sw_annotate_manual appends cell_type column", {
  mfi <- .make_mfi(3)
  map <- c("1" = "CD4+ T", "2" = "CD8+ T", "3" = "B cell")
  result <- sw_annotate_manual(mfi, map)
  expect_true("cell_type" %in% names(result))
  expect_equal(result$cell_type[1L], "CD4+ T")
  expect_equal(result$cell_type[2L], "CD8+ T")
  expect_equal(result$cell_type[3L], "B cell")
})

test_that("sw_annotate_manual places cell_type after cluster", {
  mfi    <- .make_mfi(2)
  map    <- c("1" = "CD4+ T", "2" = "CD8+ T")
  result <- sw_annotate_manual(mfi, map)
  expect_equal(names(result)[1L], "cluster")
  expect_equal(names(result)[2L], "cell_type")
})

test_that("sw_annotate_manual sets NA with warning for unmapped clusters", {
  mfi <- .make_mfi(3)
  map <- c("1" = "CD4+ T")  # clusters 2 and 3 not in map
  expect_warning(
    result <- sw_annotate_manual(mfi, map),
    "not found in 'annotation_map'"
  )
  expect_true(is.na(result$cell_type[2L]))
  expect_true(is.na(result$cell_type[3L]))
})

# ---------------------------------------------------------------------------
# sw_plot_annotation
# ---------------------------------------------------------------------------

test_that("sw_plot_annotation rejects non-data.frame annotation", {
  expect_error(sw_plot_annotation("not_df"), "data.frame or tibble")
})

test_that("sw_plot_annotation rejects annotation without cluster column", {
  df <- data.frame(cell_type = "A")
  expect_error(sw_plot_annotation(df), "'cluster' column")
})

test_that("sw_plot_annotation(type='umap') requires dimred", {
  ann <- tibble::tibble(cluster = 1L, cell_type = "CD4+ T")
  expect_error(sw_plot_annotation(ann, type = "umap"),
               "'dimred' must be provided")
})

test_that("sw_plot_annotation(type='umap') rejects dimred without dim1/dim2", {
  ann  <- tibble::tibble(cluster = 1L, cell_type = "CD4+ T")
  bad  <- tibble::tibble(cluster = 1L, x = 1, y = 2)
  expect_error(sw_plot_annotation(ann, dimred = bad, type = "umap"),
               "dim1.*dim2")
})

test_that("sw_plot_annotation(type='umap') rejects dimred without cluster", {
  ann   <- tibble::tibble(cluster = 1L, cell_type = "CD4+ T")
  bad   <- tibble::tibble(dim1 = 0, dim2 = 0)
  expect_error(sw_plot_annotation(ann, dimred = bad, type = "umap"),
               "'cluster' column")
})

test_that("sw_plot_annotation(type='heatmap') needs cluster_result when no markers in annotation", {
  ann <- tibble::tibble(cluster = 1L, cell_type = "CD4+ T",
                        score = 0.9, second_type = "Treg",
                        second_score = 0.5, n_markers_used = 5L)
  expect_error(
    sw_plot_annotation(ann, type = "heatmap"),
    "cluster_result.*must be provided"
  )
})

test_that("sw_plot_annotation(type='umap') returns ggplot when ggplot2 available", {
  skip_if_not_installed("ggplot2")
  ann <- tibble::tibble(cluster = c(1L, 1L, 2L, 2L),
                        cell_type = c("CD4+ T", "CD4+ T", "NK", "NK"))
  dr  <- tibble::tibble(cluster = c(1L, 1L, 2L, 2L),
                        dim1 = c(1, 2, 5, 6), dim2 = c(1, 2, 5, 6))
  p <- sw_plot_annotation(ann, dimred = dr, type = "umap")
  expect_s3_class(p, "ggplot")
})

test_that("sw_plot_annotation(type='heatmap') returns a plot when given MFI annotation", {
  skip_if_not_installed("ggplot2")
  mfi <- .make_mfi(3)
  map <- c("1" = "CD4+ T", "2" = "CD8+ T", "3" = "B cell")
  ann <- sw_annotate_manual(mfi, map)
  p   <- sw_plot_annotation(ann, type = "heatmap")
  # pheatmap returns an S3 pheatmap object; ggplot2 fallback returns ggplot
  expect_true(inherits(p, "pheatmap") || inherits(p, "ggplot"))
})
