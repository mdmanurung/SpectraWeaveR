# Opt-in tests: the gate_wrap.R wrappers against the real flowMagic and CytoML.
# No mocks. Each test skips when its package is missing.
# Run: testthat::test_file("tests/testthat/test-gate_wrap_real.R")

# Three samples, two blobs each (A near 1,1 and B near 4,4), as a GatingSet.
make_blob_gs <- function(seed = 1) {
  set.seed(seed)
  ffs <- lapply(1:3, function(i) {
    m <- rbind(
      cbind(rnorm(300, 1, 0.3), rnorm(300, 1, 0.3)),
      cbind(rnorm(300, 4, 0.3), rnorm(300, 4, 0.3))
    )
    colnames(m) <- c("FSC-A", "SSC-A")
    flowCore::flowFrame(m)
  })
  names(ffs) <- paste0("s", 1:3)
  fs <- flowCore::flowSet(ffs)
  flowWorkspace::GatingSet(fs)
}

# A template object as magic_template_train() would build it, with the polygon
# fixed in code instead of drawn: label "1" is blob A.
make_template <- function(gs, samples) {
  poly <- data.frame(x = c(0, 2, 2, 0), y = c(0, 0, 2, 2))
  gated <- lapply(samples, function(s) {
    df <- as.data.frame(flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[s]], "root")))
    df <- df[, c("FSC-A", "SSC-A")]
    flowMagic::magic_label_poly(df, poly, label_pol = "1")
  })
  names(gated) <- samples
  polys <- flowMagic::flowmagic_pred_to_poly_gates(
    list_df = gated, pred_label = "1", gate_label = "1", concavity_val = 10
  )
  list(all_list_out_obj = list(list(list_poly_gates = polys, list_gated_data = gated)))
}

test_that("sw_gate_flowmagic_train and sw_gate_flowmagic run on real flowMagic", {
  skip_if_not_installed("flowMagic")
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("sp")

  # flowMagic (0.99, 4c27293) calls sampleNames(), exprs(), recompute() and dplyr
  # verbs unqualified; the wrappers attach flowCore, flowWorkspace and dplyr
  # themselves, so this test does not.
  gs <- make_blob_gs()
  tmpl <- make_template(gs, c("s1", "s2"))

  trained <- suppressMessages(suppressWarnings(utils::capture.output(
    res <- sw_gate_flowmagic_train(
      gs,
      parent_node = "root", samples_id = c("s1", "s2"),
      channel_x = "FSC-A", channel_y = "SSC-A", n_gates = 1,
      train_model = "rf", n_points_per_plot = NULL, prop_down = 1,
      seed = 1, input_template_obj = tmpl, use_existing_templates = TRUE
    )
  )))
  expect_true(all(c("ref_model_info", "ref_train") %in% names(res)))

  gs2 <- flowWorkspace::gs_clone(gs)
  out <- NULL
  utils::capture.output(suppressMessages(
    out <- sw_gate_flowmagic(
      gs2,
      parent_node = "root", channel_x = "FSC-A", channel_y = "SSC-A",
      ref_train = res$ref_train, ref_model_info = res$ref_model_info,
      label_to_node = c("1" = "BlobA"), concavity_val = 10
    )
  ))
  expect_setequal(names(out), c("s1", "s2", "s3"))
  g <- out[["s3"]][["BlobA"]]
  expect_s3_class(g, "sw_gate")
  # BlobA sits near (1, 1): the predicted polygon must contain it, not blob B.
  inside <- sw_gate_apply(g, data.frame(`FSC-A` = c(1, 4), `SSC-A` = c(1, 4), check.names = FALSE))
  expect_equal(unname(as.logical(inside)), c(TRUE, FALSE))
})

test_that("sw_gate_read_flowjo reads a real FlowJo workspace", {
  skip_if_not_installed("CytoML")
  skip_if_not_installed("flowWorkspace")
  # CytoML ships no .wsp; FlowSOM ships extdata/gating.wsp (read with execute = FALSE).
  wsp <- list.files(
    c(system.file("extdata", package = "CytoML"), system.file("extdata", package = "FlowSOM")),
    pattern = "[.]wsp$", recursive = TRUE, full.names = TRUE
  )
  skip_if(length(wsp) == 0, "no .wsp file found in CytoML or FlowSOM")
  ws <- CytoML::open_flowjo_xml(wsp[[1]])
  grp <- CytoML::fj_ws_get_sample_groups(ws)
  grp_name <- unique(grp$groupName)[[1]]
  out <- sw_gate_read_flowjo(wsp[[1]], group = grp_name, transform = FALSE)
  expect_type(out, "list")
  expect_gt(length(out), 0)
})
