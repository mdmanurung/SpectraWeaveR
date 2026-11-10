# Tests for R/gate_wrap.R. Real round trips use flowCore and xml2; the
# flowGate, CytoML and flowMagic calls are mocked through .sw_gate_wrap_fn().

poly <- function() {
  sw_gate_polygon("lymph", c("FSC-A", "SSC-A"),
    rbind(c(0, 0), c(10, 0), c(10, 5), c(0, 5)),
    parent = NA
  )
}

test_that("to/from flowCore round-trips polygon and range gates", {
  skip_if_not_installed("flowCore")
  p <- poly()
  fp <- sw_gate_to_flowcore(p)
  expect_s4_class(fp, "polygonGate")
  back <- sw_gate_from_flowcore(fp, "lymph", NA)
  expect_equal(unname(back$vertices), unname(p$vertices))
  expect_equal(back$channels, p$channels)

  r <- sw_gate_range("cd3", "BV421-A", 1, 5, parent = "lymph")
  fr <- sw_gate_to_flowcore(r)
  expect_s4_class(fr, "rectangleGate")
  rb <- sw_gate_from_flowcore(fr, "cd3", "lymph")
  expect_equal(rb$type, "range")
  expect_equal(unname(rb$vertices[, 1]), c(1, 5))
  expect_equal(rb$parent, "lymph")
})

test_that("open-sided range and 2D rectangle convert", {
  skip_if_not_installed("flowCore")
  r <- sw_gate_range("pos", "X", 2, Inf, parent = NA)
  rb <- sw_gate_from_flowcore(sw_gate_to_flowcore(r), "pos", NA)
  expect_equal(unname(rb$vertices[, 1]), c(2, Inf))
  rect <- flowCore::rectangleGate(filterId = "box", "A" = c(0, 1), "B" = c(2, 3))
  g <- sw_gate_from_flowcore(rect, "box", NA)
  expect_equal(g$type, "polygon")
  expect_equal(nrow(g$vertices), 4)
  expect_error(sw_gate_from_flowcore(1, "x", NA), "polygonGate or rectangleGate")
  expect_error(sw_gate_to_flowcore("x"), "sw_gate")
})

test_that("GatingML write then read keeps gates, channels and parents", {
  skip_if_not_installed("xml2")
  gates <- list(
    lymph = poly(),
    cd3 = sw_gate_range("cd3", "BV421-A", -Inf, 5.5, parent = "lymph"),
    both = sw_gate_range("both", "A&B", 0.123456789012, 9, parent = NA)
  )
  f <- tempfile(fileext = ".xml")
  expect_identical(sw_gate_write_gatingml(gates, f), f)
  back <- sw_gate_read_gatingml(f)
  expect_named(back, c("lymph", "cd3", "both"))
  expect_equal(back$lymph$vertices, gates$lymph$vertices, ignore_attr = TRUE)
  expect_equal(back$lymph$channels, c("FSC-A", "SSC-A"))
  expect_true(is.na(back$lymph$parent))
  expect_equal(back$cd3$parent, "lymph")
  expect_equal(unname(back$cd3$vertices[, 1]), c(-Inf, 5.5))
  expect_equal(back$both$channels, "A&B")
  expect_equal(unname(back$both$vertices[1, 1]), 0.123456789012)
  expect_length(attr(back, "skipped"), 0)
})

test_that("GatingML writer validates input and reader skips unknown gates", {
  skip_if_not_installed("xml2")
  expect_error(sw_gate_write_gatingml(list(), tempfile()), "non-empty")
  expect_error(sw_gate_write_gatingml(list(poly(), poly()), tempfile()), "unique")
  f <- tempfile(fileext = ".xml")
  writeLines(c(
    "<gating:Gating-ML xmlns:gating=\"http://www.isac-net.org/std/Gating-ML/v2.0/gating\">",
    "<gating:EllipsoidGate gating:id=\"e\"/></gating:Gating-ML>"
  ), f)
  out <- sw_gate_read_gatingml(f)
  expect_length(out, 0)
  expect_equal(attr(out, "skipped"), "e")
})

test_that("missing Suggests package gives an install hint", {
  skip_if(requireNamespace("flowMagic", quietly = TRUE))
  expect_error(.sw_gate_wrap_fn("magic_template_train", "flowMagic"), "install_github")
})

fake_gs_fns <- function(gate_by_path) {
  function(name, pkg) {
    switch(name,
      sampleNames = function(x) c("s1", "s2"),
      gs_get_pop_paths = function(gs, path = "full") {
        if (path == "full") c("root", "/Lymph", "/Lymph/CD3") else c("root", "Lymph", "CD3")
      },
      gs_pop_get_gate = function(gs, y) {
        f <- gate_by_path[[y]]
        list(s1 = f, s2 = f)
      },
      stop("unmocked ", name)
    )
  }
}

test_that("sw_gate_read_flowjo passes required args and converts gates per sample", {
  skip_if_not_installed("flowCore")
  wsp <- tempfile(fileext = ".wsp")
  writeLines("x", wsp)
  seen <- list()
  gates <- list(
    "/Lymph" = sw_gate_to_flowcore(poly()),
    "/Lymph/CD3" = sw_gate_to_flowcore(sw_gate_range("cd3", "BV421-A", 1, 5, parent = "Lymph"))
  )
  base <- fake_gs_fns(gates)
  local_mocked_bindings(.sw_gate_wrap_fn = function(name, pkg) {
    if (name == "open_flowjo_xml") {
      return(function(p) {
        seen$path <<- p
        "ws"
      })
    }
    if (name == "flowjo_to_gatingset") {
      return(function(ws, ...) {
        seen$args <<- list(...)
        "gs"
      })
    }
    base(name, pkg)
  })
  out <- sw_gate_read_flowjo(wsp, group = "All Samples", transform = FALSE)
  expect_equal(seen$path, wsp)
  expect_false(seen$args$execute)
  expect_false(seen$args$transform)
  expect_equal(seen$args$name, "All Samples")
  expect_named(out, c("s1", "s2"))
  expect_named(out$s1, c("Lymph", "CD3"))
  expect_true(is.na(out$s1$Lymph$parent))
  expect_equal(out$s1$CD3$parent, "Lymph")
  expect_error(sw_gate_read_flowjo(wsp, group = "g"), "transform")
  expect_error(sw_gate_read_flowjo("nope.wsp", "g", FALSE), "does not exist")
})

test_that("sw_gate_draw calls flowGate with every choice and returns sw_gates", {
  skip_if_not_installed("flowCore")
  withr::local_options(sw.gate_draw_allow_noninteractive = TRUE)
  gs <- structure(list(), class = "GatingSet")
  seen <- NULL
  gates <- list("/Lymph" = sw_gate_to_flowcore(poly()))
  base <- function(name, pkg) {
    switch(name,
      sampleNames = function(x) c("s1", "s2"),
      gs_get_pop_paths = function(gs, path = "full") if (path == "full") c("root", "/Lymph") else c("root", "Lymph"),
      gs_pop_get_gate = function(gs, y) list(s1 = gates[[y]], s2 = gates[[y]])
    )
  }
  local_mocked_bindings(.sw_gate_wrap_fn = function(name, pkg) {
    if (name == "gs_gate_interactive") {
      return(function(gs, ...) {
        seen <<- list(...)
        invisible(gs)
      })
    }
    base(name, pkg)
  })
  out <- sw_gate_draw(gs, "Lymph", list("FSC-A", "SSC-A"),
    sample = 1, subset = "root", regate = FALSE
  )
  expect_equal(seen$filterId, "Lymph")
  expect_equal(seen$subset, "root")
  expect_false(seen$regate)
  expect_equal(seen$dims, list("FSC-A", "SSC-A"))
  expect_s3_class(out$s2$Lymph, "sw_gate")
  expect_identical(attr(out, "gs"), gs)
  expect_error(sw_gate_draw(gs, "Lymph", list("a", "b"), sample = 1, subset = "root"), "regate")
})

test_that("sw_gate_draw refuses a non-interactive session", {
  skip_if(interactive())
  expect_error(sw_gate_draw(1, "a", list("x", "y"), 1, "root", FALSE), "interactive")
})

test_that("sw_gate_flowmagic returns sw_gates per sample and node", {
  skip_if_not_installed("flowCore")
  mk <- function(off) {
    flowCore::polygonGate(
      .gate = cbind("FSC-A" = c(0, 4, 4, 0) + off, "SSC-A" = c(0, 0, 4, 4)),
      filterId = "Live"
    )
  }
  seen <- NULL
  local_mocked_bindings(.sw_gate_wrap_fn = function(name, pkg) {
    stopifnot(name == "magic_template_predict", pkg == "flowMagic")
    function(...) {
      seen <<- list(...)
      list(gs = "gs2", list_poly_gates = list(Live = list(s1 = mk(0), s2 = mk(1))))
    }
  })
  out <- sw_gate_flowmagic("gs", "root", "FSC-A", "SSC-A", ref_train = 1, ref_model_info = 2,
    label_to_node = c("1" = "Live"), concavity_val = 50
  )
  expect_equal(seen$concavity_val, 50)
  expect_true(seen$add_to_gs)
  expect_equal(seen$label_to_node, c("1" = "Live"))
  expect_named(out, c("s1", "s2"))
  expect_s3_class(out$s2$Live, "sw_gate")
  expect_true(is.na(out$s1$Live$parent))
  expect_equal(out$s2$Live$vertices[1, 1], 1, ignore_attr = TRUE)
  expect_equal(attr(out, "gs"), "gs2")
  inside <- sw_gate_apply(out$s1$Live, cbind("FSC-A" = c(1, 9), "SSC-A" = c(1, 1)))
  expect_equal(inside, c(TRUE, FALSE))
  expect_error(sw_gate_flowmagic("gs", "root", "a", "b", 1, 2, label_to_node = "Live", concavity_val = 5), "named")
  expect_error(sw_gate_flowmagic("gs", "root", "a", "b", 1, 2, label_to_node = c("1" = "L")), "concavity_val")
})

test_that("sw_gate_flowmagic_train needs exactly one downsampling choice", {
  seen <- NULL
  local_mocked_bindings(.sw_gate_wrap_fn = function(name, pkg) {
    function(...) {
      seen <<- list(...)
      list(ref_model_info = 1, ref_train = 2)
    }
  })
  args <- list(
    gs = "gs", parent_node = "root", samples_id = "s1", channel_x = "a", channel_y = "b",
    n_gates = 1, train_model = "rf", seed = 1, input_template_obj = NULL,
    use_existing_templates = FALSE
  )
  expect_error(do.call(sw_gate_flowmagic_train, c(args, list(n_points_per_plot = NULL, prop_down = NULL))), "exactly one")
  expect_error(do.call(sw_gate_flowmagic_train, c(args, list(n_points_per_plot = 100, prop_down = 0.5))), "exactly one")
  res <- do.call(sw_gate_flowmagic_train, c(args, list(n_points_per_plot = 100, prop_down = NULL)))
  expect_equal(seen$n_points_per_plot, 100)
  expect_equal(seen$seed, 1)
  expect_equal(seen$train_model, "rf")
  expect_equal(res$ref_train, 2)
  expect_error(sw_gate_flowmagic_train("gs", "root", "s1", "a", "b"), "n_gates")
})
