# tests/testthat/test-batch_scores.R
# Unit tests for R/batch_scores.R on simulated batches.

# Three batches of 5 files; marker A is shifted in batch b2, marker B is not.
make_batches <- function(shift = 1, n = 300, seed = 1) {
  set.seed(seed)
  ids <- paste0("f", 1:15)
  b <- stats::setNames(rep(c("b1", "b2", "b3"), each = 5), ids)
  d <- lapply(ids, function(i) {
    cbind(
      A = stats::rnorm(n) + if (b[[i]] == "b2") shift else 0,
      B = stats::rnorm(n)
    )
  })
  names(d) <- ids
  list(data = d, batch = b)
}

scores <- function(x, ...) {
  args <- list(
    data = x$data,
    batch = x$batch,
    markers = c("A", "B"),
    n_events = 200,
    seed = 1,
    order = NULL,
    reference_group = NULL,
    population = NULL
  )
  args[names(list(...))] <- list(...)
  do.call(sw_batch_scores, args)
}

get <- function(res, check, marker, batch = NA) {
  r <- res[res$check == check & res$marker == marker, ]
  if (!is.na(batch)) {
    r <- r[r$batch == batch, ]
  }
  r$value
}

test_that("required arguments have no default", {
  x <- make_batches()
  expect_error(sw_batch_scores(), "no default")
  expect_error(
    sw_batch_scores(
      x$data,
      x$batch,
      c("A", "B"),
      n_events = 100,
      seed = 1,
      reference_group = NULL
    ),
    "order"
  )
})

test_that("a shifted batch shows in EMD, median range and eta2 for that marker only", {
  res <- scores(make_batches(shift = 1))
  expect_s3_class(res, "sw_batch_scores")
  # EMD of a pure +1 shift between equal-shape samples is close to 1.
  expect_equal(get(res, "pair_emd", "A", "b1 ~ b2"), 1, tolerance = 0.15)
  expect_lt(get(res, "pair_emd", "A", "b1 ~ b3"), 0.3)
  expect_lt(get(res, "pair_emd", "B", "b1 ~ b2"), 0.3)
  expect_gt(get(res, "median_range", "A"), 0.8)
  expect_lt(get(res, "median_range", "B"), 0.3)
  expect_gt(get(res, "eta2", "A"), 0.8)
  expect_lt(get(res, "eta2", "B"), 0.5)
})

test_that("no batch effect gives small scores", {
  res <- scores(make_batches(shift = 0))
  expect_lt(max(res$value[res$check == "pair_emd"]), 0.3)
  expect_lt(max(res$value[res$check == "median_range"]), 0.3)
  expect_lt(max(abs(log(res$value[res$check == "spread_ratio"]))), 0.1)
})

test_that("a widened batch shows in spread_ratio for that batch and marker only", {
  x <- make_batches(shift = 0)
  for (i in names(x$batch)[x$batch == "b2"]) {
    v <- x$data[[i]][, "A"]
    x$data[[i]][, "A"] <- mean(v) + (v - mean(v)) * 1.5
  }
  res <- scores(x)
  sr <- res[res$check == "spread_ratio", ]
  a_b2 <- sr$value[sr$marker == "A" & sr$batch == "b2"]
  expect_gt(a_b2, 1.3)
  expect_lt(sr$value[sr$marker == "A" & sr$batch == "b1"], 1.1)
  expect_lt(sr$value[sr$marker == "B" & sr$batch == "b2"], 1.1)
  # The median does not move, so median_range stays small for A.
  expect_lt(get(res, "median_range", "A"), 0.3)
})

test_that("order drift inside a batch is found, and small batches are skipped", {
  x <- make_batches(shift = 0)
  # Add a steady drift to marker A across the 5 files of b1.
  ids <- names(x$batch)[x$batch == "b1"]
  for (k in seq_along(ids)) {
    x$data[[ids[k]]][, "A"] <- x$data[[ids[k]]][, "A"] + 0.5 * k
  }
  ord <- stats::setNames(as.numeric(seq_along(x$batch)), names(x$batch))
  res <- scores(x, order = ord)
  expect_equal(get(res, "order_rho", "A", "b1"), 1)
  expect_lt(abs(get(res, "order_rho", "A", "b3")), 1)
  # A batch with fewer than 4 files gets no order_rho.
  y <- make_batches()
  y$batch[names(y$batch)[y$batch == "b3"][1:2]] <- "b4"
  res2 <- scores(
    y,
    order = stats::setNames(as.numeric(seq_along(y$batch)), names(y$batch))
  )
  expect_false("b4" %in% res2$batch[res2$check == "order_rho"])
})

test_that("input problems stop with a message", {
  x <- make_batches()
  expect_error(scores(x, n_events = 5000), "fewer than `n_events`")
  expect_error(scores(x, markers = "Z"), "lacks marker")
  expect_error(scores(x, batch = x$batch[1:3]), "covering every file")
  expect_error(
    scores(x, batch = stats::setNames(rep("b1", 15), names(x$batch))),
    "two batches"
  )
})

test_that("results do not depend on file order and keep the RNG state", {
  x <- make_batches()
  set.seed(5)
  a <- stats::runif(1)
  set.seed(5)
  r1 <- scores(x)
  expect_equal(stats::runif(1), a)
  x2 <- list(
    data = x$data[rev(names(x$data))],
    batch = x$batch[rev(names(x$batch))]
  )
  r2 <- scores(x2)
  expect_equal(r1$value, r2$value)
})

# reference_group ---------------------------------------------------------------

# Four batches of three files, three donors with offsets 0, 4 and 8 on A.
# Donors are spread unevenly over batches (b1: D1 D1 D2, b2: D2 D2 D3,
# b3: D3 D3 D1, b4: D1 D2 D3), so pooled batches differ by donor mix. Batch b4
# carries a real shift of 0.5 in A for every donor.
make_donor_batches <- function(seed = 5, n = 300) {
  set.seed(seed)
  donor_of <- c(
    "D1",
    "D1",
    "D2",
    "D2",
    "D2",
    "D3",
    "D3",
    "D3",
    "D1",
    "D1",
    "D2",
    "D3"
  )
  batch_of <- rep(c("b1", "b2", "b3", "b4"), each = 3)
  ids <- paste0("f", seq_along(donor_of))
  d <- lapply(seq_along(ids), function(i) {
    off <- c(D1 = 0, D2 = 4, D3 = 8)[[donor_of[i]]]
    cbind(
      A = stats::rnorm(n) + off + if (batch_of[i] == "b4") 0.5 else 0,
      B = stats::rnorm(n)
    )
  })
  names(d) <- ids
  list(
    data = d,
    batch = stats::setNames(batch_of, ids),
    donor = stats::setNames(donor_of, ids)
  )
}

test_that("reference_group finds a batch effect hidden by donor differences", {
  x <- make_donor_batches()
  pair_a <- function(res) {
    r <- res[res$check == "pair_emd" & res$marker == "A", ]
    stats::setNames(r$value, r$batch)
  }
  all_files <- pair_a(scores(x))
  by_donor <- pair_a(scores(x, reference_group = x$donor))
  with_b4 <- grepl("b4", names(by_donor))
  # Within donors, every pair that has b4 is farther than every pair without.
  expect_gt(min(by_donor[with_b4]), max(by_donor[!with_b4]))
  # Over all files the donor mix swamps it.
  b4_all <- grepl("b4", names(all_files))
  expect_lt(min(all_files[b4_all]), max(all_files[!b4_all]))
  # All score types are still reported.
  res <- scores(x, reference_group = x$donor)
  expect_setequal(unique(res$check), unique(scores(x)$check))
})

test_that("batch reference_group is required, covers every file and needs 2 batches in a group", {
  x <- make_donor_batches()
  expect_error(
    sw_batch_scores(
      x$data,
      x$batch,
      c("A", "B"),
      n_events = 100,
      seed = 1,
      order = NULL,
      population = NULL
    ),
    "no default"
  )
  expect_error(scores(x, reference_group = x$donor[-1]), "label for every file")
  # One file per group and batch: no group has two batches.
  one_each <- stats::setNames(names(x$data), names(x$data))
  expect_error(
    scores(x, reference_group = one_each),
    "at least two batches"
  )
})

# reference_group, two-way (data frame) --------------------------------------------------

test_that("a two-way reference_group is an error for batch scores", {
  x <- make_donor_batches()
  design <- data.frame(
    file = names(x$data),
    donor = unname(x$donor),
    batch = unname(x$batch),
    stringsAsFactors = FALSE
  )
  expect_error(scores(x, reference_group = design), "remove the batch effect")
})

# 18 batches of 2 files; marker A is shifted by `shift` in batch "b07".
make_many <- function(shift = 1, n = 1000, seed = 2) {
  set.seed(seed)
  ids <- paste0("f", 1:36)
  b <- stats::setNames(rep(sprintf("b%02d", 1:18), each = 2), ids)
  d <- lapply(ids, function(i) {
    cbind(
      A = stats::rnorm(n) + if (b[[i]] == "b07") shift else 0,
      B = stats::rnorm(n)
    )
  })
  names(d) <- ids
  list(data = d, batch = b)
}

test_that("rest_emd of a shifted batch is about the shift", {
  x <- make_many(shift = 1)
  res <- scores(x, n_events = 1000)
  expect_equal(get(res, "rest_emd", "A", "b07"), 1, tolerance = 0.1)
  expect_lt(get(res, "rest_emd", "B", "b07"), 0.2)
})

test_that("rest_emd_z is finite for every marker and largest for the shifted batch", {
  x <- make_many(shift = 1)
  res <- scores(x, n_events = 1000)
  z <- res[res$check == "rest_emd_z", ]
  expect_false(anyNA(z$value))
  expect_equal(nrow(z), 2 * 18)
  za <- z[z$marker == "A", ]
  expect_equal(za$batch[which.max(za$value)], "b07")
})

test_that("rest_emd_z works with reference_group and is NA below 3 batches", {
  x <- make_many(shift = 1)
  donor <- stats::setNames(rep(c("d1", "d2"), 18), names(x$data))
  res <- scores(x, n_events = 1000, reference_group = donor)
  za <- res[res$check == "rest_emd_z" & res$marker == "A", ]
  expect_equal(za$batch[which.max(za$value)], "b07")
  two <- make_batches()
  keep <- two$batch %in% c("b1", "b2")
  small <- scores(list(data = two$data[keep], batch = two$batch[keep]))
  expect_true(all(is.na(small$value[small$check == "rest_emd_z"])))
})

# Population rows -----------------------------------------------------------------

# Two populations per file: p1 (first 60% of events) and p2. Marker A is shifted
# in batch b2 only inside p2; p2 is also rarer in b3. Order rises with file index.
make_pop <- function(n = 400, seed = 2) {
  set.seed(seed)
  ids <- paste0("f", 1:15)
  b <- stats::setNames(rep(c("b1", "b2", "b3"), each = 5), ids)
  labs <- list()
  d <- lapply(ids, function(i) {
    p2_share <- c(b1 = 0.4, b2 = 0.44, b3 = 0.15)[[b[[i]]]] +
      0.01 * (match(i, ids) %% 5 - 2)
    lab <- rep(c("p1", "p2"), c(round(n * (1 - p2_share)), round(n * p2_share)))
    shift <- if (b[[i]] == "b2") 2 else 0
    a <- stats::rnorm(n) + ifelse(lab == "p2", shift, 0)
    labs[[i]] <<- lab
    cbind(A = a, B = stats::rnorm(n))
  })
  names(d) <- ids
  list(
    data = d,
    batch = b,
    labels = labs,
    order = stats::setNames(seq_along(ids), ids)
  )
}

pop_scores <- function(x, ...) {
  args <- list(
    data = x$data,
    batch = x$batch,
    markers = c("A", "B"),
    n_events = 200,
    seed = 1,
    order = NULL,
    reference_group = NULL,
    population = list(labels = x$labels, label_source = "gate"),
    min_events = 50,
    group_summary = "median"
  )
  args[names(list(...))] <- list(...)
  do.call(sw_batch_scores, args)
}

test_that("population = NULL leaves the output identical to the all-event rows", {
  x <- make_pop()
  base <- scores(x, order = x$order)
  expect_false(any(c("population", "status") %in% names(base)))
  expect_false("batch_order_rho" %in% base$check)
  withpop <- pop_scores(x, order = x$order)
  same <- withpop[
    is.na(withpop$population) & withpop$check != "batch_order_rho",
  ]
  expect_equal(
    as.data.frame(same[, names(base)]),
    as.data.frame(base),
    ignore_attr = TRUE
  )
})

test_that("population arguments are required when population is set", {
  x <- make_pop()
  expect_error(
    scores(x, population = list(labels = x$labels, label_source = "gate")),
    "min_events"
  )
  args <- list(
    x$data,
    x$batch,
    c("A", "B"),
    n_events = 100,
    seed = 1,
    order = NULL,
    reference_group = NULL,
    population = list(labels = x$labels, label_source = "gate"),
    min_events = 50
  )
  expect_error(do.call(sw_batch_scores, args), "group_summary")
  expect_error(pop_scores(x, group_summary = "max"), "group_summary")
  expect_error(pop_scores(x, min_events = 500), "cannot exceed")
  expect_error(
    pop_scores(x, population = list(labels = x$labels)),
    "label_source"
  )
  expect_error(
    pop_scores(
      x,
      population = list(
        labels = lapply(x$labels, \(l) l[-1]),
        label_source = "gate"
      )
    ),
    "one label per event"
  )
})

test_that("a population-specific shift ranks first inside the population only", {
  x <- make_pop()
  res <- pop_scores(x)
  re <- res[res$check == "rest_emd" & res$marker == "A" & res$batch == "b2", ]
  p1 <- re$value[re$population %in% "p1"]
  p2 <- re$value[re$population %in% "p2"]
  expect_gt(p2, 1.2)
  expect_lt(p1, 0.4)
  # The shift is diluted in the all-event row.
  allv <- re$value[is.na(re$population)]
  expect_lt(allv, p2)
  # Marker B has no shift in either population.
  rb <- res[
    res$check == "rest_emd" &
      res$marker == "B" &
      res$batch == "b2" &
      !is.na(res$population),
  ]
  expect_true(all(rb$value < 0.4))
})

test_that("thin file x population cells give not_assessed rows with NA value", {
  x <- make_pop()
  # Make f1 almost free of p2.
  x$labels[["f1"]] <- c(rep("p1", 395), rep("p2", 5))
  res <- pop_scores(x, min_events = 50)
  thin <- res[res$check == "cell_events", ]
  expect_equal(thin$file, "f1")
  expect_equal(thin$population, "p2")
  expect_equal(thin$n_events, 5L)
  expect_true(is.na(thin$value))
  expect_equal(thin$status, "not_assessed")
  # A population with fewer than two batches left is not assessed as a whole.
  y <- make_pop()
  for (f in names(y$labels)[y$batch != "b1"]) {
    y$labels[[f]] <- rep("p1", 400)
  }
  r2 <- pop_scores(y)
  expect_equal(
    r2$status[r2$check == "population_scores" & r2$population == "p2"],
    "not_assessed"
  )
  expect_false("rest_emd" %in% r2$check[r2$population %in% "p2"])
})

test_that("share_rest_z finds the batch with the rare population", {
  x <- make_pop()
  res <- pop_scores(x)
  z <- res[res$check == "share_rest_z" & res$population == "p2", ]
  expect_equal(z$batch[which.max(abs(z$value))], "b3")
  sh <- res[res$check == "pop_share" & res$population == "p2", ]
  expect_equal(
    sh$value[sh$batch == "b3"],
    0.15,
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
  expect_equal(
    sh$value[sh$batch == "b1"],
    0.4,
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

test_that("batch_order_rho follows a planted drift across batches", {
  x <- make_pop()
  # Put a rising offset on marker B by batch (batch order = file index order).
  off <- c(b1 = 0, b2 = 0.5, b3 = 1)
  x$data <- Map(
    function(m, f) {
      m[, "B"] <- m[, "B"] + off[[x$batch[[f]]]]
      m
    },
    x$data,
    names(x$data)
  )
  res <- pop_scores(x, order = x$order)
  r <- res[res$check == "batch_order_rho" & is.na(res$population), ]
  expect_equal(r$value[r$marker == "B"], 1, ignore_attr = TRUE)
  expect_true(all(is.na(r$batch)))
  # No order: no batch_order_rho rows at all.
  expect_false("batch_order_rho" %in% pop_scores(x)$check)
  # Fewer than 3 batches gives NA.
  keep <- names(x$batch)[x$batch != "b3"]
  y <- list(data = x$data[keep], batch = x$batch[keep], labels = x$labels[keep])
  r2 <- pop_scores(y, order = x$order[keep])
  expect_true(all(is.na(r2$value[r2$check == "batch_order_rho"])))
})

test_that("cluster labels cap eligibility and the label source is recorded", {
  x <- make_pop()
  g <- pop_scores(x)
  cl <- pop_scores(
    x,
    population = list(labels = x$labels, label_source = "cluster")
  )
  expect_true(all(is.na(g$eligibility_cap)))
  expect_equal(
    unique(cl$eligibility_cap[!is.na(cl$population)]),
    "eligible_for_review"
  )
  expect_true(all(is.na(cl$eligibility_cap[is.na(cl$population)])))
  expect_equal(unique(cl$label_source[!is.na(cl$population)]), "cluster")
  expect_equal(attr(cl, "label_source"), "cluster")
})

test_that("population draws reproduce with a seed and restore the RNG state", {
  x <- make_pop()
  set.seed(99)
  before <- .Random.seed
  a <- pop_scores(x, seed = 5)
  expect_identical(.Random.seed, before)
  expect_equal(pop_scores(x, seed = 5), a)
  expect_false(isTRUE(all.equal(pop_scores(x, seed = 6)$value, a$value)))
})

test_that("reference groups combine population rows with group_summary", {
  x <- make_pop()
  # Group = pairs of files: each group has files from one batch only except these.
  grp <- stats::setNames(rep(paste0("d", 1:5), 3), names(x$data))
  m <- pop_scores(x, reference_group = grp, group_summary = "median")
  a <- pop_scores(x, reference_group = grp, group_summary = "mean")
  expect_equal(attr(a, "group_summary"), "mean")
  k <- m$check == "rest_emd" & m$population %in% "p2" & m$marker == "A"
  expect_gt(max(m$value[k]), 1)
  expect_false(isTRUE(all.equal(m$value, a$value)))
})

test_that("sw_batch_interaction separates a population-specific shift from a common one", {
  x <- make_pop()
  res <- pop_scores(x)
  it <- sw_batch_interaction(res, which_check = "rest_emd", min_populations = 2)
  expect_named(
    it,
    c(
      "marker",
      "batch",
      "check",
      "shift_pooled",
      "shift_sd_across_populations",
      "ratio",
      "n_populations",
      "status"
    )
  )
  ra <- it$ratio[it$marker == "A" & it$batch == "b2"]
  rb <- it$ratio[it$marker == "B" & it$batch == "b2"]
  expect_gt(ra, 0.5)
  expect_gt(ra, rb)
  # Too few populations: not_assessed with NA.
  few <- sw_batch_interaction(
    res,
    which_check = "rest_emd",
    min_populations = 3
  )
  expect_true(all(few$status == "not_assessed"))
  expect_true(all(is.na(few$ratio)))
  expect_error(sw_batch_interaction(res), "no default")
  expect_error(sw_batch_interaction(scores(x), "rest_emd", 2), "population")
  expect_error(sw_batch_interaction(res, "rest_emd", 1), "at least 2")
})
