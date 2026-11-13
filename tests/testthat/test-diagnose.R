# tests/testthat/test-diagnose.R
# Tests for R/diagnose.R: the stacked table and the per-check ranking.

make_dx <- function(seed = 1, n = 300) {
  set.seed(seed)
  ids <- paste0("f", 1:15)
  batch <- stats::setNames(rep(c("b1", "b2", "b3"), each = 5), ids)
  d <- lapply(ids, function(i) {
    cbind(
      Time = sort(stats::runif(n, 0, 100)),
      A = stats::rnorm(n) + if (batch[[i]] == "b2") 1 else 0,
      B = stats::rnorm(n),
      C = stats::rnorm(n)
    )
  })
  names(d) <- ids
  d$f3[1:90, "A"] <- 6 # saturated file
  list(data = d, batch = batch)
}

signal_set <- list(
  channels = c("A", "B", "C"),
  range_max = 6,
  range_min = NULL,
  time_channel = "Time",
  time_bin = 5,
  area_height = NULL,
  area_height_quantile = NULL,
  reference_group = NULL,
  reference = NULL,
  summary = NULL
)

test_that("every argument must be given", {
  expect_error(sw_diagnose(), "no default")
  x <- make_dx()
  expect_error(
    sw_diagnose(
      x$data,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
    "at least one module"
  )
  expect_error(
    sw_diagnose(
      NULL,
      NULL,
      signal_set,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
    "`data` is needed"
  )
  expect_error(
    sw_diagnose(
      x$data,
      "paths",
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
    "must be NULL or a list"
  )
})

test_that("batch$population is an error, not silently dropped", {
  x <- make_dx()
  expect_error(
    sw_diagnose(
      x$data,
      header = NULL,
      signal = NULL,
      anomaly = NULL,
      batch = list(
        batch = x$batch,
        markers = c("A", "B", "C"),
        n_events = 200,
        seed = 1,
        order = NULL,
        reference_group = NULL,
        population = list(labels = list(), label_source = "gate")
      ),
      cluster = NULL,
      relativise = NULL,
      otd = NULL,
      negative = NULL,
      gate = NULL,
      reference = NULL,
      batch_pop = NULL,
      flag_iqr = NULL,
      batch_otd = NULL,
      batch_mixing = NULL
    ),
    "population"
  )
})

test_that("signal and batch modules stack into one table", {
  x <- make_dx()
  dx <- sw_diagnose(
    x$data,
    header = NULL,
    signal = signal_set,
    anomaly = NULL,
    batch = list(
      batch = x$batch,
      markers = c("A", "B", "C"),
      n_events = 200,
      seed = 1,
      order = NULL,
      reference_group = NULL
    ),
    cluster = NULL,
    relativise = NULL,
    otd = NULL,
    negative = NULL,
    gate = NULL,
    reference = NULL,
    batch_pop = NULL,
    flag_iqr = NULL,
    batch_otd = NULL,
    batch_mixing = NULL
  )
  expect_s3_class(dx, "sw_diagnosis")
  expect_equal(attr(dx, "modules"), c("signal", "batch"))
  expect_setequal(
    names(dx),
    c("module", "subject", "check", "item", "value", "text")
  )
  expect_setequal(unique(dx$module), c("signal", "batch"))
  sat <- dx[dx$check == "saturation_share" & dx$item == "A", ]
  expect_equal(sat$subject[which.max(sat$value)], "f3")
  eta <- dx[dx$check == "eta2", ]
  expect_equal(eta$item[which.max(eta$value)], "A")
})

test_that("the anomaly module keeps forest scores, and novelty when asked", {
  skip_if_not_installed("isotree")
  skip_if_not_installed("mclust")
  x <- make_dx()
  an <- list(
    channels = c("A", "B", "C"),
    features = "quantile",
    quantile_step = 0.25,
    n_events = 200,
    seed = 1,
    scale = "test",
    reference_group = NULL,
    ntrees = 100,
    sample_size = 10,
    ndim = 1
  )
  dx <- sw_diagnose(
    x$data,
    NULL,
    NULL,
    an,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
  )
  expect_equal(unique(dx$check), "outlier_score")
  expect_equal(nrow(dx), 15)
  an$novelty_reference <- paste0("f", c(1:2, 4:8, 10:15)) # 13 files > 9 features
  an$novelty_quantile <- 0.05
  dx2 <- sw_diagnose(
    x$data,
    NULL,
    NULL,
    an,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
  )
  expect_setequal(unique(dx2$check), c("outlier_score", "novelty_log_density"))
  # With EMD features, one emd_z row per file and channel; a split marker tops it.
  d3 <- x$data
  d3$f7[, "B"] <- d3$f7[, "B"] + sample(c(-2, 2), nrow(d3$f7), replace = TRUE)
  an_emd <- an[c(
    "channels",
    "n_events",
    "seed",
    "scale",
    "ntrees",
    "sample_size",
    "ndim"
  )]
  an_emd <- c(
    an_emd,
    list(
      features = "emd",
      quantile_step = NULL,
      aggregate_from = "all",
      agg_size = 1500,
      reference = NULL,
      reference_group = NULL
    )
  )
  dx3 <- sw_diagnose(
    d3,
    NULL,
    NULL,
    an_emd,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
  )
  ez <- dx3[dx3$check == "emd_z", ]
  expect_equal(nrow(ez), 15 * 3)
  expect_equal(ez$subject[which.max(ez$value)], "f7")
  expect_equal(ez$item[which.max(ez$value)], "B")
})

test_that("the header module reads FCS headers", {
  skip_if_not_installed("flowCore")
  d <- withr::local_tempdir()
  mk <- function(name, cyt) {
    ff <- flowCore::flowFrame(matrix(
      stats::rnorm(40),
      ncol = 2,
      dimnames = list(NULL, c("A", "B"))
    ))
    flowCore::keyword(ff) <- list("$CYT" = cyt)
    p <- file.path(d, name)
    flowCore::write.FCS(ff, p)
    p
  }
  paths <- c(mk("r.fcs", "X"), mk("a.fcs", "X"), mk("b.fcs", "Y"))
  dx <- sw_diagnose(
    NULL,
    header = list(paths = paths, keywords = "$CYT", reference = paths[1]),
    signal = NULL,
    anomaly = NULL,
    batch = NULL,
    cluster = NULL,
    relativise = NULL,
    otd = NULL,
    negative = NULL,
    gate = NULL,
    reference = NULL,
    batch_pop = NULL,
    flag_iqr = NULL,
    batch_otd = NULL,
    batch_mixing = NULL
  )
  hit <- dx[dx$subject == "b.fcs" & dx$check == "keyword $CYT", ]
  expect_equal(hit$value, 1)
  expect_match(hit$text, "differs: Y vs X")
  expect_equal(dx$value[dx$subject == "a.fcs" & dx$check == "keyword $CYT"], 0)
})

test_that("sw_diagnosis_top ranks per check and needs n", {
  x <- make_dx()
  dx <- sw_diagnose(
    x$data,
    NULL,
    signal_set,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
  )
  expect_error(sw_diagnosis_top(dx), "n")
  expect_error(sw_diagnosis_top(data.frame(), n = 2), "must come from")
  top <- sw_diagnosis_top(dx, n = 2)
  expect_true(all(top$rank <= 2))
  sat <- top[top$check == "saturation_share" & top$item == "A", ]
  expect_equal(nrow(sat), 2)
  expect_equal(sat$subject[sat$rank == 1], "f3")
})

test_that("signed checks rank by size, so a large negative z comes first", {
  x <- make_dx()
  # Lower marker B of one file by 3 SD: its median_z is large and negative.
  x$data$f7[, "B"] <- x$data$f7[, "B"] - 3
  dx <- sw_diagnose(
    x$data,
    NULL,
    signal_set,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
  )
  z <- dx[dx$check == "median_z" & dx$item == "B", ]
  expect_lt(z$value[z$subject == "f7"], -3)
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(top$subject[top$check == "median_z" & top$item == "B"], "f7")
})

test_that("the cluster module ranks each marker and cluster as its own item", {
  x <- make_dx()
  dx <- sw_diagnose(
    x$data,
    header = NULL,
    signal = NULL,
    anomaly = NULL,
    batch = NULL,
    cluster = list(
      markers = c("A", "B", "C"),
      k = 3,
      n_events = 200,
      seed = 1,
      nstart = 3,
      min_cluster_events = 10,
      reference = NULL,
      reference_group = NULL
    ),
    relativise = NULL,
    otd = NULL,
    negative = NULL,
    gate = NULL,
    reference = NULL,
    batch_pop = NULL,
    flag_iqr = NULL,
    batch_otd = NULL,
    batch_mixing = NULL
  )
  expect_equal(attr(dx, "modules"), "cluster")
  expect_setequal(unique(dx$check), c("freq_z", "median_z"))
  expect_true(all(grepl("^n_events=", dx$text)))
  expect_setequal(unique(dx$item[dx$check == "freq_z"]), paste("cluster", 1:3))
  expect_setequal(
    unique(dx$item[dx$check == "median_z"]),
    paste(rep(c("A", "B", "C"), 3), "/ cluster", rep(1:3, each = 3))
  )
  # Each marker and cluster is ranked on its own: n = 1 gives one row per item.
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(sum(top$check == "median_z"), 9)
})

test_that("ranking keys: signed, low-first and ratio checks, without spurious warnings", {
  key <- .sw_diag_key(
    c(
      "neg_corr",
      "separation",
      "gain_batch_ratio",
      "offset",
      "novelty_log_density"
    ),
    c(-0.9, 0.5, 0.5, -3, -10)
  )
  expect_equal(key, c(0.9, -0.5, abs(log(0.5)), -3, 10))
  expect_no_warning(.sw_diag_key(c("offset", "gain_batch_ratio"), c(-5, 2)))
})

test_that("spread_ratio ranks by abs(log(value)): 0.5 and 2 are equal (DR39)", {
  key <- .sw_diag_key(rep("spread_ratio", 3), c(0.5, 2, 1))
  expect_equal(key[[1]], key[[2]])
  expect_equal(key[[3]], 0)
  expect_gt(key[[1]], key[[3]])
  expect_no_warning(.sw_diag_key(c("offset", "spread_ratio"), c(-5, 2)))
})

test_that("bimodality_z ranks by size; spread_z and rest_emd_z stay one-sided (G38)", {
  # A split blurs two populations, so the coefficient drops: a large negative z is the fault.
  key <- .sw_diag_key(
    c("bimodality_z", "spread_z", "rest_emd_z"),
    c(-4, -4, -4)
  )
  expect_equal(key, c(4, -4, -4))
})

test_that("the gate module stacks sw_gate_qc rows with gate/channel items", {
  set.seed(3)
  files <- lapply(c(0, 0, 0, 0, 1.5), function(shift) {
    cbind(A = stats::rnorm(500, 5 + shift), B = stats::rnorm(500, 5))
  })
  names(files) <- paste0("f", 1:5)
  g <- sw_gate_range("pos", "A", 4, Inf, parent = NA)
  dx <- sw_diagnose(
    files,
    header = NULL,
    signal = NULL,
    anomaly = NULL,
    batch = NULL,
    cluster = NULL,
    relativise = NULL,
    otd = NULL,
    negative = NULL,
    gate = list(gates = list(g), reference = "f1", edge_width = 0.5),
    reference = NULL,
    batch_pop = NULL,
    flag_iqr = NULL,
    batch_otd = NULL,
    batch_mixing = NULL
  )
  expect_equal(attr(dx, "modules"), "gate")
  expect_true(all(dx$module == "gate"))
  expect_true("pos/A" %in% dx$item[dx$check == "centroid_shift"])
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(top$subject[top$check == "kept_share_z"], "f5")
})

# DR33 modules: reference, batch_pop, flag_iqr, batch_otd, batch_mixing -------------------

# Every module NULL unless set, so a test names only the modules it uses.
dx_call <- function(data, ...) {
  args <- list(
    data = data,
    header = NULL,
    signal = NULL,
    anomaly = NULL,
    batch = NULL,
    cluster = NULL,
    relativise = NULL,
    otd = NULL,
    negative = NULL,
    gate = NULL,
    reference = NULL,
    batch_pop = NULL,
    flag_iqr = NULL,
    batch_otd = NULL,
    batch_mixing = NULL
  )
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(sw_diagnose, args)
}

# ref1 sits in five batches and is shifted in marker A in b2. ref2 sits in two
# batches with only 60 events per file, below the min_events of the test.
make_dx_ref <- function(shift = 2, n = 400, seed = 31) {
  set.seed(seed)
  rows <- data.frame(
    id = c(paste0("r1_b", 1:5), "r2_b1", "r2_b2"),
    batch = c(paste0("b", 1:5), "b1", "b2"),
    ref = c(rep("ref1", 5), "ref2", "ref2")
  )
  d <- lapply(seq_len(nrow(rows)), function(i) {
    m <- if (rows$ref[i] == "ref2") 60 else n
    cbind(
      A = stats::rnorm(m) + if (rows$id[i] == "r1_b2") shift else 0,
      B = stats::rnorm(m)
    )
  })
  names(d) <- rows$id
  man <- data.frame(
    sample_id = rows$id,
    batch = rows$batch,
    reference_role = "technical_replicate",
    reference_id = rows$ref
  )
  list(d = d, man = man)
}

ref_set <- function(x) {
  list(
    manifest = x$man,
    markers = c("A", "B"),
    n_events = 300,
    seed = 1,
    order = NULL,
    min_batches_per_reference = 3,
    min_events = 100,
    roles = "technical_replicate",
    transform = "asinh",
    cofactor = 6000
  )
}

test_that("every new module is a required argument", {
  x <- make_dx()
  expect_error(
    sw_diagnose(
      x$data,
      NULL,
      signal_set,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL,
      NULL
    ),
    "'reference', 'batch_pop', 'flag_iqr', 'batch_otd', 'batch_mixing'"
  )
  expect_error(
    dx_call(x$data, reference = "manifest"),
    "must be NULL or a list"
  )
  expect_error(dx_call(NULL, reference = list()), "`data` is needed")
  expect_error(dx_call(NULL, batch_otd = list()), "`data` is needed")
  expect_error(dx_call(NULL, batch_mixing = list()), "`data` is needed")
  expect_error(dx_call(NULL, batch_pop = list()), "`data` is needed")
})

test_that("the reference module finds the planted shifted batch and ranks per reference and marker", {
  x <- make_dx_ref()
  dx <- dx_call(x$d, reference = ref_set(x))
  expect_equal(attr(dx, "modules"), "reference")
  expect_true(all(dx$module == "reference"))
  expect_setequal(
    names(dx),
    c("module", "subject", "check", "item", "value", "text")
  )
  z <- dx[dx$check == "rest_emd_z" & dx$item == "ref1 / A", ]
  expect_equal(z$subject[which.max(z$value)], "b2")
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(
    top$subject[top$check == "rest_emd_z" & top$item == "ref1 / A"],
    "b2"
  )
  # count_ratio has no marker, so its item is the reference alone.
  expect_true("ref1" %in% dx$item[dx$check == "count_ratio"])
  # Same numbers as the module itself: the wrapper only restacks them.
  direct <- do.call(sw_ref_consistency, c(list(data = x$d), ref_set(x)))
  expect_equal(
    sort(dx$value[dx$check == "pair_emd"]),
    sort(direct$value[direct$check == "pair_emd"])
  )
})

test_that("a thin reference gives not_assessed rows with NA value, and nothing is dropped", {
  x <- make_dx_ref()
  dx <- dx_call(x$d, reference = ref_set(x))
  thin <- dx[grepl("^ref2", dx$item) & dx$check != "count_ratio", ]
  expect_gt(nrow(thin), 0)
  expect_true(all(is.na(thin$value)))
  expect_true(all(grepl("^not_assessed", thin$text)))
  # Every NA row says why; none is dropped without a word.
  expect_true(all(grepl("^not_assessed", dx$text[is.na(dx$value)])))
  # NA rows never enter a ranking: they follow the ranked rows with rank NA.
  top <- sw_diagnosis_top(dx, n = 5)
  thin_top <- top[grepl("^ref2", top$item) & top$check != "count_ratio", ]
  expect_true(all(is.na(thin_top$value)))
  expect_true(all(is.na(thin_top$rank)))
})

test_that("the reference module is seed reproducible and restores the RNG state", {
  x <- make_dx_ref()
  set.seed(77)
  before <- .Random.seed
  a <- dx_call(x$d, reference = ref_set(x))
  expect_identical(.Random.seed, before)
  b <- dx_call(x$d, reference = ref_set(x))
  expect_equal(a, b)
  s2 <- ref_set(x)
  s2$seed <- 2
  c2 <- dx_call(x$d, reference = s2)
  expect_false(isTRUE(all.equal(a$value, c2$value)))
})

test_that("a reference setting left out stops with the module's own message", {
  x <- make_dx_ref()
  s <- ref_set(x)
  s$min_events <- NULL
  expect_error(dx_call(x$d, reference = s), "min_events")
})

# Three populations per file. Marker A is shifted in b2 inside p2 only; p3 has 10
# events per file, below the min_events of the test.
make_dx_pop <- function(n = 400, seed = 2) {
  set.seed(seed)
  ids <- paste0("f", 1:15)
  b <- stats::setNames(rep(c("b1", "b2", "b3"), each = 5), ids)
  labs <- list()
  d <- lapply(ids, function(i) {
    cnt <- list(
      b1 = c(230, 160, 10),
      b2 = c(240, 150, 10),
      b3 = c(250, 140, 10)
    )[[b[[i]]]]
    lab <- rep(c("p1", "p2", "p3"), cnt)
    a <- stats::rnorm(n) + ifelse(lab == "p2" & b[[i]] == "b2", 2, 0)
    labs[[i]] <<- lab
    cbind(A = a, B = stats::rnorm(n))
  })
  names(d) <- ids
  list(data = d, batch = b, labels = labs)
}

pop_set <- function(x, ...) {
  args <- list(
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
  args
}

test_that("the batch_pop module ranks per marker and population and keeps thin cells visible", {
  x <- make_dx_pop()
  dx <- dx_call(x$data, batch_pop = pop_set(x))
  expect_equal(attr(dx, "modules"), "batch_pop")
  expect_true(all(dx$module == "batch_pop"))
  # The all-event rows belong to the batch module and are not stacked here.
  expect_false(any(dx$item %in% c("A", "B")))
  z <- dx[dx$check == "rest_emd_z" & dx$item == "A / p2", ]
  expect_equal(z$subject[which.max(z$value)], "b2")
  # The same marker in p1 shows no planted shift: the b2 z is far smaller.
  z1 <- dx[dx$check == "rest_emd_z" & dx$item == "A / p1", ]
  expect_gt(max(z$value), max(z1$value, na.rm = TRUE))
  # pop_share rows carry the population alone as item.
  expect_true("p2" %in% dx$item[dx$check == "pop_share"])
  # The rare population is thinner than min_events everywhere.
  rare <- dx[
    grepl("p3", dx$item) & !dx$check %in% c("pop_share", "share_rest_z"),
  ]
  expect_gt(nrow(rare), 0)
  expect_true(all(is.na(rare$value)))
  expect_true(all(grepl("^not_assessed", rare$text)))
  # Every NA row says why: thin cell, or an assessed row with no z.
  expect_true(all(grepl("^not_assessed|^no value", dx$text[is.na(dx$value)])))
  # With shares that differ by batch, share_rest_z of p2 is a real number.
  sz <- dx[dx$check == "share_rest_z" & dx$item == "p2", ]
  expect_false(anyNA(sz$value))
  # Seed reproducible, and the caller's RNG state is restored.
  set.seed(3)
  before <- .Random.seed
  a <- dx_call(x$data, batch_pop = pop_set(x))
  expect_identical(.Random.seed, before)
  expect_equal(a, dx)
})

test_that("batch_pop needs a population and a min_events, and rows match the direct call", {
  x <- make_dx_pop()
  expect_error(
    dx_call(x$data, batch_pop = pop_set(x, population = NULL)),
    "batch_pop\\$population"
  )
  s <- pop_set(x)
  s$min_events <- NULL
  expect_error(dx_call(x$data, batch_pop = s), "min_events")
  dx <- dx_call(x$data, batch_pop = pop_set(x))
  direct <- do.call(sw_batch_scores, c(list(data = x$data), pop_set(x)))
  direct <- direct[!is.na(direct$population), ]
  expect_equal(
    sum(!is.na(dx$value[dx$check == "pair_emd"])),
    sum(!is.na(direct$value[direct$check == "pair_emd"]))
  )
})

test_that("batch_order_rho is stacked when order is given and ranks by size", {
  x <- make_dx_pop()
  ord <- stats::setNames(seq_along(x$batch), names(x$batch))
  dx <- dx_call(x$data, batch_pop = pop_set(x, order = ord))
  expect_true("batch_order_rho" %in% dx$check)
  expect_equal(.sw_diag_key("batch_order_rho", -0.9), 0.9)
})

make_dx_cells <- function(seed = 5) {
  set.seed(seed)
  cells <- data.frame(
    batch = rep(paste0("b", 1:5), each = 400),
    control = "ctrl1",
    CD3 = c(stats::rnorm(1000, 0, 0.3), stats::rnorm(1000, 3, 0.3))[
      sample(2000)
    ]
  )
  shift <- cells$batch == "b5" & cells$CD3 > 1.5
  cells$CD3[shift] <- cells$CD3[shift] + 2
  cells
}

iqr_set <- function(cells, ...) {
  args <- list(
    cells = cells,
    markers = "CD3",
    thresholds = c(CD3 = 1.5),
    iqr_multiplier = 1.5,
    min_cells_per_side = 50,
    quantile_type = 7,
    min_effect = c(neg_median = 0.2, pos_median = 0.2, pct_positive = 5)
  )
  args[names(list(...))] <- list(...)
  args
}

test_that("the flag_iqr module needs no data, stacks the effect, and finds the planted batch", {
  cells <- make_dx_cells()
  dx <- dx_call(NULL, flag_iqr = iqr_set(cells))
  expect_equal(attr(dx, "modules"), "flag_iqr")
  expect_equal(unique(dx$check), "iqr_effect")
  expect_true("ctrl1 / pos_median / CD3" %in% dx$item)
  pos <- dx[dx$item == "ctrl1 / pos_median / CD3", ]
  expect_equal(pos$subject[which.max(pos$value)], "b5")
  # Only the batch outside the fence has a positive effect.
  expect_gt(pos$value[pos$subject == "b5"], 0.2)
  expect_equal(sum(pos$value > 0), 1)
  # The fence flag is not stacked; the status is.
  expect_false("flag" %in% names(dx))
  expect_true(all(grepl("^status=", dx$text)))
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(top$subject[top$item == "ctrl1 / pos_median / CD3"], "b5")
})

test_that("flag_iqr stops without min_effect, and thin panels stay as NA rows", {
  cells <- make_dx_cells()
  s <- iqr_set(cells)
  s$min_effect <- NULL
  expect_error(dx_call(NULL, flag_iqr = s), "min_effect")
  few <- iqr_set(cells, min_cells_per_side = 1e6)
  dx <- dx_call(NULL, flag_iqr = few)
  expect_true(all(is.na(dx$value)))
  expect_true(all(grepl("ineligible_min_cells", dx$text)))
})

make_dx_otd <- function(seed = 11) {
  set.seed(seed)
  shifts <- list(c(0, 0), c(0, 0), c(0, 0), c(0, 0), c(3, 0), c(3, 0))
  files <- lapply(shifts, \(s) {
    m <- matrix(stats::rnorm(600), ncol = 2, dimnames = list(NULL, c("A", "B")))
    sweep(m, 2, s, "+")
  })
  names(files) <- paste0("f", 1:6)
  list(
    data = files,
    batch = stats::setNames(c("b1", "b1", "b2", "b2", "b3", "b3"), names(files))
  )
}

otd_set <- function(x, mode) {
  list(
    batch = x$batch,
    markers = c("A", "B"),
    n_events = 120,
    blur = 1,
    seed = 3,
    max_iter = 300,
    tol = 1e-5,
    mode = mode
  )
}

test_that("the batch_otd module names the check after the mode and ranks the shifted batch first", {
  x <- make_dx_otd()
  dx <- dx_call(x$data, batch_otd = otd_set(x, "rest"))
  expect_equal(attr(dx, "modules"), "batch_otd")
  expect_equal(unique(dx$check), "otd_rest")
  expect_true(all(is.na(dx$item)))
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(top$subject, "b3")
  pw <- dx_call(x$data, batch_otd = otd_set(x, "pairwise"))
  expect_equal(unique(pw$check), "otd_pairwise")
  expect_equal(nrow(pw), 3)
  # Seed reproducible and the caller's RNG state is restored.
  set.seed(5)
  before <- .Random.seed
  a <- dx_call(x$data, batch_otd = otd_set(x, "rest"))
  expect_identical(.Random.seed, before)
  expect_equal(a, dx)
  s <- otd_set(x, "rest")
  s$tol <- NULL
  expect_error(dx_call(x$data, batch_otd = s), "tol")
})

make_dx_mix <- function(n = 300, seed = 21) {
  set.seed(seed)
  b <- rep(1:4, each = 2)
  files <- lapply(b, \(i) {
    m <- matrix(
      stats::rnorm(n * 2),
      ncol = 2,
      dimnames = list(NULL, c("A", "B"))
    )
    if (i == 4) {
      m[, 1] <- m[, 1] + 50
    }
    m
  })
  names(files) <- paste0("f", seq_along(files))
  list(data = files, batch = stats::setNames(paste0("b", b), names(files)))
}

mix_set <- function(x, ...) {
  args <- list(
    batch = x$batch,
    markers = c("A", "B"),
    k = 45,
    perplexity = 15,
    alpha = 0.05,
    n_events = 200,
    seed = 5,
    min_events = 100,
    population = NULL
  )
  args[names(list(...))] <- list(...)
  args
}

test_that("the batch_mixing module ranks the separated batch first on low LISI", {
  x <- make_dx_mix()
  dx <- dx_call(x$data, batch_mixing = mix_set(x))
  expect_equal(attr(dx, "modules"), "batch_mixing")
  expect_setequal(unique(dx$check), c("lisi_batch", "kbet_rejection"))
  expect_true("all batches" %in% dx$subject)
  lisi <- dx[dx$check == "lisi_batch", ]
  expect_equal(lisi$subject[which.min(lisi$value)], "b4")
  top <- sw_diagnosis_top(dx, n = 1)
  expect_equal(top$subject[top$check == "lisi_batch"], "b4")
  # Matches the module's own overall numbers.
  direct <- do.call(sw_batch_mixing, c(list(data = x$data), mix_set(x)))
  expect_equal(
    lisi$value[lisi$subject == "all batches"],
    direct$value[direct$check == "lisi_batch" & is.na(direct$batch)]
  )
  # Seed reproducible, and the caller's RNG state is restored.
  set.seed(9)
  before <- .Random.seed
  a <- dx_call(x$data, batch_mixing = mix_set(x))
  expect_identical(.Random.seed, before)
  expect_equal(a, dx)
  s <- mix_set(x)
  s$perplexity <- NULL
  expect_error(dx_call(x$data, batch_mixing = s), "perplexity")
})

test_that("batch_mixing units too thin for min_events stay as NA rows with their status", {
  x <- make_dx_mix()
  dx <- dx_call(
    x$data,
    batch_mixing = mix_set(x, min_events = 5000, n_events = 150)
  )
  expect_true(all(is.na(dx$value)))
  expect_true(all(grepl("status=not_assessed", dx$text)))
})

test_that("direction rules: rest_emd_z one-sided, LISI low-first, kBET high-first (DR33)", {
  # rest_emd_z: a negative z means closer to normal, so it ranks below a positive one.
  expect_equal(.sw_diag_key("rest_emd_z", c(-4, 4)), c(-4, 4))
  # LISI: the lowest value is the least mixed, so its key is the largest.
  k <- .sw_diag_key(rep("lisi_batch", 2), c(1.1, 2.9))
  expect_gt(k[[1]], k[[2]])
  # kBET rejection grows with a fault and ranks by raw value.
  expect_equal(.sw_diag_key(rep("kbet_rejection", 2), c(0.1, 0.9)), c(0.1, 0.9))
  # Signed and ratio checks of the new modules.
  expect_equal(.sw_diag_key("median_shift_z", -3), 3)
  expect_equal(.sw_diag_key("share_rest_z", -2), 2)
  r <- .sw_diag_key(rep("count_ratio", 2), c(0.5, 2))
  expect_equal(r[[1]], r[[2]])
  expect_no_warning(.sw_diag_key(c("offset", "count_ratio"), c(-5, 2)))
  # Raw checks: effect and OTD.
  expect_equal(.sw_diag_key(c("iqr_effect", "otd_rest"), c(0.5, 2)), c(0.5, 2))
})

test_that("modules from different families stack into one table with the shared columns", {
  x <- make_dx_otd()
  cells <- make_dx_cells()
  dx <- dx_call(
    x$data,
    batch_otd = otd_set(x, "rest"),
    flag_iqr = iqr_set(cells)
  )
  expect_equal(attr(dx, "modules"), c("flag_iqr", "batch_otd"))
  expect_setequal(unique(dx$module), c("flag_iqr", "batch_otd"))
})

test_that("sw_diagnosis_top keeps not-assessed rows, unranked, after the ranked ones", {
  x <- make_dx()
  # One file cut below summary$min_events: its summary_* values are NA.
  x$data$f4 <- x$data$f4[1:50, ]
  s <- signal_set
  s$summary <- list(probs = c(0.25, 0.75), min_events = 100, stat = "median")
  dx <- sw_diagnose(
    x$data,
    NULL,
    s,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
  )
  top <- sw_diagnosis_top(dx, n = 2)
  thin <- top[top$subject == "f4" & top$check == "summary_median", ]
  expect_gt(nrow(thin), 0)
  expect_true(all(is.na(thin$value)))
  expect_true(all(is.na(thin$rank)))
  expect_match(thin$text[[1]], "min_events|events")
  # Ranked rows come first, unranked rows last.
  expect_false(is.unsorted(is.na(top$rank)))
  # Every NA-value row of the diagnosis is listed.
  expect_equal(sum(is.na(top$value)), sum(is.na(dx$value)))
})
