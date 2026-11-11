# Tests for the spike-in module (R/unmix_spike.R, R/unmix_spike_helpers.R).
# Mocks only: an ordinary-least-squares unmixer on overlapping spectra, a
# coupled unmixer that shifts every cell by a file-level quantile, and a
# dropping unmixer that loses one dye.

# Fixtures ---------------------------------------------------------------------

det <- paste0("D", 1:6)
spike_spectra <- rbind(
  A = c(1, 0.6, 0.2, 0.1, 0, 0),
  B = c(0.1, 0.5, 1, 0.4, 0.1, 0),
  C = c(0, 0, 0.2, 0.6, 1, 0.3)
)
colnames(spike_spectra) <- det
key_cols <- c("Time", "FSC-A")

make_raw <- function(n = 120, seed = 1) {
  set.seed(seed)
  abund <- matrix(abs(stats::rnorm(n * 3, 500, 200)), n, 3)
  signal <- abund %*% spike_spectra + matrix(stats::rnorm(n * 6, 0, 5), n, 6)
  colnames(signal) <- det
  cbind(
    Time = seq_len(n) + seed / 1000,
    `FSC-A` = stats::runif(n, 1e4, 1e5),
    signal
  )
}

raw_mats <- list(s1 = make_raw(seed = 1), s2 = make_raw(seed = 2))
pops <- list(
  s1 = list(P1 = 1:50, P2 = 51:90),
  s2 = list(P1 = 1:50, P2 = 51:90)
)

build_mem <- function(...) {
  args <- list(
    raw = raw_mats,
    populations = pops,
    spectra = spike_spectra,
    dyes = c(mA = "A", mB = "B"),
    levels = c(0, 1, 2, 3),
    cofactor = 6000,
    spike_fraction = 0.75,
    seed = 42,
    key_cols = key_cols,
    output_dir = NULL
  )
  args[names(list(...))] <- list(...)
  do.call(sw_unmix_spike_build, args)
}

# OLS unmixer: abundances per dye, written to "<dye>-A", key columns kept.
unmix_ols <- function(x) {
  ab <- t(qr.solve(t(spike_spectra), t(x[, det])))
  colnames(ab) <- paste0(rownames(spike_spectra), "-A")
  cbind(x[, key_cols], ab)
}

# Coupled unmixer: OLS, then every cell's A channel shifted by the file median.
unmix_coupled <- function(x) {
  u <- unmix_ols(x)
  u[, "A-A"] <- u[, "A-A"] - stats::median(u[, "A-A"])
  u
}

# Dropping unmixer: dye A never comes back.
unmix_drop_a <- function(x) {
  u <- unmix_ols(x)
  u[, "A-A"] <- 0
  u
}

# Unmixed table from an in-memory design and a mock unmixer per setup.
unmixed_table <- function(design, unmixers) {
  rows <- list()
  for (setup in names(unmixers)) {
    for (s in names(design$n_events)) {
      rows[[length(rows) + 1]] <- tibble::tibble(
        setup = setup,
        sample = s,
        set = "base",
        data = list(unmixers[[setup]](design$raw[[s]]))
      )
      for (set in design$sets$set) {
        rows[[length(rows) + 1]] <- tibble::tibble(
          setup = setup,
          sample = s,
          set = set,
          data = list(unmixers[[setup]](design$data[[set]][[s]]))
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}

# Build -----------------------------------------------------------------------

test_that("build names every missing required argument", {
  expect_error(
    sw_unmix_spike_build(raw = raw_mats),
    "populations.*spectra.*dyes.*levels.*cofactor.*spike_fraction.*seed.*key_cols.*output_dir"
  )
})

test_that("levels are drawn per population with exact counts", {
  d <- build_mem()
  counts <- d$key[d$key$set == "mA", ] |>
    dplyr::count(sample, pop, level)
  p1 <- counts[counts$sample == "s1" & counts$pop == "P1", ]
  # round(0.75 * 50) = 38 spiked cells over levels 1, 2, 3; 12 unspiked
  expect_equal(p1$n[p1$level == 0], 12)
  expect_equal(sum(p1$n[p1$level > 0]), 38)
  expect_lte(max(p1$n[p1$level > 0]) - min(p1$n[p1$level > 0]), 1)
  expect_equal(
    d$fraction$fraction[
      d$fraction$set == "mA" &
        d$fraction$sample == "s1" &
        d$fraction$pop == "P1"
    ],
    38 / 50
  )
})

test_that("the same seed gives the same key and every marker spikes the same cells", {
  a <- build_mem()
  b <- build_mem()
  expect_identical(a$key, b$key)
  ka <- a$key[a$key$set == "mA", c("sample", "raw_row", "level")]
  kb <- a$key[a$key$set == "mB", c("sample", "raw_row", "level")]
  expect_identical(ka, kb)
})

test_that("build leaves the caller's random stream alone", {
  set.seed(7)
  before <- stats::runif(1)
  set.seed(7)
  build_mem()
  after <- stats::runif(1)
  expect_identical(before, after)
})

test_that("an assignment table replaces the draw", {
  a <- build_mem()
  assign_tbl <- a$key[a$key$set == "mA", c("sample", "raw_row", "pop", "level")]
  assign_tbl$level <- rev(assign_tbl$level)
  b <- build_mem(assignment = assign_tbl, spike_fraction = NULL, seed = NULL)
  expect_identical(b$key$level[b$key$set == "mA"], assign_tbl$level)
  expect_error(
    build_mem(assignment = assign_tbl),
    "set 'spike_fraction' and 'seed' to NULL"
  )
  expect_error(
    build_mem(
      assignment = assign_tbl[-1, ],
      spike_fraction = NULL,
      seed = NULL
    ),
    "match 'populations' exactly"
  )
})

test_that("build stops on overlap, duplicate keys, spectrum scale and bad levels", {
  overlap <- list(s1 = list(P1 = 1:50, P2 = 40:90), s2 = pops$s2)
  expect_error(build_mem(populations = overlap), "share 11 row")
  dup <- raw_mats
  dup$s1[2, key_cols] <- dup$s1[1, key_cols]
  expect_error(build_mem(raw = dup), "share a key")
  bad <- spike_spectra
  bad["A", ] <- bad["A", ] * 2
  expect_error(build_mem(spectra = bad), "max 1; not so for: A")
  expect_error(build_mem(levels = c(1, 2)), "contain 0")
  expect_error(
    build_mem(key_cols = c("Time", "D1")),
    "must not include spectrum detectors"
  )
  expect_error(
    build_mem(raw = as.list(c(s1 = "a.fcs", s2 = "b.fcs"))),
    "NULL only for matrix input"
  )
})

test_that("spiked values equal raw plus k times the dye spectrum", {
  d <- build_mem()
  key <- d$key[d$key$set == "mB" & d$key$sample == "s2" & d$key$k > 0, ]
  r <- key$raw_row[1]
  expect_equal(
    unname(d$data$mB$s2[r, det]),
    unname(raw_mats$s2[r, det] + key$k[1] * spike_spectra["B", ])
  )
  unspiked <- d$key$raw_row[
    d$key$set == "mB" & d$key$sample == "s2" & d$key$k == 0
  ]
  expect_equal(d$data$mB$s2[unspiked, ], raw_mats$s2[unspiked, ])
})

# Score -----------------------------------------------------------------------

test_that("OLS unmixing recovers every spike and the files are independent", {
  d <- build_mem()
  sc <- sw_unmix_spike_score(
    d,
    unmixed_table(d, list(ols = unmix_ols)),
    0.5,
    1e-6
  )
  expect_equal(sc$cells$recovery, rep(1, nrow(sc$cells)), tolerance = 1e-6)
  expect_true(all(sc$summary$share_lost == 0))
  expect_true(all(sc$file_status$file_status == "independent"))
  expect_equal(sc$summary$q50, sc$summary$median_recovery)
  expect_setequal(unique(sc$baseline_check$group), c("level0", "untouched"))
})

test_that("a coupled unmixer is flagged and its recovery is kept", {
  d <- build_mem()
  sc <- sw_unmix_spike_score(
    d,
    unmixed_table(d, list(cpl = unmix_coupled)),
    0.5,
    1e-6
  )
  fs <- sc$file_status
  expect_true(all(fs$file_status[fs$set == "mA"] == "coupled"))
  moved <- sc$baseline_check[
    sc$baseline_check$set == "mA" & sc$baseline_check$channel == "A-A",
  ]
  expect_true(all(moved$share_moved > 0))
  expect_false(anyNA(sc$cells$recovery))
  expect_true(all(sc$summary$file_status[sc$summary$set == "mA"] == "coupled"))
  expect_true(any(grepl("coupled", sc$notes)))
})

test_that("a dropped dye is lost in every spiked cell", {
  d <- build_mem()
  sc <- sw_unmix_spike_score(
    d,
    unmixed_table(d, list(drop = unmix_drop_a)),
    0.5,
    1e-6
  )
  expect_true(all(sc$summary$share_lost[sc$summary$set == "mA"] == 1))
  expect_true(all(sc$summary$share_lost[sc$summary$set == "mB"] == 0))
})

test_that("quantiles follow a known recovery vector", {
  d <- build_mem()
  sc <- sw_unmix_spike_score(
    d,
    unmixed_table(d, list(ols = unmix_ols)),
    0.5,
    1e-6
  )
  row <- sc$summary[1, ]
  rec <- sc$cells$recovery[
    sc$cells$setup == row$setup &
      sc$cells$set == row$set &
      sc$cells$sample == row$sample &
      sc$cells$pop == row$pop &
      sc$cells$level == row$level
  ]
  expect_equal(row$q05, unname(stats::quantile(rec, 0.05)))
  expect_equal(row$q95, unname(stats::quantile(rec, 0.95)))
})

test_that("score stops on missing files and on misaligned rows", {
  d <- build_mem()
  um <- unmixed_table(d, list(ols = unmix_ols))
  expect_error(sw_unmix_spike_score(d, um[-1, ], 0.5, 1e-6), "lacks 1 file")
  shuffled <- um
  shuffled$data[[2]] <- shuffled$data[[2]][
    rev(seq_len(nrow(shuffled$data[[2]]))),
  ]
  expect_error(
    sw_unmix_spike_score(d, shuffled, 0.5, 1e-6),
    "key columns differ"
  )
  expect_error(sw_unmix_spike_score(d, um), "loss_cut.*baseline_tol")
})

test_that("results give the spike cosine and an AF-passes note", {
  d <- build_mem()
  res <- structure(
    list(
      spectra = spike_spectra,
      provenance = list(args = list(n.af.passes = 5))
    ),
    class = "sw_unmix_result"
  )
  sc <- sw_unmix_spike_score(
    d,
    unmixed_table(d, list(ols = unmix_ols)),
    0.5,
    1e-6,
    results = list(ols = res)
  )
  expect_equal(sc$spike_cosine$cosine, c(1, 1))
  expect_true(any(grepl("refine.af.quantile", sc$notes)))
})

# Co-spike --------------------------------------------------------------------

test_that("co-spike spikes only level1 cells and recovers both dyes", {
  d <- build_mem()
  cs <- sw_unmix_spike_cospike(
    d,
    marker = "mA",
    dye2 = c(mC = "C"),
    levels2 = c(0, 1),
    level1 = 1,
    output_dir = NULL
  )
  expect_equal(cs$sets$set, c("cospike_mC_0", "cospike_mC_1"))
  # cells the single design gave levels 2 and 3 stay unspiked
  orig <- d$key[d$key$set == "mA", ]
  key1 <- cs$key[cs$key$set == "cospike_mC_1", ]
  expect_true(all(key1$k[orig$level %in% c(2, 3)] == 0))
  expect_true(all(key1$k2[orig$level == 1] > 0))
  expect_lt(max(cs$fraction$file_fraction), max(d$fraction$file_fraction))
  sc <- sw_unmix_spike_score(
    cs,
    unmixed_table(cs, list(ols = unmix_ols)),
    0.5,
    1e-6
  )
  expect_equal(sc$cells$recovery, rep(1, nrow(sc$cells)), tolerance = 1e-6)
  rec2 <- sc$cells$recovery_2[sc$cells$set == "cospike_mC_1"]
  expect_equal(rec2, rep(1, length(rec2)), tolerance = 1e-6)
  expect_true(all(is.na(sc$cells$recovery_2[sc$cells$set == "cospike_mC_0"])))
  expect_error(
    sw_unmix_spike_cospike(
      d,
      "mA",
      c(mC = "C"),
      1,
      level1 = 0,
      output_dir = NULL
    ),
    "non-zero level"
  )
})

# FCS round trip -----------------------------------------------------------------

test_that("FCS files carry the spike keywords and a raised range", {
  skip_if_not_installed("flowCore")
  tmp <- withr::local_tempdir()
  raw_paths <- vapply(
    names(raw_mats),
    function(s) {
      path <- file.path(tmp, paste0(s, ".fcs"))
      flowCore::write.FCS(flowCore::flowFrame(raw_mats[[s]]), path)
      path
    },
    character(1)
  )
  out <- file.path(tmp, "spiked")
  d <- build_mem(raw = as.list(raw_paths), output_dir = out)
  expect_true(all(file.exists(d$files$path)))
  expect_true(file.exists(file.path(out, "spike_key.csv")))
  expect_setequal(unique(d$files$set), c("mA", "mB", "base"))
  ff <- flowCore::read.FCS(
    d$files$path[d$files$set == "mA"][1],
    transformation = FALSE,
    truncate_max_range = FALSE
  )
  kw <- flowCore::keyword(ff)
  expect_equal(kw[["$SW_SPIKE_DYE"]], "A")
  expect_equal(kw[["$SW_SPIKE_LEVELS"]], "0,1,2,3")
  x <- flowCore::exprs(ff)
  ranges <- as.numeric(unlist(kw[paste0("$P", seq_len(ncol(x)), "R")]))
  expect_true(all(apply(x, 2, max) < ranges))
  # the round trip is float32: compare with a tolerance, never == 0
  key <- d$key[d$key$set == "mA" & d$key$sample == "s1" & d$key$k > 0, ]
  expect_equal(
    unname(x[key$raw_row[1], det]),
    unname(raw_mats$s1[key$raw_row[1], det] + key$k[1] * spike_spectra["A", ]),
    tolerance = 1e-5
  )
  expect_error(
    build_mem(raw = as.list(raw_paths), output_dir = out),
    "already exists"
  )
})

test_that("score reads unmixed FCS files listed by pairs", {
  skip_if_not_installed("flowCore")
  tmp <- withr::local_tempdir()
  raw_paths <- vapply(
    names(raw_mats),
    function(s) {
      path <- file.path(tmp, paste0(s, ".fcs"))
      flowCore::write.FCS(flowCore::flowFrame(raw_mats[[s]]), path)
      path
    },
    character(1)
  )
  d <- build_mem(
    raw = as.list(raw_paths),
    output_dir = file.path(tmp, "spiked")
  )
  jobs <- sw_unmix_spike_jobs(d, c(ols = file.path(tmp, "unmixed")))
  expect_equal(nrow(jobs), 3)
  # the user's unmixing step, here the OLS mock on each job's folder
  for (i in seq_len(nrow(jobs))) {
    dir.create(jobs$out_dir[i], recursive = TRUE)
    for (f in list.files(
      jobs$raw_dir[i],
      pattern = "[.]fcs$",
      full.names = TRUE
    )) {
      x <- flowCore::exprs(flowCore::read.FCS(
        f,
        transformation = FALSE,
        truncate_max_range = FALSE
      ))
      out <- file.path(
        jobs$out_dir[i],
        paste0(tools::file_path_sans_ext(basename(f)), " AutoSpectral.fcs")
      )
      flowCore::write.FCS(flowCore::flowFrame(unmix_ols(x)), out)
    }
  }
  pairs <- sw_unmix_spike_pairs(d, c(ols = file.path(tmp, "unmixed")))
  expect_true(all(pairs$exists))
  sc <- sw_unmix_spike_score(d, pairs, loss_cut = 0.5, baseline_tol = 1e-3)
  expect_equal(sc$cells$recovery, rep(1, nrow(sc$cells)), tolerance = 1e-3)
  expect_true(all(sc$file_status$file_status == "independent"))
})

# Helpers ----------------------------------------------------------------------

test_that("spike_rows maps gated cleaned rows to raw rows by key", {
  cleaned <- lapply(raw_mats, function(x) x[seq(2, nrow(x), by = 2), ])
  gs <- list(
    s1 = list(`/T` = c(TRUE, FALSE, rep(TRUE, 58))),
    s2 = list(`/T` = 1:10)
  )
  rows <- sw_unmix_spike_rows(gs, c(T = "/T"), cleaned, raw_mats, key_cols)
  expect_equal(rows$s1$T[1:2], c(2L, 6L))
  expect_equal(rows$s2$T, seq(2L, 20L, by = 2L))
  dup <- cleaned
  dup$s1[2, key_cols] <- dup$s1[1, key_cols]
  expect_error(
    sw_unmix_spike_rows(gs, c(T = "/T"), dup, raw_mats, key_cols),
    "share a key"
  )
  odd <- cleaned
  odd$s1[1, "Time"] <- -1
  expect_error(
    sw_unmix_spike_rows(gs, c(T = "/T"), odd, raw_mats, key_cols),
    "no raw row"
  )
})

test_that("pairs and jobs refuse in-memory designs", {
  d <- build_mem()
  expect_error(sw_unmix_spike_pairs(d, c(a = "x")), "built from FCS files")
  expect_error(sw_unmix_spike_jobs(d, c(a = "x")), "no single folder")
})
