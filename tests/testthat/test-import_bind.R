# tests/testthat/test-import_bind.R
# Unit tests for R/import_bind.R. Writes tiny FCS files to a temp folder.

# One FCS file. `nextdata = TRUE` makes it a real two-data-set file: $NEXTDATA
# gets the byte offset of a second data set that is appended after the first.
# The longer value is paid for by cutting the $FILENAME value by the same
# number of bytes, so no other offset in the file moves.
bind_fcs <- function(dir, name, nextdata = FALSE) {
  skip_if_not_installed("flowCore")
  one <- function(path) {
    m <- matrix(
      stats::rnorm(20),
      ncol = 2,
      dimnames = list(NULL, c("B1-A", "V1-A"))
    )
    flowCore::write.FCS(flowCore::flowFrame(m), path)
    path
  }
  path <- file.path(dir, name)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  one(path)
  if (nextdata) {
    r <- readBin(path, "raw", file.size(path))
    find <- \(txt) {
      pat <- charToRaw(txt)
      starts <- which(r == pat[1])
      starts[vapply(
        starts,
        \(i) all(r[i:(i + length(pat) - 1)] == pat),
        logical(1)
      )][1]
    }
    val <- charToRaw(as.character(length(r)))
    cut <- length(val) - 1
    pos_n <- find("$NEXTDATA|0|")
    pos_f <- find("|FILENAME|") + 10
    end_f <- pos_f + which(r[pos_f:length(r)] == charToRaw("|"))[1] - 1
    second <- one(tempfile(fileext = ".fcs"))
    r <- c(
      r[seq_len(pos_n + 9)],
      val,
      r[(pos_n + 11):(end_f - cut - 1)],
      r[end_f:length(r)],
      readBin(second, "raw", file.size(second))
    )
    writeBin(r, path)
  }
  path
}

bind_manifest <- function(ids) {
  sw_sample_manifest(
    data.frame(id = ids, b = rep(c("b1", "b2"), length.out = length(ids))),
    sample_id = "id",
    batch = "b"
  )
}

# Call with every required argument; `...` overrides one.
bind_call <- function(manifest, paths, ...) {
  args <- list(
    manifest = manifest,
    paths = paths,
    key = "file_stem",
    case = "exact",
    recursive = FALSE,
    on_duplicate_file = "stop",
    on_extra_file = "list",
    nextdata = "bind"
  )
  over <- list(...)
  args[names(over)] <- over
  do.call(sw_manifest_bind_files, args)
}

test_that("every impactful argument is required and has no default", {
  man <- bind_manifest(c("a", "b"))
  expect_error(sw_manifest_bind_files(), "no default")
  full <- c(
    "manifest",
    "paths",
    "key",
    "case",
    "recursive",
    "on_duplicate_file",
    "on_extra_file",
    "nextdata"
  )
  for (arg in full) {
    args <- list(
      manifest = man,
      paths = "x",
      key = "file_stem",
      case = "exact",
      recursive = FALSE,
      on_duplicate_file = "stop",
      on_extra_file = "list",
      nextdata = "bind"
    )
    args[[arg]] <- NULL
    expect_error(
      do.call(sw_manifest_bind_files, args),
      paste0("'", arg, "'"),
      info = arg
    )
  }
})

test_that("no argument has a default", {
  no_default <- vapply(
    formals(sw_manifest_bind_files),
    \(x) identical(x, quote(expr = )),
    logical(1)
  )
  expect_true(all(no_default))
})

test_that("bad values stop before any listing, header read or package check", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  man <- bind_manifest(c("a", "b"))
  # A real folder, and listing/header reads that fail if they are reached.
  testthat::local_mocked_bindings(
    .sw_bind_list = \(...) stop("listed a folder"),
    .sw_bind_nextdata = \(...) stop("read a header")
  )
  expect_error(bind_call(man, d, key = "stem"), "must be one of")
  expect_error(bind_call(man, d, case = "lower"), "must be one of")
  expect_error(bind_call(man, d, recursive = NA), "TRUE or FALSE")
  expect_error(bind_call(man, d, nextdata = "ignore"), "must be one of")
  expect_error(bind_call(man, c(d, NA)), "non-empty character")
  expect_error(
    bind_call(man, c(d, "nope_nope")),
    "Path\\(s\\) not found: nope_nope"
  )
})

test_that("bad argument values stop with a clear message", {
  man <- bind_manifest(c("a", "b"))
  expect_error(bind_call(man, "x", key = "stem"), "key")
  expect_error(bind_call(man, "x", case = "lower"), "case")
  expect_error(bind_call(man, "x", recursive = NA), "recursive")
  expect_error(bind_call(man, "x", on_duplicate_file = "first"), "on_duplicate")
  expect_error(bind_call(man, "x", on_extra_file = "drop"), "on_extra_file")
  expect_error(bind_call(man, "x", nextdata = "ignore"), "nextdata")
  expect_error(
    bind_call(data.frame(sample_id = "a"), "x"),
    "sw_sample_manifest"
  )
  expect_error(bind_call(man, character()), "paths")
})

test_that("files bind by stem and orphans are listed both ways", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  bind_fcs(d, "b.fcs")
  bind_fcs(d, "extra.fcs")
  res <- bind_call(bind_manifest(c("a", "b", "gone")), d)
  expect_s3_class(res, "sw_manifest_bind_files")
  expect_named(
    res,
    c(
      "bound",
      "manifest_no_file",
      "file_no_manifest",
      "multi_dataset",
      "not_fcs",
      "provenance"
    )
  )
  expect_setequal(res$bound$sample_id, c("a", "b"))
  expect_true(all(file.exists(res$bound$path)))
  expect_true(all(c("batch", "path", "nextdata") %in% names(res$bound)))
  expect_equal(res$manifest_no_file$sample_id, "gone")
  expect_equal(res$manifest_no_file$reason, "manifest_row_without_file")
  expect_equal(basename(res$file_no_manifest$path), "extra.fcs")
  expect_equal(res$file_no_manifest$reason, "file_without_manifest_row")
  expect_equal(nrow(res$multi_dataset), 0L)
  # Every manifest row and every file lands in exactly one table.
  expect_equal(
    nrow(res$bound) + nrow(res$manifest_no_file),
    nrow(bind_manifest(c("a", "b", "gone")))
  )
  expect_equal(nrow(res$bound) + nrow(res$file_no_manifest), 3L)
})

test_that("key = file_name needs the extension in sample_id", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  res <- bind_call(bind_manifest(c("a", "a.fcs")), d, key = "file_name")
  expect_equal(res$bound$sample_id, "a.fcs")
  expect_equal(res$manifest_no_file$sample_id, "a")
})

test_that("case rule decides whether A.fcs matches a", {
  d <- withr::local_tempdir()
  bind_fcs(d, "A.fcs")
  man <- bind_manifest(c("a", "b"))
  exact <- bind_call(man, d, case = "exact")
  expect_equal(nrow(exact$bound), 0L)
  expect_equal(nrow(exact$file_no_manifest), 1L)
  ign <- bind_call(man, d, case = "ignore")
  expect_equal(ign$bound$sample_id, "a")
  expect_equal(basename(ign$bound$path), "A.fcs")
})

test_that("case = ignore stops when manifest ids collide after folding", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  expect_error(
    bind_call(bind_manifest(c("a", "A")), d, case = "ignore"),
    "A, a|a, A"
  )
})

test_that("recursion is a choice and subfolders stay out when FALSE", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  bind_fcs(d, "sub/b.fcs")
  man <- bind_manifest(c("a", "b"))
  flat <- bind_call(man, d, recursive = FALSE)
  expect_equal(flat$bound$sample_id, "a")
  expect_equal(flat$manifest_no_file$sample_id, "b")
  deep <- bind_call(man, d, recursive = TRUE)
  expect_setequal(deep$bound$sample_id, c("a", "b"))
})

test_that("duplicate file keys stop naming the key, or leave all copies unbound", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  bind_fcs(d, "x/a.fcs")
  bind_fcs(d, "b.fcs")
  man <- bind_manifest(c("a", "b"))
  expect_error(
    bind_call(man, d, recursive = TRUE, on_duplicate_file = "stop"),
    "a"
  )
  res <- bind_call(man, d, recursive = TRUE, on_duplicate_file = "unbind")
  expect_equal(res$bound$sample_id, "b")
  expect_equal(res$manifest_no_file$sample_id, "a")
  expect_equal(res$manifest_no_file$reason, "duplicate_file_key")
  dup <- res$file_no_manifest[
    res$file_no_manifest$reason == "duplicate_file_key",
  ]
  expect_equal(nrow(dup), 2L)
})

test_that("duplicates after case folding count as duplicates", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  bind_fcs(d, "x/A.fcs")
  expect_error(
    bind_call(
      bind_manifest(c("a", "b")),
      d,
      recursive = TRUE,
      case = "ignore"
    ),
    "a"
  )
})

test_that("a repeated path is one file, not a duplicate key", {
  d <- withr::local_tempdir()
  p <- bind_fcs(d, "a.fcs")
  res <- bind_call(bind_manifest(c("a", "b")), c(p, p))
  expect_equal(nrow(res$bound), 1L)
})

test_that("on_extra_file = stop names the extra files", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  bind_fcs(d, "extra.fcs")
  expect_error(
    bind_call(bind_manifest(c("a", "b")), d, on_extra_file = "stop"),
    "extra.fcs"
  )
})

test_that("non-FCS files in a folder are listed, not silently skipped", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  writeLines("x", file.path(d, "notes.txt"))
  res <- bind_call(bind_manifest(c("a", "b")), d)
  expect_equal(basename(res$not_fcs$path), "notes.txt")
  expect_equal(res$bound$sample_id, "a")
})

test_that("a path that does not exist stops", {
  d <- withr::local_tempdir()
  expect_error(
    bind_call(bind_manifest(c("a", "b")), file.path(d, "nope.fcs")),
    "nope.fcs"
  )
})

test_that("$NEXTDATA policy: bind, unbind or stop", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs", nextdata = TRUE)
  bind_fcs(d, "b.fcs")
  man <- bind_manifest(c("a", "b"))

  bind <- bind_call(man, d, nextdata = "bind")
  expect_setequal(bind$bound$sample_id, c("a", "b"))
  expect_equal(bind$multi_dataset$sample_id, "a")
  expect_true(bind$multi_dataset$bound)
  expect_gt(as.numeric(bind$multi_dataset$nextdata), 0)
  expect_equal(bind$bound$nextdata[bind$bound$sample_id == "b"], "0")

  unb <- bind_call(man, d, nextdata = "unbind")
  expect_equal(unb$bound$sample_id, "b")
  expect_equal(unb$manifest_no_file$sample_id, "a")
  expect_equal(unb$manifest_no_file$reason, "multi_dataset_unbound")
  expect_false(unb$multi_dataset$bound)

  expect_error(bind_call(man, d, nextdata = "stop"), "a.fcs")
})

test_that("an unreadable file is a not_assessed row, not a crash", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  writeLines("junk", file.path(d, "bad.fcs"))
  res <- bind_call(bind_manifest(c("a", "bad")), d)
  expect_setequal(res$bound$sample_id, c("a", "bad"))
  bad <- res$bound[res$bound$sample_id == "bad", ]
  expect_equal(bad$header_status, "not_assessed")
  expect_true(is.na(bad$nextdata))
  expect_equal(bad$reason, "not an FCS file")
  expect_equal(res$bound$header_status[res$bound$sample_id == "a"], "assessed")
  # It also shows in multi_dataset, so the NEXTDATA check cannot hide it.
  expect_equal(res$multi_dataset$status, "not_assessed")
  expect_equal(res$multi_dataset$sample_id, "bad")
  # nextdata = stop cannot judge it, so it does not stop on it.
  expect_no_error(bind_call(bind_manifest(c("a", "bad")), d, nextdata = "stop"))
})

test_that("files on disk are not changed", {
  d <- withr::local_tempdir()
  p <- bind_fcs(d, "a.fcs", nextdata = TRUE)
  before <- unname(tools::md5sum(p))
  bind_call(bind_manifest(c("a", "b")), d)
  expect_equal(unname(tools::md5sum(p)), before)
})

test_that("provenance records the settings and package versions", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  res <- bind_call(bind_manifest(c("a", "b")), d)
  prov <- res$provenance
  expect_equal(prov$args$key, "file_stem")
  expect_equal(prov$args$nextdata, "bind")
  expect_true(all(
    c("r_version", "platform", "datetime", "packages") %in% names(prov)
  ))
  expect_true("SpectraWeaveR" %in% prov$packages$package)
  expect_equal(prov$n_manifest_rows, 2L)
  expect_equal(prov$n_files, 1L)
  expect_true("seed" %in% names(prov))
  expect_null(prov$seed)
})

test_that("a manifest column that the result adds stops, naming it", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  for (col in c("path", "file", "nextdata", "header_status", "reason")) {
    df <- data.frame(id = c("a", "b"), b = c("b1", "b2"))
    df[[col]] <- c("mine1", "mine2")
    man <- sw_sample_manifest(df, sample_id = "id", batch = "b")
    expect_error(bind_call(man, d), col, info = col)
  }
})

test_that("hidden files and folders are listed, not invisible", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  bind_fcs(d, ".hid.fcs")
  bind_fcs(d, ".sub/c.fcs")
  writeLines("x", file.path(d, ".hidden.txt"))
  man <- bind_manifest(c("a", ".hid", "c"))
  flat <- bind_call(man, d, recursive = FALSE)
  expect_setequal(flat$bound$sample_id, c("a", ".hid"))
  expect_equal(
    basename(flat$not_fcs$path[flat$not_fcs$reason == "not_fcs_extension"]),
    ".hidden.txt"
  )
  expect_equal(
    basename(flat$not_fcs$path[flat$not_fcs$reason == "directory_not_listed"]),
    ".sub"
  )
  deep <- bind_call(man, d, recursive = TRUE)
  expect_setequal(deep$bound$sample_id, c("a", ".hid", "c"))
})

test_that("a directory named like an FCS file is listed with a reason", {
  d <- withr::local_tempdir()
  dir.create(file.path(d, "a.fcs"))
  res <- bind_call(bind_manifest(c("a", "b")), d, recursive = TRUE)
  expect_equal(res$not_fcs$reason, "directory_named_fcs")
  flat <- bind_call(bind_manifest(c("a", "b")), d, recursive = FALSE)
  expect_equal(flat$not_fcs$reason, "directory_not_listed")
  expect_equal(nrow(flat$bound), 0L)
})

test_that("an unreadable sub-folder is listed with a reason", {
  skip_on_os("windows")
  skip_if(identical(Sys.info()[["user"]], "root"), "root can read anything")
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  bind_fcs(d, "locked/b.fcs")
  Sys.chmod(file.path(d, "locked"), "000")
  withr::defer(Sys.chmod(file.path(d, "locked"), "755"))
  res <- bind_call(bind_manifest(c("a", "b")), d, recursive = TRUE)
  expect_equal(res$not_fcs$reason, "directory_unreadable")
  expect_equal(res$manifest_no_file$reason, "manifest_row_without_file")
})

test_that("a broken link is a not_fcs row, not a crash", {
  skip_on_os("windows")
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  file.symlink(file.path(d, "gone_target.fcs"), file.path(d, "b.fcs"))
  res <- bind_call(bind_manifest(c("a", "b")), d)
  expect_equal(res$bound$sample_id, "a")
  expect_equal(basename(res$not_fcs$path), "b.fcs")
  expect_equal(res$not_fcs$reason, "broken_link")
})

test_that("the key comes from the link name; links to one file each bind", {
  skip_on_os("windows")
  d <- withr::local_tempdir()
  real <- bind_fcs(d, "real/run123.fcs")
  dir.create(file.path(d, "links"))
  file.symlink(real, file.path(d, "links", "a.fcs"))
  res <- bind_call(bind_manifest(c("a", "b")), file.path(d, "links"))
  expect_equal(res$bound$sample_id, "a")
  expect_equal(res$bound$file, "a.fcs")
  file.symlink(real, file.path(d, "links", "z.fcs"))
  two <- bind_call(bind_manifest(c("a", "z")), file.path(d, "links"))
  expect_equal(two$provenance$n_files, 2L)
  expect_equal(two$bound$sample_id, c("a", "z"))
  expect_equal(nrow(two$manifest_no_file), 0L)
  # Shared target is reported, not hidden.
  expect_equal(nrow(two$provenance$same_target), 2L)
  expect_equal(nrow(res$provenance$same_target), 0L)
})

test_that("many links to one file are not collapsed into one row", {
  skip_on_os("windows")
  d <- withr::local_tempdir()
  real <- bind_fcs(d, "real/run.fcs")
  dir.create(file.path(d, "many"))
  ids <- paste0("m", 1:5)
  file.symlink(real, file.path(d, "many", paste0(ids, ".fcs")))
  res <- bind_call(bind_manifest(ids[1:3]), file.path(d, "many"))
  expect_equal(res$bound$sample_id, ids[1:3])
  expect_equal(nrow(res$manifest_no_file), 0L)
  expect_equal(nrow(res$file_no_manifest), 2L)
  expect_equal(res$provenance$n_files, 5L)
})

test_that("bound is a sw_manifest and keeps its manifest columns", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  man <- bind_manifest(c("a", "b"))
  res <- bind_call(man, d)
  expect_s3_class(res$bound, "sw_manifest")
  expect_true(all(names(man) %in% names(res$bound)))
})

test_that("every listed file is in bound or file_no_manifest", {
  d <- withr::local_tempdir()
  for (n in c("a", "b", "extra", "x")) {
    bind_fcs(d, paste0(n, ".fcs"))
  }
  bind_fcs(d, "m.fcs", nextdata = TRUE)
  man <- bind_manifest(c("a", "b", "m", "gone"))
  for (nd in c("bind", "unbind")) {
    res <- bind_call(man, d, nextdata = nd)
    expect_equal(
      nrow(res$bound) + nrow(res$file_no_manifest),
      res$provenance$n_files
    )
    expect_equal(
      nrow(res$bound) + nrow(res$manifest_no_file),
      nrow(man)
    )
  }
})

test_that("an unmatched multi-dataset file is reported and never stops the call", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  bind_fcs(d, "other.fcs", nextdata = TRUE)
  res <- bind_call(bind_manifest(c("a", "b")), d, nextdata = "stop")
  expect_equal(basename(res$multi_dataset$path), "other.fcs")
  expect_true(is.na(res$multi_dataset$sample_id))
  expect_false(res$multi_dataset$bound)
  expect_equal(res$multi_dataset$status, "multi_dataset")
  expect_equal(res$file_no_manifest$reason, "file_without_manifest_row")
})

test_that("nextdata = unbind names the file in file_no_manifest", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs", nextdata = TRUE)
  res <- bind_call(bind_manifest(c("a", "b")), d, nextdata = "unbind")
  expect_equal(basename(res$file_no_manifest$path), "a.fcs")
  expect_equal(res$file_no_manifest$reason, "multi_dataset_unbound")
  expect_equal(res$file_no_manifest$sample_id, "a")
  expect_equal(nrow(res$bound), 0L)
})

test_that("unbind on duplicates works with case = ignore, spelling kept", {
  d <- withr::local_tempdir()
  bind_fcs(d, "a.fcs")
  bind_fcs(d, "x/A.fcs")
  bind_fcs(d, "B.fcs")
  res <- bind_call(
    bind_manifest(c("a", "b")),
    d,
    recursive = TRUE,
    case = "ignore",
    on_duplicate_file = "unbind"
  )
  expect_equal(res$bound$sample_id, "b")
  expect_equal(basename(res$bound$path), "B.fcs")
  expect_equal(res$manifest_no_file$reason, "duplicate_file_key")
  expect_equal(res$file_no_manifest$reason, rep("duplicate_file_key", 2))
})

# Copy an FCS file with one header text replaced by a same-length text.
bind_edit_header <- function(dir, name, from, to) {
  stopifnot(nchar(from) == nchar(to))
  p <- bind_fcs(dir, name)
  r <- readBin(p, "raw", file.size(p))
  # Search bytes: the data segment holds nul bytes that rawToChar rejects.
  head_raw <- r[seq_len(min(length(r), 4000))]
  txt <- rawToChar(head_raw[head_raw != as.raw(0)])
  pos <- regexpr(from, txt, fixed = TRUE, useBytes = TRUE)
  stopifnot(pos > 0)
  pos <- which(head_raw != as.raw(0))[pos]
  r[pos:(pos + nchar(from) - 1)] <- charToRaw(to)
  writeBin(r, p)
  p
}

test_that("a missing or non-numeric $NEXTDATA is not_assessed, kept bound", {
  d <- withr::local_tempdir()
  bind_edit_header(d, "txt.fcs", "$NEXTDATA|0|", "$NEXTDATA|x|")
  bind_edit_header(d, "gone.fcs", "$NEXTDATA|0|", "$NEXTDATE|0|")
  man <- bind_manifest(c("txt", "gone"))
  # flowCore warns about the corrupt header; that warning is let through.
  res <- suppressWarnings(bind_call(man, d, nextdata = "stop"))
  expect_setequal(res$bound$sample_id, c("txt", "gone"))
  expect_equal(res$bound$header_status, c("not_assessed", "not_assessed"))
  expect_equal(res$multi_dataset$status, c("not_assessed", "not_assessed"))
  expect_true(all(nzchar(res$bound$reason)))
  # flowCore refuses a non-number itself; the missing keyword reaches our check.
  expect_match(
    res$bound$reason[res$bound$sample_id == "gone"],
    "missing or not a number"
  )
})

test_that("a non-numeric $NEXTDATA from the header reader is not_assessed", {
  testthat::local_mocked_bindings(
    read.FCSheader = \(f) list(list("$NEXTDATA" = "abc")),
    .package = "flowCore"
  )
  out <- SpectraWeaveR:::.sw_bind_nextdata("x")
  expect_equal(out$status, "not_assessed")
  expect_equal(out$nextdata, "abc")
  expect_match(out$reason, "missing or not a number")
})

test_that("only the extra-data-segment warning is muffled", {
  testthat::local_mocked_bindings(
    read.FCSheader = \(f) {
      warning("another warning")
      list(list("$NEXTDATA" = "0"))
    },
    .package = "flowCore"
  )
  expect_warning(SpectraWeaveR:::.sw_bind_read_header("x"), "another warning")
  testthat::local_mocked_bindings(
    read.FCSheader = \(f) {
      warning("found additional data segment")
      list(list("$NEXTDATA" = "0"))
    },
    .package = "flowCore"
  )
  expect_no_warning(SpectraWeaveR:::.sw_bind_read_header("x"))
})

test_that("provenance starts with the shared keys (bind)", {
  pv <- {
    d <- withr::local_tempdir()
    bind_fcs(d, "a.fcs")
    res <- bind_call(bind_manifest(c("a", "b")), d)
    expect_identical(attr(res, "provenance"), res$provenance)
    res$provenance
  }
  expect_equal(
    head(names(pv), 8),
    c(
      "package_version",
      "flowcore_version",
      "r_version",
      "platform",
      "datetime",
      "seed",
      "session",
      "args"
    )
  )
  expect_true(all(
    c("running", "locale", "attached") %in% names(pv$session)
  ))
})
