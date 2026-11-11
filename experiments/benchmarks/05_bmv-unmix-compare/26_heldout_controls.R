# Purpose: G1 held-out controls. Split every control 50/50 at random, re-extract spectra, AF spectra and
#          variants from the TRAIN halves with each AutoSpectral code line, and unmix the TEST halves (plus
#          the full H12 unstained file) with five arms, so the controls are scored on events that did not
#          build the spectra (user plan G1, 2026-11-10).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/26_heldout_controls.R <step> [held_root]
#          from the repo root. <step> runs in its own R process because the two code lines need different
#          AutoSpectral libraries. Order:
#            split        write train/ and test/ halves (seed 42)
#            spectra_188  1.8.8 route: get.spectra.automated, get.af.spectra(refine = FALSE), get.spectral.variants
#            spectra_102  1.0.2 route: define.flow.control, clean.controls (Cells Matching Negative), spectra, AF, variants
#            as102_basic | as102_var | as188_joint | s188_c102 | s102_c188   unmix test halves + H12
#          With held_root, everything goes to that folder and only SMOKE_FILES are used (smoke test).
# Settings come from (nothing invented):
#   1.8.8 route  AutoSpectral_unmix_v1.8.8.R (n.candidates 10000, n.spectral 500, cosine.threshold 0.95, other
#                arguments default); AF from the H12 file with refine = FALSE; variants from the controls.
#   1.0.2 route  Autospectral_Setup.R (cleaning with universal negative, af.remove = FALSE; AF from H12; variants
#                with parallel = FALSE). Its own control file is used (differs from the 1.8.8 file only in
#                large.gate = TRUE for cells), reduced to the files present.
#   unmix arms   16_unmix_controls.R (same arguments per arm); flow.control is built from the train folder with
#                the control file of the line that made the spectra, as 16 does with ctrl_102 / ctrl_188.
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/*.fcs + fcs_control_file.csv
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group/unstained_AF/H12 *.fcs
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.0.2/Reference Group/fcs_control_file.csv
#          experiments/bmv_test/compare_3D/data/unmix_controls/<arm>/ (channel names of the full-data arm)
#          experiments/bmv_test/compare_3D/data/lib_as102 (AutoSpectral and AutoSpectralRcpp 1.0.2)
# Outputs: <held_root>/train/, test/             halves with the same file names and control file
#          <held_root>/split_log.csv             file order, event counts, seed
#          <held_root>/spectra_c188/, spectra_c102/  spectra.rds, af_spectra.rds, spectral_variants.rds, csv copies, elapsed.txt
#          <held_root>/unmixed/<arm>/            unmixed test halves + H12, args.txt, verify.csv, sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
LIB_102 <- file.path(DATA, "lib_as102") # AutoSpectral and AutoSpectralRcpp v1.0.2 (separate library)
SPLIT_SEED <- 42 # user plan G1: set once, files taken in sorted name order
TRAIN_FRACTION <- 0.5 # user plan G1: 50/50 split; the test half gets the extra event of an odd count
EXTRACT_SEED <- 42 # AutoSpectral subsamples events; same seed as the split so reruns match
N_CANDIDATES <- 10000L # AutoSpectral_unmix_v1.8.8.R (shipped script)
N_SPECTRAL <- 500L # same
COSINE_THRESHOLD <- 0.95 # same
WEIGHTED_102 <- FALSE # as in 06_unmix_c102.R and 16_unmix_controls.R
PIPELINE_188 <- "joint" # as in 06_unmix_c188.R and 16_unmix_controls.R
SPEED_188 <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint" (repro188 = s188_c188_slow); was "fast" to match as188_joint
SMOKE_FILES <- c(
  "Zombie NIR (Cells).fcs",
  "Unstained (Cells).fcs",
  "PE (Beads).fcs",
  "BV605 (Beads).fcs",
  "Negative_1 (Beads).fcs"
) # whole units: a cell and bead control each with its universal negative

# arm: code line, which line made the spectra, variants used (16_unmix_controls.R, ARMS)
ARMS <- list(
  as102_basic = list(code = "1.0.2", spectra = "c102", variants = FALSE),
  as102_var = list(code = "1.0.2", spectra = "c102", variants = TRUE),
  as188_joint = list(code = "1.8.8", spectra = "c188", variants = TRUE),
  s188_c102 = list(code = "1.0.2", spectra = "c188", variants = TRUE),
  s102_c188 = list(code = "1.8.8", spectra = "c102", variants = TRUE)
)
STEPS <- c("split", "spectra_188", "spectra_102", names(ARMS))

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
step <- args[[1]]
smoke <- length(args) >= 2
held_root <- if (smoke) args[[2]] else file.path(DATA, "heldout_controls")
stopifnot("step is known" = step %in% STEPS)
dir.create(held_root, recursive = TRUE, showWarnings = FALSE)
held_root <- normalizePath(held_root)
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
d102 <- file.path(DATA, "autoSpectral_v1.0.2")
ctrl_clean <- file.path(d188, "Reference Group Clean")
h12_file <- file.path(
  d188,
  "Reference Group",
  "unstained_AF",
  "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs"
)
train_dir <- file.path(held_root, "train")
test_dir <- file.path(held_root, "test")
ctrl_csv <- "fcs_control_file.csv"
code_102 <- step %in% c("spectra_102", "as102_basic", "as102_var", "s188_c102")
if (code_102) {
  .libPaths(c(LIB_102, .libPaths()))
}
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(AutoSpectral)
  library(flowCore)
})
if (code_102) {
  suppressPackageStartupMessages(library(AutoSpectralRcpp))
}
stopifnot(
  "requested AutoSpectral version loaded" = as.character(packageVersion(
    "AutoSpectral"
  )) ==
    if (code_102) "1.0.2" else "1.8.8"
)
read_ff <- \(f) read.FCS(f, truncate_max_range = FALSE, transformation = FALSE)

# Step: split -----------------------------------------------------------------------------
if (step == "split") {
  stopifnot(
    "train and test must not exist" = !dir.exists(train_dir) &&
      !dir.exists(test_dir)
  )
  dir.create(train_dir)
  dir.create(test_dir)
  files <- sort(
    list.files(ctrl_clean, pattern = "\\.fcs$"),
    method = "radix"
  )
  if (smoke) {
    files <- files[files %in% SMOKE_FILES]
  }
  set.seed(SPLIT_SEED) # once, then one sample.int per file in the order of `files`
  split_log <- map(seq_along(files), \(i) {
    ff <- read_ff(file.path(ctrl_clean, files[[i]]))
    n <- nrow(ff)
    # sorted so Time stays monotone inside each half
    train_idx <- sort(sample.int(n, floor(n * TRAIN_FRACTION)))
    test_idx <- setdiff(seq_len(n), train_idx)
    write.FCS(ff[train_idx, ], file.path(train_dir, files[[i]]))
    write.FCS(ff[test_idx, ], file.path(test_dir, files[[i]]))
    tibble(
      order = i,
      file = files[[i]],
      n = n,
      n_train = length(train_idx),
      n_test = length(test_idx),
      seed = SPLIT_SEED
    )
  }) |>
    list_rbind()
  # same control file in both halves (rows reduced to the files present in a smoke run)
  ctrl <- read.csv(file.path(ctrl_clean, ctrl_csv), check.names = FALSE)
  ctrl <- ctrl[ctrl$filename %in% files, ]
  write.csv(ctrl, file.path(train_dir, ctrl_csv), row.names = FALSE)
  write.csv(ctrl, file.path(test_dir, ctrl_csv), row.names = FALSE)
  write_csv(split_log, file.path(held_root, "split_log.csv"))
  stopifnot(
    "halves add up to the full file" = all(
      split_log$n_train + split_log$n_test == split_log$n
    )
  )
}

# Step: spectra from the train halves, 1.8.8 route ------------------------------------------
if (step == "spectra_188") {
  out_dir <- file.path(held_root, "spectra_c188")
  stopifnot("out dir must not exist" = !dir.exists(out_dir))
  dir.create(out_dir)
  setwd(out_dir) # AutoSpectral writes tables relative to the working directory
  asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
  control_file <- file.path(train_dir, ctrl_csv)
  t0 <- Sys.time()
  set.seed(EXTRACT_SEED)
  check.control.file(train_dir, control_file, asp)
  spectra <- get.spectra.automated(
    control.dir = train_dir,
    control.def.file = control_file,
    asp = asp,
    n.candidates = N_CANDIDATES,
    n.spectral = N_SPECTRAL,
    cosine.threshold = COSINE_THRESHOLD,
    figures = FALSE
  )
  af.spectra <- get.af.spectra(
    unstained.sample = h12_file,
    asp = asp,
    spectra = spectra,
    refine = FALSE,
    figures = FALSE
  )
  spectra.variants <- get.spectral.variants(
    control.dir = train_dir,
    control.def.file = control_file,
    asp = asp,
    spectra = spectra,
    figures = FALSE
  )
  saveRDS(spectra, file.path(out_dir, "spectra.rds"))
  saveRDS(af.spectra, file.path(out_dir, "af_spectra.rds"))
  saveRDS(spectra.variants, file.path(out_dir, "spectral_variants.rds"))
  write.csv(spectra, file.path(out_dir, "spectra.csv"))
  write.csv(af.spectra, file.path(out_dir, "af_spectra.csv"))
  writeLines(
    c(
      paste("elapsed_min:", round(difftime(Sys.time(), t0, units = "mins"), 1)),
      paste("seed:", EXTRACT_SEED)
    ),
    file.path(out_dir, "elapsed.txt")
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )
}

# Step: spectra from the train halves, 1.0.2 route -------------------------------------------
if (step == "spectra_102") {
  out_dir <- file.path(held_root, "spectra_c102")
  stopifnot("out dir must not exist" = !dir.exists(out_dir))
  dir.create(out_dir)
  setwd(out_dir)
  # the 1.0.2 control file (large.gate TRUE for cells), reduced to the files in the train folder
  ctrl <- read.csv(
    file.path(d102, "Reference Group", ctrl_csv),
    check.names = FALSE
  )
  ctrl <- ctrl[ctrl$filename %in% list.files(train_dir, pattern = "\\.fcs$"), ]
  control_file <- file.path(out_dir, ctrl_csv)
  write.csv(ctrl, control_file, row.names = FALSE)
  # figures = TRUE as in Autospectral_Setup.R: define.flow.control in 1.0.2 saves its gate plots and fails without a figure folder
  asp <- get.autospectral.param(cytometer = "aurora", figures = TRUE)
  t0 <- Sys.time()
  set.seed(EXTRACT_SEED)
  check.control.file(train_dir, control_file, asp)
  flow.control <- define.flow.control(train_dir, control_file, asp)
  clean.control.cells.neg <- clean.controls(
    flow.control,
    asp,
    universal.negative = TRUE,
    af.remove = FALSE,
    main.figures = FALSE
  )
  spectra <- get.fluorophore.spectra(
    clean.control.cells.neg,
    asp,
    use.clean.expr = TRUE,
    title = "Cells Matching Negative",
    figures = FALSE
  )
  af.spectra <- get.af.spectra(
    h12_file,
    asp,
    spectra,
    title = "H12 AF",
    figures = FALSE
  )
  spectra.variants <- get.spectral.variants(
    train_dir,
    control_file,
    asp,
    spectra,
    af.spectra,
    figures = FALSE,
    parallel = FALSE,
    verbose = TRUE
  )
  saveRDS(spectra, file.path(out_dir, "spectra.rds"))
  saveRDS(af.spectra, file.path(out_dir, "af_spectra.rds"))
  saveRDS(spectra.variants, file.path(out_dir, "spectral_variants.rds"))
  write.csv(spectra, file.path(out_dir, "spectra.csv"))
  write.csv(af.spectra, file.path(out_dir, "af_spectra.csv"))
  writeLines(
    c(
      paste("elapsed_min:", round(difftime(Sys.time(), t0, units = "mins"), 1)),
      paste("seed:", EXTRACT_SEED)
    ),
    file.path(out_dir, "elapsed.txt")
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )
}

# Step: unmix the test halves and H12 with one arm ------------------------------------------
if (step %in% names(ARMS)) {
  a <- ARMS[[step]]
  out_dir <- file.path(held_root, "unmixed", step)
  stopifnot("out dir must not exist" = !dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)
  src <- file.path(held_root, paste0("spectra_", a$spectra))
  # unmix.folder takes a folder: link the test halves and the full H12 file into one staging folder
  test_files <- list.files(test_dir, pattern = "\\.fcs$", full.names = TRUE)
  raw_files <- c(test_files, h12_file)
  stage_dir <- file.path(tempdir(), "controls_heldout")
  dir.create(stage_dir)
  stopifnot("file names are unique" = !anyDuplicated(basename(raw_files)))
  file.symlink(raw_files, file.path(stage_dir, basename(raw_files)))

  asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
  control_file <- if (a$spectra == "c102") {
    file.path(src, ctrl_csv)
  } else {
    file.path(train_dir, ctrl_csv)
  }
  flow.control <- reload.flow.control(train_dir, control_file, asp)
  spectra <- readRDS(file.path(src, "spectra.rds"))
  spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in the 06/07/16 scripts
  af.spectra <- readRDS(file.path(src, "af_spectra.rds"))
  spectra.variants <- if (a$variants) {
    readRDS(file.path(src, "spectral_variants.rds"))
  }

  unmix_args <- list(
    fcs.dir = stage_dir,
    spectra = spectra,
    asp = asp,
    flow.control = flow.control,
    method = "AutoSpectral",
    af.spectra = af.spectra,
    spectra.variants = spectra.variants,
    output.dir = out_dir,
    parallel = TRUE
  )
  unmix_args <- if (a$code == "1.0.2") {
    c(unmix_args, list(weighted = WEIGHTED_102, threads = NULL))
  } else {
    c(
      unmix_args,
      list(pipeline = PIPELINE_188, verbose = TRUE, speed = SPEED_188)
    )
  }
  writeLines(
    c(
      paste("arm:", step),
      paste("AutoSpectral:", a$code),
      paste("spectra from:", src),
      paste("variants:", a$variants),
      paste("flow.control from:", train_dir, "with", control_file),
      paste("weighted:", if (a$code == "1.0.2") WEIGHTED_102 else "n/a"),
      paste("pipeline:", if (a$code == "1.0.2") "n/a" else PIPELINE_188),
      paste(
        "speed:",
        if (a$code == "1.0.2") "package default (slow)" else SPEED_188
      ),
      paste("files unmixed (n):", length(raw_files)),
      paste("files:", raw_files)
    ),
    file.path(out_dir, "args.txt")
  )
  t0 <- Sys.time()
  do.call(unmix.folder, unmix_args)
  cat(
    "elapsed_min:",
    format(round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)),
    "\n",
    file = file.path(out_dir, "args.txt"),
    append = TRUE
  )

  # Verify: one output per input, same event count, same channel names as the full-data arm -----
  ref_file <- list.files(
    file.path(DATA, "unmix_controls", step),
    pattern = "AutoSpectral.*\\.fcs$",
    full.names = TRUE
  )[[1]]
  fcs_names <- \(f) {
    colnames(read.FCS(f, truncate_max_range = FALSE, which.lines = 1))
  }
  fcs_n <- \(f) as.integer(read.FCSheader(f)[[1]][["$TOT"]])
  ref_names <- fcs_names(ref_file)
  out_files <- list.files(out_dir, pattern = "\\.fcs$", full.names = TRUE)
  verify <- tibble(raw = raw_files) |>
    mutate(
      # 1.0.2 names the output after the $FIL keyword, 1.8.8 after the file on disk: accept either
      out = map_chr(raw, \(x) {
        stems <- c(
          tools::file_path_sans_ext(basename(x)),
          tools::file_path_sans_ext(read.FCSheader(x)[[1]][["$FIL"]])
        )
        paste(
          out_files[basename(out_files) %in% paste(stems, "AutoSpectral.fcs")],
          collapse = ";"
        )
      }),
      n_raw = map_int(raw, fcs_n),
      n_out = map_int(out, \(x) if (x == "") NA_integer_ else fcs_n(x)),
      same_names = map_lgl(out, \(x) {
        if (x == "") NA else identical(fcs_names(x), ref_names)
      }),
      ref_file = ref_file
    )
  write.csv(verify, file.path(out_dir, "verify.csv"), row.names = FALSE)
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )
  stopifnot(
    "one output per input with the same event count" = all(
      verify$n_raw == verify$n_out,
      na.rm = FALSE
    )
  )
}
