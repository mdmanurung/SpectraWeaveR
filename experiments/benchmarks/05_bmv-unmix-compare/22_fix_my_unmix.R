# Purpose: AutoSpectral-documented correction of bead-derived spectra applied to cells, in two steps, on the
#          as188_joint spectra: correct.unmixing.signatures() -> arm sig188; fix.my.unmix() on its output -> arm fix188.
#          Each arm unmixes the 5 raw BMV files with the corrected spectra and otherwise the as188_joint settings.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/22_fix_my_unmix.R <out_root> [smoke_pattern]
#          from the repo root. <out_root>/sig188 and <out_root>/fix188 must not exist (nothing is overwritten).
#          Full run: <out_root> = experiments/bmv_test/compare_3D/data/unmix_cells
#          Smoke test: <out_root> = a scratch folder and smoke_pattern = "^A3 " (one donor for correction and unmixing)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8 (Reference Group Clean controls, run-143803 spectra,
#          AF csv, Spectral_variants.rds, Raw, Reference Group/unstained_AF/H12, AutoSpectral_unmixed as the reference
#          for event counts and channel names) [user zip]
# Outputs: <out_root>/sig188 and <out_root>/fix188, each with 5 unmixed *.fcs, args.txt, verify.csv, sessionInfo.txt,
#          spectra.csv (spectra used for that arm) and the function results as rds/csv
#          (sig188/correct_signatures_*, fix188/fix_my_unmix_*, fix188/fix_my_unmix_files/ = files the function writes itself).
# Approved by the user on 2026-11-10 (task text). Named values come from that task; every other argument is the
# AutoSpectral 1.8.8 default. Both functions' Rd were read: fix.my.unmix() takes RAW fully stained files and the
# spectra that correct.unmixing.signatures() returns; it returns corrected spectra, not unmixed events, so the
# fix188 events come from unmix.folder() with those spectra. correct.unmixing.signatures() reads ONE fully stained
# file by path, so the 5 raw files are pooled by hand and passed as `raw.data` (the form its Rd documents; it then
# wants the unstained matrix in `unstained`, and `unstained.sample`/`fully.stained.sample` are not passed).

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
PIPELINE <- "joint" # as188_joint: repro_s188_c188_fast/args.txt
SPEED <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint" (repro188 = s188_c188_slow); was "fast" to match as188_joint
UNSTAINED_PATTERN <- "H12 Unstained" # matched unstained sample (task)
CUS_BG_MODE <- NULL # not named in the task: package default (c("global.mean", "scatter.knn", "none")[1])
FMU_BG_MODE <- "per.cell" # named in the task (also the package default)
FMU_ESTIMATOR <- "anchored" # named in the task (also the package default)

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: <out_root> [smoke_pattern]" = length(args) %in% 1:2)
out_root <- normalizePath(args[[1]], mustWork = FALSE)
smoke_pattern <- if (length(args) == 2) args[[2]] else NA_character_
dir_sig <- file.path(out_root, "sig188")
dir_fix <- file.path(out_root, "fix188")
stopifnot(
  "sig188 must not exist" = !dir.exists(dir_sig),
  "fix188 must not exist" = !dir.exists(dir_fix)
)
dir.create(dir_sig, recursive = TRUE)
dir.create(dir_fix, recursive = TRUE)
suppressPackageStartupMessages({
  library(AutoSpectral)
  library(dplyr)
})
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)
t0 <- Sys.time()

# Inputs: same as 06_unmix_c188.R cell s188 ------------------------------------------------
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
control_dir <- file.path(d188, "Reference Group Clean")
spectra_dir <- file.path(d188, "table_spectra")
spectra_file <- "Automated_autospectral_spectra_20261105_143803.csv"
af_file <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
variants_rds <- file.path(
  d188,
  "figure_spectral_variants",
  "Spectral_variants.rds"
)
raw_dir <- file.path(d188, "Raw")
ref_dir <- file.path(d188, "AutoSpectral_unmixed")
unstained_dir <- file.path(d188, "Reference Group", "unstained_AF")

asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(
  control_dir,
  file.path(control_dir, "fcs_control_file.csv"),
  asp
)
spectra <- read.spectra(spectra_file, spectra.dir = spectra_dir)
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in 06; no AF row in the 1.8.8 file
af.spectra <- read.spectra(af_file, spectra.dir = spectra_dir)
spectra.variants <- readRDS(variants_rds)
stopifnot("af.spectra has >= 2 rows" = nrow(af.spectra) >= 2)

raw_files <- sort(list.files(raw_dir, pattern = "\\.fcs$", full.names = TRUE))
stopifnot("5 raw files" = length(raw_files) == 5)
unstained_files <- list.files(
  unstained_dir,
  pattern = UNSTAINED_PATTERN,
  full.names = TRUE
)
stopifnot("one H12 unstained file" = length(unstained_files) == 1)
unstained_file <- unstained_files[[1]]
if (!is.na(smoke_pattern)) {
  raw_files <- raw_files[grepl(smoke_pattern, basename(raw_files))]
  stopifnot("smoke pattern matches one file" = length(raw_files) == 1)
}

# Raw matrices for correct.unmixing.signatures() -------------------------------------------
read_raw <- get("readFCS", envir = asNamespace("AutoSpectral"))
read_det <- function(f) {
  x <- read_raw(f, columns = flow.control$scatter.and.channel.spectral)
  as.matrix(x[, flow.control$spectral.channel, drop = FALSE])
}
raw_pool <- do.call(rbind, lapply(raw_files, read_det))
unstained_raw <- read_det(unstained_file)
stopifnot(
  "detectors match spectra" = identical(colnames(raw_pool), colnames(spectra)),
  "unstained detectors match" = identical(
    colnames(unstained_raw),
    colnames(spectra)
  )
)
message(
  "pooled stained events: ",
  nrow(raw_pool),
  "; unstained events: ",
  nrow(unstained_raw)
)

# Step 1: correct.unmixing.signatures() -----------------------------------------------------
cus_args <- list(
  spectra = spectra,
  raw.data = raw_pool,
  unstained = unstained_raw,
  asp = asp,
  verbose = TRUE
)
if (!is.null(CUS_BG_MODE)) {
  cus_args$bg.mode <- CUS_BG_MODE
}
cus <- do.call(correct.unmixing.signatures, cus_args)
spectra_sig <- cus$spectra

# Step 2: fix.my.unmix() on the corrected spectra --------------------------------------------
fmu_files_dir <- file.path(dir_fix, "fix_my_unmix_files")
asp_fmu <- asp
asp_fmu$fix.unmixing.dir <- fmu_files_dir # where the function writes its own csv and heatmap (output location only)
fmu <- fix.my.unmix(
  spectra = spectra_sig,
  unstained.sample = unstained_file,
  fully.stained.sample = raw_files,
  flow.control = flow.control,
  asp = asp_fmu,
  variants = spectra.variants,
  af.spectra = af.spectra,
  bg.mode = FMU_BG_MODE,
  estimator = FMU_ESTIMATOR
)
spectra_fix <- fmu$spectra

# Change tables ----------------------------------------------------------------------------
angle_deg <- function(a, b) {
  acos(pmin(1, sum(a * b) / sqrt(sum(a^2) * sum(b^2)))) * 180 / pi
}
spectra_change <- function(from, to) {
  stopifnot(identical(rownames(from), rownames(to)))
  tibble::tibble(
    fluorophore = rownames(from),
    angle_deg = vapply(
      seq_len(nrow(from)),
      \(i) angle_deg(from[i, ], to[i, ]),
      0
    ),
    max_abs_change = apply(abs(from - to), 1, max)
  )
}
chg_sig <- spectra_change(spectra, spectra_sig) |>
  mutate(moved = angle_deg > 0)
chg_fix <- spectra_change(spectra_sig, spectra_fix) |>
  mutate(moved = angle_deg > 0)
chg_fix_total <- spectra_change(spectra, spectra_fix) |>
  mutate(moved = angle_deg > 0)
off <- fmu$spillover
diag(off) <- 0
n_coef_nonzero <- sum(abs(off) > 0)
cat(
  sprintf(
    "signatures moved by step 1: %d of %d\n",
    sum(chg_sig$moved),
    nrow(chg_sig)
  ),
  sprintf(
    "signatures moved by step 2: %d of %d\n",
    sum(chg_fix$moved),
    nrow(chg_fix)
  ),
  sprintf(
    "fix.my.unmix off-diagonal spillover coefficients nonzero: %d of %d\n",
    n_coef_nonzero,
    length(off) - nrow(off)
  ),
  sep = ""
)

# Unmix with each set of spectra -----------------------------------------------------------
unmix_arm <- function(spectra_arm, out_dir, smoke_pattern) {
  in_dir <- raw_dir
  if (!is.na(smoke_pattern)) {
    in_dir <- file.path(
      dirname(out_dir),
      paste0("raw_smoke_", basename(out_dir))
    )
    dir.create(in_dir)
    file.symlink(raw_files, file.path(in_dir, basename(raw_files)))
  }
  unmix.folder(
    fcs.dir = in_dir,
    spectra = spectra_arm,
    asp = asp,
    flow.control = flow.control,
    method = "AutoSpectral",
    af.spectra = af.spectra,
    spectra.variants = spectra.variants,
    pipeline = PIPELINE,
    speed = SPEED,
    output.dir = out_dir,
    parallel = TRUE,
    verbose = TRUE
  )
}
unmix_arm(spectra_sig, dir_sig, smoke_pattern)
unmix_arm(spectra_fix, dir_fix, smoke_pattern)

# Verify against the as188_joint files -------------------------------------------------------
fcs_n <- function(f) as.integer(flowCore::read.FCSheader(f)[[1]][["$TOT"]])
fcs_names <- function(f) {
  colnames(flowCore::read.FCS(f, truncate_max_range = FALSE, which.lines = 1))
}
verify_arm <- function(out_dir, arm) {
  tibble::tibble(raw = raw_files) |>
    mutate(
      arm = arm,
      ref_file = file.path(
        ref_dir,
        paste(tools::file_path_sans_ext(basename(raw)), "AutoSpectral.fcs")
      ),
      out = file.path(
        out_dir,
        paste(tools::file_path_sans_ext(basename(raw)), "AutoSpectral.fcs")
      ),
      n_raw = vapply(raw, fcs_n, 0L),
      n_ref = vapply(ref_file, fcs_n, 0L),
      n_out = vapply(
        out,
        \(x) if (file.exists(x)) fcs_n(x) else NA_integer_,
        0L
      ),
      same_n = n_out == n_ref,
      same_names = vapply(
        out,
        \(x) {
          if (file.exists(x)) {
            identical(fcs_names(x), fcs_names(ref_file[match(x, out)]))
          } else {
            NA
          }
        },
        TRUE
      )
    )
}
ver_sig <- verify_arm(dir_sig, "sig188")
ver_fix <- verify_arm(dir_fix, "fix188")
write.csv(ver_sig, file.path(dir_sig, "verify.csv"), row.names = FALSE)
write.csv(ver_fix, file.path(dir_fix, "verify.csv"), row.names = FALSE)
for (v in list(ver_sig, ver_fix)) {
  stopifnot(
    "all output files exist" = !anyNA(v$n_out),
    "event counts match as188_joint" = all(v$same_n),
    "channel names match as188_joint" = all(v$same_names)
  )
}

# Save results and arguments -----------------------------------------------------------------
write.csv(spectra_sig, file.path(dir_sig, "spectra.csv"))
write.csv(spectra_fix, file.path(dir_fix, "spectra.csv"))
write.csv(spectra, file.path(dir_sig, "spectra_start_as188.csv"))
write.csv(
  cus$fit.log,
  file.path(dir_sig, "correct_signatures_fit_log.csv"),
  row.names = FALSE
)
write.csv(
  cus$proposed.spectra,
  file.path(dir_sig, "correct_signatures_proposed_spectra.csv")
)
write.csv(
  data.frame(
    fluorophore = names(cus$accepted),
    accepted = unname(cus$accepted)
  ),
  file.path(dir_sig, "correct_signatures_accepted.csv"),
  row.names = FALSE
)
write.csv(
  chg_sig,
  file.path(dir_sig, "spectra_change_vs_start.csv"),
  row.names = FALSE
)
saveRDS(
  cus[setdiff(names(cus), "dominant")],
  file.path(dir_sig, "correct_signatures_result.rds")
)
saveRDS(cus$dominant, file.path(dir_sig, "correct_signatures_dominant.rds"))

for (nm in c("spillover", "compensation", "trust")) {
  write.csv(fmu[[nm]], file.path(dir_fix, paste0("fix_my_unmix_", nm, ".csv")))
}
for (nm in c(
  "coefficient.log",
  "source.log",
  "signature.log",
  "convergence.log",
  "smooth.log",
  "af.assignment"
)) {
  if (!is.null(fmu[[nm]])) {
    write.csv(
      fmu[[nm]],
      file.path(dir_fix, paste0("fix_my_unmix_", nm, ".csv")),
      row.names = FALSE
    )
  }
}
write.csv(
  fmu$spectra.backsolved,
  file.path(dir_fix, "fix_my_unmix_spectra_backsolved.csv")
)
write.csv(
  chg_fix,
  file.path(dir_fix, "spectra_change_vs_sig188.csv"),
  row.names = FALSE
)
write.csv(
  chg_fix_total,
  file.path(dir_fix, "spectra_change_vs_start_as188.csv"),
  row.names = FALSE
)
saveRDS(
  fmu[setdiff(names(fmu), c("unmixed.final", "residual.final"))],
  file.path(dir_fix, "fix_my_unmix_result.rds")
)

fmt_formals <- function(f) {
  fm <- formals(f)
  vapply(
    names(fm),
    \(n) {
      paste0(
        n,
        " = ",
        paste(deparse(fm[[n]], width.cutoff = 500), collapse = " ")
      )
    },
    ""
  )
}
named_cus <- c(
  "spectra = starting as188_joint spectra, AF row dropped (task)",
  paste0(
    "raw.data = rbind of the detector columns of: ",
    paste(basename(raw_files), collapse = "; "),
    " (task: the 5 raw stained files; pooled by hand because the file-path argument takes one file)"
  ),
  paste0(
    "unstained = detector columns of ",
    unstained_file,
    " (task: H12; passed as a matrix because raw.data is used)"
  ),
  "asp = get.autospectral.param(cytometer = 'aurora', figures = FALSE) (used for the seed only)",
  "verbose = TRUE (package default)",
  paste0(
    "bg.mode: ",
    if (is.null(CUS_BG_MODE)) "not passed, package default" else CUS_BG_MODE
  )
)
named_fmu <- c(
  "spectra = corrected spectra returned by correct.unmixing.signatures() (Rd: 'ideally the output of')",
  paste0("unstained.sample = ", unstained_file, " (task)"),
  paste0(
    "fully.stained.sample = ",
    paste(raw_files, collapse = "; "),
    " (task)"
  ),
  "flow.control = reload.flow.control(Reference Group Clean) (as 06_unmix_c188.R)",
  paste0(
    "asp = as above with fix.unmixing.dir = ",
    fmu_files_dir,
    " (only where the function writes its own csv and heatmap)"
  ),
  paste0("variants = ", variants_rds, " (task: as188_joint variants)"),
  paste0(
    "af.spectra = ",
    file.path(spectra_dir, af_file),
    " (as188_joint AF; required by bg.mode = 'per.cell')"
  ),
  paste0("bg.mode = ", FMU_BG_MODE, " (task)"),
  paste0("estimator = ", FMU_ESTIMATOR, " (task)")
)
common <- c(
  paste("pipeline:", PIPELINE, "(as188_joint)"),
  paste("speed:", SPEED, "(as188_joint)"),
  paste("unmix.folder af.spectra:", file.path(spectra_dir, af_file)),
  paste("unmix.folder spectra.variants:", variants_rds),
  paste("unmix.folder parallel: TRUE, verbose: TRUE (as 06_unmix_c188.R)"),
  paste("start spectra:", file.path(spectra_dir, spectra_file)),
  paste("controls:", control_dir),
  paste("raw dir:", raw_dir),
  if (!is.na(smoke_pattern)) {
    paste("SMOKE TEST, raw file pattern:", smoke_pattern)
  }
)
writeLines(
  c(
    "arm: sig188 (unmixed with the spectra from correct.unmixing.signatures())",
    common,
    "",
    "## correct.unmixing.signatures arguments set here",
    named_cus,
    "",
    "## correct.unmixing.signatures all formals (package defaults unless listed above)",
    fmt_formals(correct.unmixing.signatures)
  ),
  file.path(dir_sig, "args.txt")
)
writeLines(
  c(
    "arm: fix188 (unmixed with the spectra from fix.my.unmix() run on the sig188 spectra)",
    common,
    "",
    "## fix.my.unmix arguments set here",
    named_fmu,
    "",
    "## fix.my.unmix all formals (package defaults unless listed above)",
    fmt_formals(fix.my.unmix)
  ),
  file.path(dir_fix, "args.txt")
)
for (d in c(dir_sig, dir_fix)) {
  writeLines(capture.output(sessionInfo()), file.path(d, "sessionInfo.txt"))
}
message(
  "runtime (min): ",
  round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
)
