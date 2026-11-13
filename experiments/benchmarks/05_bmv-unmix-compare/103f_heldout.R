# Purpose: stage B for one AS 1.8.8 tuning setting on 24H: held-out control error (G1) and extraction reproducibility.
#          G1 as 3D (26/27): split each zip clean control 50/50 (seed 42), extract spectra and variants from the train halves
#          (AF from the full unstained file, as 3D), unmix the test halves, and score every off-target channel on the target's
#          positive events (above q99.5 of the matching universal-negative test half; all events of a negative control):
#          spread_q = (q97.5 - q2.5) / 3.92 and bias = |median|, linear scale. Reproducibility: spectra, AF spectra and
#          variants extracted twice from the full clean controls with the same seed; largest absolute difference reported.
#          The zip unstained AF file (cell rectangle of the zip script) is unmixed with the test halves and scored on all
#          channels as control "Unstained AF (Cells)" (added 2026-11-15; it is also the AF source, so it is not held out).
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_af_tuning.sbatch heldout <refine> <n_passes>
# Inputs:  <data_dir>/autoSpectral_v1.8.8/Reference Group Clean/*.fcs + fcs_control_file.csv, Reference Group/unstained_AF/ [zip]
# Outputs: <root> = experiments/runs/<SW_RUN_DATE>_bmv24h-af-heldout (shared; the arm folder must not exist)
#          <root>/<arm>/split_log.csv       file, events, train, test
#          <root>/<arm>/heldout_scores.csv  arm, control, type, target, channel, n, spread_q, bias
#          <root>/<arm>/repro.csv           object, max abs difference between two same-seed extractions
#          <root>/<arm>/sessionInfo.txt

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ------------------------------------------------------------------------------
SEED <- 42 # split and extraction seed (3D G1, user 2026-11-10)
TRAIN_FRACTION <- 0.5 # 3D G1
POS_QUANTILE <- 0.995 # 3D G1 positive rule
QUANTILES <- c(0.025, 0.975) # 3D G1 spread
Z_WIDTH <- 3.92
NEGATIVES <- c(Beads = "Negative_1 (Beads)", Cells = "Unstained (Cells)") # zip universal negatives
UNSTAINED_AF <- "Unstained AF (Cells)" # name of the gated unstained AF file among the test files
CELL_FSC <- c(6e5, 4e6) # zip script cell rectangle (inherited: zip, as 103c)
CELL_SSC <- c(1e5, 1.6e6)
PIPELINE <- "joint"
SPEED <- "slow"
cli <- commandArgs(trailingOnly = TRUE)
REFINE <- as.logical(cli[[1]])
N_PASSES <- as.integer(cli[[2]])
stopifnot("refine TRUE/FALSE" = !is.na(REFINE), "1.8.8" = as.character(packageVersion("AutoSpectral")) == "1.8.8")
arm <- sprintf("a188_r%s_p%d", REFINE, N_PASSES)
ROOT <- file.path("experiments/runs", paste0(Sys.getenv("SW_RUN_DATE"), "_bmv", cfg("run_tag"), "-af-heldout"))
dir.create(ROOT, showWarnings = FALSE)
out <- file.path(normalizePath(ROOT), arm)
stopifnot("arm folder must not exist" = !dir.exists(out))
dir.create(out, recursive = TRUE)
d188 <- file.path(normalizePath(data_path()), "autoSpectral_v1.8.8")
clean_dir <- file.path(d188, "Reference Group Clean")
unstained <- list.files(file.path(d188, "Reference Group", "unstained_AF"), pattern = "\\.fcs$", full.names = TRUE)
stopifnot("one unstained AF file" = length(unstained) == 1)
new_asp <- \(dir) {
  a <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
  f <- names(a)[grepl("dir", names(a)) & map_lgl(a, is.character)]
  a[f] <- map(a[f], \(d) file.path(dir, d))
  a
}
extract <- \(ctrl_dir, work) {
  dir.create(work, recursive = TRUE)
  asp <- new_asp(work)
  cf <- file.path(ctrl_dir, "fcs_control_file.csv")
  set.seed(SEED)
  sp <- get.spectra.automated(control.dir = ctrl_dir, control.def.file = cf, asp = asp, n.candidates = 10000, n.spectral = 500, cosine.threshold = 0.95)
  af <- get.af.spectra(unstained.sample = unstained, asp = asp, spectra = sp, refine = REFINE)
  va <- get.spectral.variants(control.dir = ctrl_dir, control.def.file = cf, asp = asp, spectra = sp, output.dir = work)
  list(spectra = sp, af = af, variants = va, asp = asp, cf = cf)
}

# Split the clean controls 50/50 -----------------------------------------------------------------
files <- sort(list.files(clean_dir, pattern = "\\.fcs$"))
dir.create(file.path(out, "train"))
dir.create(file.path(out, "test"))
set.seed(SEED)
split_log <- map(files, \(f) {
  ff <- read.FCS(file.path(clean_dir, f), transformation = FALSE, truncate_max_range = FALSE)
  n <- nrow(ff)
  tr <- sort(sample(n, floor(n * TRAIN_FRACTION)))
  write.FCS(ff[tr, ], file.path(out, "train", f))
  write.FCS(ff[-tr, ], file.path(out, "test", f))
  tibble(file = f, events = n, train = length(tr), test = n - length(tr))
}) |> list_rbind()
stopifnot("control file copied" = all(map_lgl(c("train", "test"), \(h) file.copy(file.path(clean_dir, "fcs_control_file.csv"), file.path(out, h, "fcs_control_file.csv")))))
write_csv(split_log, file.path(out, "split_log.csv"))
uf <- Subset(read.FCS(unstained, transformation = FALSE, truncate_max_range = FALSE), rectangleGate(`FSC-A` = CELL_FSC, `SSC-A` = CELL_SSC))
write.FCS(uf, file.path(out, "test", paste0(UNSTAINED_AF, ".fcs")))

# Extract from train halves, unmix test halves -------------------------------------------------
ex <- extract(file.path(out, "train"), file.path(out, "work_train"))
dir.create(file.path(out, "unmixed"))
unmix.folder(
  fcs.dir = file.path(out, "test"), spectra = ex$spectra, asp = ex$asp,
  flow.control = reload.flow.control(file.path(out, "train"), ex$cf, ex$asp), method = "AutoSpectral",
  af.spectra = ex$af, spectra.variants = ex$variants, pipeline = PIPELINE, speed = SPEED,
  n.passes = N_PASSES, n.af.passes = N_PASSES, output.dir = file.path(out, "unmixed"), parallel = TRUE, verbose = FALSE
)

# Score off-target channels ----------------------------------------------------------------------
um <- list.files(file.path(out, "unmixed"), pattern = "\\.fcs$", full.names = TRUE)
ctl <- tibble(file = um, control = str_remove(basename(um), " AutoSpectral.*\\.fcs$")) |>
  mutate(type = str_extract(control, "(?<=\\()(Beads|Cells)(?=\\))"),
         target = if_else(control %in% c(NEGATIVES, UNSTAINED_AF), NA_character_, paste0(str_remove(control, " \\((Beads|Cells)\\)$"), "-A")))
stopifnot("one unmixed file per control + unstained AF" = nrow(ctl) == length(files) + 1 && UNSTAINED_AF %in% ctl$control, "negatives present" = all(NEGATIVES %in% ctl$control))
rd <- \(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
negs <- map(NEGATIVES, \(n) rd(ctl$file[ctl$control == n]))
scores <- pmap(ctl, \(file, control, type, target) {
  e <- rd(file)
  ch <- colnames(e)[grepl("-A$", colnames(e)) & !grepl("^(FSC|SSC)|^AF", colnames(e))]
  if (!is.na(target) && !target %in% ch) stop("103f: target channel ", target, " missing for ", control, call. = FALSE)
  rows <- if (is.na(target)) rep(TRUE, nrow(e)) else e[, target] > quantile(negs[[type]][, target], POS_QUANTILE, names = FALSE)
  map(setdiff(ch, target), \(j) {
    q <- quantile(e[rows, j], QUANTILES, names = FALSE)
    tibble(arm, control, type, target, channel = j, n = sum(rows), spread_q = (q[[2]] - q[[1]]) / Z_WIDTH, bias = abs(median(e[rows, j])))
  }) |> list_rbind()
}) |> list_rbind()
write_csv(scores, file.path(out, "heldout_scores.csv"))

# Reproducibility: two same-seed extractions from the full clean controls --------------------------------
r1 <- extract(clean_dir, file.path(out, "work_rep1"))
r2 <- extract(clean_dir, file.path(out, "work_rep2"))
vdiff <- \(a, b) if (identical(names(a), names(b))) max(map2_dbl(a, b, \(x, y) if (is.numeric(x) && identical(dim(x), dim(y))) max(abs(x - y)) else NA_real_), na.rm = FALSE) else NA_real_
tibble(arm, object = c("spectra", "af_spectra", "variants"),
       identical = c(identical(r1$spectra, r2$spectra), identical(r1$af, r2$af), identical(r1$variants, r2$variants)),
       max_abs_diff = c(if (identical(dim(r1$spectra), dim(r2$spectra))) max(abs(r1$spectra - r2$spectra)) else NA_real_,
                        if (identical(dim(r1$af), dim(r2$af))) max(abs(r1$af - r2$af)) else NA_real_,
                        tryCatch(vdiff(r1$variants, r2$variants), error = \(e) NA_real_))) |>
  write_csv(file.path(out, "repro.csv"))
saveRDS(list(rep1 = r1[c("spectra", "af", "variants")], rep2 = r2[c("spectra", "af", "variants")]), file.path(out, "repro_objects.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
