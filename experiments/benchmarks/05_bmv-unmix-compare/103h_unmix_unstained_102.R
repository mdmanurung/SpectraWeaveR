# Purpose: give the two AS 1.0.2 arms the same data-driven edges as the tuning arms (user 2026-11-14): unmix the zip's
#          unstained AF file (cell rectangle of the 1.8.8 zip script) with the zip's own 1.0.2 objects and calls
#          (ReRunUnmix.R: Cells Matching Negative spectra, AF csv, weighted FALSE; as102_var adds Spectral_variants.rds),
#          then per channel q0.5 = floor, q99.5 = edge (arcsinh, cofactor 6000).
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_af_tuning.sbatch unstained102
# Inputs:  <data_dir>/autoSpectral_v1.0.2/Reference Group/ (control file, table_spectra, figure_spectral_variants) [zip]
#          <data_dir>/autoSpectral_v1.8.8/Reference Group/unstained_AF/*.fcs [zip]
#          experiments/bmv_test/compare_3D/data/lib_as102   AutoSpectral and AutoSpectralRcpp 1.0.2
# Outputs: <out_dir> = run_out("unstained-102")
#          <out_dir>/<arm>/*.fcs, thresholds.csv (channel, n, floor, edge); <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(normalizePath("experiments/bmv_test/compare_3D/data/lib_as102"), .libPaths()))
suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
CELL_FSC <- c(6e5, 4e6) # zip 1.8.8 cell rectangle, as 103c
CELL_SSC <- c(1e5, 1.6e6)
PROBS <- c(floor = 0.005, edge = 0.995) # user 2026-11-14
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
stopifnot("1.0.2 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.0.2")
OUT_DIR <- run_out("unstained-102")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE)
OUT_DIR <- normalizePath(OUT_DIR) # absolute: the 1.0.2 calls run from the zip folder

# Zip 1.0.2 objects, as ReRunUnmix.R --------------------------------------------------------------
d102 <- file.path(normalizePath(data_path()), "autoSpectral_v1.0.2", "Reference Group")
old <- setwd(d102)
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
flow.control <- reload.flow.control(d102, "fcs_control_file.csv", asp)
spectra <- read.spectra("Cells Matching Negative_autospectral_spectra.csv")
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE]
af <- read.spectra(list.files("table_spectra", pattern = "\\_AF.*\\.csv$"))
variants <- readRDS(file.path(d102, "figure_spectral_variants", "Spectral_variants.rds"))
setwd(old)

# Unstained file in the cell rectangle ----------------------------------------------------------
unstained <- list.files(data_path("autoSpectral_v1.8.8/Reference Group/unstained_AF"), pattern = "\\.fcs$", full.names = TRUE)
stopifnot("one unstained AF file" = length(unstained) == 1)
ff <- read.FCS(unstained, transformation = FALSE, truncate_max_range = FALSE)
ff <- Subset(ff, rectangleGate(`FSC-A` = CELL_FSC, `SSC-A` = CELL_SSC))
keyword(ff)[["$FIL"]] <- basename(unstained) # 1.0.2 names its output from $FIL
gated <- file.path(tempdir(), "unstained_gated")
dir.create(gated)
write.FCS(ff, file.path(gated, basename(unstained)))

# Unmix per arm and write thresholds -------------------------------------------------------------
ARMS <- list(as102_basic = NULL, as102_var = variants)
iwalk(ARMS, \(v, arm) {
  out <- file.path(OUT_DIR, arm)
  dir.create(out)
  setwd(d102)
  args <- list(fcs.dir = gated, spectra = spectra, asp = asp, flow.control = flow.control, method = "AutoSpectral",
               af.spectra = af, weighted = FALSE, output.dir = out, parallel = TRUE, threads = NULL)
  if (!is.null(v)) args$spectra.variants <- v
  do.call(unmix.folder, args)
  setwd(old)
  um <- list.files(out, pattern = "\\.fcs$", full.names = TRUE)
  stopifnot("one unmixed unstained file" = length(um) == 1)
  m <- exprs(read.FCS(um, transformation = FALSE, truncate_max_range = FALSE))
  chans <- colnames(m)[grepl("-A$", colnames(m)) & !grepl("^(FSC|SSC)", colnames(m))]
  map(chans, \(ch) {
    x <- asinh(m[, ch] / COFACTOR)
    tibble(channel = ch, n = length(x), floor = unname(quantile(x, PROBS[["floor"]])), edge = unname(quantile(x, PROBS[["edge"]])))
  }) |>
    list_rbind() |>
    write_csv(file.path(out, "thresholds.csv"))
})
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
