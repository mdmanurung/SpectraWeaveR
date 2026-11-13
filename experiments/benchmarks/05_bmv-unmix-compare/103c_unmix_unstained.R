# Purpose: unmix the zip's unstained AF file with each tuning arm's own spectra, AF spectra and variants, so every arm gets a
#          data-driven noise floor (user 2026-11-14): per channel, q0.5 = "too negative" floor, q99.5 = positive edge.
#          The file is first cut to the zip's cell rectangle (FSC-A 6e5-4e6, SSC-A 1e5-1.6e6), as the zip cleaned its cell controls.
# Run:     SW_DATASET=24h SW_RUN_DATE=<tuning run date> SW_TUNE_SCOPE=<smoke|all> sbatch experiments/slurm/run_af_tuning.sbatch unstained
# Inputs:  experiments/runs/<SW_RUN_DATE>_bmv24h-af-tuning/<arm>/spectra.csv, af_spectra.csv, variants.rds   [103]
#          <data_dir>/autoSpectral_v1.8.8/Reference Group/unstained_AF/*.fcs, Reference Group Clean/  [zip]
# Outputs: <root>/<arm>/unstained/*.fcs        the unstained file unmixed by the arm
#          <root>/<arm>/unstained/thresholds.csv   channel, n, floor (q0.5), edge (q99.5), arcsinh cofactor 6000
#          <root>/unstained_sessionInfo.txt

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ------------------------------------------------------------------------------
CELL_FSC <- c(6e5, 4e6) # zip script cell rectangle (inherited: zip)
CELL_SSC <- c(1e5, 1.6e6)
PROBS <- c(floor = 0.005, edge = 0.995) # user 2026-11-14: q0.5 floor, q99.5 edge
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
PIPELINE <- "joint"
SPEED <- "slow"
SCOPE <- Sys.getenv("SW_TUNE_SCOPE") # smoke -> af-tuning root, all -> af-confirm root (103)
stopifnot("SW_TUNE_SCOPE is smoke or all" = SCOPE %in% c("smoke", "all"))
ROOT <- normalizePath(file.path("experiments/runs", paste0(Sys.getenv("SW_RUN_DATE"), "_bmv", cfg("run_tag"), c(smoke = "-af-tuning", all = "-af-confirm")[[SCOPE]])))

# Unstained file in the cell rectangle ----------------------------------------------------------
d188 <- file.path(normalizePath(data_path()), "autoSpectral_v1.8.8")
clean_dir <- file.path(d188, "Reference Group Clean")
control_file <- file.path(clean_dir, "fcs_control_file.csv")
unstained <- list.files(file.path(d188, "Reference Group", "unstained_AF"), pattern = "\\.fcs$", full.names = TRUE)
stopifnot("one unstained AF file" = length(unstained) == 1)
ff <- read.FCS(unstained, transformation = FALSE, truncate_max_range = FALSE)
ff <- Subset(ff, rectangleGate(`FSC-A` = CELL_FSC, `SSC-A` = CELL_SSC))
gated_dir <- file.path(tempdir(), "unstained_gated")
dir.create(gated_dir)
write.FCS(ff, file.path(gated_dir, basename(unstained)))

# Unmix per arm ---------------------------------------------------------------------------------
arms <- list.dirs(ROOT, recursive = FALSE)
arms <- arms[file.exists(file.path(arms, "variants.rds")) & !dir.exists(file.path(arms, "unstained"))]
walk(arms, \(arm_dir) {
  out <- file.path(arm_dir, "unstained")
  dir.create(file.path(out, "work"), recursive = TRUE)
  asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
  dir_fields <- names(asp)[grepl("dir", names(asp)) & map_lgl(asp, is.character)]
  asp[dir_fields] <- map(asp[dir_fields], \(d) file.path(out, "work", d))
  spectra <- as.matrix(read.csv(file.path(arm_dir, "spectra.csv"), row.names = 1, check.names = FALSE))
  af <- as.matrix(read.csv(file.path(arm_dir, "af_spectra.csv"), row.names = 1, check.names = FALSE))
  args <- readLines(file.path(arm_dir, "args.txt"))
  n_passes <- as.integer(sub(".*n.passes=n.af.passes=([0-9]+).*", "\\1", args[[2]]))
  unmix.folder(
    fcs.dir = gated_dir, spectra = spectra, asp = asp, flow.control = reload.flow.control(clean_dir, control_file, asp),
    method = "AutoSpectral", af.spectra = af, spectra.variants = readRDS(file.path(arm_dir, "variants.rds")),
    pipeline = PIPELINE, speed = SPEED, n.passes = n_passes, n.af.passes = n_passes, output.dir = out, parallel = TRUE, verbose = FALSE
  )
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
writeLines(capture.output(sessionInfo()), file.path(ROOT, "unstained_sessionInfo.txt"))
