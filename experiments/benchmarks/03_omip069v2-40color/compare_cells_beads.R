#!/usr/bin/env Rscript
# Purpose: compare reference spectra derived from cell controls with those from bead controls (OMIP-069v2).
# Run:     Rscript experiments/benchmarks/03_omip069v2-40color/compare_cells_beads.R <cells_run_dir> <beads_run_dir> <out_csv>
# Inputs:  <run_dir>/report/spectra.csv of the two runs (long form: fluorophore, detector, value)
# Outputs: <out_csv>: per fluorophore, angle between the cell and the bead spectrum
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
read_spectra <- function(run_dir) {
  long <- read.csv(file.path(run_dir, "report", "spectra.csv"), check.names = FALSE)
  fl <- unique(long$fluorophore)
  dt <- unique(long$detector)
  matrix(long$value[match(paste(rep(fl, times = length(dt)), rep(dt, each = length(fl))), paste(long$fluorophore, long$detector))],
         nrow = length(fl), dimnames = list(fl, dt))
}
cells <- read_spectra(args[[1]])
beads <- read_spectra(args[[2]])
shared <- intersect(rownames(cells), rownames(beads))
stopifnot("no shared fluorophores" = length(shared) > 0)
out <- sw_unmix_compare_spectra(cells[shared, , drop = FALSE], beads[shared, colnames(cells), drop = FALSE])
stopifnot("new file" = !file.exists(args[[3]]))
write.csv(out, args[[3]], row.names = FALSE)
cat(length(shared), "fluorophores compared;", length(setdiff(rownames(cells), rownames(beads))), "only in cells\n")
