# experiments/integration/cytobatchflagr/01_extract.R
# Purpose: read-only pass over the CytoBatchFlagR files (Zenodo 15388817, CC BY 4.0): 18 batches,
#          each with the same four control samples (tcu1, tcu2 unstimulated; tcs3, tcs4
#          stimulated), Aurora, linear unmixed values, pre-gated CD45+CD14- live single cells.
#          Keep a random subsample of events per file, the batch table and acquisition dates.
# Run:     Rscript experiments/integration/cytobatchflagr/01_extract.R <data_dir> <run_dir>
# Inputs:  <data_dir>/*.fcs, <data_dir>/Table_2.csv   from download_zenodo.sh cytobatchflagr
# Outputs: <run_dir>/events.rds        named list: one numeric matrix per file (Time + all channels)
#          <run_dir>/files.csv         file, batch, control, n_events_in_file, n_kept, date
#          <run_dir>/marker_names.rds  channel -> marker name
#          <run_dir>/sessionInfo.txt
#
# Settings are illustrative choices for a first check; not validated values.

suppressPackageStartupMessages({
  library(flowCore)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[[1]]
run_dir <- args[[2]]
stopifnot("run dir must be new" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)

# Settings (illustrative) ----------------------------------------------------
N_KEEP <- 5000 # events kept per file, drawn at random (some files hold fewer)
SEED <- 1
N_CORES <- 8

table2 <- read.csv(file.path(data_dir, "Table_2.csv"))
files <- file.path(data_dir, table2$file)
stopifnot("every listed file exists" = all(file.exists(files)))
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

# Read --------------------------------------------------------------------------
one <- function(path) {
  ff <- read.FCS(path, truncate_max_range = FALSE, transformation = FALSE)
  x <- exprs(ff)
  keep <- sample.int(nrow(x), min(N_KEEP, nrow(x)))
  desc <- pData(parameters(ff))$desc
  names(desc) <- pData(parameters(ff))$name
  list(
    events = x[keep, , drop = FALSE],
    n_in_file = nrow(x),
    date = keyword(ff)[["$DATE"]],
    marker_names = desc
  )
}
res <- mclapply(files, one, mc.cores = N_CORES, mc.set.seed = TRUE)
failed <- vapply(res, inherits, logical(1), "try-error")
if (any(failed)) {
  stop("Reading failed for ", sum(failed), " file(s).")
}
names(res) <- table2$file

# Write -------------------------------------------------------------------------
saveRDS(lapply(res, `[[`, "events"), file.path(run_dir, "events.rds"))
saveRDS(res[[1]]$marker_names, file.path(run_dir, "marker_names.rds"))
info <- data.frame(
  file = table2$file,
  batch = table2$batch,
  control = table2$control,
  n_events_in_file = vapply(res, `[[`, numeric(1), "n_in_file"),
  n_kept = vapply(res, \(r) nrow(r$events), numeric(1)),
  date = vapply(res, `[[`, character(1), "date")
)
write.csv(info, file.path(run_dir, "files.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
cat("kept", sum(info$n_kept), "events from", nrow(info), "files; fewest kept:", min(info$n_kept), "\n")
