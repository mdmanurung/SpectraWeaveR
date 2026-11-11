# experiments/integration/aurora_exvivo/01_extract.R
# Purpose: read-only pass over a folder of Aurora FCS files. Keep a random subsample of
#          events per file and the header keys needed by the diagnostics. No data is copied
#          out of the run folder and no file name is written anywhere but the run folder.
# Run:     Rscript experiments/integration/aurora_exvivo/01_extract.R <fcs_dir> <run_dir>
# Inputs:  <fcs_dir>/*.fcs  gated live single cells (marker channels already transformed)
# Outputs: <run_dir>/events.rds   named list: one numeric matrix per file (Time + marker channels)
#          <run_dir>/headers.rds  named list: header key/value vector per file
#          <run_dir>/files.csv    file, batch, sample, n_events_in_file, n_kept
#
# Settings below are illustrative choices made for a first real-data check. They are not
# validated values and not results about the donors.

suppressPackageStartupMessages({
  library(flowCore)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
fcs_dir <- args[[1]]
run_dir <- args[[2]]
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

# Settings (illustrative) ----------------------------------------------------
N_KEEP <- 5000 # events kept per file, drawn at random
SEED <- 1
N_CORES <- 8
# Columns that are not marker channels: scatter, time and the OMIQ labels.
NOT_MARKER <- "^(Time|FSC|SSC|peacoqc-passed|OmiqFilter)"
# Header keys kept for the audit (plain keys; per-parameter keys are kept too).
KEEP_KEYS <- "^(\\$(CYT|CYTSN|DATE|BTIM|ETIM|INST|OP|TOT|FIL|PAR)$|CREATOR|GUID|FILENAME|\\$P[0-9]+[NSV]$)"

files <- sort(list.files(fcs_dir, pattern = "\\.fcs$", full.names = TRUE))
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

one <- function(path) {
  hdr <- unlist(read.FCSheader(path)[[1]])
  ff <- read.FCS(path, truncate_max_range = FALSE, transformation = FALSE)
  x <- exprs(ff)
  keep <- sample.int(nrow(x), min(N_KEEP, nrow(x)))
  marker <- !grepl(NOT_MARKER, colnames(x))
  m <- x[keep, c("Time", colnames(x)[marker]), drop = FALSE]
  desc <- pData(parameters(ff))$desc
  names(desc) <- pData(parameters(ff))$name
  list(
    events = m,
    header = hdr[grepl(KEEP_KEYS, names(hdr))],
    n_in_file = nrow(x),
    marker_names = desc[colnames(m)]
  )
}

res <- mclapply(files, one, mc.cores = N_CORES, mc.set.seed = TRUE)
failed <- vapply(res, inherits, logical(1), "try-error")
if (any(failed)) {
  stop("Reading failed for ", sum(failed), " file(s).")
}
ids <- basename(files)
names(res) <- ids

saveRDS(lapply(res, `[[`, "events"), file.path(run_dir, "events.rds"))
saveRDS(lapply(res, `[[`, "header"), file.path(run_dir, "headers.rds"))
saveRDS(res[[1]]$marker_names, file.path(run_dir, "marker_names.rds"))
info <- data.frame(
  file = ids,
  batch = sub(".*(Batch[0-9]+)\\.fcs$", "\\1", ids),
  sample = sub(".*(BMV_[0-9]+[A-Z]).*", "\\1", ids),
  n_events_in_file = vapply(res, `[[`, numeric(1), "n_in_file"),
  n_kept = vapply(res, \(r) nrow(r$events), numeric(1))
)
write.csv(info, file.path(run_dir, "files.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
cat("kept", sum(info$n_kept), "events from", nrow(info), "files\n")
