# Purpose: proof chain step 5 (user, 2026-11-10): spike a known BV510 (TNF) signal into real raw events.
#          For CD4T and CD8T cells (as188_joint gates), add k x the 1.8.8 BV510 spectrum to the raw detector
#          values, with k = 6000 x sinh(level) for levels 1, 2, 3 (arcsinh units, cofactor 6000). Spectra rows
#          have max 1, so an exact unmix raises BV510 by k. Per donor, T cells are shuffled (set.seed(42),
#          user-approved) and split in quarters: unspiked, level 1, level 2, level 3. Other events untouched.
# Run:     sbatch experiments/slurm/run_af405_proof.sbatch 62 (or pixi run Rscript <this script>)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Raw/*.fcs            raw files [user zip]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed  row keys [user zip]
#          experiments/runs/2026-11-10_bmv-gated/manifest.csv, gs/as188_joint             [12_cleanup.R, 19a]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/..._143803.csv
# Outputs: experiments/bmv_test/compare_3D/data/spiked/2026-11-10_af405_proof/*.fcs  spiked raw files
#          experiments/runs/2026-11-10_bmv-af405-proof/spike_key.csv   sample, raw_row, pop, level, k
#          experiments/runs/2026-11-10_bmv-af405-proof/sessionInfo_62.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters ------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
RUN_DIR <- "experiments/runs/2026-11-10_bmv-af405-proof"
GATED <- "experiments/runs/2026-11-10_bmv-gated" # the gates 39b/45b used
OUT_FCS <- file.path(DATA, "spiked", "2026-11-10_af405_proof")
LEVELS <- c(0, 1, 2, 3) # arcsinh units added; 0 = unspiked control quarter (user-approved grid)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
SEED <- 42 # user-approved, which cells get which level
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; maps raw rows to cleaned rows
POPS <- c("CD4T", "CD8T")

stopifnot("spiked fcs dir must not exist" = !dir.exists(OUT_FCS))
dir.create(OUT_FCS, recursive = TRUE)

# BV510 spectrum on the detector channels -------------------------------------------------
s188 <- as.matrix(read.csv(
  file.path(
    DATA,
    "autoSpectral_v1.8.8/table_spectra/Automated_autospectral_spectra_20261105_143803.csv"
  ),
  row.names = 1,
  check.names = FALSE
))
b510 <- s188["BV510", ]

# Raw rows of CD4T and CD8T per donor ---------------------------------------------------------
manifest <- read_csv(
  file.path(GATED, "manifest.csv"),
  show_col_types = FALSE
) |>
  filter(arm == "as188_joint")
read_x <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
}
key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))
gs <- load_gs(file.path(GATED, "gs/as188_joint"))
raw_unmixed_dir <- file.path(DATA, "autoSpectral_v1.8.8/AutoSpectral_unmixed")
raw_dir <- file.path(DATA, "autoSpectral_v1.8.8/Raw")

# Spike each donor's raw file and record the key ------------------------------------------------
set.seed(SEED)
spike_key <- map2(manifest$sample, manifest$file, \(s, clean_file) {
  unmixed_name <- basename(clean_file)
  raw_name <- sub(" AutoSpectral\\.fcs$", ".fcs", unmixed_name)
  unmixed_key <- key_of(read_x(file.path(raw_unmixed_dir, unmixed_name)))
  keep <- unmixed_key %in% key_of(read_x(clean_file))
  kept_rows <- which(keep)
  ff <- read.FCS(
    file.path(raw_dir, raw_name),
    transformation = FALSE,
    truncate_max_range = FALSE
  )
  stopifnot(
    "raw rows equal unmixed rows" = identical(key_of(exprs(ff)), unmixed_key)
  )
  pop_rows <- map(set_names(POPS), \(p) {
    kept_rows[gh_pop_get_indices(gs[[s]], p)]
  })
  key <- imap(pop_rows, \(r, p) tibble(raw_row = r, pop = p)) |>
    list_rbind() |>
    distinct(raw_row, .keep_all = TRUE) |> # CD4T and CD8T are disjoint quadrants; keep first if not
    mutate(
      level = sample(rep_len(LEVELS, n())),
      k = COFACTOR * sinh(level),
      sample = s,
      raw_file = raw_name
    )
  x <- exprs(ff)
  x[key$raw_row, colnames(s188)] <- x[key$raw_row, colnames(s188)] +
    outer(key$k, b510)
  exprs(ff) <- x
  write.FCS(ff, file.path(OUT_FCS, raw_name))
  key
}) |>
  list_rbind()
write_csv(spike_key, file.path(RUN_DIR, "spike_key.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_62.txt")
)
