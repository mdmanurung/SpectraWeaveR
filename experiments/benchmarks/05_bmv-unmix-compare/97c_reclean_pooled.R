# Purpose: re-clean from scratch with the approved pooled SpectroFlo gates (plans/common_gates_plan.md, G10,
#          G10a; user 2026-11-11). Step a: gate SpectroFlo on the UNCLEANED files with the pooled template
#          (every row collapseDataForGating = TRUE, groupBy = 5) to get the live set. Step b: flag events
#          with the 12_cleanup.R rules, unchanged (PeacoQC RemoveMargins at -21,761 and 1e7 on every
#          fluorescence channel plus the scatter limits; below -NEG_T arcsinh in at least NEG_N fluorescence
#          channels), as a union over the 11 P arms, and drop a flagged event only when it is in the
#          SpectroFlo pooled live set (12 used the live events shared by all arms). Step c: write clean3
#          FCS for every arm (P arms and ablation arms), dropping the same rows from each. Before any file
#          is written, every arm must pair with SpectroFlo on Time, FSC-A, FSC-H and SSC-A on all rows.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/97c_reclean_pooled.R
#            from the repo root; on SLURM: sbatch experiments/slurm/run_clean3.sbatch 97c_reclean_pooled.R
# Inputs:  experiments/runs/2026-11-09_bmv-unmix-compare-final/manifest.csv   uncleaned FCS of the 11 P arms [01]
#          experiments/runs/2026-11-10_bmv-ablation-gated/manifest.csv        ablation arm names and viability
#          experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/*.fcs      uncleaned FCS of the ablation arms
#          experiments/runs/2026-11-09_bmv-unmix-compare-clean2/cleanup_counts.csv   clean2 totals, for comparison
# Outputs: <out_dir>/gs/spectroflo/            pooled GatingSet on the uncleaned SpectroFlo files
#          <out_dir>/templates/spectroflo.csv  template used
#          <out_dir>/pairing_uncleaned.csv     arm, sample, n_rows_equal, time_equal, fsca_equal, fsch_equal, ssca_equal
#          <out_dir>/manifest.csv              arm, dir, viability, file (clean3), sample, for every arm
#          <out_dir>/cleanup_counts.csv        per sample: events, SpectroFlo live, dropped by margins, by
#                                              negatives, total; clean2 totals beside them
#          <out_dir>/sessionInfo.txt
#          experiments/bmv_test/compare_3D/data/cleaned/2026-11-11_bmv-clean3/<arm>/<file>.fcs   clean3 FCS

.libPaths(c(
  "experiments/bmv_test/compare_3D/data/lib_ggplot2_35",
  "experiments/bmv_test/compare_3D/data/lib_cytoutils", # cytoUtils 0.1.1 as 19a
  .libPaths()
))
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(cytoUtils)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE) # sw_gate_run, sw_filter_margins_peacoqc

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- "experiments/runs/2026-11-11_bmv-clean3"
CLEAN_DIR <- "experiments/bmv_test/compare_3D/data/cleaned/2026-11-11_bmv-clean3"
DATA <- "experiments/bmv_test/compare_3D/data"
P_MANIFEST <- "experiments/runs/2026-11-09_bmv-unmix-compare-final/manifest.csv" # uncleaned files of the 11 P arms
A_MANIFEST <- "experiments/runs/2026-11-10_bmv-ablation-gated/manifest.csv" # ablation arms; files are clean2
CLEAN2_COUNTS <- "experiments/runs/2026-11-09_bmv-unmix-compare-clean2/cleanup_counts.csv"
REF_ARM <- "spectroflo"
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")) # FCS reads and writes run in parallel
# rules inherited unchanged from 12_cleanup.R
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
NEG_T <- 2 # below -T in arcsinh units (user, 2026-11-09)
NEG_N <- 3 # in at least N fluorescence channels (user, 2026-11-09)
MARGIN_SPECS <- list(AllFluoChannels = c(-21761, 1e7)) # PeacoQC RemoveMargins, as 12
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; finds the rows PeacoQC kept
MAX_MARGIN_FRAC <- 0.5 # stop if margins flag more than half the SpectroFlo live events
# pooled template, as 97b_spectroflo_labels_pooled.R (approved, G10a)
DEBRIS_RANGE <- "c(2e5, 1.5e6)"
LINEAGE_RANGE <- "c(0, 3)"
CD4_CD8_TMIX_ARGS <- "K=3, usePrior='no', quantile1=0.8, quantile3=0.8"
CD4 <- "Spark Blue 574-A"
CD8 <- "SBUV445-A"
CD3 <- "SBUV510-A"
TNF <- "BV510-A"
IFNG <- "BV570-A"
QUADRANTS <- c(
  CD4T = paste0(CD4, "+", CD8, "-"),
  CD8T = paste0(CD4, "-", CD8, "+"),
  DP = paste0(CD4, "+", CD8, "+"),
  DN = paste0(CD4, "-", CD8, "-")
)
LIVE_NODE <- "/nonDebris/singlets/live"

stopifnot("run dir is new" = !dir.exists(OUT_DIR))
stopifnot("clean3 dir is new" = !dir.exists(CLEAN_DIR))
dir.create(OUT_DIR, recursive = TRUE)

# Run a function over a vector in parallel; any worker error stops the script (mclapply would hide it)
par_map <- \(x, f) {
  out <- parallel::mclapply(x, f, mc.cores = N_CORES)
  failed <- map_lgl(out, \(o) inherits(o, "try-error") || is.null(o))
  if (any(failed)) {
    stop(
      "worker failed: ",
      paste(out[failed][[1]], collapse = " "),
      call. = FALSE
    )
  }
  out
}
read_x <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
}
key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))

# Uncleaned source of every arm -------------------------------------------------------------------
# P arms: the files 12_cleanup.R read. Ablation arms: unmix_cells/<arm>/<name of as188_joint's file>, the
# raw files 51_score_ablation.R masks. The three base arms in the ablation manifest are P arms and not repeated.
p_src <- read_csv(P_MANIFEST, show_col_types = FALSE) |>
  arrange(arm, sample)
a_man <- read_csv(A_MANIFEST, show_col_types = FALSE)
a_src <- a_man |>
  dplyr::filter(!arm %in% p_src$arm) |>
  mutate(file = file.path(DATA, "unmix_cells", arm, basename(file))) |>
  arrange(arm, sample)
sources <- bind_rows(p_src, a_src)
stopifnot(
  "5 files per arm" = all(count(sources, arm)$n == 5),
  "every source file exists" = all(file.exists(sources$file))
)
donors <- sort(unique(sources$sample))
ref_src <- dplyr::filter(sources, arm == REF_ARM) |>
  arrange(sample)

# a. Pooled SpectroFlo gates on the uncleaned files -------------------------------------------------
VIABILITY <- ref_src$viability[[1]]
template <- tibble(
  alias = c("nonDebris", "singlets", "live", "CD3pos", "*"),
  pop = c("+", "+", "-", "+", "+/-+/-"),
  parent = c("root", "nonDebris", "singlets", "live", "CD3pos"),
  dims = c("FSC-A", "FSC-A,FSC-H", VIABILITY, CD3, paste(CD4, CD8, sep = ",")),
  gating_method = c(
    "mindensity",
    "singletGate",
    "mindensity",
    "mindensity",
    "quadGate.tmix"
  ),
  gating_args = c(
    paste0("gate_range=", DEBRIS_RANGE),
    "wider_gate=TRUE",
    NA,
    paste0("gate_range=", LINEAGE_RANGE),
    CD4_CD8_TMIX_ARGS
  ),
  collapseDataForGating = TRUE,
  groupBy = nrow(ref_src),
  preprocessing_method = NA,
  preprocessing_args = NA
)
dir.create(file.path(OUT_DIR, "templates"))
template_file <- file.path(OUT_DIR, "templates", paste0(REF_ARM, ".csv"))
write_csv(template, template_file)
ref_gs <- sw_gate_run(
  fcs_files = ref_src$file,
  gating_template = template_file,
  transform_channels = c(VIABILITY, CD3, CD4, CD8, TNF, IFNG),
  cofactor = COFACTOR
)
sampleNames(ref_gs) <- ref_src$sample
iwalk(QUADRANTS, \(node, alias) {
  gs_pop_set_name(ref_gs, paste0(LIVE_NODE, "/CD3pos/", node), alias)
})
dir.create(file.path(OUT_DIR, "gs"))
save_gs(ref_gs, file.path(OUT_DIR, "gs", REF_ARM))
live <- map(set_names(donors), \(d) {
  gh_pop_get_indices(ref_gs[[d]], LIVE_NODE)
})

# Pairing of every arm with SpectroFlo on all uncleaned rows; stop before writing anything --------------
ref_pair <- map(set_names(donors), \(d) {
  read_x(ref_src$file[ref_src$sample == d])[, KEY_COLS, drop = FALSE]
})
stopifnot(
  "GatingSet rows are the SpectroFlo rows" = all(map_lgl(donors, \(d) {
    length(live[[d]]) == nrow(ref_pair[[d]])
  }))
)
pairing <- par_map(seq_len(nrow(sources)), \(i) {
  m <- read_x(sources$file[[i]])
  ref <- ref_pair[[sources$sample[[i]]]]
  same_n <- nrow(m) == nrow(ref)
  same_col <- \(col) same_n && identical(unname(m[, col]), unname(ref[, col]))
  tibble(
    arm = sources$arm[[i]],
    sample = sources$sample[[i]],
    n_rows_equal = same_n,
    time_equal = same_col("Time"),
    fsca_equal = same_col("FSC-A"),
    fsch_equal = same_col("FSC-H"),
    ssca_equal = same_col("SSC-A")
  )
}) |>
  list_rbind()
write_csv(pairing, file.path(OUT_DIR, "pairing_uncleaned.csv"))
stopifnot(
  "every arm pairs with SpectroFlo on Time, FSC-A, FSC-H, SSC-A (see pairing_uncleaned.csv)" = all(
    pairing$n_rows_equal &
      pairing$time_equal &
      pairing$fsca_equal &
      pairing$fsch_equal &
      pairing$ssca_equal
  )
)

# b. Flags as in 12_cleanup.R, union over the P arms ------------------------------------------------
flag_file <- \(f) {
  ff <- sw_io_read_fcs(f)
  e <- exprs(ff)
  key <- key_of(e)
  kept <- sw_filter_margins_peacoqc(
    ff,
    channel_specifications = MARGIN_SPECS
  ) |>
    exprs()
  ch <- colnames(e)[
    grepl("-A$", colnames(e)) &
      !grepl("^(FSC|SSC)|Zombie|LIVE DEAD|^AF", colnames(e))
  ]
  list(
    margin = !(key %in% key_of(kept)),
    neg = rowSums(asinh(e[, ch, drop = FALSE] / COFACTOR) < -NEG_T) >= NEG_N
  )
}
flags <- par_map(p_src$file, flag_file)
drops <- map(set_names(donors), \(d) {
  f <- flags[p_src$sample == d]
  margin <- reduce(map(f, "margin"), `|`) & live[[d]]
  neg <- reduce(map(f, "neg"), `|`) & live[[d]]
  list(margin = margin, neg = neg, drop = margin | neg)
})
clean2 <- read_csv(CLEAN2_COUNTS, show_col_types = FALSE) |>
  select(
    sample,
    n_events_clean2 = n_events,
    n_shared_live_clean2 = n_shared_live,
    n_margin_clean2 = n_margin,
    n_neg_clean2 = n_neg,
    n_dropped_clean2 = n_dropped
  )
counts <- imap(drops, \(d, s) {
  tibble(
    sample = s,
    n_events = length(d$drop),
    n_spectroflo_live = sum(live[[s]]),
    n_margin = sum(d$margin),
    n_neg = sum(d$neg),
    n_dropped = sum(d$drop)
  )
}) |>
  list_rbind() |>
  left_join(
    clean2,
    by = join_by(sample),
    relationship = "one-to-one",
    unmatched = "error"
  )
write_csv(counts, file.path(OUT_DIR, "cleanup_counts.csv"))
stopifnot(
  "margins flag more than MAX_MARGIN_FRAC of SpectroFlo live events" = all(
    counts$n_margin / counts$n_spectroflo_live < MAX_MARGIN_FRAC
  )
)

# c. Write clean3 FCS for every arm, the same rows dropped (arms pair row by row, checked above) ----------
manifest_clean3 <- sources |>
  mutate(clean_file = file.path(CLEAN_DIR, arm, basename(file)))
invisible(par_map(seq_len(nrow(manifest_clean3)), \(i) {
  out <- manifest_clean3$clean_file[[i]]
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  ff <- sw_io_read_fcs(manifest_clean3$file[[i]])
  write.FCS(ff[!drops[[manifest_clean3$sample[[i]]]]$drop, ], out)
  TRUE
}))
manifest_clean3 |>
  mutate(file = clean_file) |>
  select(arm, dir, viability, file, sample) |>
  write_csv(file.path(OUT_DIR, "manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
