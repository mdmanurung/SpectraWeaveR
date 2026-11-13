# Purpose: re-clean from scratch with the approved pooled SpectroFlo gates (plans/common_gates_plan.md, G10,
#          G10a; user 2026-11-11). Step a: gate SpectroFlo on the UNCLEANED files with the pooled template
#          (every row collapseDataForGating = TRUE, groupBy = samples per gate group, 5 in 3D) to get the live set. Step b: flag events
#          with the 12_cleanup.R rules, unchanged (PeacoQC RemoveMargins at -21,761 and 1e7 on every
#          fluorescence channel plus the scatter limits; below -NEG_T arcsinh in at least NEG_N fluorescence
#          channels), as a union over the 11 P arms, and drop a flagged event only when it is in the
#          SpectroFlo pooled live set (12 used the live events shared by all arms). Step c: write clean3
#          FCS for every arm (P arms and ablation arms), dropping the same rows from each. Before any file
#          is written, every arm must pair with SpectroFlo on Time, FSC-A, FSC-H and SSC-A on all rows.
#          One script for 3D, 24H and EXV. With gate_group set in the config (24H: condition) the pooled
#          template is fitted once per group and the GatingSets are merged; with null it is one pooled fit.
# Run:     SW_DATASET=3d SW_RUN_DATE=2026-11-11 pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/97c_reclean_pooled.R
#            from the repo root; on SLURM: sbatch experiments/slurm/run_clean3.sbatch 97c_reclean_pooled.R
#          SW_DATASET is 3d, 24h or exv and SW_RUN_DATE is YYYY-MM-DD (neither has a default).
#          SW_ABLATION_KEY (default ablation-gated) names the registry key whose manifest lists the ablation
#          arms; "none" runs the P arms only (a dataset without ablation arms yet).
# Inputs:  registry keys of datasets/<SW_DATASET>_runs.csv, files under data_path():
#          run_in("unmix-compare-final")/manifest.csv   uncleaned FCS of the P arms [01]
#          run_in("ablation-gated")/manifest.csv        ablation arm names and viability (SW_ABLATION_KEY)
#          data_path("unmix_cells", <arm>)/**/*.fcs     uncleaned FCS of the ablation arms, joined to the
#                                                       samples on the file stem
#          run_in("unmix-compare-clean2")/cleanup_counts.csv   clean2 totals, for comparison (3D only; the
#                                                       columns are left out when the key is not registered)
#          datasets/<SW_DATASET>_samples.csv, <SW_DATASET>.json   samples, channels, gate_group, lib_dirs
# Outputs: <out_dir> = run_out("clean3")
#          <out_dir>/gs/spectroflo/            pooled GatingSet on the uncleaned SpectroFlo files (all samples;
#                                              merged over the gate groups when gate_group is set)
#          <out_dir>/templates/spectroflo.csv  template used (spectroflo_<group>.csv per group when gate_group is set)
#          <out_dir>/pairing_uncleaned.csv     arm, sample, n_rows_equal, time_equal, fsca_equal, fsch_equal, ssca_equal
#          <out_dir>/manifest.csv              arm, dir, viability, file (clean3), sample, for every arm
#          <out_dir>/cleanup_counts.csv        per sample: events, SpectroFlo live, dropped by margins, by
#                                              negatives, total; clean2 totals beside them
#          <out_dir>/sessionInfo.txt
#          <data_dir>/cleaned/<basename of out_dir>/<arm>/<file>.fcs   clean3 FCS

# dataset.R first: it only loads dplyr and readr, and cfg("lib_dirs") must be known before the other libraries
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths())) # ggplot2 3.5.2 for ggcyto, cytoUtils 0.1.1 as 19a
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
OUT_DIR <- run_out("clean3")
CLEAN_DIR <- data_path("cleaned", basename(OUT_DIR))
P_MANIFEST <- file.path(run_in("unmix-compare-final"), "manifest.csv") # uncleaned files of the P arms
A_KEY <- Sys.getenv("SW_ABLATION_KEY", "ablation-gated") # registry key of the ablation manifest; "none" = P arms only
A_MANIFEST <- if (A_KEY == "none") {
  NA_character_
} else {
  file.path(run_in(A_KEY), "manifest.csv") # ablation arms; files are clean2
}
HAS_CLEAN2 <- "unmix-compare-clean2" %in% .sw05_runs$key # a clean2 run exists for 3D only
CLEAN2_COUNTS <- if (HAS_CLEAN2) {
  file.path(run_in("unmix-compare-clean2"), "cleanup_counts.csv")
}
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
CD4 <- cfg("channels")$CD4
CD8 <- cfg("channels")$CD8
CD3 <- cfg("channels")$CD3
TNF <- cfg("channels")$TNF_slot
IFNG <- cfg("channels")$IFNg_slot
VIABILITY_CHANNELS <- c(
  cfg("channels")$viability_spectroflo,
  cfg("channels")$viability_autospectral
) # left out of the negative-channel flag
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
# Ablation files are found by the file stem under unmix_cells/<arm>/ (24H: one subfolder per condition).
smp <- samples_tbl()
p_src <- read_csv(P_MANIFEST, show_col_types = FALSE) |>
  arrange(arm, sample)
if (is.na(A_MANIFEST)) {
  message(
    "97c_reclean_pooled.R: SW_ABLATION_KEY=none, the ablation arms are not cleaned"
  )
  a_src <- p_src[0, ]
} else {
  a_src <- read_csv(A_MANIFEST, show_col_types = FALSE) |>
    dplyr::filter(!arm %in% p_src$arm) |>
    distinct(arm, dir, viability) |>
    reframe(
      file = list.files(
        data_path("unmix_cells", arm),
        pattern = "\\.fcs$",
        full.names = TRUE,
        recursive = TRUE
      ),
      .by = c(arm, dir, viability)
    ) |>
    mutate(
      stem = basename(file) |>
        stringr::str_remove("\\.fcs$") |>
        stringr::str_remove("[ _]*AutoSpectral.*$")
    ) |>
    dplyr::filter(stem %in% smp$stem) |> # unmix_cells/<arm>/ also holds controls/ and objects/
    split(~arm) |>
    map(\(m) {
      inner_join(
        m,
        select(smp, sample, stem),
        by = join_by(stem),
        relationship = "one-to-one",
        unmatched = "error"
      )
    }) |>
    list_rbind() |>
    select(arm, dir, viability, file, sample) |>
    arrange(arm, sample)
}
sources <- bind_rows(p_src, a_src)
stopifnot(
  "one file per sample and arm" = all(count(sources, arm)$n == nrow(smp)),
  "manifest samples are the sample table" = setequal(
    sources$sample,
    smp$sample
  ),
  "every source file exists" = all(file.exists(sources$file))
)
donors <- sort(unique(sources$sample))
ref_src <- dplyr::filter(sources, arm == REF_ARM) |>
  arrange(sample)

# a. Pooled SpectroFlo gates on the uncleaned files -------------------------------------------------
# One pooled fit per gate group (24H: condition, user 2026-11-12); gate_group null = one fit over all samples.
VIABILITY <- ref_src$viability[[1]]
gate_group_of <- if (is.null(DS$gate_group)) {
  rep("all", nrow(smp))
} else {
  smp[[DS$gate_group]]
}
ref_src <- ref_src |>
  mutate(gate_group = gate_group_of[match(sample, smp$sample)])
# Viability search window: mindensity2 needs one (user 2026-11-12, page Q17); NULL keeps the 3D template
VIABILITY_ARGS <- if (is.null(cfg("viability_range"))) NA else paste0("gate_range=", cfg("viability_range"))
template_base <- tibble(
  alias = c("nonDebris", "singlets", "live", "CD3pos", "*"),
  pop = c("+", "+", "-", "+", "+/-+/-"),
  parent = c("root", "nonDebris", "singlets", "live", "CD3pos"),
  dims = c("FSC-A", "FSC-A,FSC-H", VIABILITY, CD3, paste(CD4, CD8, sep = ",")),
  gating_method = c(
    cfg("mindensity_method"),
    "singletGate",
    cfg("mindensity_method"),
    cfg("mindensity_method"),
    "quadGate.tmix"
  ),
  gating_args = c(
    paste0("gate_range=", DEBRIS_RANGE),
    "wider_gate=TRUE",
    VIABILITY_ARGS,
    paste0("gate_range=", LINEAGE_RANGE),
    CD4_CD8_TMIX_ARGS
  ),
  collapseDataForGating = TRUE,
  preprocessing_method = NA,
  preprocessing_args = NA
)
dir.create(file.path(OUT_DIR, "templates"))
# groupBy = number of samples in the group: openCyto pools only when groupBy is not empty
fit_group <- function(src, group) {
  template_file <- file.path(
    OUT_DIR,
    "templates",
    paste0(
      REF_ARM,
      if (is.null(DS$gate_group)) "" else paste0("_", group),
      ".csv"
    )
  )
  template_base |>
    mutate(groupBy = nrow(src), .after = collapseDataForGating) |>
    write_csv(template_file)
  gs <- sw_gate_run(
    fcs_files = src$file,
    gating_template = template_file,
    transform_channels = c(VIABILITY, CD3, CD4, CD8, TNF, IFNG),
    cofactor = COFACTOR
  )
  sampleNames(gs) <- src$sample
  gs
}
ref_fits <- ref_src |>
  split(~gate_group) |>
  imap(fit_group)
ref_gs <- if (length(ref_fits) == 1) {
  ref_fits[[1]]
} else {
  merge_list_to_gs(unname(ref_fits)) # one GatingSet over all samples, each with its group's gates
}
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
      !grepl("^(FSC|SSC)|^AF", colnames(e)) &
      !colnames(e) %in% VIABILITY_CHANNELS
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
  list_rbind()
if (HAS_CLEAN2) {
  counts <- counts |>
    left_join(
      read_csv(CLEAN2_COUNTS, show_col_types = FALSE) |>
        select(
          sample,
          n_events_clean2 = n_events,
          n_shared_live_clean2 = n_shared_live,
          n_margin_clean2 = n_margin,
          n_neg_clean2 = n_neg,
          n_dropped_clean2 = n_dropped
        ),
      by = join_by(sample),
      relationship = "one-to-one",
      unmatched = "error"
    )
} else {
  message(
    "97c_reclean_pooled.R: no clean2 run registered for dataset ",
    DS$dataset,
    ": cleanup_counts.csv has no clean2 columns"
  )
}
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
