# Purpose: copy of 92_umap_knn.R on the common SpectroFlo gates (plans/common_gates_plan.md, G2, G3). Does the
#          unmixing arm change the UMAP of CD3+ T cells? Grill 2026-11-10 (README row "UMAP kNN test"). The
#          same 10,000 CD3+ cells per donor (common CD3pos rows, seed 42; same number and seed as 92) are
#          embedded in every arm with UMAP; per cell, the share of its 30 nearest embedding neighbours that
#          are also its neighbours in the reference embedding; median per donor. References: 1.8.8 zip and
#          SpectroFlo. Noise floor: 1.8.8 zip re-embedded with seeds 1-3 vs seed 42. Two marker sets
#          (lineage; lineage + functional) and two scalings (z-score; none). Descriptive, no threshold.
#          The script stops unless FSC-A and SSC-A of every root row are identical in every arm (G8).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_script_lib.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/92_umap_knn_common.R [<run dir>]
#            (repo root; run dir must not exist; without the argument it is run_out("common-umap-knn")).
#            Figures: 92b_plot_umap_knn_common.R.
# Inputs:  datasets/<SW_DATASET>.json, <SW_DATASET>_runs.csv, samples_tbl(); run folders come from these registry keys
#          run_in("common-gated")/gs/{as188_joint,opt188_pregate_joint_p1af1,
#            af405_afrefine,as102_var,spectroflo}   common-gate GatingSets [batch 0], README.md
#          data_path("autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv")
#            dye to marker map [user zip]
#          cfg("umap_markers")   lineage, functional and colour markers of the UMAP
# Outputs: <run_dir>/cells.csv          donor, root row of each sampled CD3+ cell
#          <run_dir>/pairing.csv        arm, donor, root rows, FSC-A and SSC-A equal to 1.8.8 zip
#          <run_dir>/embeddings.rds     arm, marker set, scaling, seed, donor, cell, UMAP 1-2, colour markers (3D: TNF, Th2)
#          <run_dir>/knn_overlap.csv    comparison, reference, marker set, scaling, donor, median and mean
#                                       per-cell neighbour overlap (donor = sample id; with conditions, also
#                                       donor_id and condition)
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE) # sw_dimred_run; run from the repo root
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# dataset.R calls a bare filter(); here flowCore's filter generic comes first on the search path, so give
# its two helpers a scope where filter is dplyr's (dataset.R should say dplyr::filter itself)
.dataset_scope <- new.env(parent = globalenv())
.dataset_scope$filter <- dplyr::filter
environment(run_in) <- .dataset_scope
environment(cut_of) <- .dataset_scope

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(RUN_DIR)) {
  RUN_DIR <- run_out("common-umap-knn")
}
COMMON <- run_in("common-gated") # batch 0; gs/<arm> for all arms
# grill Q1 (user, 2026-11-10): reference first; SpectroFlo pairs row for row with the AutoSpectral arms
ARMS <- file.path(
  COMMON,
  "gs",
  c(
    "1.8.8 zip" = "as188_joint",
    opt188 = "opt188_pregate_joint_p1af1",
    "refined AF" = "af405_afrefine",
    "1.0.2 zip" = "as102_var",
    SpectroFlo = "spectroflo"
  )
) |>
  set_names(c("1.8.8 zip", "opt188", "refined AF", "1.0.2 zip", "SpectroFlo"))
REFERENCES <- c("1.8.8 zip", "SpectroFlo") # grill Q11
POP <- "CD3pos" # grill Q5, Q9: 1.8.8 zip's CD3+ rows in every arm
N_PER_DONOR <- 10000 # grill Q5
SUBSAMPLE_SEED <- 42 # grill Q10
UMAP_SEED <- 42 # grill Q10, every arm
NOISE_SEEDS <- c(1, 2, 3) # grill Q10: 1.8.8 zip re-embedded for the noise floor, as in script 89
K <- 30 # grill Q4: neighbours per cell for the overlap
N_NEIGHBORS <- 15 # grill Q6, inherited sw_dimred_run default
MIN_DIST <- 0.2 # grill Q6, inherited sw_dimred_run default
SCALINGS <- c(zscore = TRUE, none = FALSE) # grill Q8 ("Both")
COFACTOR <- 6000 # README "Transform": arcsinh, cofactor 6000, all fluorescence channels
# 19a_gate_arm.R:113 writes these channels as arcsinh already; all other fluorescence channels are linear
STORED_ARCSINH <- unlist(
  cfg("channels")[c("CD3", "CD4", "CD8", "TNF_slot", "IFNg_slot")],
  use.names = FALSE
)
# marker panel of the UMAP: a panel list is the user's choice per dataset (config), not a default here
UMAP_MARKERS <- cfg("umap_markers")
LINEAGE <- UMAP_MARKERS$lineage # README "Cluster markers": lineage only (user, 2026-11-09)
FUNCTIONAL <- UMAP_MARKERS$functional # grill Q7: every other panel marker; viability and PE (no marker) left out
COLOUR_MARKERS <- UMAP_MARKERS$colour # markers the embeddings carry for the colour figures

MARKER_SETS <- list(
  lineage = LINEAGE,
  lineage_functional = c(LINEAGE, FUNCTIONAL)
) # grill Q3
CONTROL_FILE <- data_path(
  "autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"
)

stopifnot(
  "batch 0 finished (README.md in the common dir)" = file.exists(file.path(
    COMMON,
    "README.md"
  )),
  "run dir must not exist" = !dir.exists(RUN_DIR)
)
dir.create(RUN_DIR, recursive = TRUE)
smp <- samples_tbl()

# Marker to channel ----------------------------------------------------------------------------------
channel_of <- read_csv(CONTROL_FILE, show_col_types = FALSE) |>
  dplyr::filter(marker %in% c(LINEAGE, FUNCTIONAL)) |>
  with(set_names(paste0(fluorophore, "-A"), marker))
if (!all(c(LINEAGE, FUNCTIONAL) %in% names(channel_of))) {
  stop(
    "92_umap_knn_common.R: markers missing in the control file of dataset ",
    DS$dataset,
    ": ",
    paste(setdiff(c(LINEAGE, FUNCTIONAL), names(channel_of)), collapse = ", "),
    call. = FALSE
  )
}
channel_of <- channel_of[c(LINEAGE, FUNCTIONAL)]

# Sample the same CD3+ cells per donor from the common CD3pos label ----------------------------------------
gs_list <- map(ARMS, load_gs)
donors <- sampleNames(gs_list[["1.8.8 zip"]])
stopifnot(
  "gating set samples are the sample table" = setequal(donors, smp$sample)
)
set.seed(SUBSAMPLE_SEED)
cells <- map(set_names(donors), \(d) {
  cd3_rows <- which(gh_pop_get_indices(gs_list[["1.8.8 zip"]][[d]], POP))
  tibble(row = sort(sample(cd3_rows, N_PER_DONOR)))
}) |>
  list_rbind(names_to = "donor")
write_csv(cells, file.path(RUN_DIR, "cells.csv"))

# Marker matrix per arm on arcsinh; rows paired across arms --------------------------------------------
arm_matrix <- \(gs) {
  map(donors, \(d) {
    rows <- cells$row[cells$donor == d]
    exprs(gh_pop_get_data(gs[[d]], "root"))[
      rows,
      unname(channel_of),
      drop = FALSE
    ]
  }) |>
    reduce(rbind)
}
mats <- map(gs_list, arm_matrix)
# scatter does not depend on unmixing, so equal FSC-A and SSC-A on every root row show the arms hold the same cells (G8)
root_scatter <- \(gs, d) {
  exprs(gh_pop_get_data(gs[[d]], "root"))[, c("FSC-A", "SSC-A")]
}
pairing <- gs_list |>
  imap(\(gs, arm) {
    tibble(
      donor = donors,
      n_rows = map_int(donors, \(d) nrow(root_scatter(gs, d))),
      paired = map_lgl(donors, \(d) {
        identical(root_scatter(gs, d), root_scatter(gs_list[["1.8.8 zip"]], d))
      })
    )
  }) |>
  list_rbind(names_to = "arm")
write_csv(pairing, file.path(RUN_DIR, "pairing.csv"))
stopifnot("rows pair across arms" = all(pairing$paired))
mats <- map(mats, \(m) {
  linear <- setdiff(colnames(m), STORED_ARCSINH)
  m[, linear] <- asinh(m[, linear] / COFACTOR)
  colnames(m) <- names(channel_of)
  as_tibble(m) |>
    mutate(donor = cells$donor, cell = seq_len(n()))
})

# UMAP: every arm at seed 42, 1.8.8 zip also at the noise seeds --------------------------------------
runs <- tidyr::expand_grid(
  arm = names(ARMS),
  marker_set = names(MARKER_SETS),
  scaling = names(SCALINGS),
  seed = UMAP_SEED
) |>
  bind_rows(tidyr::expand_grid(
    arm = "1.8.8 zip",
    marker_set = names(MARKER_SETS),
    scaling = names(SCALINGS),
    seed = NOISE_SEEDS
  ))
embeddings <- runs |>
  mutate(
    umap = pmap(list(arm, marker_set, scaling, seed), \(a, ms, sc, sd) {
      sw_dimred_run(
        mats[[a]],
        markers = MARKER_SETS[[ms]],
        method = "umap",
        n_neighbors = N_NEIGHBORS,
        min_dist = MIN_DIST,
        max_cells = NULL,
        scale = SCALINGS[[sc]],
        seed = sd,
        n_threads = 8
      ) |>
        select(dim1, dim2) |>
        # no subsampling inside sw_dimred_run, so rows stay in mats order
        bind_cols(select(mats[[a]], donor, cell, all_of(COLOUR_MARKERS)))
    })
  )
saveRDS(embeddings, file.path(RUN_DIR, "embeddings.rds"))

# Neighbour overlap per cell, median per donor --------------------------------------------------------
neighbours <- embeddings |>
  mutate(
    nn = map(umap, \(u) {
      FNN::get.knn(as.matrix(u[, c("dim1", "dim2")]), k = K)$nn.index
    })
  )
overlap <- \(nn_a, nn_ref) {
  map_dbl(seq_len(nrow(nn_a)), \(i) {
    length(intersect(nn_a[i, ], nn_ref[i, ])) / K
  })
}
ref_nn <- neighbours |>
  dplyr::filter(arm %in% REFERENCES, seed == UMAP_SEED) |>
  select(reference = arm, marker_set, scaling, nn_ref = nn)
comparisons <- neighbours |>
  mutate(
    comparison = if_else(
      seed == UMAP_SEED,
      arm,
      paste0("1.8.8 zip, seed ", seed)
    )
  ) |>
  select(comparison, arm, seed, marker_set, scaling, nn) |>
  inner_join(
    ref_nn,
    by = join_by(marker_set, scaling),
    relationship = "many-to-many"
  ) |>
  # an arm is not scored against itself; the noise seeds are scored only against 1.8.8 zip seed 42
  dplyr::filter(
    comparison != reference,
    seed == UMAP_SEED | reference == "1.8.8 zip"
  )
knn_overlap <- comparisons |>
  mutate(
    per_cell = map2(nn, nn_ref, \(a, r) {
      tibble(donor = cells$donor, overlap = overlap(a, r))
    })
  ) |>
  select(comparison, reference, marker_set, scaling, per_cell) |>
  tidyr::unnest(per_cell) |>
  summarise(
    median_overlap = median(overlap),
    mean_overlap = mean(overlap),
    .by = c(comparison, reference, marker_set, scaling, donor)
  )
# with conditions, the sample's donor id and condition join the table; 3D and EXV keep their columns
if (!is.null(DS$conditions)) {
  knn_overlap <- knn_overlap |>
    left_join(
      select(smp, donor = sample, donor_id = donor, condition),
      by = join_by(donor),
      relationship = "many-to-one",
      unmatched = "error"
    )
}
write_csv(knn_overlap, file.path(RUN_DIR, "knn_overlap.csv"))

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
