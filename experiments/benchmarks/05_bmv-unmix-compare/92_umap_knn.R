# Purpose: does the unmixing arm change the UMAP of CD3+ T cells? Grill 2026-11-10 (README row "UMAP kNN
#          test"). The same 10,000 CD3+ cells per donor (1.8.8 zip CD3pos rows, seed 42) are embedded in
#          every arm with UMAP; per cell, the share of its 30 nearest embedding neighbours that are also
#          its neighbours in the reference embedding; median per donor. References: 1.8.8 zip and
#          SpectroFlo. Noise floor: 1.8.8 zip re-embedded with seeds 1-3 vs seed 42. Two marker sets
#          (lineage; lineage + functional) and two scalings (z-score; none). Descriptive, no threshold.
# Run:     sbatch experiments/slurm/run_script_lib.sbatch \
#            experiments/benchmarks/05_bmv-unmix-compare/92_umap_knn.R experiments/runs/2026-11-10_bmv-umap-knn
#            (repo root; run dir must not exist). Figures: 92b_plot_umap_knn.R.
# Inputs:  experiments/runs/2026-11-10_bmv-ablation-gated/gs/{as188_joint,opt188_pregate_joint_p1af1,
#            af405_afrefine,as102_var}   GatingSets [19a_gate_arm.R via 51_score_ablation.R]
#          experiments/runs/2026-11-10_bmv-gated-pooled/gs/spectroflo   GatingSet [19a_gate_arm.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv
#            dye to marker map [user zip]
# Outputs: <run_dir>/cells.csv          donor, root row of each sampled CD3+ cell
#          <run_dir>/embeddings.rds     arm, marker set, scaling, seed, donor, cell, UMAP 1-2, TNF and Th2
#          <run_dir>/knn_overlap.csv    comparison, reference, marker set, scaling, donor, median and mean
#                                       per-cell neighbour overlap
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

# Parameters ---------------------------------------------------------------------------------------
RUN_DIR <- commandArgs(trailingOnly = TRUE)[1]
GATED <- "experiments/runs/2026-11-10_bmv-ablation-gated/gs"
GATED_POOLED <- "experiments/runs/2026-11-10_bmv-gated-pooled/gs"
# grill Q1 (user, 2026-11-10): reference first; SpectroFlo pairs row for row with the AutoSpectral arms
ARMS <- c(
  "1.8.8 zip" = file.path(GATED, "as188_joint"),
  opt188 = file.path(GATED, "opt188_pregate_joint_p1af1"),
  "refined AF" = file.path(GATED, "af405_afrefine"),
  "1.0.2 zip" = file.path(GATED, "as102_var"),
  SpectroFlo = file.path(GATED_POOLED, "spectroflo")
)
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
STORED_ARCSINH <- c(
  "SBUV510-A",
  "Spark Blue 574-A",
  "SBUV445-A",
  "BV510-A",
  "BV570-A"
)
# README "Cluster markers": lineage only (user, 2026-11-09)
LINEAGE <- c(
  "CD45RA",
  "TCRgd",
  "CD8",
  "CD3",
  "CD161",
  "CD25",
  "CD38",
  "CD7",
  "CD4",
  "CD56",
  "CD57",
  "CXCR5",
  "CCR7",
  "CD127"
)
# grill Q7: every other panel marker; viability and PE (no marker) left out
FUNCTIONAL <- c(
  "GLUT1",
  "Th2",
  "TNF",
  "IFNg",
  "IL-2",
  "CTLA-4",
  "ATP5A",
  "TCRvd2",
  "IL-10",
  "FoxP3",
  "CD71",
  "CS",
  "PD-1",
  "T-bet",
  "CPT1A",
  "ATGL",
  "G6PD"
)
MARKER_SETS <- list(
  lineage = LINEAGE,
  lineage_functional = c(LINEAGE, FUNCTIONAL)
) # grill Q3
CONTROL_FILE <- "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)

# Marker to channel ----------------------------------------------------------------------------------
channel_of <- read_csv(CONTROL_FILE, show_col_types = FALSE) |>
  dplyr::filter(marker %in% c(LINEAGE, FUNCTIONAL)) |>
  with(set_names(paste0(fluorophore, "-A"), marker))
channel_of <- channel_of[c(LINEAGE, FUNCTIONAL)]

# Sample the same CD3+ cells per donor from 1.8.8 zip's gates ----------------------------------------
gs_list <- map(ARMS, load_gs)
donors <- sampleNames(gs_list[["1.8.8 zip"]])
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
# scatter does not depend on unmixing, so equal FSC-A per row shows the arms hold the same cells
stopifnot(
  "rows pair across arms" = map_lgl(gs_list, \(gs) {
    identical(
      unlist(map(donors, \(d) {
        exprs(gh_pop_get_data(gs[[d]], "root"))[
          cells$row[cells$donor == d],
          "FSC-A"
        ]
      })),
      unlist(map(donors, \(d) {
        exprs(gh_pop_get_data(gs_list[["1.8.8 zip"]][[d]], "root"))[
          cells$row[cells$donor == d],
          "FSC-A"
        ]
      }))
    )
  }) |>
    all()
)
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
        bind_cols(select(mats[[a]], donor, cell, TNF, Th2))
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
write_csv(knn_overlap, file.path(RUN_DIR, "knn_overlap.csv"))

writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
