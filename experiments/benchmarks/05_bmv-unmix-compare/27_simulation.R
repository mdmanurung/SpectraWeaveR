# Purpose: G2 simulation. Ground-truth accuracy of unmixing with known abundances. Simulate raw Aurora events
#          with AutoSpectral::sim.flow.data (1.8.8 spectra, AF spectra, variants, real scatter), unmix them
#          with three code lines, and score unmixed minus truth against plain OLS.
# Run:     from the repo root, one stage at a time (each stage fails loudly; never overwrites):
#            pixi run Rscript <this> sim   <sim_dir> <n_cells>   # simulate, write raw FCS + truth
#            pixi run Rscript <this> unmix <sim_dir> <arm>       # arm = as188_joint | as102_var | as102_basic
#            pixi run Rscript <this> score <sim_dir> <run_dir>   # tables, OLS floor
#          <sim_dir> = experiments/bmv_test/compare_3D/data/simulation/<tag> (full: tag full_1e5, n_cells 100000;
#          smoke: tag smoke_5000, n_cells 5000, run_dir outside experiments/runs).
#          Full run_dir = experiments/runs/2026-11-10_bmv-simulation
# Arms (same 1.8.8 spectra, AF spectra and controls; only code and variants differ):
#   as188_joint  = 1.8.8 unmix.folder, pipeline joint, speed slow (fast in as188_joint; no effect under joint), AF spectra, variants (as 06_unmix_c188.R cell s188)
#   as102_var    = 1.0.2 unmix.folder, weighted FALSE, AF spectra, variants (as 06_unmix_c102.R cell s188_143803)
#   as102_basic  = as as102_var without variants (as the as102_basic arm of 16_unmix_controls.R)
#   ols          = plain least squares in base R, no AF, no variants (floor; computed in the score stage)
# Inputs:  experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/table_spectra/ (spectra, AF csv) [user zip]
#          .../autoSpectral_v1.8.8/figure_spectral_variants/Spectral_variants.rds [user zip]
#          .../autoSpectral_v1.8.8/Reference Group Clean (flow.control for unmix.folder) [user zip]
#          .../autoSpectral_v1.8.8/Raw/*.fcs (5 stained samples: scatter and FCS template) [user zip]
#          .../lib_as102 (AutoSpectral 1.0.2) [06_unmix_c102.R]
# Outputs: <sim_dir>/raw/sim_seed<seed>.fcs simulated raw events, detector names as the instrument
#          <sim_dir>/truth/sim_seed<seed>.rds full sim.flow.data truth and params
#          <sim_dir>/truth/sim_seed<seed>_truth.csv.gz per event: abundance and layer per fluorophore, AF row and abundance
#          <sim_dir>/unmixed/<arm>/*.fcs and args.txt (what was run, runtime)
#          <sim_dir>/sim_args.txt every sim.flow.data argument and its source
#          <run_dir>/sim_scores.csv per seed, arm, fluorophore, stratum: n, rmse, bias, spread, pileup
#          <run_dir>/sim_summary.csv per arm: median over fluorophores of the metrics, ratio to as188_joint, difference to OLS
#          <run_dir>/sessionInfo.txt

# Parameters ------------------------------------------------------------------------------
DATA <- normalizePath("experiments/bmv_test/compare_3D/data")
LIB_102 <- file.path(DATA, "lib_as102") # AutoSpectral and AutoSpectralRcpp v1.0.2 (separate library)
SEEDS <- c(42L, 43L, 44L) # 3 independent draws (user, 2026-11-10)
DIM_ABUNDANCE <- 0.005 # relative abundance of the dim layer = layer.centroids[3], the package default (user, 2026-11-10)
DIM_LAYER <- "low" # sim.flow.data names the third centroid layer "low"
SPECTRA_FILE <- "Automated_autospectral_spectra_20261105_143803.csv" # 1.8.8 run used for as188_joint
AF_FILE <- "H12 Unstained_BMV_Aurora_3dCD328_Batch1.fcs Autofluorescence spectra 20261105_143933.csv"
WEIGHTED_102 <- FALSE # as in 06_unmix_c102.R
PIPELINE_188 <- "joint" # as in 06_unmix_c188.R
SPEED_188 <- "slow" # always slow (user, 2026-11-10); no effect under pipeline "joint" (repro188 = s188_c188_slow); was "fast" to match as188_joint
PILEUP_CUT <- sinh(0.2) * 6000 # |x| below this is "at zero" on arcsinh cofactor 6000 (user, 2026-11-10; cofactor: CLAUDE.md)
NEG_QUANTILES <- c(0.025, 0.975) # negative spread = (q97.5 - q2.5) / 3.92, the SD of a normal (user, 2026-11-10)
REF_ARM <- "as188_joint" # reference arm (README)
ARMS <- list(
  as188_joint = list(code = "1.8.8", variants = TRUE),
  as102_var = list(code = "1.0.2", variants = TRUE),
  as102_basic = list(code = "1.0.2", variants = FALSE)
)

# Setup -----------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stage <- args[[1]]
sim_dir <- args[[2]]
stopifnot(
  "stage is sim, unmix or score" = stage %in% c("sim", "unmix", "score")
)
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
ctrl_dir <- file.path(d188, "Reference Group Clean")
tab_dir <- file.path(d188, "table_spectra")
if (stage == "unmix" && ARMS[[args[[3]]]]$code == "1.0.2") {
  .libPaths(c(LIB_102, .libPaths()))
}
suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})
if (stage == "unmix" && ARMS[[args[[3]]]]$code == "1.0.2") {
  suppressPackageStartupMessages(library(AutoSpectralRcpp))
}
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
spectra <- read.spectra(SPECTRA_FILE, spectra.dir = tab_dir)
spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in the 06 scripts
af.spectra <- read.spectra(AF_FILE, spectra.dir = tab_dir)
fluors <- rownames(spectra)

# Stage sim -------------------------------------------------------------------------------
if (stage == "sim") {
  n_cells <- as.integer(args[[3]])
  stopifnot("sim dir must not exist" = !dir.exists(sim_dir))
  stopifnot(
    "dim abundance is the package default" = formals(
      sim.flow.data
    )$layer.centroids[[4]] ==
      DIM_ABUNDANCE
  )
  walk(c("raw", "truth"), \(x) {
    dir.create(file.path(sim_dir, x), recursive = TRUE)
  })
  variants <- readRDS(file.path(
    d188,
    "figure_spectral_variants",
    "Spectral_variants.rds"
  ))

  # real stained samples: scatter pool for sim.flow.data and the FCS template
  real <- list.files(
    file.path(d188, "Raw"),
    pattern = "\\.fcs$",
    full.names = TRUE
  ) |>
    map(\(f) read.FCS(f, truncate_max_range = FALSE, transformation = FALSE))
  real_exprs <- map(real, exprs)
  scatter_pool <- real_exprs |>
    map(\(x) x[, c("FSC-A", "SSC-A")]) |>
    reduce(rbind)
  # height columns are not simulated: scale area by the pooled median H/A ratio (not used by unmixing)
  pool <- reduce(real_exprs, rbind)
  h_ratio <- c(
    "FSC-H" = median(pool[, "FSC-H"] / pool[, "FSC-A"]),
    "SSC-H" = median(pool[, "SSC-H"] / pool[, "SSC-A"]),
    "SSC-B-H" = median(pool[, "SSC-B-H"] / pool[, "SSC-B-A"])
  )
  template <- real[[1]]
  time_range <- range(pool[, "Time"])
  rm(pool, real_exprs)

  sim_one <- function(seed) {
    sim <- sim.flow.data(
      spectra = spectra,
      asp = asp,
      n.cells = n_cells,
      af.spectra = af.spectra,
      variants = variants,
      scatter.data = scatter_pool,
      seed = seed
    )
    stopifnot(
      "sim detectors are FCS detectors" = all(
        colnames(sim$raw) %in% colnames(template)
      )
    )
    m <- matrix(
      NA_real_,
      n_cells,
      ncol(template),
      dimnames = list(NULL, colnames(template))
    )
    m[, colnames(sim$raw)] <- sim$raw
    m[, "FSC-A"] <- sim$truth$scatter[, 1]
    m[, "SSC-A"] <- sim$truth$scatter[, 2]
    m[, "SSC-B-A"] <- sim$truth$scatter[, 2]
    m[, "FSC-H"] <- m[, "FSC-A"] * h_ratio[["FSC-H"]]
    m[, "SSC-H"] <- m[, "SSC-A"] * h_ratio[["SSC-H"]]
    m[, "SSC-B-H"] <- m[, "SSC-B-A"] * h_ratio[["SSC-B-H"]]
    m[, "Time"] <- seq(time_range[[1]], time_range[[2]], length.out = n_cells)
    fn <- sprintf("sim_seed%d", seed)
    ff <- new(
      "flowFrame",
      exprs = m,
      parameters = parameters(template),
      description = keyword(template)
    )
    keyword(ff) <- list(
      "$FIL" = paste0(fn, ".fcs"),
      "$TOT" = as.character(n_cells)
    )
    write.FCS(ff, file.path(sim_dir, "raw", paste0(fn, ".fcs")))
    saveRDS(
      sim[c("truth", "params")],
      file.path(sim_dir, "truth", paste0(fn, ".rds"))
    )
    tr <- sim$truth
    bind_cols(
      tibble(
        event = seq_len(n_cells),
        af_row = tr$af.row,
        af_abundance = tr$af.abundance,
        fsc_a = tr$scatter[, 1],
        ssc_a = tr$scatter[, 2]
      ),
      as_tibble(tr$abundances[, fluors]) |>
        rename_with(\(x) paste0("abundance__", x)),
      as_tibble(tr$expression.layer) |> rename_with(\(x) paste0("layer__", x))
    ) |>
      write_csv(file.path(sim_dir, "truth", paste0(fn, "_truth.csv.gz")))
  }
  t0 <- Sys.time()
  walk(SEEDS, sim_one)
  writeLines(
    c(
      "sim.flow.data arguments (source: ?sim.flow.data in AutoSpectral 1.8.8; unlisted arguments are package defaults)",
      paste("spectra:", file.path(tab_dir, SPECTRA_FILE), "[user, 2026-11-10]"),
      paste("asp: get.autospectral.param(cytometer = 'aurora') [06 scripts]"),
      paste("n.cells:", n_cells, "[user, 2026-11-10; smoke test uses fewer]"),
      paste("af.spectra:", file.path(tab_dir, AF_FILE), "[user, 2026-11-10]"),
      paste(
        "variants:",
        file.path(d188, "figure_spectral_variants", "Spectral_variants.rds"),
        "[user, 2026-11-10]"
      ),
      paste(
        "scatter.data: FSC-A and SSC-A pooled from",
        nrow(scatter_pool),
        "events of the 5 Raw stained samples [user, 2026-11-10]"
      ),
      paste("seed:", paste(SEEDS, collapse = ", "), "[user, 2026-11-10]"),
      paste(
        "layer.centroids (default):",
        paste(formals(sim.flow.data)$layer.centroids[-1], collapse = ", "),
        "; dim layer =",
        DIM_ABUNDANCE
      ),
      "all other arguments: package default (complexity 1/3, layer.cv 0.2, af.scale.range c(0.001, 0.05), af.scatter.slope 0.1, shot.noise, spillover.noise, detector.noise, af.variation, spectral.variation all TRUE)",
      "FCS: detector columns = sim$raw; FSC-A, SSC-A (and SSC-B-A = SSC-A) = sim scatter; H columns = area x pooled median H/A; Time = even grid over real range",
      paste(
        "elapsed_min:",
        format(round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
      )
    ),
    file.path(sim_dir, "sim_args.txt")
  )
}

# Stage unmix -----------------------------------------------------------------------------
if (stage == "unmix") {
  arm <- args[[3]]
  a <- ARMS[[arm]]
  stopifnot(
    "requested AutoSpectral version loaded" = as.character(packageVersion(
      "AutoSpectral"
    )) ==
      a$code
  )
  out_dir <- file.path(sim_dir, "unmixed", arm)
  stopifnot("out dir must not exist" = !dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)
  flow.control <- reload.flow.control(
    ctrl_dir,
    file.path(ctrl_dir, "fcs_control_file.csv"),
    asp
  )
  variants <- if (a$variants) {
    readRDS(file.path(
      d188,
      "figure_spectral_variants",
      "Spectral_variants.rds"
    ))
  } else {
    NULL
  }
  unmix_args <- list(
    fcs.dir = file.path(sim_dir, "raw"),
    spectra = spectra,
    asp = asp,
    flow.control = flow.control,
    method = "AutoSpectral",
    af.spectra = af.spectra,
    spectra.variants = variants,
    output.dir = out_dir,
    parallel = TRUE
  )
  if (a$code == "1.0.2") {
    unmix_args <- c(unmix_args, list(weighted = WEIGHTED_102, threads = NULL))
  } else {
    unmix_args <- c(
      unmix_args,
      list(pipeline = PIPELINE_188, speed = SPEED_188, verbose = TRUE)
    )
  }
  writeLines(
    c(
      paste("arm:", arm),
      paste("AutoSpectral:", a$code),
      paste("raw:", file.path(sim_dir, "raw")),
      paste("spectra:", file.path(tab_dir, SPECTRA_FILE)),
      paste("af:", file.path(tab_dir, AF_FILE)),
      paste(
        "variants:",
        if (a$variants) "Spectral_variants.rds (1.8.8)" else "none"
      ),
      paste("flow.control from:", ctrl_dir),
      paste(
        "settings:",
        if (a$code == "1.0.2") {
          paste("weighted", WEIGHTED_102)
        } else {
          paste("pipeline", PIPELINE_188, "speed", SPEED_188)
        }
      )
    ),
    file.path(out_dir, "args.txt")
  )
  t0 <- Sys.time()
  do.call(unmix.folder, unmix_args)
  cat(
    "elapsed_min:",
    format(round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)),
    "\n",
    file = file.path(out_dir, "args.txt"),
    append = TRUE
  )
}

# Stage score -----------------------------------------------------------------------------
if (stage == "score") {
  run_dir <- args[[3]]
  stopifnot("run dir must not exist" = !dir.exists(run_dir))
  dir.create(run_dir, recursive = TRUE)

  # per fluorophore and stratum: all events, dim layer, truly negative events
  score_one <- function(x, tru, lyr) {
    err <- x - tru
    is_dim <- !is.na(lyr) & lyr == DIM_LAYER
    is_neg <- tru == 0
    tibble(
      stratum = c("all", "dim", "neg"),
      n = c(length(x), sum(is_dim), sum(is_neg)),
      rmse = c(
        sqrt(mean(err^2)),
        sqrt(mean(err[is_dim]^2)),
        sqrt(mean(err[is_neg]^2))
      ),
      bias = c(mean(err), mean(err[is_dim]), mean(err[is_neg])),
      spread = c(
        NA,
        NA,
        diff(quantile(x[is_neg], NEG_QUANTILES, names = FALSE)) / 3.92
      ),
      pileup = c(NA, NA, mean(abs(x[is_neg]) < PILEUP_CUT))
    )
  }
  scores <- map(set_names(SEEDS), \(seed) {
    fn <- sprintf("sim_seed%d", seed)
    tr <- readRDS(file.path(sim_dir, "truth", paste0(fn, ".rds")))$truth
    raw_fcs <- exprs(read.FCS(
      file.path(sim_dir, "raw", paste0(fn, ".fcs")),
      truncate_max_range = FALSE,
      transformation = FALSE
    ))
    # plain OLS on the file the arms read (float32 as written)
    ols <- t(solve(
      tcrossprod(spectra),
      spectra %*% t(raw_fcs[, colnames(spectra)])
    ))
    colnames(ols) <- fluors
    read_arm <- function(arm) {
      f <- list.files(
        file.path(sim_dir, "unmixed", arm),
        pattern = fn,
        full.names = TRUE
      )
      stopifnot("one unmixed file per arm and draw" = length(f) == 1)
      u <- exprs(read.FCS(
        f,
        truncate_max_range = FALSE,
        transformation = FALSE
      ))
      stopifnot(
        "unmixed events align with raw events" = isTRUE(all.equal(
          u[, "FSC-A"],
          raw_fcs[, "FSC-A"],
          tolerance = 1e-6
        ))
      )
      x <- u[, paste0(fluors, "-A")]
      colnames(x) <- fluors
      x
    }
    mats <- c(list(ols = ols), map(set_names(names(ARMS)), read_arm))
    imap(mats, \(x, arm) {
      map(set_names(fluors), \(f) {
        score_one(x[, f], tr$abundances[, f], tr$expression.layer[, f])
      }) |>
        list_rbind(names_to = "fluor")
    }) |>
      list_rbind(names_to = "arm")
  }) |>
    list_rbind(names_to = "seed")
  write_csv(scores, file.path(run_dir, "sim_scores.csv"))

  # mean over the 3 draws per arm, fluorophore and stratum; then median over fluorophores
  by_fluor <- scores |>
    summarise(
      across(c(rmse, bias, spread, pileup), mean),
      .by = c(arm, fluor, stratum)
    ) |>
    mutate(abs_bias = abs(bias), .keep = "unused") |>
    pivot_longer(
      c(rmse, abs_bias, spread, pileup),
      names_to = "metric",
      values_drop_na = TRUE
    )
  ref <- by_fluor |>
    filter(arm == REF_ARM) |>
    select(fluor, stratum, metric, ref_value = value)
  ols_ref <- by_fluor |>
    filter(arm == "ols") |>
    select(fluor, stratum, metric, ols_value = value)
  summary_tbl <- by_fluor |>
    left_join(
      ref,
      by = join_by(fluor, stratum, metric),
      relationship = "many-to-one",
      unmatched = "error"
    ) |>
    left_join(
      ols_ref,
      by = join_by(fluor, stratum, metric),
      relationship = "many-to-one",
      unmatched = "error"
    ) |>
    summarise(
      n_fluor = n(),
      median_value = median(value),
      median_ratio_to_ref = median(value / ref_value),
      median_diff_to_ols = median(value - ols_value),
      .by = c(arm, stratum, metric)
    ) |>
    arrange(metric, stratum, arm)
  write_csv(summary_tbl, file.path(run_dir, "sim_summary.csv"))
  writeLines(
    capture.output(sessionInfo()),
    file.path(run_dir, "sessionInfo.txt")
  )
}
