# experiments/benchmarks/06_anomaly-features/run.R
# Purpose: frozen benchmark v1 for sw_anomaly_features() + sw_anomaly_flag_outlier().
#          Scores how well the forest score ranks injected-anomaly files (AUROC).
# Run:     Rscript experiments/benchmarks/06_anomaly-features/run.R <candidate> [shuffle]
#          candidate: both | quantile | emd | emd_reference
#          shuffle:   shuffle labels to check the evaluator rewards no noise
# Inputs:  none (simulated here; card in the tusoskill run folder)
# Outputs: one JSON line on stdout: mean_auroc, per-effect AUROC, runtime_s
#
# All values below are benchmark-only and illustrative. They are not package defaults.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source("R/corrector.R")
source("R/anomaly_features.R")
source("R/anomaly_flag.R")

args <- commandArgs(trailingOnly = TRUE)
CANDIDATE <- args[[1]]
SHUFFLE <- length(args) > 1 && args[[2]] == "shuffle"

# Benchmark settings (frozen) ------------------------------------------------
# Each value can be overridden by an environment variable (SW_<NAME>) for a
# sensitivity run. The defaults are the frozen v1 values; the primary metric is
# only ever reported at the defaults.
env_num <- function(name, default) {
  as.numeric(Sys.getenv(paste0("SW_", name), as.character(default)))
}
SEED_START <- env_num("SEED_START", 1)
SEEDS <- SEED_START + 0:(env_num("N_SEEDS", 10) - 1) # simulated cohorts
N_CLEAN <- env_num("N_CLEAN", 60) # clean files per cohort
N_BAD <- env_num("N_BAD", 12) # injected files per cohort
N_CH <- env_num("N_CH", 8) # channels
N_EV <- env_num("N_EV", 1000) # events per file
FILE_JITTER <- env_num("FILE_JITTER", 0.1) # SD of file-to-file population-mean jitter
EFFECT_SCALE <- env_num("EFFECT_SCALE", 1) # multiplies every injected effect size
N_CH_AFFECTED <- env_num("N_CH_AFFECTED", 1) # channels hit in an injected file
POP_ONLY <- env_num("POP_ONLY", 0) == 1 # effect only on the rarest population
SCORER <- Sys.getenv("SW_SCORER", "iforest") # iforest | novelty | robustz
EFFECTS <- c("shift", "stretch", "bimodal") # injected effects, cycled over bad files
# Feature and flag settings shared by every candidate (illustrative, not tuned):
QSTEP <- env_num("QSTEP", 0.1)
NEV <- env_num("NEV", 500)
AGG <- env_num("AGG", 5000)
NTREES <- env_num("NTREES", 300)
FOREST_SAMPLE <- env_num("FOREST_SAMPLE", 32)
FOREST_NDIM <- env_num("FOREST_NDIM", 1)
CHANNELS <- paste0("C", seq_len(N_CH))

# Simulation -----------------------------------------------------------------
simulate_cohort <- function(seed) {
  set.seed(seed)
  mu <- matrix(rep(c(-2, 0, 2), N_CH), nrow = 3)
  make_file <- function() {
    pop <- sample(1:3, N_EV, replace = TRUE, prob = c(0.5, 0.3, 0.2))
    jit <- matrix(stats::rnorm(3 * N_CH, 0, FILE_JITTER), 3)
    x <- mu[pop, ] + jit[pop, ] + matrix(stats::rnorm(N_EV * N_CH), N_EV)
    colnames(x) <- CHANNELS
    attr(x, "pop") <- pop
    x
  }
  files <- lapply(seq_len(N_CLEAN + N_BAD), function(i) make_file())
  effect <- c(rep("clean", N_CLEAN), rep(EFFECTS, length.out = N_BAD))
  for (i in which(effect != "clean")) {
    x <- files[[i]]
    rows <- if (POP_ONLY) attr(x, "pop") == 3 else rep(TRUE, nrow(x))
    for (ch in sample(CHANNELS, N_CH_AFFECTED)) {
      v <- x[rows, ch]
      x[rows, ch] <- switch(
        effect[i],
        shift = v + 0.5 * EFFECT_SCALE,
        stretch = mean(v) + (v - mean(v)) * (1 + 0.5 * EFFECT_SCALE),
        bimodal = v +
          sample(c(-0.5, 0.5) * EFFECT_SCALE, length(v), replace = TRUE)
      )
    }
    files[[i]] <- x
  }
  ord <- sample(length(files)) # file order carries no label
  files <- files[ord]
  names(files) <- paste0("f", seq_along(files))
  list(files = files, effect = effect[ord])
}

auroc <- function(score, pos) {
  r <- rank(score)
  n1 <- sum(pos)
  n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Candidates -----------------------------------------------------------------
cand <- list(
  both = list(features = c("quantile", "emd"), from = "all"),
  quantile = list(features = "quantile", from = "all"),
  emd = list(features = "emd", from = "all"),
  emd_reference = list(features = c("quantile", "emd"), from = "reference")
)[[CANDIDATE]]

t0 <- Sys.time()
res <- lapply(SEEDS, function(s) {
  sim <- simulate_cohort(s)
  ids <- names(sim$files)
  # Reference = 30 known-clean files, as a lab would hold known-good controls.
  ref <- head(ids[sim$effect == "clean"], 30)
  feat <- sw_anomaly_features(
    sim$files,
    channels = CHANNELS,
    features = cand$features,
    quantile_step = QSTEP,
    n_events = NEV,
    seed = s,
    scale = "simulated",
    aggregate_from = cand$from,
    agg_size = AGG,
    reference = ref
  )
  fl <- if (SCORER == "iforest") {
    sw_anomaly_flag_outlier(
      feat,
      ntrees = NTREES,
      sample_size = FOREST_SAMPLE,
      ndim = FOREST_NDIM,
      cutoff = 0.5,
      seed = s
    )
  } else if (SCORER == "novelty") {
    # Score = minus log-density; reference files are left out of the evaluation.
    nv <- sw_anomaly_flag_novelty(feat, reference = ref, quantile = 0.05)
    nv <- nv[nv$role == "test", ]
    tibble(file = nv$file, score = -nv$log_density)
  } else {
    # Comparator: sum of robust z-scores over all features.
    m <- as.matrix(feat[-1])
    z <- sweep(
      abs(sweep(m, 2, apply(m, 2, median))),
      2,
      apply(m, 2, mad) + 1e-12,
      "/"
    )
    tibble(file = feat$file, score = rowSums(z))
  }
  lab <- sim$effect[match(fl$file, ids)]
  if (SHUFFLE) {
    set.seed(1000 + s)
    lab <- sample(lab)
  }
  pair_auroc <- function(e) {
    keep <- lab %in% c("clean", e)
    auroc(fl$score[keep], lab[keep] == e)
  }
  tibble(
    seed = s,
    effect = c("all", EFFECTS),
    auroc = c(
      auroc(fl$score, lab != "clean"),
      vapply(EFFECTS, pair_auroc, numeric(1))
    )
  )
}) |>
  bind_rows()

by_effect <- res |> summarise(auroc = mean(auroc), .by = effect)
out <- list(
  candidate = CANDIDATE,
  shuffle = SHUFFLE,
  mean_auroc = by_effect$auroc[by_effect$effect == "all"],
  # 95% interval of the mean over cohorts (normal approximation, few seeds).
  ci95_half = 1.96 * sd(res$auroc[res$effect == "all"]) / sqrt(length(SEEDS)),
  by_effect = as.list(setNames(by_effect$auroc, by_effect$effect)),
  runtime_s = as.numeric(difftime(Sys.time(), t0, units = "secs"))
)
cat(jsonlite::toJSON(out, auto_unbox = TRUE, digits = 4), "\n")
