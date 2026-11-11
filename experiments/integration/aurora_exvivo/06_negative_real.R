# experiments/integration/aurora_exvivo/06_negative_real.R
# Purpose: run sw_negative_audit() on real files. (1) How bimodal is each real channel (GMM
#          separation), so how far the negative checks can be trusted per channel. (2) Inject
#          negative-population faults into a random channel of 25% of the real donor files of one
#          batch and score the hit channel's row against the whole-file checks.
# Run:     Rscript experiments/integration/aurora_exvivo/06_negative_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/negative_separation.csv     per channel: median GMM separation over files
#          <run_dir>/negative_injected_long.csv  per arm, seed, fault, check: AUROC
#          <run_dir>/fig_negative_real.png (fault arm only), figure_source_data_negative_real.csv
#            (both arms)
# G40: null arm (no injection, same draws). arm = "fault" or "null"; the null scores the clean
#      base files against the same bad labels and hit channels, after the fault arm of each
#      iteration, so the fault arm's random stream is unchanged.
#
# Every value below is illustrative; none is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
BATCH <- "Batch1"
SEEDS <- 1:5
FRACTION_BAD <- 0.25
N_EVENTS <- 2000
HYPER_K <- 4
FAULTS <- tibble::tribble(
  ~fault            , ~size ,
  "neg_shift"       , 0.3   , # negative SD units
  "neg_shift"       , 1     ,
  "neg_widen"       , 1.5   , # factor on the negative spread around its median
  "hyper_negative"  , 0.03  , # share of negative events moved to median - 8 negative SDs
  "angled_negative" , 0.3 # slope added from a second channel, joint negatives only
)

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
base <- ev[info$file[info$batch == BATCH & grepl("BMV_[0-9]+[A-Z]", info$file)]]
ch <- setdiff(colnames(base[[1]]), "Time")
names(base) <- paste0("f", seq_along(base)) # file names are not needed and stay out of outputs

# (1) Bimodality per channel on the real files ---------------------------------
clean <- sw_negative_audit(
  base,
  ch,
  N_EVENTS,
  seed = 1,
  hyper_k = HYPER_K,
  pairs = NULL,
  spread_quantile = NULL,
  min_separation = 3
) # G24
neg_info <- clean |>
  filter(check %in% c("separation", "neg_median", "neg_sd")) |>
  summarise(value = median(value, na.rm = TRUE), .by = c(channel, check)) |>
  tidyr::pivot_wider(names_from = check, values_from = value)
write.csv(
  neg_info,
  file.path(run_dir, "negative_separation.csv"),
  row.names = FALSE
)
cat(
  "channels with median GMM separation >= 3:",
  sum(neg_info$separation >= 3),
  "of",
  nrow(neg_info),
  "\n"
)

# (2) Injected faults ------------------------------------------------------------
auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}
# Negative events for injection: below the midpoint between the cohort negative median and the
# positive side, per channel (a fixed cut, so the fault does not depend on a per-file fit).
cut_of <- setNames(neg_info$neg_median + 3 * neg_info$neg_sd, neg_info$channel)
inject <- function(x, fault, size, a, b) {
  neg_b <- x[, b] < cut_of[[b]]
  med <- neg_info$neg_median[neg_info$channel == b]
  s <- neg_info$neg_sd[neg_info$channel == b]
  switch(
    fault,
    neg_shift = {
      x[neg_b, b] <- x[neg_b, b] + size * s
    },
    neg_widen = {
      x[neg_b, b] <- med + (x[neg_b, b] - med) * size
    },
    hyper_negative = {
      idx <- which(neg_b)[seq_len(round(size * sum(neg_b)))]
      x[idx, b] <- med - 8 * s
    },
    angled_negative = {
      both <- neg_b & x[, a] < cut_of[[a]]
      x[both, b] <- x[both, b] + size * (x[both, a] - stats::median(x[both, a]))
    }
  )
  x
}

# Scores one set of files against the bad-file labels, for the fault and the null arm.
score_arm <- function(files, bad, pick, seed) {
  na <- sw_negative_audit(
    files,
    ch,
    N_EVENTS,
    seed = seed,
    hyper_k = HYPER_K,
    pairs = matrix(pick, 1),
    spread_quantile = NULL,
    min_separation = 3 # G24
  )
  sa <- sw_fcs_signal_audit(
    files,
    ch,
    NULL,
    NULL,
    time_channel = NULL,
    time_bin = NULL,
    area_height = NULL,
    area_height_quantile = NULL,
    reference_group = NULL,
    reference = NULL
  )
  row_of <- function(tab, check, item) {
    d <- tab[tab$check == check & tab$channel == item, ]
    abs(d$value)[match(names(files), d$file)]
  }
  pair <- paste(pick, collapse = " ~ ")
  scores <- list(
    "neg_median_z" = row_of(na, "neg_median_z", pick[2]),
    "neg_sd_z" = row_of(na, "neg_sd_z", pick[2]),
    "hyper_negative_share" = row_of(na, "hyper_negative_share", pick[2]),
    "neg_corr_z" = row_of(na, "neg_corr_z", pair),
    "median_z (whole file)" = row_of(sa, "median_z", pick[2]),
    "sd_z (whole file)" = row_of(sa, "sd_z", pick[2])
  )
  tibble(
    check = names(scores),
    auroc = vapply(scores, \(v) auroc(v, bad), numeric(1))
  )
}

long <- lapply(SEEDS, function(seed) {
  lapply(seq_len(nrow(FAULTS)), function(i) {
    fault <- FAULTS$fault[[i]]
    size <- FAULTS$size[[i]]
    set.seed(seed * 1000 + i)
    bad <- seq_along(base) %in%
      sample(seq_along(base), round(FRACTION_BAD * length(base)))
    # Hit channels drawn from channels the GMM separates (median separation >= 3).
    pick <- sample(neg_info$channel[neg_info$separation >= 3], 2)
    files <- base
    files[bad] <- lapply(
      files[bad],
      inject,
      fault = fault,
      size = size,
      a = pick[1],
      b = pick[2]
    )
    fault_arm <- score_arm(files, bad, pick, seed)
    # G40 null arm: same bad labels and hit channels, clean files. It runs last in the iteration
    # (the package calls may set the seed); the next iteration starts with its own set.seed().
    null_arm <- score_arm(base, bad, pick, seed)
    bind_rows(
      mutate(
        fault_arm,
        arm = "fault",
        seed = seed,
        fault = fault,
        size = size,
        .before = 1
      ),
      mutate(
        null_arm,
        arm = "null",
        seed = seed,
        fault = fault,
        size = size,
        .before = 1
      )
    )
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(
  long,
  file.path(run_dir, "negative_injected_long.csv"),
  row.names = FALSE
)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(arm, fault, size, check)
  ) |>
  mutate(case = paste(fault, size))
write.csv(
  mat,
  file.path(run_dir, "figure_source_data_negative_real.csv"),
  row.names = FALSE
)
print(as.data.frame(mutate(mat, auroc = round(auroc, 2))[c(
  "arm",
  "case",
  "check",
  "auroc"
)]))
p <- ggplot(filter(mat, arm == "fault"), aes(auroc, check)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(
    aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)),
    size = 0.2
  ) +
  facet_wrap(vars(case), nrow = 1) +
  labs(
    x = "AUROC on real files (mean, 95% CI over 5 seeds; hit channel's row)",
    y = NULL
  )
ggsave(
  file.path(run_dir, "fig_negative_real.png"),
  p,
  width = 220,
  height = 70,
  units = "mm",
  dpi = 200
)
