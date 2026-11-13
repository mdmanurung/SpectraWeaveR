# experiments/integration/aurora_exvivo/08_spill_real.R
# Purpose: inject a leak from a bright channel into a negative one (tandem breakdown, dye
#          interaction) into 25% of real donor files of one batch, and score spill_corr_z for that
#          directed pair against the whole-file checks.
# Run:     Rscript experiments/integration/aurora_exvivo/08_spill_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv, negative_separation.csv   from 01 and 06
# Outputs: <run_dir>/spill_injected_long.csv  per arm, seed, leak, check: AUROC
#          <run_dir>/fig_spill_real.png (fault arm only), figure_source_data_spill_real.csv (both arms)
# G40: null arm (no injection, same draws). arm = "fault" or "null"; the null scores the clean
#      base files against the same bad labels and source/target pair, after the fault arm of each
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
LEAK <- c(0.05, 0.15) # share of the source signal added to the target, source-positive events only

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
base <- ev[info$file[info$batch == BATCH & grepl("BMV_[0-9]+[A-Z]", info$file)]]
names(base) <- paste0("f", seq_along(base))
ch <- setdiff(colnames(base[[1]]), "Time")
neg_info <- read.csv(file.path(run_dir, "negative_separation.csv"))
cut_of <- setNames(neg_info$neg_median + 3 * neg_info$neg_sd, neg_info$channel)
separated <- neg_info$channel[neg_info$separation >= 3]

auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}

# Scores one set of files against the bad-file labels, for the fault and the null arm.
score_arm <- function(files, bad, pick, seed) {
  na <- sw_negative_audit(
    files,
    pick,
    N_EVENTS,
    seed,
    HYPER_K,
    pairs = matrix(pick, 1),
    spread_quantile = NULL,
    min_separation = 3
  ) # G24
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
    reference = NULL,
    summary = NULL
  )
  row_of <- function(tab, check, item) {
    d <- tab[tab$check == check & tab$channel %in% item, ]
    # A file whose rows are all NA (too few events) stays NA, not -Inf.
    top <- \(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    as.numeric(tapply(abs(d$value), d$file, top)[names(files)])
  }
  scores <- list(
    "spill_corr_z (source > target)" = row_of(
      na,
      "spill_corr_z",
      paste(pick, collapse = " > ")
    ),
    "neg_corr_z (pair)" = row_of(
      na,
      "neg_corr_z",
      paste(pick, collapse = " ~ ")
    ),
    "neg_sd_z (target)" = row_of(na, "neg_sd_z", pick[2]),
    "corr_shift (whole file, max over pairs)" = row_of(
      sa[sa$check == "corr_shift", ],
      "corr_shift",
      NA
    ),
    "sd_z (whole file, target)" = row_of(sa, "sd_z", pick[2])
  )
  tibble(
    check = names(scores),
    auroc = vapply(scores, \(v) auroc(v, bad), numeric(1))
  )
}

long <- lapply(SEEDS, function(seed) {
  lapply(seq_along(LEAK), function(i) {
    set.seed(seed * 1000 + i)
    bad <- seq_along(base) %in%
      sample(seq_along(base), round(FRACTION_BAD * length(base)))
    pick <- sample(separated, 2) # source, target
    files <- base
    files[bad] <- lapply(files[bad], function(x) {
      # Leak only where the source is positive and the target negative; centred source signal.
      leak <- x[, pick[1]] >= cut_of[[pick[1]]] &
        x[, pick[2]] < cut_of[[pick[2]]]
      x[leak, pick[2]] <- x[leak, pick[2]] +
        LEAK[[i]] * (x[leak, pick[1]] - stats::median(x[leak, pick[1]]))
      x
    })
    fault_arm <- score_arm(files, bad, pick, seed)
    # G40 null arm: same bad labels and pair, clean files. It runs last in the iteration (the
    # package calls may set the seed); the next iteration starts with its own set.seed().
    null_arm <- score_arm(base, bad, pick, seed)
    bind_rows(
      mutate(
        fault_arm,
        arm = "fault",
        seed = seed,
        leak = LEAK[[i]],
        .before = 1
      ),
      mutate(null_arm, arm = "null", seed = seed, leak = LEAK[[i]], .before = 1)
    )
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(
  long,
  file.path(run_dir, "spill_injected_long.csv"),
  row.names = FALSE
)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(arm, leak, check)
  )
write.csv(
  mat,
  file.path(run_dir, "figure_source_data_spill_real.csv"),
  row.names = FALSE
)
print(as.data.frame(mutate(mat, auroc = round(auroc, 2))))
p <- ggplot(filter(mat, arm == "fault"), aes(auroc, check)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(
    aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)),
    size = 0.2
  ) +
  facet_wrap(vars(leak), labeller = label_both) +
  labs(x = "AUROC on real files (mean, 95% CI over 5 seeds)", y = NULL)
ggsave(
  file.path(run_dir, "fig_spill_real.png"),
  p,
  width = 170,
  height = 60,
  units = "mm",
  dpi = 200
)
