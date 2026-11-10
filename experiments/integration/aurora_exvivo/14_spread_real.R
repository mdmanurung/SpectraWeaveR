# experiments/integration/aurora_exvivo/14_spread_real.R
# Purpose: positive control for the spread check. Inject spillover spread from a bright channel
#          into the negative of another into 25% of real donor files of one batch: target-negative
#          events get zero-mean noise whose SD grows with the source signal, so the target median
#          and the pair correlation barely move. Score spread_z for that directed pair against
#          neg_sd_z (target) and spill_corr_z.
# Run:     Rscript experiments/integration/aurora_exvivo/14_spread_real.R <extract_dir> <run_dir>
# Inputs:  <extract_dir>/events.rds, files.csv, negative_separation.csv   from 01 and 06
# Outputs: <run_dir>/spread_injected_long.csv
#          <run_dir>/fig_spread_real.png, figure_source_data_spread_real.csv
#
# Fault sizes, events, hyper_k and seeds are illustrative, as in 08_spill_real.R. The two
# spread_quantile values are reported as variants.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
extract_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) stop("run dir exists: ", run_dir, call. = FALSE)
dir.create(run_dir, recursive = TRUE)

BATCH <- "Batch1"
SEEDS <- 1:5
FRACTION_BAD <- 0.25
N_EVENTS <- 2000
HYPER_K <- 4
# Noise SD added to the target negative, in target-negative SDs, at the source's 90th percentile.
SPREAD <- c(0.5, 1)
QUANTILES <- c(0.9, 0.95)

ev <- readRDS(file.path(extract_dir, "events.rds"))
info <- read.csv(file.path(extract_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
base <- ev[info$file[info$batch == BATCH & grepl("BMV_[0-9]+[A-Z]", info$file)]]
names(base) <- paste0("f", seq_along(base))
neg_info <- read.csv(file.path(extract_dir, "negative_separation.csv"))
cut_of <- setNames(neg_info$neg_median + 3 * neg_info$neg_sd, neg_info$channel)
sd_of <- setNames(neg_info$neg_sd, neg_info$channel)
separated <- neg_info$channel[neg_info$separation >= 3]

auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}

long <- lapply(SEEDS, function(seed) {
  lapply(seq_along(SPREAD), function(i) {
    set.seed(seed * 1000 + i)
    bad <- seq_along(base) %in% sample(seq_along(base), round(FRACTION_BAD * length(base)))
    pick <- sample(separated, 2) # source, target
    files <- base
    files[bad] <- lapply(files[bad], function(x) {
      s <- x[, pick[1]]
      tgt_neg <- x[, pick[2]] < cut_of[[pick[2]]]
      # Source signal above its cut, scaled to 1 at the source's 90th percentile.
      k <- pmax(s - cut_of[[pick[1]]], 0) / (stats::quantile(s, 0.9) - cut_of[[pick[1]]])
      x[tgt_neg, pick[2]] <- x[tgt_neg, pick[2]] +
        stats::rnorm(sum(tgt_neg)) * SPREAD[[i]] * sd_of[[pick[2]]] * k[tgt_neg]
      x
    })
    lapply(QUANTILES, function(q) {
      na <- sw_negative_audit(files, pick, N_EVENTS, seed, HYPER_K, pairs = matrix(pick, 1), spread_quantile = q, min_separation = 3) # G24
      row_of <- function(check, item) {
        d <- na[na$check == check & na$channel == item, ]
        d$value[match(names(files), d$file)]
      }
      scores <- list(
        "spread_z (source > target)" = row_of("spread_z", paste(pick, collapse = " > ")),
        "neg_sd_z (target)" = abs(row_of("neg_sd_z", pick[2])),
        "spill_corr_z (source > target)" = abs(row_of("spill_corr_z", paste(pick, collapse = " > ")))
      )
      tibble(
        seed = seed,
        spread = SPREAD[[i]],
        spread_quantile = q,
        check = names(scores),
        auroc = vapply(scores, \(v) auroc(v, bad), numeric(1))
      )
    }) |>
      bind_rows()
  }) |>
    bind_rows()
}) |>
  bind_rows()
write.csv(long, file.path(run_dir, "spread_injected_long.csv"), row.names = FALSE)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(spread, spread_quantile, check)
  )
write.csv(mat, file.path(run_dir, "figure_source_data_spread_real.csv"), row.names = FALSE)
print(as.data.frame(mutate(mat, across(c(auroc, ci95_half), \(x) round(x, 2)))))
p <- ggplot(mat, aes(auroc, check)) +
  geom_vline(xintercept = c(0.5, 0.9), linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)), size = 0.2) +
  facet_grid(rows = vars(spread_quantile), cols = vars(spread), labeller = label_both) +
  labs(x = "AUROC on real files (mean, 95% CI over 5 seeds)", y = NULL)
ggsave(file.path(run_dir, "fig_spread_real.png"), p, width = 170, height = 70, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
