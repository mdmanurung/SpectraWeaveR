# experiments/integration/aurora_exvivo/09_batch_real.R
# Purpose: inject a batch-level fault (one marker in every donor file of one batch) into the real
#          multi-batch files and ask whether sw_batch_scores() ranks the hit marker above the 39
#          others. The real batches already differ (batch nested with donor), so the injected
#          marker must stand out above real batch differences.
# Run:     Rscript experiments/integration/aurora_exvivo/09_batch_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/batch_injected_long.csv  per arm, seed, fault, check: AUROC
#          <run_dir>/fig_batch_real.png (fault arm only), figure_source_data_batch_real.csv (both arms)
# G41: SW_SEEDS, SW_NULL_ONLY, SW_ONLY_CASE (defaults keep G40 behaviour)
# G40: null arm (no injection, same draws). arm = "fault" or "null"; the null scores the clean
#      base files against the same hit marker (hit batch drawn but unused), after the fault arm
#      of each iteration, so the fault arm's random stream is unchanged.
#
# Every value below is illustrative; none is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
# G41 switches (defaults keep G40 behaviour): seeds, null arm only, restricted grid.
SEEDS <- eval(parse(text = Sys.getenv("SW_SEEDS", "1:10")))
NULL_ONLY <- Sys.getenv("SW_NULL_ONLY") == "1" # injection draws still run; fault scoring is skipped
ONLY_CASE <- Sys.getenv("SW_ONLY_CASE") # fault:size pairs, comma separated, e.g. batch_shift:0.5; empty = all
N_EVENTS <- 1000
FAULTS <- tibble::tribble(
  ~fault          , ~size ,
  "batch_shift"   , 0.25  , # marker SD units, added in every file of the hit batch
  "batch_shift"   , 0.5   ,
  "batch_stretch" , 1.3   , # width factor around each file's median
  "batch_drift"   , 0.5   , # SD units, rising linearly over the hit batch's acquisition order
  "batch_bimodal" , 0.5   , # SD units: each event moved by +size or -size at random
  "batch_bimodal" , 1
)

# Grid rows kept; i stays the row of FAULTS, so the draws match the full run.
FAULT_ROWS <- if (nzchar(ONLY_CASE)) {
  which(
    paste0(FAULTS$fault, ":", FAULTS$size) %in% strsplit(ONLY_CASE, ",")[[1]]
  )
} else {
  seq_len(nrow(FAULTS))
}

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
donor <- grepl("BMV_[0-9]+[A-Z]", info$file)
base <- ev[donor]
batch <- setNames(info$batch[donor], paste0("f", seq_along(base)))
names(base) <- names(batch)
ch <- setdiff(colnames(base[[1]]), "Time")
set.seed(1) # draw for the per-marker SD used to size the faults
marker_sd <- apply(
  do.call(rbind, lapply(base, \(m) m[sample.int(nrow(m), 500), ch])),
  2,
  sd
)
order_in_batch <- stats::ave(seq_along(batch), batch, FUN = seq_along)
names(order_in_batch) <- names(batch)

auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}

# Scores one set of files for the hit marker among all markers, for the fault and the null arm.
score_arm <- function(files, mk, seed) {
  sc <- sw_batch_scores(
    files,
    batch,
    ch,
    N_EVENTS,
    seed,
    order = order_in_batch,
    reference_group = NULL
  )
  per_marker <- list(
    pair_emd = sc[sc$check == "pair_emd", ] |>
      summarise(v = max(value), .by = marker),
    median_range = sc[sc$check == "median_range", ] |>
      transmute(marker, v = value),
    eta2 = sc[sc$check == "eta2", ] |> transmute(marker, v = value),
    spread_ratio = sc[sc$check == "spread_ratio", ] |>
      summarise(v = max(abs(log(value))), .by = marker),
    order_rho = sc[sc$check == "order_rho", ] |>
      summarise(v = max(abs(value)), .by = marker)
  )
  tibble(
    check = names(per_marker),
    auroc = vapply(
      per_marker,
      \(d) auroc(d$v[match(ch, d$marker)], ch == mk),
      numeric(1)
    )
  )
}

long <- lapply(SEEDS, function(seed) {
  lapply(FAULT_ROWS, function(i) {
    fault <- FAULTS$fault[[i]]
    size <- FAULTS$size[[i]]
    set.seed(seed * 1000 + i)
    hit_batch <- sample(unique(batch), 1)
    mk <- sample(ch, 1)
    hit_files <- names(batch)[batch == hit_batch]
    files <- base
    for (f in hit_files) {
      x <- files[[f]]
      x[, mk] <- switch(
        fault,
        batch_shift = x[, mk] + size * marker_sd[[mk]],
        batch_stretch = stats::median(x[, mk]) +
          (x[, mk] - stats::median(x[, mk])) * size,
        batch_drift = x[, mk] +
          size *
            marker_sd[[mk]] *
            (order_in_batch[[f]] - 1) /
            (length(hit_files) - 1),
        batch_bimodal = x[, mk] +
          sample(c(-size, size), nrow(x), replace = TRUE) * marker_sd[[mk]]
      )
      files[[f]] <- x
    }
    # The injection above and its draws always run; only the fault scoring is skipped.
    fault_arm <- if (NULL_ONLY) NULL else score_arm(files, mk, seed)
    # G40 null arm: same hit marker, clean files. It runs last in the iteration (sw_batch_scores
    # may set the seed); the next iteration starts with its own set.seed().
    null_arm <- score_arm(base, mk, seed)
    bind_rows(
      if (!is.null(fault_arm)) {
        mutate(
          fault_arm,
          arm = "fault",
          seed = seed,
          fault = fault,
          size = size,
          .before = 1
        )
      },
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
  file.path(run_dir, "batch_injected_long.csv"),
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
  file.path(run_dir, "figure_source_data_batch_real.csv"),
  row.names = FALSE
)
print(as.data.frame(mutate(mat, auroc = round(auroc, 2))[c(
  "arm",
  "case",
  "check",
  "auroc"
)]))
if (!NULL_ONLY) {
  p <- ggplot(filter(mat, arm == "fault"), aes(auroc, check)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
    geom_pointrange(
      aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)),
      size = 0.2
    ) +
    facet_wrap(vars(case), nrow = 2) +
    labs(
      x = "AUROC of the hit marker among 40 (mean, 95% CI over 10 seeds; real batches)",
      y = NULL
    )
  ggsave(
    file.path(run_dir, "fig_batch_real.png"),
    p,
    width = 200,
    height = 60,
    units = "mm",
    dpi = 200
  )
}
