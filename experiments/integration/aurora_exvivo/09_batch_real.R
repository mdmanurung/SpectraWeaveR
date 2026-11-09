# experiments/integration/aurora_exvivo/09_batch_real.R
# Purpose: inject a batch-level fault (one marker in every donor file of one batch) into the real
#          multi-batch files and ask whether sw_batch_scores() ranks the hit marker above the 39
#          others. The real batches already differ (batch nested with donor), so the injected
#          marker must stand out above real batch differences.
# Run:     Rscript experiments/integration/aurora_exvivo/09_batch_real.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/batch_injected_long.csv
#          <run_dir>/fig_batch_real.png, figure_source_data_batch_real.csv
#
# Every value below is illustrative; none is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
SEEDS <- 1:10
N_EVENTS <- 1000
FAULTS <- tibble::tribble(
  ~fault, ~size,
  "batch_shift", 0.25, # marker SD units, added in every file of the hit batch
  "batch_shift", 0.5,
  "batch_stretch", 1.3, # width factor around each file's median
  "batch_drift", 0.5 # SD units, rising linearly over the hit batch's acquisition order
)

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
donor <- grepl("BMV_[0-9]+[A-Z]", info$file)
base <- ev[donor]
batch <- setNames(info$batch[donor], paste0("f", seq_along(base)))
names(base) <- names(batch)
ch <- setdiff(colnames(base[[1]]), "Time")
set.seed(1) # draw for the per-marker SD used to size the faults
marker_sd <- apply(do.call(rbind, lapply(base, \(m) m[sample.int(nrow(m), 500), ch])), 2, sd)
order_in_batch <- stats::ave(seq_along(batch), batch, FUN = seq_along)
names(order_in_batch) <- names(batch)

auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}

long <- lapply(SEEDS, function(seed) {
  lapply(seq_len(nrow(FAULTS)), function(i) {
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
        batch_stretch = stats::median(x[, mk]) + (x[, mk] - stats::median(x[, mk])) * size,
        batch_drift = x[, mk] +
          size * marker_sd[[mk]] * (order_in_batch[[f]] - 1) / (length(hit_files) - 1)
      )
      files[[f]] <- x
    }
    sc <- sw_batch_scores(files, batch, ch, N_EVENTS, seed, order = order_in_batch)
    per_marker <- list(
      pair_emd = sc[sc$check == "pair_emd", ] |> summarise(v = max(value), .by = marker),
      median_range = sc[sc$check == "median_range", ] |> transmute(marker, v = value),
      eta2 = sc[sc$check == "eta2", ] |> transmute(marker, v = value),
      spread_ratio = sc[sc$check == "spread_ratio", ] |>
        summarise(v = max(abs(log(value))), .by = marker),
      order_rho = sc[sc$check == "order_rho", ] |> summarise(v = max(abs(value)), .by = marker)
    )
    tibble(
      seed = seed,
      fault = fault,
      size = size,
      check = names(per_marker),
      auroc = vapply(per_marker, \(d) auroc(d$v[match(ch, d$marker)], ch == mk), numeric(1))
    )
  })
}) |>
  unlist(recursive = FALSE) |>
  unname() |>
  bind_rows()
write.csv(long, file.path(run_dir, "batch_injected_long.csv"), row.names = FALSE)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    .by = c(fault, size, check)
  ) |>
  mutate(case = paste(fault, size))
write.csv(mat, file.path(run_dir, "figure_source_data_batch_real.csv"), row.names = FALSE)
print(as.data.frame(mutate(mat, auroc = round(auroc, 2))[c("case", "check", "auroc")]))
p <- ggplot(mat, aes(auroc, check)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)), size = 0.2) +
  facet_wrap(vars(case), nrow = 1) +
  labs(x = "AUROC of the hit marker among 40 (mean, 95% CI over 10 seeds; real batches)", y = NULL)
ggsave(file.path(run_dir, "fig_batch_real.png"), p, width = 200, height = 60, units = "mm", dpi = 200)
