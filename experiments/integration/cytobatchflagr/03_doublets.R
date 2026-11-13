# experiments/integration/cytobatchflagr/03_doublets.R
# Purpose: real-data test of the doublet check (area_height_z in sw_fcs_signal_audit). The
#          CytoBatchFlagR files are pre-gated live singlets (CD45+CD14-), 72 files in 18 batches.
#          Doublets are injected into 25% of files: an injected event is two events passing
#          together, so every area channel is their sum and every height channel their maximum.
#          AUROC ranks faulty files above clean ones; the real batch effects stay in.
# G40: null arm (no injection, same draws)
# Run:     Rscript experiments/integration/cytobatchflagr/03_doublets.R <extract_run_dir> <run_dir>
# Inputs:  <extract_run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/doublets_long.csv, figure_source_data_doublets.csv   (both have a column
#          `arm`: "fault" and "null"; null = same bad labels, clean files, same scoring)
#          <run_dir>/fig_doublets.png   (arm == "fault" only), session_info.txt
#
# Fault sizes and seeds are illustrative; the two quantiles are the variants the user approved.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
extract_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) {
  stop("run dir exists: ", run_dir, call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

SEEDS <- 1:5
FRACTION_BAD <- 0.25
DOUBLET_SHARE <- c(0.02, 0.05, 0.10)
QUANTILES <- c(0.95, 0.99)
AREA_HEIGHT <- c(`FSC-A` = "FSC-H")

ev <- readRDS(file.path(extract_dir, "events.rds"))
names(ev) <- paste0("f", seq_along(ev))
shared <- Reduce(intersect, lapply(ev, colnames))
area_ch <- grep("-A$", shared, value = TRUE)
height_ch <- grep("-H$", shared, value = TRUE)
markers <- setdiff(area_ch, c("FSC-A", "SSC-A", "SSC-B-A", "AF-A"))

auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}
add_doublets <- function(x, share) {
  k <- sample.int(nrow(x), round(share * nrow(x)))
  partner <- sample.int(nrow(x), length(k))
  x[k, area_ch] <- x[k, area_ch] + x[partner, area_ch]
  x[k, height_ch] <- pmax(x[k, height_ch], x[partner, height_ch])
  x
}

long <- lapply(SEEDS, function(seed) {
  lapply(DOUBLET_SHARE, function(share) {
    set.seed(seed * 1000 + round(share * 100))
    bad <- seq_along(ev) %in%
      sample(seq_along(ev), round(FRACTION_BAD * length(ev)))
    files <- ev
    files[bad] <- lapply(files[bad], add_doublets, share = share)
    # G40: null arm = the same clean files and bad labels, no doublets. lapply runs "fault" for
    # every quantile before "null", and sw_fcs_signal_audit draws no random numbers.
    arms <- list(fault = files, null = ev)
    lapply(names(arms), function(which_arm) {
      files <- arms[[which_arm]]
      lapply(QUANTILES, function(q) {
        sa <- sw_fcs_signal_audit(
          files,
          markers,
          NULL,
          NULL,
          time_channel = NULL,
          time_bin = NULL,
          area_height = AREA_HEIGHT,
          area_height_quantile = q,
          reference_group = NULL,
          reference = NULL,
          summary = NULL
        )
        z <- sa[sa$check == "area_height_z", ]
        fsc_sd <- sa[sa$check == "sd_z" & sa$channel == markers[1], ]
        gz <- sa[sa$check == "global_sd_z", ]
        tibble(
          seed = seed,
          arm = which_arm,
          share = share,
          quantile = q,
          check = c("area_height_z", "global_sd_z", "sd_z_one_marker"),
          auroc = c(
            auroc(z$value[match(names(files), z$file)], bad),
            auroc(abs(gz$value[match(names(files), gz$file)]), bad),
            auroc(abs(fsc_sd$value[match(names(files), fsc_sd$file)]), bad)
          )
        )
      }) |>
        bind_rows()
    }) |>
      bind_rows()
  }) |>
    bind_rows()
}) |>
  bind_rows()
write.csv(long, file.path(run_dir, "doublets_long.csv"), row.names = FALSE)

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc) / sqrt(n()),
    auroc = mean(auroc),
    .by = c(arm, share, quantile, check)
  )
write.csv(
  mat,
  file.path(run_dir, "figure_source_data_doublets.csv"),
  row.names = FALSE
)
print(as.data.frame(mutate(mat, across(c(auroc, ci95_half), \(x) round(x, 3)))))
p <- ggplot(filter(mat, arm == "fault"), aes(share, auroc, colour = check)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = 2, colour = "grey60") +
  geom_pointrange(
    aes(ymin = auroc - ci95_half, ymax = pmin(auroc + ci95_half, 1)),
    size = 0.2
  ) +
  geom_line() +
  facet_wrap(vars(quantile), labeller = label_both) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    x = "Injected doublets (share of events)",
    y = "AUROC, 72 real files (mean, 95% CI, 5 seeds)",
    colour = NULL
  )
ggsave(
  file.path(run_dir, "fig_doublets.png"),
  p,
  width = 150,
  height = 70,
  units = "mm",
  dpi = 200
)
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
