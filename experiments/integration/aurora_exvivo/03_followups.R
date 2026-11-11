# experiments/integration/aurora_exvivo/03_followups.R
# Purpose: two follow-ups to 02_diagnose.R on the real files.
#          (1) Size of the largest time gap in the file with the largest max_gap, against
#              all other files, as a share of each file's span.
#          (2) Exact EMD per marker between reference files of different batches, against
#              donor-vs-donor EMD within and across batches. Reference files are those whose
#              name contains "BMV_REF". If they are the same material, their difference is a
#              batch effect free of donor effects; the script cannot check that.
# Run:     Rscript experiments/integration/aurora_exvivo/03_followups.R <run_dir>
# Inputs:  <run_dir>/events.rds, files.csv   from 01_extract.R
# Outputs: <run_dir>/figure_source_data_ref_emd.csv  EMD per pair and marker
#          <run_dir>/fig_ref_vs_donor_emd.png         median EMD per marker and pair type
#
# Pair counts (40 same-batch pairs per batch, 160 cross-batch pairs) and the seed are
# illustrative; they only set how well the donor yardstick is estimated.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
PAIRS_PER_BATCH <- 40
CROSS_PAIRS <- 160
SEED <- 1

ev <- readRDS(file.path(run_dir, "events.rds"))
info <- read.csv(file.path(run_dir, "files.csv"))
stopifnot("files match events" = identical(info$file, names(ev)))
ch <- setdiff(colnames(ev[[1]]), "Time")

# (1) Largest time gap as a share of the file span ----------------------------
gap_share <- vapply(
  ev,
  \(m) {
    t <- sort(m[, "Time"])
    max(diff(t)) / diff(range(t))
  },
  numeric(1)
)
top_gap <- names(which.max(gap_share))
cat(
  "largest gap share:",
  round(max(gap_share), 3),
  "; next largest:",
  round(sort(gap_share, decreasing = TRUE)[2], 3),
  "; 99th percentile:",
  round(quantile(gap_share, 0.99), 3),
  "\n"
)

# (2) Reference vs donor EMD ---------------------------------------------------
refs <- info$file[grepl("BMV_REF", info$file)]
donors <- info[grepl("BMV_[0-9]+[A-Z]", info$file), ]
pair_emd <- function(file_a, file_b, pair_kind) {
  # Computed before the table is built, so no argument shares a column's name.
  emd_values <- vapply(
    ch,
    \(m) sw_emd_1d(ev[[file_a]][, m], ev[[file_b]][, m]),
    numeric(1)
  )
  tibble(
    a = file_a,
    b = file_b,
    marker = ch,
    kind = pair_kind,
    emd = emd_values
  )
}
set.seed(SEED)
ref_pairs <- combn(refs, 2, simplify = FALSE)
same_pairs <- lapply(split(donors$file, donors$batch), \(fs) {
  replicate(PAIRS_PER_BATCH, sample(fs, 2), simplify = FALSE)
}) |>
  unlist(recursive = FALSE) |>
  unname() # names such as "Batch11" would turn the pair tables into packed columns
cross_pairs <- replicate(
  CROSS_PAIRS,
  {
    b <- sample(unique(donors$batch), 2)
    c(
      sample(donors$file[donors$batch == b[1]], 1),
      sample(donors$file[donors$batch == b[2]], 1)
    )
  },
  simplify = FALSE
)
all_emd <- bind_rows(
  lapply(ref_pairs, \(p) {
    pair_emd(p[1], p[2], "reference vs reference (different batches)")
  }),
  lapply(same_pairs, \(p) pair_emd(p[1], p[2], "donor vs donor (same batch)")),
  lapply(cross_pairs, \(p) {
    pair_emd(p[1], p[2], "donor vs donor (different batches)")
  })
)
write.csv(
  all_emd,
  file.path(run_dir, "figure_source_data_ref_emd.csv"),
  row.names = FALSE
)

med <- all_emd |>
  summarise(median = median(emd), .by = c(marker, kind)) |>
  pivot_wider(names_from = kind, values_from = median)
ref_col <- "reference vs reference (different batches)"
same_col <- "donor vs donor (same batch)"
cat(
  "markers where reference-vs-reference EMD exceeds the same-batch donor median:",
  sum(med[[ref_col]] > med[[same_col]]),
  "of",
  nrow(med),
  "\n"
)

fig <- all_emd |>
  mutate(marker = factor(marker, levels = med$marker[order(med[[ref_col]])])) |>
  ggplot(aes(emd, marker, colour = kind)) +
  stat_summary(fun = median, geom = "point", size = 1.3) +
  labs(
    x = "Median exact EMD between two files (marker scale)",
    y = NULL,
    colour = NULL
  ) +
  theme(legend.position = "top", legend.direction = "vertical")
ggsave(
  file.path(run_dir, "fig_ref_vs_donor_emd.png"),
  fig,
  width = 110,
  height = 150,
  units = "mm",
  dpi = 200
)
