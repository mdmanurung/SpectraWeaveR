# experiments/integration/glehr2026/12_af700_cd8dim.R
# Purpose: test whether CD8-dim cells widen the Fortessa AF700 (CD8) panel negative (0.53 against
#          0.13 for the AF700 reagent alone, 11). The authors' panel files hold CD3+ cells only
#          (lympho-single_cd3 gate), so subsets come from CD4 alone. CD4+ T cells are truly
#          CD8-negative: if their
#          AF700 SD is far below the panel negative's, the width comes from other cells the negative
#          takes in; if near it, the true negatives are wide themselves (then summed spread from
#          CD4-APC and CD3-AA750 is the likelier cause).
# Run:     Rscript experiments/integration/glehr2026/12_af700_cd8dim.R <data_dir> <run_dir>
# Inputs:  <data_dir>/BAD/<instrument>/<donor>/*relativised*12-panel.fcs  authors' gated panel
#          experiments/runs/2026-11-10_glehr-af700-pb-negative/negative_audit.csv  panel neg_sd, 08
#          experiments/runs/2026-11-10_glehr-unstained-noise-v2/unstained_sd.csv   unstained SD, 10
#          experiments/runs/2026-11-10_glehr-af700-stain/single_stain_sd.csv      AF700-stain SD, 11
# Outputs: <run_dir>/subset_sd.csv              donor, instrument, subset, n, AF700 SD
#          <run_dir>/fig_af700_subset_sd.png, figure_source_data_af700_subset_sd.csv
#          <run_dir>/fig_af700_subset_density.png, figure_source_data_af700_subset_density.csv
#          <run_dir>/session_info.txt
#
# Decisions G16, G18 (plans/2026-11-10_diagnostics-next.md).

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[[1]]
run_dir <- args[[2]]
if (dir.exists(run_dir)) {
  stop("run dir exists: ", run_dir, call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

TARGET <- "AF700-A" # CD8
CD4 <- "APC-A"
SPLIT_SEED <- 1 # G18 (revised: CD4 split only)
SPLIT_NSTART <- 10
MIN_EVENTS <- 200 # a subset with fewer events gets NA
LATER_SESSION <- paste0("d", 6:9) # Aurora donors 6-9 acquired later with other voltages
INSTRUMENTS <- c("aurora", "fortessa")
SUBSETS <- c("CD4+ T", "CD4- T")

robust_sd <- function(v) {
  if (length(v) < MIN_EVENTS) NA_real_ else 1.4826 * mad(v, constant = 1)
}
# k-means, k = 2; TRUE for the cluster with the higher centre.
split_high <- function(v) {
  set.seed(SPLIT_SEED)
  km <- kmeans(v, centers = 2, nstart = SPLIT_NSTART)
  km$cluster == which.max(km$centers[, 1])
}

# Subsets and AF700 SD per donor and instrument ----------------------------------------------
donor_dirs <- sort(basename(list.dirs(
  file.path(data_dir, "BAD", "aurora"),
  recursive = FALSE
)))
grid <- expand.grid(
  donor_i = seq_along(donor_dirs),
  instrument = INSTRUMENTS,
  stringsAsFactors = FALSE
)
per_file <- lapply(seq_len(nrow(grid)), function(j) {
  dir <- file.path(
    data_dir,
    "BAD",
    grid$instrument[j],
    donor_dirs[grid$donor_i[j]]
  )
  f <- list.files(
    dir,
    pattern = "relativised.*12-panel\\.fcs$",
    full.names = TRUE
  )
  e <- exprs(read.FCS(
    f,
    truncate_max_range = FALSE,
    transformation = FALSE,
    emptyValue = FALSE
  ))
  x <- asinh(e[, c(TARGET, CD4)] / 1000)
  tibble(
    donor = paste0("d", grid$donor_i[j]),
    instrument = grid$instrument[j],
    subset = if_else(split_high(x[, CD4]), "CD4+ T", "CD4- T"),
    af700 = x[, TARGET],
    cd4 = x[, CD4]
  )
})
events <- bind_rows(per_file)
subset_sd <- events |>
  summarise(n = n(), sd = robust_sd(af700), .by = c(donor, instrument, subset))
write.csv(subset_sd, file.path(run_dir, "subset_sd.csv"), row.names = FALSE)
print(as.data.frame(
  subset_sd |>
    summarise(n_min = min(n), n_median = median(n), .by = c(instrument, subset))
))

# AF700 SD per sample type and subset -----------------------------------------------------------
# Comparators from 10 (unstained), 11 (AF700 stain, CD3- cells) and 08 (panel negative).
comparators <- bind_rows(
  read.csv(
    "experiments/runs/2026-11-10_glehr-unstained-noise-v2/unstained_sd.csv"
  ) |>
    filter(channel == TARGET) |>
    transmute(donor, instrument, sample_type = "unstained", sd = unstained_sd),
  read.csv(
    "experiments/runs/2026-11-10_glehr-af700-stain/single_stain_sd.csv"
  ) |>
    filter(stain == "AF700", side == "CD3-") |>
    transmute(
      donor,
      instrument,
      sample_type = "AF700 stain, CD3-",
      sd = af700_sd
    ),
  read.csv(
    "experiments/runs/2026-11-10_glehr-af700-pb-negative/negative_audit.csv"
  ) |>
    filter(check == "neg_sd", channel == TARGET) |>
    separate_wider_delim(file, "_", names = c("donor", "instrument")) |>
    transmute(donor, instrument, sample_type = "panel negative", sd = value)
)
sd_frame <- bind_rows(
  comparators,
  subset_sd |>
    filter(subset %in% SUBSETS) |>
    transmute(donor, instrument, sample_type = paste("panel", subset), sd)
) |>
  mutate(
    sample_type = factor(
      sample_type,
      levels = c(
        "unstained",
        "AF700 stain, CD3-",
        "panel CD4+ T",
        "panel CD4- T",
        "panel negative"
      )
    ),
    session = if_else(
      donor %in% LATER_SESSION,
      "donors d6-d9 (later Aurora session)",
      "d1-d5"
    )
  )
stopifnot(
  "one SD per donor, instrument and sample type" = nrow(sd_frame) == 9 * 2 * 5
)
write.csv(
  sd_frame,
  file.path(run_dir, "figure_source_data_af700_subset_sd.csv"),
  row.names = FALSE
)
print(as.data.frame(
  sd_frame |>
    summarise(
      median_sd = round(median(sd, na.rm = TRUE), 3),
      min = round(min(sd, na.rm = TRUE), 3),
      max = round(max(sd, na.rm = TRUE), 3),
      .by = c(instrument, sample_type)
    ) |>
    arrange(instrument, sample_type)
))
p_sd <- ggplot(sd_frame, aes(sample_type, sd, group = donor)) +
  geom_line(alpha = 0.4) +
  geom_point(aes(shape = session), size = 1.2) +
  facet_wrap(vars(instrument)) +
  scale_y_log10() +
  labs(
    x = NULL,
    y = "AF700 (CD8) robust SD, asinh(x / 1000) of relativised values",
    caption = "Lines: donors; robust SD = 1.4826 x MAD; panel files are CD3+; subsets split by k-means on CD4"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(
  file.path(run_dir, "fig_af700_subset_sd.png"),
  p_sd,
  width = 170,
  height = 100,
  units = "mm",
  dpi = 200
)

# AF700 density per subset ------------------------------------------------------------------------
density_frame <- events |>
  filter(subset %in% SUBSETS) |>
  mutate(subset = factor(subset, levels = SUBSETS))
set.seed(1) # display subsample only
density_frame <- slice_sample(
  density_frame,
  n = 30000,
  by = c(instrument, subset)
)
write.csv(
  density_frame,
  file.path(run_dir, "figure_source_data_af700_subset_density.csv"),
  row.names = FALSE
)
p_den <- ggplot(density_frame, aes(af700, colour = instrument)) +
  geom_density(adjust = 0.7) +
  facet_wrap(vars(subset), ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-2, 5)) +
  labs(
    x = "AF700 (CD8), asinh(x / 1000) of relativised values",
    y = "Density",
    caption = "Donors pooled; up to 30,000 events per instrument and subset"
  )
ggsave(
  file.path(run_dir, "fig_af700_subset_density.png"),
  p_den,
  width = 120,
  height = 120,
  units = "mm",
  dpi = 200
)

writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
