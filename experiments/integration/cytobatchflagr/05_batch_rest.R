# experiments/integration/cytobatchflagr/05_batch_rest.R
# Purpose: does `rest_emd_z` (EMD of one batch to all other batches pooled, robust z over batches
#          within a marker) find a batch fault better than the pairwise `pair_emd` maximum? Same
#          real files as 04: four control samples acquired in 18 batches. One marker is split into
#          two modes (+/- size SD) in every file of one batch. Each marker gets one score (the
#          maximum over batches of rest_emd_z, or over batch pairs of pair_emd); AUROC is the rank
#          of the hit marker among all markers. The null arm scores the same files with no fault
#          injected against the same randomly drawn hit marker (as the null in 07); its AUROC
#          should sit near 0.5.
#          Decisions: G28 (new check), G29 (AUROC unit as in 04), G30 (sizes, seeds), G31 (stop rule).
# Run:     Rscript experiments/integration/cytobatchflagr/05_batch_rest.R <extract_run_dir> <run_dir> <seeds>
#          <seeds> is a comma list or range, for example 1 (timing) or 1:20. <run_dir> must not exist.
# Inputs:  <extract_run_dir>/events.rds, files.csv, marker_names.rds   from 01_extract.R
#                (experiments/runs/2026-11-10_cytobatchflagr/)
# Outputs: <run_dir>/batch_rest_long.csv                one row per seed, size, check, arm
#          <run_dir>/figure_source_data_batch_rest.csv  mean AUROC and 95% CI per size, check, arm
#          <run_dir>/fig_batch_rest.png
#          <run_dir>/session_info.txt
#
# Choices reused from 04_reference_group.R (approved by the user): 25 shared named markers (G9 says 30; the code keeps the 25 named ones),
# asinh(x / 6000), 2000 events per file for the batch scores, the batch fault and its injection
# lines, the seed rule set.seed(seed * 1000 + g), the AUROC function and the CI formula
# (1.96 * SD / sqrt(n seeds)). Marker SD for the fault size is taken over all files.
# The bar from the user: mean AUROC >= 0.90 and lower CI >= 0.80.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
extract_dir <- args[[1]]
run_dir <- args[[2]]
SEEDS <- eval(parse(text = paste0("c(", args[[3]], ")"))) # G30: 20 seeds after the timing run
if (dir.exists(run_dir)) {
  stop("run dir exists: ", run_dir, call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# Settings --------------------------------------------------------------------------------------
SIZES <- c(0.5, 1) # two-mode split, in marker SD (G30)
COFACTOR <- 6000 # arcsinh scale, as in 04
BATCH_EVENTS <- 2000 # events per file for the batch scores, as in 04

# Data ------------------------------------------------------------------------------------------

ev <- readRDS(file.path(extract_dir, "events.rds"))
info <- read.csv(file.path(extract_dir, "files.csv"))
marker_names <- readRDS(file.path(extract_dir, "marker_names.rds"))
stopifnot("files match events" = identical(info$file, names(ev)))
ids <- paste0("f", seq_along(ev))
names(ev) <- ids
shared <- Reduce(intersect, lapply(ev, colnames))
markers <- intersect(names(marker_names)[!is.na(marker_names)], shared)
base <- lapply(ev, \(m) asinh(m[, markers] / COFACTOR))
batch <- setNames(as.character(info$batch), ids)
marker_sd <- apply(do.call(rbind, base), 2, sd)

auroc <- function(score, pos) {
  ok <- !is.na(score)
  r <- rank(score[ok])
  p <- pos[ok]
  n1 <- sum(p)
  (sum(r[p]) - n1 * (n1 + 1) / 2) / (n1 * sum(!p))
}

# Fault and scores ------------------------------------------------------------------------------

long <- lapply(SEEDS, function(seed) {
  lapply(seq_along(SIZES), function(g) {
    size <- SIZES[g]
    set.seed(seed * 1000 + g)
    mk <- sample(markers, 1) # hit marker
    hit_batch <- sample(unique(batch), 1)
    files <- base
    bad <- batch == hit_batch
    files[bad] <- lapply(files[bad], function(x) {
      x[, mk] <- x[, mk] +
        sample(c(-size, size), nrow(x), replace = TRUE) * marker_sd[[mk]]
      x
    })
    # Null: the clean files, scored against the same randomly drawn hit marker.
    arms <- list(true = files, null = base)
    lapply(names(arms), function(arm) {
      sc <- sw_batch_scores(
        arms[[arm]],
        batch,
        markers,
        BATCH_EVENTS,
        seed,
        order = NULL,
        reference_group = NULL
      )
      per_marker <- list(
        pair_emd = sc[sc$check == "pair_emd", ] |>
          summarise(v = max(value), .by = marker),
        rest_emd_z = sc[sc$check == "rest_emd_z", ] |>
          summarise(v = max(value), .by = marker)
      )
      res <- vapply(
        per_marker,
        \(d) auroc(d$v[match(markers, d$marker)], markers == mk),
        numeric(1)
      )
      tibble(
        seed = seed,
        size = size,
        check = names(res),
        arm = arm,
        auroc = unname(res)
      )
    }) |>
      bind_rows()
  }) |>
    bind_rows()
}) |>
  bind_rows()
write.csv(long, file.path(run_dir, "batch_rest_long.csv"), row.names = FALSE)

# Summary ---------------------------------------------------------------------------------------

mat <- long |>
  summarise(
    ci95_half = 1.96 * sd(auroc, na.rm = TRUE) / sqrt(sum(!is.na(auroc))),
    auroc = mean(auroc, na.rm = TRUE),
    n_seeds = n(),
    .by = c(size, check, arm)
  ) |>
  mutate(
    case = paste("batch_split", size, "SD"),
    arm = factor(arm, levels = c("true", "null")),
    meets_bar = auroc >= 0.9 & auroc - ci95_half >= 0.8
  )
write.csv(
  mat,
  file.path(run_dir, "figure_source_data_batch_rest.csv"),
  row.names = FALSE
)
print(as.data.frame(mutate(
  mat,
  across(c(auroc, ci95_half), \(x) round(x, 3))
)[c("case", "check", "arm", "auroc", "ci95_half", "meets_bar")]))

p <- ggplot(mat, aes(auroc, check, colour = arm)) +
  geom_vline(xintercept = c(0.5, 0.9), linetype = 2, colour = "grey60") +
  geom_point(
    data = mutate(
      long,
      case = paste("batch_split", size, "SD"),
      arm = factor(arm, levels = c("true", "null"))
    ),
    position = position_dodge(width = 0.6),
    size = 0.4,
    alpha = 0.5
  ) +
  geom_pointrange(
    aes(xmin = auroc - ci95_half, xmax = pmin(auroc + ci95_half, 1)),
    position = position_dodge(width = 0.6),
    size = 0.2
  ) +
  facet_wrap(vars(case), ncol = 1) +
  scale_colour_manual(values = c(true = "#D55E00", null = "grey40")) +
  labs(
    x = sprintf(
      "AUROC of the hit marker among %d markers (mean and 95%% CI over %d seeds; dots = seeds)",
      length(markers),
      length(SEEDS)
    ),
    y = NULL,
    colour = "batch labels"
  )
ggsave(
  file.path(run_dir, "fig_batch_rest.png"),
  p,
  width = 120,
  height = 90,
  units = "mm",
  dpi = 200
)
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
