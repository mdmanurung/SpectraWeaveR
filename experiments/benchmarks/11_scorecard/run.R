# experiments/benchmarks/11_scorecard/run.R
# Purpose: one scorecard of the real-file evidence. For each fault class and dataset, read the
#          designed check's AUROC at the agreed size from saved source data (no reruns), apply the
#          user's bar (D1) and give a status: met, met (null untested), not met, or observed on
#          real data with no AUROC. The best check at the same size is shown beside it.
# Run:     Rscript experiments/benchmarks/11_scorecard/run.R <run_dir>
# Inputs:  experiments/benchmarks/11_scorecard/scorecard_map.csv   one row per class and dataset
#          every source file the map names (reference/figures/*, experiments/runs/*)
# Outputs: <run_dir>/scorecard.csv                 class, dataset, size, check, AUROC, CI, status
#          <run_dir>/fig_scorecard.png, figure_source_data_scorecard.csv
#          <run_dir>/session_info.txt
#
# Decisions: G39 (this scorecard), D1 (bar), D2 (sizes).

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
if (dir.exists(run_dir)) {
  stop("run dir exists: ", run_dir, call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# The bar (D1) ------------------------------------------------------------------------------------
BAR_MEAN <- 0.90
BAR_LOWER <- 0.80
BAR_SEEDS <- 5
BAR_NULL <- 0.55

map <- read.csv(
  "experiments/benchmarks/11_scorecard/scorecard_map.csv",
  na.strings = ""
)

# One map row -> its designed row, the best row in its group and the null ------------------------
pick <- function(d, expr) d[eval(parse(text = expr), d), , drop = FALSE]
read_row <- function(i) {
  m <- map[i, ]
  if (!is.na(m$status_override)) {
    return(tibble(i = i, check = NA_character_))
  }
  d <- read.csv(m$source)
  # The check's name sits in `check`, `score` (05_marker_level) or `arm` (06_frozen_reference).
  check_col <- intersect(c("check", "score", "arm"), names(d))[1]
  hit <- pick(d, m$row_filter)
  if (nrow(hit) != 1) {
    stop(
      sprintf(
        "map row %d (%s, %s) matched %d rows",
        i,
        m$class,
        m$dataset,
        nrow(hit)
      ),
      call. = FALSE
    )
  }
  grp <- pick(d, m$group_filter)
  best <- grp[which.max(grp[[m$auroc_col]]), ]
  null_auroc <- NA_real_
  if (!is.na(m$null_filter)) {
    nul <- pick(d, m$null_filter)
    stopifnot("one null row" = nrow(nul) == 1)
    null_auroc <- nul[[m$auroc_col]]
  }
  # Seeds: compare with a fixed SEEDS line in the producing script when there is one.
  src <- readLines(m$seeds_from)
  fixed <- regmatches(src, regexpr("SEEDS <- 1:[0-9]+", src))
  seeds_checked <- if (length(fixed) == 1) {
    stopifnot(
      "seeds match the script" = as.integer(sub(".*1:", "", fixed)) == m$n_seeds
    )
    "script"
  } else {
    "run argument"
  }
  tibble(
    i = i,
    check = hit[[check_col]],
    auroc = hit[[m$auroc_col]],
    ci95_half = hit$ci95_half,
    best_check = best[[check_col]],
    best_auroc = best[[m$auroc_col]],
    null_auroc = null_auroc,
    seeds_checked = seeds_checked
  )
}
rows <- bind_rows(lapply(seq_len(nrow(map)), read_row))

# Status per class and dataset -------------------------------------------------------------------
scorecard <- bind_cols(
  map[c("class", "dataset", "size", "n_seeds", "status_override", "note")],
  rows[-1]
) |>
  mutate(
    lower = auroc - ci95_half,
    status = case_when(
      !is.na(status_override) ~ "observed on real data, no AUROC",
      auroc >= BAR_MEAN &
        lower >= BAR_LOWER &
        n_seeds >= BAR_SEEDS &
        !is.na(null_auroc) &
        null_auroc <= BAR_NULL ~ "met",
      auroc >= BAR_MEAN &
        lower >= BAR_LOWER &
        n_seeds >= BAR_SEEDS &
        is.na(null_auroc) ~ "met (null untested)",
      .default = "not met"
    )
  ) |>
  select(
    class,
    dataset,
    size,
    check,
    auroc,
    ci95_half,
    lower,
    n_seeds,
    seeds_checked,
    null_auroc,
    best_check,
    best_auroc,
    status,
    note
  )
write.csv(scorecard, file.path(run_dir, "scorecard.csv"), row.names = FALSE)
print(as.data.frame(
  scorecard |>
    mutate(across(c(auroc, lower, best_auroc, null_auroc), \(x) round(x, 2))) |>
    select(class, dataset, check, auroc, lower, null_auroc, status)
))
print(count(scorecard, status))

# Figure ------------------------------------------------------------------------------------------
status_levels <- c("met", "met (null untested)", "not met")
plot_frame <- scorecard |>
  filter(status %in% status_levels) |>
  mutate(
    row = paste0(class, " (", dataset, ", ", size, ")"),
    row = reorder(row, auroc),
    status = factor(status, levels = status_levels)
  )
write.csv(
  plot_frame,
  file.path(run_dir, "figure_source_data_scorecard.csv"),
  row.names = FALSE
)
p <- ggplot(plot_frame, aes(auroc, row, colour = status)) +
  geom_vline(
    xintercept = c(0.5, BAR_LOWER, BAR_MEAN),
    linetype = c(3, 2, 2),
    colour = "grey60"
  ) +
  geom_pointrange(
    aes(xmin = pmax(auroc - ci95_half, 0), xmax = pmin(auroc + ci95_half, 1)),
    size = 0.25
  ) +
  scale_colour_manual(
    values = c(
      met = "#009E73",
      `met (null untested)` = "#56B4E9",
      `not met` = "#D55E00"
    )
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    x = "AUROC, designed check, real files (mean, 95% CI)",
    y = NULL,
    colour = NULL,
    caption = "Bar (D1): mean >= 0.90, lower CI >= 0.80, >= 5 seeds, null <= 0.55. Dashed: 0.80 and 0.90; dotted: chance"
  )
ggsave(
  file.path(run_dir, "fig_scorecard.png"),
  p,
  width = 180,
  height = 150,
  units = "mm",
  dpi = 200
)

writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
