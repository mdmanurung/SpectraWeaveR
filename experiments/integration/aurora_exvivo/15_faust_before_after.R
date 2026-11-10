# experiments/integration/aurora_exvivo/15_faust_before_after.R
# Purpose: annotate the same events with FAUST before and after cyCombine batch correction and ask
#          whether correction kept the biology (label agreement, phenotype shares by batch). No
#          fault is injected. Donors are nested in batches; the four reference files (one per
#          batch) are the only unconfounded batch contrast and are drawn apart.
# Run:     Rscript experiments/integration/aurora_exvivo/15_faust_before_after.R <events_dir> <run_dir>
#          [<labels_dir>]
#          <run_dir> must not exist. With <labels_dir>, FAUST is not run: its labels_before.rds and
#          labels_after.rds are read, and the correction, compare and figures run again.
# Inputs:  <events_dir>/events.rds, files.csv, marker_names.rds   from 01_extract.R
#          <labels_dir>/labels_before.rds, labels_after.rds      from an earlier run (optional)
# Outputs: <run_dir>/faust_compare_long.csv          sw_faust_compare() table
#          <run_dir>/labels_before.rds, labels_after.rds   sw_faust_labels() results
#          <run_dir>/fig_faust_prop_by_batch.png, figure_source_data_faust_prop_by_batch.csv
#          <run_dir>/fig_faust_ari.png, figure_source_data_faust_ari.csv
#          <run_dir>/run_info.txt, sessionInfo.txt
#          <run_dir>/faust_before/, faust_after/        FAUST project folders
#
# Every setting below is a FAUST convention or an earlier user decision (G3, plans/2026-11-10_
# diagnostics-next.md); none is a package default.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
events_dir <- args[[1]]
run_dir <- args[[2]]
labels_dir <- if (length(args) >= 3) args[[3]] else NA_character_
stopifnot("run dir is new" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)

N_EVENTS <- 5000 # events per file (G3); the files hold 5000, so all rows are kept
SUBSAMPLE_SEED <- 1
FAUST_SEED <- 1
THREADS <- 8
DEPTH_SCORE_THRESHOLD <- 0.01 # FAUST convention
SELECTION_QUANTILE <- 0.5 # FAUST convention
NAME_OCCURRENCE_NUM <- 0 # FAUST convention
DENSITY_SUBSAMPLE_THRESHOLD <- 1e6 # FAUST defaults (G3)
DENSITY_SUBSAMPLE_SIZE <- 1e6
DENSITY_SUBSAMPLE_ITERATIONS <- 1
ANNOTATION_FOREST_DEPTH <- 3 # FAUST default, user decision G11
CORRECT_SEED <- 42
CORRECT_NORM_METHOD <- "rank" # user decision Q7; sw_correct_run() default is "scale". Other settings stay at its defaults: 8 x 8 SOM, rlen 10
DROP_MARKERS <- c("AF1-A", "AF-A", "LIVE DEAD NIR-A") # user decision Q8: autofluorescence and viability dye are not phenotype markers; FAUST and the medians use the other 37
MATCH_BY <- "marker_state" # user decision 2026-11-10: FAUST orders markers per run, so exact strings fail to match
N_SHOWN <- 8 # phenotypes drawn in the share figure, largest mean share before correction

ev <- readRDS(file.path(events_dir, "events.rds"))
info <- read.csv(file.path(events_dir, "files.csv"))
mk <- readRDS(file.path(events_dir, "marker_names.rds"))
stopifnot("files match events" = identical(info$file, names(ev)))

# Inputs ---------------------------------------------------------------------------------------
# Anonymous file ids; the mapping stays in files.csv. Channels are renamed to marker names.
ch <- setdiff(colnames(ev[[1]]), "Time")
marker_names <- unname(mk[ch])
is_donor <- grepl("BMV_[0-9]+[A-Z]", info$file)
tbl <- data.frame(
  sample = sprintf("f%03d", seq_len(nrow(info))),
  batch = info$batch,
  donor = ifelse(
    is_donor,
    sub(".*(BMV_[0-9]+[A-Z]).*", "\\1", info$file),
    "reference"
  ),
  role = ifelse(is_donor, "donor", "reference")
)
before <- lapply(ev, \(m) {
  m <- m[, ch]
  colnames(m) <- marker_names
  m
}) |>
  setNames(tbl$sample)
faust_markers <- setdiff(marker_names, DROP_MARKERS)
markers <- setNames(faust_markers, faust_markers)
batch <- setNames(tbl$batch, tbl$sample)
role <- setNames(tbl$role, tbl$sample)

# After data: cyCombine through sw_correct_run() --------------------------------------------------
uncorrected <- bind_rows(lapply(tbl$sample, \(s) {
  as_tibble(before[[s]]) |>
    mutate(sample = s, batch = .env$batch[[s]], row_id = row_number())
}))
t_correct <- system.time(
  corrected <- sw_correct_run(
    uncorrected,
    markers = faust_markers, # user decision G13: correct on the same 37 markers FAUST uses
    seed = CORRECT_SEED,
    norm_method = CORRECT_NORM_METHOD
  )
)
stopifnot(
  "correction kept row order" = identical(
    corrected$row_id,
    uncorrected$row_id
  ) &&
    identical(corrected$sample, uncorrected$sample)
)
after <- lapply(
  split(
    seq_len(nrow(corrected)),
    factor(corrected$sample, levels = tbl$sample)
  ),
  \(i) {
    as.matrix(corrected[i, faust_markers])
  }
)

# FAUST, before and after, same settings -------------------------------------------------------
faust_one <- function(mats, folder) {
  gs <- sw_faust_gatingset(
    mats,
    markers,
    tbl,
    n_events = N_EVENTS,
    subsample_seed = SUBSAMPLE_SEED
  )
  elapsed <- system.time(
    run <- sw_faust_run(
      gs,
      project_path = file.path(run_dir, folder),
      starting_pop = "root",
      markers = faust_markers,
      depth_score_threshold = DEPTH_SCORE_THRESHOLD,
      selection_quantile = SELECTION_QUANTILE,
      experimental_unit = "",
      imputation_hierarchy = "",
      channel_bounds = "",
      name_occurrence_num = NAME_OCCURRENCE_NUM,
      supervised_list = NA,
      seed = FAUST_SEED,
      threads = THREADS,
      density_subsample_threshold = DENSITY_SUBSAMPLE_THRESHOLD,
      density_subsample_size = DENSITY_SUBSAMPLE_SIZE,
      density_subsample_iterations = DENSITY_SUBSAMPLE_ITERATIONS,
      annotation_forest_depth = ANNOTATION_FOREST_DEPTH,
      annotations_approved = TRUE,
      draw_histograms = FALSE
    )
  )
  list(
    labels = sw_faust_labels(run, mats),
    seconds = unname(elapsed[["elapsed"]])
  )
}
if (is.na(labels_dir)) {
  res_before <- faust_one(before, "faust_before")
  res_after <- faust_one(after, "faust_after")
} else {
  res_before <- list(
    labels = readRDS(file.path(labels_dir, "labels_before.rds")),
    seconds = NA
  )
  res_after <- list(
    labels = readRDS(file.path(labels_dir, "labels_after.rds")),
    seconds = NA
  )
}
saveRDS(res_before$labels, file.path(run_dir, "labels_before.rds"))
saveRDS(res_after$labels, file.path(run_dir, "labels_after.rds"))

# Compare -------------------------------------------------------------------------------------
cmp <- sw_faust_compare(
  res_before$labels,
  res_after$labels,
  before,
  after,
  batch = batch,
  role = role,
  markers = faust_markers,
  match_by = MATCH_BY,
  ari_scope = "file",
  median_scope = "file"
)
write.csv(cmp, file.path(run_dir, "faust_compare_long.csv"), row.names = FALSE)

# Figures -------------------------------------------------------------------------------------
# FAUST token M~i~n~ (two levels): the label lists the high markers as M+, then "rest -"; the
# all-zero bin stays "unclassified". Each phenotype is wrapped on its own.
prettify <- function(x) {
  vapply(
    x,
    \(p) {
      if (grepl("^0(_0)*$", p)) {
        return("unclassified")
      }
      tok <- regmatches(p, gregexpr("[^~]+~[0-9]+~[0-9]+~", p))[[1]]
      parts <- strsplit(tok, "~", fixed = TRUE)
      pos <- vapply(
        Filter(\(z) z[[2]] == "2", parts),
        \(z) paste0(z[[1]], "+"),
        ""
      )
      paste(
        strwrap(paste(c(pos, "rest -"), collapse = " "), 30),
        collapse = "\n"
      )
    },
    ""
  )
}
shown <- cmp |>
  filter(check == "prop_before", !is.na(value)) |>
  summarise(mean_share = mean(value), .by = phenotype) |>
  slice_max(mean_share, n = N_SHOWN, with_ties = FALSE) |>
  pull(phenotype)
prop <- cmp |>
  filter(check %in% c("prop_before", "prop_after"), phenotype %in% shown) |>
  transmute(
    sample = subject,
    phenotype = factor(phenotype, levels = shown),
    version = fct(sub("prop_", "", check), levels = c("before", "after")),
    share = value,
    batch = .env$batch[subject],
    role = .env$role[subject]
  )
lab <- setNames(
  unname(prettify(shown)),
  shown
)
stopifnot("panel labels are unique" = !anyDuplicated(lab))
write.csv(
  mutate(prop, phenotype_label = lab[as.character(phenotype)]),
  file.path(run_dir, "figure_source_data_faust_prop_by_batch.csv"),
  row.names = FALSE
)
med <- summarise(
  prop,
  share = median(share, na.rm = TRUE),
  .by = c(phenotype, version, batch)
)
p_prop <- ggplot(prop, aes(batch, share, colour = version)) +
  geom_point(
    data = filter(prop, role == "donor"),
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.7,
      seed = 1
    ),
    size = 0.3,
    alpha = 0.5
  ) +
  geom_point(
    data = filter(prop, role == "reference"),
    aes(shape = "reference file"),
    position = position_dodge(width = 0.7),
    size = 1.8,
    colour = "black",
    fill = NA
  ) +
  geom_crossbar(
    data = med,
    aes(ymin = share, ymax = share),
    position = position_dodge(width = 0.7),
    width = 0.6,
    linewidth = 0.3
  ) +
  scale_shape_manual(values = c("reference file" = 5), name = NULL) +
  scale_colour_manual(values = c(before = "#4C78A8", after = "#E45756")) +
  facet_wrap(
    vars(phenotype),
    ncol = 4,
    scales = "free_y",
    labeller = as_labeller(lab)
  ) +
  scale_x_discrete(labels = \(x) sub("Batch", "B", x)) +
  theme(strip.text = element_text(size = 4.5)) +
  expand_limits(y = 0) +
  labs(
    x = NULL,
    y = "Share of the file's events",
    colour = "Median per batch (bar)"
  )
ggsave(
  file.path(run_dir, "fig_faust_prop_by_batch.png"),
  p_prop,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 200
)

ari <- cmp |>
  filter(check == "ari") |>
  transmute(
    sample = subject,
    ari = value,
    batch = .env$batch[subject],
    role = .env$role[subject]
  )
write.csv(
  ari,
  file.path(run_dir, "figure_source_data_faust_ari.csv"),
  row.names = FALSE
)
p_ari <- ggplot(ari, aes(batch, ari)) +
  geom_point(
    data = filter(ari, role == "donor"),
    position = position_jitter(width = 0.2, height = 0, seed = 1),
    size = 0.6,
    alpha = 0.6
  ) +
  geom_point(
    data = filter(ari, role == "reference"),
    aes(shape = "reference file"),
    size = 2,
    colour = "#E45756"
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    aes(ymin = after_stat(y), ymax = after_stat(y)),
    width = 0.5,
    linewidth = 0.3
  ) +
  scale_shape_manual(values = c("reference file" = 5), name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(labels = \(x) sub("Batch", "B", x)) +
  labs(
    x = NULL,
    y = "ARI before vs after correction",
    caption = "One point per file; bar = batch median"
  )
ggsave(
  file.path(run_dir, "fig_faust_ari.png"),
  p_ari,
  width = 89,
  height = 70,
  units = "mm",
  dpi = 200
)

# Run info ------------------------------------------------------------------------------------
status <- readLines("/proc/self/status")
writeLines(
  c(
    paste("correct_seconds", round(unname(t_correct[["elapsed"]]))),
    paste("faust_before_seconds", round(res_before$seconds)),
    paste("faust_after_seconds", round(res_after$seconds)),
    paste("labels_from", if (is.na(labels_dir)) "this run" else labels_dir),
    paste("match_by", MATCH_BY),
    paste(
      "phenotypes_shared",
      length(intersect(
        unique(cmp$phenotype[cmp$check == "prop_before" & !is.na(cmp$value)]),
        unique(cmp$phenotype[cmp$check == "prop_after" & !is.na(cmp$value)])
      ))
    ),
    paste("n_files", nrow(tbl)),
    paste("n_events_total", sum(lengths(res_before$labels$kept_rows))),
    paste("phenotypes_before", n_distinct(res_before$labels$cells$phenotype)),
    paste("phenotypes_after", n_distinct(res_after$labels$cells$phenotype)),
    paste(
      "unclassified_share_before",
      round(mean(grepl("^0(_0)*$", res_before$labels$cells$phenotype)), 3)
    ),
    paste(
      "unclassified_share_after",
      round(mean(grepl("^0(_0)*$", res_after$labels$cells$phenotype)), 3)
    ),
    grep("^(VmHWM|VmPeak)", status, value = TRUE)
  ),
  file.path(run_dir, "run_info.txt")
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
