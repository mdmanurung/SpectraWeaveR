# Purpose: old vs new spike-in tables (plans/common_gates_plan.md G5). The old numbers (spiked cells picked from the
#          per-arm as188_joint gates) are the INDEX section 18, 20, 22, 23 and 25 tables; the new ones come from the
#          88cg family (spiked cells picked from the SpectroFlo labels, clean3 pooled gates). Same readout as the
#          INDEX: per table cell, the median over the 5 donors of the share of +1 arcsinh cells losing more than
#          half of the spike (loss cut recovery < 0.5), the donor range, and the median of the donor median
#          recoveries. Also compares the 88f GLUT1 quartiles, the 62/64 af405 proof, the 88h drift counts and the
#          unspiked-cell checks, and how many cells each population holds.
# Run:     sbatch experiments/slurm/run_spike_cg.sbatch 88icg_compare_old_new.R
# Inputs:  old and new score tables under experiments/runs/ (paths in OLD_NEW below) [88 family, 88cg family]
#          old and new spike keys (population sizes and overlap)
# Outputs: <OUT_DIR>/compare_pct_lost.csv     one row per table cell: section, setup or setting, marker, pop, arm,
#                                             old_ and new_ median % lost, min, max, median recovery, delta
#          <OUT_DIR>/compare_af405.csv        62/64 proof: recovery and squeezed share by setup, pop, level
#          <OUT_DIR>/compare_drift.csv        88h unspiked cells moved, old vs new
#          <OUT_DIR>/compare_baseline.csv     unspiked-cell max |spiked - base|, old vs new
#          <OUT_DIR>/compare_population.csv   spiked cells per donor and population, old vs new, and overlap
#          <OUT_DIR>/fig_compare_pct_lost.{pdf,png,csv}  old against new % lost per table cell
#          <OUT_DIR>/sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")

# Parameters ---------------------------------------------------------------------------------------
RUNS <- "experiments/runs"
OUT_DIR <- file.path(RUNS, "2026-11-11_bmv-spike-cg-compare")
FIG_LEVEL <- 1 # the level of the INDEX tables (+1 arcsinh)
SPIKE_FILES <- tribble(
  ~section , ~old                             , ~new                                ,
  "s18"    , "2026-11-10_bmv-spike-in-v2"     , "2026-11-11_bmv-spike-in-cg"        ,
  "s20"    , "2026-11-10_bmv-spike-in-panel2" , "2026-11-11_bmv-spike-in-panel2-cg" ,
  "s22"    , "2026-11-11_bmv-spike-in-cd3neg" , "2026-11-11_bmv-spike-in-cd3neg-cg" ,
  "s25"    , "2026-11-11_bmv-spike-passes"    , "2026-11-11_bmv-spike-passes-cg"
)
KEY_FILES <- tribble(
  ~section , ~old                                   , ~new                                      ,
  "s18"    , "2026-11-10_bmv-spike-in"              , "2026-11-11_bmv-spike-in-cg-build"        ,
  "s20"    , "2026-11-10_bmv-spike-in-panel2-build" , "2026-11-11_bmv-spike-in-panel2-cg-build" ,
  "s22"    , "2026-11-11_bmv-spike-in-cd3neg-build" , "2026-11-11_bmv-spike-in-cd3neg-cg-build"
)
KEY_MARKER <- c(s18 = "TNF", s20 = "IFNg", s22 = "TNF") # one marker per key set: all markers share the cells

stopifnot("compare dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Spike scores: per donor share lost, then median and range over donors ------------------------------------
donor_stats <- \(df, arm_col) {
  df |>
    filter(level == FIG_LEVEL) |>
    rename(arm = all_of(arm_col)) |>
    summarise(
      n_donors = n(),
      median_pct_lost = 100 * median(share_lost),
      min_pct_lost = 100 * min(share_lost),
      max_pct_lost = 100 * max(share_lost),
      median_recovery = median(median_recovery),
      .by = c(arm, marker, pop)
    )
}
spike <- pmap(SPIKE_FILES, \(section, old, new) {
  arm_col <- if (section == "s25") "setting" else "setup"
  read_one <- \(run) {
    read_csv(
      file.path(RUNS, run, "spike_summary.csv"),
      show_col_types = FALSE
    ) |>
      donor_stats(arm_col)
  }
  full_join(
    read_one(old) |> rename_with(\(n) paste0("old_", n), -c(arm, marker, pop)),
    read_one(new) |> rename_with(\(n) paste0("new_", n), -c(arm, marker, pop)),
    by = c("arm", "marker", "pop"),
    relationship = "one-to-one"
  ) |>
    mutate(section = section)
}) |>
  list_rbind()

# Co-spike (section 23): per setup and GLUT1 level ----------------------------------------------------------
cospike_stats <- \(run) {
  read_csv(file.path(run, "cospike_summary.csv"), show_col_types = FALSE) |>
    summarise(
      n_donors = n(),
      median_pct_lost = median(pct_tnf_lost),
      min_pct_lost = min(pct_tnf_lost),
      max_pct_lost = max(pct_tnf_lost),
      median_recovery = median(median_tnf_recovery),
      median_glut1_recovery = median(median_glut1_recovery),
      .by = c(setup, glut1_level, pop)
    )
}
cospike <- full_join(
  cospike_stats(file.path(RUNS, "2026-11-11_bmv-cospike-tnf-glut1")) |>
    rename_with(\(n) paste0("old_", n), -c(setup, glut1_level, pop)),
  cospike_stats(file.path(RUNS, "2026-11-11_bmv-cospike-tnf-glut1-cg")) |>
    rename_with(\(n) paste0("new_", n), -c(setup, glut1_level, pop)),
  by = c("setup", "glut1_level", "pop"),
  relationship = "one-to-one"
) |>
  mutate(
    section = "s23_cospike",
    arm = paste0(setup, " GLUT1 +", glut1_level),
    marker = "TNF"
  ) |>
  select(-setup, -glut1_level)

# TNF loss by natural GLUT1 quartile (section 23, script 88f) -----------------------------------------------
quartile_stats <- \(run) {
  read_csv(
    file.path(RUNS, run, "summary_by_glut1.csv"),
    show_col_types = FALSE
  ) |>
    summarise(
      n_donors = n(),
      median_pct_lost = median(pct_lost),
      min_pct_lost = min(pct_lost),
      max_pct_lost = max(pct_lost),
      median_recovery = median(median_recovery),
      .by = c(setup, glut1_bin, pop)
    )
}
quartile <- full_join(
  quartile_stats("2026-11-11_bmv-tnf-spike-by-glut1") |>
    rename_with(\(n) paste0("old_", n), -c(setup, glut1_bin, pop)),
  quartile_stats("2026-11-11_bmv-tnf-spike-by-glut1-cg") |>
    rename_with(\(n) paste0("new_", n), -c(setup, glut1_bin, pop)),
  by = c("setup", "glut1_bin", "pop"),
  relationship = "one-to-one"
) |>
  mutate(
    section = "s23_glut1_quartile",
    arm = paste(setup, glut1_bin),
    marker = "TNF"
  ) |>
  select(-setup, -glut1_bin)

compare_pct <- bind_rows(spike, cospike, quartile) |>
  mutate(
    delta_pct_lost = new_median_pct_lost - old_median_pct_lost,
    delta_recovery = new_median_recovery - old_median_recovery
  ) |>
  relocate(section, arm, marker, pop)
write_csv(compare_pct, file.path(OUT_DIR, "compare_pct_lost.csv"))

# 62/64 af405 proof ------------------------------------------------------------------------------------------
proof_stats <- \(run) {
  read_csv(file.path(RUNS, run, "spike_summary.csv"), show_col_types = FALSE) |>
    summarise(
      n_donors = n(),
      median_recovery = median(median_recovery),
      median_squeezed_share = median(squeezed_share),
      max_squeezed_share = max(squeezed_share),
      .by = c(setup, pop, level)
    )
}
compare_af405 <- full_join(
  proof_stats("2026-11-10_bmv-af405-proof") |>
    rename_with(\(n) paste0("old_", n), -c(setup, pop, level)),
  proof_stats("2026-11-11_bmv-af405-proof-cg") |>
    rename_with(\(n) paste0("new_", n), -c(setup, pop, level)),
  by = c("setup", "pop", "level"),
  relationship = "one-to-one"
)
write_csv(compare_af405, file.path(OUT_DIR, "compare_af405.csv"))

# 88h drift and baseline checks -----------------------------------------------------------------------------
compare_drift <- full_join(
  read_csv(
    file.path(RUNS, "2026-11-11_bmv-spike-passes/drift_summary.csv"),
    show_col_types = FALSE
  ) |>
    rename_with(\(n) paste0("old_", n), -c(setting, marker, pop)),
  read_csv(
    file.path(RUNS, "2026-11-11_bmv-spike-passes-cg/drift_summary.csv"),
    show_col_types = FALSE
  ) |>
    rename_with(\(n) paste0("new_", n), -c(setting, marker, pop)),
  by = c("setting", "marker", "pop"),
  relationship = "one-to-one"
)
write_csv(compare_drift, file.path(OUT_DIR, "compare_drift.csv"))
baseline <- \(run, arm_col) {
  read_csv(
    file.path(RUNS, run, "baseline_check.csv"),
    show_col_types = FALSE
  ) |>
    rename(arm = all_of(arm_col)) |>
    summarise(
      max_abs_diff = max(max_abs_diff),
      n = sum(n),
      .by = any_of(c("arm", "marker", "pop"))
    )
}
compare_baseline <- pmap(SPIKE_FILES, \(section, old, new) {
  arm_col <- if (section == "s25") "setting" else "setup"
  full_join(
    baseline(old, arm_col) |>
      rename(old_max_abs_diff = max_abs_diff, old_n = n),
    baseline(new, arm_col) |>
      rename(new_max_abs_diff = max_abs_diff, new_n = n),
    by = intersect(c("arm", "marker", "pop"), names(baseline(old, arm_col))),
    relationship = "one-to-one"
  ) |>
    mutate(section = section)
}) |>
  list_rbind()
write_csv(compare_baseline, file.path(OUT_DIR, "compare_baseline.csv"))

# Population sizes and overlap of old and new spiked cells (key sets) ----------------------------------------
compare_population <- pmap(KEY_FILES, \(section, old, new) {
  read_key <- \(run) {
    read_csv(
      file.path(RUNS, run, paste0("spike_key_", KEY_MARKER[[section]], ".csv")),
      show_col_types = FALSE
    ) |>
      select(sample, pop, raw_row, level)
  }
  old_key <- read_key(old)
  new_key <- read_key(new)
  full_join(
    old_key |> summarise(old_n = n(), .by = c(sample, pop)),
    new_key |> summarise(new_n = n(), .by = c(sample, pop)),
    by = c("sample", "pop"),
    relationship = "one-to-one"
  ) |>
    left_join(
      inner_join(
        old_key,
        new_key,
        by = c("sample", "pop", "raw_row"),
        relationship = "one-to-one"
      ) |>
        summarise(n_in_both = n(), .by = c(sample, pop)),
      by = c("sample", "pop"),
      relationship = "one-to-one"
    ) |>
    mutate(section = section)
}) |>
  list_rbind()
write_csv(compare_population, file.path(OUT_DIR, "compare_population.csv"))

# Figure: old against new % lost per table cell --------------------------------------------------------------
plot_frame <- compare_pct |>
  filter(!is.na(old_median_pct_lost), !is.na(new_median_pct_lost)) |>
  mutate(
    section = factor(
      section,
      levels = c(
        "s18",
        "s20",
        "s22",
        "s23_cospike",
        "s23_glut1_quartile",
        "s25"
      )
    )
  )
section_pal <- categorical_palette(levels(plot_frame$section))
p <- ggplot(plot_frame, aes(old_median_pct_lost, new_median_pct_lost)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "grey60",
    linetype = "dashed"
  ) +
  geom_point(aes(colour = section), size = 1) +
  scale_colour_manual(values = section_pal, name = "INDEX section") +
  coord_equal() +
  labs(
    x = "Old: cells losing >50% of spike (%, donor median)",
    y = "New: cells losing >50% of spike (%, donor median)",
    title = "Spike loss, per-arm gates vs SpectroFlo labels"
  )
save_panel(p, plot_frame, OUT_DIR, "fig_compare_pct_lost", 89, 89)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
