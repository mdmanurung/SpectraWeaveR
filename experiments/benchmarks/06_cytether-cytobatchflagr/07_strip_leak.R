#!/usr/bin/env Rscript
# Purpose: anchor plan A8 leak check: in the stripped batches, do reference-based effects (unstripped heb_kv_bap)
#          fall inside the anchor-based 90% intervals of heb_kv_bap_anc?
# Run:     Rscript experiments/benchmarks/06_cytether-cytobatchflagr/07_strip_leak.R <strip_results> <full_results> <out_dir>
# Inputs:  <strip_results>/heb_kv_bap_anc/fold_<c>/{effects.csv.gz,status.csv} (02c_correct_strip.R)
#          <full_results>/heb_kv_bap/fold_<c>/effects.csv.gz (02_correct.R, no batch stripped)
# Outputs (in <out_dir>, which must not exist):
#   leak_pairs.csv     fold, batch, population, marker, anchor-based and reference-based effect, se, df, inside
#   leak_summary.csv   per fold and pooled: n, share inside, median |difference|, pass (pooled share >= 0.85)
#   fig_leak.png, fig_leak.pdf   anchor-based vs reference-based effect, one panel per fold
#   sessionInfo.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(ggplot2)
})

# Parameters (approved, anchor plan A8) ------------------------------------------------
LEVEL <- 0.90 # 90% t interval
PASS_SHARE <- 0.85 # user 2026-11-10

args <- commandArgs(trailingOnly = TRUE)
strip_dir <- args[[1]]
full_dir <- args[[2]]
out_dir <- args[[3]]
if (dir.exists(out_dir)) stop(sprintf("%s exists; use a new out_dir.", out_dir), call. = FALSE)

# Pair the effects ---------------------------------------------------------------------
folds <- c("tcu1", "tcu2", "tcs3", "tcs4")
pairs <- map(folds, function(f) {
  stripped <- strsplit(read_csv(file.path(strip_dir, "heb_kv_bap_anc", paste0("fold_", f), "status.csv"),
                                show_col_types = FALSE)$stripped, ";")[[1]]
  anc <- read_csv(file.path(strip_dir, "heb_kv_bap_anc", paste0("fold_", f), "effects.csv.gz"),
                  show_col_types = FALSE) |>
    filter(batch %in% stripped) |>
    select(batch, population, marker, est_anchor = estimate, se_anchor = se_donor, df_anchor = df, support)
  ref <- read_csv(file.path(full_dir, "heb_kv_bap", paste0("fold_", f), "effects.csv.gz"),
                  show_col_types = FALSE) |>
    filter(batch %in% stripped) |>
    select(batch, population, marker, est_reference = estimate)
  inner_join(anc, ref, by = join_by(batch, population, marker), relationship = "one-to-one") |>
    mutate(fold = f, .before = 1)
}) |>
  list_rbind() |>
  mutate(
    half_width = qt(1 - (1 - LEVEL) / 2, df_anchor) * se_anchor,
    inside = abs(est_anchor - est_reference) <= half_width
  )

summary <- bind_rows(
  summarise(pairs, n = n(), share_inside = mean(inside), median_abs_diff = median(abs(est_anchor - est_reference)),
            .by = fold),
  summarise(pairs, fold = "pooled", n = n(), share_inside = mean(inside),
            median_abs_diff = median(abs(est_anchor - est_reference)))
) |>
  mutate(pass = if_else(fold == "pooled", share_inside >= PASS_SHARE, NA))

# Write ------------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE)
write_csv(pairs, file.path(out_dir, "leak_pairs.csv"))
write_csv(summary, file.path(out_dir, "leak_summary.csv"))

# Figure: anchor-based vs reference-based effect; points outside the interval in a second colour.
fig_frame <- mutate(pairs, interval = if_else(inside, "inside 90% interval", "outside"))
p <- ggplot(fig_frame, aes(est_reference, est_anchor, colour = interval)) +
  geom_abline(linewidth = 0.3, linetype = 2) +
  geom_point(size = 0.6, alpha = 0.6) +
  facet_wrap(~fold, nrow = 1) +
  coord_equal() +
  scale_colour_manual(values = c("inside 90% interval" = "grey40", outside = "#d95f02")) +
  labs(x = "effect, reference donors (unstripped fit, asinh)", y = "effect, internal anchors (asinh)",
       colour = NULL, title = "Stripped batches: anchor-based vs reference-based HEB effects") +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig_leak.png"), p, width = 180, height = 70, units = "mm", dpi = 200)
ggsave(file.path(out_dir, "fig_leak.pdf"), p, width = 180, height = 70, units = "mm")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
