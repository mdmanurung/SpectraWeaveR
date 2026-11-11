# experiments/integration/glehr2026/01_invert_authors.R
# Purpose: recover the authors' per-channel CD3 offset m and gain M from their own relativised
#          files (Glehr et al. 2026, BAD study, Aurora), with no gating or choice of ours. Each
#          event of the authors' file is matched to the original panel on scatter and time, which
#          relativisation leaves unchanged; a line fit of relativised on original value per channel
#          gives slope (high - low) / (M - m) and intercept low - m * slope, with low 0, high 1e4.
# Run:     Rscript experiments/integration/glehr2026/01_invert_authors.R <data_dir> <run_dir> <instrument>
# Inputs:  <data_dir>/BAD/aurora/<donor>/*_12-panel.fcs   original and relativised panel per donor
# Outputs: <run_dir>/authors_m_M.csv   donor (d1..), channel, n_matched, slope, intercept, r2, m, M
#          <run_dir>/fig_authors_fit.png, figure_source_data_authors_fit.csv

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(ggplot2)
})
theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[[1]]
run_dir <- args[[2]]
instrument <- args[[3]] # "aurora" or "fortessa"
if (dir.exists(run_dir)) stop("run dir exists: ", run_dir, call. = FALSE)
dir.create(run_dir, recursive = TRUE)

# Paper targets (Box 1): CD3- median to 0, CD3+ median to 1e4.
LOW <- 0
HIGH <- 1e4
KEY <- c("TIME", "FSC-A", "FSC-H", "SSC-A", "SSC-H")

read_exprs <- function(path) {
  exprs(read.FCS(path, truncate_max_range = FALSE, transformation = FALSE, emptyValue = FALSE))
}

donors <- sort(list.dirs(file.path(data_dir, "BAD", instrument), recursive = FALSE))
fits <- lapply(seq_along(donors), function(i) {
  fs <- list.files(donors[i], pattern = "12-panel\\.fcs$", full.names = TRUE)
  orig <- read_exprs(fs[grepl("_ungated_", fs)])
  rel <- read_exprs(fs[grepl("relativised", fs)])
  # Match events on the unchanged columns; drop keys that are not unique in the original.
  k_orig <- do.call(paste, as.data.frame(orig[, KEY]))
  k_rel <- do.call(paste, as.data.frame(rel[, KEY]))
  uniq <- !(k_orig %in% k_orig[duplicated(k_orig)])
  hit <- match(k_rel, k_orig[uniq])
  o <- orig[uniq, , drop = FALSE][hit[!is.na(hit)], ]
  r <- rel[!is.na(hit), ]
  chans <- setdiff(colnames(orig), KEY)
  lapply(chans, function(ch) {
    f <- stats::lm(r[, ch] ~ o[, ch])
    s <- unname(coef(f)[2])
    c0 <- unname(coef(f)[1])
    tibble(
      donor = paste0("d", i),
      channel = ch,
      n_rel = nrow(rel),
      n_matched = nrow(r),
      slope = s,
      intercept = c0,
      r2 = summary(f)$r.squared,
      m = (LOW - c0) / s,
      M = (LOW - c0) / s + (HIGH - LOW) / s
    )
  }) |>
    bind_rows()
}) |>
  bind_rows()
write.csv(fits, file.path(run_dir, "authors_m_M.csv"), row.names = FALSE)
print(as.data.frame(mutate(fits, across(where(is.numeric), \(x) signif(x, 4)))))

# Figure: implied m and M per channel, one point per donor.
plot_frame <- fits |>
  tidyr::pivot_longer(c(m, M), names_to = "parameter", values_to = "value")
write.csv(plot_frame, file.path(run_dir, "figure_source_data_authors_fit.csv"), row.names = FALSE)
p <- ggplot(plot_frame, aes(value, channel)) +
  geom_point(size = 0.8, alpha = 0.7) +
  facet_wrap(vars(parameter), scales = "free_x") +
  labs(x = "Implied value on the original linear scale (one point per donor)", y = NULL)
ggsave(file.path(run_dir, "fig_authors_fit.png"), p, width = 160, height = 70, units = "mm", dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
