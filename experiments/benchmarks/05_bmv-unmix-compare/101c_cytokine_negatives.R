# Purpose: count cells below -1 (arcsinh, cofactor 6000) in the five 24H cytokine channels, per arm, sample and condition,
#          and show IL-10 and Th2 as ggcyto hexbins after the transform they lacked in the GatingSet.
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_common_gs.sbatch 101c_cytokine_negatives.R
# Inputs:  run_in("common-gated")/gs/<arm>, run_in("common-cuts")/cuts.csv   [97e]
# Outputs: <out_dir> = run_out("cytokine-negatives")
#          <out_dir>/negatives.csv                 arm, sample, donor, condition, channel, marker, n, share below -1, median, pos cut, % above cut
#          <out_dir>/figures/negatives_share.png   share below -1 per donor, channel x arm, coloured by condition
#          <out_dir>/figures/il10_th2_hexbin.png   IL-10 and Th2 vs CD8, smoke donor, rows = arm, columns = condition
#          <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
suppressPackageStartupMessages({
  library(flowCore); library(flowWorkspace); library(ggcyto); library(ggplot2)
  library(dplyr, exclude = "filter"); library(purrr); library(readr); library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters -----------------------------------------------------------------------------
ARMS <- c("spectroflo", "as102_basic", "as102_var", "as188_joint")
PARENT <- "CD3pos" # the parent of the cytokine readouts
FLOOR <- -1 # "negative" threshold on arcsinh values; looked at, not a cut
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
RAW_MIN <- 100 # a channel whose max abs value exceeds this is raw, so it is transformed here
MARKERS <- c(`BV510-A` = "TNF", `BV570-A` = "IFNg", `BV605-A` = "IL-2", `BB700-A` = "IL-10", `eFluor 450-A` = "Th2")
CD8 <- cfg("channels")$CD8
N_PLOT <- 100000 # display only
SHOW_ARMS <- c("as102_var", "as188_joint") # hexbin arms; display only
OUT_DIR <- run_out("cytokine-negatives")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Values per arm and sample -----------------------------------------------------------------
smp <- samples_tbl()
cuts <- read_csv(file.path(run_in("common-cuts"), "cuts.csv"), show_col_types = FALSE) |>
  dplyr::filter(channel %in% names(MARKERS)) |>
  select(arm, sample, channel, pos_cut)
gss <- map(set_names(ARMS), \(a) load_gs(file.path(run_in("common-gated"), "gs", a)))
vals <- \(a, s) {
  m <- exprs(gh_pop_get_data(gss[[a]][[s]], PARENT))[, c(names(MARKERS), CD8), drop = FALSE]
  # the GatingSet holds arcsinh values for some channels and raw ones for the rest; arcsinh values stay below RAW_MIN
  raw <- colnames(m)[apply(abs(m), 2, max) > RAW_MIN]
  m[, raw] <- asinh(m[, raw, drop = FALSE] / COFACTOR)
  m
}
grid <- expand.grid(arm = ARMS, sample = smp$sample, stringsAsFactors = FALSE) |> as_tibble()
res <- pmap(grid, \(arm, sample) {
  m <- vals(arm, sample)
  map(names(MARKERS), \(ch) {
    ct <- cuts$pos_cut[cuts$arm == arm & cuts$sample == sample & cuts$channel == ch]
    stopifnot("one cut" = length(ct) == 1)
    tibble(arm, sample, channel = ch, marker = MARKERS[[ch]], n = nrow(m), share_below = mean(m[, ch] < FLOOR),
           median = median(m[, ch]), pos_cut = ct, pct_above = 100 * mean(m[, ch] >= ct))
  }) |> list_rbind()
}) |>
  list_rbind() |>
  left_join(select(smp, sample, donor, condition), by = join_by(sample), relationship = "many-to-one")
write_csv(res, file.path(OUT_DIR, "negatives.csv"))

# Figures ------------------------------------------------------------------------------------
res <- mutate(res, arm = factor(arm, levels = ARMS), marker = factor(marker, levels = MARKERS),
              condition = factor(condition, levels = c("uRBC", "iRBC", "SEB")))
p <- ggplot(res, aes(arm, 100 * share_below, colour = condition)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6, seed = 1), size = 0.8) +
  facet_wrap(vars(marker), scales = "free_y", nrow = 1) +
  labs(x = NULL, y = sprintf("%% of CD3+ T cells below %g (arcsinh)", FLOOR), colour = NULL,
       caption = "One point per donor and condition (10 each).") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(OUT_DIR, "figures", "negatives_share.png"), p, width = 190, height = 80, units = "mm", dpi = 150, bg = "white")

set.seed(1)
smoke <- paste0(cfg("smoke_donor"), "_", c("uRBC", "iRBC", "SEB"))
g2 <- expand.grid(ch = c("BB700-A", "eFluor 450-A"), a = SHOW_ARMS, s = smoke, stringsAsFactors = FALSE) |> as_tibble()
panels <- pmap(g2, \(ch, a, s) {
  m <- vals(a, s)[, c(ch, CD8), drop = FALSE]
  r <- res[res$arm == a & res$sample == s & res$channel == ch, ]
  ff <- flowFrame(if (nrow(m) > N_PLOT) m[sample(nrow(m), N_PLOT), , drop = FALSE] else m)
  identifier(ff) <- "pooled"
  ff@parameters@data$desc <- c(MARKERS[[ch]], "CD8")
  as.ggplot(ggcyto(ff, aes(x = !!rlang::sym(ch), y = !!rlang::sym(CD8))) +
    geom_hex(bins = 60) +
    geom_gate(rectangleGate(set_names(list(c(r$pos_cut, Inf)), ch)), colour = "red") +
    geom_vline(xintercept = FLOOR, colour = "grey40", linetype = 2) +
    scale_fill_viridis_c(trans = "log10", guide = "none") +
    labs(x = sprintf("%s (%s)", MARKERS[[ch]], ch), y = sprintf("CD8 (%s)", CD8),
         title = sprintf("%s | %s | CD3+ T", a, s),
         subtitle = sprintf("%.1f%% below %g; %.1f%% above cut %.2f", 100 * r$share_below, FLOOR, r$pct_above, r$pos_cut)) +
    theme(plot.title = element_text(size = 6), plot.subtitle = element_text(size = 6)))
})
ggsave(file.path(OUT_DIR, "figures", "il10_th2_hexbin.png"), cowplot::plot_grid(plotlist = panels, ncol = 6),
       width = 290, height = 120, units = "mm", dpi = 150, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
