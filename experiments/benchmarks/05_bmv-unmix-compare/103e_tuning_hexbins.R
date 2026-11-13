# Purpose: ggcyto hexbins of the four distinct AS 1.8.8 tuning settings (threshold 0.01) on the smoke donor: CD8 vs TNF in
#          CD3+ T cells (uRBC, SEB; the TNF-CD8 streak) and AF vs TNF in live cells (SEB; the AF link of 101g). Red = the arm's
#          TNF edge (unstained q99.5) as a gate, dashed = its floor (unstained q0.5).
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_common_gs.sbatch 103e_tuning_hexbins.R <tuning root>
# Inputs:  <tuning root>/<arm>/*.fcs, unstained/thresholds.csv [103, 103c]; run_in("pop-gates")/labels.csv.gz [99a]
# Outputs: <out_dir> = run_out("tuning-hexbins")
#          <out_dir>/figures/tuning_tnf_hexbin.png   rows = setting, columns = uRBC CD3+, SEB CD3+, SEB live vs AF
#          <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
suppressPackageStartupMessages({
  library(flowCore); library(ggcyto); library(ggplot2)
  library(dplyr, exclude = "filter"); library(purrr); library(readr); library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters -----------------------------------------------------------------------------
ROOT <- commandArgs(trailingOnly = TRUE)[[1]]
ARMS <- c("a188_rFALSE_p1_t0.01", "a188_rTRUE_p1_t0.01", "a188_rFALSE_p5_t0.01", "a188_rTRUE_p5_t0.01") # threshold has no effect
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
TNF <- "BV510-A"
CD8 <- cfg("channels")$CD8
AF <- "AF-A"
N_PLOT <- 100000 # display only
OUT_DIR <- run_out("tuning-hexbins")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Panels ---------------------------------------------------------------------------------------
smoke <- samples_tbl() |> dplyr::filter(.data$donor == cfg("smoke_donor"))
labels <- read_csv(file.path(run_in("pop-gates"), "labels.csv.gz"), show_col_types = FALSE) |>
  dplyr::filter(sample %in% smoke$sample)
set.seed(1)
PANELS <- tribble(~cond, ~pop, ~y, ~pop_name,
  "uRBC", "CD3pos", CD8, "CD3+ T", "SEB", "CD3pos", CD8, "CD3+ T", "SEB", "live", AF, "live")
panels <- pmap(expand.grid(i = seq_len(nrow(PANELS)), arm = ARMS, stringsAsFactors = FALSE), \(i, arm) {
  pn <- PANELS[i, ]
  s <- smoke$sample[smoke$condition == pn$cond]
  m <- exprs(read.FCS(file.path(ROOT, arm, paste0(smoke$stem[smoke$sample == s], " AutoSpectral.fcs")),
                      transformation = FALSE, truncate_max_range = FALSE))
  lb <- labels[labels$sample == s, ]
  stopifnot("rows pair with labels" = nrow(m) == nrow(lb) && isTRUE(all.equal(unname(m[, "FSC-A"]), lb$FSC_A)))
  z <- asinh(m[lb[[pn$pop]], c(TNF, pn$y), drop = FALSE] / COFACTOR)
  th <- read_csv(file.path(ROOT, arm, "unstained", "thresholds.csv"), show_col_types = FALSE)
  e <- th$edge[th$channel == TNF]
  fl <- th$floor[th$channel == TNF]
  pct <- c(100 * mean(z[, TNF] < fl), 100 * mean(z[, TNF] >= e))
  if (nrow(z) > N_PLOT) z <- z[sample(nrow(z), N_PLOT), , drop = FALSE]
  ff <- flowFrame(z)
  identifier(ff) <- "pooled"
  as.ggplot(ggcyto(ff, aes(x = !!rlang::sym(TNF), y = !!rlang::sym(pn$y))) +
    geom_hex(bins = 60) +
    geom_gate(rectangleGate(set_names(list(c(e, Inf)), TNF)), colour = "red") +
    geom_vline(xintercept = fl, linetype = 2, colour = "red") +
    scale_fill_viridis_c(trans = "log10", guide = "none") +
    ggcyto_par_set(limits = list(x = c(-4, 5), y = c(-3, 5))) +
    labs(x = "TNF (BV510-A)", y = if (pn$y == AF) "AF (AF-A)" else "CD8 (SBUV445-A)",
         title = sprintf("%s | %s | %s", arm, s, pn$pop_name),
         subtitle = sprintf("%.2f%% below floor %.2f; %.2f%% above edge %.2f", pct[[1]], fl, pct[[2]], e)) +
    theme(plot.title = element_text(size = 5.5), plot.subtitle = element_text(size = 5.5)))
})
ggsave(file.path(OUT_DIR, "figures", "tuning_tnf_hexbin.png"), cowplot::plot_grid(plotlist = panels, ncol = nrow(PANELS)),
       width = 200, height = 250, units = "mm", dpi = 150, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
