# Purpose: compare the four 24H unmixing arms on TNF and IFNg in CD4 T cells (same AS 1.8.8 labels, same cuts).
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_common_gs.sbatch 101b_tnf_cd4_arms.R
# Inputs:  run_in("common-gated")/gs/<arm>, run_in("common-cuts")/cuts.csv   labels, cuts [97e]
# Outputs: <out_dir> = run_out("tnf-cd4-arms")
#          <out_dir>/tnf_cd4_arms.csv            arm, sample, donor, condition, n CD4T, cuts, % TNF+, % IFNg+
#          <out_dir>/figures/tnf_cd4_arms_pct.png  % TNF+ per donor, by arm and condition
#          <out_dir>/figures/tnf_ifng_cd4_arms_hexbin.png  smoke donor, TNF vs IFNg; rows = arms, columns = uRBC, iRBC, SEB
#          <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
suppressPackageStartupMessages({
  library(flowCore); library(flowWorkspace); library(ggcyto); library(ggplot2)
  library(dplyr, exclude = "filter"); library(purrr); library(readr); library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters -----------------------------------------------------------------------------
ARMS <- c("spectroflo", "as102_basic", "as102_var", "as188_joint") # the four zip arms
POP <- "CD4T"
TNF <- "BV510-A"
IFNG <- "BV570-A"
N_PLOT <- 100000 # events per panel; display only
OUT_DIR <- run_out("tnf-cd4-arms")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Per-arm percentages --------------------------------------------------------------------
smp <- samples_tbl()
cuts <- read_csv(file.path(run_in("common-cuts"), "cuts.csv"), show_col_types = FALSE) |>
  dplyr::filter(channel %in% c(TNF, IFNG)) |>
  tidyr::pivot_wider(id_cols = c(arm, sample), names_from = channel, values_from = pos_cut)
gss <- map(set_names(ARMS), \(a) load_gs(file.path(run_in("common-gated"), "gs", a)))
ev <- \(a, s) exprs(gh_pop_get_data(gss[[a]][[s]], POP))[, c(TNF, IFNG), drop = FALSE]
res <- expand.grid(arm = ARMS, sample = smp$sample, stringsAsFactors = FALSE) |>
  as_tibble() |>
  pmap(\(arm, sample) {
    m <- ev(arm, sample)
    ct <- cuts[cuts$arm == arm & cuts$sample == sample, ]
    stopifnot("one cut row" = nrow(ct) == 1)
    tibble(arm, sample, n_cd4 = nrow(m), cut_tnf = ct[[TNF]], cut_ifng = ct[[IFNG]],
           pct_tnf = 100 * mean(m[, TNF] >= ct[[TNF]]), pct_ifng = 100 * mean(m[, IFNG] >= ct[[IFNG]]))
  }) |>
  list_rbind() |>
  left_join(select(smp, sample, donor, condition), by = join_by(sample), relationship = "many-to-one")
write_csv(res, file.path(OUT_DIR, "tnf_cd4_arms.csv"))

# Figures ------------------------------------------------------------------------------------
res <- mutate(res, arm = factor(arm, levels = ARMS), condition = factor(condition, levels = c("uRBC", "iRBC", "SEB")))
p_pct <- ggplot(res, aes(arm, pct_tnf)) +
  geom_line(aes(group = donor), colour = "grey70", linewidth = 0.2) +
  geom_point(size = 1) +
  facet_wrap(vars(condition), scales = "free_y") +
  labs(x = NULL, y = "% of CD4 T cells above the TNF cut",
       caption = "One point per donor (10), lines join a donor across arms; same labels and cut in every arm.") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "figures", "tnf_cd4_arms_pct.png"), p_pct, width = 130, height = 70, units = "mm", dpi = 150, bg = "white")

set.seed(1)
SHOW <- paste0(cfg("smoke_donor"), "_", c("uRBC", "iRBC", "SEB")) # one donor, three conditions; display only
grid <- expand.grid(s = SHOW, a = ARMS, stringsAsFactors = FALSE) |> as_tibble()
panels <- pmap(grid, \(s, a) {
  m <- ev(a, s)
  r <- res[res$arm == a & res$sample == s, ]
  ff <- flowFrame(if (nrow(m) > N_PLOT) m[sample(nrow(m), N_PLOT), , drop = FALSE] else m)
  identifier(ff) <- "pooled"
  ff@parameters@data$desc <- c("TNF", "IFNg")
  p <- ggcyto(ff, aes(x = !!rlang::sym(TNF), y = !!rlang::sym(IFNG))) +
    geom_hex(bins = 60) +
    geom_gate(rectangleGate(set_names(list(c(r$cut_tnf, Inf)), TNF)), colour = "red") +
    geom_gate(rectangleGate(set_names(list(c(r$cut_ifng, Inf)), IFNG)), colour = "red") +
    scale_fill_viridis_c(trans = "log10", guide = "none") +
    ggcyto_par_set(limits = list(x = c(-2, 5), y = c(-2, 5))) +
    labs(x = "TNF (BV510-A)", y = "IFNg (BV570-A)", title = sprintf("%s | %s | CD4+ T", a, s),
         subtitle = sprintf("%.1f%% TNF+, %.1f%% IFNg+", r$pct_tnf, r$pct_ifng)) +
    theme(plot.title = element_text(size = 6), plot.subtitle = element_text(size = 6))
  as.ggplot(p)
})
ggsave(file.path(OUT_DIR, "figures", "tnf_ifng_cd4_arms_hexbin.png"), cowplot::plot_grid(plotlist = panels, ncol = 3),
       width = 190, height = 240, units = "mm", dpi = 150, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
