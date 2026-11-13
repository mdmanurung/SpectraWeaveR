# Purpose: find the markers with the most cells at negative values. Every unmixed fluorescence channel, every arm, every
#          24H sample, live cells (arcsinh, cofactor 6000). Reports the share of cells below three floors; nothing is cut.
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch --mem=64G experiments/slurm/run_common_gs.sbatch 101f_negatives_all_markers.R
# Inputs:  run_in("common-gated")/gs/<arm>   GatingSets with the approved labels [97e]
#          <data_dir>/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv   marker of each dye [zip]
# Outputs: <out_dir> = run_out("negatives-all-markers")
#          <out_dir>/negatives_all.csv                  arm, sample, condition, channel, marker, n, share below each floor, q01, median
#          <out_dir>/negatives_ranked.csv               per arm and channel: median and max over samples of share below RANK_FLOOR
#          <out_dir>/figures/negatives_heatmap.png      channel x arm, median % below RANK_FLOOR
#          <out_dir>/figures/negatives_top_hexbin.png   top markers of the reference arm, marker vs SSC-A, smoke donor, SEB
#          <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
suppressPackageStartupMessages({
  library(flowCore); library(flowWorkspace); library(ggcyto); library(ggplot2)
  library(dplyr, exclude = "filter"); library(purrr); library(readr); library(stringr); library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters -----------------------------------------------------------------------------
ARMS <- c("spectroflo", "as102_basic", "as102_var", "as188_joint")
PARENT <- "live" # all live singlets, so every marker is judged on every cell type
FLOORS <- c(-0.5, -1, -2) # arcsinh values below which a cell counts as negative; looked at, never applied as a cut
RANK_FLOOR <- -1 # floor used to rank markers
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
RAW_MIN <- 100 # a channel whose max abs value exceeds this is raw in the GatingSet, so it is transformed here
N_TOP <- 6 # markers drawn as hexbins
N_PLOT <- 100000 # display only
OUT_DIR <- run_out("negatives-all-markers")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Share below each floor ----------------------------------------------------------------------
smp <- samples_tbl()
ctrl <- read_csv(data_path("autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"), show_col_types = FALSE) |>
  dplyr::filter(!fluorophore %in% c("AF", "Negative")) |>
  distinct(fluorophore, marker)
marker_of <- set_names(ctrl$marker, paste0(ctrl$fluorophore, "-A"))
gss <- map(set_names(ARMS), \(a) load_gs(file.path(run_in("common-gated"), "gs", a)))
one <- \(a, s) {
  m <- exprs(gh_pop_get_data(gss[[a]][[s]], PARENT))
  keep <- colnames(m)[str_detect(colnames(m), "-A$") & !str_detect(colnames(m), "^(FSC|SSC)")]
  m <- m[, keep, drop = FALSE]
  raw <- apply(abs(m), 2, max) > RAW_MIN
  m[, raw] <- asinh(m[, raw, drop = FALSE] / COFACTOR)
  map(keep, \(ch) {
    x <- m[, ch]
    tibble(arm = a, sample = s, channel = ch, n = length(x), !!!set_names(map(FLOORS, \(f) mean(x < f)), paste0("share_below_", FLOORS)),
           q01 = unname(quantile(x, 0.01)), median = median(x))
  }) |> list_rbind()
}
res <- expand.grid(arm = ARMS, sample = smp$sample, stringsAsFactors = FALSE) |>
  as_tibble() |>
  pmap(\(arm, sample) one(arm, sample)) |>
  list_rbind() |>
  mutate(marker = unname(marker_of[channel])) |>
  left_join(select(smp, sample, condition), by = join_by(sample), relationship = "many-to-one")
write_csv(res, file.path(OUT_DIR, "negatives_all.csv"))
rank_col <- paste0("share_below_", RANK_FLOOR)
ranked <- res |>
  summarise(median_pct = 100 * median(.data[[rank_col]]), max_pct = 100 * max(.data[[rank_col]]), .by = c(arm, channel, marker)) |>
  arrange(arm, desc(median_pct))
write_csv(ranked, file.path(OUT_DIR, "negatives_ranked.csv"))

# Figures ------------------------------------------------------------------------------------
lab <- \(ch) if_else(is.na(marker_of[ch]), ch, sprintf("%s (%s)", marker_of[ch], str_remove(ch, "-A$")))
ord <- ranked |> dplyr::filter(arm == "as188_joint") |> arrange(median_pct) |> pull(channel)
ph <- ranked |>
  mutate(label = factor(lab(channel), levels = lab(ord)), arm = factor(arm, levels = ARMS)) |>
  ggplot(aes(arm, label, fill = median_pct)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", median_pct)), size = 1.8) +
  scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt") +
  labs(x = NULL, y = NULL, fill = sprintf("%% live cells\nbelow %g\n(median of 30)", RANK_FLOOR),
       caption = "Rows ordered by the AS 1.8.8 arm; arcsinh (cofactor 6000).") +
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "figures", "negatives_heatmap.png"), ph, width = 120, height = 190, units = "mm", dpi = 150, bg = "white")

top <- ranked |> dplyr::filter(arm == "as188_joint") |> slice_max(median_pct, n = N_TOP) |> pull(channel)
s <- paste0(cfg("smoke_donor"), "_SEB")
set.seed(1)
m <- exprs(gh_pop_get_data(gss[["as188_joint"]][[s]], PARENT))
panels <- map(top, \(ch) {
  z <- m[, c(ch, "SSC-A")]
  if (max(abs(z[, ch])) > RAW_MIN) z[, ch] <- asinh(z[, ch] / COFACTOR)
  if (nrow(z) > N_PLOT) z <- z[sample(nrow(z), N_PLOT), , drop = FALSE]
  ff <- flowFrame(z)
  identifier(ff) <- "pooled"
  r <- ranked[ranked$arm == "as188_joint" & ranked$channel == ch, ]
  as.ggplot(ggcyto(ff, aes(x = !!rlang::sym(ch), y = `SSC-A`)) +
    geom_hex(bins = 60) +
    geom_vline(xintercept = RANK_FLOOR, linetype = 2, colour = "red") +
    scale_fill_viridis_c(trans = "log10", guide = "none") +
    labs(x = lab(ch), title = sprintf("%s | as188_joint | parent: %s", s, PARENT),
         subtitle = sprintf("median over 30 samples: %.2f%% below %g", r$median_pct, RANK_FLOOR)) +
    theme(plot.title = element_text(size = 6), plot.subtitle = element_text(size = 6)))
})
ggsave(file.path(OUT_DIR, "figures", "negatives_top_hexbin.png"), cowplot::plot_grid(plotlist = panels, ncol = 3),
       width = 190, height = 120, units = "mm", dpi = 150, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
