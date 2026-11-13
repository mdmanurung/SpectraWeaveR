# Purpose: find which channels go with the TNF (BV510-A) negatives. For live cells of every 24H sample, compare cells below
#          the floor with the rest on every other fluorescence channel, and rank channels by Spearman correlation with TNF.
#          A channel that is high in the negative cells, or tightly anti-correlated with TNF, points to the source of the artifact.
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch --mem=64G experiments/slurm/run_common_gs.sbatch 101g_tnf_negative_cause.R
# Inputs:  run_in("common-gated")/gs/<arm>   [97e]
#          <data_dir>/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv   marker of each dye [zip]
# Outputs: <out_dir> = run_out("tnf-negative-cause")
#          <out_dir>/negative_vs_rest.csv     arm, channel, marker, median in negative cells, median in the rest, difference
#          <out_dir>/correlation.csv          arm, channel, marker, Spearman rho with TNF over all sampled cells
#          <out_dir>/figures/tnf_vs_top.png   TNF against the five channels with the largest difference, as188_joint and spectroflo
#          <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
suppressPackageStartupMessages({
  library(flowCore); library(flowWorkspace); library(ggcyto); library(ggplot2)
  library(dplyr, exclude = "filter"); library(purrr); library(readr); library(stringr); library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters -----------------------------------------------------------------------------
ARMS <- c("spectroflo", "as188_joint", "as102_var") # worst, the reference, and a clean arm
TNF <- "BV510-A"
PARENT <- "live"
FLOOR <- -1 # arcsinh; the floor of 101f, looked at, not a cut
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
RAW_MIN <- 100 # a channel whose max abs value exceeds this is raw in the GatingSet, so it is transformed here
N_PER_SAMPLE <- 20000 # live cells sampled per sample; keeps all 30 samples equal
N_TOP <- 5 # channels drawn
N_PLOT <- 100000 # display only
SEED <- 1
OUT_DIR <- run_out("tnf-negative-cause")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Sampled live cells per arm ------------------------------------------------------------------
smp <- samples_tbl()
ctrl <- read_csv(data_path("autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"), show_col_types = FALSE) |>
  dplyr::filter(!fluorophore %in% c("AF", "Negative")) |>
  distinct(fluorophore, marker)
marker_of <- set_names(ctrl$marker, paste0(ctrl$fluorophore, "-A"))
lab <- \(ch) if_else(is.na(marker_of[ch]), ch, sprintf("%s (%s)", marker_of[ch], str_remove(ch, "-A$")))
set.seed(SEED)
cells <- map(set_names(ARMS), \(a) {
  gs <- load_gs(file.path(run_in("common-gated"), "gs", a))
  map(smp$sample, \(s) {
    m <- exprs(gh_pop_get_data(gs[[s]], PARENT))
    m <- m[sample(nrow(m), min(N_PER_SAMPLE, nrow(m))), , drop = FALSE]
    keep <- colnames(m)[str_detect(colnames(m), "-A$") & !str_detect(colnames(m), "^(FSC|SSC)")]
    raw <- keep[apply(abs(m[, keep, drop = FALSE]), 2, max) > RAW_MIN]
    m[, raw] <- asinh(m[, raw, drop = FALSE] / COFACTOR)
    m[, c(keep, "SSC-A")]
  }) |> reduce(rbind)
})

# Negative cells against the rest, and correlation --------------------------------------------------
cmp <- map(set_names(ARMS), \(a) {
  m <- cells[[a]]
  neg <- m[, TNF] < FLOOR
  others <- setdiff(colnames(m), c(TNF, "SSC-A"))
  tibble(arm = a, channel = others, marker = unname(marker_of[others]), n_neg = sum(neg),
         median_neg = map_dbl(others, \(ch) median(m[neg, ch])), median_rest = map_dbl(others, \(ch) median(m[!neg, ch]))) |>
    mutate(difference = median_neg - median_rest)
}) |> list_rbind()
write_csv(cmp, file.path(OUT_DIR, "negative_vs_rest.csv"))
cor_tbl <- map(set_names(ARMS), \(a) {
  m <- cells[[a]]
  others <- setdiff(colnames(m), c(TNF, "SSC-A"))
  tibble(arm = a, channel = others, marker = unname(marker_of[others]),
         rho = map_dbl(others, \(ch) suppressWarnings(cor(m[, TNF], m[, ch], method = "spearman"))))
}) |> list_rbind()
write_csv(cor_tbl, file.path(OUT_DIR, "correlation.csv"))

# Figure: TNF against the channels that differ most, per arm ------------------------------------------
set.seed(SEED)
panels <- map(c("as188_joint", "spectroflo"), \(a) {
  tops <- cmp |> dplyr::filter(arm == a) |> slice_max(abs(difference), n = N_TOP) |> pull(channel)
  map(tops, \(ch) {
    z <- cells[[a]][, c(ch, TNF)]
    z <- z[sample(nrow(z), min(N_PLOT, nrow(z))), , drop = FALSE]
    ff <- flowFrame(z)
    identifier(ff) <- "pooled"
    r <- cmp[cmp$arm == a & cmp$channel == ch, ]
    as.ggplot(ggcyto(ff, aes(x = !!rlang::sym(ch), y = !!rlang::sym(TNF))) +
      geom_hex(bins = 60) +
      geom_hline(yintercept = FLOOR, linetype = 2, colour = "red") +
      scale_fill_viridis_c(trans = "log10", guide = "none") +
      labs(x = lab(ch), y = "TNF (BV510)", title = sprintf("%s | live cells, 30 samples", a),
           subtitle = sprintf("median %.2f in TNF<%g cells vs %.2f in the rest", r$median_neg, FLOOR, r$median_rest)) +
      theme(plot.title = element_text(size = 6), plot.subtitle = element_text(size = 6)))
  })
}) |> list_flatten()
ggsave(file.path(OUT_DIR, "figures", "tnf_vs_top.png"), cowplot::plot_grid(plotlist = panels, ncol = N_TOP, byrow = TRUE),
       width = 290, height = 100, units = "mm", dpi = 150, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
