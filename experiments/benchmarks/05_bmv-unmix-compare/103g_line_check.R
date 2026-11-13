# Purpose: line check across 7 arms on the 24H smoke donor (user 2026-11-14): for every marker pair (i, j), in live cells
#          positive for i and negative for j, does j lean with i? Adequate unmixing keeps the i+ population flat along the
#          i axis. Measures: robust slope of j on i (Theil-Sen over the medians of 10 quantile bins of i, arcsinh units per
#          unit) and Spearman rho; report only. Also ggcyto hexbins of the worst pairs and of TNF vs CD8 in CD3+ T cells.
#          Edges (data-driven, user 2026-11-14): AS 1.8.8 tuning arms and AS 1.0.2 arms = q0.5 / q99.5 of their own unmixed
#          unstained control; SpectroFlo = per file and channel, q0.5 / q99.5 of the negative component of a 2-component 1D
#          Gaussian mixture with unequal variances (mclust model "V" fitted by EM here; mclust is not installed).
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch --mem=96G experiments/slurm/run_common_gs.sbatch 103g_line_check.R <tuning root> <unstained-102 dir>
# Inputs:  <tuning root>/<arm>/*.fcs, unstained/thresholds.csv [103, 103c]; <unstained-102 dir>/<arm>/thresholds.csv [103h]
#          run_in("common-gated")/gs/<arm> for spectroflo, as102_basic, as102_var [97e]; run_in("pop-gates")/labels.csv.gz [99a]
# Outputs: <out_dir> = run_out("line-check")
#          <out_dir>/edges.csv          arm, sample, channel, floor, edge, source
#          <out_dir>/line_check.csv     arm, sample, channel_i, channel_j, markers, n, slope, rho
#          <out_dir>/line_summary.csv   arm, sample: median and q90 of |slope| and |rho| over pairs
#          <out_dir>/figures/line_summary.png      |slope| per pair, per arm (summary, not events)
#          <out_dir>/figures/line_worst_pairs.png  ggcyto hexbins of the 6 worst pairs, SEB, every arm, edges as quadrants
#          <out_dir>/figures/tnf_cd8_arms.png      ggcyto TNF vs CD8 in CD3+ T cells, uRBC and SEB, every arm
#          <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
suppressPackageStartupMessages({
  library(flowCore); library(flowWorkspace); library(ggcyto); library(ggplot2)
  library(dplyr, exclude = "filter"); library(purrr); library(readr); library(stringr); library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters -----------------------------------------------------------------------------
cli <- commandArgs(trailingOnly = TRUE)
TUNE_ROOT <- cli[[1]]
U102 <- cli[[2]]
TUNE_ARMS <- c("a188_rFALSE_p1_t0.01", "a188_rTRUE_p1_t0.01", "a188_rFALSE_p5_t0.01", "a188_rTRUE_p5_t0.01")
GS_ARMS <- c("spectroflo", "as102_basic", "as102_var")
ARMS <- c(GS_ARMS, TUNE_ARMS)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
RAW_MIN <- 100 # a channel whose max abs value exceeds this is raw, so it is transformed here
PROBS <- c(0.005, 0.995) # floor and edge (user 2026-11-14)
N_BINS <- 10 # quantile bins of i for the slope
MIN_CELLS <- 200 # i+ j- cells needed for a pair; fewer gives NA
N_SUB <- 20000 # cells per pair and per GMM fit; seeded
SEED <- 1
TNF <- "BV510-A"
CD8 <- cfg("channels")$CD8
N_WORST <- 6
N_PLOT <- 60000 # display only
OUT_DIR <- run_out("line-check")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Events of the smoke donor, paired with the labels --------------------------------------------------
smoke <- samples_tbl() |> dplyr::filter(.data$donor == cfg("smoke_donor"))
labels <- read_csv(file.path(run_in("pop-gates"), "labels.csv.gz"), show_col_types = FALSE) |>
  dplyr::filter(sample %in% smoke$sample)
ctrl <- read_csv(data_path("autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"), show_col_types = FALSE) |>
  dplyr::filter(!fluorophore %in% c("AF", "Negative")) |>
  distinct(fluorophore, marker)
marker_of <- set_names(ctrl$marker, paste0(ctrl$fluorophore, "-A"))
lab <- \(ch) if_else(is.na(marker_of[ch]), ch, sprintf("%s (%s)", marker_of[ch], str_remove(ch, "-A$")))
gss <- map(set_names(GS_ARMS), \(a) load_gs(file.path(run_in("common-gated"), "gs", a)))
events <- \(arm, s) {
  m <- if (arm %in% GS_ARMS) {
    exprs(gh_pop_get_data(gss[[arm]][[s]], "root"))
  } else {
    exprs(read.FCS(file.path(TUNE_ROOT, arm, paste0(smoke$stem[smoke$sample == s], " AutoSpectral.fcs")),
                   transformation = FALSE, truncate_max_range = FALSE))
  }
  lb <- labels[labels$sample == s, ]
  stopifnot("rows pair with labels" = nrow(m) == nrow(lb) && isTRUE(all.equal(unname(m[, "FSC-A"]), lb$FSC_A, tolerance = 1e-6)))
  fl <- colnames(m)[str_detect(colnames(m), "-A$") & !str_detect(colnames(m), "^(FSC|SSC|AF)")]
  raw <- fl[apply(abs(m[, fl, drop = FALSE]), 2, max) > RAW_MIN]
  m[, raw] <- asinh(m[, raw, drop = FALSE] / COFACTOR)
  list(m = m[, c(fl, "SSC-A")], lb = lb, fl = fl)
}
ev <- map(set_names(ARMS), \(a) map(set_names(smoke$sample), \(s) events(a, s)))

# Edges ------------------------------------------------------------------------------------------
gmm_negative <- \(x) {
  # 2-component 1D Gaussian mixture, unequal variances, EM from the 25th and 75th percentiles
  mu <- unname(quantile(x, c(0.25, 0.75))); s2 <- rep(var(x), 2); p <- c(0.5, 0.5)
  for (k in 1:500) {
    d <- cbind(p[1] * dnorm(x, mu[1], sqrt(s2[1])), p[2] * dnorm(x, mu[2], sqrt(s2[2])))
    w <- d / pmax(rowSums(d), .Machine$double.xmin)
    nk <- colSums(w); mu_new <- colSums(w * x) / nk
    s2 <- pmax(colSums(w * (x - rep(mu_new, each = length(x)))^2) / nk, 1e-6); p <- nk / length(x)
    if (max(abs(mu_new - mu)) < 1e-6) { mu <- mu_new; break }
    mu <- mu_new
  }
  i <- which.min(mu)
  qnorm(PROBS, mu[i], sqrt(s2[i]))
}
set.seed(SEED)
edges <- map(ARMS, \(a) {
  if (a %in% TUNE_ARMS || a %in% c("as102_basic", "as102_var")) {
    f <- if (a %in% TUNE_ARMS) file.path(TUNE_ROOT, a, "unstained", "thresholds.csv") else file.path(U102, a, "thresholds.csv")
    th <- read_csv(f, show_col_types = FALSE)
    map(smoke$sample, \(s) mutate(select(th, channel, floor, edge), arm = a, sample = s, source = "unstained control")) |> list_rbind()
  } else {
    map(smoke$sample, \(s) {
      e <- ev[[a]][[s]]
      live <- e$m[e$lb$live, e$fl, drop = FALSE]
      map(e$fl, \(ch) {
        x <- live[, ch]
        q <- gmm_negative(x[sample(length(x), min(N_SUB, length(x)))])
        tibble(channel = ch, floor = q[[1]], edge = q[[2]], arm = a, sample = s, source = "GMM negative")
      }) |> list_rbind()
    }) |> list_rbind()
  }
}) |> list_rbind()
write_csv(edges, file.path(OUT_DIR, "edges.csv"))

# Line check -------------------------------------------------------------------------------------------
theil_sen <- \(x, y) {
  ij <- combn(length(x), 2)
  dx <- x[ij[2, ]] - x[ij[1, ]]
  median((y[ij[2, ]] - y[ij[1, ]])[dx != 0] / dx[dx != 0])
}
set.seed(SEED)
lines <- map(ARMS, \(a) map(smoke$sample, \(s) {
  e <- ev[[a]][[s]]
  live <- e$m[e$lb$live, e$fl, drop = FALSE]
  ed <- edges |> dplyr::filter(arm == a, sample == s)
  edge <- set_names(ed$edge, ed$channel)
  chans <- intersect(e$fl, names(edge))
  pos <- map(set_names(chans), \(ch) live[, ch] >= edge[[ch]])
  map(chans, \(i) map(setdiff(chans, i), \(j) {
    idx <- which(pos[[i]] & !pos[[j]])
    if (length(idx) < MIN_CELLS) return(tibble(arm = a, sample = s, channel_i = i, channel_j = j, n = length(idx), slope = NA_real_, rho = NA_real_))
    if (length(idx) > N_SUB) idx <- sample(idx, N_SUB)
    xi <- live[idx, i]; xj <- live[idx, j]
    b <- cut(xi, unique(quantile(xi, seq(0, 1, length.out = N_BINS + 1))), include.lowest = TRUE)
    mx <- tapply(xi, b, median); my <- tapply(xj, b, median)
    tibble(arm = a, sample = s, channel_i = i, channel_j = j, n = length(idx),
           slope = if (length(mx) >= 3) theil_sen(as.numeric(mx), as.numeric(my)) else NA_real_,
           rho = suppressWarnings(cor(xi, xj, method = "spearman")))
  }) |> list_rbind()) |> list_rbind()
}) |> list_rbind()) |> list_rbind() |>
  mutate(marker_i = unname(marker_of[channel_i]), marker_j = unname(marker_of[channel_j]))
write_csv(lines, file.path(OUT_DIR, "line_check.csv"))
summ <- lines |>
  summarise(n_pairs = sum(!is.na(slope)), median_abs_slope = median(abs(slope), na.rm = TRUE), q90_abs_slope = quantile(abs(slope), 0.9, na.rm = TRUE),
            median_abs_rho = median(abs(rho), na.rm = TRUE), .by = c(arm, sample))
write_csv(summ, file.path(OUT_DIR, "line_summary.csv"))

# Figures ---------------------------------------------------------------------------------------------
lines_f <- lines |> mutate(arm = factor(arm, levels = ARMS))
p1 <- ggplot(dplyr::filter(lines_f, !is.na(slope)), aes(abs(slope), arm)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(data = mutate(summ, arm = factor(arm, levels = ARMS)), aes(median_abs_slope, arm, shape = sample), colour = "red", size = 1.2) +
  coord_cartesian(xlim = c(0, quantile(abs(lines$slope), 0.98, na.rm = TRUE))) +
  labs(x = "|slope| of j on i in i+ j- live cells (arcsinh/arcsinh)", y = NULL, shape = "median, sample",
       caption = "Box: all pairs and samples; 0 = population flat along the i axis.")
ggsave(file.path(OUT_DIR, "figures", "line_summary.png"), p1, width = 160, height = 80, units = "mm", dpi = 150, bg = "white")

quad_panel <- \(a, s, pop, xi, yj, title) {
  e <- ev[[a]][[s]]
  z <- e$m[e$lb[[pop]], c(xi, yj), drop = FALSE]
  if (nrow(z) > N_PLOT) z <- z[sample(nrow(z), N_PLOT), , drop = FALSE]
  ed <- edges |> dplyr::filter(arm == a, sample == s)
  ff <- flowFrame(z)
  identifier(ff) <- "pooled"
  l <- lines |> dplyr::filter(arm == a, sample == s, channel_i == xi, channel_j == yj)
  as.ggplot(ggcyto(ff, aes(x = !!rlang::sym(xi), y = !!rlang::sym(yj))) +
    geom_hex(bins = 60) +
    geom_gate(quadGate(.gate = set_names(c(ed$edge[ed$channel == xi], ed$edge[ed$channel == yj]), c(xi, yj))), colour = "red") +
    scale_fill_viridis_c(trans = "log10", guide = "none") +
    ggcyto_par_set(limits = list(x = c(-3, 5), y = c(-3, 5))) +
    labs(x = lab(xi), y = lab(yj), title = title,
         subtitle = if (nrow(l) == 1) sprintf("i+ j- cells: slope %.3f, rho %.2f (n %d)", l$slope, l$rho, l$n) else "") +
    theme(plot.title = element_text(size = 5), plot.subtitle = element_text(size = 5)))
}
set.seed(SEED)
seb <- smoke$sample[smoke$condition == "SEB"]
common <- reduce(map(ARMS, \(a) ev[[a]][[seb]]$fl), intersect) # channels every arm has
worst <- lines |>
  dplyr::filter(sample == seb, channel_i %in% common, channel_j %in% common) |>
  summarise(m = median(abs(slope), na.rm = TRUE), .by = c(channel_i, channel_j)) |>
  slice_max(m, n = N_WORST)
p2 <- pmap(expand.grid(k = seq_len(nrow(worst)), a = ARMS, stringsAsFactors = FALSE), \(k, a) {
  quad_panel(a, seb, "live", worst$channel_i[[k]], worst$channel_j[[k]], sprintf("%s | %s | live", a, seb))
})
ggsave(file.path(OUT_DIR, "figures", "line_worst_pairs.png"), cowplot::plot_grid(plotlist = p2, ncol = nrow(worst)),
       width = 50 * nrow(worst), height = 45 * length(ARMS), units = "mm", dpi = 130, bg = "white", limitsize = FALSE)
p3 <- pmap(expand.grid(s = smoke$sample[match(c("uRBC", "SEB"), smoke$condition)], a = ARMS, stringsAsFactors = FALSE), \(s, a) {
  quad_panel(a, s, "CD3pos", TNF, CD8, sprintf("%s | %s | CD3+ T", a, s))
})
ggsave(file.path(OUT_DIR, "figures", "tnf_cd8_arms.png"), cowplot::plot_grid(plotlist = p3, ncol = 2),
       width = 110, height = 45 * length(ARMS), units = "mm", dpi = 140, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
