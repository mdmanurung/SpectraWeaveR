# Purpose: re-run the three mindensity2 cuts of the 24H pooled AS 1.8.8 gates (viability, CD3, CD56 on CD3-) in logicle space
#          and compare them with the arcsinh cuts, on the same events. The arcsinh search window of each gate is mapped
#          through raw intensity into logicle, so the window covers the same raw range. Logicle is estimated on all events
#          of the pooled AS 1.8.8 files, one fit per channel.
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch --mem=96G experiments/slurm/run_common_gs.sbatch 101e_mindensity_logicle.R
# Inputs:  run_in("common-gated")/gs/<gate_ref_arm>   GatingSet with the approved parents [97e]
# Outputs: <out_dir> = run_out("mindensity-logicle")
#          <out_dir>/cuts_by_transform.csv        gate, channel, parent, transform, window, cut (own scale), cut (raw), % of parent above
#          <out_dir>/logicle_params.csv           channel, w, t, m, a
#          <out_dir>/figures/mindensity_<gate>.png   FSC-A vs channel in the parent, arcsinh and logicle, red = cut, dashed = window
#          <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
suppressPackageStartupMessages({
  library(flowCore); library(flowWorkspace); library(openCyto); library(ggcyto); library(ggplot2)
  library(dplyr, exclude = "filter"); library(purrr); library(readr); library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters -----------------------------------------------------------------------------
REF_ARM <- cfg("gate_ref_arm")
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
RAW_MIN <- 100 # a channel whose max abs value exceeds this is raw intensity in the GatingSet
LOGICLE_M <- 4.5 # decades, flowCore default
GATES <- tibble::tribble( # arcsinh search windows are the approved ones (viability and CD3 0 to 3, CD56 0.3 to 2.5)
  ~gate, ~channel, ~parent, ~lo, ~hi,
  "live", "viability", "singlets", 0, 3,
  "CD3pos", cfg("channels")$CD3, "live", 0, 3,
  "NKpos", "RB744-A", "CD3neg_live", 0.3, 2.5
)
N_PLOT <- 100000 # display only
OUT_DIR <- run_out("mindensity-logicle")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Events ---------------------------------------------------------------------------------------
smp <- samples_tbl()
GATES$channel[GATES$gate == "live"] <- cfg("channels")$viability_autospectral
gs <- load_gs(file.path(run_in("common-gated"), "gs", REF_ARM))
chs <- unique(c(GATES$channel, "FSC-A"))
to_raw <- \(m) {
  a <- apply(abs(m), 2, max) <= RAW_MIN
  m[, a] <- sinh(m[, a, drop = FALSE]) * COFACTOR
  m
}
get_raw <- \(node) map(smp$sample, \(s) to_raw(exprs(gh_pop_get_data(gs[[s]], node))[, chs, drop = FALSE])) |> reduce(rbind)
all_ev <- get_raw("root")
parents <- map(set_names(unique(GATES$parent)), get_raw)

# Logicle on all events, one fit per channel ------------------------------------------------------
tl <- estimateLogicle(flowFrame(all_ev), GATES$channel, m = LOGICLE_M)
lfun <- map(set_names(GATES$channel), \(ch) tl@transforms[[ch]]@f)
params <- map(GATES$channel, \(ch) {
  e <- environment(lfun[[ch]])
  tibble(channel = ch, w = as.numeric(get("w", e)), t = as.numeric(get("t", e)), m = as.numeric(get("m", e)), a = as.numeric(get("a", e)))
}) |> list_rbind()
write_csv(params, file.path(OUT_DIR, "logicle_params.csv"))

# mindensity2 in both spaces ------------------------------------------------------------------
fwd <- list(arcsinh = \(ch, x) asinh(x / COFACTOR), logicle = \(ch, x) lfun[[ch]](x))
cut_of <- \(x, ch, win) {
  m <- cbind(x)
  colnames(m) <- ch
  g <- gate_mindensity2(flowFrame(m), channel = ch, gate_range = win)
  c(g@min[[ch]], g@max[[ch]])[is.finite(c(g@min[[ch]], g@max[[ch]]))][[1]]
}
res <- pmap(GATES, \(gate, channel, parent, lo, hi) {
  raw <- parents[[parent]][, channel]
  raw_win <- sinh(c(lo, hi)) * COFACTOR # same raw window in both spaces
  imap(fwd, \(f, tr) {
    win <- f(channel, raw_win)
    cut <- cut_of(f(channel, raw), channel, win)
    inv <- if (tr == "arcsinh") sinh(cut) * COFACTOR else approx(f(channel, sort(raw)), sort(raw), xout = cut, ties = "ordered")$y
    tibble(gate, channel, parent, transform = tr, win_lo = win[[1]], win_hi = win[[2]], cut = cut, cut_raw = inv,
           pct_above = 100 * mean(f(channel, raw) >= cut), n_parent = length(raw))
  }) |> list_rbind()
}) |> list_rbind()
write_csv(res, file.path(OUT_DIR, "cuts_by_transform.csv"))

# Figures: the channel against FSC-A in the parent, both spaces ----------------------------------------
set.seed(1)
walk(seq_len(nrow(GATES)), \(i) {
  g <- GATES[i, ]
  idx <- sample(nrow(parents[[g$parent]]), min(N_PLOT, nrow(parents[[g$parent]])))
  panels <- imap(fwd, \(f, tr) {
    r <- res[res$gate == g$gate & res$transform == tr, ]
    z <- cbind(f(g$channel, parents[[g$parent]][idx, g$channel]), parents[[g$parent]][idx, "FSC-A"])
    colnames(z) <- c(g$channel, "FSC-A")
    ff <- flowFrame(z)
    identifier(ff) <- "pooled"
    as.ggplot(ggcyto(ff, aes(x = !!rlang::sym(g$channel), y = `FSC-A`)) +
      geom_hex(bins = 60) +
      geom_gate(rectangleGate(set_names(list(c(r$cut, Inf)), g$channel)), colour = "red") +
      geom_vline(xintercept = c(r$win_lo, r$win_hi), linetype = 2, colour = "grey40") +
      scale_fill_viridis_c(trans = "log10", guide = "none") +
      labs(x = g$channel, title = sprintf("%s | %s | parent: %s | AS 1.8.8 pooled", g$gate, tr, g$parent),
           subtitle = sprintf("cut %.3f (raw %.0f): %.1f%% of parent above; dashed = search window", r$cut, r$cut_raw, r$pct_above)) +
      theme(plot.title = element_text(size = 6), plot.subtitle = element_text(size = 6)))
  })
  ggsave(file.path(OUT_DIR, "figures", paste0("mindensity_", g$gate, ".png")), cowplot::plot_grid(plotlist = panels, ncol = 2),
         width = 190, height = 80, units = "mm", dpi = 150, bg = "white")
})
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
