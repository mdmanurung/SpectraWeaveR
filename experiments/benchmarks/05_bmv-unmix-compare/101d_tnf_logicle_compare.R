# Purpose: compare arcsinh (cofactor 6000) with logicle for TNF and IFNg in 24H CD4 T cells, per arm. Logicle is estimated
#          on all events (root node, all 30 samples pooled) two ways: per arm, and once on the reference arm and shared.
#          The cut of each sample is the one of 97e (donor uRBC q99.5), turned back to raw intensity and passed through each
#          transform, so % positive can only differ if a transform is not monotone; the figures show the cloud shape.
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> sbatch experiments/slurm/run_common_gs.sbatch 101d_tnf_logicle_compare.R
# Inputs:  run_in("common-gated")/gs/<arm>, run_in("common-cuts")/cuts.csv   [97e]
# Outputs: <out_dir> = run_out("tnf-logicle")
#          <out_dir>/logicle_params.csv                 arm or shared, channel, w, t, m, a, n events used
#          <out_dir>/tnf_pct_by_transform.csv           arm, sample, transform, % TNF+, % IFNg+ in CD4T
#          <out_dir>/figures/tnf_ifng_transforms_<cond>.png   rows = transform, columns = arm, smoke donor, cond = SEB, uRBC
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
REF_ARM <- "as188_joint" # arm whose logicle parameters are shared in the "logicle_shared" variant
POP <- "CD4T"
CH <- c(TNF = "BV510-A", IFNg = "BV570-A")
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
RAW_MIN <- 100 # a channel whose max abs value exceeds this is raw intensity; arcsinh values stay below it
LOGICLE_M <- 4.5 # decades, flowCore default; the only logicle value fixed by hand, w and t are estimated
N_PLOT <- 100000 # display only
LIM_PROBS <- c(0.001, 0.999) # display axis limits per transform; display only
OUT_DIR <- run_out("tnf-logicle")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Raw intensities ------------------------------------------------------------------------
smp <- samples_tbl()
cuts <- read_csv(file.path(run_in("common-cuts"), "cuts.csv"), show_col_types = FALSE) |>
  dplyr::filter(channel %in% CH) |>
  select(arm, sample, channel, pos_cut)
gss <- map(set_names(ARMS), \(a) load_gs(file.path(run_in("common-gated"), "gs", a)))
to_raw <- \(m) {
  asinhed <- apply(abs(m), 2, max) <= RAW_MIN
  m[, asinhed] <- sinh(m[, asinhed, drop = FALSE]) * COFACTOR
  m
}
get_raw <- \(a, s, node) to_raw(exprs(gh_pop_get_data(gss[[a]][[s]], node))[, CH, drop = FALSE])
all_ev <- map(set_names(ARMS), \(a) map(smp$sample, \(s) get_raw(a, s, "root")) |> reduce(rbind))
cd4 <- map(set_names(ARMS), \(a) map(set_names(smp$sample), \(s) get_raw(a, s, POP)))

# Logicle on all events ------------------------------------------------------------------------
fit <- \(m, who) {
  n_used <- nrow(m)
  ff <- flowFrame(m)
  tl <- estimateLogicle(ff, unname(CH), m = LOGICLE_M)
  fun <- map(set_names(unname(CH)), \(ch) tl@transforms[[ch]]@f) # logicle function per channel
  list(
    fun = fun,
    par = map(set_names(unname(CH)), \(ch) {
      e <- environment(fun[[ch]])
      tibble(who = who, channel = ch, w = as.numeric(get("w", e)), t = as.numeric(get("t", e)), m = as.numeric(get("m", e)), a = as.numeric(get("a", e)), n_events = n_used)
    }) |> list_rbind()
  )
}
lg_arm <- map2(all_ev, ARMS, fit)
lg_ref <- fit(all_ev[[REF_ARM]], "shared")
write_csv(list_rbind(c(map(lg_arm, "par"), list(lg_ref$par))), file.path(OUT_DIR, "logicle_params.csv"))

# Transforms of the events and the cuts -------------------------------------------------------------
TRANSFORMS <- c("arcsinh", "logicle_per_arm", "logicle_shared")
tf <- \(variant, a, ch, x) {
  switch(variant,
    arcsinh = asinh(x / COFACTOR),
    logicle_per_arm = lg_arm[[a]]$fun[[ch]](x),
    logicle_shared = lg_ref$fun[[ch]](x))
}
res <- expand.grid(arm = ARMS, sample = smp$sample, variant = TRANSFORMS, stringsAsFactors = FALSE) |>
  as_tibble() |>
  pmap(\(arm, sample, variant) {
    m <- cd4[[arm]][[sample]]
    ct <- map_dbl(CH, \(ch) {
      v <- cuts$pos_cut[cuts$arm == arm & cuts$sample == sample & cuts$channel == ch]
      stopifnot("one cut" = length(v) == 1)
      tf(variant, arm, ch, sinh(v) * COFACTOR) # cut back to raw intensity, then through this transform
    })
    tibble(arm, sample, variant, pct_tnf = 100 * mean(tf(variant, arm, CH[["TNF"]], m[, CH[["TNF"]]]) >= ct[["TNF"]]),
           pct_ifng = 100 * mean(tf(variant, arm, CH[["IFNg"]], m[, CH[["IFNg"]]]) >= ct[["IFNg"]]))
  }) |>
  list_rbind() |>
  left_join(select(smp, sample, donor, condition), by = join_by(sample), relationship = "many-to-one")
write_csv(res, file.path(OUT_DIR, "tnf_pct_by_transform.csv"))

# Figures ------------------------------------------------------------------------------------
set.seed(1)
lims <- map(set_names(TRANSFORMS), \(v) map(set_names(names(CH)), \(k) {
  x <- unlist(map(ARMS, \(a) tf(v, a, CH[[k]], all_ev[[a]][sample.int(nrow(all_ev[[a]]), 20000), CH[[k]]])))
  quantile(x, LIM_PROBS)
}))
for (cond in c("SEB", "uRBC")) {
  s <- paste0(cfg("smoke_donor"), "_", cond)
  g <- expand.grid(variant = TRANSFORMS, arm = ARMS, stringsAsFactors = FALSE) |> as_tibble()
  panels <- pmap(g, \(variant, arm) {
    m <- cd4[[arm]][[s]]
    if (nrow(m) > N_PLOT) m <- m[sample(nrow(m), N_PLOT), , drop = FALSE]
    z <- cbind(tf(variant, arm, CH[["TNF"]], m[, CH[["TNF"]]]), tf(variant, arm, CH[["IFNg"]], m[, CH[["IFNg"]]]))
    colnames(z) <- unname(CH)
    ff <- flowFrame(z)
    identifier(ff) <- "pooled"
    ff@parameters@data$desc <- names(CH)
    r <- res[res$arm == arm & res$sample == s & res$variant == variant, ]
    ct <- map_dbl(CH, \(ch) {
      tf(variant, arm, ch, sinh(cuts$pos_cut[cuts$arm == arm & cuts$sample == s & cuts$channel == ch]) * COFACTOR)
    })
    as.ggplot(ggcyto(ff, aes(x = !!rlang::sym(CH[["TNF"]]), y = !!rlang::sym(CH[["IFNg"]]))) +
      geom_hex(bins = 60) +
      geom_gate(rectangleGate(set_names(list(c(ct[["TNF"]], Inf)), CH[["TNF"]])), colour = "red") +
      geom_gate(rectangleGate(set_names(list(c(ct[["IFNg"]], Inf)), CH[["IFNg"]])), colour = "red") +
      scale_fill_viridis_c(trans = "log10", guide = "none") +
      ggcyto_par_set(limits = list(x = unname(lims[[variant]]$TNF), y = unname(lims[[variant]]$IFNg))) +
      labs(x = "TNF (BV510-A)", y = "IFNg (BV570-A)", title = sprintf("%s | %s | %s | CD4+ T", variant, arm, s),
           subtitle = sprintf("%.2f%% TNF+, %.2f%% IFNg+", r$pct_tnf, r$pct_ifng)) +
      theme(plot.title = element_text(size = 5), plot.subtitle = element_text(size = 5)))
  })
  ggsave(file.path(OUT_DIR, "figures", paste0("tnf_ifng_transforms_", cond, ".png")),
         cowplot::plot_grid(plotlist = panels, ncol = length(ARMS), byrow = FALSE),
         width = 260, height = 190, units = "mm", dpi = 150, bg = "white")
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
