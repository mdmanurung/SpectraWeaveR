# Purpose: show the positive cut of each focus marker on real events, as ggcyto bivariate hexbins, before any
#          readout is computed from it. 24H: the five cytokines in CD3+ T cells, one donor, uRBC / iRBC / SEB, with
#          that donor's uRBC-derived cut. Ex vivo: CD1c and HLA-DR in CD3-negative live cells, three donors, with the
#          cut of the config. Reference arm only (AS 1.8.8); other arms use the same cut.
# Run:     SW_DATASET=<24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/100a_cut_figures.R
#          from the repo root; on SLURM: sbatch experiments/slurm/run_common_gs.sbatch 100a_cut_figures.R
# Inputs:  run_in("common-gated")/gs/<ref arm>   GatingSet with the approved labels and cut gates [97e]
#          run_in("common-cuts")/cuts.csv        arm, sample, channel, pos_cut, neg_cut* [97e]
#          <data_dir>/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv   marker of each dye [zip]
#          datasets/<SW_DATASET>.json            channels, smoke_donor, gate_ref_arm
# Outputs: <out_dir> = run_out("cut-figures")
#          <out_dir>/figures/cuts_<ds>.png   rows = markers, columns = samples; red line = the cut
#          <out_dir>/cut_summary.csv         sample, marker, channel, parent, cut, % of parent above the cut
#          <out_dir>/sessionInfo.txt

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters ---------------------------------------------------------------------------------------
REF_ARM <- cfg("gate_ref_arm")
PARENT <- list(`24h` = "CD3pos", exv = "CD3neg_live")[[cfg("dataset")]] # populations the readouts are scored in
PARENT_NAME <- c(CD3pos = "CD3+ T cells", CD3neg_live = "CD3- live cells")[[
  PARENT
]]
Y_CHANNEL <- list(`24h` = cfg("channels")$CD8, exv = "SSC-A")[[cfg("dataset")]] # second axis: CD8 splits T cells; SSC-A for CD3-
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
RAW_MIN <- 100 # a channel whose max abs value exceeds this is raw in the GatingSet, so it is transformed here
N_PLOT <- 100000 # events per panel; display only
OUT_DIR <- run_out("cut-figures")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Samples, markers and cuts ---------------------------------------------------------------------------
smp <- samples_tbl()
cuts <- read_csv(
  file.path(run_in("common-cuts"), "cuts.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(arm == REF_ARM)
channels <- unique(cuts$channel)
ctrl <- read_csv(
  data_path("autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(!fluorophore %in% c("AF", "Negative")) |>
  distinct(fluorophore, marker)
marker_of <- set_names(ctrl$marker, paste0(ctrl$fluorophore, "-A"))
stopifnot(
  "every cut channel has a marker in the control file" = all(
    channels %in% names(marker_of)
  )
)
shown <- if (identical(cfg("dataset"), "24h")) {
  smp |>
    dplyr::filter(donor == cfg("smoke_donor")) |>
    arrange(match(condition, c("uRBC", "iRBC", "SEB")))
} else {
  smp |>
    arrange(sample) |>
    slice(c(1, ceiling(n() / 2), n()))
}
cond_txt <- \(s) {
  if ("condition" %in% names(smp)) smp$condition[smp$sample == s] else "ex vivo"
}

# Events and gates ----------------------------------------------------------------------------------------
gs <- load_gs(file.path(run_in("common-gated"), "gs", REF_ARM))
set.seed(1)
panel <- \(s, ch) {
  m <- exprs(gh_pop_get_data(gs[[s]], PARENT))[, c(ch, Y_CHANNEL), drop = FALSE]
  raw <- colnames(m)[apply(abs(m), 2, max) > RAW_MIN & colnames(m) != "SSC-A"] # cuts are arcsinh; SSC-A stays raw
  m[, raw] <- asinh(m[, raw, drop = FALSE] / COFACTOR)
  cut <- cuts$pos_cut[cuts$sample == s & cuts$channel == ch]
  stopifnot("one cut per sample and channel" = length(cut) == 1)
  pct <- 100 * mean(m[, ch] >= cut)
  ff <- flowFrame(
    if (nrow(m) > N_PLOT) m[sample(nrow(m), N_PLOT), , drop = FALSE] else m
  )
  identifier(ff) <- "pooled"
  p <- ggcyto(ff, aes(x = !!rlang::sym(ch), y = !!rlang::sym(Y_CHANNEL))) +
    geom_hex(bins = 60) +
    geom_gate(rectangleGate(set_names(list(c(cut, Inf)), ch)), colour = "red") +
    scale_fill_viridis_c(trans = "log10", guide = "none") +
    labs(
      x = sprintf("%s (%s)", marker_of[[ch]], ch),
      y = if (Y_CHANNEL == "SSC-A") "SSC-A" else sprintf("CD8 (%s)", Y_CHANNEL),
      title = sprintf("%s | %s | parent: %s", s, cond_txt(s), PARENT_NAME),
      subtitle = sprintf("cut %.2f: %.1f%% of parent above the cut", cut, pct)
    ) +
    theme(
      plot.title = element_text(size = 6),
      plot.subtitle = element_text(size = 6)
    )
  list(
    plot = as.ggplot(p),
    row = tibble(
      sample = s,
      marker = marker_of[[ch]],
      channel = ch,
      parent = PARENT,
      cut = cut,
      pct_above = pct
    )
  )
}
grid <- expand.grid(
  ch = channels,
  s = shown$sample,
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  arrange(match(ch, channels), match(s, shown$sample))
out <- pmap(grid, \(ch, s) panel(s, ch))
write_csv(list_rbind(map(out, "row")), file.path(OUT_DIR, "cut_summary.csv"))
ggsave(
  file.path(OUT_DIR, "figures", paste0("cuts_", cfg("dataset"), ".png")),
  cowplot::plot_grid(plotlist = map(out, "plot"), ncol = nrow(shown)),
  width = 60 * nrow(shown) + 10,
  height = 62 * length(channels) + 10,
  units = "mm",
  dpi = 150,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
