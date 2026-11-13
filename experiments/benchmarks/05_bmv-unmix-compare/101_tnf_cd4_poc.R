# Purpose: proof of concept. TNF (BV510-A) staining in CD4 T cells of the 24H dataset after medium (uRBC), iRBC and SEB,
#          next to the 3D dataset, both in AS 1.8.8 and on each dataset's own CD4T label. TNF is on BV510 in both panels.
#          24H cut: the donor's uRBC q99.5 (user decision); 3D cut: fixed arcsinh 1.5 (3D decision G7a). The two cut
#          rules differ, and 3D has no unstimulated control, so the percentages are shown but not equated.
# Run:     SW_DATASET=24h SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/101_tnf_cd4_poc.R
#          from the repo root; on SLURM: sbatch experiments/slurm/run_common_gs.sbatch 101_tnf_cd4_poc.R
# Inputs:  run_in("common-gated")/gs/as188_joint, run_in("common-cuts")/cuts.csv   24H labels and cuts [97e]
#          3D: key common-gated of datasets/3d_runs.csv, gs/as188_joint   3D labels and gates [97e, 2026-11-11]
# Outputs: <out_dir> = run_out("tnf-cd4-poc")
#          <out_dir>/figures/tnf_cd4_hexbin.png    row 1: 24H donor, uRBC / iRBC / SEB; row 2: three 3D donors
#          <out_dir>/figures/tnf_cd4_pct.png       % TNF+ in CD4T per donor: 24H by condition, 3D
#          <out_dir>/tnf_cd4_summary.csv           dataset, sample, donor, condition, n CD4T, cuts, % TNF+, % IFNg+
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
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters ---------------------------------------------------------------------------------------
stopifnot("run with SW_DATASET=24h" = identical(cfg("dataset"), "24h"))
ARM <- "as188_joint" # reference arm, both datasets
POP <- "CD4T" # node alias in both GatingSets
TNF <- "BV510-A" # TNF on BV510 in 3D and 24H (checked against both control files)
IFNG <- "BV570-A" # IFNg on BV570 in 3D and 24H
CUT_3D <- 1.5 # 3D cut for every cytokine (3D decision G7a)
SHOW_3D <- c("A3", "A6", "B11") # 3D donors drawn as hexbins; display only
N_PLOT <- 100000 # events per panel; display only
COND_ORDER <- c("uRBC", "iRBC", "SEB")
OUT_DIR <- run_out("tnf-cd4-poc")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# GatingSets and cuts -------------------------------------------------------------------------------------
smp <- samples_tbl()
gs24 <- load_gs(file.path(run_in("common-gated"), "gs", ARM))
cuts24 <- read_csv(
  file.path(run_in("common-cuts"), "cuts.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(arm == ARM, channel %in% c(TNF, IFNG)) |>
  tidyr::pivot_wider(
    id_cols = sample,
    names_from = channel,
    values_from = pos_cut
  )
path3d <- read_csv(
  "experiments/benchmarks/05_bmv-unmix-compare/datasets/3d_runs.csv",
  show_col_types = FALSE
) |>
  dplyr::filter(key == "common-gated") |>
  pull(path)
stopifnot("one 3D common-gated run" = length(path3d) == 1)
gs3d <- load_gs(file.path(path3d, "gs", ARM))

events <- \(gs, s) {
  exprs(gh_pop_get_data(gs[[s]], POP))[, c(TNF, IFNG), drop = FALSE]
}
summ <- bind_rows(
  map(smp$sample, \(s) {
    m <- events(gs24, s)
    ct <- cuts24[cuts24$sample == s, ]
    tibble(
      dataset = "24H",
      sample = s,
      n_cd4 = nrow(m),
      cut_tnf = ct[[TNF]],
      cut_ifng = ct[[IFNG]],
      pct_tnf = 100 * mean(m[, TNF] >= ct[[TNF]]),
      pct_ifng = 100 * mean(m[, IFNG] >= ct[[IFNG]])
    )
  }) |>
    list_rbind() |>
    left_join(
      select(smp, sample, donor, condition),
      by = join_by(sample),
      relationship = "one-to-one"
    ),
  map(sampleNames(gs3d), \(s) {
    m <- events(gs3d, s)
    tibble(
      dataset = "3D",
      sample = s,
      n_cd4 = nrow(m),
      cut_tnf = CUT_3D,
      cut_ifng = CUT_3D,
      pct_tnf = 100 * mean(m[, TNF] >= CUT_3D),
      pct_ifng = 100 * mean(m[, IFNG] >= CUT_3D),
      donor = s,
      condition = "3D stimulated (no control)"
    )
  }) |>
    list_rbind()
)
write_csv(summ, file.path(OUT_DIR, "tnf_cd4_summary.csv"))

# Hexbins: TNF against IFNg in CD4 T cells, the cuts as one gate -----------------------------------------------------
set.seed(1)
panel <- \(gs, s, title, subtitle, cut_t, cut_i) {
  m <- events(gs, s)
  ff <- flowFrame(
    if (nrow(m) > N_PLOT) m[sample(nrow(m), N_PLOT), , drop = FALSE] else m
  )
  identifier(ff) <- "pooled"
  ff@parameters@data$desc <- c("TNF", "IFNg")
  p <- ggcyto(ff, aes(x = !!rlang::sym(TNF), y = !!rlang::sym(IFNG))) +
    geom_hex(bins = 60) +
    geom_gate(
      rectangleGate(set_names(list(c(cut_t, Inf)), TNF)),
      colour = "red"
    ) +
    geom_gate(
      rectangleGate(set_names(list(c(cut_i, Inf)), IFNG)),
      colour = "red"
    ) +
    scale_fill_viridis_c(trans = "log10", guide = "none") +
    ggcyto_par_set(limits = list(x = c(-2, 5), y = c(-2, 5))) +
    labs(
      x = sprintf("TNF (%s)", TNF),
      y = sprintf("IFNg (%s)", IFNG),
      title = title,
      subtitle = subtitle
    ) +
    theme(
      plot.title = element_text(size = 6),
      plot.subtitle = element_text(size = 6)
    )
  as.ggplot(p)
}
donor24 <- cfg("smoke_donor")
p24 <- map(COND_ORDER, \(cn) {
  s <- smp$sample[smp$donor == donor24 & smp$condition == cn]
  r <- summ[summ$sample == s, ]
  panel(
    gs24,
    s,
    sprintf(
      "24H | donor %s | %s | AS 1.8.8 | parent: CD4+ T cells",
      donor24,
      cn
    ),
    sprintf(
      "TNF cut %.2f (donor uRBC q99.5): %.1f%% TNF+; %d CD4T events",
      r$cut_tnf,
      r$pct_tnf,
      r$n_cd4
    ),
    r$cut_tnf,
    r$cut_ifng
  )
})
p3d <- map(SHOW_3D, \(s) {
  r <- summ[summ$sample == s, ]
  panel(
    gs3d,
    s,
    sprintf(
      "3D | donor %s | stimulated, no unstimulated control | AS 1.8.8 | parent: CD4+ T cells",
      s
    ),
    sprintf(
      "TNF cut %.2f (fixed): %.1f%% TNF+; %d CD4T events",
      CUT_3D,
      r$pct_tnf,
      r$n_cd4
    ),
    CUT_3D,
    CUT_3D
  )
})
ggsave(
  file.path(OUT_DIR, "figures", "tnf_cd4_hexbin.png"),
  cowplot::plot_grid(plotlist = c(p24, p3d), ncol = 3),
  width = 190,
  height = 135,
  units = "mm",
  dpi = 150,
  bg = "white"
)

# % TNF+ per donor (summary of the same counts; not a gate figure) ----------------------------------------------------
pct_plot <- summ |>
  mutate(
    group = if_else(dataset == "24H", condition, "3D"),
    group = factor(group, levels = c(COND_ORDER, "3D"))
  )
med <- pct_plot |>
  summarise(med = median(pct_tnf), .by = group)
p_pct <- ggplot(pct_plot, aes(group, pct_tnf)) +
  geom_point(
    aes(colour = group),
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 1.4
  ) +
  geom_crossbar(
    data = med,
    aes(y = med, ymin = med, ymax = med),
    width = 0.5,
    linewidth = 0.3
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = NULL,
    y = "% of CD4 T cells above the TNF cut",
    colour = NULL,
    caption = "One point per donor (24H: 10; 3D: 5); bar = median. 24H cut = donor's uRBC q99.5; 3D cut = fixed 1.5."
  ) +
  theme(legend.position = "none")
ggsave(
  file.path(OUT_DIR, "figures", "tnf_cd4_pct.png"),
  p_pct,
  width = 89,
  height = 75,
  units = "mm",
  dpi = 150,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
