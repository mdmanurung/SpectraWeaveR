# Purpose: show per-cell CD8T TNF (BV510) against IFNg (BV570) for arms at the ends of the two coherence
#          readouts (user, 2026-11-10): TNF+ among IFNg+ highest (swap_spec102, as188_joint) and lowest
#          (as102_var, s102_c188); the TNF-IFNg odds ratio ranks them the other way. Every donor shown.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/68_tnf_ifng_scatter.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-ablation-gated/gs/<arm>/   GatingSets, BV510/BV570 stored arcsinh
#            [19a_gate_arm.R via 51_score_ablation.R]
#          experiments/runs/2026-11-10_bmv-coherence-or/coherence_or.csv   readouts per donor and arm [66]
# Outputs: experiments/runs/2026-11-10_bmv-tnf-ifng-scatter/fig_tnf_ifng_scatter.png
#          experiments/runs/2026-11-10_bmv-tnf-ifng-scatter/quadrant_labels.csv   panel text, per donor and arm
#          experiments/runs/2026-11-10_bmv-tnf-ifng-scatter/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
})

# Parameters ---------------------------------------------------------------------------------------
GATED <- "experiments/runs/2026-11-10_bmv-ablation-gated/gs"
ARMS <- c("swap_spec102", "as188_joint", "as102_var", "s102_c188") # high to low TNF+ among IFNg+
POP <- "CD8T"
CUT <- 1 # fixed cytokine cut, arcsinh (user, 2026-11-10)
AXIS <- c(-1.5, 4) # display window only; points outside are clamped to the edge
OUT_DIR <- "experiments/runs/2026-11-10_bmv-tnf-ifng-scatter"

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# CD8T events per arm and donor ------------------------------------------------------------------
cells <- ARMS |>
  set_names() |>
  map(\(a) {
    gs <- load_gs(file.path(GATED, a))
    sampleNames(gs) |>
      set_names() |>
      map(\(s) {
        x <- exprs(gh_pop_get_data(gs[[s]], POP))
        tibble(TNF = x[, "BV510-A"], IFNg = x[, "BV570-A"])
      }) |>
      list_rbind(names_to = "donor")
  }) |>
  list_rbind(names_to = "arm") |>
  mutate(
    arm = factor(arm, levels = ARMS),
    TNF = pmin(pmax(TNF, AXIS[1]), AXIS[2]),
    IFNg = pmin(pmax(IFNg, AXIS[1]), AXIS[2])
  )

# Panel text: the readouts from 66 for the same arm and donor ---------------------------------------
labels <- read_csv(
  "experiments/runs/2026-11-10_bmv-coherence-or/coherence_or.csv",
  show_col_types = FALSE
) |>
  filter(arm %in% ARMS) |>
  mutate(
    arm = factor(arm, levels = ARMS),
    text = sprintf(
      "TNF+|IFNg+ %.0f%%\nTNF+|IFNg- %.0f%%\nOR %.1f",
      100 * p_pos,
      100 * p_neg,
      odds_ratio
    )
  ) |>
  select(arm, donor, text)
stopifnot(
  "one label per panel" = nrow(labels) == length(ARMS) * n_distinct(cells$donor)
)
write_csv(labels, file.path(OUT_DIR, "quadrant_labels.csv"))

# Figure ---------------------------------------------------------------------------------------------
p <- ggplot(cells, aes(IFNg, TNF)) +
  geom_point(size = 0.05, alpha = 0.15) +
  geom_hline(yintercept = CUT, colour = "red", linewidth = 0.3) +
  geom_vline(xintercept = CUT, colour = "red", linewidth = 0.3) +
  geom_text(
    data = labels,
    aes(x = AXIS[1] + 0.1, y = AXIS[2] - 0.1, label = text),
    hjust = 0,
    vjust = 1,
    size = 2.3
  ) +
  facet_grid(arm ~ donor) +
  coord_fixed(xlim = AXIS, ylim = AXIS) +
  labs(
    x = "IFNg, BV570 (arcsinh, cofactor 6000)",
    y = "TNF, BV510 (arcsinh, cofactor 6000)",
    title = "CD8T cells per arm (rows) and donor (columns); red lines = cut at arcsinh 1"
  ) +
  theme_bw(base_size = 8)
ggsave(
  file.path(OUT_DIR, "fig_tnf_ifng_scatter.png"),
  p,
  width = 11,
  height = 9.5,
  dpi = 150,
  bg = "white"
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
