# Purpose: show per-cell CD8T TNF (BV510) for the per-dye variant swap (45_variant_per_dye.R), against
#          IFNg (BV570) and against GLUT1 (Alexa Fluor 405), so the squeeze is visible (user, 2026-11-10).
#          Arms: 1.8.8 baseline, SBUV510-only and Alexa Fluor 405-only 1.0.2 variants (the two dyes that
#          moved TNF), and all 1.0.2 variants. Each arm on its own gates. Every donor shown.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/71_per_dye_scatter.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-ablation-gated/gs/<arm>/   GatingSets; BV510/BV570 stored arcsinh,
#            Alexa Fluor 405 stored linear  [19a_gate_arm.R via 51_score_ablation.R]
# Outputs: experiments/runs/2026-11-10_bmv-per-dye-scatter/fig_per_dye_scatter.png
#          experiments/runs/2026-11-10_bmv-per-dye-scatter/panel_labels.csv   TNF+ % per arm and donor (panel text)
#          experiments/runs/2026-11-10_bmv-per-dye-scatter/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters ---------------------------------------------------------------------------------------
GATED <- "experiments/runs/2026-11-10_bmv-ablation-gated/gs"
ARMS <- c(
  "1.8.8 variants" = "as188_joint",
  "+ SBUV510 1.0.2" = "var102_only_SBUV510",
  "+ AF405 (GLUT1) 1.0.2" = "var102_only_Alexa_Fluor_405",
  "all 1.0.2 variants" = "swap_var102"
)
POP <- "CD8T"
CUT <- 1 # fixed cytokine cut, arcsinh (user, 2026-11-10)
COFACTOR <- 6000 # arcsinh cofactor for the linear GLUT1 channel, as for the cytokines (CLAUDE.md)
AXIS <- c(-1.5, 4) # display window only; points outside are clamped to the edge
OUT_DIR <- "experiments/runs/2026-11-10_bmv-per-dye-scatter"

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# CD8T events per arm and donor ------------------------------------------------------------------
cells <- ARMS |>
  map(\(a) {
    gs <- load_gs(file.path(GATED, a))
    sampleNames(gs) |>
      set_names() |>
      map(\(s) {
        x <- exprs(gh_pop_get_data(gs[[s]], POP))
        tibble(
          TNF = x[, "BV510-A"],
          IFNg = x[, "BV570-A"],
          GLUT1 = asinh(x[, "Alexa Fluor 405-A"] / COFACTOR)
        )
      }) |>
      list_rbind(names_to = "donor")
  }) |>
  list_rbind(names_to = "arm") |>
  mutate(arm = factor(arm, levels = names(ARMS)))

# Panel text: TNF+ % at the cut, from the plotted cells -------------------------------------------
labels <- cells |>
  summarise(pct_tnf = 100 * mean(TNF > CUT), n = n(), .by = c(arm, donor)) |>
  mutate(text = sprintf("TNF+ %.0f%%", pct_tnf))
write_csv(labels, file.path(OUT_DIR, "panel_labels.csv"))

# Plot frame: one row block per x marker, TNF on y, clamped to the window ----------------------------
plot_df <- cells |>
  pivot_longer(c(IFNg, GLUT1), names_to = "x_marker", values_to = "x") |>
  mutate(
    x_marker = factor(x_marker, levels = c("IFNg", "GLUT1")),
    x = pmin(pmax(x, AXIS[1]), AXIS[2]),
    TNF = pmin(pmax(TNF, AXIS[1]), AXIS[2])
  )
label_df <- labels |>
  cross_join(tibble(
    x_marker = factor(c("IFNg", "GLUT1"), levels = c("IFNg", "GLUT1"))
  ))

# Figure ---------------------------------------------------------------------------------------------
p <- ggplot(plot_df, aes(x, TNF)) +
  geom_point(size = 0.05, alpha = 0.15) +
  geom_hline(yintercept = CUT, colour = "red", linewidth = 0.3) +
  geom_text(
    data = label_df,
    aes(x = AXIS[1] + 0.1, y = AXIS[2] - 0.1, label = text),
    hjust = 0,
    vjust = 1,
    size = 2.3
  ) +
  facet_grid(
    x_marker + arm ~ donor,
    labeller = label_wrap_gen(width = 14, multi_line = FALSE)
  ) +
  coord_fixed(xlim = AXIS, ylim = AXIS) +
  labs(
    x = "IFNg, BV570 (top block) or GLUT1, Alexa Fluor 405 (bottom block); arcsinh, cofactor 6000",
    y = "TNF, BV510 (arcsinh, cofactor 6000)",
    title = "CD8T cells: 1.0.2 variants swapped in one dye at a time (rows) per donor (columns)",
    subtitle = "Red line = TNF cut at arcsinh 1. Each arm on its own gates."
  ) +
  theme_bw(base_size = 8)
ggsave(
  file.path(OUT_DIR, "fig_per_dye_scatter.png"),
  p,
  width = 9,
  height = 15,
  dpi = 150,
  bg = "white"
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
