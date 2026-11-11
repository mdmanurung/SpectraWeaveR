#!/usr/bin/env Rscript
# Purpose: check whether each arm's TNF signal hangs together with the cell's other activation readouts,
#          with no unstimulated control. Two readouts per arm, donor and parent (CD4T, CD8T):
#          (1) % TNF+ among IFNg+ cells, at the approved fixed cut (arcsinh 1 for both cytokines);
#          (2) Spearman rho of TNF, and of IFNg as the reference cytokine, with each activation or effector
#          marker in the panel. Rank correlation needs no further cut.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/46_tnf_coherence.R \
#            experiments/runs/2026-11-10_bmv-gated-pooled experiments/runs/2026-11-10_bmv-tnf-coherence
#          from the repo root.
# Inputs:  <gated_dir>/manifest.csv [01_manifest-gate.R]
#          <gated_dir>/gs/<arm>/ saved GatingSets with CD4T and CD8T [19a_gate_arm.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  ARM_LABELS, VERSION_ROWS, ZIP_ARMS
# Outputs: <run_dir>/tnf_among_ifng.csv           % TNF+ of IFNg+ cells per arm, donor, parent
#          <run_dir>/rho_activation.csv            Spearman rho per arm, donor, parent, cytokine, marker
#          <run_dir>/fig_a_tnf_among_ifng.png      panel a
#          <run_dir>/fig_a_source.csv              data plotted in panel a
#          <run_dir>/fig_b_rho_activation.png      panel b
#          <run_dir>/fig_b_source.csv              data plotted in panel b
#          <run_dir>/sessionInfo.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # same ggplot2 as 38/40

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS, VERSION_ROWS, ZIP_ARMS, ZIP_FILL

# Parameters ------------------------------------------------------------------------------
CUT <- 1 # fixed cytokine cut, arcsinh (cofactor 6000), both cytokines (user, 2026-11-10)
CYTOKINES <- c(TNF = "BV510-A", IFNg = "BV570-A") # as in 19a_gate_arm.R; arcsinh in the GatingSet
# Every activation or effector marker in the panel (CD69, CD154, CD137 and GzmB are not in it). Rank
# correlation is invariant to the monotone transform, so these channels are used as stored (linear).
ACTIVATION <- c(
  `IL-2` = "BV605-A",
  CD25 = "SBUV605-A",
  CD71 = "Spark YG 593-A",
  CD38 = "SBUV740-A",
  `CTLA-4` = "BV650-A",
  `PD-1` = "PE-Fire 640-A"
)
POPULATIONS <- c("CD4T", "CD8T")
ARM_ORDER <- unlist(VERSION_ROWS, use.names = FALSE)
GROUP_BREAKS <- c(1.5, 5.5, 9.5) # dashed lines between the version groups

cli_args <- commandArgs(trailingOnly = TRUE)
gated_dir <- cli_args[[1]]
run_dir <- cli_args[[2]]
stopifnot("46 run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir)
manifest <- read_csv(
  file.path(gated_dir, "manifest.csv"),
  show_col_types = FALSE
)

# Per arm, donor and parent: % TNF+ of IFNg+ and Spearman rho -----------------------------------
results <- ARM_ORDER |>
  set_names() |>
  map(\(a) {
    gs <- load_gs(file.path(gated_dir, "gs", a))
    sampleNames(gs) |>
      set_names() |>
      map(\(s) {
        POPULATIONS |>
          set_names() |>
          map(\(pop) {
            x <- exprs(gh_pop_get_data(gs[[s]], pop))
            ifng_pos <- x[, CYTOKINES[["IFNg"]]] > CUT
            rho <- CYTOKINES |>
              map(\(cy_ch) {
                map_dbl(ACTIVATION, \(m) {
                  cor(x[, cy_ch], x[, m], method = "spearman")
                })
              })
            list(
              pct = tibble(
                n_ifng_pos = sum(ifng_pos),
                pct_tnf_of_ifng = 100 *
                  mean(x[ifng_pos, CYTOKINES[["TNF"]]] > CUT),
                pct_tnf_of_ifng_neg = 100 *
                  mean(x[!ifng_pos, CYTOKINES[["TNF"]]] > CUT)
              ),
              rho = rho |>
                map(\(r) enframe(r, name = "marker", value = "rho")) |>
                list_rbind(names_to = "cytokine")
            )
          })
      })
  })

# Long tables ------------------------------------------------------------------------------
# flatten arm -> donor -> parent into rows
tnf_among_ifng <- results |>
  map(\(by_s) {
    map(by_s, \(by_p) list_rbind(map(by_p, "pct"), names_to = "parent")) |>
      list_rbind(names_to = "donor")
  }) |>
  list_rbind(names_to = "arm")
rho_activation <- results |>
  map(\(by_s) {
    map(by_s, \(by_p) list_rbind(map(by_p, "rho"), names_to = "parent")) |>
      list_rbind(names_to = "donor")
  }) |>
  list_rbind(names_to = "arm")
write_csv(tnf_among_ifng, file.path(run_dir, "tnf_among_ifng.csv"))
write_csv(rho_activation, file.path(run_dir, "rho_activation.csv"))

# Figure a: % TNF+ among IFNg+ cells, donors joined by lines ---------------------------------
fig_a <- tnf_among_ifng |>
  mutate(arm = factor(arm, levels = ARM_ORDER))
write_csv(fig_a, file.path(run_dir, "fig_a_source.csv"))
p_a <- ggplot(fig_a, aes(arm, pct_tnf_of_ifng)) +
  annotate(
    "rect",
    xmin = match(ZIP_ARMS, ARM_ORDER) - 0.5,
    xmax = match(ZIP_ARMS, ARM_ORDER) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = ZIP_FILL
  ) +
  geom_vline(
    xintercept = GROUP_BREAKS,
    linetype = "dashed",
    colour = "grey60"
  ) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey40") +
  geom_line(aes(group = donor, colour = donor), alpha = 0.6, linewidth = 0.4) +
  geom_point(aes(colour = donor), size = 1.8) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  scale_x_discrete(labels = ARM_LABELS) +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(~parent) +
  labs(x = NULL, y = "% TNF+ of IFNg+ cells") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
ggsave(
  file.path(run_dir, "fig_a_tnf_among_ifng.png"),
  p_a,
  width = 11,
  height = 5.5,
  dpi = 130,
  bg = "white"
)

# Figure b: Spearman rho with activation markers, TNF next to IFNg ---------------------------
# donor median per arm, parent, cytokine and marker; donors drawn as small points behind it
fig_b <- rho_activation |>
  mutate(
    arm = factor(arm, levels = rev(ARM_ORDER)),
    cytokine = factor(cytokine, levels = names(CYTOKINES)),
    marker = factor(marker, levels = names(ACTIVATION))
  )
write_csv(fig_b, file.path(run_dir, "fig_b_source.csv"))
fig_b_median <- fig_b |>
  summarise(rho = median(rho), .by = c(arm, parent, cytokine, marker))
p_b <- ggplot(fig_b, aes(rho, arm, colour = cytokine)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(size = 0.9, alpha = 0.4, position = position_dodge(width = 0.6)) +
  geom_point(
    data = fig_b_median,
    shape = 18,
    size = 3,
    position = position_dodge(width = 0.6)
  ) +
  scale_colour_manual(
    values = c(TNF = "#D55E00", IFNg = "#0072B2"),
    name = NULL
  ) +
  scale_y_discrete(labels = ARM_LABELS) +
  facet_grid(parent ~ marker) +
  labs(x = "Spearman rho (diamond = donor median)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(
      size = 9,
      face = if_else(rev(ARM_ORDER) %in% ZIP_ARMS, "bold", "plain")
    ),
    axis.text.x = element_text(size = 8),
    legend.position = "top"
  )
ggsave(
  file.path(run_dir, "fig_b_rho_activation.png"),
  p_b,
  width = 14,
  height = 8,
  dpi = 130,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
