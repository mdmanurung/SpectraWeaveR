# Purpose: trace the IFNg+ TNF- T cells of AutoSpectral 1.0.2 (as102_var) (user, 2026-11-10). The same events
#          are followed into 1.8.8 (as188_joint), the AF405-only swap (var102_only_Alexa_Fluor_405: 1.8.8 with
#          only the Alexa Fluor 405 variants from 1.0.2) and the knock-out (af405proof_ko_bv510: those variants
#          minus their BV510 part). Asks: are these cells TNF+ elsewhere, and are they GLUT1-high?
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/69_ifng_tnfneg_source.R   (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-ablation-gated/gs/<arm>/  GatingSets; rows align across arms
#            (51_score_ablation.R masks every arm with as188_joint's cleanup rows) [19a_gate_arm.R via 51]
# Outputs: experiments/runs/2026-11-10_bmv-ifng-tnfneg-v2/cells.csv.gz          per IFNg+ cell: TNF per arm, GLUT1, AF
#          experiments/runs/2026-11-10_bmv-ifng-tnfneg-v2/share_tnfpos.csv      per donor, pop, arm: share TNF+
#          experiments/runs/2026-11-10_bmv-ifng-tnfneg-v2/fig_a_tnf_188_vs_102.png  TNF 1.8.8 vs 1.0.2, IFNg+ cells
#          experiments/runs/2026-11-10_bmv-ifng-tnfneg-v2/fig_b_share_tnfpos.png   1.0.2 IFNg+TNF- cells TNF+ elsewhere
#          experiments/runs/2026-11-10_bmv-ifng-tnfneg-v2/fig_c_glut1_af.png       GLUT1 and AF, TNF- vs TNF+ (1.0.2)
#          experiments/runs/2026-11-10_bmv-ifng-tnfneg-v2/sessionInfo.txt

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
REF_ARM <- "as188_joint" # its CD4T/CD8T gates define the cells in every arm (as in 51)
ARM_102 <- "as102_var" # AutoSpectral 1.0.2 zip output
COMPARE_ARMS <- c(
  "as188_joint" = "1.8.8",
  "var102_only_Alexa_Fluor_405" = "1.8.8 + 1.0.2 AF405 variants",
  "af405proof_ko_bv510" = "same, BV510 part removed"
)
POPS <- c("CD4T", "CD8T") # user, 2026-11-10
CUT <- 1 # fixed cytokine cut, arcsinh (user, 2026-11-10)
COFACTOR <- 6000 # arcsinh cofactor for the linear GLUT1 and AF channels (CLAUDE.md)
OUT_DIR <- "experiments/runs/2026-11-10_bmv-ifng-tnfneg-v2" # v2: GLUT1 colour scale -1 to 4 (user)

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# Events: same rows in every arm ---------------------------------------------------------------------
# BV510/BV570 are stored as arcsinh in the GatingSets; GLUT1 (Alexa Fluor 405) and AF-A are linear
arms <- c(ARM_102, names(COMPARE_ARMS))
gating_sets <- arms |>
  set_names() |>
  map(\(a) load_gs(file.path(GATED, a)))
donors <- sampleNames(gating_sets[[REF_ARM]])
stopifnot(
  "same samples in every arm" = all(map_lgl(gating_sets, \(g) {
    setequal(sampleNames(g), donors)
  }))
)

cells <- donors |>
  set_names() |>
  map(\(s) {
    ev <- map(gating_sets, \(g) exprs(gh_pop_get_data(g[[s]], "root")))
    stopifnot(
      "rows align across arms" = all(map_int(ev, nrow) == nrow(ev[[REF_ARM]]))
    )
    POPS |>
      set_names() |>
      map(\(pop) {
        rows <- gh_pop_get_indices(gating_sets[[REF_ARM]][[s]], pop)
        tibble(
          ifng_102 = ev[[ARM_102]][rows, "BV570-A"],
          tnf_102 = ev[[ARM_102]][rows, "BV510-A"],
          glut1 = asinh(ev[[REF_ARM]][rows, "Alexa Fluor 405-A"] / COFACTOR),
          af = asinh(ev[[REF_ARM]][rows, "AF-A"] / COFACTOR),
          !!!map(names(COMPARE_ARMS) |> set_names(), \(a) {
            ev[[a]][rows, "BV510-A"]
          })
        )
      }) |>
      list_rbind(names_to = "pop")
  }) |>
  list_rbind(names_to = "donor") |>
  filter(ifng_102 > CUT) |>
  mutate(
    tnf_call_102 = if_else(tnf_102 > CUT, "TNF+ in 1.0.2", "TNF- in 1.0.2")
  )
write_csv(cells, file.path(OUT_DIR, "cells.csv.gz"))

# Share of 1.0.2 IFNg+ TNF- cells that are TNF+ in each comparison arm ---------------------------
share_tnfpos <- cells |>
  filter(tnf_call_102 == "TNF- in 1.0.2") |>
  pivot_longer(
    all_of(names(COMPARE_ARMS)),
    names_to = "arm",
    values_to = "tnf"
  ) |>
  summarise(
    n_cells = n(),
    share_tnfpos = mean(tnf > CUT),
    .by = c(donor, pop, arm)
  ) |>
  mutate(arm_label = factor(COMPARE_ARMS[arm], levels = unname(COMPARE_ARMS)))
write_csv(share_tnfpos, file.path(OUT_DIR, "share_tnfpos.csv"))

# Figure a: TNF of the same IFNg+ cells, 1.8.8 vs 1.0.2, coloured by GLUT1 ---------------------------
frame_a <- cells |>
  mutate(across(c(as188_joint, tnf_102), \(x) pmin(pmax(x, -1.5), 4))) |>
  arrange(glut1)
p_a <- ggplot(frame_a, aes(as188_joint, tnf_102, colour = glut1)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_hline(yintercept = CUT, colour = "red", linewidth = 0.3) +
  geom_vline(xintercept = CUT, colour = "red", linewidth = 0.3) +
  scale_colour_viridis_c(
    name = "GLUT1\n(arcsinh)",
    limits = c(-1, 4), # display range (user, 2026-11-10); values outside are squished to the ends
    oob = scales::squish
  ) +
  facet_grid(pop ~ donor) +
  coord_fixed(xlim = c(-1.5, 4), ylim = c(-1.5, 4)) +
  labs(
    x = "TNF in 1.8.8 (as188_joint), arcsinh",
    y = "TNF in 1.0.2 (as102_var), arcsinh",
    title = "IFNg+ T cells (IFNg+ in 1.0.2): the same cells' TNF in both versions"
  ) +
  theme_bw(base_size = 8)
ggsave(
  file.path(OUT_DIR, "fig_a_tnf_188_vs_102.png"),
  p_a,
  width = 11,
  height = 5.2,
  dpi = 150,
  bg = "white"
)

# Figure b: where the 1.0.2 TNF- cells are TNF+ --------------------------------------------------------
p_b <- ggplot(share_tnfpos, aes(arm_label, share_tnfpos)) +
  geom_line(aes(group = donor, colour = donor), alpha = 0.6) +
  geom_point(aes(colour = donor), size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, linewidth = 0.3) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~pop) +
  labs(
    x = NULL,
    y = "Share TNF+ (arcsinh > 1)",
    title = "IFNg+ cells called TNF- by 1.0.2: share TNF+ in other arms (bar = donor median)"
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(
  file.path(OUT_DIR, "fig_b_share_tnfpos.png"),
  p_b,
  width = 8,
  height = 4.2,
  dpi = 150,
  bg = "white"
)

# Figure c: GLUT1 and AF of IFNg+ cells, TNF- vs TNF+ in 1.0.2, donor medians -----------------------
frame_c <- cells |>
  pivot_longer(c(glut1, af), names_to = "marker", values_to = "value") |>
  summarise(
    median_value = median(value),
    n_cells = n(),
    .by = c(donor, pop, marker, tnf_call_102)
  ) |>
  mutate(
    marker = factor(
      marker,
      levels = c("glut1", "af"),
      labels = c("GLUT1 (AF405, arcsinh)", "AF-A (arcsinh)")
    )
  )
p_c <- ggplot(frame_c, aes(tnf_call_102, median_value)) +
  geom_line(aes(group = donor, colour = donor)) +
  geom_point(aes(colour = donor), size = 2) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  facet_grid(marker ~ pop, scales = "free_y") +
  labs(
    x = NULL,
    y = "Donor median",
    title = "IFNg+ cells: GLUT1 and autofluorescence by 1.0.2 TNF call"
  ) +
  theme_bw(base_size = 9)
ggsave(
  file.path(OUT_DIR, "fig_c_glut1_af.png"),
  p_c,
  width = 6.5,
  height = 5,
  dpi = 150,
  bg = "white"
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
