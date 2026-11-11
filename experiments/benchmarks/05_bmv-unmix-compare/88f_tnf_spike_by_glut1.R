# Purpose: grill 2026-11-11 Q6 free check (README row "Co-spike (grill Q6, Q7; script 88d)"). Does TNF
#          spike loss track the cell's natural GLUT1 (Alexa Fluor 405) level? The TNF +1 arcsinh cells of
#          the five-setup spike-in are split into GLUT1 quartiles per donor (arcsinh, cofactor 6000) of the
#          1.8.8 zip unspiked unmix, so a cell has the same bin in every setup. Per bin, setup and T cell
#          population: % of cells losing >50% of the spike and median recovery, donors as points; plus a
#          ggcyto hexbin of TNF before vs after per GLUT1 bin. Natural GLUT1 is not randomised: bins also
#          differ in cell state, so a trend is an association.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/88f_tnf_spike_by_glut1.R  (repo root)
# Inputs:  experiments/runs/2026-11-10_bmv-spike-in-v2/spike_cells.csv   TNF level-1 cells, 5 setups [88]
#          experiments/bmv_test/compare_3D/data/unmix_cells/spk_c188_base/*.fcs  1.8.8 zip unspiked unmix [88]
# Outputs: <out_dir>/cells_by_glut1.csv          per TNF level-1 cell: GLUT1 bin, recovery
#          <out_dir>/summary_by_glut1.csv        per setup, pop, donor, bin: % lost >50%, median recovery
#          <out_dir>/fig_loss_by_glut1.{pdf,png,csv}   % lost per bin, donors as points
#          <out_dir>/fig_raw_by_glut1_<pop>.{pdf,png}  ggcyto hexbin TNF before vs after, bins x setups
#          <out_dir>/sessionInfo.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto
suppressPackageStartupMessages({
  library(flowCore)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lib/plots.R")

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- "experiments/runs/2026-11-11_bmv-tnf-spike-by-glut1"
SPIKE_CELLS <- "experiments/runs/2026-11-10_bmv-spike-in-v2/spike_cells.csv"
BASE_DIR <- "experiments/bmv_test/compare_3D/data/unmix_cells/spk_c188_base"
GLUT1_CH <- "Alexa Fluor 405-A"
LEVEL <- 1 # the +1 arcsinh cells, as the summary figures (user, 2026-11-10)
N_BINS <- 4 # quartiles of natural GLUT1 per donor (user, 2026-11-11: the default)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md)
LOSS_FRAC <- 0.5 # "losing most of the spike", as 88
SETUP_LABELS <- c(
  c188 = "1.8.8 zip",
  afrefine = "refined AF",
  opt188 = "opt188",
  s188_c102 = "1.0.2 code,\n1.8.8 spectra",
  s102_c102 = "1.0.2 zip"
)
AXIS <- c(-1.5, 4) # display window only
BINS <- 50 # hexagons per panel; display only

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)

# TNF level-1 cells with their natural GLUT1 from the 1.8.8 zip unspiked unmix ---------------------------
tnf <- read_csv(SPIKE_CELLS, show_col_types = FALSE) |>
  dplyr::filter(marker == "TNF", level == LEVEL)
glut1 <- distinct(tnf, raw_file) |>
  mutate(
    glut1_raw = map(raw_file, \(f) {
      ff <- read.FCS(
        file.path(BASE_DIR, sub("\\.fcs$", " AutoSpectral.fcs", f)),
        transformation = FALSE,
        truncate_max_range = FALSE
      )
      exprs(ff)[, GLUT1_CH]
    })
  )
cells <- tnf |>
  left_join(
    glut1,
    by = "raw_file",
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(
    glut1 = asinh(map2_dbl(glut1_raw, raw_row, \(g, r) g[[r]]) / COFACTOR)
  ) |>
  select(-glut1_raw) |>
  # quartile bins per donor, set on the cells of one setup so every setup shares the same cut points
  mutate(
    glut1_bin = cut(
      glut1,
      quantile(
        glut1[setup == "c188"],
        probs = seq(0, 1, length.out = N_BINS + 1)
      ),
      labels = paste0("GLUT1 Q", seq_len(N_BINS)),
      include.lowest = TRUE
    ),
    .by = sample
  ) |>
  mutate(
    before = asinh(base_value / COFACTOR),
    after = asinh((base_value + recovered) / COFACTOR),
    setup_label = factor(SETUP_LABELS[setup], levels = SETUP_LABELS)
  )
write_csv(
  select(cells, sample, pop, raw_row, setup, glut1, glut1_bin, recovery),
  file.path(OUT_DIR, "cells_by_glut1.csv")
)

summary_by_glut1 <- cells |>
  summarise(
    pct_lost = 100 * mean(recovery < LOSS_FRAC),
    median_recovery = median(recovery),
    glut1_median = median(glut1),
    n = n(),
    .by = c(setup, setup_label, pop, sample, glut1_bin)
  )
write_csv(
  select(summary_by_glut1, -setup_label),
  file.path(OUT_DIR, "summary_by_glut1.csv")
)

# Figure 1: % losing >50% per GLUT1 bin, donors as points, setups as panels -------------------------------
loss_df <- summary_by_glut1 |>
  mutate(donor = sample)
pal_setup <- categorical_palette(levels(loss_df$setup_label))
p_loss <- ggplot(loss_df, aes(glut1_bin, pct_lost, colour = setup_label)) +
  stat_summary(
    aes(group = setup_label),
    fun = median,
    geom = "line",
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_point(
    aes(shape = donor),
    position = position_jitter(width = 0.12, height = 0, seed = 1)
  ) +
  scale_colour_manual(values = pal_setup, guide = "none") +
  facet_grid(pop ~ setup_label) +
  labs(
    x = "Natural GLUT1 (Alexa Fluor 405), quartile per donor",
    y = "Cells losing >50% of the TNF spike (%)",
    shape = "Donor",
    title = "TNF +1 arcsinh spike by natural GLUT1"
  ) +
  rotate_x()
save_panel(
  p_loss,
  loss_df,
  OUT_DIR,
  "fig_loss_by_glut1",
  width_mm = 183,
  height_mm = 110
)

# Figure 2 per population: ggcyto hexbin of TNF before vs after, GLUT1 bins x setups ----------------------
grid_x <- seq(AXIS[1], AXIS[2], length.out = 200)
curves <- bind_rows(
  tibble(
    x = grid_x,
    y = asinh(sinh(grid_x) + sinh(LEVEL)),
    curve = "full recovery"
  ),
  tibble(x = grid_x, y = grid_x, curve = "spike lost")
)
raw_figure <- \(which_pop) {
  d <- dplyr::filter(cells, pop == which_pop)
  frames <- d |>
    mutate(name = paste(setup_label, glut1_bin, sep = "|")) |>
    split(~name) |>
    map(\(x) flowFrame(cbind(before = x$before, after = x$after)))
  fs <- flowSet(frames)
  pd <- tibble(name = sampleNames(fs)) |>
    tidyr::separate_wider_delim(
      name,
      "|",
      names = c("setup_label", "glut1_bin"),
      cols_remove = FALSE
    ) |>
    mutate(setup_label = factor(setup_label, levels = SETUP_LABELS))
  pData(fs) <- as.data.frame(pd) |>
    `rownames<-`(pd$name)
  labels <- summary_by_glut1 |>
    dplyr::filter(pop == which_pop) |>
    summarise(pct_lost = median(pct_lost), .by = c(setup_label, glut1_bin))
  p <- ggcyto(fs, aes(x = before, y = after)) +
    geom_hex(bins = BINS) +
    geom_line(
      data = curves,
      aes(x, y, linetype = curve),
      inherit.aes = FALSE,
      colour = "grey30",
      linewidth = 0.4
    ) +
    geom_label(
      data = labels,
      aes(x = -Inf, y = Inf, label = sprintf("lost >50%%: %.1f%%", pct_lost)),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 7,
      size.unit = "pt",
      fontface = "bold",
      fill = "white",
      alpha = 0.85,
      label.size = 0
    ) +
    facet_grid(glut1_bin ~ setup_label, switch = "y") +
    ggcyto_par_set(limits = list(x = AXIS, y = AXIS)) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      trans = "sqrt",
      guide = "none"
    ) +
    scale_linetype_manual(
      values = c("full recovery" = "dashed", "spike lost" = "solid"),
      name = NULL
    ) +
    labs(
      x = "TNF, unspiked unmix (arcsinh)",
      y = "TNF, +1 arcsinh spiked unmix",
      title = sprintf(
        "%s: TNF spike per cell by natural GLUT1 quartile (5 donors pooled; %% = donor median)",
        which_pop
      )
    ) +
    ggpubr::theme_pubr(base_size = 8, border = TRUE) +
    theme(
      strip.text.x = element_text(size = 7),
      strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(1, "mm")
    )
  stem <- file.path(OUT_DIR, sprintf("fig_raw_by_glut1_%s", which_pop))
  ggsave(
    paste0(stem, ".pdf"),
    as.ggplot(p),
    width = 183,
    height = 150,
    units = "mm",
    device = cairo_pdf
  )
  ggsave(
    paste0(stem, ".png"),
    as.ggplot(p),
    width = 183,
    height = 150,
    units = "mm",
    dpi = 200,
    bg = "white"
  )
}
raw_figure("CD4T")
raw_figure("CD8T")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
