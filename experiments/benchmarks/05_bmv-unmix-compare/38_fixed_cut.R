# Purpose: compare % TNF+ and % IFNg+ across unmixing arms with ONE fixed cut per cytokine for every arm
#          and donor, so differences come from unmixing, not from gate placement. Adds the fixed gates to
#          the saved GatingSets in memory (not written back), then draws (1) boxplots per arm with donors
#          joined by lines and (2) ggcyto TNF vs IFNg hexbins with the fixed gates (geom_gate).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/38_fixed_cut.R <run_dir> <cut> [out_root]
#          from the repo root, after 19a_gate_arm.R has run for every arm in <run_dir>/manifest.csv.
#          <cut> is the arcsinh cut used for both cytokines, e.g. 1.5 or 1.
#          [out_root] is where figures_fixed_cut_<cut>/ goes; default <run_dir>.
# Inputs:  <run_dir>/manifest.csv [01_manifest-gate.R]
#          <run_dir>/gs/<arm>/ saved GatingSets with CD4T and CD8T [19a_gate_arm.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  labels, VERSION_ROWS, version_grid()
# Outputs: <out_root>/figures_fixed_cut_<cut>/pct_fixed_paired.png       boxplots with paired donor lines
#          <out_root>/figures_fixed_cut_<cut>/pct_fixed_source.csv       % above the cut per arm, donor, parent
#          <out_root>/figures_fixed_cut_<cut>/tnf_vs_ifng_fixed_<population>_<sample>.png  hexbins, fixed gates
#          <out_root>/figures_fixed_cut_<cut>/sessionInfo_38.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto

# Bioconductor first, tidyverse parts last; dplyr::filter() is called by name.
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
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS, VERSION_ROWS, version_grid()

# Parameters ------------------------------------------------------------------------------
# One fixed cut per cytokine, arcsinh(x / 6000), same in every arm and donor (user, 2026-11-10): automatic
# gates landed in arbitrary places (TNF has no valley; IFNg tailgate failed in several arms).
# The cut comes from the command line (1.5 approved; 1 asked as a check by the user, 2026-11-10).
cli_args <- commandArgs(trailingOnly = TRUE)
CUT <- as.numeric(cli_args[[2]])
CUTS <- c(TNF = CUT, IFNg = CUT)
CHANNELS <- c(TNF = "BV510-A", IFNg = "BV570-A") # as in 19a_gate_arm.R
POPULATIONS <- c("CD4T", "CD8T")
BINS <- 100 # hexagon bins; display only
AXIS_LIMITS <- c(-0.5, 5) # arcsinh display range (user, 2026-11-09: trims are plots only)
ARM_ORDER <- unlist(VERSION_ROWS, use.names = FALSE) # the four code-version groups, in figure order
GROUP_BREAKS <- c(1.5, 5.5, 9.5) # dashed lines between the version groups

run_dir <- cli_args[[1]]
out_root <- if (length(cli_args) > 2) cli_args[[3]] else run_dir
# the first run (cut 1.5) wrote figures_fixed_cut/; later cuts get their own folder
fig_dir <- file.path(out_root, paste0("figures_fixed_cut_", CUT))
stopifnot("38 outputs must not exist" = !dir.exists(fig_dir))
dir.create(fig_dir, recursive = TRUE)
manifest <- read_csv(file.path(run_dir, "manifest.csv"), show_col_types = FALSE)
arms <- unique(manifest$arm)
stopifnot("every arm sits in one version row" = setequal(ARM_ORDER, arms))

# Load the GatingSets and add the fixed cytokine gates (in memory) -----------------------------
# Each gate is a 1D rectangle from the cut to +Inf, named <cytokine>fix, under CD4T and CD8T.
gating_sets <- arms |>
  set_names() |>
  map(\(a) {
    gs <- load_gs(file.path(run_dir, "gs", a))
    for (pop in POPULATIONS) {
      iwalk(CUTS, \(cut, cy) {
        bounds <- list(c(cut, Inf)) |> set_names(CHANNELS[[cy]])
        gs_pop_add(
          gs,
          rectangleGate(filterId = paste0(cy, "fix"), .gate = bounds),
          parent = pop
        )
      })
    }
    recompute(gs)
    gs
  })

# % above the fixed cut per arm, donor, parent and cytokine -----------------------------------
pct <- gating_sets |>
  imap(\(gs, arm) {
    gs_pop_get_count_fast(gs) |> as_tibble() |> mutate(arm = arm)
  }) |>
  list_rbind() |>
  dplyr::filter(str_detect(Population, "/(CD4T|CD8T)/(TNF|IFNg)fix$")) |>
  transmute(
    arm = factor(arm, levels = ARM_ORDER),
    donor = name,
    parent = str_extract(Population, "CD[48]T"),
    cytokine = factor(
      str_extract(Population, "(TNF|IFNg)(?=fix$)"),
      levels = names(CUTS)
    ),
    n_parent = ParentCount,
    pct_pos = 100 * Count / ParentCount
  )
stopifnot(
  "one value per arm, donor, parent, cytokine" = nrow(pct) ==
    length(ARM_ORDER) * 5 * 2 * 2
)
write_csv(pct, file.path(fig_dir, "pct_fixed_source.csv"))

# Figure 1: boxplots per arm, donors joined by lines ------------------------------------------
p_box <- ggplot(pct, aes(arm, pct_pos)) +
  annotate(
    "rect",
    xmin = match(ZIP_ARMS, ARM_ORDER) - 0.5,
    xmax = match(ZIP_ARMS, ARM_ORDER) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = ZIP_FILL
  ) + # zip arms shaded
  geom_vline(
    xintercept = GROUP_BREAKS,
    linetype = "dashed",
    colour = "grey70"
  ) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey40") +
  geom_line(aes(group = donor, colour = donor), alpha = 0.6, linewidth = 0.4) +
  geom_point(aes(colour = donor), size = 1.8) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  scale_x_discrete(labels = ARM_LABELS) +
  facet_grid(parent ~ cytokine, scales = "free_y") +
  labs(
    x = NULL,
    y = "% positive of parent",
    title = paste0("Cut at arcsinh ", CUT)
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(fig_dir, "pct_fixed_paired.png"),
  p_box,
  width = 12,
  height = 10,
  dpi = 130,
  bg = "white"
)

# Figure 2: TNF vs IFNg hexbin per donor and parent, fixed gates drawn ---------------------------
pct_key <- pct |>
  mutate(key = paste(arm, donor, parent, cytokine)) |>
  (\(d) set_names(d$pct_pos, d$key))()
for (s in sort(unique(manifest$sample))) {
  for (pop in POPULATIONS) {
    panels <- gating_sets |>
      imap(\(gs, arm) {
        ggcyto(
          gs[s],
          aes(
            x = !!rlang::sym(CHANNELS[["TNF"]]),
            y = !!rlang::sym(CHANNELS[["IFNg"]])
          ),
          subset = pop
        ) +
          geom_hex(bins = BINS) +
          geom_gate(paste0(pop, "/", c("TNFfix", "IFNgfix"))) +
          ggcyto_par_set(limits = list(x = AXIS_LIMITS, y = AXIS_LIMITS)) +
          labs(
            title = ARM_LABELS[[arm]],
            subtitle = sprintf(
              "TNF+ %.1f%%  IFNg+ %.1f%%",
              pct_key[[paste(arm, s, pop, "TNF")]],
              pct_key[[paste(arm, s, pop, "IFNg")]]
            ),
            x = "TNF (arcsinh)",
            y = "IFNg (arcsinh)"
          ) +
          PANEL_THEME
      })
    ggsave(
      file.path(fig_dir, paste0("tnf_vs_ifng_fixed_", pop, "_", s, ".png")),
      version_grid(panels),
      width = 14,
      height = 10,
      dpi = 120,
      bg = "white"
    )
  }
}
writeLines(
  capture.output(sessionInfo()),
  file.path(fig_dir, "sessionInfo_38.txt")
)
