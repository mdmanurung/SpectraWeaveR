# Purpose: compare % cytokine-positive per parent population (TNF+ and IFNg+ in CD4 T and CD8 T) across
#          unmixing arms: boxplot per arm over donors, donors as points joined by lines across arms
#          (paired, same raw events). Each arm uses its own gates from 19a_gate_arm.R.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/32_pct_pos_paired.R <run_dir> [out_root]
#          from the repo root. [out_root] is where the figure folder goes; default <run_dir>.
# Inputs:  <run_dir>/fsc_cyto_gate_counts.csv  events per gate, arm and donor [19_fsc_vs_cytokines.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  arm display labels
# Outputs: <out_root>/figures_pct_pos/pct_pos_paired.png  boxplots with paired donor lines
#          <out_root>/figures_pct_pos/pct_pos_paired_source.csv  plotted values
#          <out_root>/figures_pct_pos/sessionInfo_32.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R") # ARM_LABELS

# Parameters ------------------------------------------------------------------------------
# x order: the four code-version groups used in every figure of this benchmark
ARM_ORDER <- c(
  "spectroflo",
  "as102_basic",
  "as102_var",
  "best102_neg",
  "best102_af",
  "as188_joint",
  "as188_fast_csv",
  "best188_p1",
  "best188_p5",
  "s188_c102",
  "s102_c188"
)
GROUP_BREAKS <- c(1.5, 5.5, 9.5) # dashed lines between the version groups
# Failed own gate, display only: below 0.5% positive (the failed B11 CD8T TNF gates sit at 0.02%), or
# IFNg+ above 50% (IFNg+ is 15-30% where the tailgate sits right; TNF+ is often above 50% legitimately)
FAIL_PCT_LOW <- 0.5
IFNG_FAIL_PCT <- 50

cli_args <- commandArgs(trailingOnly = TRUE)
run_dir <- cli_args[[1]]
out_root <- if (length(cli_args) > 1) cli_args[[2]] else run_dir
fig_dir <- file.path(out_root, "figures_pct_pos")
stopifnot("32 outputs must not exist" = !dir.exists(fig_dir))
dir.create(fig_dir, recursive = TRUE)

# Plot frame: % positive per arm, donor, parent and cytokine -----------------------------------
pct <- read_csv(
  file.path(run_dir, "fsc_cyto_gate_counts.csv"),
  show_col_types = FALSE
) |>
  filter(str_detect(Population, "/(CD4T|CD8T)/(TNF|IFNg)pos$")) |>
  transmute(
    arm = factor(arm, levels = ARM_ORDER),
    donor = name,
    parent = str_extract(Population, "CD[48]T"),
    cytokine = factor(
      str_extract(Population, "(TNF|IFNg)(?=pos$)"),
      levels = c("TNF", "IFNg")
    ),
    pct_pos = 100 * Count / ParentCount
  ) |>
  mutate(
    failed = pct_pos < FAIL_PCT_LOW |
      (cytokine == "IFNg" & pct_pos > IFNG_FAIL_PCT)
  )
stopifnot(
  "every arm is in ARM_ORDER" = !anyNA(pct$arm),
  "one value per arm, donor, parent, cytokine" = nrow(pct) ==
    length(ARM_ORDER) * 5 * 2 * 2
)
write_csv(pct, file.path(fig_dir, "pct_pos_paired_source.csv"))

# Figure ----------------------------------------------------------------------------------
p <- ggplot(pct, aes(arm, pct_pos)) +
  geom_vline(
    xintercept = GROUP_BREAKS,
    linetype = "dashed",
    colour = "grey70"
  ) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey40") +
  geom_line(aes(group = donor, colour = donor), alpha = 0.6, linewidth = 0.4) +
  geom_point(aes(colour = donor, shape = failed), size = 1.8, stroke = 0.7) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 1),
    name = "Failed gate"
  ) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  scale_x_discrete(labels = ARM_LABELS) +
  facet_grid(parent ~ cytokine, scales = "free_y") +
  labs(
    x = NULL,
    y = "% positive of parent",
    title = "Automatic gates, one per arm"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(fig_dir, "pct_pos_paired.png"),
  p,
  width = 12,
  height = 10,
  dpi = 130,
  bg = "white"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(fig_dir, "sessionInfo_32.txt")
)
