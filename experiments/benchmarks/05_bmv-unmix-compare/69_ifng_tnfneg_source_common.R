# Purpose: Common-gates copy (plans/common_gates_plan.md, batch 3): gates are the common SpectroFlo labels and the
#          cytokine cut is arcsinh 1.5 (G7a). Traces the IFNg+ TNF- T cells of AutoSpectral 1.0.2 (as102_var) (user, 2026-11-10). The same events
#          are followed into 1.8.8 (as188_joint), the AF405-only swap (var102_only_Alexa_Fluor_405: 1.8.8 with
#          only the Alexa Fluor 405 variants from 1.0.2) and the knock-out (af405proof_ko_bv510: those variants
#          minus their BV510 part). Asks: are these cells TNF+ elsewhere, and are they GLUT1-high?
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/69_ifng_tnfneg_source_common.R   (repo root)
#          Needs the TNF and IFNg readouts: stops when the dataset's targets are not TNF and IFNg (EXV).
# Inputs:  run_in("common-gated")/gs/<arm>/  common-gate GatingSets; rows align across arms [batch 0, 97d]
#          datasets/<SW_DATASET>.json (channels, targets, conditions, pos_rule), <ds>_samples.csv
#          run_in("common-cuts")/cuts.csv when the dataset's positive rule is not fixed [97e]
# Outputs: run_out("ifng-tnfneg-common")/cells.csv.gz          per IFNg+ cell: TNF per arm, GLUT1, AF (donor = sample id)
#          run_out("ifng-tnfneg-common")/share_tnfpos.csv      per sample, pop, arm: share TNF+
#          run_out("ifng-tnfneg-common")/fig_a_tnf_188_vs_102.{pdf,png,csv}  TNF 1.8.8 vs 1.0.2, IFNg+ cells
#          run_out("ifng-tnfneg-common")/fig_b_share_tnfpos.{pdf,png,csv}   1.0.2 IFNg+TNF- cells TNF+ elsewhere
#          run_out("ifng-tnfneg-common")/fig_c_glut1_af.{pdf,png,csv}       GLUT1 and AF, TNF- vs TNF+ (1.0.2)
#          run_out("ifng-tnfneg-common")/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ---------------------------------------------------------------------------------------
GATED <- file.path(run_in("common-gated"), "gs") # batch 0: SpectroFlo labels as gates
TNF_LABEL <- cfg("targets")$tnf_label
IFNG_LABEL <- cfg("targets")$ifng_label
if (!identical(c(TNF_LABEL, IFNG_LABEL), c("TNF", "IFNg"))) {
  stop(
    "69_ifng_tnfneg_source_common.R: needs TNF and IFNg, missing in dataset ",
    DS$dataset,
    " (the two dye slots are ",
    TNF_LABEL,
    " and ",
    IFNG_LABEL,
    ")",
    call. = FALSE
  )
}
TNF_SLOT <- cfg("channels")$TNF_slot # stored as arcsinh in the GatingSets
IFNG_SLOT <- cfg("channels")$IFNg_slot
GLUT1_SLOT <- paste0(cfg("targets")$glut1_dye, "-A") # linear in the GatingSets
REF_ARM <- "as188_joint" # its CD4T/CD8T gates define the cells in every arm (as in 51)
ARM_102 <- "as102_var" # AutoSpectral 1.0.2 zip output
COMPARE_ARMS <- c(
  "as188_joint" = "1.8.8",
  "var102_only_Alexa_Fluor_405" = "1.8.8 + 1.0.2 AF405 variants",
  "af405proof_ko_bv510" = "same, BV510 part removed"
)
POPS <- c("CD4T", "CD8T") # user, 2026-11-10
# cytokine cuts: cut_of(arm, sample, channel, "pos"); 3D is the fixed 1.5 in every arm (G7a, user 2026-11-11; was 1)
# the positive cut in the axis label: the fixed value, or the per-donor rule
POS_TEXT <- if (identical(cfg("pos_rule")$type, "fixed")) {
  sprintf("arcsinh > %s", cfg("pos_rule")$value)
} else {
  "above the positive cut"
}
COFACTOR <- 6000 # arcsinh cofactor for the linear GLUT1 and AF channels (CLAUDE.md)
OUT_DIR <- run_out("ifng-tnfneg-common") # v2: GLUT1 colour scale -1 to 4 (user)
SAMPLES <- samples_tbl()
COND_KEY <- if (is.null(DS$conditions)) character() else "condition"
# 24H figures plot the donor over conditions (columns plot_donor, condition); 3D and EXV plot the sample
UNIT_COLS <- if (is.null(DS$conditions)) {
  character()
} else {
  c("plot_donor", "condition")
}
COLOUR_VAR <- if (is.null(DS$conditions)) "donor" else "plot_donor"
add_plot_unit <- \(d) {
  if (is.null(DS$conditions)) {
    d
  } else {
    left_join(
      d,
      select(SAMPLES, donor = sample, plot_donor = donor, condition),
      by = "donor",
      relationship = "many-to-one",
      unmatched = c("error", "drop")
    )
  }
}

stopifnot("output dir must not exist" = !dir.exists(OUT_DIR))
dir.create(OUT_DIR)
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R") # theme_set(), save_panel()
source("experiments/benchmarks/05_bmv-unmix-compare/lib/palettes.R") # categorical_palette() for more than 8 donors

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
  })),
  "gating sets hold the samples of samples_tbl()" = setequal(
    donors,
    SAMPLES$sample
  )
)
# positive cuts per arm and sample, read once (cut_of re-reads cuts.csv when the rule is not fixed)
cut_tbl <- tidyr::expand_grid(
  arm = arms,
  sample = donors,
  channel = c(TNF_SLOT, IFNG_SLOT)
) |>
  mutate(
    pos_cut = pmap_dbl(
      list(arm, sample, channel),
      \(a, s, ch) cut_of(a, s, ch, "pos")
    )
  )
# cut of one arm and channel by sample id
cut_by_sample <- \(a, ch) {
  hit <- cut_tbl[cut_tbl$arm == a & cut_tbl$channel == ch, ]
  set_names(hit$pos_cut, hit$sample)
}

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
          ifng_102 = ev[[ARM_102]][rows, IFNG_SLOT],
          tnf_102 = ev[[ARM_102]][rows, TNF_SLOT],
          glut1 = asinh(ev[[REF_ARM]][rows, GLUT1_SLOT] / COFACTOR),
          af = asinh(ev[[REF_ARM]][rows, "AF-A"] / COFACTOR),
          !!!map(names(COMPARE_ARMS) |> set_names(), \(a) {
            ev[[a]][rows, TNF_SLOT]
          })
        )
      }) |>
      list_rbind(names_to = "pop")
  }) |>
  list_rbind(names_to = "donor") |>
  filter(ifng_102 > cut_by_sample(ARM_102, IFNG_SLOT)[donor]) |>
  mutate(
    tnf_call_102 = if_else(
      tnf_102 > cut_by_sample(ARM_102, TNF_SLOT)[donor],
      "TNF+ in 1.0.2",
      "TNF- in 1.0.2"
    )
  ) |>
  add_plot_unit()
write_csv(cells, file.path(OUT_DIR, "cells.csv.gz"))

# Share of 1.0.2 IFNg+ TNF- cells that are TNF+ in each comparison arm ---------------------------
share_tnfpos <- cells |>
  filter(tnf_call_102 == "TNF- in 1.0.2") |>
  pivot_longer(
    all_of(names(COMPARE_ARMS)),
    names_to = "arm",
    values_to = "tnf"
  ) |>
  left_join(
    filter(cut_tbl, channel == TNF_SLOT),
    join_by(donor == sample, arm),
    relationship = "many-to-one",
    unmatched = c("error", "drop")
  ) |>
  summarise(
    n_cells = n(),
    share_tnfpos = mean(tnf > pos_cut),
    .by = c(donor, pop, arm, any_of(UNIT_COLS))
  ) |>
  mutate(arm_label = factor(COMPARE_ARMS[arm], levels = unname(COMPARE_ARMS)))
write_csv(share_tnfpos, file.path(OUT_DIR, "share_tnfpos.csv"))

# Figure a: TNF of the same IFNg+ cells, 1.8.8 vs 1.0.2, coloured by GLUT1 ---------------------------
frame_a <- cells |>
  mutate(across(c(as188_joint, tnf_102), \(x) pmin(pmax(x, -1.5), 4))) |>
  arrange(glut1)
cut_lines <- tibble(
  donor = donors,
  x_cut = cut_by_sample(REF_ARM, TNF_SLOT)[donors],
  y_cut = cut_by_sample(ARM_102, TNF_SLOT)[donors]
) |>
  add_plot_unit()
p_a <- ggplot(frame_a, aes(as188_joint, tnf_102, colour = glut1)) +
  geom_point(size = 0.1, alpha = 0.5, stroke = 0) +
  geom_hline(
    data = cut_lines,
    aes(yintercept = y_cut),
    colour = "red",
    linewidth = 0.3
  ) +
  geom_vline(
    data = cut_lines,
    aes(xintercept = x_cut),
    colour = "red",
    linewidth = 0.3
  ) +
  scale_colour_viridis_c(
    name = "GLUT1\n(arcsinh)",
    limits = c(-1, 4), # display range (user, 2026-11-10); values outside are squished to the ends
    oob = scales::squish
  ) +
  facet_grid(as.formula(paste(
    paste(c("pop", COND_KEY), collapse = " + "),
    "~",
    COLOUR_VAR
  ))) +
  coord_fixed(xlim = c(-1.5, 4), ylim = c(-1.5, 4)) +
  labs(
    x = sprintf("%s in 1.8.8 (as188_joint), arcsinh", TNF_LABEL),
    y = sprintf("%s in 1.0.2 (as102_var), arcsinh", TNF_LABEL)
  )
save_panel(p_a, frame_a, OUT_DIR, "fig_a_tnf_188_vs_102", 183, 90)

# donor colours: Dark2 holds 8, 24H has more donors
donor_scale <- \(d) {
  if (n_distinct(d[[COLOUR_VAR]]) <= 8) {
    scale_colour_brewer(palette = "Dark2", name = "Donor")
  } else {
    scale_colour_manual(
      values = categorical_palette(sort(unique(d[[COLOUR_VAR]]))),
      name = "Donor"
    )
  }
}

# Figure b: where the 1.0.2 TNF- cells are TNF+ --------------------------------------------------------
p_b <- ggplot(share_tnfpos, aes(arm_label, share_tnfpos)) +
  geom_line(
    aes(group = donor, colour = .data[[COLOUR_VAR]]),
    alpha = 0.6,
    linewidth = 0.3
  ) +
  geom_point(aes(colour = .data[[COLOUR_VAR]]), size = 1) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, linewidth = 0.3) +
  donor_scale(share_tnfpos) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(c("pop", COND_KEY)) +
  labs(x = NULL, y = sprintf("Share %s+ (%s)", TNF_LABEL, POS_TEXT)) +
  scale_x_discrete(labels = scales::label_wrap(12))
save_panel(p_b, share_tnfpos, OUT_DIR, "fig_b_share_tnfpos", 183, 80)

# Figure c: GLUT1 and AF of IFNg+ cells, TNF- vs TNF+ in 1.0.2, donor medians -----------------------
frame_c <- cells |>
  pivot_longer(c(glut1, af), names_to = "marker", values_to = "value") |>
  summarise(
    median_value = median(value),
    n_cells = n(),
    .by = c(donor, pop, marker, tnf_call_102, any_of(UNIT_COLS))
  ) |>
  mutate(
    marker = factor(
      marker,
      levels = c("glut1", "af"),
      labels = c("GLUT1 (AF405, arcsinh)", "AF-A (arcsinh)")
    )
  )
p_c <- ggplot(frame_c, aes(tnf_call_102, median_value)) +
  geom_line(aes(group = donor, colour = .data[[COLOUR_VAR]]), linewidth = 0.3) +
  geom_point(aes(colour = .data[[COLOUR_VAR]]), size = 1) +
  donor_scale(frame_c) +
  facet_grid(
    as.formula(paste("marker ~", paste(c("pop", COND_KEY), collapse = " + "))),
    scales = "free_y"
  ) +
  labs(x = NULL, y = "Donor median") +
  scale_x_discrete(labels = scales::label_wrap(12))
save_panel(p_c, frame_c, OUT_DIR, "fig_c_glut1_af", 183, 90)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
