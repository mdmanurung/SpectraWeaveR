# Purpose: % TNF+ and % IFNg+ of CD4T and CD8T per arm and sample at the positive cut of the dataset (3D: one
#          fixed cut, arcsinh 1.5, cofactor 6000) on the common SpectroFlo gates (plan common_gates_plan.md:
#          G1a, G7a, G10b). Every arm uses the same cell labels, so differences come from unmixing only.
#          Copy of 38_fixed_cut.R: new gated dir, cut 1.5, new run dir; figures redrawn in the Nature style
#          with donors as points and a median bar. The two cytokine slots and their labels come from the
#          dataset config (TNF and IFNg in 3D and 24H; CD1c and HLA-DR in EXV).
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> \
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/38c_fixed_cut_common.R [<run_dir>]
#          from the repo root, after batch 0 has written GATED_DIR/README.md. <run_dir> overrides the registry
#          output run_out("common-cut15-38c"); it must not exist.
# Inputs:  run_in("common-gated")/manifest.csv, gs/<arm>/  common-gate GatingSets [97_common_gs.R]
#          cut_of(arm, sample, channel, "pos")  positive cut [config pos_rule; cuts.csv from 97e when not fixed]
#          experiments/benchmarks/05_bmv-unmix-compare/datasets/  config, run registry, sample table
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R, report_figures.R  labels, builders
# Outputs: <run_dir>/pct_fixed_source.csv                 % above the cut per arm, sample, parent, cytokine
#          <run_dir>/fig_pct_fixed.{pdf,png,csv}          samples as points, median bar
#          <run_dir>/<tnf>_vs_<ifng>_fixed_<population>_<sample>.{pdf,png,csv}  slot 1 against slot 2 hexbins,
#                                                         cut drawn (tnf_vs_ifng_fixed_... in 3D)
#          <run_dir>/sessionInfo_38c.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), run_out(), cut_of()
source("experiments/benchmarks/05_bmv-unmix-compare/report_figures.R") # theme, save_panel, builders, ARM_LABELS

# Parameters ------------------------------------------------------------------------------
GATED_DIR <- run_in("common-gated") # batch 0 output (G9: read-only here)
# These T-cell nodes score TNF and IFNg; a dataset that labels the slots otherwise has no such target
if (
  !identical(
    c(cfg("targets")$tnf_label, cfg("targets")$ifng_label),
    c("TNF", "IFNg")
  )
) {
  stop(
    "38c_fixed_cut_common.R: needs TNF and IFNg, missing in dataset ",
    DS$dataset,
    " (slots are labelled ",
    cfg("targets")$tnf_label,
    " and ",
    cfg("targets")$ifng_label,
    ")",
    call. = FALSE
  )
}
# Cytokine slots, as in 19a_gate_arm.R; stored as arcsinh. Names are the dataset's marker labels.
TNF_LAB <- cfg("targets")$tnf_label
IFNG_LAB <- cfg("targets")$ifng_label
CHANNELS <- set_names(
  c(cfg("channels")$TNF_slot, cfg("channels")$IFNg_slot),
  c(TNF_LAB, IFNG_LAB)
)
# Cut text for titles: the fixed value when the rule is fixed (3D: 1.5), else the per-sample cut
CUT_IS_FIXED <- identical(cfg("pos_rule")$type, "fixed")
CUT_VALUE_TXT <- if (CUT_IS_FIXED) cfg("pos_rule")$value else "per sample"
CUT_TXT <- if (CUT_IS_FIXED) {
  paste("arcsinh", CUT_VALUE_TXT)
} else {
  "the per-sample cut"
}
NODES <- c(
  CD4T = "/nonDebris/singlets/live/CD3pos/CD4T",
  CD8T = "/nonDebris/singlets/live/CD3pos/CD8T"
)
BINS <- 60 # hexagon bins per axis; display only
AXIS_LIMITS <- c(-0.5, 5) # arcsinh display range (user, 2026-11-09: trims are plots only)
ARM_ORDER <- unlist(VERSION_ROWS, use.names = FALSE) # the 11 figure arms, in figure order

run_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(run_dir)) {
  run_dir <- run_out("common-cut15-38c")
}
stopifnot("38c run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)
manifest <- read_csv(
  file.path(GATED_DIR, "manifest.csv"),
  show_col_types = FALSE
)
stopifnot(
  "all 11 figure arms are in the gated manifest" = all(
    ARM_ORDER %in% manifest$arm
  )
)

# Events of each arm, donor and population: the two cytokines --------------------------------
events <- ARM_ORDER |>
  set_names() |>
  map(\(a) {
    gs <- load_gs(file.path(GATED_DIR, "gs", a))
    sampleNames(gs) |>
      set_names() |>
      map(\(s) {
        NODES |>
          map(\(node) {
            x <- exprs(gh_pop_get_data(gs[[s]], node))[, CHANNELS, drop = FALSE]
            colnames(x) <- names(CHANNELS)
            as_tibble(x)
          }) |>
          list_rbind(names_to = "parent")
      }) |>
      list_rbind(names_to = "donor")
  }) |>
  list_rbind(names_to = "arm")

# Positive cut per arm, sample and cytokine, joined onto the events ----------------------------
cut_tbl <- expand_grid(
  arm = ARM_ORDER,
  donor = samples_tbl()$sample,
  cytokine = names(CHANNELS)
) |>
  mutate(
    cut = pmap_dbl(list(arm, donor, cytokine), \(a, s, k) {
      cut_of(a, s, CHANNELS[[k]], "pos")
    })
  )
events <- events |>
  left_join(
    pivot_wider(
      cut_tbl,
      names_from = cytokine,
      values_from = cut,
      names_prefix = "cut_"
    ),
    by = join_by(arm, donor),
    relationship = "many-to-one",
    unmatched = "error"
  )

# % above the cut per arm, sample, parent and cytokine ------------------------------------------
pct <- names(CHANNELS) |>
  set_names() |>
  map(\(k) {
    events |>
      summarise(
        n_parent = n(),
        pct_pos = 100 * mean(.data[[k]] >= .data[[paste0("cut_", k)]]),
        .by = c(arm, donor, parent)
      )
  }) |>
  list_rbind(names_to = "cytokine") |>
  mutate(
    arm = factor(arm, levels = ARM_ORDER),
    cytokine = factor(cytokine, levels = names(CHANNELS))
  ) |>
  arrange(arm, donor, parent, cytokine) |>
  mutate(arm = as.character(arm), cytokine = as.character(cytokine)) |>
  select(arm, donor, parent, cytokine, n_parent, pct_pos)
stopifnot(
  "gated samples equal the sample table" = setequal(
    events$donor,
    samples_tbl()$sample
  ),
  "one value per arm, sample, parent, cytokine" = nrow(pct) ==
    length(ARM_ORDER) * nrow(samples_tbl()) * length(NODES) * length(CHANNELS)
)
# Identity check: my counts above the cut equal the saved TNFpos / IFNgpos gates (fixed 1.5, README)
node_counts <- ARM_ORDER |>
  set_names() |>
  map(\(a) {
    gs_pop_get_count_fast(load_gs(file.path(GATED_DIR, "gs", a))) |>
      as_tibble() |>
      dplyr::filter(str_detect(Population, "/(CD4T|CD8T)/(TNF|IFNg)pos$"))
  }) |>
  list_rbind(names_to = "arm") |>
  transmute(
    arm,
    donor = name,
    parent = str_extract(Population, "CD[48]T"),
    cytokine = unname(c(TNF = TNF_LAB, IFNg = IFNG_LAB)[str_extract(
      Population,
      "(TNF|IFNg)(?=pos$)"
    )]),
    n_node = Count
  )
stopifnot(
  "counts above the cut equal the TNFpos / IFNgpos node counts" = nrow(
    inner_join(
      mutate(pct, n_above = round(pct_pos / 100 * n_parent)),
      node_counts,
      by = join_by(arm, donor, parent, cytokine),
      relationship = "one-to-one",
      unmatched = "error"
    ) |>
      dplyr::filter(n_above != n_node)
  ) ==
    0
)
write_csv(pct, file.path(run_dir, "pct_fixed_source.csv"))

# Figure 1: % positive per arm, donors as points over a median bar ------------------------------
redraw_panel(
  file.path(run_dir, "pct_fixed_source.csv"),
  "pct_pos",
  paste0("% positive at ", CUT_TXT, ", common SpectroFlo gates"),
  run_dir,
  "fig_pct_fixed"
)

# Figure 2: TNF vs IFNg hexbins, one file per donor and parent, cut drawn ------------------------
pct_label <- pct |>
  pivot_wider(
    names_from = cytokine,
    values_from = pct_pos,
    names_prefix = "pct_"
  ) |>
  mutate(
    panel = sprintf(
      "%s\n%s+ %.1f%%  %s+ %.1f%%",
      ARM_LABELS[arm],
      TNF_LAB,
      .data[[paste0("pct_", TNF_LAB)]],
      IFNG_LAB,
      .data[[paste0("pct_", IFNG_LAB)]]
    ),
    panel = factor(
      as.character(panel),
      levels = unique(as.character(panel))
    )
  )
for (s in sort(unique(events$donor))) {
  for (pop in names(NODES)) {
    plot_df <- events |>
      dplyr::filter(donor == s, parent == pop) |>
      left_join(
        dplyr::filter(pct_label, donor == s, parent == pop),
        by = join_by(arm, donor, parent),
        relationship = "many-to-one"
      ) |>
      mutate(
        panel = factor(
          as.character(panel),
          levels = levels(pct_label$panel)
        ) |>
          droplevels()
      )
    cut_df <- plot_df |>
      distinct(
        panel,
        cut_x = .data[[paste0("cut_", TNF_LAB)]],
        cut_y = .data[[paste0("cut_", IFNG_LAB)]]
      )
    p <- ggplot(plot_df, aes(.data[[TNF_LAB]], .data[[IFNG_LAB]])) +
      geom_hex(bins = BINS) +
      geom_vline(
        data = cut_df,
        aes(xintercept = cut_x),
        colour = "red",
        linewidth = 0.3
      ) +
      geom_hline(
        data = cut_df,
        aes(yintercept = cut_y),
        colour = "red",
        linewidth = 0.3
      ) +
      scale_fill_viridis_c(
        option = "magma",
        transform = "log10",
        name = "Cells"
      ) +
      coord_cartesian(xlim = AXIS_LIMITS, ylim = AXIS_LIMITS) +
      facet_wrap(vars(panel), ncol = 4, axes = "all", axis.labels = "margins") +
      labs(
        x = paste0(TNF_LAB, " (arcsinh)"),
        y = paste0(IFNG_LAB, " (arcsinh)"),
        title = paste(pop, s, "- cut", CUT_VALUE_TXT)
      )
    hex_counts <- layer_data(p, 1) |>
      transmute(panel = levels(plot_df$panel)[PANEL], x, y, cells = count)
    save_panel(
      p,
      hex_counts,
      run_dir,
      paste0(
        tolower(TNF_LAB),
        "_vs_",
        tolower(IFNG_LAB),
        "_fixed_",
        pop,
        "_",
        s
      ),
      183,
      130
    )
  }
}
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_38c.txt")
)
