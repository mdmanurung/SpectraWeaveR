# Purpose: hexbin of the two cytokine slots (y; TNF and IFNg in 3D and 24H, CD1c and HLA-DR in EXV) against
#          SSC-A (x) on CD4T and CD8T with the positive cut (3D: fixed arcsinh 1.5, G7a) drawn, one figure per
#          sample, population and cytokine with one panel per arm; plus % positive per arm per population.
#          Same cut and common SpectroFlo gates as 38c_fixed_cut_common.R.
#          Copy of 40_ssc_fixed_cut.R: new gated dir, cut 1.5, new run dir, Nature-style figures.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> \
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/40c_ssc_fixed_cut_common.R [<run_38c> [<run_dir>]]
#          from the repo root, after 38c_fixed_cut_common.R. <run_38c> defaults to run_in("common-cut15-38c")
#          and <run_dir> to run_out("common-cut15-40c"); <run_dir> must not exist.
# Inputs:  run_in("common-gated")/manifest.csv, gs/<arm>/  common-gate GatingSets [97_common_gs.R]
#          <run_38c>/pct_fixed_source.csv  % above the cut [38c_fixed_cut_common.R; registry key
#            common-cut15-38c is the output key of 38c]
#          cut_of(arm, sample, channel, "pos")  positive cut [config pos_rule; cuts.csv from 97e when not fixed]
#          experiments/benchmarks/05_bmv-unmix-compare/datasets/  config, run registry, sample table
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R, report_figures.R  labels, builders
# Outputs: <run_dir>/ssc_vs_<cytokine>_<population>_<sample>.{pdf,png,csv}  hexbins, cut drawn
#          <run_dir>/pct_pos_<population>.{pdf,png,csv}  % positive, donors as points, median bar
#          <run_dir>/sessionInfo_40c.txt

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
    "40c_ssc_fixed_cut_common.R: needs TNF and IFNg, missing in dataset ",
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
CHANNELS <- set_names(
  c(cfg("channels")$TNF_slot, cfg("channels")$IFNg_slot),
  c(cfg("targets")$tnf_label, cfg("targets")$ifng_label)
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
Y_LIMITS <- c(-0.5, 6) # arcsinh display range (user, 2026-11-09: trims are plots only)
ARM_ORDER <- unlist(VERSION_ROWS, use.names = FALSE)

args <- commandArgs(trailingOnly = TRUE)
run_38c <- if (length(args) >= 1) args[[1]] else run_in("common-cut15-38c")
run_dir <- if (length(args) >= 2) args[[2]] else run_out("common-cut15-40c")
stopifnot("40c run dir must not exist" = !dir.exists(run_dir))
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
pct <- read_csv(
  file.path(run_38c, "pct_fixed_source.csv"),
  show_col_types = FALSE
)

# Events of each arm, donor and population: SSC-A and the two cytokines -----------------------
events <- ARM_ORDER |>
  set_names() |>
  map(\(a) {
    gs <- load_gs(file.path(GATED_DIR, "gs", a))
    sampleNames(gs) |>
      set_names() |>
      map(\(s) {
        NODES |>
          map(\(node) {
            x <- exprs(gh_pop_get_data(gs[[s]], node))[,
              c("SSC-A", CHANNELS),
              drop = FALSE
            ]
            colnames(x) <- c("SSC", names(CHANNELS))
            as_tibble(x)
          }) |>
          list_rbind(names_to = "parent")
      }) |>
      list_rbind(names_to = "donor")
  }) |>
  list_rbind(names_to = "arm")
# Positive cut per arm, sample and cytokine, joined onto the events
events <- events |>
  left_join(
    expand_grid(
      arm = ARM_ORDER,
      donor = samples_tbl()$sample,
      cytokine = names(CHANNELS)
    ) |>
      mutate(
        cut = pmap_dbl(list(arm, donor, cytokine), \(a, s, k) {
          cut_of(a, s, CHANNELS[[k]], "pos")
        })
      ) |>
      pivot_wider(
        names_from = cytokine,
        values_from = cut,
        names_prefix = "cut_"
      ),
    by = join_by(arm, donor),
    relationship = "many-to-one",
    unmatched = "error"
  )
# the % in each panel title must equal 38c's, recomputed here from the same events
check <- names(CHANNELS) |>
  set_names() |>
  map(\(k) {
    events |>
      summarise(
        pct_now = 100 * mean(.data[[k]] >= .data[[paste0("cut_", k)]]),
        .by = c(arm, donor, parent)
      )
  }) |>
  list_rbind(names_to = "cytokine") |>
  inner_join(
    pct,
    by = join_by(arm, donor, parent, cytokine),
    relationship = "one-to-one"
  )
stopifnot(
  "40c % equals 38c %" = nrow(check) == nrow(pct) &&
    all(abs(check$pct_now - check$pct_pos) < 1e-9)
)

# SSC-A vs cytokine hexbins, cut drawn ------------------------------------------------------------
for (s in sort(unique(events$donor))) {
  for (pop in names(NODES)) {
    for (cy in names(CHANNELS)) {
      labels <- pct |>
        dplyr::filter(donor == s, parent == pop, cytokine == cy) |>
        mutate(
          panel = sprintf("%s\n%s+ %.1f%%", ARM_LABELS[arm], cy, pct_pos)
        ) |>
        (\(d) set_names(d$panel, d$arm))()
      plot_df <- events |>
        dplyr::filter(donor == s, parent == pop) |>
        mutate(
          panel = factor(
            unname(labels[arm]),
            levels = unname(labels[ARM_ORDER])
          )
        )
      cut_df <- plot_df |>
        distinct(panel, cut_y = .data[[paste0("cut_", cy)]])
      p <- ggplot(plot_df, aes(SSC, .data[[cy]])) +
        geom_hex(bins = BINS) +
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
        scale_x_continuous(
          labels = scales::label_number(scale = 1e-6, suffix = "M")
        ) +
        coord_cartesian(ylim = Y_LIMITS) +
        facet_wrap(
          vars(panel),
          ncol = 4,
          axes = "all",
          axis.labels = "margins"
        ) +
        labs(
          x = "SSC-A",
          y = paste0(cy, " (arcsinh)"),
          title = paste(pop, s, "- cut", CUT_VALUE_TXT)
        )
      hex_counts <- layer_data(p, 1) |>
        transmute(panel = levels(plot_df$panel)[PANEL], x, y, cells = count)
      save_panel(
        p,
        hex_counts,
        run_dir,
        paste0("ssc_vs_", cy, "_", pop, "_", s),
        183,
        130
      )
    }
  }
}

# % positive per arm, one figure per population, donors as points over a median bar -----------
for (pop in names(NODES)) {
  pop_csv <- file.path(run_dir, paste0("pct_pos_", pop, "_source.csv"))
  write_csv(dplyr::filter(pct, parent == pop), pop_csv)
  redraw_panel(
    pop_csv,
    "pct_pos",
    paste0(pop, ", cut at ", CUT_TXT, ", common SpectroFlo gates"),
    run_dir,
    paste0("pct_pos_", pop)
  )
}
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "sessionInfo_40c.txt")
)
