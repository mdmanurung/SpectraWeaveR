# Purpose: redraw the report's table-driven figures (s00 to s15) from the common-gate outputs (SpectroFlo labels,
#          clean3, cut arcsinh 1.5; plans/common_gates_plan.md G7a, G13, G15). Same builders as report_figures.R
#          (donor colour and shape, median bar behind the points, no lines across unordered arms); only the sources,
#          titles and the cut wording change. Figures whose source does not read a gate (controls, spectra,
#          simulation) keep their original table and are marked "unchanged source" in the manifest.
#          s10_spike is a spike figure and is left to batch 5.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/95c_restyle_report_tables.R   (repo root)
# Inputs:  see FIGURES below; every path is under experiments/runs/
# Outputs: experiments/runs/2026-11-11_bmv-restyle-common/report_figures/<key>.{pdf,png,csv}
#          experiments/runs/2026-11-11_bmv-restyle-common/report_figures/captions.csv
#          experiments/runs/2026-11-11_bmv-restyle-common/report_figures/af_swap_gap_closed_source.csv  derived table
#          experiments/runs/2026-11-11_bmv-restyle-common/manifest_parts/95c_restyle_report_tables.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")

# Parameters ---------------------------------------------------------------------------------------
SCRIPT <- "95c_restyle_report_tables.R"
OUT <- fresh_dir("report_figures")
OLD_FIG <- "2026-11-10_bmv-report/figures-v7" # where the report cited these figures before
SUMMARY_DIR <- "experiments/runs/2026-11-10_bmv-summary"
READOUT_POS <- paste("CD8T TNF above", CUT_TXT) # readout label in the common gap tables

# Gap panel with the 1.5 wording -------------------------------------------------------------------
# panel_gap() in report_figures.R keys the positive-share readout to "above arcsinh 1"; the common
# tables say "above arcsinh 1.5", so the label map is redefined here and BUILDERS$gap picks it up.
panel_gap <- function(src, title) {
  readout_labels <- c(
    squeezed_share = "CD8T TNF squeezed near 0",
    pos_share = READOUT_POS
  )
  gap <- read_csv(src, show_col_types = FALSE) |>
    mutate(
      readout = if_else(
        .data$readout %in% readout_labels,
        names(readout_labels)[match(.data$readout, readout_labels)],
        .data$readout
      )
    ) |>
    filter(.data$readout %in% names(readout_labels))
  stopifnot("both gap readouts present" = n_distinct(gap$readout) == 2)
  arm_order <- gap |>
    filter(.data$readout == "squeezed_share") |>
    summarise(med = median(.data$frac_closed), .by = arm) |>
    arrange(.data$med) |>
    pull(arm)
  plot_df <- gap |>
    mutate(
      donor = factor_by_pal(.data$donor, DONOR_PAL),
      arm = factor(.data$arm, levels = arm_order),
      readout = factor(readout_labels[.data$readout], levels = readout_labels)
    )
  p <- ggplot(plot_df, aes(frac_closed, arm)) +
    value_grid("x") +
    geom_vline(xintercept = c(0, 1), colour = "grey50") +
    geom_vline(xintercept = GAP_LINES, linetype = "dashed", colour = "grey60") +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.5,
      linewidth = 0.2
    ) +
    geom_point(aes(colour = donor, shape = donor), size = 1.2) +
    donor_scales() +
    facet_wrap(~readout, axes = "all", axis.labels = "margins") +
    labs(
      x = "Fraction of the as188_joint to s102_c188 gap closed",
      y = NULL,
      title = title
    )
  list(
    plot = p,
    data = plot_df,
    width_mm = 183,
    height_mm = 30 + ROW_MM * length(arm_order)
  )
}
BUILDERS$gap <- \(d, src, title) {
  out <- panel_gap(src, title)
  out$caption <- paste(
    "Fraction of the as188_joint to s102_c188 gap each arm closes, for the CD8T TNF squeezed and positive",
    paste0("(cut ", CUT_TXT, ") shares."),
    "Each point is one donor; the bar is the median; dashed lines mark the D3 criterion (every donor 0.5, median 0.8)."
  )
  out
}

# Visual-pass overrides (2026-11-11) -----------------------------------------------------------------
# Builders from report_figures.R, changed here so the library file stays as certified.
# Plain-text titles and x labels of the four section 03 ranking figures.
RANK_TEXT <- list(
  fig_streak_sym = list(
    title = "Streak: TNF–CD8 Spearman rho",
    x = "TNF-CD8 Spearman rho (BV510-A vs SBUV445-A)"
  ),
  fig_ssm_sym = list(
    title = "Spillover spread: quadrature SS",
    x = "Median quadrature SS over source i, target j"
  ),
  fig_tilt_sym = list(
    title = "Tilt: median |rho|",
    x = "Median |Spearman rho| over channel pairs"
  ),
  fig_spread_sym = list(
    title = "Spread: median Nguyen rSD",
    x = "Median Nguyen rSD (q84 − q50)"
  )
)
BUILDERS$rank <- \(d, src, title) {
  d <- rename(d, donor = sample)
  txt <- RANK_TEXT[[title]]
  out <- arm_donor(d, "value", str_wrap(txt$x, 45), txt$title)
  out$caption <- paste("Ranking metric per arm.", PAIRED)
  out
}

# Score tables: strips named for what they show, no x title "Value" (strips carry the quantity).
# pos_share is a fraction in the table; the strip says "% positive", so it is shown in percent.
READOUT_STRIPS <- c(
  pos_share = "% positive",
  rsd_ratio = "neg. rSD ratio",
  squeezed_share = "squeezed share"
)
BUILDERS$scores <- \(d, src, title) {
  long <- d |>
    rename(donor = sample) |>
    mutate(
      panel = paste(
        .data$pop,
        c(`BV510-A` = "TNF", `BV570-A` = "IFNg")[.data$channel]
      ),
      pos_share = 100 * .data$pos_share
    ) |>
    pivot_longer(
      c(squeezed_share, pos_share, rsd_ratio),
      names_to = "readout",
      values_to = "value"
    ) |>
    mutate(readout = factor(READOUT_STRIPS[.data$readout], READOUT_STRIPS))
  out <- arm_donor(long, "value", "", title, rows = "readout", cols = "panel")
  out$plot <- out$plot + labs(x = NULL)
  out$width_mm <- 183
  out$caption <- paste(
    "TNF and IFNg readouts (strips: % positive, neg. rSD ratio, squeezed share) per arm, population and cytokine (rows of panels).",
    PAIRED
  )
  out
}

BUILDERS$interventions <- \(d, src, title) {
  d <- mutate(d, arm = .data$arm_label)
  out <- arm_donor(
    d,
    "value",
    "Share of population cells",
    title,
    rows = "readout",
    cols = "pop"
  )
  out$caption <- paste(
    "Share of cells squeezed near 0 and above the cut, per intervention arm (rows) and population (columns).",
    PAIRED
  )
  out
}

# Passes grid: cell text to 3 significant digits, so small differences show; one column per readout,
# so the grid has no empty slot.
BUILDERS$passes <- \(d, src, title) {
  p <- ggplot(
    d,
    aes(factor(n_passes), factor(n_af_passes), fill = fill_scaled)
  ) +
    geom_tile(colour = "white") +
    geom_text(
      aes(label = signif(median, 3), colour = fill_scaled > 0.5),
      size = 5,
      size.unit = "pt",
      show.legend = FALSE
    ) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
    scale_fill_viridis_c(name = "Scaled median") +
    facet_wrap(
      ~readout,
      nrow = 1,
      labeller = labeller(readout = label_wrap_gen(14)),
      axes = "all",
      axis.labels = "margins"
    ) +
    labs(x = "Unmixing passes", y = "AF passes", title = title)
  list(
    plot = p,
    data = d,
    width_mm = 183,
    height_mm = 55,
    caption = "Median over 5 donors per pass grid cell (text, 3 significant digits); the fill is that median scaled within each readout."
  )
}

# Offsets: a small vertical dodge per set so sets with the same offset (chain4_spec, chain5_bead at
# Alexa Fluor 405) both show.
BUILDERS$offsets <- \(d, src, title) {
  dyes <- rev(unique(d$dye))
  sets <- unique(d$set)
  step <- 0.16
  plot_df <- d |>
    mutate(
      y = match(.data$dye, dyes) +
        (match(.data$set, sets) - (length(sets) + 1) / 2) * step
    )
  p <- ggplot(plot_df, aes(offset, y, colour = set)) +
    value_grid("x") +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(size = 1.2) +
    scale_colour_manual(values = categorical_palette(d$set), name = "Set") +
    scale_y_continuous(breaks = seq_along(dyes), labels = dyes) +
    labs(x = "Offset", y = NULL, title = title)
  list(
    plot = p,
    data = plot_df,
    width_mm = 89,
    height_mm = 30 + ROW_MM * n_distinct(d$dye),
    caption = "Offset per dye for each control set; sets are dodged vertically by 0.16 of a row so equal offsets stay visible."
  )
}

# Screens: tile cells with no score are filled grey80 ("not tested"). Rows that were never written
# (marker not scored in that population) are added as NA so the grey shows. PE has no marker in the
# control file (marker NA in fcs_control_file.csv), so its row reads "PE (no marker)".
fill_untested <- function(d, label_col, x_col) {
  d |>
    mutate(
      "{label_col}" := if_else(
        .data[[label_col]] == "NA (PE)",
        "PE (no marker)",
        .data[[label_col]]
      )
    ) |>
    mutate(
      "{label_col}" := factor(.data[[label_col]], unique(.data[[label_col]])),
      "{x_col}" := factor(.data[[x_col]], unique(.data[[x_col]])),
      pop = factor(.data$pop, unique(.data$pop))
    ) |>
    complete(.data[[x_col]], .data[[label_col]], .data$pop) |>
    mutate(across(where(is.factor), as.character))
}
screen_old <- BUILDERS$screen
BUILDERS$screen <- \(d, src, title) {
  out <- screen_old(fill_untested(d, "label", "arm"), src, title)
  out$plot <- out$plot + scale_x_discrete(labels = arm_label)
  out$caption <- paste(
    "Change per marker (rows) and arm (columns) against the reference, median of 5 donors, per population.",
    "Blue is lower, red higher; grey is not tested. PE has no marker in the control file."
  )
  out
}
map_old <- BUILDERS$dye_channel_map
BUILDERS$dye_channel_map <- \(d, src, title) {
  out <- map_old(
    fill_untested(d, "channel_label", "swapped_dye"),
    src,
    title
  )
  out$caption <- paste(
    "Median change in squeezed share per channel (rows) when one dye (columns) takes 1.0.2 variants,",
    "for CD8T and live cells (the table has no CD4T). Grey is not tested. PE has no marker in the control file."
  )
  out
}

# Derived table: AF set swap gap closed ----------------------------------------------------------------
# The AF-swap gap (as102_var and swap_af102) has no run of its own on the common gates. Its arms are scored in
# common-ablation (as102_var) and swap-test-common (as188_joint, s102_c188, swap_af102), so the gap-closed
# table is rebuilt from those saved scores with the definition the old table used:
# (value - as188_joint) / (s102_c188 - as188_joint) for CD8T BV510-A squeezed_share and pos_share.
gap_from_scores <- function(scores_csv, which_arms) {
  scores <- read_csv(scores_csv, show_col_types = FALSE) |>
    filter(.data$pop == "CD8T", .data$channel == "BV510-A") |>
    select(donor = sample, arm, squeezed_share, pos_share) |>
    pivot_longer(
      c(squeezed_share, pos_share),
      names_to = "readout",
      values_to = "value"
    )
  ref <- scores |>
    filter(.data$arm %in% c("as188_joint", "s102_c188")) |>
    pivot_wider(names_from = arm, values_from = value) |>
    rename(start = as188_joint, target = s102_c188)
  scores |>
    filter(.data$arm %in% which_arms) |>
    left_join(
      ref,
      by = join_by(donor, readout),
      relationship = "many-to-one",
      unmatched = "error"
    ) |>
    mutate(frac_closed = (value - start) / (target - start))
}
# The rebuild must reproduce the original table from the original scores (one check).
old_check <- gap_from_scores(
  file.path(RUNS, "2026-11-10_bmv-swap-test/scores.csv"),
  "swap_af102"
) |>
  left_join(
    read_csv(
      file.path(RUNS, "2026-11-10_bmv-af-swap/fig_gap_closed_source.csv"),
      show_col_types = FALSE
    ) |>
      filter(.data$arm == "swap_af102") |>
      mutate(
        readout = if_else(
          str_detect(readout, "squeezed"),
          "squeezed_share",
          "pos_share"
        )
      ) |>
      select(donor, arm, readout, old = frac_closed),
    by = join_by(donor, arm, readout),
    relationship = "one-to-one",
    unmatched = "error"
  )
stopifnot(
  "derivation reproduces the old af-swap table" = max(abs(
    old_check$frac_closed - old_check$old
  )) <
    1e-9
)
readout_key <- c(
  squeezed_share = "CD8T TNF squeezed near 0",
  pos_share = READOUT_POS
)
af_swap_gap <- bind_rows(
  gap_from_scores(
    file.path(RUNS, "2026-11-11_bmv-swap-test-common/scores.csv"),
    "swap_af102"
  ),
  read_csv(
    file.path(RUNS, "2026-11-11_bmv-common-ablation/fig_gap_closed.csv"),
    show_col_types = FALSE
  ) |>
    filter(.data$arm == "as102_var")
) |>
  mutate(
    readout = if_else(
      .data$readout %in% names(readout_key),
      readout_key[.data$readout],
      .data$readout
    )
  ) |>
  select(donor, arm, readout, value, start, target, frac_closed)
AF_SWAP_SRC <- file.path(OUT, "af_swap_gap_closed_source.csv")
write_csv(af_swap_gap, AF_SWAP_SRC)

# Figures ---------------------------------------------------------------------------------------------
# key, source (under experiments/runs/, or the derived table), layout, title, source status
C <- "2026-11-11_bmv-" # prefix of the common-gate runs
C34 <- paste0(C, "common-cut15-34e") # G16 rerun of the 34c streak, tilt, spread and SSM rows
FIGURES <- tribble(
  ~key                 , ~src                                                                      , ~layout           , ~title                               , ~status            ,
  "s00_per_dye"        , paste0(C, "variant-per-dye-common/fig_cd8t_tnf_squeezed.csv")             , "dye_share"       , "CD8T TNF squeezed, one dye swapped" , "redrawn"          ,
  "s01_pct_fixed"      , paste0(C, "common-cut15-38c/pct_fixed_source.csv")                        , "pct_pos"         , paste("% positive at", CUT_TXT)      , "redrawn"          ,
  "s03_streak"         , paste0(C34, "/fig_streak_sym_source.csv")                                 , "rank"            , "fig_streak_sym"                     , "redrawn"          ,
  "s03_ssm"            , paste0(C34, "/fig_ssm_sym_source.csv")                                    , "rank"            , "fig_ssm_sym"                        , "redrawn"          ,
  "s03_tilt"           , paste0(C34, "/fig_tilt_sym_source.csv")                                   , "rank"            , "fig_tilt_sym"                       , "redrawn"          ,
  "s03_spread"         , paste0(C34, "/fig_spread_sym_source.csv")                                 , "rank"            , "fig_spread_sym"                     , "redrawn"          ,
  "s04_g1"             , "2026-11-10_bmv-g1g2-figures/fig_g1_heldout_source.csv"                   , "g1"              , "G1 held-out controls"               , "unchanged source" ,
  "s04_g2"             , "2026-11-10_bmv-g1g2-figures/fig_g2_simulation_source.csv"                , "g2"              , "G2 simulation"                      , "unchanged source" ,
  "s05_swap"           , paste0(C, "swap-test-common/scores.csv")                                  , "scores"          , "Swap test"                          , "redrawn"          ,
  "s05_count"          , paste0(C, "variant-count-common/scores.csv")                              , "scores"          , "Variant count"                      , "redrawn"          ,
  "s05_af_swap"        , AF_SWAP_SRC                                                               , "gap"             , "AF set swap"                        , "redrawn, derived" ,
  "s05_bv510_by_set"   , "2026-11-10_bmv-variant-bv510-by-set/fig_variant_bv510_by_set_source.csv" , "variant_set"     , "BV510 part by set"                  , "unchanged source" ,
  "s05_coherence_or"   , paste0(C, "common-coherence-or/coherence_or.csv")                         , "coherence"       , "TNF-IFNg odds ratio"                , "redrawn"          ,
  "s07_neg_spread"     , paste0(C, "common-cut15-47c/fig_f_neg_spread.csv")                        , "neg_spread"      , "Negative spread"                    , "redrawn"          ,
  "s07_pct_pos"        , paste0(C, "common-cut15-47c/fig_a_pct_pos.csv")                           , "pct_pos"         , "IL-2, IL-10 and Th2"                , "redrawn"          ,
  "s09_block_a"        , paste0(C, "common-ablation/fig_gap_closed.csv")                           , "gap"             , "Block A swaps"                      , "redrawn"          ,
  "s09_chain"          , paste0(C, "common-chain-all/fig_gap_closed.csv")                          , "gap"             , "Rebuild chain"                      , "redrawn"          ,
  "s09_chain7a8"       , paste0(C, "common-chain7a8/fig_gap_closed.csv")                           , "gap"             , "Chain 7a and 8"                     , "redrawn"          ,
  "s10_interventions"  , paste0(C, "common-af405-proof/fig2_interventions.csv")                    , "interventions"   , "AF405 interventions"                , "redrawn"          ,
  "s10_dose"           , paste0(C, "common-af405-proof/fig4_dose.csv")                             , "dose"            , "AF405 spread dose"                  , "redrawn"          ,
  "s11_gls"            , paste0(C, "gls-noise-common/scores/scores.csv")                           , "scores"          , "GLS noise model"                    , "redrawn"          ,
  "s12_passes"         , paste0(C, "common-passes-grid/fig_passes_grid.csv")                       , "passes"          , "Passes grid"                        , "redrawn"          ,
  "s13_opt188"         , paste0(C, "common-opt188/report/opt188_readouts.csv")                     , "readouts"        , "opt188 readouts"                    , "redrawn"          ,
  "s15_102code_spread" , "2026-11-10_bmv-af405-102code-clean/fig_af405_spread_source.csv"          , "spread_coef"     , "AF405 spread, clean controls"       , "unchanged source" ,
  "s15_102code_gap"    , paste0(C, "common-af405-102code-clean-score/fig_gap_closed.csv")          , "gap"             , "1.0.2-code AF405 variants"          , "redrawn"          ,
  "s15_origin_spread"  , "2026-11-10_bmv-af405-af-origin/fig_af_origin_spread_source.csv"          , "spread_coef"     , "AF405 spread by AF half"            , "unchanged source" ,
  "s15_origin_shape"   , "2026-11-10_bmv-af405-af-origin/fig_af_origin_shape_source.csv"           , "af_shape"        , "Spread shape"                       , "unchanged source" ,
  "s15_offsets"        , "2026-11-10_bmv-automated-raw-controls/fig_offsets_source.csv"            , "offsets"         , "Reference offsets"                  , "unchanged source" ,
  "s15_sbuv510"        , paste0(C, "common-sbuv510-decompose-score/fig_gap_closed.csv")            , "gap"             , "SBUV510 decomposition"              , "redrawn"          ,
  "s15_ssm_zip"        , "2026-11-10_bmv-own-ssm/fig_ssm_vs_zip_source.csv"                        , "ssm_zip"         , "Own SSM vs zip"                     , "unchanged source" ,
  "s15_screen"         , paste0(C, "all-marker-screen-common/fig_screen_squeezed.csv")             , "screen"          , "All-marker screen"                  , "redrawn"          ,
  "s15_channel_map"    , paste0(C, "per-dye-channel-map-common/fig_per_dye_channel_map.csv")       , "dye_channel_map" , "Per-dye channel map"                , "redrawn"
)

# Absolute sources start with "experiments/"; the others sit under RUNS.
src_path <- \(s) if_else(str_starts(s, "experiments/"), s, file.path(RUNS, s))
drawn <- FIGURES |>
  filter(only_key(.data$key)) |>
  mutate(
    res = pmap(
      list(src_path(.data$src), .data$layout, .data$title, .data$key),
      \(src, layout, title, key) redraw_panel(src, layout, title, OUT, key)
    )
  ) |>
  unnest(res)

# Section 06: angle between the 1.0.2 and 1.8.8 dye spectra; spectra do not read a gate (unchanged source).
# The AF row has no 1.0.2 spectrum (angle NA); it is dropped and the caption says so.
angle_src <- file.path(SUMMARY_DIR, "06_spectra_variants/spectra_compare.csv")
angles <- read_csv(angle_src, show_col_types = FALSE) |>
  filter(!is.na(.data$angle_deg)) |>
  mutate(fluorophore = reorder(.data$fluorophore, .data$angle_deg))
p_angle <- ggplot(angles, aes(angle_deg, fluorophore)) +
  geom_point(size = 1.2) +
  labs(
    x = "Angle, 1.0.2 vs 1.8.8 (degrees)",
    y = NULL,
    title = "Dye spectra, 1.0.2 vs 1.8.8"
  ) +
  value_grid("x")
if (only_key("s06_spectra_angle")) {
  save_panel(p_angle, angles, OUT, "s06_spectra_angle", 89, 120)
  drawn <- bind_rows(
    drawn,
    tibble(
      key = "s06_spectra_angle",
      src = angle_src,
      status = "unchanged source",
      redrawn = file.path(OUT, "s06_spectra_angle.png"),
      caption = "Angle between each dye's 1.0.2 and 1.8.8 spectrum; each point is one dye. AF: no 1.0.2 spectrum."
    )
  )
}

# Captions and manifest -----------------------------------------------------------------------------
# A partial redraw replaces only the rows of the figures it drew.
cap_file <- file.path(OUT, "captions.csv")
caps <- select(drawn, key, caption)
if (!is.null(ONLY)) {
  caps <- bind_rows(
    filter(read_csv(cap_file, show_col_types = FALSE), !key %in% caps$key),
    caps
  )
}
write_csv(caps, cap_file)
pwalk(
  drawn,
  \(key, src, status, redrawn, ...) {
    add_manifest(
      SCRIPT,
      file.path(OLD_FIG, paste0(key, ".png")),
      redrawn,
      src,
      status = status
    )
  }
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT, "sessionInfo_95c_report_tables.txt")
)
