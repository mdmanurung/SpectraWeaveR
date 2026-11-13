# Purpose: shared figure builders for the BMV report notebooks (90_report.qmd, 91_all_figures.qmd), in
#          the lucid-coder figure style (user, 2026-11-10). Each builder reads a SAVED table only and returns
#          list(plot, data, width_mm, height_mm) for save_panel(); nothing is recomputed from events.
# Run:     SW_DATASET=<3d|24h|exv> source("experiments/benchmarks/05_bmv-unmix-compare/report_figures.R") from the
#          repo root (SW_RUN_DATE only where the calling script writes a run folder). The dataset sets the
#          donor list and palette (samples_tbl()), the cytokine labels and channels (cfg("targets"), cfg("channels")).
#          Nature spec (user, 2026-11-11): widths 89 or 183 mm, text 5-7 pt; save_panel() enforces it.
# Inputs:  experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R   DS, cfg(), samples_tbl() (registry keys:
#            none; the saved tables come from the calling script)
#          experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R   Nature theme, value_grid(), save_panel()
#          experiments/benchmarks/05_bmv-unmix-compare/lucid/palettes.R, plots.R
#            [lucid-coder 0.5.0 templates, copied 2026-11-10]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R   ARM_LABELS for the 11 main arms
# Outputs: none; defines DONORS, N_DONORS, DONOR_PAL, DONOR_SHAPES, arm_label(), panel_gap(), panel_arm_donor(),
#          panel_signed_tile(), panel_set_points(), LAYOUTS, BUILDERS and redraw_panel()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lucid/palettes.R")
source("experiments/benchmarks/05_bmv-unmix-compare/lucid/plots.R")
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R")
# Every facet draws its own axis lines and ticks, labels on the outer panels only (r-0008).
FACET_AXES <- list(axes = "all", axis.labels = "margins")
# Height per category row in a flipped panel: one 7 pt label plus spacing.
ROW_MM <- 4.5

# Dataset labels and channels -------------------------------------------------------------------
# The two cytokine dye slots keep their channel; the marker label comes from the dataset config.
TNF_LABEL <- cfg("targets")$tnf_label
IFNG_LABEL <- cfg("targets")$ifng_label
TNF_DYE <- cfg("targets")$tnf_dye
TNF_SLOT <- cfg("channels")$TNF_slot
IFNG_SLOT <- cfg("channels")$IFNg_slot
CD8_CHANNEL <- cfg("channels")$CD8

# Donor encoding -------------------------------------------------------------------------------
# Donors are the replicates: one colour and one shape per donor, the same in every figure. Without
# conditions the table's `donor` column holds the sample id (the well in 3D and EXV); with conditions
# (24H) the samples are donor x condition, so the replicate is the donor id and `condition` facets.
DONORS <- if (is.null(DS$conditions)) {
  samples_tbl()$sample
} else {
  unique(samples_tbl()$donor)
}
N_DONORS <- length(DONORS)
DONOR_PAL <- categorical_palette(DONORS)
# the first five are the 3D shapes; the rest serve datasets with more donors
DONOR_SHAPE_POOL <- c(16, 17, 15, 18, 25, 8, 3, 4, 7, 9)
stopifnot(
  "a point shape for every donor" = N_DONORS <= length(DONOR_SHAPE_POOL)
)
DONOR_SHAPES <- setNames(DONOR_SHAPE_POOL[seq_len(N_DONORS)], DONORS)
# Table with a `sample` column to the replicate columns: `donor` (and `condition` when the dataset has them)
with_donor <- function(d) {
  if (is.null(DS$conditions)) {
    return(rename(d, donor = sample))
  }
  d |>
    left_join(
      select(samples_tbl(), sample, donor, condition),
      by = join_by(sample),
      relationship = "many-to-one",
      unmatched = "error"
    ) |>
    select(-sample)
}
# Gap criterion from the README decision table (D3): median >= 0.8 and every donor >= 0.5.
GAP_LINES <- c(0.5, 0.8)

# Arm IDs keep their table keys; the 11 main arms read as their figure labels.
arm_label <- function(x) {
  if_else(x %in% names(ARM_LABELS), unname(ARM_LABELS[x]), x)
}

donor_scales <- function() {
  list(
    scale_colour_manual(values = DONOR_PAL, name = "Donor"),
    scale_shape_manual(values = DONOR_SHAPES, name = "Donor")
  )
}

# Gap closed per arm -------------------------------------------------------------------------
# A gap_closed table (donor, arm, readout, frac_closed): fraction of the as188_joint -> s102_c188 gap
# each arm closes. Arms sorted by median squeeze; donors as points over a median crossbar.
panel_gap <- function(src, title) {
  readout_labels <- c(
    squeezed_share = sprintf("CD8T %s squeezed near 0", TNF_LABEL),
    pos_share = sprintf("CD8T %s above arcsinh 1", TNF_LABEL)
  )
  # Some runs saved the readout key, others its label; both map to the key.
  gap <- read_csv(src, show_col_types = FALSE) |>
    mutate(
      readout = if_else(
        .data$readout %in% readout_labels,
        names(readout_labels)[match(.data$readout, readout_labels)],
        .data$readout
      )
    ) |>
    filter(.data$readout %in% names(readout_labels))
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

# One value per arm and donor ------------------------------------------------------------------
# Paired design: every arm unmixes the same donors; colour and shape mark the donor. No lines join a
# donor across arms, because the arms have no order (user, 2026-11-11).
# Arms run down the y axis with horizontal labels (r-0012), first arm on top; the value is on x with
# a thin grid (r-0014). The table's `rows` (readouts with their own units) become columns with
# their own x scale; its `cols` (populations) become rows. With conditions, `condition` is one more facet
# row (24H: one row block per condition), and the many-readout wrap is off.
panel_arm_donor <- function(
  df,
  value,
  value_label,
  title,
  arm_levels,
  rows = NULL,
  cols = NULL
) {
  plot_df <- df |>
    mutate(
      donor = factor_by_pal(.data$donor, DONOR_PAL),
      arm = factor(.data$arm, levels = arm_levels)
    )
  stopifnot("every arm has a level" = !anyNA(plot_df$arm))
  by_condition <- !is.null(DS$conditions)
  if (by_condition && !"condition" %in% names(plot_df)) {
    stop(
      "report_figures.R panel_arm_donor: dataset ",
      DS$dataset,
      " has conditions but the table has no condition column",
      call. = FALSE
    )
  }
  p <- ggplot(plot_df, aes(.data[[value]], arm)) +
    value_grid("x") +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.5,
      linewidth = 0.2,
      orientation = "y"
    ) +
    geom_point(aes(colour = donor, shape = donor), size = 1.1) +
    donor_scales() +
    scale_y_discrete(limits = rev, labels = arm_label) +
    labs(x = value_label, y = NULL, title = title) +
    # room between panels so the end tick labels of neighbouring free x axes do not touch
    theme(panel.spacing.x = unit(4, "mm"))
  n_readouts <- if (is.null(rows)) 1 else n_distinct(plot_df[[rows]])
  # many readouts and no population split: wrap them three to a row so strips and ticks fit
  wrap <- !by_condition && is.null(cols) && n_readouts > 4
  if (by_condition) {
    p <- p +
      facet_grid(
        rows = vars(!!!rlang::syms(c(cols, "condition"))),
        cols = if (is.null(rows)) NULL else vars(.data[[rows]]),
        scales = "free_x",
        axes = "all",
        axis.labels = "margins"
      )
  } else if (wrap) {
    p <- p +
      facet_wrap(
        vars(.data[[rows]]),
        ncol = 3,
        scales = "free_x",
        axes = "all",
        axis.labels = "margins"
      )
  } else if (!is.null(rows) || !is.null(cols)) {
    p <- p +
      facet_grid(
        rows = if (is.null(cols)) NULL else vars(.data[[cols]]),
        cols = if (is.null(rows)) NULL else vars(.data[[rows]]),
        scales = "free_x",
        axes = "all",
        axis.labels = "margins"
      )
  }
  n_facet_cols <- if (wrap) 3 else n_readouts
  n_facet_rows <- if (wrap) {
    ceiling(n_readouts / 3)
  } else if (is.null(cols)) {
    1
  } else {
    n_distinct(plot_df[[cols]])
  }
  if (by_condition) {
    n_facet_rows <- n_facet_rows * n_distinct(plot_df$condition)
  }
  list(
    plot = p,
    data = plot_df,
    width_mm = if (n_facet_cols > 1) 183 else 89,
    height_mm = min(240, 30 + n_facet_rows * (6 + ROW_MM * length(arm_levels)))
  )
}

# Signed tile map ------------------------------------------------------------------------------
# Signed differences on blue-white-red centred at 0; `limit` is the symmetric saturation, values beyond
# are squished. Cells with no value stay grey ("not tested"), never 0.
panel_signed_tile <- function(
  df,
  x,
  y,
  fill,
  limit,
  fill_label,
  title,
  cols = NULL,
  label = NULL
) {
  p <- ggplot(df, aes(.data[[x]], .data[[y]], fill = .data[[fill]])) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-limit, limit),
      oob = scales::squish,
      na.value = "grey80",
      name = fill_label
    ) +
    labs(x = NULL, y = NULL, title = title) +
    rotate_x()
  if (!is.null(label)) {
    p <- p + geom_text(aes(label = .data[[label]]), size = 6, size.unit = "pt")
  }
  if (!is.null(cols)) {
    # more than 15 columns per panel: stack the panels so each gets the full width for its labels
    stack <- n_distinct(df[[x]]) > 15
    p <- p +
      facet_wrap(
        vars(.data[[cols]]),
        nrow = if (stack) NULL else 1,
        ncol = if (stack) 1 else NULL,
        axes = "all",
        axis.labels = "margins"
      )
  }
  n_stacked <- if (!is.null(cols) && n_distinct(df[[x]]) > 15) {
    n_distinct(df[[cols]])
  } else {
    1
  }
  list(
    plot = p,
    data = df,
    width_mm = 183,
    height_mm = min(240, 45 + n_stacked * (6 + 2.6 * n_distinct(df[[y]])))
  )
}

# Many points per set (variants, controls) ------------------------------------------------------
# Not donors: each point is a variant or control. A box from 10 points, jittered points (seeded) on
# top. Sets run down the y axis with their n (r-0012); the value is on x with a thin grid (r-0014).
panel_set_points <- function(
  df,
  x,
  value,
  value_label,
  title,
  x_levels,
  cols = NULL,
  seed = 1
) {
  plot_df <- df |>
    mutate(set = factor(.data[[x]], levels = x_levels))
  stopifnot("every set has a level" = !anyNA(plot_df$set))
  # n counts points per set in one panel: every panel holds the same variants or controls
  n_df <- if (is.null(cols)) {
    plot_df
  } else {
    filter(plot_df, .data[[cols]] == first(.data[[cols]]))
  }
  p <- ggplot(plot_df, aes(.data[[value]], set)) +
    value_grid("x") +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_boxplot(outlier.shape = NA, fill = NA, width = 0.55, linewidth = 0.3) +
    geom_point(
      position = position_jitter(width = 0, height = 0.18, seed = seed),
      size = 0.5,
      alpha = 0.5
    ) +
    scale_y_discrete(limits = rev, labels = n_labels(n_df$set)) +
    labs(x = value_label, y = NULL, title = title)
  if (!is.null(cols)) {
    p <- p +
      facet_wrap(
        vars(.data[[cols]]),
        scales = "free_x",
        axes = "all",
        axis.labels = "margins"
      )
  }
  n_panels <- if (is.null(cols)) 1 else n_distinct(plot_df[[cols]])
  n_wrap_rows <- ceiling(n_panels / min(3, n_panels))
  list(
    plot = p,
    data = plot_df,
    width_mm = if (n_panels > 1) 183 else 89,
    height_mm = min(240, 30 + n_wrap_rows * (8 + 7 * length(x_levels)))
  )
}

# Event panels (F3/F4) --------------------------------------------------------------------
# Display transform of the unmix event panels: pseudo_log_sigma = 1000, biplot_bins = 200 and
# biplot_pair BUV395-A x BV421-A in every benchmark config (01, 02, 03); mirrored, not chosen here.
PSEUDO_LOG_SIGMA <- 1000
BIPLOT_BINS <- 200
BIPLOT_PAIR <- c("BUV395-A", "BV421-A")

# Layout per column signature ---------------------------------------------------
# One builder per table layout. Each takes the table and the figure stem and returns
# list(plot, data, width_mm, height_mm, caption).
LAYOUTS <- c(
  "donor,arm,readout,value,start,target,frac_closed" = "gap",
  "arm,donor,value,open,row,col" = "cyto",
  "arm,sample,value" = "rank",
  "arm,sample,value,median_spread_q,n_channels" = "rank",
  "arm,sample,value,n_cells,n_na,n_low_n" = "rank",
  "arm,sample,value,n_pairs,n_na" = "rank",
  "arm,donor,parent,cytokine,pct_pos,failed" = "pct_pos",
  "donor,arm,parent,cytokine,pct_pos" = "pct_pos",
  "arm,donor,parent,cytokine,n_parent,pct_pos" = "pct_pos",
  "donor,arm,parent,cytokine,value,density,version" = "density_curves",
  "set,variant,spread_coef" = "spread_coef",
  "arm,pop,channel,marker,median_d_squeezed,min_d_squeezed,median_rsd_ratio,max_rsd_ratio,median_d_pos,label" = "screen",
  "donor,arm,parent,pair,rho,kind" = "rho_pairs",
  "dye,arm,sample,squeezed_share,median_share" = "dye_share",
  "fluorophore_x,fluorophore_y,cosine" = "cosine",
  "fluorophore_x,fluorophore_y,cosine,cosine_raw" = "cosine",
  "arm,control,group,n_controls,metric,log2_ratio" = "g1",
  "arm,donor,parent,cytokine,marker,rho" = "rho_marker",
  "arm,donor,parent,n_ifng_pos,pct_tnf_of_ifng,pct_tnf_of_ifng_neg" = "tnf_of_ifng",
  "arm,fluor,metric,ratio,ratio_min_seed,ratio_max_seed,mean_value,mean_ref" = "g2",
  "arm,n_passes,n_af_passes,readout,median,min,max,n_donors,fill_scaled" = "passes",
  "cell,detector,series,value" = "anatomy",
  "detector,value,set,line" = "af_shape",
  "donor,arm,arm_label,pop,readout,value" = "interventions",
  "donor,arm,parent,among,n_among,pct_il2" = "il2",
  "donor,arm,pop,cytokine,rsd_neg,ref_rsd,pos_share,squeezed_share,rsd_ratio,arm_label,alpha" = "dose",
  "donor,parent,cytokine,arm,n_shared,rsd,rsd_ratio" = "neg_spread",
  "sample,arm,bv510_asinh" = "rawaf_density",
  "sample,pop,glut1_decile,n,share_newly_squeezed,median_tnf_change,median_glut1,median_af,share_glut1_above_thr" = "glut1",
  "set,dye,offset" = "offsets",
  "set,dye,variant,spread_coef,offset_coef" = "variant_set",
  "setup,pop,level,glut1_decile,sample,median_recovery" = "recovery_glut1",
  "setup,sample,pop,level,n,median_recovery,q25,q75,squeezed_share" = "spike",
  "source,target,ours,zip,zip_reading" = "ssm_zip",
  "swapped_dye,channel,marker,pop,median_d_squeezed,min_d_squeezed,median_rsd_ratio,median_d_pos,channel_label" = "dye_channel_map",
  "variant,detector,delta,version" = "delta_overlap",
  "variant_dye,part,channel_dye,coef,n_variants,coef_shown" = "overlap_matrix",
  "version,row,measure,value" = "variant_cosines",
  "method,file,event_id,channel,value" = "events_dist",
  "method,file,event_id,x,y" = "events_biplot",
  "method,fluorophore,detector,value_raw,value_norm,angle_deg,panel_label" = "spectra_overlay",
  "arm,sample,pop,channel,n,rsd_neg,pos_share,squeezed_share,ref_rsd,rsd_ratio" = "scores",
  "donor,arm,readout,value" = "readouts",
  "experiment,arm,donor,pct_tnf,pct_ifng,n_ifng_pos,p_pos,p_neg,odds_ratio" = "coherence",
  "row_dye,col_dye,ssm,label" = "ssm_long",
  "source,target,n_pos,n_neg,ss" = "own_ssm",
  "donor,pop,arm,n_cells,share_tnfpos,arm_label" = "share_tnfpos"
)

# Axis labels of the ranking and cytokine metrics, copied from scripts 30, 31 and 34.
METRIC_LABELS <- c(
  fig_spread_sym = "Median over channels of Nguyen rSD (q84 - q50, linear)",
  fig_stain_index_sym = "Median over channels of stain index",
  fig_hyperneg_sym = "Median over channels of hypernegatives (%)",
  fig_streak_sym = sprintf(
    "%s-CD8 Spearman rho (%s vs %s)",
    TNF_LABEL,
    TNF_SLOT,
    CD8_CHANNEL
  ),
  fig_tilt_sym = "Median |Spearman rho| over channel pairs",
  fig_ssm_sym = "Median quadrature SS over source i, target j",
  fig_spread = "Median over channels of Nguyen rSD (q84 - q50, linear)",
  fig_ssm = "Median SS over source i, target j",
  fig_tilt = "Median |Spearman rho| over channel pairs",
  fig_pct_pos_own = "% positive, own gate",
  fig_stain_index = "Stain index (linear)",
  fig_neg_rsd = "Negative Nguyen rSD (linear, MFI)",
  fig_pileup = "Fraction of negatives with |x| < 0.5 x reference rSD",
  fig_emd = "EMD vs as188_joint (arcsinh units)",
  fig_biology = "Value"
)

PAIRED <- sprintf(
  "Each point is one donor (colour and shape); the bar is the median of %d donors.",
  N_DONORS
)
FLAGGED <- "A grey ring marks a value from a failed own gate or a pinned reference gate."

# Paired donors per arm, with an optional ring on flagged values.
arm_donor <- function(
  d,
  value,
  label,
  title,
  rows = NULL,
  cols = NULL,
  flag = NULL
) {
  out <- panel_arm_donor(
    d,
    value,
    str_wrap(label, 45),
    title,
    unique(d$arm),
    rows = rows,
    cols = cols
  )
  if (!is.null(flag)) {
    out$plot <- out$plot +
      geom_point(
        data = filter(out$data, .data[[flag]]),
        shape = 1,
        size = 2.2,
        colour = "grey30"
      )
  }
  out
}

# Fill scale for unsigned matrices (cosine, spillover spread): white to red from 0.
seq_fill <- function(name, limits = NULL) {
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    limits = limits,
    name = name,
    na.value = "grey80"
  )
}

pseudo_log <- scales::transform_pseudo_log(sigma = PSEUDO_LOG_SIGMA)
EVENT_BREAKS <- c(-1e5, -1e4, 0, 1e4, 1e5, 1e6)
event_axis <- function(scale_fn, name) {
  scale_fn(
    name = name,
    transform = pseudo_log,
    breaks = EVENT_BREAKS,
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  )
}

BUILDERS <- list(
  gap = \(d, src, title) {
    out <- panel_gap(src, title)
    out$caption <- paste(
      sprintf(
        "Fraction of the as188_joint to s102_c188 gap each arm closes, for the CD8T %s squeezed and positive shares.",
        TNF_LABEL
      ),
      "Each point is one donor; the bar is the median; dashed lines mark the D3 criterion (every donor 0.5, median 0.8)."
    )
    out
  },
  cyto = \(d, src, title) {
    d <- mutate(d, open = as.logical(.data$open))
    out <- arm_donor(
      d,
      "value",
      METRIC_LABELS[[title]],
      title,
      rows = "row",
      cols = "col",
      flag = "open"
    )
    out$caption <- paste(
      "Cytokine metric per arm, population (rows) and cytokine (columns).",
      PAIRED,
      FLAGGED
    )
    out
  },
  rank = \(d, src, title) {
    d <- with_donor(d)
    out <- arm_donor(d, "value", METRIC_LABELS[[title]], title)
    out$caption <- paste("Ranking metric per arm.", PAIRED)
    out
  },
  pct_pos = \(d, src, title) {
    flag <- if ("failed" %in% names(d)) "failed" else NULL
    out <- arm_donor(
      d,
      "pct_pos",
      "% positive of parent",
      title,
      rows = "cytokine",
      cols = "parent",
      flag = flag
    )
    out$caption <- paste(
      "Share of cytokine-positive cells per arm, cytokine (rows) and parent population (columns).",
      PAIRED,
      if (is.null(flag)) "" else "A grey ring marks a failed gate."
    )
    out
  },
  density_curves = \(d, src, title) {
    # The table holds every cytokine; the figure name says which one it shows.
    plot_df <- d |>
      filter(.data$cytokine == str_remove(.env$title, "^fig_b_density_")) |>
      mutate(
        donor = factor_by_pal(.data$donor, DONOR_PAL),
        arm = factor(.data$arm, levels = unique(.data$arm))
      )
    p <- ggplot(plot_df, aes(value, density, colour = donor)) +
      geom_line(linewidth = 0.3) +
      scale_colour_manual(values = DONOR_PAL, name = "Donor") +
      facet_grid(
        parent ~ arm,
        labeller = labeller(arm = \(x) str_wrap(arm_label(x), 12)),
        axes = "all",
        axis.labels = "margins"
      ) +
      labs(x = "Cytokine (arcsinh)", y = "Density", title = title)
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 85,
      caption = "Saved density curves of the cytokine channel per arm (columns) and parent (rows); one line per donor."
    )
  },
  spread_coef = \(d, src, title) {
    out <- panel_set_points(
      d,
      "set",
      "spread_coef",
      "AF405 spread coefficient",
      title,
      unique(d$set)
    )
    out$caption <- "AF405 spread coefficient per variant set; each point is one spectral variant, the box gives quartiles."
    out
  },
  screen = \(d, src, title) {
    rsd <- str_detect(title, "rsd")
    plot_df <- d |>
      mutate(
        fill = if (rsd) {
          log2(.data$median_rsd_ratio)
        } else {
          .data$median_d_squeezed
        },
        arm = factor(.data$arm, levels = unique(.data$arm))
      )
    limit <- max(abs(plot_df$fill), na.rm = TRUE)
    out <- panel_signed_tile(
      plot_df,
      "arm",
      "label",
      "fill",
      limit,
      if (rsd) "log2 median rSD ratio" else "Median change in squeezed share",
      title,
      cols = "pop"
    )
    out$caption <- paste(
      sprintf(
        "Change per marker (rows) and arm (columns) against the reference, median of %d donors, per population.",
        N_DONORS
      ),
      "Blue is lower, red higher; grey is not tested."
    )
    out
  },
  rho_pairs = \(d, src, title) {
    out <- arm_donor(
      d,
      "rho",
      "Spearman rho",
      title,
      rows = "pair",
      cols = "parent"
    )
    out$caption <- paste(
      "Spearman rho between cytokine pairs (rows) per parent population (columns).",
      PAIRED
    )
    out
  },
  dye_share = \(d, src, title) {
    d <- d |> with_donor() |> mutate(arm = .data$dye)
    out <- arm_donor(
      d,
      "squeezed_share",
      sprintf("CD8T %s squeezed share", TNF_LABEL),
      title
    )
    out$caption <- paste(
      sprintf(
        "CD8T %s squeezed share when only one dye takes 1.0.2 variants (x).",
        TNF_LABEL
      ),
      PAIRED
    )
    out
  },
  cosine = \(d, src, title) {
    lv <- unique(d$fluorophore_x)
    plot_df <- d |>
      mutate(
        fluorophore_x = factor(.data$fluorophore_x, levels = lv),
        fluorophore_y = factor(.data$fluorophore_y, levels = lv)
      )
    p <- ggplot(plot_df, aes(fluorophore_x, fluorophore_y, fill = cosine)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      seq_fill("Cosine", c(0, 1)) +
      coord_fixed() +
      labs(x = NULL, y = NULL, title = title) +
      rotate_x()
    list(
      plot = p,
      data = plot_df,
      width_mm = if (length(lv) <= 12) 89 else 183,
      height_mm = min(183, 40 + 3 * length(lv)),
      caption = "Cosine similarity between fluorophore spectra; 1 means identical shape."
    )
  },
  g1 = \(d, src, title) {
    out <- panel_set_points(
      d,
      "arm",
      "log2_ratio",
      "log2 ratio to reference",
      title,
      unique(d$arm),
      cols = "metric"
    )
    out$caption <- "Held-out controls (G1): log2 ratio of each metric to the reference per arm; each point is one control."
    out
  },
  rho_marker = \(d, src, title) {
    d <- mutate(d, pair = paste(.data$cytokine, "-", .data$marker))
    out <- arm_donor(
      d,
      "rho",
      "Spearman rho",
      title,
      rows = "pair",
      cols = "parent"
    )
    out$caption <- paste(
      "Spearman rho between a cytokine and an activation marker (rows) per parent (columns).",
      PAIRED
    )
    out
  },
  tnf_of_ifng = \(d, src, title) {
    long <- d |>
      pivot_longer(
        c(pct_tnf_of_ifng, pct_tnf_of_ifng_neg),
        names_to = "among",
        values_to = "pct"
      ) |>
      mutate(
        among = if_else(
          .data$among == "pct_tnf_of_ifng",
          paste0(IFNG_LABEL, "+"),
          paste0(IFNG_LABEL, "-")
        )
      )
    out <- arm_donor(
      long,
      "pct",
      paste0("% ", TNF_LABEL, "+"),
      title,
      rows = "among",
      cols = "parent"
    )
    out$caption <- paste(
      sprintf(
        "Share of %s+ cells among %s+ and %s- cells (rows) per parent (columns).",
        TNF_LABEL,
        IFNG_LABEL,
        IFNG_LABEL
      ),
      PAIRED
    )
    out
  },
  g2 = \(d, src, title) {
    d <- mutate(d, log2_ratio = log2(.data$ratio))
    out <- panel_set_points(
      d,
      "arm",
      "log2_ratio",
      "log2 ratio to reference",
      title,
      unique(d$arm),
      cols = "metric"
    )
    out$caption <- "Simulation (G2): log2 ratio of each metric to the reference per arm; each point is one fluor, mean over seeds."
    out
  },
  passes = \(d, src, title) {
    p <- ggplot(
      d,
      aes(factor(n_passes), factor(n_af_passes), fill = fill_scaled)
    ) +
      geom_tile(colour = "white") +
      geom_text(
        aes(label = signif(median, 2), colour = fill_scaled > 0.5),
        size = 6,
        size.unit = "pt",
        show.legend = FALSE
      ) +
      scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
      scale_fill_viridis_c(name = "Scaled median") +
      facet_wrap(~readout, axes = "all", axis.labels = "margins") +
      labs(x = "Unmixing passes", y = "AF passes", title = title)
    list(
      plot = p,
      data = d,
      width_mm = 183,
      height_mm = 80,
      caption = sprintf(
        "Median over %d donors per pass grid cell (text); the fill is that median scaled within each readout.",
        N_DONORS
      )
    )
  },
  anatomy = \(d, src, title) {
    plot_df <- mutate(
      d,
      detector = factor(.data$detector, levels = unique(.data$detector))
    )
    p <- ggplot(
      plot_df,
      aes(detector, value, colour = series, group = series)
    ) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_line(linewidth = 0.3) +
      scale_colour_manual(
        values = categorical_palette(plot_df$series),
        name = NULL
      ) +
      facet_wrap(~cell, axes = "all", axis.labels = "margins") +
      labs(x = "Detector", y = "Signal", title = title) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 120,
      caption = "Per-detector anatomy of single spiked cells (panels): each line is one saved series."
    )
  },
  af_shape = \(d, src, title) {
    # Rows carry either a variant set or a named reference line (set NA).
    plot_df <- d |>
      mutate(
        detector = factor(.data$detector, levels = unique(.data$detector)),
        series = coalesce(.data$set, .data$line)
      )
    p <- ggplot(
      plot_df,
      aes(detector, value, colour = series, group = series)
    ) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_line(linewidth = 0.4) +
      scale_colour_manual(
        values = categorical_palette(plot_df$series),
        name = NULL
      ) +
      labs(x = "Detector", y = "Normalised signal", title = title) +
      # detectors stay on x: a spectrum reads left to right, so it is not flipped
      rotate_x()
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 85,
      caption = sprintf(
        "Mean AF405 variant spread per detector for each variant set, with the AF and %s (1.8.8) reference shapes.",
        TNF_DYE
      )
    )
  },
  interventions = \(d, src, title) {
    d <- mutate(d, arm = .data$arm_label)
    out <- arm_donor(d, "value", "Value", title, rows = "readout", cols = "pop")
    out$caption <- paste(
      "Readouts (rows) per intervention arm and population (columns).",
      PAIRED
    )
    out
  },
  il2 = \(d, src, title) {
    out <- arm_donor(
      d,
      "pct_il2",
      "% IL-2+",
      title,
      rows = "among",
      cols = "parent"
    )
    out$caption <- paste(
      "Share of IL-2+ cells among each cytokine subset (rows) per parent (columns).",
      PAIRED
    )
    out
  },
  dose = \(d, src, title) {
    long <- d |>
      pivot_longer(
        c(pos_share, squeezed_share, rsd_ratio),
        names_to = "readout",
        values_to = "value"
      ) |>
      mutate(donor = factor_by_pal(.data$donor, DONOR_PAL))
    p <- ggplot(long, aes(alpha, value)) +
      geom_line(aes(group = donor), colour = "grey75", linewidth = 0.3) +
      geom_point(aes(colour = donor, shape = donor), size = 1.4) +
      donor_scales() +
      facet_grid(
        readout ~ pop,
        scales = "free_y",
        axes = "all",
        axis.labels = "margins"
      ) +
      scale_x_continuous(breaks = c(0, 0.5, 1)) +
      labs(x = "Variant dose (alpha)", y = NULL, title = title)
    list(
      plot = p,
      data = long,
      width_mm = 89,
      height_mm = 100,
      caption = "Readouts (rows) per variant dose alpha and population (columns); each point is one donor, lines join a donor."
    )
  },
  neg_spread = \(d, src, title) {
    out <- arm_donor(
      d,
      "rsd_ratio",
      "Negative rSD ratio to reference",
      title,
      rows = "cytokine",
      cols = "parent"
    )
    out$caption <- paste(
      "Negative spread (rSD) relative to the reference per cytokine (rows) and parent (columns).",
      PAIRED
    )
    out
  },
  rawaf_density = \(d, src, title) {
    plot_df <- d |>
      rename(donor = sample) |>
      mutate(arm = factor(.data$arm, levels = unique(.data$arm)))
    p <- ggplot(plot_df, aes(bv510_asinh, colour = arm)) +
      geom_density(linewidth = 0.3) +
      scale_colour_manual(
        values = categorical_palette(levels(plot_df$arm)),
        name = "Arm"
      ) +
      facet_wrap(~donor, nrow = 1, axes = "all", axis.labels = "margins") +
      labs(x = paste0(TNF_DYE, " (arcsinh)"), y = "Density", title = title)
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 60,
      caption = sprintf(
        "%s (%s) distribution per arm and donor (panels), drawn from the saved event values.",
        TNF_DYE,
        TNF_LABEL
      )
    )
  },
  glut1 = \(d, src, title) {
    plot_df <- d |>
      rename(donor = sample) |>
      mutate(donor = factor_by_pal(.data$donor, DONOR_PAL))
    p <- ggplot(
      plot_df,
      aes(glut1_decile, share_newly_squeezed, colour = donor, shape = donor)
    ) +
      geom_line(aes(group = donor), linewidth = 0.3) +
      geom_point(size = 1.2) +
      donor_scales() +
      facet_wrap(~pop, axes = "all", axis.labels = "margins") +
      scale_x_continuous(breaks = 1:10) +
      labs(x = "GLUT1 decile", y = "Share newly squeezed", title = title)
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 65,
      caption = "Share of cells newly squeezed near 0 per GLUT1 decile and population; one line per donor."
    )
  },
  offsets = \(d, src, title) {
    p <- ggplot(d, aes(offset, dye, colour = set)) +
      value_grid("x") +
      geom_vline(xintercept = 0, colour = "grey60") +
      geom_point(size = 1.2) +
      scale_colour_manual(values = categorical_palette(d$set), name = "Set") +
      scale_y_discrete(limits = rev) +
      labs(x = "Offset", y = NULL, title = title)
    list(
      plot = p,
      data = d,
      width_mm = 89,
      height_mm = 30 + ROW_MM * n_distinct(d$dye),
      caption = "Offset per dye for each control set."
    )
  },
  variant_set = \(d, src, title) {
    long <- d |>
      pivot_longer(
        c(spread_coef, offset_coef),
        names_to = "coef",
        values_to = "value"
      ) |>
      mutate(panel = paste(.data$dye, .data$coef))
    out <- panel_set_points(
      long,
      "set",
      "value",
      "Coefficient",
      title,
      unique(long$set),
      cols = "panel"
    )
    out$caption <- "Spread and offset coefficients of each spectral variant (points) per control set, by dye (panels)."
    out
  },
  recovery_glut1 = \(d, src, title) {
    plot_df <- d |>
      rename(donor = sample) |>
      mutate(donor = factor_by_pal(.data$donor, DONOR_PAL))
    p <- ggplot(plot_df, aes(glut1_decile, median_recovery, colour = setup)) +
      geom_line(aes(group = interaction(donor, setup)), linewidth = 0.3) +
      scale_colour_manual(
        values = categorical_palette(plot_df$setup),
        name = "Setup"
      ) +
      facet_grid(
        level ~ pop,
        labeller = label_both,
        axes = "all",
        axis.labels = "margins"
      ) +
      scale_x_continuous(breaks = 1:10) +
      labs(x = "GLUT1 decile", y = "Median spike recovery", title = title)
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 110,
      caption = "Median recovery of the spiked signal per GLUT1 decile, spike level (rows) and population; one line per donor and setup."
    )
  },
  spike = \(d, src, title) {
    plot_df <- d |>
      rename(donor = sample) |>
      mutate(donor = factor_by_pal(.data$donor, DONOR_PAL))
    p <- ggplot(plot_df, aes(level, median_recovery)) +
      geom_line(aes(group = donor), colour = "grey75", linewidth = 0.3) +
      geom_point(aes(colour = donor, shape = donor), size = 1.4) +
      donor_scales() +
      facet_grid(pop ~ setup, axes = "all", axis.labels = "margins") +
      labs(x = "Spike level", y = "Median spike recovery", title = title)
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 80,
      caption = "Median recovery of the spiked signal per spike level, population (rows) and setup (columns); each point is one donor."
    )
  },
  ssm_zip = \(d, src, title) {
    p <- ggplot(d, aes(zip, ours)) +
      geom_point(size = 0.6, alpha = 0.5) +
      scale_x_log10() +
      scale_y_log10() +
      facet_wrap(~zip_reading, axes = "all", axis.labels = "margins") +
      # both axes are the same quantity, so 1:1
      theme(aspect.ratio = 1) +
      labs(
        x = "SSM from the zip (log scale)",
        y = "SSM, own computation (log scale)",
        title = title
      )
    list(
      plot = p,
      data = d,
      width_mm = 183,
      height_mm = 95,
      caption = "Own spillover spread per dye pair against the zip value, with the zip read as stored and transposed; both axes log scale, so pairs at 0 sit on the axis edge."
    )
  },
  dye_channel_map = \(d, src, title) {
    limit <- max(abs(d$median_d_squeezed), na.rm = TRUE)
    out <- panel_signed_tile(
      d,
      "swapped_dye",
      "channel_label",
      "median_d_squeezed",
      limit,
      "Median change in squeezed share",
      title,
      cols = "pop"
    )
    out$height_mm <- min(230, out$height_mm)
    out$caption <- "Median change in squeezed share per channel (rows) when one dye (columns) takes 1.0.2 variants, per population."
    out
  },
  delta_overlap = \(d, src, title) {
    plot_df <- mutate(
      d,
      detector = factor(.data$detector, levels = unique(.data$detector))
    )
    p <- ggplot(plot_df, aes(detector, delta, group = variant)) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_line(linewidth = 0.2, alpha = 0.3) +
      facet_wrap(~version, ncol = 1, axes = "all", axis.labels = "margins") +
      labs(x = "Detector", y = "Variant minus mean spectrum", title = title) +
      rotate_x()
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 100,
      caption = "Deviation of each AF405 spectral variant (lines) from its mean spectrum, per detector and code version."
    )
  },
  overlap_matrix = \(d, src, title) {
    limit <- max(abs(d$coef_shown), na.rm = TRUE)
    out <- panel_signed_tile(
      d,
      "channel_dye",
      "variant_dye",
      "coef_shown",
      limit,
      "Coefficient",
      title,
      cols = "part"
    )
    out$caption <- "Overlap coefficient of each dye's variants (rows) onto each channel dye (columns), split by part."
    out
  },
  variant_cosines = \(d, src, title) {
    out <- panel_set_points(
      d,
      "version",
      "value",
      "Value",
      title,
      unique(d$version),
      cols = "measure"
    )
    out$caption <- "Geometry of each AF405 variant delta (points) per code version, for four measures (panels)."
    out
  },
  events_dist = \(d, src, title) {
    plot_df <- mutate(
      d,
      method = factor(.data$method, levels = unique(.data$method))
    )
    n_ch <- n_distinct(plot_df$channel)
    p <- ggplot(plot_df, aes(value, colour = method)) +
      geom_density(linewidth = 0.3, key_glyph = "path") +
      event_axis(scale_x_continuous, "Unmixed value (pseudo-log, sigma 1000)") +
      scale_colour_manual(
        values = categorical_palette(levels(plot_df$method)),
        name = "Method"
      ) +
      facet_wrap(
        ~channel,
        scales = "free_y",
        ncol = min(6, n_ch),
        axes = "all",
        axis.labels = "margins"
      ) +
      labs(y = "Density", title = title) +
      rotate_x()
    list(
      plot = p,
      data = slice_sample(plot_df, n = min(nrow(plot_df), 20000)),
      width_mm = 183,
      height_mm = 30 + 30 * ceiling(n_ch / 6),
      caption = "Distribution of unmixed values per channel (panels) and method, from the saved event sample."
    )
  },
  events_biplot = \(d, src, title) {
    plot_df <- mutate(
      d,
      method = factor(.data$method, levels = unique(.data$method))
    )
    p <- ggplot(plot_df, aes(x, y)) +
      geom_bin2d(bins = BIPLOT_BINS) +
      event_axis(scale_x_continuous, BIPLOT_PAIR[[1]]) +
      event_axis(scale_y_continuous, BIPLOT_PAIR[[2]]) +
      scale_fill_viridis_c(transform = "log10", name = "Events") +
      facet_wrap(~method, nrow = 1, axes = "all", axis.labels = "margins") +
      labs(title = title) +
      rotate_x()
    list(
      plot = p,
      data = slice_sample(plot_df, n = min(nrow(plot_df), 20000)),
      width_mm = 183,
      height_mm = 75,
      caption = "Event density of two unmixed channels per method (panels), from the saved event sample; both axes pseudo-log (sigma 1000)."
    )
  },
  spectra_overlay = \(d, src, title) {
    plot_df <- mutate(
      d,
      detector = factor(.data$detector, levels = unique(.data$detector))
    )
    n_fl <- n_distinct(plot_df$panel_label)
    p <- ggplot(
      plot_df,
      aes(detector, value_norm, colour = method, group = method)
    ) +
      geom_line(linewidth = 0.4, alpha = 0.7) +
      scale_colour_manual(
        values = categorical_palette(plot_df$method),
        name = "Method"
      ) +
      facet_wrap(
        ~panel_label,
        ncol = min(6, n_fl),
        axes = "all",
        axis.labels = "margins"
      ) +
      labs(x = "Detector", y = "Normalised signal", title = title) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = min(230, 30 + 25 * ceiling(n_fl / 6)),
      caption = "Fluorophore spectra per method (panels give the angle between methods in degrees)."
    )
  },
  scores = \(d, src, title) {
    # Two channels per table: the TNF and IFNg dye slots, named by the dataset's labels.
    long <- d |>
      with_donor() |>
      mutate(
        panel = paste(
          .data$pop,
          setNames(c(TNF_LABEL, IFNG_LABEL), c(TNF_SLOT, IFNG_SLOT))[
            .data$channel
          ]
        )
      ) |>
      pivot_longer(
        c(squeezed_share, pos_share, rsd_ratio),
        names_to = "readout",
        values_to = "value"
      )
    out <- arm_donor(
      long,
      "value",
      "Value",
      title,
      rows = "readout",
      cols = "panel"
    )
    out$width_mm <- 183
    out$caption <- paste(
      sprintf(
        "%s and %s readouts (rows) per arm, population and cytokine (columns).",
        TNF_LABEL,
        IFNG_LABEL
      ),
      PAIRED
    )
    out
  },
  readouts = \(d, src, title) {
    out <- arm_donor(d, "value", "Value", title, rows = "readout")
    out$caption <- paste("Readouts (columns) per arm.", PAIRED)
    out
  },
  coherence = \(d, src, title) {
    plot_df <- d |>
      mutate(
        donor = factor_by_pal(.data$donor, DONOR_PAL),
        experiment = factor(
          .data$experiment,
          levels = unique(.data$experiment)
        ),
        arm = factor(.data$arm, levels = unique(.data$arm))
      )
    p <- ggplot(plot_df, aes(odds_ratio, arm)) +
      value_grid("x") +
      geom_vline(xintercept = 1, colour = "grey60") +
      stat_summary(
        fun = median,
        fun.min = median,
        fun.max = median,
        geom = "crossbar",
        width = 0.5,
        linewidth = 0.2,
        orientation = "y"
      ) +
      geom_point(aes(colour = donor, shape = donor), size = 1.1) +
      donor_scales() +
      scale_y_discrete(limits = rev, labels = arm_label) +
      scale_x_log10() +
      facet_grid(
        experiment ~ .,
        scales = "free_y",
        space = "free_y",
        axes = "all",
        axis.labels = "margins"
      ) +
      labs(
        x = sprintf(
          "Odds ratio %s+ given %s+ (log scale)",
          TNF_LABEL,
          IFNG_LABEL
        ),
        y = NULL,
        title = title
      )
    list(
      plot = p,
      data = plot_df,
      width_mm = 89,
      height_mm = min(
        240,
        30 +
          8 * n_distinct(plot_df$experiment) +
          ROW_MM * nrow(distinct(plot_df, experiment, arm))
      ),
      caption = paste(
        sprintf(
          "Odds of %s+ in %s+ versus %s- cells per arm, grouped by experiment.",
          TNF_LABEL,
          IFNG_LABEL,
          IFNG_LABEL
        ),
        PAIRED
      )
    )
  },
  ssm_long = \(d, src, title) {
    lv <- unique(d$row_dye)
    plot_df <- mutate(
      d,
      row_dye = factor(.data$row_dye, levels = rev(lv)),
      col_dye = factor(.data$col_dye, levels = lv)
    )
    p <- ggplot(plot_df, aes(col_dye, row_dye, fill = ssm)) +
      geom_tile() +
      seq_fill("SSM") +
      coord_fixed() +
      labs(x = "Receiving dye", y = "Spreading dye", title = title) +
      rotate_x() +
      theme(axis.text = element_text(size = 6))
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 183,
      caption = "Spillover spread matrix of the 1.8.8 unmix; empty cells (diagonal) are grey."
    )
  },
  own_ssm = \(d, src, title) {
    lv <- unique(d$source)
    plot_df <- mutate(
      d,
      source = factor(.data$source, levels = rev(lv)),
      target = factor(.data$target, levels = lv)
    )
    p <- ggplot(plot_df, aes(target, source, fill = ss)) +
      geom_tile() +
      seq_fill("SS") +
      coord_fixed() +
      labs(x = "Target j", y = "Source i", title = title) +
      rotate_x() +
      theme(axis.text = element_text(size = 6))
    list(
      plot = p,
      data = plot_df,
      width_mm = 183,
      height_mm = 183,
      caption = "Own quadrature spillover spread from source dye i into target channel j."
    )
  },
  share_tnfpos = \(d, src, title) {
    d <- mutate(d, arm = .data$arm_label)
    out <- arm_donor(
      d,
      "share_tnfpos",
      sprintf(
        "Share %s+ among %s+ %s- source",
        TNF_LABEL,
        IFNG_LABEL,
        TNF_LABEL
      ),
      title,
      cols = "pop"
    )
    out$caption <- paste(
      sprintf(
        "Share of %s+ cells per arm and population (columns).",
        TNF_LABEL
      ),
      PAIRED
    )
    out
  }
)

# Redraw one figure -----------------------------------------------------------------------------
# Reads the saved table, builds the layout's panel, writes it with save_panel(); returns the PNG
# path and the caption. Event panels write a seeded 20,000-event subsample as their source CSV.
# Builders that draw their own panels do not split by condition; a dataset with conditions stops there.
# 91_all_figures.qmd should key LAYOUTS with table_signature(), and stop on a mapped table with no layout.
CONDITION_LAYOUTS <- c(
  "cyto",
  "rank",
  "pct_pos",
  "rho_pairs",
  "dye_share",
  "rho_marker",
  "tnf_of_ifng",
  "interventions",
  "il2",
  "neg_spread",
  "scores",
  "readouts",
  "share_tnfpos"
)
# Tables the scripts write with `donor` = sample id (24H: BMV_006C_SEB): put the donor id into `donor`,
# which keys DONOR_PAL and DONOR_SHAPES, and add `condition`. Tables with a `sample` column go through
# with_donor() in their builder; tables that already carry `condition` are left as they are.
with_donor_ids <- function(d) {
  if (
    is.null(DS$conditions) ||
      !"donor" %in% names(d) ||
      "condition" %in% names(d) ||
      !all(d$donor %in% samples_tbl()$sample)
  ) {
    return(d)
  }
  d |>
    left_join(
      select(samples_tbl(), donor = sample, donor_id = donor, condition),
      by = join_by(donor),
      relationship = "many-to-one",
      unmatched = "error"
    ) |>
    mutate(donor = .data$donor_id) |>
    select(-donor_id)
}
# Column signature of a saved table without `condition`, the key into LAYOUTS (a dataset with conditions
# may save the same table with or without that column)
table_signature <- function(path) {
  paste(
    setdiff(
      names(read_csv(path, n_max = 1, show_col_types = FALSE)),
      "condition"
    ),
    collapse = ","
  )
}
redraw_panel <- function(src, layout, title, out_dir, name) {
  if (!is.null(DS$conditions) && !layout %in% CONDITION_LAYOUTS) {
    stop(
      "report_figures.R redraw_panel: layout '",
      layout,
      "' does not split by condition (dataset ",
      DS$dataset,
      "); figure ",
      name,
      call. = FALSE
    )
  }
  d <- read_csv(src, show_col_types = FALSE) |>
    with_donor_ids()
  set.seed(1)
  out <- BUILDERS[[layout]](d, src, title)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  save_panel(out$plot, out$data, out_dir, name, out$width_mm, out$height_mm)
  tibble(
    redrawn = file.path(out_dir, paste0(name, ".png")),
    caption = out$caption
  )
}
