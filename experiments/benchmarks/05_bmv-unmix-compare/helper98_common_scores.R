# Purpose: the two figures shared by 37b, 39b and 43 on the common gates (plans/common_gates_plan.md, batch 3):
#          BV510 density per arm (donor x population) and the three per-donor scores per arm. Replaces
#          helper41_plot-scores.R for the *_common.R copies: lib/theme.R sizes, density curves precomputed so
#          the saved source csv is small, "above cut" labelled with the cut that was used.
# Run:     source("experiments/benchmarks/05_bmv-unmix-compare/helper98_common_scores.R") after lib/theme.R.
#          SW_DATASET (3d, 24h or exv) is required; this file sources datasets/dataset.R itself.
# Inputs:  none; the callers pass the data. Channel names and the marker labels come from datasets/<ds>.json
#          (channels TNF_slot and IFNg_slot, targets tnf_label and ifng_label).
# Outputs: density_curves(), plot_density_curves(), plot_arm_scores(): data frame / ggplot objects.

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg()

# Density curve per arm, donor and population, same bandwidth rule as geom_density(adjust = 0.5) -------
# dens: one row per event with arm, sample, pop, arcsinh
density_curves <- function(dens, from, to, adjust = 0.5, n = 256) {
  dens |>
    dplyr::reframe(
      {
        d <- stats::density(
          arcsinh,
          adjust = adjust,
          from = from,
          to = to,
          n = n
        )
        tibble::tibble(x = d$x, density = d$y)
      },
      .by = c(arm, sample, pop)
    )
}

# neg_cut: one number (same cut everywhere) or a data frame with columns sample and cut (one cut per sample)
plot_density_curves <- function(curves, arm_colors, neg_cut) {
  cut_layer <- if (is.data.frame(neg_cut)) {
    ggplot2::geom_vline(
      data = neg_cut,
      ggplot2::aes(xintercept = cut),
      linetype = "dotted",
      linewidth = 0.3
    )
  } else {
    ggplot2::geom_vline(
      xintercept = neg_cut,
      linetype = "dotted",
      linewidth = 0.3
    )
  }
  ggplot2::ggplot(curves, ggplot2::aes(x, density, colour = arm)) +
    ggplot2::geom_line(linewidth = 0.3) +
    cut_layer +
    ggplot2::facet_grid(pop ~ sample, scales = "free_y") +
    ggplot2::scale_colour_manual(values = arm_colors) +
    ggplot2::scale_x_continuous(breaks = c(0, 2, 4)) +
    ggplot2::coord_cartesian(xlim = c(-1.5, 4)) +
    ggplot2::labs(
      x = sprintf(
        "%s (%s), arcsinh",
        cfg("channels")$TNF_slot,
        cfg("targets")$tnf_label
      ),
      y = "Density",
      colour = NULL
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::theme(legend.position = "top")
}

# Three scores per arm, donors as points joined by lines -----------------------------------------------
# scores: scores.csv of 37b, 39b or 43 (arm, sample, pop, channel, rsd_ratio, squeezed_share, pos_share)
# pos_cut: one number (labelled in the panel title) or a vector that differs by sample (labelled "per-sample cut")
plot_arm_scores <- function(scores, arm_colors, pos_cut) {
  metric_labels <- c(
    rsd_ratio = "rSD of negatives\n(ratio to as188_joint)",
    squeezed_share = "Share squeezed\nnear 0",
    pos_share = if (length(unique(pos_cut)) == 1) {
      sprintf("Share above\narcsinh %s", unique(pos_cut))
    } else {
      "Share above\nper-sample cut"
    }
  )
  channel_labels <- c(
    cfg("targets")$tnf_label,
    cfg("targets")$ifng_label
  ) |>
    setNames(c(cfg("channels")$TNF_slot, cfg("channels")$IFNg_slot))
  plot_frame <- scores |>
    dplyr::select(
      arm,
      sample,
      pop,
      channel,
      dplyr::all_of(names(metric_labels))
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(names(metric_labels)),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric = factor(metric_labels[metric], levels = metric_labels),
      row = paste(channel_labels[as.character(channel)], pop)
    )
  p <- ggplot2::ggplot(plot_frame, ggplot2::aes(arm, value, group = sample)) +
    ggplot2::geom_line(colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = arm), size = 1) +
    ggplot2::facet_grid(row ~ metric, scales = "free_y") +
    ggplot2::scale_colour_manual(values = arm_colors, guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = ggplot2::margin(5, 5, 14, 5) # room for the rotated arm names
    )
  list(plot = p, frame = plot_frame)
}
