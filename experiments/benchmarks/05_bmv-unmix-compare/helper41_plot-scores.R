# Purpose: the two figures shared by 37b_score_swap.R, 39b_score_variant_count.R and 41_replot_scores.R:
#          BV510 density per arm (donor x population) and the three per-donor scores per arm.
# Run:     source("experiments/benchmarks/05_bmv-unmix-compare/helper41_plot-scores.R") from the repo root.
# Inputs:  none; the callers pass the data.
# Outputs: plot_bv510_density(), plot_scores(): ggplot objects.

# BV510 density per arm --------------------------------------------------------------------
# dens: one row per event with arm, sample, pop, arcsinh; neg_cut: dotted line, arcsinh
plot_bv510_density <- function(dens, arm_colors, neg_cut) {
  ggplot2::ggplot(dens, ggplot2::aes(arcsinh, colour = arm)) +
    ggplot2::geom_density(linewidth = 0.5, adjust = 0.5) +
    ggplot2::geom_vline(
      xintercept = neg_cut,
      linetype = "dotted",
      linewidth = 0.3
    ) +
    ggplot2::facet_grid(pop ~ sample, scales = "free_y") +
    ggplot2::scale_colour_manual(values = arm_colors) +
    ggplot2::scale_x_continuous(breaks = c(0, 2, 4)) +
    ggplot2::coord_cartesian(xlim = c(-1.5, 4)) +
    ggplot2::labs(x = "BV510-A (TNF), arcsinh", y = "Density", colour = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::theme(legend.position = "top")
}

# Three scores per arm, donors as points joined by lines ------------------------------------
# scores: scores.csv of 37b or 39b (arm, sample, pop, channel, rsd_ratio, squeezed_share, pos_share)
plot_scores <- function(scores, arm_colors) {
  metric_labels <- c(
    rsd_ratio = "rSD of negatives\n(ratio to as188_joint)",
    squeezed_share = "Share squeezed\nnear 0",
    pos_share = "Share above\narcsinh 1"
  )
  channel_labels <- c(`BV510-A` = "TNF", `BV570-A` = "IFNg")
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
  ggplot2::ggplot(plot_frame, ggplot2::aes(arm, value, group = sample)) +
    ggplot2::geom_line(colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = arm), size = 1.8) +
    ggplot2::facet_grid(row ~ metric, scales = "free_y") +
    ggplot2::scale_colour_manual(values = arm_colors, guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.margin = ggplot2::margin(5, 5, 5, 30) # room for the first rotated arm name
    )
}
