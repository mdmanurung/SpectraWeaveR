# The common plot forms, drawn the same way by every figure script; plots.py is
# the plotnine twin, layer for layer.
# Source theme.R and palettes.R first; each function returns a ggplot object and
# never saves it, so the caller passes it to save_panel().
# `title` takes one short descriptive title from the caller ("CD14", "CD14+
# monocyte pseudobulks (top 500 HVGs)"); no function invents one, and no
# subtitle or caption is drawn: the caption goes in the report (fig.text).
# Columns are passed as strings. Statistics stop at the lm fit in
# plot_scatter_fit() and the plot_box() brackets the caller asks for.
# Colours, theme and dashed threshold lines are fixed house style. Every tunable
# value is a required argument; a grill asks the user for it per project:
#   padj_max:       which adjusted-p cutoff marks a gene significant?
#   lfc_min:        which absolute log2 fold change marks a gene significant?
#   limit:          at which absolute value should the heatmap colour scale saturate?
#   jitter_width:   how wide should the jitter spread points within a box?
#   seed:           which seed fixes the jitter so reruns draw the same panel?
#   max_size:       what is the largest dot size in the dotplot, in mm?
#   test, p_adjust: which test and adjustment for the box-plot brackets (fig.stats)?
#                   Both or neither, passed to ggpubr::stat_pwc(): test
#                   "wilcox_test" or "t_test"; p_adjust any p.adjust() method
#                   ("BH", "holm", ...). plots.py takes "wilcox.test" or "t.test".

# fig.wilke: under this many points per group a box misleads; draw points + median.
FEW_POINTS <- 10
# fig.volcano: blue down, red up, grey not significant.
DIRECTION <- c(down = "#2166AC", up = "#B2182B", ns = "grey60")
# fig.heatmap: rotate column labels past this many columns.
ROTATE_PAST <- 8
# Functions, not constants: ggplot2 may load after this file is sourced.
rotate_x <- function() theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# fig.embedding: no frame, ticks, tick labels or axis titles; corner arrows instead.
no_axes <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
}

# Stop with the function name when a required column is absent.
check_cols <- function(data, cols, fn) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf("%s: missing columns: %s", fn, toString(missing)), call. = FALSE)
  }
  invisible(data)
}

# Factor in palette order; palette_for() stops first on any label the palette lacks.
factor_by_pal <- function(x, pal) {
  palette_for(x, pal)
  factor(x, levels = names(pal))
}

# fig.n: each group label over its count, "ctrl\n(n = 4)".
n_labels <- function(groups) {
  n <- table(groups)
  \(g) sprintf("%s\n(n = %d)", g, as.integer(n[g]))
}

# fig.embedding: two short arrows at the bottom-left, "UMAP 1" and "UMAP 2".
corner_arrows <- function(x, y, axis_name) {
  x0 <- min(x)
  y0 <- min(y)
  dx <- 0.2 * (max(x) - x0)
  dy <- 0.2 * (max(y) - y0)
  head <- arrow(length = unit(1.3, "mm"), type = "closed")
  list(
    annotate("segment", x = x0, y = y0, xend = x0 + dx, yend = y0, arrow = head),
    annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + dy, arrow = head),
    annotate("text", x = x0 + dx / 2, y = y0 - dy / 4, label = paste(axis_name, 1), size = 2.5),
    annotate(
      "text",
      x = x0 - dx / 4,
      y = y0 + dy / 2,
      label = paste(axis_name, 2),
      angle = 90,
      size = 2.5
    )
  )
}

# Points, scattermore when installed, for dense point clouds.
embedding_points <- function() {
  if (requireNamespace("scattermore", quietly = TRUE)) {
    return(scattermore::geom_scattermore(pointsize = 2))
  }
  geom_point(size = 0.3)
}

# Significant = padj < padj_max and |lfc| > lfc_min, coloured by direction. With
# `label`, significant genes are labelled: all, or the top `n_label` per direction.
plot_volcano <- function(
  data,
  lfc,
  padj,
  padj_max,
  lfc_min,
  label = NULL,
  n_label = NULL,
  title = NULL
) {
  check_cols(data, c(lfc, padj, label), "plot_volcano")
  significant <- data[[padj]] < padj_max & abs(data[[lfc]]) > lfc_min
  data$direction <- ifelse(significant, ifelse(data[[lfc]] > 0, "up", "down"), "ns")
  p <- ggplot(data, aes(.data[[lfc]], -log10(.data[[padj]]), colour = direction)) +
    geom_point(size = 0.5) +
    geom_hline(yintercept = -log10(padj_max), linetype = "dashed") +
    geom_vline(xintercept = c(-lfc_min, lfc_min), linetype = "dashed") +
    scale_colour_manual(values = DIRECTION, guide = "none") +
    labs(
      x = expression(log[2] ~ "fold change"),
      y = expression(-log[10] ~ "adjusted p"),
      title = title
    )
  if (is.null(label)) {
    return(p)
  }
  shown <- data[data$direction != "ns", ]
  if (!is.null(n_label)) {
    shown <- shown[order(shown[[padj]]), ]
    shown <- shown[stats::ave(shown[[padj]], shown$direction, FUN = seq_along) <= n_label, ]
  }
  p + geom_text(aes(label = .data[[label]]), data = shown, vjust = -0.5, size = 2.5)
}

# Box per group with the rows on top; pair joins the same id across groups.
# `test` and `p_adjust` add adjusted-p brackets (fig.stats); with `pair` the test
# is paired on the pair id.
plot_box <- function(
  data,
  group,
  value,
  pal,
  jitter_width,
  seed,
  value_label,
  pair = NULL,
  test = NULL,
  p_adjust = NULL,
  title = NULL
) {
  check_cols(data, c(group, value, pair), "plot_box")
  if (is.null(test) != is.null(p_adjust)) {
    stop("plot_box: pass test and p_adjust together, or neither", call. = FALSE)
  }
  data[[group]] <- factor_by_pal(data[[group]], pal)
  if (!is.null(pair)) {
    # stat_pwc pairs values by row order within each group: sort by pair id.
    data <- data[order(data[[pair]]), ]
    if (!is.null(test) && length(unique(table(data[[group]]))) > 1) {
      stop("plot_box: a paired test needs every pair in every group", call. = FALSE)
    }
  }
  p <- ggplot(data, aes(.data[[group]], .data[[value]])) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_colour_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = n_labels(data[[group]])) +
    labs(x = NULL, y = value_label, title = title)
  # fig.wilke: under FEW_POINTS per group a box misleads; draw points and a median.
  p <- if (min(table(data[[group]])) < FEW_POINTS) {
    p +
      stat_summary(
        aes(colour = .data[[group]]),
        fun = median,
        fun.min = median,
        fun.max = median,
        geom = "crossbar",
        width = 0.5
      )
  } else {
    p + geom_boxplot(aes(fill = .data[[group]]), outlier.shape = NA)
  }
  p <- if (is.null(pair)) {
    p + geom_point(position = position_jitter(width = jitter_width, height = 0, seed = seed))
  } else {
    # No jitter with pairs, so each line meets its points.
    p + geom_line(aes(group = .data[[pair]]), colour = "grey60") + geom_point()
  }
  if (is.null(test)) {
    return(p)
  }
  p +
    ggpubr::stat_pwc(
      method = test,
      method.args = list(paired = !is.null(pair)),
      p.adjust.method = p_adjust,
      label = "p.adj = {p.adj.format}",
      # ggpubr's default is 3.88 mm (about 11 pt); match the 8 pt theme.
      label.size = 7 / .pt
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
}

# UMAP or similar. `axis_name` names the embedding ("UMAP"); the corner arrows
# read "UMAP 1" and "UMAP 2".
plot_embedding <- function(data, x, y, colour, pal, axis_name, title = NULL) {
  check_cols(data, c(x, y, colour), "plot_embedding")
  data[[colour]] <- factor_by_pal(data[[colour]], pal)
  ggplot(data, aes(.data[[x]], .data[[y]], colour = .data[[colour]])) +
    embedding_points() +
    scale_colour_manual(values = pal) +
    coord_equal() +
    corner_arrows(data[[x]], data[[y]], axis_name) +
    labs(colour = colour, title = title) +
    no_axes()
}

# fig.continuous: viridis colour, highest values drawn last (on top).
plot_feature <- function(data, x, y, value, axis_name, value_label, title = NULL) {
  check_cols(data, c(x, y, value), "plot_feature")
  data <- data[order(data[[value]]), ]
  ggplot(data, aes(.data[[x]], .data[[y]], colour = .data[[value]])) +
    geom_point(size = 0.3) +
    scale_colour_viridis_c() +
    coord_equal() +
    corner_arrows(data[[x]], data[[y]], axis_name) +
    labs(colour = value_label, title = title) +
    no_axes()
}

# Diverging fill centred at 0 with symmetric limits; values beyond limit are squished.
# Row and column order are the caller's factor levels, first row on top: no clustering.
# `col_annotation` names a column of `data` that maps each column to a group,
# drawn as a strip above the heatmap in `annotation_pal` colours.
plot_heatmap <- function(
  data,
  row,
  col,
  value,
  limit,
  col_annotation = NULL,
  annotation_pal = NULL,
  title = NULL,
  x_label = col,
  y_label = row,
  fill_label = value
) {
  check_cols(data, c(row, col, value, col_annotation), "plot_heatmap")
  # Caller factor levels win; otherwise first-seen order, as in plots.py.
  for (key in c(row, col)) {
    if (!is.factor(data[[key]])) data[[key]] <- factor(data[[key]], unique(data[[key]]))
  }
  p <- ggplot(data, aes(.data[[col]], .data[[row]], fill = .data[[value]])) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-limit, limit),
      oob = scales::squish
    ) +
    scale_y_discrete(limits = rev) +
    labs(x = x_label, y = y_label, fill = fill_label, title = title)
  if (nlevels(data[[col]]) > ROTATE_PAST) {
    p <- p + rotate_x()
  }
  if (is.null(col_annotation)) {
    return(p)
  }
  strip <- data[!duplicated(data[[col]]), c(col, col_annotation)]
  ann <- data.frame(
    x = as.integer(strip[[col]]),
    y = nlevels(data[[row]]) + 1,
    group = factor_by_pal(strip[[col_annotation]], annotation_pal)
  )
  # fig.heatmap: one butt-ended segment per column reads as a solid strip.
  p +
    geom_segment(
      aes(x = x - 0.5, xend = x + 0.5, y = y, yend = y, colour = group),
      data = ann,
      linewidth = 3,
      lineend = "butt",
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = annotation_pal) +
    labs(colour = col_annotation)
}

# fig.dotplot: size = fraction expressing (0-1), viridis colour = mean.
# `mean_label` titles the colourbar with its unit ("Mean expression (log1p CP10K)");
# the size legend is always "Fraction expressing". First group on top.
plot_dotplot <- function(data, feature, group, mean, frac, max_size, mean_label, title = NULL) {
  check_cols(data, c(feature, group, mean, frac), "plot_dotplot")
  ggplot(data, aes(.data[[feature]], .data[[group]])) +
    geom_point(aes(size = .data[[frac]], colour = .data[[mean]])) +
    scale_size_area(
      max_size = max_size,
      limits = c(0, 1),
      breaks = c(0.25, 0.5, 0.75, 1),
      labels = scales::percent
    ) +
    scale_colour_viridis_c() +
    scale_y_discrete(limits = rev) +
    guides(size = guide_legend(override.aes = list(colour = "grey60"))) +
    labs(
      x = feature,
      y = group,
      size = "Fraction expressing",
      colour = mean_label,
      title = title
    ) +
    rotate_x()
}

# Stacked bar of each label's share of the count per sample.
plot_composition <- function(
  data,
  sample,
  label,
  count,
  pal,
  title = NULL,
  x_label = sample,
  y_label = sprintf("%% of %s", count),
  fill_label = label
) {
  check_cols(data, c(sample, label, count), "plot_composition")
  data[[label]] <- factor_by_pal(data[[label]], pal)
  ggplot(data, aes(.data[[sample]], .data[[count]], fill = .data[[label]])) +
    geom_col(position = "fill") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = x_label, y = y_label, fill = fill_label, title = title)
}

# Points with a linear fit and its 95% confidence band.
plot_scatter_fit <- function(
  data,
  x,
  y,
  label = NULL,
  title = NULL,
  x_label = x,
  y_label = y
) {
  check_cols(data, c(x, y, label), "plot_scatter_fit")
  p <- ggplot(data, aes(.data[[x]], .data[[y]])) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95) +
    labs(x = x_label, y = y_label, title = title)
  if (is.null(label)) {
    return(p)
  }
  p + geom_text(aes(label = .data[[label]]), vjust = -0.5, size = 2.5)
}

if (sys.nframe() == 0) {
  source("theme.R")
  source("palettes.R")
  # A null device, so no stray Rplots.pdf lands in templates/.
  pdf(NULL)
  donors <- data.frame(
    donor_id = rep(sprintf("D%02d", 1:4), times = 2),
    condition = rep(names(pal_condition), each = 4),
    log2_expr = c(5.1, 4.8, 5.4, 5.0, 8.2, 7.9, 8.5, 8.4)
  )
  de <- data.frame(
    gene = c("IFI27", "ACTB", "CD3E", "MX1"),
    log2_fc = c(3.1, 0.1, -1.2, 2.0),
    q_value = c(0.001, 0.9, 0.04, 0.2)
  )
  cells <- data.frame(
    umap_1 = c(0.1, 0.5, 1.2, 2.0, 2.4, 3.1),
    umap_2 = c(1.0, 0.2, 2.2, 1.5, 0.4, 2.8),
    cell_type = rep(names(pal_cell_type), times = 2)
  )
  grid <- tidyr::expand_grid(
    gene = c("IFI27", "ACTB"),
    cell_type = names(pal_cell_type)
  )
  grid$z <- c(-5, -1, 0, 1, 2, 5)
  grid$frac <- c(0.1, 0.5, 0.9, 0.3, 0.7, 1.0)
  # Ten samples, so the heatmap rotates its labels and draws the strip.
  wide <- data.frame(
    gene = rep(c("IFI27", "ACTB"), each = 10),
    sample = rep(sprintf("S%02d", 1:10), times = 2),
    condition = rep(rep(names(pal_condition), each = 5), times = 2),
    z = seq(-2, 2, length.out = 20)
  )
  cells$cd3e <- c(0, 0.5, 1.2, 3.1, 0.2, 2.2)
  many <- data.frame(
    cell_type = rep(names(pal_cell_type), each = 12),
    cd3e = rep(seq(0, 3, length.out = 12), times = 3)
  )
  counts <- data.frame(
    sample = rep(c("S1", "S2"), each = 3),
    cell_type = rep(names(pal_cell_type), times = 2),
    n_cells = c(10, 20, 30, 5, 5, 40)
  )
  panels <- list(
    volcano = list(plot_volcano(de, "log2_fc", "q_value", 0.05, 1, label = "gene"), de),
    box = list(
      plot_box(donors, "condition", "log2_expr", pal_condition, 0.2, 1, "log2 expression"),
      donors
    ),
    box_paired = list(
      plot_box(
        donors,
        "condition",
        "log2_expr",
        pal_condition,
        0.2,
        1,
        "log2 expression",
        pair = "donor_id"
      ),
      donors
    ),
    volcano_top = list(
      plot_volcano(
        de,
        "log2_fc",
        "q_value",
        0.5,
        1,
        label = "gene",
        n_label = 1,
        title = "Disease vs healthy"
      ),
      de
    ),
    box_many = list(
      plot_box(many, "cell_type", "cd3e", pal_cell_type, 0.2, 1, "CD3E (log1p)"),
      many
    ),
    feature = list(
      plot_feature(
        cells,
        "umap_1",
        "umap_2",
        "cd3e",
        "UMAP",
        "CD3E (log1p CP10K)",
        title = "CD3E"
      ),
      cells
    ),
    heatmap_annotated = list(
      plot_heatmap(
        wide,
        "gene",
        "sample",
        "z",
        2,
        col_annotation = "condition",
        annotation_pal = pal_condition
      ),
      wide
    ),
    embedding = list(
      plot_embedding(cells, "umap_1", "umap_2", "cell_type", pal_cell_type, "UMAP"),
      cells
    ),
    heatmap = list(plot_heatmap(grid, "cell_type", "gene", "z", 4), grid),
    dotplot = list(
      plot_dotplot(grid, "gene", "cell_type", "z", "frac", 4, "Mean expression (z-score)"),
      grid
    ),
    composition = list(
      plot_composition(counts, "sample", "cell_type", "n_cells", pal_cell_type),
      counts
    ),
    scatter_fit = list(plot_scatter_fit(de, "log2_fc", "q_value", label = "gene"), de)
  )
  # ggpubr draws the brackets (theme.R already needs it).
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    panels$box_test <- list(
      plot_box(
        donors,
        "condition",
        "log2_expr",
        pal_condition,
        0.2,
        1,
        "log2 expression",
        test = "wilcox_test",
        p_adjust = "BH",
        title = "CD14"
      ),
      donors
    )
    panels$box_paired_test <- list(
      plot_box(
        donors,
        "condition",
        "log2_expr",
        pal_condition,
        0.2,
        1,
        "log2 expression",
        pair = "donor_id",
        test = "t_test",
        p_adjust = "holm"
      ),
      donors
    )
  }
  # A failed stat_pwc only warns and drops the bracket; make that fail the check.
  for (name in intersect(c("box_test", "box_paired_test"), names(panels))) {
    withCallingHandlers(
      ggplot_build(panels[[name]][[1]]),
      warning = \(w) stop(sprintf("%s: %s", name, conditionMessage(w)), call. = FALSE)
    )
  }
  run_dir <- tempfile()
  dir.create(run_dir)
  for (name in names(panels)) {
    save_panel(panels[[name]][[1]], panels[[name]][[2]], run_dir, name, 89, 60)
    stopifnot(file.exists(file.path(run_dir, paste0(name, ".pdf"))))
  }
  stopifnot(inherits(try(plot_volcano(de, "lfc", "q_value", 0.05, 1), silent = TRUE), "try-error"))
  stopifnot(inherits(
    try(
      plot_box(donors, "condition", "log2_expr", pal_condition, 0.2, 1, "x", test = "t_test"),
      silent = TRUE
    ),
    "try-error"
  ))
  donors$condition[[1]] <- "unknown"
  stopifnot(inherits(
    try(
      plot_box(donors, "condition", "log2_expr", pal_condition, 0.2, 1, "log2 expression"),
      silent = TRUE
    ),
    "try-error"
  ))
  invisible(dev.off())
}
