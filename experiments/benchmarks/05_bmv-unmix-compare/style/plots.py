"""The common plot forms, drawn one fixed way in plotnine; plots.R is the R twin.

Each function takes a long DataFrame and column names, returns a plotnine ggplot
and never saves it: the caller imports theme and passes the plot to save_panel().
Each function copies the columns it needs into a frame with fixed names, so the
layers below mirror plots.R layer for layer.
`title=` takes one short descriptive title from the caller ("CD14", "CD14+
monocyte pseudobulks (top 500 HVGs)"); no function invents one, and no subtitle
or caption is drawn: the caption goes in the report (fig.text).
Colours and line styles are fixed house style. Every tunable value is a
required keyword argument; a grill asks the user for each one per project:
  padj_max:     which adjusted-p cutoff marks a gene significant?
  lfc_min:      which absolute log2 fold change must a significant gene exceed?
  limit:        at which absolute value should the heatmap colour scale saturate?
  seed:         which seed fixes the box-plot jitter?
  jitter_width: how far (in box widths, each side) may box-plot points jitter?
  max_size:     what dot size (mm) does frac = 1 get?
  test, p_adjust: which test and which adjustment for the box-plot brackets
                (fig.stats)? Both or neither. Python passes them to
                plotnine_extra.stat_pwc: test "wilcox.test" or "t.test";
                p_adjust "holm", "bonferroni", "BH", "fdr", "BY", "hommel" or
                "none". plots.R takes the ggpubr names ("wilcox_test", "t_test").
"""

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from mizani.bounds import squish
from mizani.labels import label_percent
from plotnine import (
    aes,
    after_stat,
    annotate,
    arrow,
    coord_equal,
    element_blank,
    element_text,
    geom_boxplot,
    geom_col,
    geom_hline,
    geom_line,
    geom_point,
    geom_segment,
    geom_smooth,
    geom_text,
    geom_tile,
    geom_vline,
    ggplot,
    guide_legend,
    guides,
    labs,
    position_jitter,
    scale_colour_cmap,
    scale_colour_manual,
    scale_fill_gradient2,
    scale_fill_manual,
    scale_size_area,
    scale_x_discrete,
    scale_y_continuous,
    scale_y_discrete,
    stat_summary,
    theme,
)

from palettes import palette_for

# Matplotlib colormaps for skills that still draw their own matplotlib panels.
diverging = LinearSegmentedColormap.from_list(
    "diverging", ["#2166AC", "white", "#B2182B"]
)
grey_to_red = LinearSegmentedColormap.from_list("grey_to_red", ["#E5E5E5", "#B2182B"])
# fig.wilke: under this many points per group a box misleads; draw points + median.
FEW_POINTS = 10
# fig.volcano: blue down, red up, grey not significant.
DIRECTION = {"down": "#2166AC", "up": "#B2182B", "ns": "#999999"}
# fig.heatmap: rotate column labels past this many columns.
ROTATE_PAST = 8
ROTATE_X = theme(axis_text_x=element_text(rotation=90, ha="center", va="top"))
# fig.embedding: no frame, ticks, tick labels or axis titles; corner arrows instead.
NO_AXES = theme(
    panel_border=element_blank(),
    panel_grid=element_blank(),
    axis_text=element_blank(),
    axis_ticks=element_blank(),
    axis_title=element_blank(),
)


def check_cols(data, cols, fn):
    missing = [col for col in cols if col not in data.columns]
    if missing:
        raise KeyError(f"{fn}: missing columns: {', '.join(missing)}")


def levels_of(values, pal=None):
    """Palette key order if pal is given, else Categorical order, else first seen."""
    if pal is not None:
        palette_for(values.unique(), pal)
        return [key for key in pal if key in set(values.astype(str))]
    if isinstance(values.dtype, pd.CategoricalDtype):
        return list(values.cat.categories)
    return list(values.unique())


def as_factor(values, levels):
    return pd.Categorical(values.astype(str), categories=[str(lv) for lv in levels])


def n_labels(groups):
    """fig.n: each group label over its count, `ctrl\\n(n = 4)`."""
    n = groups.value_counts()
    return lambda gs: [f"{g}\n(n = {n[g]})" for g in gs]


def corner_arrows(df, axis_name):
    """fig.embedding: two short arrows at the bottom-left, "UMAP 1" and "UMAP 2"."""
    x0, y0 = df["x"].min(), df["y"].min()
    dx, dy = 0.2 * (df["x"].max() - x0), 0.2 * (df["y"].max() - y0)
    head = arrow(length=0.05, type="closed")
    return [
        annotate("segment", x=x0, y=y0, xend=x0 + dx, yend=y0, arrow=head),
        annotate("segment", x=x0, y=y0, xend=x0, yend=y0 + dy, arrow=head),
        annotate("text", x=x0 + dx / 2, y=y0 - dy / 4, label=f"{axis_name} 1", size=7),
        annotate(
            "text",
            x=x0 - dx / 4,
            y=y0 + dy / 2,
            label=f"{axis_name} 2",
            angle=90,
            size=7,
        ),
    ]


def plot_volcano(
    data, lfc, padj, *, padj_max, lfc_min, label=None, n_label=None, title=None
):
    """Significant = padj < padj_max and |lfc| > lfc_min, coloured by direction.

    With `label`, significant genes are labelled: all of them, or the top
    `n_label` per direction by adjusted p.
    """
    check_cols(data, [lfc, padj] + ([label] if label else []), "plot_volcano")
    significant = (data[padj] < padj_max) & (data[lfc].abs() > lfc_min)
    df = pd.DataFrame(
        {
            "lfc": data[lfc],
            "padj": data[padj],
            "neg_log_p": -np.log10(data[padj]),
            "direction": np.select(
                [significant & (data[lfc] > 0), significant & (data[lfc] < 0)],
                ["up", "down"],
                "ns",
            ),
        }
    )
    p = (
        ggplot(df, aes("lfc", "neg_log_p", colour="direction"))
        + geom_point(size=0.5)
        + geom_hline(yintercept=-np.log10(padj_max), linetype="dashed")
        + geom_vline(xintercept=[-lfc_min, lfc_min], linetype="dashed")
        + scale_colour_manual(values=DIRECTION, guide=None)
        + labs(x=r"$\log_2$ fold change", y=r"$-\log_{10}$ adjusted p", title=title)
    )
    if label is None:
        return p
    df["label"] = data[label]
    shown = df[df["direction"] != "ns"]
    if n_label is not None:
        shown = shown.sort_values("padj").groupby("direction").head(n_label)
    return p + geom_text(aes(label="label"), data=shown, va="bottom", size=6)


def plot_box(
    data,
    group,
    value,
    pal,
    *,
    jitter_width,
    seed,
    value_label,
    pair=None,
    test=None,
    p_adjust=None,
    title=None,
):
    """Points per group over a box; below FEW_POINTS in any group, a median tick only.

    `value_label` is the y-axis title with its unit (fig.labels), never a column
    name. `test` and `p_adjust` add adjusted-p brackets (fig.stats); with `pair`
    the test is paired.
    """
    check_cols(data, [group, value] + ([pair] if pair else []), "plot_box")
    if (test is None) != (p_adjust is None):
        raise ValueError("plot_box: pass test and p_adjust together, or neither")
    levels = levels_of(data[group], pal)
    df = pd.DataFrame({"group": as_factor(data[group], levels), "value": data[value]})
    if pair is not None:
        # stat_pwc pairs values by row order within each group: sort by pair id.
        df = df.assign(pair=data[pair]).sort_values("pair")
    p = (
        ggplot(df, aes("group", "value"))
        + scale_fill_manual(values=pal, guide=None)
        + scale_colour_manual(values=pal, guide=None)
        + scale_x_discrete(labels=n_labels(df["group"]))
        + labs(x="", y=value_label, title=title)
    )
    # fig.wilke: under FEW_POINTS per group a box misleads; draw points and a median.
    if df["group"].value_counts().min() < FEW_POINTS:
        p += stat_summary(
            aes(colour="group"),
            fun_y=np.median,
            fun_ymin=np.median,
            fun_ymax=np.median,
            geom="crossbar",
            width=0.5,
        )
    else:
        p += geom_boxplot(aes(fill="group"), outlier_shape="")
    if pair is None:
        p += geom_point(
            position=position_jitter(width=jitter_width, height=0, random_state=seed)
        )
    else:
        # No jitter with pairs, so each line meets its points.
        p += geom_line(aes(group="pair"), colour="#999999")
        p += geom_point()
    if test is None:
        return p
    from plotnine_extra import stat_pwc

    if pair is not None and df.groupby("group", observed=True).size().nunique() > 1:
        raise ValueError("plot_box: a paired test needs every pair in every group")
    return (
        p
        + stat_pwc(
            aes(label=after_stat("'p.adj = ' + p_adj.map('{:.2g}'.format)")),
            method=test,
            p_adjust_method=p_adjust,
            paired=pair is not None,
            label_size=7,
        )
        + scale_y_continuous(expand=(0.05, 0, 0.15, 0))
    )


def plot_embedding(data, x, y, colour, pal, *, axis_name, title=None):
    """`axis_name` names the embedding ("UMAP"); arrows read "UMAP 1" and "UMAP 2"."""
    check_cols(data, [x, y, colour], "plot_embedding")
    df = pd.DataFrame(
        {
            "x": data[x],
            "y": data[y],
            "colour": as_factor(data[colour], levels_of(data[colour], pal)),
        }
    )
    return (
        ggplot(df, aes("x", "y", colour="colour"))
        + geom_point(size=0.3, raster=True)
        + scale_colour_manual(values=pal)
        + coord_equal()
        + corner_arrows(df, axis_name)
        + labs(colour=colour, title=title)
        + NO_AXES
    )


def plot_feature(data, x, y, value, *, axis_name, value_label, title=None):
    """fig.continuous: viridis colour, highest values drawn last (on top)."""
    check_cols(data, [x, y, value], "plot_feature")
    df = pd.DataFrame({"x": data[x], "y": data[y], "value": data[value]})
    df = df.sort_values("value")
    return (
        ggplot(df, aes("x", "y", colour="value"))
        + geom_point(size=0.3, raster=True)
        + scale_colour_cmap("viridis")
        + coord_equal()
        + corner_arrows(df, axis_name)
        + labs(colour=value_label, title=title)
        + NO_AXES
    )


def plot_heatmap(
    data,
    row,
    col,
    value,
    *,
    limit,
    col_annotation=None,
    annotation_pal=None,
    title=None,
    x_label=None,
    y_label=None,
    fill_label=None,
):
    """Diverging fill centred at 0, squished at +-limit; first row on top.

    `col_annotation` names a column of `data` that maps each column to a group,
    drawn as a strip above the heatmap in `annotation_pal` colours.
    """
    annotated = [col_annotation] if col_annotation else []
    check_cols(data, [row, col, value] + annotated, "plot_heatmap")
    rows, cols = levels_of(data[row]), levels_of(data[col])
    df = pd.DataFrame(
        {
            "row": as_factor(data[row], rows),
            "col": as_factor(data[col], cols),
            "value": data[value],
        }
    )
    p = (
        ggplot(df, aes("col", "row", fill="value"))
        + geom_tile()
        + scale_fill_gradient2(
            low="#2166AC",
            mid="white",
            high="#B2182B",
            midpoint=0,
            limits=(-limit, limit),
            oob=squish,
        )
        + scale_y_discrete(limits=[str(r) for r in rows[::-1]])
        + labs(
            x=x_label or col,
            y=y_label or row,
            fill=fill_label or value,
            title=title,
        )
    )
    if len(cols) > ROTATE_PAST:
        p += ROTATE_X
    if col_annotation is None:
        return p
    strip = data[[col, col_annotation]].drop_duplicates(col).set_index(col)
    groups = strip[col_annotation].reindex(cols)
    palette_for(groups.unique(), annotation_pal)
    ann = pd.DataFrame(
        {
            "x": np.arange(1, len(cols) + 1),
            "group": as_factor(groups, levels_of(groups, annotation_pal)),
        }
    )
    # fig.heatmap: one butt-ended segment per column reads as a solid strip.
    return (
        p
        + geom_segment(
            aes(x="x - 0.5", xend="x + 0.5", y="y", yend="y", colour="group"),
            data=ann.assign(y=len(rows) + 1),
            size=3,
            lineend="butt",
            inherit_aes=False,
        )
        + scale_colour_manual(values=annotation_pal)
        + labs(colour=col_annotation)
    )


def plot_dotplot(data, feature, group, mean, frac, *, max_size, mean_label, title=None):
    """fig.dotplot: size = fraction expressing (0-1), viridis colour = mean.

    `mean_label` titles the colourbar with its unit ("Mean expression (log1p
    CP10K)"); the size legend is always "Fraction expressing". First group on top.
    """
    check_cols(data, [feature, group, mean, frac], "plot_dotplot")
    groups = levels_of(data[group])
    df = pd.DataFrame(
        {
            "feature": as_factor(data[feature], levels_of(data[feature])),
            "group": as_factor(data[group], groups),
            "mean": data[mean],
            "frac": data[frac],
        }
    )
    return (
        ggplot(df, aes("feature", "group"))
        + geom_point(aes(size="frac", colour="mean"))
        + scale_size_area(
            max_size=max_size,
            limits=(0, 1),
            breaks=[0.25, 0.5, 0.75, 1],
            labels=label_percent(),
        )
        + scale_colour_cmap("viridis")
        + scale_y_discrete(limits=[str(g) for g in groups[::-1]])
        + guides(size=guide_legend(override_aes={"colour": "#999999"}))
        + labs(
            x=feature,
            y=group,
            size="Fraction expressing",
            colour=mean_label,
            title=title,
        )
        + ROTATE_X
    )


def plot_composition(
    data,
    sample,
    label,
    count,
    pal,
    *,
    title=None,
    x_label=None,
    y_label=None,
    fill_label=None,
):
    """Stacked bar of each label's share of the count per sample."""
    check_cols(data, [sample, label, count], "plot_composition")
    df = pd.DataFrame(
        {
            "sample": as_factor(data[sample], levels_of(data[sample])),
            "label": as_factor(data[label], levels_of(data[label], pal)),
            "count": data[count],
        }
    )
    return (
        ggplot(df, aes("sample", "count", fill="label"))
        + geom_col(position="fill")
        + scale_fill_manual(values=pal)
        + scale_y_continuous(labels=label_percent())
        + labs(
            x=x_label or sample,
            y=y_label or f"% of {count}",
            fill=fill_label or label,
            title=title,
        )
    )


def plot_scatter_fit(data, x, y, label=None, *, title=None, x_label=None, y_label=None):
    """Points with a linear fit and its 95% confidence band."""
    check_cols(data, [x, y] + ([label] if label else []), "plot_scatter_fit")
    df = pd.DataFrame({"x": data[x], "y": data[y]})
    p = (
        ggplot(df, aes("x", "y"))
        + geom_point(size=0.5)
        + geom_smooth(method="lm", se=True, level=0.95)
        + labs(x=x_label or x, y=y_label or y, title=title)
    )
    if label is None:
        return p
    df["label"] = data[label]
    return p + geom_text(aes(label="label"), va="bottom", size=6)


if __name__ == "__main__":
    import tempfile
    from pathlib import Path

    from palettes import pal_cell_type, pal_condition
    from theme import save_panel

    rng = np.random.default_rng(1)
    genes = pd.DataFrame(
        {
            "gene": [f"G{i}" for i in range(20)],
            "log2_fc": rng.normal(0, 2, 20),
            "padj": rng.uniform(0.0001, 1, 20),
        }
    )
    donors = pd.DataFrame(
        {
            "donor": [f"D{i}" for i in range(6)] * 2,
            "condition": ["healthy"] * 6 + ["disease"] * 6,
            "score": rng.normal(0, 1, 12),
        }
    )
    cells = pd.DataFrame(
        {
            "umap_1": rng.normal(0, 1, 60),
            "umap_2": rng.normal(0, 1, 60),
            "cell_type": ["t_cell", "b_cell", "monocyte"] * 20,
            "cd3e": rng.gamma(1, 1, 60),
        }
    )
    grid = pd.DataFrame(
        {
            "gene": ["CD3E", "MS4A1", "CD14"] * 3,
            "cell_type": np.repeat(["t_cell", "b_cell", "monocyte"], 3),
            "z": rng.normal(0, 3, 9),
            "frac": rng.uniform(0, 1, 9),
        }
    )
    # Ten samples, so the heatmap rotates its labels and draws the strip.
    wide = pd.DataFrame(
        {
            "gene": np.repeat(["CD3E", "MS4A1"], 10),
            "sample": [f"S{i:02d}" for i in range(10)] * 2,
            "condition": (["healthy"] * 5 + ["disease"] * 5) * 2,
            "z": rng.normal(0, 1, 20),
        }
    )
    counts = pd.DataFrame(
        {
            "sample": np.repeat(["S1", "S2"], 3),
            "cell_type": ["t_cell", "b_cell", "monocyte"] * 2,
            "n_cells": [50, 30, 20, 10, 60, 30],
        }
    )
    box = {
        "pal": pal_condition,
        "seed": 1,
        "jitter_width": 0.2,
        "value_label": "score (a.u.)",
    }
    figures = {
        "volcano": (
            plot_volcano(
                genes, "log2_fc", "padj", padj_max=0.1, lfc_min=0.5, label="gene"
            ),
            genes,
        ),
        "volcano-top": (
            plot_volcano(
                genes,
                "log2_fc",
                "padj",
                padj_max=0.5,
                lfc_min=0.5,
                label="gene",
                n_label=2,
                title="Disease vs healthy",
            ),
            genes,
        ),
        "box": (plot_box(donors, "condition", "score", **box), donors),
        "box-test": (
            plot_box(
                donors,
                "condition",
                "score",
                test="wilcox.test",
                p_adjust="BH",
                title="CD14",
                **box,
            ),
            donors,
        ),
        "box-paired": (
            plot_box(
                donors,
                "condition",
                "score",
                pair="donor",
                test="t.test",
                p_adjust="holm",
                **box,
            ),
            donors,
        ),
        "box_many": (
            plot_box(
                cells,
                "cell_type",
                "cd3e",
                pal_cell_type,
                jitter_width=0.2,
                seed=1,
                value_label="CD3E (log1p)",
            ),
            cells,
        ),
        "embedding": (
            plot_embedding(
                cells, "umap_1", "umap_2", "cell_type", pal_cell_type, axis_name="UMAP"
            ),
            cells,
        ),
        "feature": (
            plot_feature(
                cells,
                "umap_1",
                "umap_2",
                "cd3e",
                axis_name="UMAP",
                value_label="CD3E (log1p CP10K)",
                title="CD3E",
            ),
            cells,
        ),
        "heatmap": (plot_heatmap(grid, "cell_type", "gene", "z", limit=2), grid),
        "heatmap-annotated": (
            plot_heatmap(
                wide,
                "gene",
                "sample",
                "z",
                limit=2,
                col_annotation="condition",
                annotation_pal=pal_condition,
            ),
            wide,
        ),
        "dotplot": (
            plot_dotplot(
                grid,
                "gene",
                "cell_type",
                "z",
                "frac",
                max_size=4,
                mean_label="Mean expression (z-score)",
            ),
            grid,
        ),
        "composition": (
            plot_composition(counts, "sample", "cell_type", "n_cells", pal_cell_type),
            counts,
        ),
        "scatter-fit": (plot_scatter_fit(genes, "log2_fc", "padj"), genes),
    }
    with tempfile.TemporaryDirectory() as run_dir:
        for name, (plot, data) in figures.items():
            save_panel(plot, data, run_dir, name, 89, 60)
            assert (Path(run_dir) / f"{name}.pdf").exists(), name
    for bad_call in (
        lambda: plot_volcano(genes, "lfc", "padj", padj_max=0.1, lfc_min=0.5),
        lambda: plot_box(donors, "condition", "score", **{**box, "pal": pal_cell_type}),
        lambda: plot_box(donors, "condition", "score", test="t.test", **box),
    ):
        try:
            bad_call()
        except (KeyError, ValueError):
            pass
        else:
            raise SystemExit("a missing column, label or p_adjust did not raise")
