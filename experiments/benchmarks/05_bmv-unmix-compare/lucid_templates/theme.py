"""Figure theme and panel export for every figure script; theme.R is the twin."""

from pathlib import Path

import matplotlib
from plotnine import ggplot, theme, theme_classic, theme_set

# Type 42 embeds editable TrueType text instead of outlines. plotnine has no
# theme_pubr; theme_classic at 8 pt with the legend on top is its nearest twin
# (R: theme_pubr(base_size = 8)).
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8})
theme_set(theme_classic(base_size=8) + theme(legend_position="top"))


def save_panel(plot, data, run_dir, name, width_mm, height_mm):
    """Write <run_dir>/<name>.pdf (vector), .png (review) and .csv (source data).

    `plot` is a plotnine ggplot or a matplotlib Figure. `data` is the same frame
    the panel was drawn from, so the table and the panel cannot drift apart.
    """
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(
            f"save_panel: run_dir {run_dir} does not exist; create it before saving."
        )
    if isinstance(plot, ggplot):
        size = {"width": width_mm, "height": height_mm, "units": "mm", "verbose": False}
        plot.save(run_dir / f"{name}.pdf", **size)
        plot.save(run_dir / f"{name}.png", dpi=300, **size)
    else:
        plot.set_size_inches(width_mm / 25.4, height_mm / 25.4)
        plot.savefig(run_dir / f"{name}.pdf")
        plot.savefig(run_dir / f"{name}.png", dpi=300)
    data.to_csv(run_dir / f"{name}.csv", index=False)
    return run_dir / name
