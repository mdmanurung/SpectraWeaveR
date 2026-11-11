# Purpose: build the hand-written artifact page of the BMV unmix comparison on the common SpectroFlo gates
#          (rerun, cut arcsinh 1.5, 2026-11-11). Text verbatim from the new INDEX.md; every cited figure is
#          resolved through MANIFEST_restyle_common.csv (restyled copy first, cited run-dir path otherwise) by
#          dashboard/figs_common.py. Figures are the 300-dpi PNGs; the Conclusions bullets carry their figures.
#          A cited figure with no file stops the build. Nothing is published.
# Run:     pixi run python3 experiments/benchmarks/05_bmv-unmix-compare/dashboard/artifact/build_common.py   (repo root)
#          then publish OUT_DIR/index.html with img/ and data/ as Artifact files (by hand).
# Inputs:  experiments/runs/2026-11-11_bmv-summary-common/INDEX.md      text, verbatim (via pandoc)
#          experiments/runs/2026-11-11_bmv-summary-common/CHANGES.md    the "What changed" table
#          experiments/runs/2026-11-11_bmv-restyle-common/MANIFEST_restyle_common.csv   old -> new figure paths
#          <figure>.csv beside each cited PNG                            plotted data, loaded when opened
#          experiments/benchmarks/05_bmv-unmix-compare/dashboard/artifact/{page.html, cards_common.json}
# Outputs: OUT_DIR/index.html     the page
#          OUT_DIR/img/*.png      one 300-dpi PNG per figure
#          OUT_DIR/data/*.json    plotted data per figure (first MAX_ROWS rows)
#          OUT_DIR/figures.csv    section, figure, kind, caption, image, data file

import csv
import html
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path("experiments/benchmarks/05_bmv-unmix-compare/dashboard/artifact")
sys.path.insert(0, str(HERE.parent))
import figs_common as fc  # noqa: E402

# Parameters ---------------------------------------------------------------------------------------
RUNS = fc.RUNS
CHANGES_MD = RUNS / "2026-11-11_bmv-summary-common/CHANGES.md"
OUT_DIR = RUNS / os.environ.get("SW_OUT_NAME", "2026-11-11_bmv-dashboard-artifact-common")  # set SW_OUT_NAME to build into a new dir
MAX_ROWS = 500  # rows per data table; the full CSV path is named under the table
MAX_FILE_MB = 15  # Artifact limit per file
SKIP_SECTIONS = {"08"}  # 08 points at the superseded offline dashboard this page replaces
SPIKE_FOLDERS = ["18_spike_in/", "20_spike_in_panel2/", "22_spike_cd3neg/", "23_cospike/", "25_spike_passes/"]
BANNER = "Rerun on common SpectroFlo gates, cut arcsinh 1.5 (2026-11-11)"

assert not OUT_DIR.exists(), f"output dir must not exist: {OUT_DIR}"
manifest_text = fc.MANIFEST.read_text()
no_spike = [d for d in SPIKE_FOLDERS if f"{fc.RESTYLE}{d}" not in manifest_text]
assert not no_spike, f"manifest has no spike rows yet for {no_spike}"

# Figures cited by the INDEX, resolved ---------------------------------------------------------------------
sections, items, missing = fc.collect()
assert not missing, "cited figures with no file:\n" + "\n".join(missing)
items = fc.dedupe(items)

# INDEX text, converted once with pandoc ------------------------------------------------------------------
# Chunks: the intro, each section, and, in Conclusions, each bullet (so its figures sit under it).
lines = fc.INDEX_MD.read_text().split("\n")
chunks = []  # (section index or -1 for the intro, first line, last line)
chunks.append((-1, 1, sections[0][1] - 1))  # line 0 is the page title
for s_idx, (heading, first, last) in enumerate(sections):
    if heading == "Conclusions":
        starts = [i for i in range(first + 1, last + 1) if lines[i].startswith("- ")]
        chunks.append((s_idx, first + 1, starts[0] - 1))
        chunks += [(s_idx, a, b - 1) for a, b in zip(starts, starts[1:] + [last + 1])]
    else:
        chunks.append((s_idx, first + 1, last))
marker = "\n\nSWMARKER{}SWMARKER\n\n"
joined = "".join(marker.format(k) + "\n".join(lines[a : b + 1]) for k, (_, a, b) in enumerate(chunks))
converted = subprocess.run(
    ["quarto", "pandoc", "-f", "markdown-smart", "-t", "html", "--wrap=none"],
    input=joined, capture_output=True, text=True, check=True,
).stdout
pieces = re.split(r"<p>SWMARKER(\d+)SWMARKER</p>", converted)
chunk_html = {int(pieces[i]): pieces[i + 1] for i in range(1, len(pieces), 2)}
assert len(chunk_html) == len(chunks), "every INDEX chunk converted"


# Page parts ---------------------------------------------------------------------------------------------
def slug(text):
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")[:48]


figures_log, made = [], {}  # made: figure path -> (image name, data json name or "")


def figure_html(it, n):
    _, heading, _, token, _, path, kind, caption = it
    if path not in made:
        stem = f"f{len(made) + 1:03d}"
        src = RUNS / path
        assert src.stat().st_size <= MAX_FILE_MB * 1024**2, f"over {MAX_FILE_MB} MB: {src}"
        shutil.copyfile(src, OUT_DIR / "img" / f"{stem}.png")
        data = ""
        if src.with_suffix(".csv").exists():
            with src.with_suffix(".csv").open(newline="") as f:
                rows = list(csv.reader(f))
            json.dump({"columns": rows[0], "rows": rows[1:][:MAX_ROWS], "n": len(rows) - 1}, (OUT_DIR / "data" / f"{stem}.json").open("w"))
            data = f"data/{stem}.json"
        made[path] = (f"img/{stem}.png", data, len(rows) - 1 if data else 0)
    img, data, n_rows = made[path]
    note = f"{min(n_rows, MAX_ROWS)} of {n_rows} rows · {RUNS / Path(path).with_suffix('.csv')}" if data else "No plotted table was saved beside this figure; it is drawn from per-event data."
    figures_log.append([heading, path, kind, caption, img, data])
    return f"""
<figure class="fig">
  <a href="{img}" target="_blank" rel="noopener"><img src="{img}" alt="{html.escape(caption)}" loading="lazy"></a>
  <figcaption>{html.escape(caption)} <code>{html.escape(path)}</code></figcaption>
  <details class="data"{f' data-src="{data}"' if data else ''}>
    <summary>Plotted data</summary>
    <p class="data-note">{html.escape(note)}</p>
    <div class="table-wrap"></div>
  </details>
</figure>"""


def figures_of(s_idx, a, b):
    return "".join(figure_html(it, 0) for it in items if it[0] == s_idx and a <= it[2] <= b)


(OUT_DIR / "img").mkdir(parents=True)
(OUT_DIR / "data").mkdir()
nav, body = [], []
for s_idx, (heading, first, last) in enumerate(sections):
    number = heading[:2]
    if number in SKIP_SECTIONS:
        continue
    anchor = slug(heading)
    if number.isdigit():
        name = heading.split(":", 1)[1].strip() if ":" in heading else heading
        label = f'<span class="num">{number}</span>{html.escape(name)}'
        nav.append(f'<a href="#{anchor}"><span class="num">{number}</span>{html.escape(name.split(" (")[0])}</a>')
    else:
        label = html.escape(heading)
        nav.append(f'<a href="#{anchor}">{html.escape(heading)}</a>')
    parts = ""
    for k, (c_idx, a, b) in enumerate(chunks):
        if c_idx != s_idx:
            continue
        parts += chunk_html[k]
        if heading == "Conclusions":  # a bullet's figures sit directly under it
            figs = figures_of(s_idx, a, b)
            parts += f'<div class="figs">{figs}</div>' if figs else ""
    if heading != "Conclusions":
        figs = figures_of(s_idx, first, last)
        parts += f'<div class="figs">{figs}</div>' if figs else ""
    body.append(f'<section id="{anchor}"><h2>{label}</h2><div class="prose">{parts}</div></section>')

# What changed: CHANGES.md as a table, last tab of the rail
changes_html = subprocess.run(
    ["quarto", "pandoc", "-f", "markdown-smart", "-t", "html", "--wrap=none"],
    input="\n".join(CHANGES_MD.read_text().split("\n")[1:]), capture_output=True, text=True, check=True,
).stdout
nav.append('<a href="#what-changed">What changed</a>')
body.append(f'<section id="what-changed"><h2>What changed</h2><div class="prose">{changes_html}</div></section>')

cards = json.loads((HERE / "cards_common.json").read_text())
cards_html = "".join(
    f'<div class="card {c["state"]}"><div class="c-label">{html.escape(c["label"])}</div>'
    f'<div class="c-value">{html.escape(c["value"])}</div><div class="c-note">{html.escape(c["note"])}</div></div>'
    for c in cards
)

banner_css = """
.banner { background: var(--accent-soft); border-bottom: 1px solid var(--rule); }
.banner .wrap { padding-block: 10px; font-size: 14px; display: flex; gap: 8px 16px; flex-wrap: wrap; align-items: baseline; }
.banner strong { font-weight: 600; }
.banner a { font-weight: 600; }
</style>"""
banner_html = f'<div class="banner"><div class="wrap"><strong>{BANNER}</strong><a href="#what-changed">What changed: every old claim against its new value</a></div></div>'

page = (HERE / "page.html").read_text()
swaps = [
    ("</style>", banner_css),
    ("</header>", "</header>\n" + banner_html),
    ("Text from INDEX.md as of 2026-11-11. Figures are the certified current set at Nature spec (300-dpi images, shown scaled to the column); open a figure for the full-resolution image and expand \"Plotted data\" for the table behind it.",
     "Text from the common-gate INDEX.md as of 2026-11-11. Figures are 300-dpi images at Nature spec, shown scaled to the column; open one for full resolution and expand \"Plotted data\" for the table behind it."),
    ("Built by <code>dashboard/artifact/build.py</code> from INDEX.md, 90_report.qmd and audit_freshness.csv.",
     "Built by <code>dashboard/artifact/build_common.py</code> from the common-gate INDEX.md, CHANGES.md and MANIFEST_restyle_common.csv."),
]
for old, new in swaps:
    assert page.count(old) == 1, f"page.html no longer holds: {old[:50]}"
    page = page.replace(old, new)
page = (
    page.replace("{{CARDS}}", cards_html)
    .replace("{{NAV}}", "".join(nav))
    .replace("{{INTRO}}", chunk_html[0])
    .replace("{{SECTIONS}}", "".join(body))
    .replace("{{N_FIG}}", str(len(made)))
)
(OUT_DIR / "index.html").write_text(page)
with (OUT_DIR / "figures.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["section", "figure", "kind", "caption", "image", "data"])
    w.writerows(figures_log)
total = sum(p.stat().st_size for p in OUT_DIR.rglob("*") if p.is_file())
print(f"{len(made)} figures ({len(figures_log)} placements), {len(body)} sections, {total / 1024**2:.1f} MB -> {OUT_DIR}")
