# Purpose: build the hand-written artifact page of the BMV unmix comparison (user, 2026-11-11): text verbatim
#          from INDEX.md, figures only from the certified list (audit_freshness.csv, current = TRUE), in the
#          section order and figure order of 90_report.qmd. Images are the 300-dpi copies of the local
#          Quarto dashboard (user, 2026-11-11: "figures should be 300 dpi"); each figure's plotted data is written as JSON and loaded only when opened.
# Run:     pixi run python3 experiments/benchmarks/05_bmv-unmix-compare/dashboard/artifact/build.py   (repo root)
#          then publish OUT_DIR/index.html with img/ and data/ as Artifact files.
# Inputs:  experiments/runs/2026-11-10_bmv-summary/INDEX.md                 text, verbatim (via pandoc)
#          experiments/benchmarks/05_bmv-unmix-compare/90_report.qmd         section order, figure calls, captions
#          experiments/runs/2026-11-11_bmv-restyle/audit_freshness.csv       certified figures
#          experiments/runs/2026-11-11_bmv-dashboard-v6/img_log.csv, img/    300-dpi images [96_dashboard.qmd]
#          experiments/benchmarks/05_bmv-unmix-compare/dashboard/artifact/{page.html, cards.json}
# Outputs: OUT_DIR/index.html        the page
#          OUT_DIR/img/*.png         web images, one per figure
#          OUT_DIR/data/*.json       plotted data per figure (first MAX_ROWS rows)
#          OUT_DIR/figures.csv       section, figure, caption, image, data file

import csv
import html
import json
import re
import shutil
import subprocess
from pathlib import Path

# Parameters ---------------------------------------------------------------------------------------
RUNS = Path("experiments/runs")
HERE = Path("experiments/benchmarks/05_bmv-unmix-compare/dashboard/artifact")
INDEX_MD = RUNS / "2026-11-10_bmv-summary/INDEX.md"
REPORT_QMD = Path("experiments/benchmarks/05_bmv-unmix-compare/90_report.qmd")
AUDIT_CSV = RUNS / "2026-11-11_bmv-restyle/audit_freshness.csv"
WEB_DIR = RUNS / "2026-11-11_bmv-dashboard-v6"  # 300-dpi images, certified set (local Quarto dashboard v6)
FIG_DIR = "2026-11-10_bmv-report/figures-v7"  # s* redraws the certified list names (tables equal v8, v9)
KEPT_DIR = "2026-11-11_bmv-restyle"
OUT_DIR = RUNS / "2026-11-11_bmv-dashboard-artifact-v7"  # v7: §17 caption fix; v6 revised prose
MAX_ROWS = 500  # rows per data table; the full CSV path is named under the table
SKIP_SECTIONS = {"08"}  # 08 points at the superseded offline dashboard this page replaces

assert not OUT_DIR.exists(), f"output dir must not exist: {OUT_DIR}"

# Figures per report section, in report order ----------------------------------------------------------
report = REPORT_QMD.read_text()
titles = dict(re.findall(r'^\s*"(s\d\d[^"]*)"\s*,\s*"[^"]*"\s*,\s*"[^"]*"\s*,\s*"([^"]*)"', report, re.M))
titles["s06_spectra_angle"] = "Angle between each dye's 1.0.2 and 1.8.8 spectrum; each point is one dye."
call = re.compile(
    r'show\(\s*"([^"]+)"\s*\)'
    r'|show_kept\(\s*"([^"]+)"\s*,\s*"((?:[^"\\]|\\.)*)"\s*\)'
    r'|show_run\(\s*"([^"]+)"\s*,\s*"((?:[^"\\]|\\.)*)"\s*\)'
    # walk(c("IFNg", ...), \(x) show_kept(sprintf("...%s...", x), sprintf("...%s...", x))): split figures
    r'|walk\(\s*c\(([^)]*)\)\s*,\s*\\\(x\)\s*show_kept\(\s*sprintf\(\s*"([^"]+)"\s*,\s*x\s*\)\s*,'
    r'\s*sprintf\(\s*"((?:[^"\\]|\\.)*)"\s*,\s*x\s*\)\s*\)\s*\)',
    re.S,
)
report_sections = {}
for block in re.split(r"^## ", report, flags=re.M)[1:]:
    number = block[:2]
    if not number.isdigit():
        continue
    figs = []
    for m in call.finditer(block):
        if m.group(1):
            figs.append((f"{FIG_DIR}/{m.group(1)}.png", titles.get(m.group(1), m.group(1))))
        elif m.group(2):
            figs.append((f"{KEPT_DIR}/{m.group(2)}", m.group(3)))
        elif m.group(4):
            figs.append((m.group(4), m.group(5)))
        else:
            for x in re.findall(r'"([^"]+)"', m.group(6)):
                figs.append((f"{KEPT_DIR}/{m.group(7) % x}", m.group(8) % x))
    report_sections[number] = figs

with AUDIT_CSV.open() as f:
    certified = {r["figure"] for r in csv.DictReader(f) if r["current"] == "TRUE"}
with (WEB_DIR / "img_log.csv").open() as f:
    web = {r["rel"]: Path(r["web"]) for r in csv.DictReader(f)}
cited = [rel for figs in report_sections.values() for rel, _ in figs]
missing = [rel for rel in cited if rel not in certified or rel not in web]
assert not missing, f"figures not certified or without a web image: {missing}"

# INDEX text, converted once with pandoc ------------------------------------------------------------------
# Sections are split on "## " headings; a marker paragraph between them lets one pandoc run convert all.
index = INDEX_MD.read_text()
parts = re.split(r"^## ", index, flags=re.M)
intro, chunks = parts[0], parts[1:]
headings = [c.split("\n", 1)[0].strip() for c in chunks]
bodies = [c.split("\n", 1)[1] if "\n" in c else "" for c in chunks]
intro_body = intro.split("\n", 1)[1]  # drop the "# ..." page title
marker = "\n\nSWMARKER{}SWMARKER\n\n"
joined = intro_body + "".join(marker.format(i) + b for i, b in enumerate(bodies))
converted = subprocess.run(
    ["quarto", "pandoc", "-f", "markdown-smart", "-t", "html", "--wrap=none"],
    input=joined,
    capture_output=True,
    text=True,
    check=True,
).stdout
pieces = re.split(r"<p>SWMARKER(\d+)SWMARKER</p>", converted)
intro_html = pieces[0]
section_html = {int(pieces[i]): pieces[i + 1] for i in range(1, len(pieces), 2)}
assert len(section_html) == len(headings), "every INDEX section converted"


# Page parts ---------------------------------------------------------------------------------------------
def slug(text):
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")[:48]


def figure_html(rel, caption, n):
    stem = f"f{n:03d}"
    img_name = f"{stem}.png"
    shutil.copyfile(web[rel], OUT_DIR / "img" / img_name)
    src_csv = RUNS / rel.replace(".png", ".csv")
    data_note = "No plotted table was saved beside this figure; it is drawn from per-event data."
    data_attr = ""
    if src_csv.exists():
        with src_csv.open(newline="") as f:
            rows = list(csv.reader(f))
        header, body = rows[0], rows[1:]
        json.dump(
            {"columns": header, "rows": body[:MAX_ROWS], "n": len(body)},
            (OUT_DIR / "data" / f"{stem}.json").open("w"),
        )
        shown = min(len(body), MAX_ROWS)
        data_note = f"{shown} of {len(body)} rows · {src_csv}"
        data_attr = f' data-src="data/{stem}.json"'
    figures_log.append([current_section, rel, caption, f"img/{img_name}", data_attr and f"data/{stem}.json"])
    return f"""
<figure class="fig">
  <a href="img/{img_name}" target="_blank" rel="noopener"><img src="img/{img_name}" alt="{html.escape(caption)}" loading="lazy"></a>
  <figcaption>{html.escape(caption)} <code>{html.escape(rel)}</code></figcaption>
  <details class="data"{data_attr}>
    <summary>Plotted data</summary>
    <p class="data-note">{html.escape(data_note)}</p>
    <div class="table-wrap"></div>
  </details>
</figure>"""


(OUT_DIR / "img").mkdir(parents=True)
(OUT_DIR / "data").mkdir()
figures_log = []
nav, sections = [], []
n_fig = 0
for i, heading in enumerate(headings):
    number = heading[:2]
    if number in SKIP_SECTIONS:
        continue
    anchor = slug(heading)
    current_section = heading
    if number.isdigit():
        name = heading.split(":", 1)[1].strip() if ":" in heading else heading
        label = f'<span class="num">{number}</span>{html.escape(name)}'
        nav.append(f'<a href="#{anchor}"><span class="num">{number}</span>{html.escape(name.split(" (")[0])}</a>')
    else:
        label = html.escape(heading)
        nav.append(f'<a href="#{anchor}">{html.escape(heading)}</a>')
    figs_html = ""
    for rel, caption in report_sections.get(number, []):
        n_fig += 1
        figs_html += figure_html(rel, caption, n_fig)
    if figs_html:
        figs_html = f'<div class="figs">{figs_html}</div>'
    sections.append(
        f'<section id="{anchor}"><h2>{label}</h2><div class="prose">{section_html[i]}</div>{figs_html}</section>'
    )

cards = json.loads((HERE / "cards.json").read_text())
cards_html = "".join(
    f'<div class="card {c["state"]}"><div class="c-label">{html.escape(c["label"])}</div>'
    f'<div class="c-value">{html.escape(c["value"])}</div><div class="c-note">{html.escape(c["note"])}</div></div>'
    for c in cards
)

page = (HERE / "page.html").read_text()
page = (
    page.replace("{{CARDS}}", cards_html)
    .replace("{{NAV}}", "".join(nav))
    .replace("{{INTRO}}", intro_html)
    .replace("{{SECTIONS}}", "".join(sections))
    .replace("{{N_FIG}}", str(n_fig))
)
(OUT_DIR / "index.html").write_text(page)
with (OUT_DIR / "figures.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["section", "figure", "caption", "image", "data"])
    w.writerows(figures_log)
print(f"{n_fig} figures, {len(sections)} sections -> {OUT_DIR}")
