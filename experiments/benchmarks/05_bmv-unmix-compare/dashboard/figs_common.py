# Purpose: list every figure the common-gate INDEX cites, resolved to a file that exists. Shared by the
#          artifact builder (artifact/build_common.py) and the Quarto dashboard (96c_dashboard.qmd).
#          A cited run-dir figure is swapped for its restyled copy when MANIFEST_restyle_common.csv has one;
#          otherwise the cited path is used. A figure that resolves to no file is an error.
# Run:     python3 experiments/benchmarks/05_bmv-unmix-compare/dashboard/figs_common.py OUT.csv   (repo root)
# Inputs:  experiments/runs/2026-11-11_bmv-summary-common/INDEX.md              cited figures, in text order
#          experiments/runs/2026-11-11_bmv-restyle-common/MANIFEST_restyle_common.csv   old -> new paths
# Outputs: OUT.csv  section_idx, heading, line_no, token, cited, path (relative to experiments/runs),
#                   kind (restyled, cited, run-dir), caption

import csv
import fnmatch
import re
import sys
from pathlib import Path, PurePosixPath

RUNS = Path("experiments/runs")
INDEX_MD = RUNS / "2026-11-11_bmv-summary-common/INDEX.md"
MANIFEST = RUNS / "2026-11-11_bmv-restyle-common/MANIFEST_restyle_common.csv"
RESTYLE = "2026-11-11_bmv-restyle-common/"
CAPTION_MAX = 240  # characters of the citing sentence kept as the figure caption

# The INDEX names these without a run folder (§02); the figures that match the common labels are in the clean3 gate run (G10b).
OVERRIDES = {"figures/strategy_<donor>.png": "2026-11-11_bmv-common-gates-clean3/figures/strategy_<donor>.png"}

PREFIX = "experiments/runs/"
strip = lambda s: s[len(PREFIX) :] if s.startswith(PREFIX) else s
NA = {"", "NA"}


def load_manifest():
    by_old, by_src = {}, {}
    with MANIFEST.open(newline="") as f:
        for r in csv.DictReader(f):
            if r["new_figure"] in NA:
                continue
            new = strip(r["new_figure"])
            if r["old_figure"] not in NA:
                by_old[r["old_figure"]] = new
            for src in re.split(r"[;,]\s*", r["source"]):
                if src not in NA:
                    by_src[(str(PurePosixPath(strip(src)).parent), PurePosixPath(new).name)] = new
    return by_old, by_src


BY_OLD, BY_SRC = load_manifest()


def resolve_one(rel):
    """Return (path, kind) for one cited path without wildcards, or None when no file exists."""
    p = PurePosixPath(rel)
    for cand in (BY_OLD.get(rel), BY_SRC.get((str(p.parent), p.name)), rel):
        if cand and (RUNS / cand).is_file():
            return cand, ("restyled" if cand != rel else "cited" if rel.startswith(RESTYLE) else "run-dir")
    return None


def resolve(rel):
    """All (path, kind) matches of a cited path; a <placeholder> becomes a wildcard."""
    if "<" not in rel and "*" not in rel:
        hit = resolve_one(rel)
        return [hit] if hit else []
    p = PurePosixPath(re.sub(r"<[^>]+>", "*", rel))  # '*' in a cited path is a wildcard too
    names = {f.name for f in (RUNS / p.parent).glob(p.name)}
    # restyled copies of a run-dir folder: manifest names whose source sits in this folder
    names |= {n for (d, n) in BY_SRC if d == str(p.parent) and fnmatch.fnmatch(n, p.name)}
    hits = {resolve_one(f"{p.parent}/{n}") for n in sorted(names)}
    return sorted(h for h in hits if h)


def trim_sentence(line, pos):
    """The sentence of the line that holds the citation at pos; a '- `file`: text' bullet gives its text."""
    bullet = re.match(r"^\s*- `[^`]+\.png[^`]*`(?:, `[^`]+`)*:\s*(.*)$", line)
    if bullet and pos < line.index(":") + 1:
        text = bullet.group(1)
    else:
        cut = line.rfind(". ", 0, pos)
        start = cut + 2 if cut >= 0 else 0  # no earlier sentence: start at the line start (not at character 1)
        end = line.find(". ", pos)
        text = line[start : end + 1 if end >= 0 else len(line)]
    text = re.sub(r"^\s*-\s+", "", text)

    def drop_or_name(m):  # a trailing citation is dropped (the path shows under the figure); a path used as a subject keeps its file name
        rest = text[m.end():]
        if re.match(r"\s*([.;:,)—]|$)", rest):
            return ""
        return " " + m.group(1).rsplit("/", 1)[-1]

    text = re.sub(r"\s*[—(]?\s*`([^`]*\.png[^`]*)`\)?", drop_or_name, text).strip()
    text = re.sub(r"(?:\s*[,;]){2,}", ",", text)  # several dropped paths in a list leave stacked separators
    text = re.sub(r"[`*]", "", text).strip(" —;,:")  # a dropped path leaves a dangling comma or colon
    text = re.sub(r"[,;:]+\.$", ".", text)  # ",." where the dropped path sat before the sentence end
    if text.count(")") > text.count("("):  # sentence split inside a parenthesis: keep the part after it closes
        m = re.search(r"\)\.\s+(\S.*)$", text)
        if m:
            text = m.group(1)
    if text.count("(") > text.count(")"):  # the dropped path closed a parenthesis opened earlier in the sentence
        text = re.sub(r"[;,:]+\.\s", "). ", text, count=1)
    if text.count("(") > text.count(")"):
        text = text.rstrip(".;:, ") + ")."
    if re.match(r"[a-z]+\b(?![_.-])", text):  # a '- `file`: text' bullet gives a lower-case start; identifiers keep their case
        text = text[0].upper() + text[1:]
    if text and text[-1] not in ".!?)":
        text += "."
    return text


TITLES = {}  # figure stem -> title, from report_figures/captions.csv (the fallback when no sentence of the INDEX fits)
_cap_file = RUNS / "2026-11-11_bmv-restyle-common/report_figures/captions.csv"
if _cap_file.exists():
    with _cap_file.open(newline="") as _f:
        TITLES = {r["key"]: r["caption"].strip() for r in csv.DictReader(_f)}

BAD_CAPTION = re.compile(r",,|:,|;,|,\.|;\.|\(\s|\s\)|…|\b(and|are|in|of|the|to|from|is)[,.]$|^[,;:)]|\(\)|\bare\s*\.|:\s*[;.]")


def caption_ok(text, limit=CAPTION_MAX):
    """A caption is one whole sentence: short, balanced, ended, with no dropped-path residue."""
    return (
        0 < len(text) <= limit
        and text.count("(") == text.count(")")
        and text[-1] in ".!?)"
        and re.match(r"[A-Z0-9]|[a-z0-9]+[_.-]", text) is not None
        and not BAD_CAPTION.search(text)
    )


def caption_of(line, pos, token=""):
    """Whole-sentence caption: the citing sentence; else the first sentence of the line naming the figure; else its title."""
    text = trim_sentence(line, pos)
    if caption_ok(text):
        return text
    name = token.rsplit("/", 1)[-1]
    for sent in re.split(r"(?<=[.!?])\s+(?=[A-Z`*])", line):
        if name and name in sent:
            at = line.find(sent)
            cand = trim_sentence(line, at + sent.find(name))
            if caption_ok(cand):
                return cand
    stem = name[:-4] if name.endswith(".png") else name
    if stem in TITLES and caption_ok(TITLES[stem], 400):
        return TITLES[stem]
    return stem.replace("_", " ").capitalize() + "."


def collect():
    lines = INDEX_MD.read_text().split("\n")
    sections = []  # (heading, first line, last line), line numbers are 0-based indices into lines
    for i, l in enumerate(lines):
        if l.startswith("## "):
            if sections:
                sections[-1][2] = i - 1
            sections.append([l[3:].strip(), i, len(lines) - 1])
    items, missing = [], []
    for s_idx, (heading, first, last) in enumerate(sections):
        ctx = []  # folders named earlier in this section, most recent first
        for i in range(first, last + 1):
            line = lines[i]
            toks = [(m.group(1), m.start()) for m in re.finditer(r"`([^`]+)`", line)]
            folders = []  # folders named on this line, from full paths and folder tokens
            for t, _ in toks:
                if t.startswith("20") and t.endswith("/"):
                    folders.append(t.rstrip("/"))
                elif t.startswith("20") and "/" in t and t.endswith(".png"):
                    folders.append(str(PurePosixPath(t).parent))
            for t, pos in toks:
                if not t.endswith(".png"):
                    continue
                t = OVERRIDES.get(t, t)
                full = t.startswith("20") and "/" in t
                dirs = [""] if full else folders + ctx
                hits = []
                for d in dirs:
                    hits = resolve(f"{d}/{t}" if d else t)
                    if hits:
                        break
                if not hits:
                    missing.append(f"line {i + 1} (section {heading[:30]}): {t}")
                    continue
                for path, kind in hits:
                    items.append([s_idx, heading, i, t, d + "/" + t if d else t, path, kind, caption_of(line, pos, t)])
            ctx = [d for d in dict.fromkeys(folders + ctx)]
    return sections, items, missing


def dedupe(items):
    """Keep the first citation of each file within a section."""
    seen, out = set(), []
    for it in items:
        if (it[0], it[5]) not in seen:
            seen.add((it[0], it[5]))
            out.append(it)
    return out


if __name__ == "__main__":
    sections, items, missing = collect()
    if missing:
        sys.exit("missing figures:\n" + "\n".join(missing))
    items = dedupe(items)
    with open(sys.argv[1], "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["section_idx", "heading", "line_no", "token", "cited", "path", "kind", "caption"])
        w.writerows(items)
    print(f"{len(items)} figures in {len(sections)} sections -> {sys.argv[1]}")
