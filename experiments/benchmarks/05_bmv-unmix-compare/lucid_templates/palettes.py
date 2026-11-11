"""Project palettes and figure choices the house style leaves open.

Fill each field once per project; every figure script reads this file.
  dotplot defaults:    [fill in: size range, colour scale, dot.min]
  volcano thresholds:  [fill in; one seen pair is Q < 0.1 and |log2FC| > 0.5]
  heatmap anchors:     [fill in; one seen example is -4 / 0 / 4]
  font family:         [fill in: e.g. Helvetica]

Palettes are dicts keyed by biological label, never by position.
The two below are placeholders; the project replaces them with its own labels.
"""

pal_condition = {"healthy": "#4477AA", "disease": "#EE6677"}
pal_cell_type = {"t_cell": "#228833", "b_cell": "#CCBB44", "monocyte": "#AA3377"}

# Colour-blind-safe categorical sets (fig.colourblind). Hex values from the
# sources: Paul Tol (tol-colors 2.2.0 colors.json), Okabe-Ito (R
# grDevices::palette.colors), Tableau "Color Blind" (ggthemes); standard
# Tableau 10 is not safe (red against green) and is not listed.
TOL_BRIGHT = [
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
]
OKABE_ITO = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000",
]
TOL_MUTED = [
    "#CC6677",
    "#332288",
    "#DDCC77",
    "#117733",
    "#88CCEE",
    "#882255",
    "#44AA99",
    "#999933",
    "#AA4499",
]
TABLEAU_CB10 = [
    "#1170AA",
    "#FC7D0B",
    "#A3ACB9",
    "#57606C",
    "#5FA2CE",
    "#C85200",
    "#7B848F",
    "#A3CCE9",
    "#FFBC79",
    "#C8D0D9",
]
# Tol discrete rainbow: his 29-colour ramp and, for each n (11-23), the
# 0-based indexes of the colours to use.
TOL_RAINBOW_HEX = [
    "#E8ECFB",
    "#D9CCE3",
    "#D1BBD7",
    "#CAACCB",
    "#BA8DB4",
    "#AE76A3",
    "#AA6F9E",
    "#994F88",
    "#882E72",
    "#1965B0",
    "#437DBF",
    "#5289C7",
    "#6195CF",
    "#7BAFDE",
    "#4EB265",
    "#90C987",
    "#CAE0AB",
    "#F7F056",
    "#F7CB45",
    "#F6C141",
    "#F4A736",
    "#F1932D",
    "#EE8026",
    "#E8601C",
    "#E65518",
    "#DC050C",
    "#A5170E",
    "#72190E",
    "#42150A",
]
TOL_RAINBOW_INDEX = {
    11: [8, 9, 11, 13, 14, 16, 17, 20, 23, 25, 27],
    12: [2, 5, 8, 9, 11, 13, 14, 16, 17, 20, 23, 25],
    13: [2, 5, 8, 9, 11, 13, 14, 15, 16, 17, 20, 23, 25],
    14: [2, 5, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25],
    15: [2, 5, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 27],
    16: [2, 4, 6, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 27],
    17: [2, 4, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 27],
    18: [2, 4, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 26, 27],
    19: [1, 3, 4, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 26, 27],
    20: [1, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 19, 21, 23, 25, 26, 27],
    21: [1, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24, 25, 26, 27],
    22: [
        1,
        3,
        4,
        6,
        7,
        8,
        9,
        10,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        20,
        22,
        24,
        25,
        26,
        27,
        28,
    ],
    23: [
        0,
        1,
        3,
        4,
        6,
        7,
        8,
        9,
        10,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        20,
        22,
        24,
        25,
        26,
        27,
        28,
    ],
}


def categorical_palette(labels):
    """A colour-blind-safe dict keyed by label, from the smallest set that fits.

    Up to 7 Tol bright, 8 Okabe-Ito, 9 Tol muted, 10 Tableau Color Blind,
    11-23 Tol discrete rainbow; past 23 no set stays distinguishable, so it raises.
    """
    labels = list(dict.fromkeys(str(label) for label in labels))
    n = len(labels)
    for colours in (TOL_BRIGHT, OKABE_ITO, TOL_MUTED, TABLEAU_CB10):
        if n <= len(colours):
            return dict(zip(labels, colours, strict=False))
    if n in TOL_RAINBOW_INDEX:
        return {
            lab: TOL_RAINBOW_HEX[i]
            for lab, i in zip(labels, TOL_RAINBOW_INDEX[n], strict=True)
        }
    raise ValueError(
        f"categorical_palette: {n} categories; no colour set stays colour-blind safe "
        "past 23. Facet, label directly or lump rare categories."
    )


def text_colour(fill):
    """Black or white text for a hex fill, whichever contrasts more (WCAG luminance)."""
    rgb = [int(fill.lstrip("#")[i : i + 2], 16) / 255 for i in (0, 2, 4)]
    lin = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
    lum = 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    return "#000000" if (lum + 0.05) / 0.05 >= 1.05 / (lum + 0.05) else "#FFFFFF"


def palette_for(labels, pal):
    """Colours for a sequence of labels; an unknown label raises.

    Failing here means a mislabelled or renamed group can never fall through
    as an unplanned colour.
    """
    missing = sorted({str(label) for label in labels} - pal.keys())
    if missing:
        raise KeyError(
            f"palette_for: labels not in the palette: {missing}; "
            "add them to palettes.py or fix the labels."
        )
    return [pal[str(label)] for label in labels]


if __name__ == "__main__":
    assert palette_for(["disease"], pal_condition) == ["#EE6677"]
    try:
        palette_for(["unknown"], pal_condition)
    except KeyError:
        pass
    else:
        raise SystemExit("palette_for accepted an unknown label")
    assert list(categorical_palette(["a", "b"]).values()) == TOL_BRIGHT[:2]
    assert len(set(categorical_palette(range(23)).values())) == 23
    assert all(len(v) == n for n, v in TOL_RAINBOW_INDEX.items())
    assert text_colour("#000000") == "#FFFFFF" and text_colour("#F0E442") == "#000000"
