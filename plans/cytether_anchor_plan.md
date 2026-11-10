# Plan: internal anchors in HEB (CD3 and clean negative modes), 2026-11-10

Source idea: Glehr et al., eBioMedicine 2026, CD3-relativisation (DOI 10.1016/j.ebiom.2026.106427, PMC13488037),
Box 1. The published method needs CD3 single stains in every channel; CytoBatchFlagR (Zenodo 15388817) has only
fully stained files, so the exact method is not reproducible here. Review (Fable, 2026-11-10; numbers measured on
T6, 10,000 cells per sample, scripts in the session scratchpad): CD3 = BUV563-A, clean 2-means split (4.5 to 5.3
pooled SD); CD3+ median SD over batches 0.27, over donors 0.08, donor x batch residual 0.14 to 0.16; negative
modes of CD16, CD19, CD4 donor-stable and batch-sensitive; positive modes of CD25, CD56, CD69, HLA-DR, FOXP3
donor-dominated; stimulated donors wobble more (within-donor SD 0.31 to 0.40 vs 0.23 to 0.26).

Aim: let HEB correct batches without reference donors (now abstained) from internal anchors, with honest SEs,
without hurting batches it already corrects. This relaxes "reference-only" to "reference donors plus
donor-invariant internal anchors"; reports must say so.

## Approved rows

| ID | Decision | Source |
|---|---|---|
| A0 | Option (b): anchors enter HEB's mixed model as extra observations (not pre-normalisation, not affine-only) | user 2026-11-10 |
| A1 | Anchors: CD3+ median in the CD3 channel, plus the negative mode of each marker passing A2; other positive modes excluded | user 2026-11-10 |
| A2 | A marker's split counts when its two k-means modes are separated by >= 3 pooled SDs | user 2026-11-10 |
| A3 | An anchor informs its own channel only (no latent gain, no cross-channel loadings) | user 2026-11-10 |
| A4 | Synthetic S18 (shared per batch x channel affine effect on the linear scale) and S18w (S18 plus one stimulated-like donor whose anchor wobbles, SD 0.3 asinh per batch); dev seeds 441-460 only | user 2026-11-10 |
| A5 | S18: 3 batches per seed without reference donors, seeded draw | user 2026-11-10 |
| A6 | Anchor rows: residual = known sampling variance var(cells)/n + their own estimated extra variance (separate from donor rows); no donor term | user 2026-11-10 |
| A7 | Modes: k-means k = 2 per sample and marker, seed 1; mode value = cluster median | user 2026-11-10 |
| A8 | T6 check: remove the reference role of all 4 controls in 3 batches (seed 1); compare anchor-corrected vs abstained EMD on surface markers; guardrail: stimulated-minus-unstimulated contrast moves <= 0.05 asinh | user 2026-11-10 |
| A9 | An anchor joins the designs of populations whose anchor-batch median for that marker falls in the anchor's k-means cluster | user 2026-11-10 |
| A10 | S18 effect sizes: offset SD 0.10 and log-gain SD 0.20 per batch x channel (linear scale) | user 2026-11-10 (unaudited judgement) |
| inherited | min_cells 10 (G1), anchor batch01 (Q8), asinh cofactor 6000 (Q4), heb_pool_pop_kv_bap settings (F6, F8) | earlier plans |

## Pass / fail

- Synthetic (S18, S18w): pass = held-out reference RMSE no worse than heb_pool_pop_kv_bap on batches with reference
  donors, coverage >= 0.85 on anchor-only effects, and anchor-only batches better than abstaining. Fail = coverage
  < 0.85, or worse than abstaining under the S18w wobble.
- T6 (A8): pass = EMD below abstaining on surface markers in the stripped batches and the contrast guardrail holds.

## Steps

1. Generator: S18 and S18w in gen_scenarios.py (anchor populations and affine effect), dev seeds 441-460.
2. Engine: anchor summaries (A7) and anchor rows in the HEB design (A6, A9) behind a required argument.
3. Unit tests: anchors never use held-out cells; a batch with anchors and no reference donors gets finite effects.
4. SLURM: run lead vs lead+anchors on S18, S18w and the existing dev sets (regression check), score, figure.
5. T6 check (A8) on SLURM, then figure and ledger record.
