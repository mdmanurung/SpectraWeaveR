# Unmix event metrics: one-stop diagnostics suite (approved 2026-11-15)

Goal: port the event-level unmixing metrics of benchmark 05 (101f, 103c, 103d, 103f, 103g) into the package and wire them
into `sw_unmix_diagnose()` and `sw_unmix_report()`, next to the spectra, SIF, condition number, spike-in and screen checks.

## Decisions (user, 2026-11-15)

| Item | Decision |
|---|---|
| Metric groups | negatives + spread; stain index; pairwise line check; held-out G1 + extraction reproducibility |
| API | new exported functions + `event_metrics =` and `heldout =` slots in `sw_unmix_diagnose()`, report sections |
| Positive / negative | AutoSpectral 1.8.8 flat edges (unstained q0.005 floor, q0.995 edge, `variants$neg.thresholds` / `$thresholds`) plus the spread-aware edge of `AutoSpectral::get.spread.thresholds()` (kappa required) |
| Stain index | AutoSpectral `calculate.ssi()` convention (median pos - median neg) / (2 MAD neg) plus the 3D quantile SI (spread_q) |
| `sw_unmix_quality()` | removed now (breaking); `.sw_channel_stats()` keeps mean, sd, cv, pct_negative in diagnose |

## Functions

- `sw_unmix_edges(variants, unstained, probs)`: per-channel floor and edge (linear), from an AutoSpectral variants object or
  from an unmixed unstained matrix at quantiles the caller passes.
- `sw_unmix_event_metrics(data, edges, cofactor, negative_cut, spread_probs, min_pos, hyper_k, spread, kappa, pairs, n_bins,
  min_cells)`: per file x channel `pct_below_cut` (asinh < negative_cut), `pct_below_floor`, `neg_median`, `neg_rsd`, `spread_q`,
  `pileup`, `hyper_share`, `n_pos`, `si_mad`, `si_q`; per file x pair (i+ j-) Theil-Sen `slope` over `n_bins` quantile-bin
  medians and Spearman `rho`. Every knob required; `spread`/`kappa` NULL means flat edges; `pairs` NULL skips the line check.
- `sw_unmix_heldout(...)`: wraps AutoSpectral (Suggests): split controls, extract on train, unmix test + gated unstained,
  score off-target spread_q and bias; optional same-seed repeat extraction. Every knob required.
- Diagnose: check ids `event_channel` (item `file|channel`), `event_pair` (item `file|i|j`, metrics `abs_slope`, `abs_rho`),
  `heldout` (item `control|channel`). New limit names map to these metrics; nothing judged without a limit.
- Report: `event_channels.csv`, `event_pairs.csv`, `heldout_scores.csv`, `heldout_repro.csv`; HTML sections with per-channel
  medians over files, the 25 pairs with the largest |slope|, and held-out medians by control type. Pair rows are left out
  of the full HTML checks table (CSV only).
