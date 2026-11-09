# 04: unmixing-error correction on the Mendeley 7-colour lung file

Status: **finished** (`runs/2026-11-09_mendeley-correction`). Two of seven spectra changed and were accepted. Results are in `../RESULTS.md`.

## Question and claim

The user stated the question and claim for the unmix report on 2026-11-09 (`plans/beta_unmix_report_figure_plan.md`). This benchmark's own draft question: does AutoSpectral's two-step correction (`correct.unmixing.signatures`, then `fix.my.unmix`), applied through `sw_unmix_correct()`, change the reference spectra, and does unmixing with the corrected spectra differ from unmixing with the original ones?

No improvement is claimed or tested for. AutoSpectral article 19 calls its own error-correction idea speculative, and the correction functions are an extension judged against the unmodified baseline (ROADMAP "Gold standard").

## Design

- Same data and controls as benchmark 01 (settings read from its `config.R`).
- Correction steps: both, in AutoSpectral's required order. All tuning arguments keep AutoSpectral's defaults.
- Unmixing: OLS and AutoSpectral fast, each with original and corrected spectra. Variants are re-extracted from the corrected spectra, as `fix.my.unmix` documents.
- Comparator: the SpectroFlo unmixed file from the dataset.
- Figures F1 (original vs corrected spectra), F3 and F4 (the four runs plus SpectroFlo). Sample size, seed, sigma and bins are the benchmark-01 values.

## Outputs

`runs/<date>_mendeley-correction/`: `correction_changes.csv` (per fluorophore angle and change), `report/`, `figures/`, `unmixed/<run>/`, `provenance`, `sessionInfo.txt`.
