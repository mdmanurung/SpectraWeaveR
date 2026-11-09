# 01: Mendeley 7-colour mouse tissue (AutoSpectral article 01 dataset)

Status: **v2 run finished** (`runs/2026-11-09_mendeley-7color-v2`). Results are in `../RESULTS.md`.

## Question and claim

The user stated the question and claim on 2026-11-09. Their words are in `plans/beta_unmix_report_figure_plan.md`. The draft below is the earlier wording:
- Question: does a SpectraWeaveR run of the AutoSpectral article-01 workflow use sensible reference spectra and give unmixed data consistent with unmodified AutoSpectral and with SpectroFlo?
- Claim: the spectra and unmixed output agree with the comparators within the differences listed in the report tables.

## Data

`autospectral_mendeley` in `../../datasets/registry.csv` (Mendeley Data `ch5dnspd79`, 2.1 GB, SHA-256 checked).
Raw: `Raw/Set1/{Reference Group, Unstained, Stained}`; comparators: `Unmixed/Set1/...` (SpectroFlo). Control file: `fcs_control_file.csv`, shipped with the dataset.

## Approved choices (user, 2026-11-08)

| Choice | Decision | Source |
|---|---|---|
| Meaning of "working" | Run completes; report builds; our spectra are compared with spectracle's and our unmixed files with SpectroFlo's. No pass/fail cutoffs. | user, Q1 |
| Where it runs | SLURM, partition `all`, account `div2-para` | user, Q2 and SLURM question |
| Figures | F1 spectra overlay, F2 similarity matrix, F3 unmixed distributions, F4 biplot vs SpectroFlo; F5 stays a table | user, `plans/beta_unmix_report_figure_plan.md` |
| Verdicts | Numbers only; pass/warn/fail only where a limit is supplied | user |

| F3/F4 event sample | 50,000 events per file, seed 1 | user, 2026-11-09 |
| Symmetric-log sigma (F3, F4) | 1000 | user, 2026-11-09 |
| Biplot bins (F4) | 200 per axis | user, 2026-11-09 |
| spectracle AF row (F1) | kept as a separate table; F1 compares fluorophores only | user, 2026-11-09 |
| Methods on the stained lung file | OLS, AutoSpectral fast, AutoSpectral slow (three runs, same file) | user, 2026-11-09 |

## Still to be set by the user

Nothing open. The verdict limits are SIF warn 4, fail 6 (`sw_unmix_limits_autospectral()`, user choice 2026-11-09). Other metrics are not judged.

## Steps (workflow of article 01, wrapped)

1. `sw_unmix_setup()` with the dataset control file.
2. `sw_unmix_prepare()` (reference spectra).
3. `sw_unmix_extract_af()` on the matching unstained tissue file; `sw_unmix_extract_variants()`.
4. `sw_unmix_run()` on the stained files: an OLS baseline and the AutoSpectral method.
5. `sw_unmix_diagnose()`, `sw_unmix_report()`.
6. Compare against `Unmixed/Set1/Stained/*.fcs`.

## Outputs

`../../runs/<date>_mendeley-7color/` with `report/` (HTML, CSV tables, `provenance.csv`), figures and `sessionInfo.txt`.
