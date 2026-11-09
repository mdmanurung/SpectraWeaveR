# 02: OMIP-109 45-colour

Status: **run completed**. Results are in `../RESULTS.md`.

## Settings carried over from benchmark 01 (user confirmed, 2026-11-09)

| Choice | Decision |
|---|---|
| Methods | OLS, AutoSpectral fast, AutoSpectral slow |
| Seed | 1 |
| F3/F4 event sample | 50,000 events |
| Symmetric-log sigma (F3, F4) | 1000 |
| Biplot bins (F4) | 200 per axis |
| Verdict limits | `limits = sw_unmix_limits_autospectral()`: `max_sif` warn 4, fail 6. All other metrics are not judged. |

Run output: `../../runs/2026-11-09_omip109-45color/`. Settings: `config.R`. Entry script: `run.R`.
