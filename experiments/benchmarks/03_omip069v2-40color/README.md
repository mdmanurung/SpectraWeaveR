# 03: OMIP-069v2 40-colour PBMC

Status: **all runs completed**. The `-cells-set-v2` and `-beads-set-v2` runs are the current ones; the first `-cells-set` and `-beads-set` attempts failed on the first runner version and stay as history. Results are in `../RESULTS.md`.

## Settings carried over from benchmark 01 (user confirmed, 2026-11-09)

| Choice | Decision |
|---|---|
| Methods | OLS, AutoSpectral fast, AutoSpectral slow |
| Seed | 1 |
| F3/F4 event sample | 50,000 events |
| Symmetric-log sigma (F3, F4) | 1000 |
| Biplot bins (F4) | 200 per axis |
| Verdict limits | `limits = sw_unmix_limits_autospectral()`: `max_sif` warn 4, fail 6. All other metrics are not judged. |

## Control file and the two spectra sets

The rerun follows the AutoSpectral article `02b_Control_File_Complex_Example`. One control file holds all controls: `control_combined.csv`. It has cell rows and bead rows named `<dye> beads`. The unstained cells row is the AF row. The bead negative row is `Negative`.

Two runs unmix with that file, one per spectra set:

| Spectra set | Config | Entry script | Run folder |
|---|---|---|---|
| cells | `config_cells_set.R` | `run_cells_set.R` | `../../runs/2026-11-09_omip069v2-cells-set-v2` |
| beads | `config_beads_set.R` | `run_beads_set.R` | `../../runs/2026-11-09_omip069v2-beads-set-v2` |

The cells set uses cell spectra only. The beads set uses the bead spectrum where one exists, else the cell spectrum.

## History

The older separate-file runs stay as history: `-cells` and `-beads-v2` (configs `config_cells.R`, `config_beads.R`).
