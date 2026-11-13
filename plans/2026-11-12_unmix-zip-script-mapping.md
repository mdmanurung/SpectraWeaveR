# Mapping: zip script `AutoSpectral_unmix_v1.8.8.R` to SpectraWeaveR wrappers (2026-11-12)

Source: `experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmix_v1.8.8.R` (read in full). Wrapper signatures from `R/unmix.R`. AutoSpectral internals were not read. Nothing here was run.

## Step map

| # | Zip script call | Wrapper | How the zip's values pass through |
|---|---|---|---|
| 0 | `get.autospectral.param(cytometer = "aurora", figures = TRUE)` | `sw_unmix_setup(control_dir, cytometer = "aurora", control_file, output_dir, figures = TRUE)` | `cytometer` and `figures` match. |
| 1 | `read.flowSet` + `rectangleGate` + `Subset` + `write.FCS` into `Reference Group Clean/` | **no wrapper** | Do it with flowCore before step 2. Cells gate: FSC-A 6e5 to 4e6, SSC-A 1e5 to 1.6e6. Beads gate: FSC-A 3e5 to 8e5, SSC-A 1.3e6 to 1.9e6. These are the zip's values, kept for reference only. |
| 2 | `create.control.file()` | `sw_unmix_setup(..., control_file = NULL)` | Writes the control file. |
| 3 | Hand edits: peak channel (`V7-A`, `UV9-A`, `R6-A`), marker join, `universal.negative` | **no wrapper** | Edit the CSV, then pass it as `control_file`. Needed columns: filename, fluorophore, marker, channel, control.type, universal.negative. |
| 4 | `check.control.file()` | `sw_unmix_setup(..., control_file = <edited csv>)` | The setup docs name `check.control.file()` for a supplied file. Verify at run time. |
| 5 | `get.spectra.automated(n.candidates = 10000, n.spectral = 500, cosine.threshold = 0.95)` | `sw_unmix_prepare(setup, workflow = "automated", n.candidates = 10000, n.spectral = 500, cosine.threshold = 0.95)` | The three arguments go through `...` under AutoSpectral's names. |
| 6 | `reload.flow.control()` | inside `sw_unmix_prepare()` (automated workflow) | Returned as `$flow_control`. |
| 7 | `get.af.spectra(unstained.sample, spectra, refine = FALSE)` | `sw_unmix_extract_af(unstained_fcs, setup, spectra, refine = FALSE)` | The script stops unless `Reference Group/unstained_AF/` holds exactly one FCS. |
| 8 | `get.spectral.variants(control.dir, control.def.file, asp, spectra)` | `sw_unmix_extract_variants(setup, spectra)` | `n_cells` and `som_dim` stay NULL (AutoSpectral defaults), as in the zip. |
| 9 | `unmix.folder(method = "AutoSpectral", pipeline = "joint", parallel = TRUE, ...)` per subfolder of `Raw/` | `sw_unmix_run(input = <folder>, spectra, setup, flow_control, method = "AutoSpectral", af_spectra, spectra_variants, speed = "slow", pipeline = "joint", output_dir, parallel = TRUE)` | `pipeline` goes through `...`; check it against the 1.8.8 formals. The zip sets no `speed`; the project rule is `speed = "slow"`, with no effect under the joint pipeline. |
| all | steps 0 to 9 in one call | `sw_unmix_pipeline(control_dir, sample_input, unstained_fcs, ...)` | Cannot do steps 1 and 3 (gating, hand edits). Use the stepwise calls for this dataset. |

## Files the zip run consumed

- `Raw/`: 5 sample FCS files.
- `Reference Group/`: single-stain controls (cells and beads), `Unstained (Cells)`, `Negative_1 (Beads)`, and `unstained_AF/` with one FCS.
- `Reference Group Clean/`: made from `Reference Group/` in step 1, plus the edited control CSV.
- The fluorophore-to-marker table the script reads from a `W:\` path. It only labels markers, so the unmixing does not depend on it.

Outputs, not inputs: `figure_*`, `table_spectra/`, `AutoSpectral_unmixed/`. Editor leftovers: `.RDataTmp`, `.Rproj.user`, `Project.Rproj`.

## Gaps to settle before an example script exists

1. Debris gating and control-file editing have no wrapper. The example would call flowCore and base R for them.
2. Whether `sw_unmix_run()` passes `pipeline = "joint"` through `...` is untested here. The parity test in `tests/testthat/` checks argument names against the 1.8.8 formals.
3. The marker table is outside the zip. The example needs its own fluorophore-to-marker CSV.
