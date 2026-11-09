# Beta unmixing benchmarks: results

Runs of 2026-11-09 (folders in `experiments/runs/`, git-ignored). Every run used the shared runner `_common/run_unmix.R`. Figures come from `_common/render_figures.R`. Finished run folders are never edited, so re-rendered figures sit in `figures_v2/`.

Settings carried over from the figure plan (`plans/beta_unmix_report_figure_plan.md`): 50,000 events, seed 1, pseudo-log sigma 1000, 200 bins, methods OLS, AutoSpectral fast and AutoSpectral slow. The user approved F1 to F4 on 2026-11-08. They have not yet confirmed these settings for the OMIP datasets, nor the rename maps in the configs.

## Status per dataset

| Run folder | Result |
|---|---|
| `2026-11-09_mendeley-7color-v2` | Completed. 7 fluorophores, 64 detectors. |
| `2026-11-09_omip109-45color` | Completed. 37 fluorophores matched to spectracle, 6 unmatched (`spectra_unmatched.csv`). |
| `2026-11-09_omip069v2-cells` | Completed. |
| `2026-11-09_omip069v2-cells-set-v2`, `-beads-set-v2` | Completed. One combined control file (`03_omip069v2-40color/control_combined.csv`), one run per spectra set, 40 fluorophores each. `-cells-set` and `-beads-set` failed on the first runner version (variants names); the v2 runner derives each set's control file from the combined file. |
| `2026-11-09_omip069v2-beads-v2` | Completed (earlier separate run). 39 fluorophores. `-beads` failed (no unstained AF row); `-40color` failed (draft control file). |
| `2026-11-09_mendeley-correction` | Completed. Two of seven spectra changed and were accepted. |

## Findings

**Mendeley spectra match spectracle.** Six of seven angles are 0.1 to 2.2 degrees. BUV805 is the largest at 4.5 degrees, with small differences in the flat region between detectors. Figure inspected: `experiments/runs/2026-11-09_mendeley-7color-v2/figures_v2/panels/spectra-overlay.png`. Source: `spectra_compare.csv`.

**OLS differs from the AutoSpectral methods in the negative population.** In BUV395, BUV805, BV421 and PE the OLS negative peak sits right of zero. Both AutoSpectral methods and SpectroFlo centre near zero. The two AutoSpectral methods overlap closely. SpectroFlo has a wider negative population in BUV395 and BUV805 than AutoSpectral. Figure inspected: `experiments/runs/2026-11-09_mendeley-7color-v2/figures_v2/panels/unmixed-distributions.png`. This is a visual read on 50,000 events from one file, so it supports a description, not a ranking.

**OMIP-069v2 beads have close pairs inside dye families.** The highest off-diagonal cosine is 0.975 (Super Bright 436 vs BV421). Other high pairs sit in the BV421, BV480/BV510, BB515, PE, PerCP and APC families. Figure inspected: `experiments/runs/2026-11-09_omip069v2-beads-v2/figures_v2/panels/similarity-matrix.png`.

**OMIP-109 populations separate the same way across methods.** In the BUV395 vs BV421 biplot the four methods (OLS, AutoSpectral fast, AutoSpectral slow, SpectroFlo) show the same occupied regions and no extra diagonal smear in one method only. Figure inspected: `experiments/runs/2026-11-09_omip109-45color/figures_v2/panels/biplot.png`.

**Cells vs beads spectra differ most in PerCP-Cy5.5.** Of 39 fluorophores, the median angle is 1.37 degrees and the maximum 19.29 degrees (PerCP-Cy5.5).
`unsupported: no inspected figure` for the cells-vs-beads comparison. `table fallback: one scalar per fluorophore from a pairwise comparison, source experiments/runs/2026-11-09_omip069v2-beads-v2/cells_vs_beads_spectra.csv`.

**Correction run.** `sw_unmix_correct()` changed BUV805 (7.05 degrees) and PE (2.87 degrees) and left five spectra unchanged. In PE-Cy7 and APC, whose spectra did not change, the corrected and uncorrected curves overlay. The BUV805 and PE panels show small differences in the negative peak height. The figure cannot say whether the changes are improvements. Figure inspected: `experiments/runs/2026-11-09_mendeley-correction/figures/unmixed-distributions.png`. Source: `correction_changes.csv`. Not inspected: the biplot and spectra overlay of that run.

**SIF flags the close dye pairs (OMIP-069v2, SIF limits 4 and 6).** In the cells set, 4 of 40 fluorophores fail and 6 warn. In the beads set, 5 fail and 5 warn. The largest values are BV510 (13.6 cells, 13.3 beads) and Spark Violet 538 (11.2, 11.3), partners of each other. Super Bright 436 (7.4, 8.0) and Alexa Fluor 647 (7.4, 7.1) also fail. The similarity matrix shows the same pairs as bright off-diagonal cells: BV510/Spark Violet 538, Super Bright 436/BV421 (0.975) and APC/Alexa Fluor 647. Figure inspected: `experiments/runs/2026-11-09_omip069v2-beads-set-v2/figures/panels/similarity-matrix.png`. Source: `report/checks.csv`, metric `max_sif`. The cutoffs 4 and 6 come from AutoSpectral article 01 and have no further source; Mage et al. use 3 and 8 for the same quantity, so a different choice changes how many dyes fail.

**A flagged channel shows a wider negative population in OLS and SpectroFlo.** In Super Bright 436 (SIF 7.4) the AutoSpectral fast and slow curves have a narrow peak at 0. The OLS and SpectroFlo curves are wide and flat across about -10^4 to 10^4. BUV395, with a passing SIF, shows all four methods overlaid. Figure inspected: `experiments/runs/2026-11-09_omip069v2-cells-set-v2/figures/panels/unmixed-distributions.png`. This is one file and four channels, so it shows an association, not a rule.

**Cells vs beads, set runs.** The comparison reproduces the earlier numbers: 40 fluorophores, median 1.37 degrees, maximum 19.29 degrees (PerCP-Cy5.5). `table fallback: one number per fluorophore from a pairwise comparison, source experiments/runs/2026-11-09_omip069v2-beads-set-v2/cells_vs_beads_spectra.csv`. Each set unmixes the same stained file, so the two sets can be compared channel by channel.

**Method note.** The set runs derive a control file per set from `control_combined.csv` and extract each row on its own, instead of subsetting the spectra of one run. Each row is extracted on its own, so the spectra match what a subset of one run gives. I did not test that equality.

**Metric note.** A first version of the SIF check took the largest entry of each row of the hotspot matrix. Mage et al. define SIF as the diagonal only, so the metric now uses the diagonal and names the largest off-diagonal partner in the reason column.

## Dashboards

Interactive HTML dashboards for these runs are in `experiments/runs/2026-11-09_dashboard/` (`omip069v2_cells_vs_beads.html`, `mendeley_7color.html`, `omip109_45color.html`). They are built by `sw_unmix_dashboard()` with `pseudo_log_sigma = 1000` and `bins = 200`, the approved figure settings. Tabs inspected in a headless browser: SIF, Spectra, Similarity, Distributions, Biplot and Compare for OMIP-069v2; SIF for Mendeley. The OMIP-109 file was built but its tabs were not inspected.

Rebuilt dashboards from the same run folders are in `experiments/runs/2026-11-09_dashboard_v2/` (`*_before.html` with the first front end, `*_after.html` with the current one, plus a two-run and an eight-run file). The embedded payloads (`*.json`) decode to the same numbers; the biplot bins travel as delta-coded keys. The runner now writes `report/dashboard.html` for new runs.

The Mendeley and OMIP-109 runs predate the `max_sif` check. Their SIF values are computed from `report/spectra.csv` and labelled "computed", with no verdict. For Mendeley the values are 1.00 to 1.08. This is a display of a saved-spectra calculation, not a saved check: `table fallback: one number per dye, source report/spectra.csv`.

## Open items

- No condition-number, cosine or negative-spread limit has a sourced cutoff, so those metrics stay `not_judged`.
- The user set the SIF cutoffs to 4 and 6 (AutoSpectral). Mage et al. use 3 and 8.
- Which spectra set (cells or beads) a study should unmix with is not decided. The user asked to run both.
- The F3 channel subset and the rename maps in the OMIP configs are still the first-run choices.

## Test status (2026-11-09, pixi env)

Full `devtools::test()`: 1768 expectations, 11 failures, in 11 files, all of them older files (annotate, assistant, batch_correct, cluster, composable_pipeline, differential, gate, gating_utils, pipeline, qc, targets_backend). The count matches the baseline of 11 failing files at the original commit. I did not diff the failure lists by name. All new test files pass. `devtools::document()` changed NAMESPACE by one blank line, so the hand-merge was faithful.
