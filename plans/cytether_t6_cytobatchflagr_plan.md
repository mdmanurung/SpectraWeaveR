# T6: cyTether on CytoBatchFlagR (real spectral data, repeated controls)

Question: on real Aurora data with the same control donors in every batch, are `paired_shift` and HEB not inferior to cyCombine, CytoNorm and CytofIn on a held-out control donor?

## Data

Zenodo 15388817 (CC-BY-4.0), at `$SW_DATA_DIR/cytobatchflagr_zenodo/`, per-file md5 verified (`md5_verified.tsv`).

- 72 FCS files: 18 batches (batch01 to batch18) x 4 control aliquots. `tcu1` and `tcu2` are unstimulated; `tcs3` and `tcs4` were stimulated 7 days with anti-CD3/CD28 and IL-2.
- The files were pre-gated in FlowJo on CD45+CD14- live single cells and are already unmixed (Aurora, `$SPILLOVER` present). There are 28 marker channels plus AF, with 25,055 to 970,597 events per file.
- Metadata: `Table_2.csv` (file, batch, control, sample_id); panel in `Table_1.pdf`.
- No study level, no patient samples.

## Approved choices

| # | Choice | Source |
|---|---|---|
| N1 | Stand-in primary: takes Z247's slot in the Q11 family (2 datasets x 3 tools) until Z247 arrives; Z247 is then added as an extra | user 2026-11-09 |
| N2b | 20 markers with the same antibody and fluorophore in all 18 batches. Replaces N2 (26 markers) after the panel was found to change: APC GITR to CCR4 at batch05; HELIOS and CCR7 swap Alexa 488 and 700 at batch08; CD27, CD28 and CD163 are absent in later batches. Also dropped: CD45, CD14 (gating markers), AF and scatter | user 2026-11-09 |
| N3 | 50,000 cells per file, subsampling seed 1. The two smaller files (25,055 and 36,969) contribute all their cells | user 2026-11-09 |
| N4 | HEB hierarchy flat: study = batch (the `heb_flat` structure), `pool = "population"` | user 2026-11-09 |
| N5 | FlowSOM 10x10 grid, 20 metaclusters, fit on the anchor batch's reference cells, then frozen | user 2026-11-09 (grid) and Q5 |
| N6 | Seeds: tool defaults (cyCombine 473; FlowSOM and CytoNorm 1; subsampling 1), recorded in the run manifest | user 2026-11-09 |
| N7 | CytofIn mode `meanshift` (the README's worked example; no documented default). Pooled-anchor target, so our anchor has no effect on it | user 2026-11-09 |
| O1 | OpenProblems scores use separate 2-donor folds: one unstimulated and one stimulated control held out (4 pairings), the other two as references. Groups = stimulated vs unstimulated, 18 samples each. The EMD primary keeps the 1-donor folds | user 2026-11-09 |
| O2 | Scaling over all three task controls (no, perfect, shuffle integration): min-max, NaN to 0, clip to [0, 1], as on the website | user 2026-11-09 |
| O3 | Replicate split: seeded random 50/50 of each evaluated sample's cells (seed 1); each method runs once per split | user 2026-11-09 |
| O4 | Marker types: lineage = CD3 CD4 CD8 CD19 CD56 CD16 CD45RA CD103 CCR6; functional = GZMB GNLY PFN CD69 HLADR CD25 PD1 CD39 FOXP3 TIGIT CD95 | user 2026-11-09 |
| O5 | perfect_integration batch = batch01; metric FlowSOM 10x10, 20 metaclusters (as N5) | user 2026-11-09 |
| O7 | perfect_integration for T6 (the task's version needs one batch per donor per split): each evaluated (donor, batch) sample in split k is replaced by that donor's batch01 cells of split k, labelled with that batch | user 2026-11-10 |
| S1-S4 | Scale diagnostic (after fair cyCombine beat heb_kv_bap in 3 of 4 folds, emd_scores_v5): split held-out EMD into total and shape (each file centred on its median); prototype heb_kv_bap_scale = HEB shift then x' = m_b + s (x - m_b), s = MAD(anchor-batch calibrating cells) / MAD(batch calibrating cells) per batch x population x marker, donors pooled, >= 10 cells; methods none, heb_kv_bap, cycombine_ref, cytonorm, 4 folds, EMD max over batch pairs; stop when shape explains < 50% of the HEB - cycombine_ref gap. Scripts 08_emd_decompose.R, 09_fig_emd_decompose.py | user 2026-11-10 (Q1-Q4) |
| O6 | OpenProblems thresholds inherited from task_cyto_batch_integration v0.0.1-rc3 (`6cdd9e8`), as fixed in `experiments/benchmarks/_common/openproblems_cyto_metrics.py`; LISI uses exact kNN | user 2026-11-09 (confirm) |
| Q1, D1 | Endpoint: CytoNorm 2020 EMD (absolute, 0.1 bins, maximum over batch pairs) per population and marker on the held-out donor; median over population x marker | inherited: `plans/cytether_benchmark_plan.md`, `plans/cytether_benchmark_strategy.md` |
| Q4 | arcsinh cofactor 6000 (spectral) | inherited: `plans/cytether_benchmark_plan.md` |
| Q7, Q12, D7 | Non-inferior if within 10% of the best competitor's median, on EMD and on each OpenProblems score, direction-aware | inherited |
| Q8 | Anchor: most reference donors, ties to earliest, so batch01. Leave-one-donor-out: 4 folds, each holding out one control, with the other 3 as references | inherited |
| Q10 | Competitors in documented reference mode, pinned commits (cyCombine `fcd42db`, CytoNorm `278f352`, cytofin `cb3e24b`) | inherited |
| Q11 | Holm over the family; resample batches and donors, never cells | inherited |
| G1 | `min_cells` 10 | inherited: `experiments/benchmarks/05_cytether-synthetic/run_method.R` |

## Methods

`none`, `paired_shift`, `heb` (flat, pool population), cyCombine, CytoNorm and CytofIn. The CytofIn adapter does not exist yet and must be built under Q10.

## Steps

1. Read and subsample: `flowCore::read.FCS(truncate_max_range = FALSE)`, 50,000 cells per file with seed 1, keep the 26 markers, then arcsinh with cofactor 6000. Write one parquet with sample metadata to a fresh run directory.
2. Per fold: run FlowSOM on the anchor's 3 reference files (10x10 grid, 20 metaclusters, seed 1) and map every cell to the frozen metaclusters.
3. Per fold and method: fit on the 3 reference donors in all batches, then apply to all cells. The held-out donor never enters fitting.
4. Score the held-out donor: EMD per metacluster x marker over batch pairs (primary), the cyCombine MAD score, and the 9 OpenProblems scores (D5). Safety gates: frozen invariance, and the protected contrast (within-batch stimulated minus unstimulated difference, change against raw).
5. Summarise: median EMD per method and fold, bootstrap over batches, the Holm family per Q11.

## Figures (plot-first)

- Held-out EMD per method, one point per fold over the median (primary).
- Per-batch held-out EMD heatmap, method x batch.
- Protected-contrast change per method.
- Marker density per batch for 3 markers, before and after correction, for the held-out donor.

## Open work before running

- OpenProblems' 9 scores have to be implemented from the task repo's metric definitions, because the pipeline needs its S3 data. D7 makes them part of the claim.
- Build the CytofIn adapter.
- Check whether `tcs3` and `tcs4` come from different donors from `tcu1` and `tcu2` (Table_1.pdf). This changes what the protected contrast means.
