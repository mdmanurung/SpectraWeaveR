# Plan: gating scenario vignettes, systematic QC and cleaning, and a summary report per step

Date: 2026-11-13. Status: **decisions answered 2026-11-13, awaiting final confirm**. No code written.
Sources: `experiments/repo_crawls/spectraweaver_incorporation_proposals.md` (ids M01 to M75),
`experiments/repo_crawls/cytometry_github_studies_and_pipelines_catalog.md`, `plans/2026-11-10_gating.md`,
`vignettes/cleaning-gate.qmd`, `vignettes/tutorial-02-gating-signal-qc.qmd`, `R/` listing.

## Where the package stands

- Gating is largely built (plan 2026-11-10): `sw_gate_density_cut/apply`, `sw_gate_transfer`, `sw_gate_qc`,
  `sw_gate_plot`, `sw_gate_read_flowjo`, `sw_gate_read_gatingml`, `sw_gate_draw`, flowMagic wrappers.
  The FlowJo reader exists, so the manual-vs-auto agreement item (M62) has its input.
- Vignettes show one gating route on synthetic data (`cleaning-gate.qmd`, `tutorial-02`). No vignette
  covers controls, FMO, longitudinal reuse, failed gates or agreement with manual gates.
- Only the unmix step has a report function (`sw_unmix_report()`). Import, clean, QC, gate, batch,
  cluster, annotate and differential return objects and plots, with no run folder of tables and figures.
- `sw_qc_batch()` keeps a 5-column summary and drops the PeacoQC object (R/qc.R:170-195).

## Part A. Report contract (all steps)

One shape for every step, copied from `sw_unmix_report()`:

- `sw_report_<step>(result, out_dir, ...)` writes to a **new** run folder (stops if it exists):
  `tables/*.csv`, `figures/*.png` (+ `source_data/<panel>.csv`), `provenance.json`, `sessionInfo.txt`.
- Returns `list(tables, figures, source_data, dir)`. Numbers only: no pass or fail, no per-sample
  roll-up score (user decision Q14). Thin rows are `not_assessed` with a reason.
- Every removal step reports events in, events out, events removed. A silent filter is a bug.
- Figures show every sample as a point or line over any summary. Axes carry units. Panels share axes.
- A notebook (not the function) writes the prose. Quarto HTML stays a local file.

## Part B. Steps, wrappers to add, vignette scenarios, report panels

### B1. Import and harmonise (new report: `sw_report_import`)

- Add: `sw_fcs_integrity_audit` (M22), `sw_fcs_compensation_state` (M30), `sw_fcs_marker_dye_map` (M38),
  `sw_marker_names` (M39), `sw_fcs_panel_match` (M21), `sw_fcs_harmonize` (M17).
- Report: file x check table (keyword and channel differences against a reference file); marker x batch
  dye map; panel-match table. Figure: file x marker presence tile map.

### B2. Signal QC and cleaning, systematic (new report: `sw_report_clean`)

Order follows the validated pipeline: pre-flight, per-file signal audit, margins, transform, PeacoQC,
cleaning gates, sample ledger.

| Stage | Wrapper (source id) | Reports |
|---|---|---|
| Pre-flight | B1 functions; `sw_fcs_signal_audit()` | saturation, dead channel, non-finite, rate, correlation per file |
| Drift and instability | `sw_signal_segment_drift` (M60), `sw_fcs_rate_sections` (M49), `sw_signal_changepoints` (M63, wraps `changepoint` from Suggests) | per file and channel: segment range, largest jump, changepoints |
| Cofactor evidence | `sw_cofactor_scan` (M20) | bimodality and negative spread over a required cofactor grid; no pick |
| PeacoQC | `sw_qc_breakdown` (M59, needs `sw_qc_batch()` to keep the result object) | IT, MAD and consecutive-bin split; channel contribution; events per bin |
| Sensitivity | `sw_qc_sweep` (M19), `sw_qc_removal_audit` (M51), `sw_qc_mask_validate` (M56, only with a truth mask) | removed share across a required IT_limit and MAD grid; kept vs removed divergence over time |
| Cleaning gates | `sw_gate_density_cut(method = "tail_slope")` (M61), `sw_gate_threshold_sweep` (M18) | debris, singlet, live cut candidates and the share kept at each |
| Ledger | `sw_sample_ledger` (M11), `sw_count_floor` (M32), `sw_retention_association` (M24) | one row per sample and removal step with a written reason; count floor; retention vs batch and condition |

Report tables: `event_waterfall.csv` (sample x step), `qc_breakdown.csv`, `drift.csv`, `sweep.csv`,
`ledger.csv`. Report figures:

1. Event retention per sample across steps (lines per sample, steps on x).
2. PeacoQC removal over acquisition time per sample (ribbon of kept vs removed, one facet per sample group).
3. Drift tile map, sample x channel, colour = range-normalised drift, rows ordered by largest drift.
4. Removed share across the sweep grid, one line per sample.
5. Cleaning-gate overlays on raw FSC-A vs SSC-A and viability, via `sw_gate_plot()`.

### B3. Gating (`sw_report_gate`; vignette `gating-scenarios.qmd`)

Wrappers: `sw_gate_control_cut` (M04), `sw_fmo_cutoff` (M03), `sw_pop_control_subtract` (M23),
`sw_control_fpr_limit` (M08), `sw_gate_membership` (M05), `sw_gate_agreement` (M62),
`sw_gate_impute` (M42), `sw_gate_cut_consensus` (M31), `sw_bool_combinations` (M02),
`sw_gate_qc` new checks (M47).

Vignette scenarios (each uses one repo as its model; values in examples are illustrative and labelled so):

| # | Scenario | Model repo | Functions shown |
|---|---|---|---|
| S1 | One sample: template CSV for margins, debris, singlets, live, lymphocytes; `groupBy` and `collapseDataForGating` | FredHutch immune-flow-gating, HIPCCyto, flowFun | `sw_gate_template`, `sw_gate_run` |
| S2 | Stimulated vs unstimulated control per donor (ICS) | TB-TASA, paper-opencyto | `sw_gate_control_cut`, `sw_pop_control_subtract`, `sw_control_fpr_limit` |
| S3 | FMO-guided positivity | cytotypr, Fox et al. | `sw_fmo_cutoff`, `sw_bool_combinations` |
| S4 | Reuse one template over timepoints or batches; move gates with landmarks | Chuong, McCarthy, DuraClone | `sw_gate_transfer`, `sw_gate_cut_consensus`, `sw_gate_qc` |
| S5 | Spectral specifics: gate unmixed data on arcsinh 6000, scatter-dependent markers | den Braanker, SoegaardLab | `sw_gate_density_cut`, M09 later |
| S6 | Agreement with manual FlowJo gates | everything-flow-cytometry, flowWorkspace | `sw_gate_read_flowjo`, `sw_gate_agreement` |
| S7 | Failure review: a gate that fails, a sample that needs imputing | HIPCCyto, Bin He Lab | `sw_gate_qc`, `sw_gate_impute`, `sw_gate_plot` |

Report tables: gate x sample stats (events, share of parent, edge share, centroid shift, kept-share z),
agreement table, control-limit table. Figures: gate overlays per sample (reference vs adapted), gate-shift
dots per sample and gate, parent-share dots per sample.

### B4. Unmixing (extend existing `sw_unmix_report`)

Add `sw_unmix_secondary_index` (M65), `sw_unmix_control_residual` (M46), `sw_unmix_residual` (M73, after
checking AutoSpectral 1.8.8 exports), `sw_control_peak_check` (M35). Only the report check rows and
dashboard tab change. Blocked on matched evaluation tubes (ROADMAP B3 to B10); build the functions, run on
the existing BMV data only.

### B5. Batch (`sw_report_batch`) — frozen in ROADMAP

Needs an approved plan row first: `sw_correct_node_audit` (M67), `sw_correct_metric_anchors` (M52),
`sw_correct_peak_consistency` (M53), `sw_replicate_spread` (M10), `sw_plate_position_scores` (M12).
Report: stack existing `sw_batch_scores`, `sw_batch_otd`, `sw_batch_mixing`, `sw_correct_metrics`,
`sw_batch_decide` outputs through `sw_batch_matrix()`.

### B6. Clustering (`sw_report_cluster`)

Wrappers: `sw_cluster_map_outliers` (M29), `sw_marker_resolution` (M72), `sw_cluster_label_agreement`
(M70), `sw_cluster_overlap_scores` (M28), `sw_cluster_replicate_agreement` (M54, frozen side).
Report: cluster size and share per sample, marker medians per cluster, map-distance outliers, stability
(`sw_cluster_stability`). Figures: cluster share dots per sample, marker x cluster tile map, distance density.

### B7. Annotation and differential (`sw_report_annotate`, `sw_report_diff`)

Wrappers: `sw_diff_abundance_glmm` (M37), `sw_variance_share` (M33), `sw_marker_nrs` (M69),
`sw_stim_response` (M01). Report: annotation table with score per cluster, abundance model table with the
multiple-testing family named, volcano and per-donor box plots over points. Smallest evidence base in the
survey; the critic asked for a separate iidx / cyCONDOR / CATALYST pass first.

## Build order (proposed)

1. Part A contract, `sw_report_clean` with M59 (small `sw_qc_batch()` change) and M60. One vignette
   `qc-and-cleaning.qmd`.
2. `sw_report_gate`, gating wrappers M03, M04, M23, M08, M62, then `gating-scenarios.qmd` S1 to S7.
3. `sw_report_import` and B1 wrappers.
4. `sw_report_cluster`, then annotation and differential after the extra survey.
5. Unmix additions, then batch (after you lift the freeze).

## Decisions for the user (each needs approval)

| # | Choice | Recommended | Why |
|---|---|---|---|
| RD1 | Report contract in Part A | Approve as written | One shape across steps; matches `sw_unmix_report()` |
| RD2 | Vignette data | Synthetic frames for every scenario, plus one real public spectral set (FR-FCM-Z4KT, den Braanker) in a non-run (`eval: false`) section | Package ships no data; synthetic runs in CI |
| RD3 | Build order | As above | QC and cleaning is your first ask; gating mostly built |
| RD4 | Comparator cleaners (flowAI, flowCut) | Wrap flowCut/flowAI only for the `sw_qc_sweep` comparison, from Suggests | They are the comparators in three repos; GPL, so wrap only |
| RD5 | Batch items (B5) | Keep frozen | ROADMAP A2 freeze |
| RD6 | Section B3 example values (cuts, quantiles, grids) | Mark every value "illustrative, choose your own"; each function keeps it a required argument | impactful-choices rule |
| RD7 | Extra survey for annotate, differential, report/app, and the repos not yet read (B02, B03, B13, C03, C07, C09, D05, D14) | Run it before step 4 | Critic found these uncovered |

Every impactful value in the proposals file (section 4 there, over 100 rows) stays a required argument.
Before each build, the row list for that function is collapsed against `ROADMAP.md` choices and marked
settled, inherited-confirm or new.

## Decisions answered (user, 2026-11-13)

| Row | Answer |
|---|---|
| RD1 / Q1 report contract | Run folder, numbers only, no roll-up score. |
| RD2 / Q2 vignette data | Synthetic frames everywhere, plus one real public spectral set (FR-FCM-Z4KT, availability unverified) in a non-run section. |
| RD4 / Q3 cleaners | Wrap flowAI and flowCut from Suggests, for the QC sweep comparison only. |
| RD7 / Q4 extra survey | Run **now, before any build** (changed from "before step 4"). |
| RD3 / Q5 build order | **Import first** (changed from "QC and cleaning first"): B1 import, then B2 QC and cleaning, then B3 gating, B6 cluster, B7 annotate and differential. Pipeline order. |
| RD5 / Q6 batch | Stay frozen (ROADMAP A2). |
| RD6 / Q7 example values | Labelled illustrative; every impactful value stays a required argument. |

Order now: extra survey, then B1, B2, B3, B6, B7, then B4 additions. The "Build order" list above is superseded.
