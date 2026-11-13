# SpectraWeaveR proposals from external repos

## 1. Summary

1. Build six items first: M60, M62, M03, M29, M59 and M72. Each scores value 4 with fit 4 or 5.
2. M60 `sw_signal_segment_drift` adds range-normalised drift to the signal audit. M62 `sw_gate_agreement` scores manual against automated gates per event.
3. M03 `sw_fmo_cutoff`, M29 `sw_cluster_map_outliers`, M59 `sw_qc_breakdown` and M72 `sw_marker_resolution` are small and need no new Imports. M59 needs one change to `sw_qc_batch()`.
4. Next come five value-4 items at effort M (M65, M46, M73, M67, M20) and three value-3 items with fit 5 (M04, M30, M61). M66 (quantile spline engine) is later: effort L, and it changes the corrector contract.
5. Every item is a clean rewrite from the idea, except M63, which would wrap `changepoint` from Suggests. Every impactful value is a required argument. Section 4 lists what you must approve first.

Scope note: ROADMAP.md freezes new batch-side work. M52, M53, M54, M66, M67 and M41 touch that side and need an approved plan row first.

## 2. Tier table

Licence action is `clean-rewrite-from-idea` unless stated. Effort: S small, M medium, L large. Value and fit are 1 to 5. Repo URLs are GitHub paths.

### do_first

| id | sw_name | Source repo and file:lines | What it does | Gap evidence | Effort | Fit | Value |
|---|---|---|---|---|---|---|---|
| M60 | `sw_signal_segment_drift` | https://github.com/jmeskas/flowCut R/flowCut.R:535-553, 577-581, 946-1000 (Artistic-2.0, idea only) | Per file and channel: range of segment means and largest adjacent jump, both over the 2-98 percentile range, plus worst channel. No cutoff. | `sw_fcs_signal_audit()` has only `time_trend`, `rate_cv`, `max_gap` (R/fcs_signal_audit.R:330-342). No segment-drift hit in R/, NAMESPACE or ROADMAP.md. | S | 5 | 4 |
| M62 | `sw_gate_agreement` | https://github.com/RGLab/flowWorkspace R/getStats.R:331-366 (AGPL, idea only); https://github.com/gfinak/OpenCytoTutorial OpenCytoTutorial.Rmd:306-338 (CC BY-NC, idea only) | Per sample and population: event-level Jaccard, precision, recall of automated vs manual gates, plus per-population count correlation. No pooled correlation, no 0.80 target. | Closest are `sw_gate_qc()` (geometry shift only) and `sw_faust_compare()` (ARI). No gate-agreement function in NAMESPACE or ROADMAP.md. | M | 5 | 4 |
| M03 | `sw_fmo_cutoff` | https://github.com/aef1004/cytotypr R/get_99.R:13-20; R/filter_FMO.R:13-24; R/fe.R:21-29 (MIT) | Per-marker cutoff at a quantile of the FMO tube, exact-name map, labels 0/1 per stained event, exceedance shares. | No FMO code in R/, NAMESPACE or ROADMAP.md. `sw_flag_iqr(thresholds=)` consumes analyst cutoffs. ROADMAP B10 is a parked, different positivity idea. | S | 4 | 4 |
| M29 | `sw_cluster_map_outliers` | https://github.com/saeyslab/FlowSOM R/3_buildMST.R:165-246, 297-400 (GPL >= 2, idea only) | Per-cell distance to the assigned SOM node vs a per-node median + k x MAD limit from reference events. Per cell, per node, per file tables. | `sw_cluster_predict()` calls `kohonen::map` (R/cluster.R:379) and drops `$distances`. No outlier function. | S | 4 | 4 |
| M59 | `sw_qc_breakdown` | https://github.com/saeyslab/PeacoQC R/PeacoQC.R:425-484; R/PeacoQC_helper_functions.R:491-557, 648-679 (GPL >= 3, reads result fields only) | Table of PeacoQC step split (IT, MAD, consecutive), per-channel MAD contribution, split and trending channels, events per bin. | `sw_qc_batch()` keeps only a 5-column summary and drops the PeacoQC object (R/qc.R:170-195). No hit for these fields. | S | 4 | 4 |
| M72 | `sw_marker_resolution` | https://github.com/hammerlab/t-cell-data Flow.Rmd:335-348 (no licence, idea only) | Per file, population, marker: (median_pop - median_rest) / (MAD_pop + MAD_rest), then a donor summary. NA with a reason when MAD sum is zero. | No resolution or separation function. Nearest is `sw_cluster_scores()` (per-file robust z). | S | 4 | 4 |

### next

| id | sw_name | Source repo and file:lines | What it does | Gap evidence | Effort | Fit | Value |
|---|---|---|---|---|---|---|---|
| M65 | `sw_unmix_secondary_index` | https://github.com/Olsman/CompensAID R/CalculateSSI.R:12-32; R/GetSegment.R:11-59; R/UpdateMatrixInfo.R:97-101 (GPL >= 3, idea only) | Secondary-over-primary stain index in segments of the positive primary, per declared pair, on a stained sample. | No SSI code in R/, NAMESPACE or ROADMAP.md. `sw_unmix_negative_shape` needs paired setups. | M | 4 | 4 |
| M46 | `sw_unmix_control_residual` | https://github.com/camellialambert/CARROT pipeline_cytometrie.R:484-542 (no licence, idea only) | Per single-stain control and channel: positive and negative medians before and after unmixing, and their difference. | No pos-minus-neg matrix per control. ROADMAP.md:155 records the gap. | M | 4 | 4 |
| M73 | `sw_unmix_residual` | https://github.com/UMGCCCFCSR/CytometryInR course/14_Unmixing/index.qmd:936-937, 1005-1029, 1128-1196 (AGPL-3, idea only) | Per-detector RMS residual of events fitted to a signature set; compares full panel, drop-each and add-extra sets. | No residual getter in R/ or NAMESPACE. AutoSpectral 1.8.8 not checked. | M | 4 | 4 |
| M67 | `sw_correct_node_audit` | https://github.com/biosurf/cyCombine R/02_batch_correct.R:446-531; R/utils_helper.R:108-141 (MIT; helper adapted from sva) | Per SOM node: batches present, skipped or thin, covariate and anchor status, values at the input bound. | `sw_design_audit()` checks confounding globally only (R/design_audit.R:421-427). | M | 4 | 4 |
| M20 | `sw_cofactor_scan` | https://github.com/mBiocoder/Nature_immunology_code_repository CD8_Tm_FSFC/Preprocessing_functions.R:79-117 (MIT) | Per channel and cofactor grid: bimodality, Bartlett and negative-spread curves. No pick. | Only a single required cofactor exists (R/clean_gate.R:144-157). `.sw_bimodality()` (R/fcs_signal_audit.R:489-502) is reusable. | M | 4 | 4 |
| M04 | `sw_gate_control_cut` | https://github.com/SATVILab/DataTidyMusvosviPredictTBTASA script:181-190, 256 (MIT); https://github.com/mehrnoushmalek/flowDensity R/openCyto_wrappers.R:~50-135; https://github.com/RGLab/cytogallery neg_control.Rmd:27-52 | Pairs each tube group with its control; cut from `sw_gate_density_cut()` on the control, applied to the group. | `sw_gate_density_cut/apply` exist but have no control pairing. | S | 5 | 3 |
| M30 | `sw_fcs_compensation_state` | https://github.com/wguesdon/everything-flow-cytometry R/panels.R:132-162; R/harmonisation.R:213-235 (no licence found) | Header-only state per file: no matrix, identity, matrix present, already compensated. | `sw_fcs_header_audit()` does not parse the matrix or apply flag. | S | 5 | 3 |
| M61 | `sw_gate_density_cut(method = "tail_slope")` | https://github.com/mehrnoushmalek/flowDensity R/helper_functions.R:611-709, 90-122 (Artistic-2.0) | Adds a one-peak tail-slope cut to the existing cut function. | `sw_gate_density_cut` has only valley, percentile, mad (R/gate_density.R:45-75). | M | 5 | 3 |

### later

Ordered by value, fit, then effort. Items to fold into other work are marked.

| id | sw_name | Source repo and file:lines | What it does | Gap evidence | Effort | Fit | Value |
|---|---|---|---|---|---|---|---|
| M66 | `sw_corrector_engine_quantile_spline` | https://github.com/saeyslab/CytoNorm R/QuantileNorm.R:308-410, 538-547 (GPL >= 2) | Per-population monotone quantile spline engine. Needs a contract change (nonlinear maps). | No spline code; registry holds paired_shift and two heb engines. | L | 3 | 4 |
| M18 | `sw_gate_threshold_sweep` | https://github.com/DavideCapozzi/pbmc-fcs-batch-correction scripts/diagnose_fcs.R:298-355 (no licence); https://github.com/RGLab/paper-opencyto generate-GatingTemplate-HVTN065.r:24-60 | Share above each candidate cut over a required grid. | One cut per call today. | S | 4 | 3 |
| M19 | `sw_qc_sweep` | https://github.com/HdBraanker/Spectral_Flow_Workflow script.Rmd:493-505 (GPL-3); pbmc repo scripts/diagnose_fcs.R:185-222 | PeacoQC over an IT_limit x MAD grid; percent removed per sample. | R/clean_gate.R:118-143 tells users to loop by hand. | S | 4 | 3 |
| M22 | `sw_fcs_integrity_audit` | https://github.com/mBiocoder/Nature_immunology_code_repository CD8_Tm_FSFC/Compare.R:104-115 (MIT); https://github.com/codeneeded/NK_Manuscript FLOWSOM/Scripts/QC_TARA.R:53-66 | `$TOT` vs rows, time span, unreadable files, channel set. Never edits `$TOT`. | No `$TOT`/`$ETIM` code in R/. Time-span part rests on an unopened file. | S | 4 | 3 |
| M23 | `sw_pop_control_subtract` | https://github.com/RGLab/paper-opencyto lib/helpers.R:114-137 (no licence) | Test minus control per pair key; stops when a control is missing. | No stimulation or control subtraction. | S | 4 | 3 |
| M32 | `sw_count_floor` | https://github.com/wguesdon/everything-flow-cytometry R/counting.R:32-117 (no licence) | Poisson or binomial CV floor, events needed, excess CV. May fit better as columns in `sw_ref_composition()`. | `sw_ref_composition()` has the floor but not events-needed. | S | 4 | 3 |
| M10 | `sw_replicate_spread` | https://github.com/binhe-lab/E013-Pho4-evolution 03-Pho4-chimera/03-flow-full-data-lfs/analysis/20231229-pool-qc-PHO5.Rmd:307-336, 441-523 (no root licence) | n, mean, CV or SD of per-file medians per replicate group, with leave-one-out. | No replicate CV function. | S | 4 | 3 |
| M12 | `sw_plate_position_scores` | same repo, same Rmd:~338-440 (plots only in source) | Spearman vs column, eta2 vs row per run and control. | No well parsing. `order_rho` covers drift only. | S | 4 | 3 |
| M24 | `sw_retention_association` | https://github.com/maryellenlynall/2019_bcell_stress bcellstress07.Rmd:150-197 (no licence) | Events kept vs batch and condition as effect sizes. | No retention-vs-design function. | S | 4 | 3 |
| M05 | `sw_gate_membership` | https://github.com/SoegaardLab/2026_Fisher_Garcia_Frattari_Naasz_et_al spectral_flow_ics/0_source.R:297-370 (no licence) | Per-event gate-membership matrix from a GatingSet with alignment checks. | Only single-node internal helper (R/unmix_spike_helpers.R:8). | S | 4 | 3 |
| M08 | `sw_control_fpr_limit` | https://github.com/chuongj/ee_gap1_archmuts newGates.R:652-690 (no licence) | Limit from control gate frequencies (mean+SD or median+IQR, rule required). | No limit-from-controls function. | S | 4 | 3 |
| M34 | `sw_replicate_distance` | https://github.com/saeyslab/Liechti_et_al_TwinStudy 210708_figures.R:251-310 (GPL-3) | Within-unit distance vs a label-shuffled null. | `sw_ref_noise_floor()` has no label permutation. | S | 4 | 3 |
| M38 | `sw_fcs_marker_dye_map` | pbmc-fcs-batch-correction scripts/diagnose_fcs_panels.R:143-375 (no licence) | Marker x batch dye table from headers; unsplit names reported. | No dye-by-batch view. Feeds ROADMAP D5. | S | 4 | 3 |
| M39 | `sw_marker_names` | pbmc-fcs-batch-correction src/functions/fcs_io.R:20-87 (no licence) | Exact alias map, collision stop, shared marker set. Better built with ROADMAP D5. | `sw_channel_set_markers()` is single-frame, warns only. | S | 4 | 3 |
| M40 | `sw_fcs_header_signatures` | https://github.com/FredHutch/immune-flow-gating bin/process.R:28-112, 166-241 (GPL-3); https://github.com/paiva-s-lab/FlowCT R/unify.FCSheaders.R:~41-53 | Group files by header signature with counts; no auto reference. | `sw_fcs_header_audit()` needs a named reference. | S | 4 | 3 |
| M52 | `sw_correct_metric_anchors` | https://github.com/openproblems-bio/task_cyto_batch_integration src/control_methods/utils.py:3-16 (MIT) | Seeded shuffle candidate for `sw_correct_metrics()`. Could be a doc example. | Identity anchor already exists via `identity_id`. | S | 4 | 3 |
| M58 | `sw_anomaly_feature_stability` | https://github.com/AUMC-HEMA/CytoScan-manuscript 02-Liechti_28-experiments.Rmd:~112-180 (no licence) | Feature agreement across an n_events grid with seeded draws. | Only `sw_cluster_stability()` has this pattern. | S | 4 | 3 |
| M70 | `sw_cluster_label_agreement` | https://github.com/saeyslab/FlowSOM R/2_buildSOM.R:332-362; R/4_metaClustering.R:199-215 (GPL >= 2) | Cluster purity and per-label F vs gate labels. | Only `.sw_ari` in R/faust_wrap.R. | S | 4 | 3 |
| M02 | `sw_bool_combinations` | https://github.com/SoegaardLab/2026_Fisher_Garcia_Frattari_Naasz_et_al spectral_flow_ics/0_source.R:532-595, 652-682; https://github.com/aef1004/cytotypr R/fe.R, count_calc.R, identified_pop_perc.R | Combination labels, zero fill, optional background subtraction. Absorbs `sw_phenotype_share`. | No Boolean combination code. | M | 4 | 3 |
| M09 | `sw_scatter_dependence` | https://github.com/binhe-lab/E013-Pho4-evolution ...20230216-gate-normalize-data.Rmd:235-294 | Log-log slope of marker on scatter per file. Source draws a loess only. Defer normalise step. | No marker-vs-scatter function. | M | 4 | 3 |
| M17 | `sw_fcs_harmonize` | https://github.com/saeyslab/BonteEtAl_AML_Dx_EuroFlow Rscripts/0_HelperFunctions.R:13-28; Rscripts/1_DataPreprocessing.R:44-197 (no licence) | Logged channel, marker and spillover renames. | Header audit changes nothing. | M | 4 | 3 |
| M35 | `sw_control_peak_check` | https://github.com/wguesdon/everything-flow-cytometry R/spillover_compute.R:505-588 | Positive share and brightest detector per control. No verdict. | `sw_unmix_check_controls()` is file-level only. | M | 4 | 3 |
| M41 | `sw_cluster_match` | pbmc-fcs-batch-correction src/functions/cluster_matching.R:16-219; task_cyto_batch_integration src/metrics/abundance_preservation/helper.R:58-486 | Match clusters across batches; per-pair marker deltas. | Nothing matches clusters across batches. | M | 4 | 3 |
| M42 | `sw_gate_impute` | https://github.com/RGLab/HIPCCyto R/impute.R:77-198 (GPL-3) | Replace failed-sample gates by nearest donor. Ship nearest only. | No donor-selection function. | M | 4 | 3 |
| M47 | `sw_gate_qc` new checks | https://github.com/RGLab/HIPCCyto R/save.R:570-672 (GPL-3) | Gate centre d2, size ratio, in-gate dip p. Needs `diptest` in Suggests. | `sw_gate_qc()` compares one reference sample only. | M | 4 | 3 |
| M49 | `sw_fcs_rate_sections` | https://github.com/prybakowska/CytoQP R/clean_flow_rate_and_signal.R:427-836; https://github.com/giannimonaco/flowAI R/anomaly-detection.R:85-146; flowCut R/flowCut.R:810-840 | Anomalous time bins and event share; removes nothing. | Only file-level `rate_cv`, `max_gap`. | M | 4 | 3 |
| M51 | `sw_qc_removal_audit` | https://github.com/1ordinateur/FlowMOP_paper benchmarks/analyse_fluorescence_instability.R:23-205 | Removed vs kept divergence per time bin, circular-shift null. | No fluorescence-based QC audit. | M | 4 | 3 |
| M53 | `sw_correct_peak_consistency` | task_cyto_batch_integration src/metrics/ratio_consistent_peaks/helper.R:6-101; script.py:63-226 (MIT) | Peak-count case table before vs after correction across two splits. | No before/after peak count. | M | 4 | 3 |
| M54 | `sw_cluster_replicate_agreement` | task_cyto_batch_integration src/metrics/flowsom_mapping_similarity/helper.R:35-125; script.R:54-139 (MIT) | Cross-replicate label-composition difference, both directions. | `sw_cluster_stability()` is across seeds only. | M | 4 | 3 |
| M63 | `sw_signal_changepoints` | flowAI R/signal.R:7-125; R/auto-qc.R:138, 254 (GPL >= 2) | Binary-segmentation changepoints in bin medians. Wrap `changepoint` from Suggests (the only item with that licence action). | No changepoint search. | M | 4 | 3 |
| M68 | `sw_signal_above_noise` | https://github.com/saeyslab/iidx 03_OutlierAndNoiseDetection.R:~152-232; 02_FlowSOMFeatureExtraction.R:349-385; 99_Aux.R:838-920 (Apache-2.0) | Cluster MFI vs unstained median + n x MAD per batch and marker. | No unstained-floor scoring. ROADMAP B10 is related. | M | 4 | 3 |
| M69 | `sw_marker_nrs` | https://github.com/HelenaLC/CATALYST R/utils-differential.R:73-83; R/plotNRS.R:43-95 (GPL >= 2) | Per-sample marker non-redundancy score from PCA loadings. | No NRS or per-sample loading score. | S | 4 | 3 |
| M01 | `sw_stim_response` | https://github.com/SATVILab/DataTidyMusvosviPredictTBTASA Script to perform TB-TASA.R:~307-331 (MIT) | Per-tube Fisher test and fold change vs control. | None; `sw_diff_abundance()` is across donors. Not on ROADMAP. | M | 3 | 3 |
| M25 | `sw_dimred_project` | Spectral_Flow_Workflow script.Rmd:235-390 (GPL-3) | UMAP fit on train files, held-out projection. | `sw_dimred_run()` has `ret_model = FALSE`. | M | 3 | 3 |
| M28 | `sw_cluster_overlap_scores` | https://github.com/prybakowska/CyTOF_analysis_Pipeline1 functions.R:670-930; https://github.com/prybakowska/CytoQP R/files_quality_check.R:352-547 | AOF-like overlap per marker. The cytutils formula is unread. | No AOF code. | M | 3 | 3 |
| M31 | `sw_gate_cut_consensus` | https://github.com/mehrnoushmalek/DuraClone-gating HelperFunctions.R:5-31 (academic licence) | Flag or replace outlier cuts. Default should be flag. | Cuts table exists; no cross-sample step. | S | 3 | 3 |
| M33 | `sw_variance_share` | everything-flow-cytometry R/harmonisation.R:~380-418 | Donor variance share via `lme4`. ICC was out of scope in plans/2026-11-11_reference-and-batch-diagnostics.md:204. | `eta2` only. | M | 3 | 3 |
| M37 | `sw_diff_abundance_glmm` | https://github.com/masashi-CU/ISC-flow-cytometry functions/masc.R:33-139 (no licence); iidx 04a_TestsDA.R:28-426 (Apache-2.0) | Per-cluster GLMM with LRT; NA with reason on failure. | `sw_diff_abundance()` is edgeR/voom only. | M | 3 | 3 |
| M07 | `sw_gate_control_check` | https://github.com/chuongj/ee_gap1_archmuts newGates.R:286-314 (no licence; repo owner is chuongj) | Control share in expected gate over `sw_gate_qc()` output. Thin layer; may be an argument. | `sw_gate_qc()` has no control concept. | S | 3 | 3 |
| M11 | `sw_sample_ledger` | binhe-lab E013-Pho4-evolution ...20231229-pool-qc-PHO5.Rmd:52, 117-120, 243, 258, 275 | Reasoned exclusion log. Source stores no reasons. | Event-level audit only. | S | 3 | 2 |
| M21 | `sw_fcs_panel_match` | BonteEtAl_AML_Dx_EuroFlow Rscripts/0_HelperFunctions.R:~30-36; 1_DataPreprocessing.R:92-197 | Marker-set match per file. May be a multi-reference mode of the header audit. | No panel-set matching. | S | 4 | 2 |
| M26 | `sw_cluster_control_cv` | https://github.com/saeyslab/CytoNorm R/evaluation.R:150-215; Spectral_Flow_Workflow script.Rmd:~161-181 | CV of cluster shares across control files over a k grid (kmeans rewrite). | `sw_cluster_scores()` has no CV or k sweep. | S | 4 | 2 |
| M48 | `sw_io_event_budget` | https://github.com/peterleary/MARMOT inst/pipeline/MARMOT_Pipeline.qmd:636-737 (MIT) | Per-sample event cap from a declared RAM budget. | `sw_io_subsample()` takes a count. | S | 3 | 2 |
| M56 | `sw_qc_mask_validate` | https://github.com/saeyslab/PeacoQC_evaluation R_Scripts/FlowCAP4/6_Quantitative_analysis.R:9-117 (no licence) | Sensitivity and specificity vs bad-event truth. Needs event-level labels. | No QC-mask scorer. | S | 4 | 2 |
| M55 | `sw_gate_compare_proportions` | https://github.com/RGLab/flowcap3 reports/SeraCareComparison/GatingComparison.Rmd:79-195 (no licence) | Mixed-model method comparison on logit proportions. | None; ROADMAP C6 mentions lme4 only. | L | 3 | 2 |

## 3. Sketches for do_first and next

Convention for all: no defaults for listed arguments, caller transforms (cofactor recorded, never applied), no pass or fail, thin rows are `not_assessed`, dropped rows are reported in an attribute. Each plugs into `sw_diagnose()` as a module only if you approve it; that is a separate plan row.

**M60 `sw_signal_segment_drift`**
- Signature: `(data, channels, time_channel, events_per_segment, range_probs, statistic, tail, min_segments, min_events)`. Required choices: all nine.
- Output: long tibble `file, channel, check (segment_range_frac | segment_max_jump_frac), value, n_segments, n_events_used, n_dropped, status`, plus worst-channel row.
- Plugs in: standalone, stackable in `sw_diagnose()`. Use first-pass (pre-removal) numbers, not flowCut's post-cleaning values. `tail` chooses fold-into-last-segment (flowCut) or drop.
- First test: a synthetic file with a step in one channel gives a jump near step/range at the right segment, and a flat channel gives about 0.

**M62 `sw_gate_agreement`**
- Signature: `(manual_member, auto_member, name_map, reference, min_events, rare_events, min_samples, cor_method, statistic)`. Required: all.
- Output: sample x population `n_manual, n_auto, n_both, jaccard, precision, recall, rare, status`; per-population correlation rows.
- Plugs in: after `sw_gate_apply()` on a shared scale; validates `sw_clean_run` and `sw_gate_transfer` against `sw_gate_read_flowjo` gates. Stops on unmapped names or event-count mismatch.
- First test: identical membership matrices give Jaccard 1; a swapped name map errors.

**M03 `sw_fmo_cutoff`**
- Signature: `(fmo_data, stained_data, fmo_map, probs, min_events, scale, cofactor, parent_population, pooling, batch)`. Required: all but `batch` (NULL means one pooled cutoff, stated).
- Output: `list(cutoffs, labels, dropped, exceedance)`. Exact-name map; stops on missing or duplicate match or a thin FMO tube.
- Plugs in: after gating, before `sw_flag_iqr(thresholds=)`. FMO tubes need a manifest role or map table.
- First test: map CD4 and CD44 tubes; assert no double match and that the minimum stained event is not NA.

**M29 `sw_cluster_map_outliers`**
- Signature: `(cluster_result, newdata, reference_data, mad_allowed, min_events_per_node, by, seed)`. Required: all.
- Output: `list(per_cell, per_node, per_file)`; events in `not_assessed` nodes counted separately.
- Plugs in: after `sw_cluster_run()` (SOM only); uses `kohonen::map()$distances`. Reference events should be held-out files.
- First test: held-out events from a shifted marker produce more outliers than unshifted ones; a node with fewer than `min_events_per_node` reference events gives NA.

**M59 `sw_qc_breakdown`**
- Signature: `(qc_results)`; no impactful argument. Reads fields PeacoQC already computes.
- Output: `per_sample` (IT, MAD, consecutive percent, status) and `channels` tibbles. Missing `ITPercentage` gives NA plus status.
- Plugs in: needs `sw_qc_batch()` to keep `qc_result$PeacoQC` per sample (existing function changes; needs its own test).
- First test: a MAD-only object gives `it_status = "not_requested"`; fewer than 150 bins gives `skipped_too_few_bins`.

**M72 `sw_marker_resolution`**
- Signature: `(x, file_id, population, markers, reference, min_events, n_events, seed, min_files, mad_constant, transform_note, label_source)`. Required: all.
- Output: `per_file`, `summary`, `removed`. Zero or non-finite MAD sum gives NA, status `mad_zero`.
- Plugs in: before and after unmixing, relativisation or correction on the same populations. Cluster labels cap the claim, as in `sw_batch_scores`.
- First test: two well-separated populations give a large positive value; identical medians give 0; zero MAD gives NA.

**M65 `sw_unmix_secondary_index`**
- Signature: `(data, markers, pairs, cutoffs, negative_gap, n_segments, min_events, neg_spread, segment_scheme, merge_thin)`. Required: all.
- Output: long tibble per sample, pair, segment with `ssi, status`. Returns all segments (upstream keeps only the highest merged one).
- Plugs in: after unmixing, on stained samples; stackable in `sw_diagnose()`. Doc must say a rising secondary can be true co-expression.
- First test: a planted linear spill of 0.05 gives SSI rising with segment; a thin segment gives NA, not a merge, when `merge_thin = FALSE`.

**M46 `sw_unmix_control_residual`**
- Signature: `(before, after, controls, pos_gate, neg_gate, channels, scale, cofactor, min_events, centre)`. Required: all (cofactor NULL only for linear).
- Output: control, stage, channel, `n_pos, n_neg, pos_median, neg_median, pos_minus_neg`; diagonal NA. Source used compensation, not unmixing; "after" must be unmixed data.
- First test: a control with planted off-diagonal spill shows it before and near 0 after.

**M73 `sw_unmix_residual`**
- Signature: `(x, signatures, detectors, weights, n_events, seed, scale, clip_negative, compare, extra_signatures, event_filter)`. Required: all.
- Output: `file, signature_set, detector, rms_residual, n_events_used, n_events_dropped` plus delta vs the full panel. Reuse `sw_unmix_sample_events()` for sampling. Check first whether AutoSpectral 1.8.8 exports an equivalent (call, do not rebuild).
- First test: drop one true signature from a synthetic mix; residual rises at that signature's peak detector.

**M67 `sw_correct_node_audit`**
- Signature: `(df, label, batch_col, covar, anchor, min_cells_per_batch, skew_margin, covar_confounded, max_uncorrected_share, corrected)`. Required: all (NULL for covar, anchor, corrected, max_uncorrected_share is explicit).
- Output: one row per node with action and covariate status. Runs between `sw_correct_som()` and `sw_correct_apply()`; bound counts are inferred, not exact.
- First test: a node with one batch gives `single_batch_skipped`; a perfectly confounded covariate gives `dropped_confounded`.

**M20 `sw_cofactor_scan`**
- Signature: `(x, channels, cofactors, curves, n_events, seed, min_events, by_file)`. Required: all. Raw linear input only. Source `b` is a multiplier; cofactor = 1/b.
- Output: long tibble `file, channel, cofactor, curve, value, n_events, status`. Reuses `.sw_bimodality()`; GMM negative needs mclust (Suggests). Read flowVS source before the Bartlett curve.
- First test: a mixture with a known negative gives a bimodality curve peaking near the true scale; a non-positive cofactor stops.

**M04 `sw_gate_control_cut`**
- Signature: `(data, manifest, sample_col, group_col, control_col, control_value, channels, method, bandwidth, min_peak_height, min_peak_distance, prob, k, mad_side, min_control_events, on_missing_control)`. Plus `pool_controls`. Required: all.
- Output: long tibble of cut, control false-positive share, per-sample share above. A group with no control stops and names the group.
- First test: group without a control tube errors by name.

**M30 `sw_fcs_compensation_state`**
- Signature: `(paths, matrix_keywords, apply_keywords, zero_tol, expected_state)`. Required: all.
- Output: file, keyword used, matrix size, largest off-diagonal, state. Header only; never applies a matrix. Read keyword text and parse it; stop on unparsable text.
- First test: identity `$SPILLOVER` with `APPLY COMPENSATION TRUE` gives a stated state, with both facts in separate columns.

**M61 `sw_gate_density_cut(method = "tail_slope")`**
- New required arguments: `side, slope_alpha, height_cap, slope_half_window, step, min_events, min_peak_height, bandwidth`. Returns NA with a reason, never Inf. Index the tested grid point directly (flowDensity maps back with the wrong step).
- First test: a one-peak density with a known tail returns a cut right of the peak; fewer than `min_events` returns NA with a reason.

## 4. Needs your decision

Each row is an impactful choice to approve before building. None may become a default. Choices an earlier stage made (marked inherited) need a one-line confirm. Later-tier items list their choices in the verified candidate records; they are asked for when an item is promoted.

| Item | Choice |
|---|---|
| M60 | events_per_segment (flowCut 500) |
| M60 | range_probs (flowCut 0.02, 0.98) |
| M60 | statistic per segment (mean vs median) |
| M60 | channels included |
| M60 | scale of input and cofactor (inherited: 6000) |
| M60 | min_segments |
| M60 | time ordering: time_channel vs row order |
| M60 | tail handling: fold vs drop |
| M62 | reference side (manual or auto) |
| M62 | name_map between manual and auto populations |
| M62 | min_events per sample x population |
| M62 | rare_events threshold |
| M62 | min_samples for correlation |
| M62 | correlation type (Pearson or Spearman) |
| M62 | statistic returned (count, freq, events) |
| M62 | data scale the manual gates were drawn on |
| M03 | probs (cytotypr 0.99) |
| M03 | FMO-tube to marker map |
| M03 | scale and cofactor before the quantile |
| M03 | min events in the FMO tube |
| M03 | parent population |
| M03 | pooled vs per-batch cutoff |
| M03 | tie rule (> vs >=, include.lowest) |
| M29 | mad_allowed (FlowSOM 4) |
| M29 | reference events: training vs held-out files |
| M29 | min events per node |
| M29 | marker set and scale (kohonen does not scale) |
| M29 | distance metric (Euclidean) |
| M29 | spread statistic (median + MAD vs quantile) |
| M59 | none (reads PeacoQC fields). Confirm the `sw_qc_batch()` change to keep the PeacoQC object. |
| M72 | transform and cofactor (inherited) |
| M72 | population labels and label_source |
| M72 | reference group: rest or named set |
| M72 | min_events per file x population |
| M72 | unit: per file then donor summary |
| M72 | n_events and seed |
| M72 | mad_constant (1.4826 vs 1) |
| M72 | MAD floor for near-empty negatives (possible new value) |
| M72 | min_files for the donor summary |
| M65 | n_segments (CompensAID 4) |
| M65 | min events per segment and negative (50) |
| M65 | positivity cutoff per primary channel |
| M65 | negative_gap band |
| M65 | neg_spread: sd vs MAD |
| M65 | segment_scheme: equal width vs equal count |
| M65 | declared pair table |
| M65 | thin segments: merge vs NA |
| M65 | transform and cofactor (inherited) |
| M46 | pos and neg gate per control |
| M46 | scale and cofactor |
| M46 | min_events per gate |
| M46 | centre statistic |
| M46 | channel-to-fluorophore map after unmixing |
| M46 | what "before" is (raw detectors or prior setup) |
| M73 | weights: OLS vs weighted, matching the real unmixing |
| M73 | detector set and event gate |
| M73 | n_events and seed |
| M73 | residual scale: raw vs relative |
| M73 | signature sets to compare |
| M73 | clip_negative |
| M73 | signature normalisation and intercept |
| M67 | min_cells_per_batch (new threshold) |
| M67 | skew_margin (cyCombine hard-codes 5) |
| M67 | covar_confounded: record or stop |
| M67 | max_uncorrected_share (NULL = no stop) |
| M67 | SOM grid and seed (inherited: ROADMAP.md:245 lists 8, 8, 42 as unconfirmed) |
| M20 | cofactor grid per channel |
| M20 | channels |
| M20 | n_events and seed |
| M20 | pooled vs per-file scan |
| M20 | min_events |
| M20 | which curves |
| M20 | negative-peak definition (GMM model and components) |
| M20 | Bartlett binning |
| M04 | control_col and control_value |
| M04 | group_col |
| M04 | cut method and its knobs |
| M04 | min_control_events |
| M04 | one vs pooled control tubes |
| M04 | channels to gate and scale |
| M30 | zero_tol (source: exactly 0) |
| M30 | matrix_keywords and apply_keywords |
| M30 | expected_state per file |
| M30 | precedence when identity matrix and apply flag TRUE |
| M61 | bandwidth, side, slope_alpha, height_cap |
| M61 | slope window and step |
| M61 | min_peak_height and min_events (flowDensity hard-codes 20) |
| M61 | fallback to inflection or none |
| All | Whether to take on non-spectral assay utilities at all (M01, M23) and batch-side items despite the ROADMAP freeze |
| All | Whether each new check also gets a `sw_diagnose()` module (adds a settings list per module) |
| All | Plot-first rule: each new function that feeds a conclusion needs a planned, inspected figure (M20 and M46 especially) |

## 5. Rejected or already covered

| Idea | Reason |
|---|---|
| M57 anomaly spike recall | Gap is narrow. `sw_synref_shift()` plants shift, stretch_median and two_mode. Repo benchmarks `04_injected_faults.R`, `08_detection-limits` and `06_anomaly-features` already score flags. Missing: a depletion family and a TP/FP/FN tally at a cutoff, both small additions to existing code. Source stretch multiplies about zero, which is wrong on an asinh scale. Source separate effect needs `flowDensity::deGate`. |
| PIDgeon outlier share (second half of M29) | Same rule as FlowSOM; its `outliers[, 4]` indexing fails on current FlowSOM. Merged into M29. |
| `sw_phenotype_share` | Subset of M02 (no background). Merged. |
| `sw_gate_concordance` | Same idea as M62. Merged. |
| Gate fit-on-control names (three claims) | One idea. Merged into M04. |
| CytoNorm `testCV` wrap (M26 claim 1) | Needs a FlowSOM object and FlowSOM/CytoNorm are not in Suggests. Replaced by a kmeans rewrite (M26, later). |
| Source-copied verdict fields | Not ported anywhere: responder flag defaults, `delta_pass`, `recommended_threshold`, `it_recoverable`, the 0.80 agreement target, the 1.5 or 2 CV cutoff. They break Q14 (no pass or fail) or are unplanned impactful choices. |
| Source behaviours not to copy | Silent NA-to-0 or NA-to-1, substring name matching, first-match tie rules, overwriting `$TOT` or PnS, unseeded draws. |
| M52 perfect-integration anchor | The source "perfect" control uses the same split on both sides; it has no meaning for before/after metrics. `sw_ref_noise_floor()` already gives the ceiling idea. Only the seeded shuffle part stays (later). |
| Separate dim-red stability metric (M25) | Source scan is visual only. A stability metric would be new design, not a port. |

## 6. Repos that could not be read

None failed to open. Reading limits that affect trust:

- Nothing was executed; all sources were read only.
- `GreshamLab/ee_gap1_archmuts` does not exist. The real repo is `chuongj/ee_gap1_archmuts` (M07, M08).
- Not opened or unverified: `lynall/bcellstress06.Rmd` (M17), `ismmshimc/cytutils` AOF formula and licence (M28), `nacl` `compare_panels` (M21), the ISC `try/catch` block (M22), `saeyslab/ComputationalCytometry` (M26), and the saeys `helperFunctions.R` time-span check (M22).
- M73: R was not on PATH, so AutoSpectral 1.8.8 exports could not be checked.
- M20: the flowVS Bartlett method was not read; the Bartlett and negative-peak curves are new design.
- M65, M72 and several others rely on repos with no licence file, so only the idea is reused.


## 7. Critic follow-ups

- Differential, report and app stages have almost no proposals. Differential has only M37 (GLMM, value 3, 'later'). No proposal covers sw_report or sw_qc_app or sw_unmix_app/dashboard features. Catalog repos iidx (D22, differential and app layout), CATALYST (D13, only M69 NRS), cyCONDOR (B13) and FlowCT (B10, only header unify) were mined for one idea each at most.
- Catalog repos with no proposal citing them: B02 FredHutch/nf-flow-cytometry-qc (explicitly requested in the catalog, 'explicitly requested' tag), B03 DrCytometer/ImmunePhenotyping (same author as AutoSpectral, the repo's gold standard), B06 ppreshant/flow_cytometry, B07 flowFun, B13 cyCONDOR, C03 GatingMethod_evalutation, C05 CytoNorm_Figures, C07 faust_flowcap4_analysis, C09 UNITO, D05/D06 CytoExploreR and CytoGate, D14 FlowSOM_protocol, E09 cytoverse-bioc-2023, E10 R_flowcytometry_course, and A04, A05, A06, A07, A10, A11, A20, A24, A25, A26. Section 6 says 'None failed to open', but the doc never states which repos were skipped, so silent non-coverage is not reported. B02 and B03 are the biggest omissions given they were explicitly requested or are closest to AutoSpectral.
- The verify cap dropped several items that map directly to existing ROADMAP rows or CLAUDE.md fit, so the cap looks arbitrary. M45 (shape-only JSON run log) matches ROADMAP D6 (provenance block, 'planned', nothing calls sw_io_event_audit). M74 (unmixed FCS keyword consistency check) matches ROADMAP D7 (FCS write exporter) and the unmix stage. M43 (output hash baseline) is a refactor-safety tool for a repo with 139 exports. M27 (5th-percentile robust z for over-compensation) is an unmix-stage check with no unmix candidate besides M65/M46/M73. M75 (KDE-mode centring) is adjacent to sw_relativise_*. M71 (observed/expected cluster share, entropy) extends sw_cluster_scores. M64 (counting beads) is the only absolute-count idea. M14 duplicates and would improve M19.
- Unmix stage is thin and no unmix item is in do_first. Only M65, M46, M73 (all 'next') plus dropped M27 and M74. Nothing covers unmixing quality across controls beyond M35 (later). Check whether M73 duplicates ROADMAP B11/AutoSpectral 1.8.8 getters before spending effort; the proposal admits R was not on PATH, so this was never verified.
- Stage coverage is uneven in do_first: four of six items sit in QC and clustering; gate has M62 only; there is nothing for batch except frozen-side items, nothing for dimred, annotation (R/annotate.R) or pipeline orchestration (ROADMAP D8, D9). Annotation has zero proposals although cell-type annotation is a named pipeline step.
- Ideas that may have been verified away too quickly. M57 is rejected whole, yet the doc itself says two small additions are missing (a depletion family in sw_synref_shift and a TP/FP/FN tally at a cutoff); that is a concrete S-effort follow-up, not a rejection. M52 keeps only a doc example. M25 is rejected for lacking a stability metric, but sw_cluster_stability already has the pattern to reuse for a UMAP-projection stability check. M26's CytoNorm testCV idea was replaced by a kmeans rewrite and placed in 'later' at value 2 with no justification for the downgrade.
- Tiering conflicts with ROADMAP. M67 (batch node audit, SOM-based) is in 'next' with value 4 although batch work is frozen; the doc notes the freeze in a scope note but still ranks it above unfrozen items. M46 and M65 depend on paired or planted data that the roadmap says is not available (B3 to B10 parked until independent evaluation tubes exist). M62 depends on manual FlowJo gates and sw_gate_read_flowjo; confirm that reader exists and has tests before ranking it first.
- Section 4 lists more than 100 decision rows with no priority and no sense of which are already settled in AGENTS.md or ROADMAP.md's choice register (lines 188-260). Inherited choices (cofactor 6000, SOM 8/8/42) are tagged on some items but not checked against the register, so confirms may repeat work the user already did.
- Evidence is mostly from reading, nothing run. Several 'gap evidence' claims rest on grep of R/ and NAMESPACE only. M60 and M59 claims about flowCut and PeacoQC result fields were not executed against the installed versions. M20's Bartlett curve and M28's AOF formula were never read, so those two rows are new design, not ports.
- Licence handling is uneven. Many sources are 'no licence', which means idea only; yet the doc lists file:line ranges for them as if portable. M63 is the only wrap, but flowAI-derived M49 and M63 overlap and are listed separately.
- Open and triage the uncovered catalog repos in this order: B02 nf-flow-cytometry-qc, B03 ImmunePhenotyping, B13 cyCONDOR, C03 GatingMethod_evalutation, C07 faust_flowcap4_analysis, C09 UNITO, D14 FlowSOM_protocol, D05 CytoExploreR, B07 flowFun, B06 ppreshant. Add the explicit list of skipped repos to section 6.
- Add a differential and report and app survey: read iidx 04a_TestsDA.R and the iidx Shiny layout, cyCONDOR differential functions, and CATALYST diff-abundance wrappers; propose at most three items for the differential stage and one for sw_qc_app or sw_unmix_app.
- Re-verify the dropped candidates M45, M74, M43, M27, M75 and M71 with the same checks used for the others. Fold M45 into ROADMAP D6 and M74 into D7 as proposals; M43 as a repo test or experiments tool.
- Add a unmix-stage do_first candidate: check AutoSpectral 1.8.8 exports (R is needed on PATH; use pixi), then decide M73 call-or-build; pair M27 with M46.
- Add a small M57 follow-up: depletion family in sw_synref_shift plus a TP/FP/FN tally at a user cutoff, with the required-argument list.
- Add an annotation-stage survey (R/annotate.R vs cell-curator or celltypist-style scoring in cyCONDOR, FlowCT, CytoPipeline) and a cheap hash-baseline regression test.
- Cross-check Section 4 rows against ROADMAP.md lines 188-260 and any AGENTS.md, collapse repeated choices into one row per shared parameter (cofactor, min_events, seed, n_events), and tag each as settled, inherited-confirm or new.
- Move M67 and any other batch-side items into a separate 'frozen, needs approval' tier so 'next' only holds unfrozen work; reconsider M62 as the first build only after checking sw_gate_read_flowjo tests.

Run stats: {"surveyed": 10, "candidates": 95, "merged": 75, "verified": 60, "confirmed": 59, "unreachable": []}


---

# Addendum 2026-11-13: second pass on differential, annotation, gating, report and app stages

Source of this pass: 44 new candidates (ids N01 to N47) checked against the cloned repos and this repo. Nothing was executed; all sources were read only. Plan-step labels refer to the sections of `plans/2026-11-13_scenarios-and-step-reports.md`, not to ROADMAP.md rows (ROADMAP B3 is the parked study manifest and B7 is lock and audit).

## A1. Summary

- 39 items are confirmed with corrections and worth a plan row; 1 (N04) failed confirmation, 4 duplicate older items (N25, N27, N29, N42).
- Value 4 items are all fit 4: N13 (sample ordination), N14 (label transfer), N16 (hit card), N34 (population statistics), N35 (back-gating numbers), N36 (template pre-flight).
- No item has value 5. Most are small rewrites; 23 of the 39 are effort S.
- Every source with a GPL, CC BY-SA or no licence is idea-only. No code is copied. Apache-2.0 and MIT sources are also rewritten, so no new NOTICE entry is needed beyond the existing iidx one.
- Three items change an existing exported function or its tests: N15 (`sw_annotate_manual` stops instead of warns, breaks `tests/testthat/test-annotate.R:253-262`), N03 (`sw_qc_batch` names the failing sample), N32 (`sw_clean_run` keeps per-step event ids). Each needs your sign-off.

## A2. Ranked table

Ranked by value, then fit, then low effort. Licence column: `rewrite` = clean rewrite from the idea. Repo URLs: iidx = https://github.com/saeyslab/iidx; CytoExploreR = https://github.com/DillonHammill/CytoExploreR; cyCONDOR = https://github.com/lorenzobonaguro/cyCONDOR; CATALYST = https://github.com/HelenaLC/CATALYST; MARMOT = https://github.com/peterleary/MARMOT; FlowSOM_protocol = https://github.com/saeyslab/FlowSOM_protocol; cytoverse = https://github.com/OzetteTech/cytoverse-bioc-2023; B02 = https://github.com/FredHutch/nf-flow-cytometry-qc; B03 = https://github.com/DrCytometer/ImmunePhenotyping; B06 = https://github.com/ppreshant/flow_cytometry; EFC = https://github.com/wguesdon/everything-flow-cytometry; GME = https://github.com/hung-ching-chang/GatingMethod_evalutation; A04 = https://github.com/kauralasoo/flow_cytomtery_genetics; ISC = https://github.com/masashi-CU/ISC-flow-cytometry; Pepper = https://github.com/jackspepper/Pepper-2026; BEV = https://github.com/edfadeev/Marine_BEVs_across_South_Pacific; fibroblast = https://github.com/davismcc/fibroblast-clonality.

| Rank | id | sw_name | Plan step | Source and file:lines | What it does | Gap evidence | Effort | Licence action | Fit | Value |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | N13 | `sw_sample_ordination` | B7 report, also B2 | CATALYST `R/clrDR.R:111-140`, `R/pbMDS.R:70-93`; FlowSOM_protocol `R/FlowSOM_protocol.R:239-255`; MARMOT `inst/pipeline/MARMOT_Pipeline.qmd:2220-2243` (only calls `pbMDS`) | One point per sample from CLR cluster shares or marker medians, PCA or MDS; scores, loadings, variance tables | `R/dimred.R:126` is cell level; `R/batch_correct.R:628` is per batch; no CLR anywhere | S | rewrite (CATALYST GPL >= 2; FlowSOM_protocol no licence; MARMOT MIT) | 4 | 4 |
| 2 | N14 | `sw_annotate_transfer` | B7 annotate | cyCONDOR `R/data_projection.R:~100-181` | kNN label transfer with CV folds by file or donor; per-class recall, precision; `below_min_vote` NA. The source reports only accuracy by k, so vote share, recall and precision are new design | `R/annotate.R:75-110,266` labels clusters, no CV; no `class::` in DESCRIPTION | M | rewrite (GPL-3) | 4 | 4 |
| 3 | N16 | `sw_diff_hit_card` | B7 report | iidx `04a_TestsDA.R:205-262`, `06b_Plotting.R:284-400,828-930` | Donor points per group, coloured by batch, for each named hit; group effect per batch as a dot row. The per-batch table is our design | `sw_plot_abundance` (`R/differential.R:837`) has no donor points and no batch | M | rewrite (Apache-2.0) | 4 | 4 |
| 4 | N34 | `sw_pop_stats` | B3 report | CytoExploreR `R/cyto_stats_compute.R:196-206`, `R/cyto_stats-internal.R:113,357-364,448-455`; cytoverse `vignettes/7_Reporting.Rmd:58-82,164-167` | Long sample x node table: share of a named ancestor, linear median, CV; NA with reason | Only `gh_pop_get_count` at `R/clean_gate.R:778`; `R/gate_qc.R:15-26` is parent share per gate, no marker stats; the plan lists the table at `plans/2026-11-13_scenarios-and-step-reports.md:83` with no producing function | M | rewrite (CytoExploreR GPL-2); use `gh_pop_get_count`, already in Imports | 4 | 4 |
| 5 | N35 | `sw_gate_backgate` | B3 | CytoExploreR `R/cyto_plot_gating_scheme.R:183-195,399-433` | Edge share, distance bins and centroid of a descendant population against an ancestor gate, plus overlay. The source only draws; the numbers are our design | `R/gate_qc.R:76,105-107` measures parent events only; `R/gate_plot.R:52-60` has no population overlay | M | rewrite (GPL-2) | 4 | 4 |
| 6 | N36 | `sw_gate_template_check` | B3 | CytoExploreR `R/cyto_gate-helpers.R:519`, `R/cyto_gatingTemplate-helpers.R:614-640`; fibroblast `src/Rmd/facs_analysis_template.Rmd:129-130,202-203,275-277,354-358` | Pre-flight on the template CSV: dims resolve, alias unique, parent earlier, ref and bool gates exist, `groupBy` columns exist; explicit group-to-channel map | `R/gate.R:340-362` checks only `file.exists`, then openCyto fails late at `:424` | M | rewrite (GPL-2; CC BY-SA 4.0) | 4 | 4 |
| 7 | N02 | `sw_fcs_batch_date_check` | B1 | B03 `source/perform_batch_correction.R:30-48,79-82`, `source/preprocess_flow_files.R:67-76` | Crosstab of declared batch (and manifest `run`) against `$DATE`, `$BTIM`, `$CYTSN`; tile map. B03 builds batches from `$DATE` and never checks a declared batch | `R/fcs_header_audit.R:20-21,58` reads the keywords against one reference file, no manifest join; `R/manifest.R:9-20` | S | rewrite (GPL-3) | 4 | 3 |
| 8 | N07 | `.sw_report_readback_check`, `.sw_report_collect` | report, Part A | B02 `bin/summarize.py:67-97`, `bin/flow_qc_widget.py:932-987`; B06 `analyze_fcs_flowcal.py:243-254` | Re-read written CSVs and FCS; stop if counts differ; stop on metadata column collision. B02 checks in-memory rows only, so the disk re-read is our extension | `R/clean_gate.R:807` derives `n_removed` by diff, so in minus out is true by construction; no `sw_report_*` writer except unmix | S | rewrite (B02 no licence) | 4 | 3 |
| 9 | N11 | `sw_diff_ledger` | B7 diff | iidx `99_Aux.R:735-798`, `04a_TestsDA.R:296-323` | Sample x test included flag with reason; cluster-filter removals table; tile map | `R/differential.R:247-257` drops clusters silently; `~474` drops samples silently; attributes at `311-313` hold no support | S | rewrite (Apache-2.0) | 4 | 3 |
| 10 | N15 | `sw_annotate_manual` (stricter), `sw_annotation_table_write/read`, `sw_annotate_merge` | B7 annotate | CATALYST `R/mergeClusters.R:62-80`; MARMOT `inst/app/server-relabel.R:229-262` | Stop on missing, unknown or repeated keys; editable CSV; merge table. CATALYST does not stop on missing clusters; the table columns are our design | `R/annotate.R:436-443` warns, `:450` adds the column; `tests/testthat/test-annotate.R:253-262` locks the warning | S | rewrite (GPL >= 2; MIT) | 4 | 3 |
| 11 | N19 | `sw_plot_diff_heatmap` | B7 report | CATALYST `R/plotDiffHeatmap.R:115-278`, `R/utils-differential.R:256-265` | Top-n hits x sample heatmap, row z-score, source CSV; the yes/no significance bar is left out | Only volcano, boxplot and abundance plots exist (`R/differential.R:578,729,837`) | S | rewrite (GPL >= 2) | 4 | 3 |
| 12 | N21 | `.sw_run_log` (+ `step_ledger.csv`) | report, Part A | MARMOT `inst/pipeline/MARMOT_Pipeline.qmd:378-478,151,800,812,2081-2108` | `run_log.json` with status and last step, rethrowing `tryCatch`, per-step ledger. In MARMOT the crash fields are never written | `ROADMAP.md:140` (D6, planned); M45 dropped earlier; fold into D6 | S | rewrite (MIT) | 4 | 3 |
| 13 | N22 | `sw_cluster_split_table`, `sw_plot_cluster_tree` | B6 | MARMOT `inst/app/helpers/plot_helpers.R:604-695` | Full parent-child edge table across a k grid; tree figure. The source draws core edges only | `R/cluster_stability.R:55` works at one k; ggraph and tidygraph not in DESCRIPTION | S | rewrite (MIT; cite clustree) | 4 | 3 |
| 14 | N28 | `sw_event_labels_align` | B3, serves M62 | FlowSOM_protocol `R/FlowSOM_protocol.R:333-341,460-472,123-130` | Map labels between raw and cleaned events by (file, id); NA with reason for removed events; stop on duplicate or missing ids | `R/utils.R:585` sets `Original_ID`; no align function; `sw_cluster_predict` returns input rows only. Whether `sw_clean_run` keeps `Original_ID` is unverified | S | rewrite (no licence) | 4 | 3 |
| 15 | N30 | `sw_cluster_file_contribution` | B6 | FlowSOM_protocol `R/FlowSOM_protocol.R:352-370` | Observed over expected share per cluster x unit, with counts; thin cells `not_assessed` | `R/cluster_scores.R` `freq_z` already gives the robust z of the same shares; only the ratio is new | S | rewrite (no licence) | 4 | 3 |
| 16 | N37 | `sw_gate_edit_effect` | B3 | CytoExploreR `R/cyto_gate-helpers.R:828,1128-1132,1310-1313,1336-1359` | Before and after counts, shares, gained and lost events for an edited gate and its descendants. The source takes no snapshot | `R/gate_wrap.R:312-340` (`regate`); `R/clean_gate.R:778` counts per step only | S | rewrite (GPL-2) | 4 | 3 |
| 17 | N39 | `sw_gate_partition_check` | B3 | cytoverse `vignettes/3_Gating_1.Rmd:507-522,534-545` | Per parent: events in no declared child, in 2 or more, sum minus parent. The vignette example is a tautology; the check is our design | `R/gate_qc.R:254-256` masks are reusable; no sibling check anywhere | S | rewrite (MIT) | 4 | 3 |
| 18 | N44 | `sw_range_loss` | B2 | Pepper `Analysis Pipeline.Rmd:559,604-612,686-742` | Per file and channel: events below and above fixed bounds, plus outside-any counted once. Pepper counts an event twice if it is out on two channels | `R/ref_coverage.R` uses quantile bounds per batch; `R/clean_gate.R` margin step removes and reports one total | S | rewrite (MIT) | 4 | 3 |
| 19 | N03 | `tests/testthat/test-fail-loud-import-qc.R` | B1 | B02 `bin/clean.Rscript:45-50,63-68`, `bin/compensate.Rscript:25-31,37-56`; B06 `scripts_general_fns/g15_beads_functions.py:58,83` | Failure-injection tests: each stops and names file and step | Gaps: `R/qc.R:148-190` (`sw_qc_batch` no sample name), `R/utils.R:85-89`, `R/clean_gate.R:577` (duplicate basenames). Already covered: `tests/testthat/test-clean_gate.R:299,456,642,964`, `test-manifest.R:35`. No `tryCatch` in `R/` swallows errors | M | rewrite (B02 no licence) | 4 | 3 |
| 20 | N12 | `sw_diff_share_rank` | B7 diff | cyCONDOR `R/differential_testing.R:200-368,731-945` | Paired Wilcoxon and Friedman on shares with a hard pairing check; base stats. Kendall W and r computed by hand | `R/differential.R:194,260-265` handles two groups, no pairing; no `wilcox`/`friedman` in `R/` | M | rewrite (GPL-3) | 4 | 3 |
| 21 | N24 | `sw_plot_cluster_biaxial` | B6, app | iidx `06b_Plotting.R:1530-1839` | One cluster over the rest on shared bins for a marker pair. The unstained path in the source reads an undefined `agg_u` (`:1574,1583`) | `R/unmix_figures.R:842-929` covers unmix only; `R/cluster.R:299` is a heatmap; no cluster explorer app exists | M | rewrite (Apache-2.0) | 4 | 3 |
| 22 | N33 | `sw_gate_lda_check` | B3 | GME `CyTOF Linear Classifier/CyTOF_LDApredict.R:19-23,44-54`, `autogating.R:9-18` | LDA trained on gate labels from training samples; rejection threshold; unknown share and agreement on held-out samples. The agreement metric is our design | No LDA in `R/`; MASS is used (`R/gate_density.R:186`) but not in DESCRIPTION | M | rewrite (MIT); `MASS::lda` from Suggests | 4 | 3 |
| 23 | N41 | `sw_gate_operator_spread` | B3 | EFC `R/harmonisation.R:~321-418`, `R/comparison.R:51-120` | CV, relative bias and z across operators gating the same vials, with a population map. Donor variance share is dropped (M33) | No operator function; `plans/2026-11-11_reference-and-batch-diagnostics.md:204` puts variance components out of scope | M | rewrite (no licence) | 4 | 3 |
| 24 | N43 | `sw_gate_density_cut(method = "mixture")` | B3, S5 | A04 `analysis/functions/flow_functions.R:1-41`; EFC `R/naive_memory.R:135-162,172-247` | mclust mixture with BIC component choice; cut from the highest-mean component; NA with reason on failure, no fallback | `R/gate_density.R:30-76` has valley, percentile, mad only; mclust is in Suggests (`DESCRIPTION:54`) | M | rewrite (A04 Apache-2.0; EFC no licence) | 4 | 3 |
| 25 | N01 | `sw_manifest_bind_files` | B1 | B03 `source/run_analysis_pipeline.R:55-95`; MARMOT `MARMOT_Pipeline.qmd:640-664` | Exact-key join of manifest to FCS files; lists orphans both ways; stops on duplicate keys; reports `$NEXTDATA` | No `NEXTDATA` in `R/`; `R/unmix.R:1565-1600` stops only on zero files. `read.FCS` behaviour on `$NEXTDATA` not checked | S | rewrite (GPL-3) | 3 | 3 |
| 26 | N08 | `sw_fcs_channel_summary` | B1 | B06 `analyze_fcs.R:57-100`, `analyze_fcs_flowcal.py:197-218` | Per file x channel n, non-finite count, mean, quantiles in caller units | `R/anomaly_features.R:55-58` gives a feature matrix only; `sw_fcs_signal_audit` gives z only | S | rewrite (GPL-3) | 3 | 3 |
| 27 | N09 | `sw_plate_layout_read` | B1 | B06 `scripts_general_fns/12-read_layout_fns.R:9-96`, `analyze_fcs.R:75-85` | Plate grid to long well table; 1:1 well-to-file join; fixes the padded vs unpadded well mismatch | No well parser; unlocks M12 (`proposals.md:54`), which is not built | S | rewrite (GPL-3) | 3 | 3 |
| 28 | N17 | `sw_cluster_profile` | B6 | iidx `05_Profiling.R:115-217` | Cluster x marker quartiles with an All row; gate-label composition with an Unlabeled column | `R/differential.R:103-126` and `R/cluster.R:247-269` give medians only; composition half overlaps M70 | S | rewrite (Apache-2.0) | 3 | 3 |
| 29 | N18 | `sw_diff_signal_support` | B7 diff | iidx `99_Aux.R:840-923`, `06b_Plotting.R:24-60,167-185` | Join the per-sample above-floor flag (M68) to DE rows as `n_test`, `n_above_floor`, `rate`. Source's `match()` at `:170` has swapped arguments | `sw_diff_expression` returns no sample support; blocked until M68 exists (`proposals.md:77`) | S | rewrite (Apache-2.0) | 3 | 3 |
| 30 | N23 | `.sw_app_export_ui`, `.sw_app_export_server` | app | iidx `06_Dashboard.R:1025-1100,716-727`; MARMOT `inst/app/server-download.R:8-22` | Download source-data CSV and `app_state.csv` from `sw_qc_app` and `sw_unmix_app` | `grep -c download` gives 0 in `R/qc_app.R` and `R/unmix_app.R` | S | rewrite (Apache-2.0; MIT) | 3 | 3 |
| 31 | N38 | `sw_plot_gate_tree` | B3 | CytoExploreR `R/cyto_plot_gating_tree.R:64-123,126-150` | Static tree with median share and range across samples; node table. Helper and negated flags need a stated source | No tree function; numbers come from `sw_gate_qc` kept share | S | rewrite (GPL-2) | 3 | 3 |
| 32 | N45 | `sw_duplicate_tubes` | manifest or cleaning | ISC `immunosenescence.Rmd:~601-627,349-350` | Side-by-side table of files per tube plus a caller keep table with reasons. The source keeps most cells silently | Only duplicate `sample_id` stops (`R/manifest.R:116-118`); M11 not built | S | rewrite (no licence) | 3 | 3 |
| 33 | N47 | `gating-scenarios.qmd` S7b | B3 | BEV `R/process_FCM.R:37-46,84-104,~150-185` | Scenario where a gate fails and no substitute gate is drawn; no new function | `R/gate_density.R:236-325` already returns status `failed`; S7 is planned | S | rewrite (no licence) | 3 | 3 |
| 34 | N06 | `sw_qc_cleaner_overlap` | B2 | B02 `bin/clean.Rscript:12-68` | Run PeacoQC, flowAI, flowClean, flowCut on the same raw frame; removed share, kept-set overlap, removals over time | `plans/2026-11-13_scenarios-and-step-reports.md:131,146` approved flowAI and flowCut only; flowClean appears nowhere; M19, M51 cover one cleaner | M | wrap from Suggests (flowClean, flowCut Artistic-2.0; flowAI GPL >= 2; PeacoQC GPL >= 3); B02 idea only | 3 | 3 |
| 35 | N31 | `vignettes/failure-gallery.qmd` | scenarios (RD7), not B3 | FlowSOM_protocol `R/Figures_and_timing.Rmd:90-151,153-245,247-257` | Three planted breaks, each with the check that shows it. Needs N28, N30 and a new scatter-dominance diagnostic | `sw_synref_shift` (`R/synref_shift.R:9-10,208,270`) plants gain only; no existing check fires for the unscaled SSC-A case | M | rewrite (no licence) | 3 | 3 |
| 36 | N32 | `sw_clean_step_plot` | B2 | FlowSOM_protocol `R/FlowSOM_protocol.R:119-169,189-221` | Kept vs removed overlay per cleaning step and sample; caption counts from `clean$audit` | `R/clean_gate.R:800-815` audit holds counts; per-step event ids are not kept, so `sw_clean_run` must change | M | rewrite (no licence) | 3 | 3 |
| 37 | N20 | `sw_report_diff` | B7 report | iidx `06_Dashboard.R:195-330,983-989`, `06c_InfoText.R:100-155` | Seven fixed sections. iidx supplies only selectors, volcano, hit table, single-hit view and method text | `plans/2026-11-13_scenarios-and-step-reports.md:108-113` names it with no layout; needs N13 and the diff wrappers first | L | rewrite (Apache-2.0) | 3 | 3 |
| 38 | N05 | `sw_bead_peak_audit` | B1 | B06 `scripts_general_fns/g15_beads_functions.py:58,73-83,89-92,111-121` | Per-peak median, CV and fit against a supplied MEF table. The source returns a transform and reports none of these | No bead code in `R/`; M64 is absolute counts only; MEF tables cover few of ~64 spectral detectors | M | rewrite (GPL-3; FlowCal is Python) | 2 | 3 |
| 39 | N46 | `sw_cut_scope_check` | other | ISC `immunosenescence.Rmd:312-362` | Scope label and overlap warning for a quantile cut; pooled share equals `prob` by construction | `R/clean_gate.R:69,249,319` (`pool_samples`); `sw_gate_density_apply` `cuts` attribute | S | rewrite (no licence) | 3 | 2 |

## A3. Report-panel and app ideas by plan step

Every figure below is unmade. Under the plot-first rule, any conclusion drawn from one is `unsupported: no inspected figure` until it is saved and inspected, and making it needs a grill plan. Tables carry their own `table fallback` reason where they stand in for a figure.

**B1 import (`sw_report_import`)**
- N01: table of bound files, manifest rows without file, files without manifest row, multi-dataset files. `table fallback: names and counts, nothing to plot`.
- N02: batch x date tile map (fill = file count) with its source table; optional batch x serial map.
- N08: `channel_summary.csv`; per-file median and IQR dots, one panel per channel, files in acquisition or batch order.
- N05: peak table plus log median vs log MEF panels with a residual strip, one panel per channel.
- N09: plate map tile grid (wells without file marked).
- N03: none.

**B2 qc-clean (`sw_report_clean`)**
- N44: file x channel below and above counts; dot matrix via `sw_plot_batch_matrix()`.
- N32: up to five kept-vs-removed panels per sample (margins, qc, non_debris, singlets, live) with an `n_in / n_kept / n_removed` caption.
- N06: removed share and pairwise overlap tables; kept/removed ribbon per method along acquisition time.
- N07 and N21: `run_log.json` and `step_ledger.csv`; `table fallback: counts and skip reasons`.

**B3 gate (`sw_report_gate`)**
- N34: sample x node table (share of named ancestor, linear median, CV, NA reasons); per-donor points by node.
- N35: per-sample target-vs-ancestor overlay and edge-share table.
- N36: row x check tibble of template rows; optional tile figure.
- N37: before and after table; paired dots per node.
- N38: static gate tree beside a node table.
- N39: share-in-no-child and share-in-several-children dots per parent.
- N33: sample x population dot matrix (colour = unknown share, size = agreement), held-out samples only.
- N41: operator dots over donor mean; operator x population matrix of relative bias and z, no cutoff.
- N43: density with fitted components and the cut, plus a BIC panel.
- N46: per-sample dead-share dots with the pooled share as a reference line.

**B6 cluster (`sw_report_cluster`)**
- N30: cluster x unit dot matrix (colour = log2 ratio, size = n).
- N22: cluster tree across k from the full edge table.
- N24: cluster biaxial panel on shared bins, with an app tab once a cluster explorer exists.
- N17: marker x cluster quartile dots; cluster x gate-label composition bars with a grey Unlabeled segment.

**B7 annotate and diff (`sw_report_annotate`, `sw_report_diff`)**
- N13: donor ordination (colour = group, shape = batch, size = cells) with loadings and variance tables.
- N11: sample x test tile map plus `model_support.csv` and `cluster_removals.csv`.
- N12: paired-donor slope plot and a table of the smallest attainable p per population.
- N14: per-class recall and precision dots per held-out file; confusion heatmap; vote-share density.
- N15: merge table (original, label, n_clusters, n_events).
- N16: one hit card per named cluster; per-batch effect dot row.
- N18: volcano with fill = support rate; per-sample above-floor strip.
- N19: top-n hit heatmap with a source CSV.
- N20: the seven-section frame that hosts the above.

**App**
- N23: two downloads in `sw_qc_app` and `sw_unmix_app` (`*_source_data_<time>.csv`, `*_app_state.csv`). `sw_qc_app` builds its plot from `sw_diagnosis_top()` output, so its source data is that table.

## A4. Vignette scenario additions for `gating-scenarios.qmd`

Existing list: `plans/2026-11-13_scenarios-and-step-reports.md:75-81` (S1 to S7). All example values stay labelled illustrative. Each addition needs its function built first, except S7b.

| Scenario | Addition | Item | Functions shown |
|---|---|---|---|
| S1 | Run the template pre-flight before `sw_gate_run`; show a deliberately broken row (unresolved dim, parent listed later, missing `groupBy` column) | N36 | `sw_gate_template_check`, `sw_gate_run` |
| S1 | Gate tree and node table for the one-sample template | N38 | `sw_plot_gate_tree` |
| S1 | Check the declared child sets of a parent for events in no child or in two | N39 | `sw_gate_partition_check` |
| S1 | Count events outside fixed bounds before choosing `margin_range` | N44 | `sw_range_loss` |
| S2, S3 | Show the scope of a quantile viability or positivity cut and the share of the events it came from | N46 | `sw_cut_scope_check` |
| S4 | Counts and shares of a named ancestor per sample after moving gates; descendant population against the moved gate | N34, N35 | `sw_pop_stats`, `sw_gate_backgate` |
| S5 | Mixture cut for a marker with one broad positive shoulder, with the BIC panel and the failure reason | N43 | `sw_gate_density_cut(method = "mixture")` |
| S6 | Align manual labels on raw events to auto labels on cleaned events before scoring | N28 | `sw_event_labels_align`, `sw_gate_agreement` |
| S6 | Classifier cross-check of gate labels on held-out samples | N33 | `sw_gate_lda_check` |
| S6 | Several operators gating the same vials | N41 | `sw_gate_operator_spread` |
| S7 (sub-case S7b) | A sample where the density gate fails: status `failed`, no substitute gate, sample absent from the gate list. State that SpectraWeaveR has no flowClust wrapper, so this is not a like-for-like swap. Parent share per sample as points with a within-group robust z, failed samples open | N47 | `sw_gate_density_apply`, `sw_gate_qc`, `sw_gate_plot` |
| S7 | Compare counts before and after a hand edit of one gate and its descendants | N37 | `sw_gate_edit_effect` |

N31 (`failure-gallery.qmd`) is a separate vignette in the RD7 line of the plan, not a gating scenario. It depends on N28 and N30.

## A5. Impactful choices that need your approval

No value is proposed. Every parameter is a required argument with no default. Rows marked `inherited` are decided elsewhere and need a one-line confirm. The `min_events` row stands for one value per function, set separately.

| # | Parameter | Items | Status |
|---|---|---|---|
| 1 | Abundance cluster filter `min_prop` | N11 | inherited: `ROADMAP.md:254-255` (0.01, default at `R/differential.R:197`) |
| 2 | `min_cells` per sample x cluster | N11, N13, N30 | inherited: `ROADMAP.md:254-255` (3L, `R/differential.R:392`) |
| 3 | Minimum donors per group or per batch x group cell | N11, N12, N16 | new |
| 4 | What counts as a dropped sample (NA only, or NA and zero) | N11 | new |
| 5 | Share denominator and scale (sample events, parent, named ancestor, root) | N11, N12, N34, N38 | new |
| 6 | `min_events` floors for `not_assessed` | N17, N30, N33, N34, N35, N39, N08, N43, N46 | new, per function |
| 7 | Scale and cofactor of caller-transformed data (arcsinh 6000) | N13, N14, N17, N31, N32, N33, N34, N44, N46 | inherited: `CLAUDE.md` design constraints |
| 8 | Seed and event draw size (`n_events`, `max_events`) | N14, N24, N35, N32, N43 | new, per function |
| 9 | Quantile `probs` | N08, N17, N46 | new |
| 10 | Multiple-testing family and p-adjust method (cyCONDOR default is Bonferroni) | N12, N20 | new |
| 11 | Paired vs unpaired, groups compared, level order | N12, N16 | new |
| 12 | Rank test, post-hoc test and its gate threshold, exact vs normal, ties | N12 | new |
| 13 | Effect-size definition (r, Kendall W, log2 ratio vs difference) | N12, N16 | new |
| 14 | Hits to show (named list or top n by a column), `top_n` | N16, N19, N20 | new |
| 15 | Row scaling, z clip (CATALYST fixes 2.5) and abundance transform (CATALYST uses `asin(sqrt())`) | N19 | new |
| 16 | Sort key, column order, row and column clustering, colour limits | N19 | new |
| 17 | Summary over donor points and its interval method | N16 | new |
| 18 | Ordination input (CLR shares, marker medians, both), CLR pseudocount, log base, scaling, PCA vs MDS, dimensions | N13 | new |
| 19 | Marker medians: median vs mean | N13, N08 | new |
| 20 | kNN `k`, CV folds, `fold_by` (file or donor), events per class, class balance, tie rule | N14 | new |
| 21 | `min_vote` for a label | N14 | new |
| 22 | Training vs held-out files | N14, N33 | new |
| 23 | `label_source` (gate or cluster; cluster labels cap evidence) | N14, N17, N30 | new, same rule as `sw_batch_scores` |
| 24 | Cluster resolution: `k` grid, method per k, node vs metacluster | N22, N30 | inherited for SOM grid and seed: `ROADMAP.md` choice register; new for the k grid |
| 25 | Expected-share rule (proportional to events vs equal) and unit (file, batch, donor) | N30 | new |
| 26 | Edge draw threshold, core-edge rule and edge width statistic | N22 | new |
| 27 | Bins, `n_cluster`, `n_rest`, axis `limit_probs`, background type | N24 | new |
| 28 | Pooling of events (pooled vs file then summarise); gates defining labels; Unlabeled name | N17 | new |
| 29 | `unmapped` policy, `overwrite`, label vocabulary | N15 | new; changes an exported signature |
| 30 | Match key, case rule, recursion, duplicate and extra-file policy, `$NEXTDATA` policy | N01 | new |
| 31 | Header keywords for date, time, serial; date format; time zone; `run_gap_hours`; declared batch column | N02 | new |
| 32 | Policy for zero-event files, channel-set mismatch, and whether `sw_qc_batch` is changed to name the sample | N03 | new; code change |
| 33 | Edge width per channel and the target-ancestor pairs | N35 | new; edge distance rule inherited from `sw_gate_qc` |
| 34 | Group-to-channel map, name-match rule, missing-dim policy, helper gates shown | N36 | new |
| 35 | LDA rejection threshold and priors | N33 | new |
| 36 | Mixture: component candidates, BIC-gain rule, tail rule, model family, positive component, min events per component | N43 | new |
| 37 | Edit comparison: `group_by`, descendants (all or children), `count_floor` | N37 | new |
| 38 | Declared child sets per parent and whether overlap is allowed | N39 | new |
| 39 | Lower and upper bounds per channel, strict vs inclusive, missing-bound policy | N44 | new; relates to `margin_range` in `ROADMAP.md:217-225` |
| 40 | Which cleaners run, their settings, time bin width, kept vs removed overlap | N06 | flowAI and flowCut approved in RD4 and Q3; flowClean new |
| 41 | Operator reference, sigma, variance estimator, exclusion handling, population map | N41 | new |
| 42 | Keep table and the dead-share definition | N45 | new |
| 43 | Cut scope (pooled, per sample, control-derived), control tubes, tie rule at the cut | N46 | new |
| 44 | Which steps and channel pairs are drawn, samples shown | N32 | new; needs per-step ids in `sw_clean_run` |
| 45 | Floor definition and `min_above_floor`; `min_cells` must equal the DE call | N18 | floor set in M68; others new |
| 46 | Peaks, MEF table per detector, scatter floor, density fraction, bandwidth, fit model | N05 | new |
| 47 | Well regex, plate size, grid orientation, replicate rule, empty-well policy | N09 | new |
| 48 | Report frame: predictor, confounders, model family, effect threshold | N20 | new |
| 49 | Label alignment: direction, id and file columns, `removed_reason` | N28 | new |
| 50 | Report collect: join key and metadata columns | N07 | new |
| 51 | Skip rule for a step in the run log | N21 | must be user supplied; the log only records it |
| 52 | Export file name format and `shinyjs` dependency | N23 | new, small |
| 53 | Gate-tree node-size denominator, edge statistic, palette keyed by population | N38 | new |

## A6. Rejected or already covered

- **N04** (`sw_gate_density_fraction`) did not confirm. `sw_gate_density_2d()` (`R/gate_density.R:137-216`) already does 2-D kernel density gating with a contour. The narrow gap is per-event density rank and the discarded share per sample; extend `sw_gate_density_2d` instead. Value drops to 2.
- **N25** equals M37, **N27** and **N29** equal M62, **N42** equals M18. No new row.
- **Merged names:** N01 absorbs the duplicate `sw_manifest_match_files`. N13 absorbs `sw_sample_mds` and `sw_file_median_pca`. N11 merges the two iidx claims (diagnostic and report-app). N36 merges two survey entries.
- **Parts dropped inside kept items:** donor variance share in N41 (M33 covers it); the significance bar in N19; the rate-as-filter slider in N18; threshold lines in N24 (M68); any separation test or permanova in N13; mode in N08; numeric fill codes in N28; `JudgeAgreement` from N41.
- **Fold into existing plan rows:** N21 into ROADMAP D6 (with M45); N07 into the Part A contract; N47 into S7; N45 feeds M11; N38 overlaps the planned parent-share dots.
- **N03 is narrower than the survey said.** `sw_clean_run` and the manifest already have name-the-sample tests. MARMOT-style silent failures (B02, B06, B03) do not occur in `R/`, so the tests guard against regressions. Real gaps: `sw_qc_batch`, `sw_io_read_fcs`, duplicate basenames.
- **N18 is blocked** until M68 exists. **N19, N20 and N24** have no host function or app yet. **N09** only matters if M12 is built.
- **N31 is weaker than claimed.** No existing check fires for the unscaled SSC-A case, and N28 and N30 do not exist yet.
- **Wrong labels in the survey:** the "B7 annotate-diff" and "B3 gating" tags for N11, N12, N18, N31 and N45 do not match ROADMAP rows; I could not find "annotate-diff" in `plans/` or `ROADMAP.md`. Treat the plan-step column here as the scenarios-plan section.

## A7. Repos still unread or unreachable

- Unreachable: none reported by this pass. Two URLs failed: `RGLab/cytoverse-bioc-2023` (real repo is OzetteTech; the `seandavi` URL returned not found).
- Parts of read repos not opened: B02 `parse_samplesheet.py` and flowFun `R/utility_functions.R:856-877` (N01 claims, unverified); iidx `99_CheckInputs.R` (N01); iidx `06_Dashboard.R` beyond line 330, found by grep only (N20); cytoverse vignette `5_Gating_gatingTemplate.Rmd` (N36, N38 half unverified); `CyTOF_LDAtrain` (N33); `scater::calculateMDS` internals (N13); clustree package (N22); `gating_template_fitc.csv` swap in fibroblast (N36).
- Catalog repos from the earlier section 7 list that no confirmed item in this pass cites, so they remain unread as far as this evidence shows: C05 CytoNorm_Figures, C07 faust_flowcap4_analysis, C09 UNITO, D06 CytoGate, E10 R_flowcytometry_course, FlowCT (B10) differential functions, A05, A07, A11, A20, A24, A26.
- Not run: R was not on PATH in earlier passes, so AutoSpectral 1.8.8 exports for M73 are still unchecked. This pass executed nothing.
