# Plan: next round of technical diagnostics (2026-11-10)

Branch: `worktree-cytoscan-diagnostics`. Status: **awaiting approval**. No code written for this plan.
Goal (user, 2026-11-09): "work autonomously until we can consistently diagnose for any technical
problems including but not limited to batch effects extensively". "Consistently" has no numeric bar
yet; D1 and D2 below propose one. Every number here is quoted from a file, with its path.

## Where things stand

Simulation finds most faults. 15 of 17 fault classes have a check with mean AUROC 0.98 or higher
(`reference/figures/fault-detection-matrix.png`, cited in `reference/diagnostics.qmd`). The weak class is
a split into two modes in one marker: 0.67 in a file, 0.86 for a batch with seeds from 0.57 to 1.0
(same figure). Limits in simulation (`reference/figures/detection-limits.png`): a shift is found at
0.5 SD (`median_z` 0.99) and is near chance at 0.2 SD; a stretch at x1.2 (1.00) but not x1.05 (0.57);
a rare-population shift at 10% of events (0.95) but not 2% (0.49).

Real files carry the same limits at larger sizes. On 36 Aurora donor files the hit marker's own z
finds a 1 SD shift at 0.95 and a x1.5 stretch at 0.94; the max over 40 markers drops these to 0.74 and
0.69 (`reference/figures/real-injected-marker-level.png`), which is why Q14 keeps per-marker rows and
no per-file score. Negative-population checks score 0.93 to 1.00 on real files, 0.84 at a 0.3 SD shift,
on the 10 of 40 channels whose mixture separation was 3 or more (`reference/figures/negative-faults-real.png`).
A proportional leak is found at 15% (0.94) and not well at 5% (0.68) (`reference/figures/spill-leak-real.png`).
Batch faults on 144 real files in 4 batches: `pair_emd` 0.92 at 0.5 SD and 0.75 at 0.25 SD,
`spread_ratio` 0.87 at x1.3, `order_rho` 0.86 for drift (`reference/figures/batch-faults-real.png`).
k-means at k = 6 is stable (seed-pair ARI 1.00); k = 12 moves a median 12% of events
(`reference/figures/cluster-stability-real.png`). Shuffled labels give 0.50 to 0.54 in simulation
(`reference/diagnostics.qmd`, line 130).

Relativisation (E1) reaches the no-instrument floor in simulation only: at donor CD3 log SD 0.3 the batch EMD goes from 0.164
observed to 0.034 per batch, truth 0.035 (`experiments/runs/2026-11-09_relativisation-v2/fig_relativisation.png`).
It cannot be tested on the Aurora extract (OMIQ-transformed, no CD3 or single-stain controls).

Coverage table (`reference/diagnostics.qmd`, lines 72 to 97): 14 of 24 rows say "not tested" on real
files. `experiments/integration/aurora_exvivo/10_file_faults_real.R` (running at the time of writing,
seeds 1:5, 25% faulty, batch 1) injects 9 of them: saturation, dead channel, non-finite, clog, drift,
spillover, debris, bimodal, rare population. After it, 5 rows stay untested: header mismatch, all
markers widened, composition change, batch split into two modes, relativisation gain.

Dataset facts to reconcile. `reference/diagnostics.qmd` line 214 says 160 files and ledger iteration 60 "160+";
`experiments/runs/2026-11-09_aurora-exvivo-real/files.csv` has 149 rows and
`plans/2026-11-09_faust-wrapper.md` says 149 files (37, 37, 37, 38 per batch). Donor is nested in batch.
Only the four `BMV_REF` files (one per batch) give a batch contrast without donor; on 11 of 40 markers
their EMD exceeds the same-batch donor median, BV510-A 0.307 vs 0.057
(`experiments/runs/2026-11-09_aurora-exvivo-real/fig_ref_vs_donor_emd.png`). All three entries in
`experiments/datasets/registry.csv` are single-batch control sets, so none tests batch or E1.

Ledger: `experiments/runs/2026-11-09_tusoskill/experiments/results.jsonl` has 73 lines; `reports/status.json`
says iteration 73 of 100 on the task "file-level anomaly features and flags". `final_method/final_report.md`
still says "Iterations: 50 / 50" with TODO sections, so the ledger has outgrown its benchmark definition.

Rule gaps found: `experiments/benchmarks/07_fault-detection/run.R` and `08_detection-limits/run.R` `source()`
files from `R/`, against rule 2 in `experiments/README.md`; script 10 uses `pkgload::load_all()`.
The faust plan names its script `05_faust_before_after.R`, but `05_marker_level.R` exists; next free is 11.

## Work items, ranked

| # | Item | Purpose | Evidence produced | Hours | Depends on |
|---|---|---|---|---|---|
| W1 | Fill the coverage table from script 10 | Close 9 "not tested" rows | `fig_file_faults_real.png` inspected, copied to `reference/figures/` with its source CSV; table rows updated; ledger line 74 | 2 | job finishing |
| W2 | Scorecard against the bar | Make "consistently" measurable | `experiments/benchmarks/11_scorecard/run.R` reads every `reference/figures/*-source.csv`, applies D1/D2, writes one table and figure: per fault class, met / not met / untested; fills `final_report.md`; fixes the `source()` rule in 07 and 08 | 5 | W1, D1, D2, D9 |
| W3 | Cover the 4 cheap residual rows | Header, all widened, composition, batch split on real files | Extend scripts 10 and 09 with these faults; add rows to the table | 3 | W1 |
| W4 | Bimodal split check | Lift the weakest class (0.67 / 0.86) | New per-file, per-marker check in `sw_fcs_signal_audit()`; benchmark 07 rerun and script 10 rerun; AUROC before vs after | 6 | D3 |
| W5 | Separation gate for negative checks | Say on which channels the negative checks mean anything | `min_separation` required argument; rows below it get `NA` and a reason; count of covered channels per file in the table; app shows it | 4 | D4 |
| W6 | Events-per-file sweep | Test whether the small-fault limits are sampling-bound | Benchmark 08 with the D6 grid; `detection-limits-by-events.png`; rare-population share on the same grid | 3 | D6 |
| W7 | Batch scores by role | A batch contrast free of donor on the real files | `sw_batch_scores()` gains an optional `role` split; reference-only EMD per marker (4 files, no AUROC); figure beside `fig_ref_vs_donor_emd.png` | 4 | D7 |
| W8 | Second dataset, then faust | Test E1, the unmixing-control checks and batch without donor; FAUST before/after | Rerun scripts 02, 06, 09, 10 and the E1 fit on the new data; `11_faust_before_after.R` per the faust plan | 10 + 8 | D5, D8, D10 |

Per-marker reading aid: `sw_qc_app()` already shows the largest values per check and all values per subject
(`R/qc_app.R`). A file x marker tile view per check (display only, no cutoff) goes into W5 as one tab.

Parallel: W1 runs alone until the job ends. W4, W5, W6 and W7 are independent and can run together once
their decision rows are answered. W2 waits for W1 and D1, D2. W3 follows W1. W8 waits for data and installs.

## Decision rows (each needs user approval)

| # | Choice | Recommended | Reason | Status |
|---|---|---|---|---|
| D1 | The consistency bar: metric, tier, seeds, CI, null | Per fault class: the named check's mean AUROC on real files with the fault injected is 0.90 or more over at least 5 seeds, with the lower 95% CI at 0.80 or more; shuffled-label null at 0.55 or less | 0.90 sits above the weak results (simulation: split 0.67 and 0.86; real: batch shift 0.75 at 0.25 SD, leak 0.68 at 5%). At 0.90 two real batch rows fail today: stretch `spread_ratio` 0.87 at x1.3 and drift `order_rho` 0.86 at 0.5 SD (`batch-faults-real.png`). Fork: keep 0.90 and carry those two as open rows (retest at a larger size or after W7), or set 0.85, which passes them. The CI condition is unchecked yet; `batch-faults-real-source.csv` shows `ci95_half` 0.16 at 0.25 SD, so W2 is where it is first evaluated | needs user approval |
| D2 | Fault size the bar applies to, per class | The smallest size already tested on real files: shift 1 SD, stretch x1.5, negative shift 1 SD, leak 15%, batch shift 0.5 SD, batch stretch x1.3, drift 0.5 SD; rare population 10% | A bar without a size means nothing; smaller sizes (0.3 SD, 5% leak, 0.25 SD batch) are below 0.90 now and would make the goal unreachable without new methods | needs user approval |
| D3 | Bimodal check method | Own-code bimodality coefficient from the published formula first; `diptest` (GPL) wrapped from Suggests if the coefficient fails benchmark 07 | No new dependency for the first try; GPL code is wrapped, never copied | needs user approval |
| D4 | `min_separation` for negative checks | Required argument; the experiments used 3 (`reference/diagnostics.qmd`, line 189), which is the experiment's value, not validated | A rows-below-cutoff rule changes which channels report, so it is the user's | needs user approval |
| D5 | Second dataset | Ask the user for a lab Aurora run first; else search FlowRepository and Zenodo. Required: multi-batch, same material in every batch, per-batch single stains, linear export, batch not nested in donor. No accession is proposed | None of the three registry entries qualifies; E1 and the control checks have no real-file evidence | needs user approval |
| D6 | Events-per-file grid for W6 | 500, 2000, 5000 events; rare-population shares 2, 5, 10% | 500 was the simulation setting and 5000 the real-file extract; the ledger (iteration 7) says 500-event sampling noise dominated a 0.1 SD jitter | needs user approval |
| D7 | `role` split in `sw_batch_scores()` | Optional `role` vector; scores reported per role; no pooling across roles | Reference files are the only donor-free contrast; 4 files give EMD, not AUROC | needs user approval |
| D8 | Install faust and scamp (GitHub RGLab, GPL-3) into the pixi env | Yes, into `.pixi` of the main checkout; Suggests + Remotes as in `plans/2026-11-09_faust-wrapper.md` | Route "wrap now" is approved; the real-data check is blocked on the install; F1 to F19 in that plan still need answers | needs user approval |
| D9 | `k` for scorecard reruns of the cluster checks | 6 | ARI 1.00 at k = 6; k = 12 (used in `02_diagnose.R`) moves 12% of events (`cluster-stability-real.png`) | needs user approval |
| D10 | Ledger: close at 73 or extend | Extend to 100: log W1, W3 to W7 as iterations 74 onward, then fill `final_report.md` with the scorecard | The 50-iteration floor is met; the report's benchmark definition must be updated either way | needs user approval |
| D11 | Residual size rows below D2 (0.3 SD, 5% leak, 0.25 SD batch, 2% population) | Keep them in "What the checks cannot see" with their measured AUROC; do not chase them this round | Simulation shows effect size, not k or cohort contamination, sets these limits (`detection-limits.png`, ledger iteration 63) | needs user approval |

Impactful choices stay required arguments with no default (`CLAUDE.md`). Every new check (W4, W5, W7)
gets its own ledger line, benchmark row and coverage-table row before it is called covered.

## Goal met under the proposed bar

The goal is met when every row of the coverage table has a real-file AUROC at or above D1 at its D2
size, or sits in "What the checks cannot see" with its measured limit; no row says "not tested" except
those that need data the user declined to supply (relativisation gain, if D5 is declined); the
shuffled-label null stays at chance; the scorecard reruns from the `*-source.csv` files and prints the
met / not met / untested counts; `sw_qc_app()` lets a reader browse per-marker rows per check. Batch
effects count as covered when shift, stretch, split and drift each meet the bar on real multi-batch
files and the reference-only view (W7) is in the report.

## Order of work

1. Decision rows D1 to D11 (you).
2. W1 as soon as the job ends; W3 right after.
3. W4, W5, W6, W7 in parallel.
4. W2 scorecard; ledger and `final_report.md`; docs and CLAUDE.md rows; tests; atomic commits.
5. W8 when data and installs arrive.

## Decisions answered (user, 2026-11-10)

| Row | Answer |
|---|---|
| D1 consistency bar | Named check's mean real-file AUROC >= 0.90 over >= 5 seeds, lower 95% CI >= 0.80, shuffled-label null <= 0.55. |
| D2 size per class | The smallest size already tested on real files; smaller sizes are listed as limits. |
| D5 second dataset | Search public repositories; download only after the user approves a specific dataset. |
| D8 faust | Install faust and scamp (GPL-3) now; then ask the FAUST choice rows F1 to F19. |

## Decisions answered, round "fix all gaps" (user, 2026-11-10)

| Row | Answer |
|---|---|
| G1 spillover-spread check | Per channel pair i -> j: SD of j-negative events where i is above a required quantile, over SD of j-negative events where i is below its median. Log ratio, robust z across files. One required knob (`spread_quantile`). |
| G2 same-donor reference | New required argument `reference_group` (NULL = today). Robust z taken within each group. Re-benchmark the four donor-limited gaps on CytoBatchFlagR. |
| G3 FAUST rows F1-F19 | FAUST conventions: all shared fluorescent markers; `root`; depth 0.01; quantile 0.5; experimental unit and imputation `""`; empirical bounds; name occurrence 0; supervised `NA`; seed 1; FAUST density-subsample defaults; annotations approved; 5000 events, subsample seed 1; 8 threads; arcsinh 6000; after = `sw_correct_run()` as tested; exact phenotype match; ARI per file; median per phenotype per file. |
| G4 ledger | Keep target 150. |
| G5 batch scores with groups | Per donor across batches, combined by the median over donors. |
| G6 emd_z with groups | EMD to the pooled events of the file's own group; z within the group. |
| G7 FAUST correction | `sw_correct_run()` with `norm_method = "rank"`; other cyCombine defaults (8x8 SOM, rlen 10, seed 42). |
| G8 FAUST markers | 37 markers: drop AF1-A, AF-A and LIVE DEAD NIR-A. |
| G9 reference_group details | Confirmed: global_sd_z cohort SD within group; agg_size per group; batch scores skip groups with < 2 batches, order_rho within group, median over groups per row; benchmark on asinh(x/6000), 30 shared markers. |
| G10 next for the four limits | Two-way residual: robust z of the value after removing donor and batch medians (median polish); test on the CBF donor-by-batch design. |
| G11 FAUST annotationForestDepth | 3 (FAUST default), passed explicitly. |
| G12 FAUST prop_by_batch | Range over batches of the per-batch median share, before and after correction. |
| G13 FAUST correction markers | Batch correction uses the same 37 markers as FAUST (AF and dead dye excluded). |
| G14 FAUST phenotype match (F17) | `match_by = "marker_state"`: sort each phenotype's marker tokens before matching. FAUST orders markers per run, so exact strings matched only the unclassified bin. |
| G15 unstained gate (detector-noise test) | Box from the same donor x instrument gated panel: FSC-A, SSC-A and FSC-H/FSC-A within its 5th-95th percentiles. Revised from 1st-99th: the authors' Fortessa gate holds a 2.4% SSC tail, so the 99th-percentile box kept granulocytes. |
| G16 detector-noise reading | Report the unstained / panel-negative Fortessa-over-Aurora SD ratio per donor (median, range) with no cutoff; word the claim by how close it is to 1. Inherited and confirmed: 6 channels, authors' maps + asinh(x/1000), 1.4826 x MAD, panel negatives from the 08 run, all 9 donors. |
| G17 AF700 stain test | Target AF700 only; single stains gated with the G15 box (5th-95th, FSC-A, SSC-A, FSC-H/FSC-A); CD3 split as in 02 (k-means k = 2 on asinh(x/6000), seed 1, nstart 10, < 200 events per side gives NA); reading as G16. |
| G18 CD8-dim test | Split CD4 (APC) per panel file with k-means, k = 2 on asinh(x/1000) of relativised values, seed 1, nstart 10; subsets CD4+ T (primary true CD8-negative) and CD4- T; < 200 events gives NA; reading as G16. Revised: the panel files are already CD3+ (lympho-single_cd3), so a CD3 split only halves the T cells. |
| G19 AF700 compensation test | (E1) MASS::rlm (Huber) of AF700 on APC + AA750 in panel CD4+ T cells (k-means CD4 split, G18), linear relativised values, plus AF700 SD on asinh(x/1000) before and after removing the fit; (E2) AF700 median CD3+ minus CD3- in the APC and AA750 single stains (G15 gate, G17 split, authors map), linear scale. Reading as G16. |

## Decisions answered, round "plan all gaps" (user, 2026-11-10)

D3 and D4 above are answered by G20 and G24. D4's text reference is `reference/diagnostics.qmd` near line 229 (not 189).

| Row | Answer |
|---|---|
| G20 bimodal statistic (D3, A1) | Own-code bimodality coefficient (BC); wrap `diptest` from Suggests only if BC misses the bar in 07. |
| G21 bimodal input and events (A2, A3) | Values as given (as `median_z`/`sd_z`); all events in the file. |
| G22 bimodal score unit (A4) | Inherited per benchmark: 07 max \|z\| per file; 11 the hit marker's z. Report the new check's own AUROC too. |
| G23 bimodal real-file run (A7) | Rerun 11 at 0.5, 1 and 1.5 SD, 5 seeds, with a `bimodality_z` arm. |
| G24 negative gate cutoff (D4, B1) | `min_separation` required, no default; benchmarks use 3 (inherited from 06). |
| G25 negative gate unit (B2) | Per channel: median separation across files, not per file. |
| G26 negative gate action (B3) | Gated rows keep the row with value `NA` and the reason in `detail`; `separation` always kept; new per-channel check `split_ok`. |
| G27 negative gate scope (B4) | `neg_*`, `hyper_negative_share`, their z, and pair checks where either channel is gated. |
| G28 CBF batch score (C1) | New check `rest_emd_z`: EMD of each batch to all other batches pooled, robust z over batches within each marker. |
| G29 CBF AUROC unit (C2) | Inherited from 04: hit marker among markers. |
| G30 CBF seeds and sizes (C3, C4) | Time 1 seed, then 20 seeds; sizes 0.5 SD and 1 SD. |
| G31 CBF stop rule (C7) | If `rest_emd_z` stays below the bar at 0.5 SD, record a limit and close the gap. |
| G32 frozen reference scope (D0, D1) | CBF counts as the real multi-batch dataset (lifts the ROADMAP freeze for this test); reference B01-B09, scored B10-B18. |
| G33 frozen reference content (D2, D3) | Per-marker median and MAD of file-level medians and SDs, S3 list with digest (as `sw_frozen_model`); `median_z` and `sd_z` only, via a new `reference` argument on `sw_fcs_signal_audit`. |
| G34 frozen reference fault (D4) | 1 SD shift on one marker in 80% of scored files; 10 seeds. Runs after G20-G23. |
| G35 G20 fallback triggered | BC missed the bar in 07 (0.58). User: run the `diptest` fallback anyway ("Do all"). Hartigan dip per file and channel, robust z over files, as benchmark arms in 07 and 11 (`dip_z`); moved into the package only if it meets the bar. |
| G36 frozen reference calibration | Approved as planned (2026-11-11): reference batch01-09, scored batch10-18 (G32); asinh 6000, 25 named markers, all kept events (G9); median_z and sd_z, frozen vs cohort (G33); clean null with per-batch \|z\| and a Spearman trend on batch date, reference files scored in-sample as baseline; share of \|z\| >= 2, 3, 4, 5 with no single cutoff; G34 fault rerun saving per-file z; stop rule: if clean later batches often pass \|z\| >= 3, record a limit and name rolling or bridge references as a later fix. |
