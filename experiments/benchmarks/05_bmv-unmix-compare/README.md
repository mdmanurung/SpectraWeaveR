# 05: four unmixing methods on the BMV 3-day-stim set

Status: **finished** (`runs/2026-11-09_bmv-unmix-compare-final`, 11 arms; reports `report_unpinned/report.html` (primary) and `report_all/report.html`). Tusoskill ledger: 22 evaluations in `runs/2026-11-09_bmv-unmix-compare/tusoskill`. Earlier runs `...-6arm` and `...-9arm` are superseded.

## Question and claim

Question (user, 2026-11-09): how much does the unmixing method change data quality on the same 5 raw Aurora files, so one algorithm can be chosen? No improvement is claimed. Differences are reported against the reference arm and the user decides from the figures. ROADMAP:186 says no cutoff for spread or condition number is sourced, so nothing is graded pass or fail.

## Data

Zip: `experiments/bmv_test/compare_3D/compare_3Dstim.zip`, unpacked to `experiments/bmv_test/compare_3D/data/` (git-ignored, checksums in `data/MD5SUMS`). Samples: A3, A6, B11, B2, D1 (40-colour panel, 3-day CD3/CD28 stimulation, Batch1).

| Arm | Folder under `data/` | What it is |
|---|---|---|
| `spectroflo` | `manual/` | SpectroFlo 3.3.0 export. Viability channel `LIVE DEAD NIR-A`. |
| `as102_basic` | `autoSpectral_v1.0.2/AutospectralUnmix/` | AutoSpectral 1.0.2, `weighted = F`, AF spectra, no variants. |
| `as102_var` | `autoSpectral_v1.0.2/AutospectralUnmix_optimized/` | 1.0.2 with per-cell spectral variants. |
| `as188_joint` | `autoSpectral_v1.8.8/AutoSpectral_unmixed/` | 1.8.8, `pipeline = "joint"`, variants, cleaned controls. **Reference.** |

Checked on 2026-11-09: all four arms have the same event count per sample, and `Time`, `FSC-A`, `FSC-H`, `SSC-A` are identical across arms. Events are paired by row. Columns are joined by detector name (`$PnN`); marker labels differ between arms. The viability channel is `Zombie NIR-A` in the AutoSpectral arms.

## Approved choices

| Choice | Decision | Source |
|---|---|---|
| Events | Live singlet non-debris cells, gated with openCyto | user, 2026-11-09 |
| Live gate | Gate per arm, keep events live in all arms (shared set) | user, 2026-11-09 |
| Reference arm | `as188_joint` | user, 2026-11-09 |
| Transform | arcsinh, cofactor 6000, all fluorescence channels | inherited: `CLAUDE.md`, confirmed by user |
| Separation markers | CD3, CD4, CD8, CD45RA, CCR7, CD127, CD25, CD56, TCRgd, CD161, TNF, IFNg, IL-2, IL-10 | user, 2026-11-09 |
| P1 nonDebris | `mindensity` on `FSC-A`, `gate_range = c(2e5, 1.5e6)`, per sample | plan 2026-11-09 |
| P2 singlets | `singletGate` on `FSC-A,FSC-H`, `wider_gate = TRUE` | plan 2026-11-09 |
| P3 live | `mindensity` on arcsinh viability, pop `-`, per sample and arm | plan 2026-11-09 |
| P4 negatives for channel j | below the per-arm per-sample `mindensity` cut, in all arms | plan 2026-11-09 |
| P5 spread | 1.4826 x MAD of negatives, linear scale | plan 2026-11-09 |
| P6 spread matrix | (rSD_j in i+j- minus rSD_j in i-j-) / sqrt(median_i+ minus median_i-) | plan 2026-11-09 |
| P7 lower tail | fraction of negatives below median - 3 rSD | plan 2026-11-09 |
| P8 leaning | Spearman rho among events negative for both; median abs and 10 worst pairs | plan 2026-11-09 |
| P9 stain index | (median pos - median neg) / (2 rSD neg) | plan 2026-11-09 |
| P10 clustering | `sw_cluster_run(method = "som", xdim = 10, ydim = 10, n_metaclusters = 20, seed = 42)` | plan 2026-11-09 |
| P11 agreement | ARI vs `as188_joint` per donor; cytokine+ % in CD4+ and CD8+ | plan 2026-11-09 |
| P12 summary | per-donor log2(arm / as188), donors as points, no cutoff | plan 2026-11-09 |
| P13 spectra | `sw_unmix_compare_spectra()`, `sw_unmix_diagnose()` | plan 2026-11-09 |
| Cut range | mindensity `gate_range = c(0, 3)` on arcsinh, all channels and arms (P4, P9, cytokine and lineage gates); default put cuts in the tails | user, 2026-11-09 |
| Cluster markers | lineage only: CD3, CD4, CD8, CD45RA, CCR7, CD127, CD25, CD56, TCRgd, CD161, CD57, CD38, CD7, CXCR5 | user, 2026-11-09 |
| P6 medians | median_i+ over all i+ events (above the arm's own cut), median_i- over shared negatives; rSD_j over i+j- and i-j- | user, 2026-11-09 (grill Q2) |
| P9 events | positives = above the cut in all arms (reference arm's cut where fewer than 100 shared positives); negatives = shared negative set | user, 2026-11-09 (grill Q3) |
| P13 1.8.8 spectra | report all three runs (`..._142441`, `..._143455`, `..._143803`); which one unmixed the data is not proven | user, 2026-11-09 (grill Q1) |
| Framing | pipeline comparison: arms differ in spectra, control cleaning and algorithm. Controlled arms (same spectra, both algorithms) only if results are close | user, 2026-11-09 (grill Q4) |
| Decision rule | readout priority fixed before results: (1) negative spread, (2) stain index, (3) artefacts, (4) downstream agreement. Higher ranks decide; lower ranks break ties. No cutoffs | user, 2026-11-09 (grill Q5) |
| Decision rule (revised) | guardrails reject first (pile-up index, TNF-CD8 streak, hypernegatives, reproducibility); then rank by G1 held-out control error (cell controls + H12; beads secondary); G2 simulation confirms; quantile spread then quantile stain index break ties. Replaces the 2026-11-09 order because the 1.0.2 zero pile-up fools MAD-based spread | user, 2026-11-10 (Fable plan) |
| Shared negatives | below the reference arm's (`as188_joint`) cut only, so adding arms never changes earlier numbers | user, 2026-11-10 |
| G1 held-out controls | random 50/50 event split per control, seed 42, `Reference Group Clean` files; spectra re-extracted from one half per code line, off-target error scored on the other | user, 2026-11-10 |
| G2 simulation | `sim.flow.data()`, 1e5 cells per draw, seeds 42/43/44, dim layer at 0.005 relative abundance, 1.8.8 spectra, AF and variants with real scatter; both code lines; RMSE and bias per dye | user, 2026-11-10 |
| Dropped as ranking metrics | MAD spread, centring, roundness, negative-cell fraction (pile-up), exclusivity 18 and 21 (gate placement), TNF/IL-2 positivity, ARI (inside noise floor) | user, 2026-11-10 |
| Spread and SSM (30) | Nguyen robust SD (84th - 50th percentile) primary, 2.5-97.5 spread as check; quadrature SSM; all-pair tilt | user, 2026-11-10 |
| Minimum population (30, 31) | metric is NA when a population has < 100 events (`MIN_N = 100`) | user, 2026-11-10 ("do all above") |
| Pinned cuts (30) | excluded for spread, stain index and SSM, both source and target; kept for tilt and the TNF-CD8 streak | user, 2026-11-10 ("do all above") |
| CD4/CD8 gate (19a) | openCyto `quadGate.tmix`, K = 3, quantile 0.8, pooled over the 5 donors per arm | user, 2026-11-10 |
| Cytokine gates (19a) | TNF: mindensity `gate_range = c(1, 4)`; IFNg: cytoUtils tailgate `tol = 0.01, num_peaks = 2, ref_peak = 1`; no unstimulated control exists | user, 2026-11-10 |
| Cytokine % positive (38, 40) | one fixed cut for every arm and donor, TNF (BV510-A) and IFNg (BV570-A) above arcsinh 1; replaces the automatic gates for % positive. 1.5 was used first (`figures_fixed_cut/`) and is kept as a check | user, 2026-11-10 |
| Variant subsample seed (39) | seed 42 for `var102_n26` | user, 2026-11-10 |
| Noise-model GLS arms (42) | `gls188_noise`: AutoSpectral 1.8.8 `estimate.noise.model()` + GLS unmix, package defaults; `gls188_noise_spill`: same with `include.spillover = TRUE`; noise-fit seed `asp$bird.seed` | user, 2026-11-10 |
| Other cytokines (47) | IL-2, IL-10, Th2: % positive at a fixed cut of arcsinh 1 in every arm and donor, plus cut-free rank metrics (density, Spearman rho, shared-negative rSD); linear channels transformed arcsinh(x / 6000) in memory | user, 2026-11-10 ("Both") |
| Ablation baseline (48+) | start `as188_joint`; chain target `s102_c188`; final target `as102_var`; `opt188` extra endpoint | user, 2026-11-10 (D1) |
| Ablation criterion | fraction of the per-donor gap `as188_joint` → `s102_c188` closed on CD8T TNF squeezed share and share above arcsinh 1; "resembles 1.0.2" = median ≥ 0.8 and every donor ≥ 0.5 on both; "contributes" = ≥ 0.1 in ≥ 4/5 donors | user, 2026-11-10 (D3) |
| Ablation readouts | full gating for every arm (19a gates, TNF+/IFNg+ at cut 1, TNF+ among IFNg+, 39b squeezed/pos share, Th2/IL-10/TNF/IFNg negative rSD) | user, 2026-11-10 (D4) |
| Ablation order | Block A object swaps (spectra × AF × variants cube + 1.0.2-thresholds-only) first; then regenerate chain output side first: variants → filter → AF → spectra → bead gate → cell negative → cell gating → speed | user, 2026-11-10 (D8/D12) |
| Bead pre-gate fix | FSC 2e5–9e5, SSC 1.3e6–2.25e6; also a single arm `bead_wide` off `as188_joint` | user, 2026-11-10 (O1) |
| Ablation seeds | `asp$bird.seed` inside AutoSpectral; `set.seed(42)` around extraction; one replicate of chain1 at seed 43 | user, 2026-11-10 (D5) |
| Zombie NIR control | reuse the 84,528-event `Reference Group Clean` file in regenerate arms (shipped raw has 1,578 events) | user, 2026-11-10 (D13) |
| Cell universal negative | two chain arms: `chain6a` = same `Unstained (Cells)`, cell rectangle FSC ≥ 4e5 (SSC 1e5–1.6e6, FSC ≤ 4e6); `chain6b` = H12 unstained sample; `opt188` uses whichever the chain result supports (to confirm) | user, 2026-11-10 (O3, after plots in `runs/2026-11-10_bmv-cell-gating/`) |
| Cell landmark gating | `define.gate.landmarks(n.cells = 2000, percentile = 70)` for the six metabolic cell controls (`gate.name` lymphocyte; control-file copy in `runs/2026-11-10_bmv-cell-gating/`); density gating for beads | user, 2026-11-10 (O2, from the tune.gate grid) |
| `opt188` settings | AF `refine = TRUE`; `clean.controls` `af.remove = FALSE`; `n.passes = 1`; `get.spectral.variants(stained.sample = A3 raw)`, necessity threshold 0.01 | user, 2026-11-10 |
| Unmix speed | `speed = "slow"` in every AutoSpectral unmix call; under `pipeline = "joint"` it has no effect (`repro188` fast = `s188_c188_slow` on A3), so earlier joint arms stand | user, 2026-11-10 |
| Legacy pipeline in ablation | `chain8_legacy` (chain7 + `pipeline = "legacy"`) and Block A `legacy188`, `legacy102in`; 1.8.8 legacy is not 1.0.2 code | user, 2026-11-10 |
| Passes grid (52) | `n.passes` × `n.af.passes` over {1, 5, 10} on the zip 1.8.8 objects, joint, slow; no winner chosen until a selection criterion is set | user, 2026-11-10 |
| Control-file gate fix (49, chain4 on) | copy `experiments/runs/2026-11-10_bmv-ablation-chain/fcs_control_file_gatefix.csv`: `Unstained (Cells).fcs` moved from `smallGate_density_orphan` to `smallGate_2`; `Negative_1 (Beads).fcs` alone in `smallGate_density_orphan`; original untouched; chain7 on uses the landmark copy. Needed because 1.8.8 `define.flow.control` rejects mixed cells/beads in one gate ("Inconsistent values for control.type") | user, 2026-11-10 |
| Chain7 split (49) | `chain7a_pregate`: `define.gate.landmarks(n.cells = 2000, percentile = 70)` on the landmark control-file copy, polygon applied to the lymphocyte-group files (cell controls and Unstained), then the chain route on the gated copies; `chain7b_legacy`: the same gate passed to `define.flow.control(gating.system = "landmarks", gate.list)` for the spectra only; `chain8_legacy` builds on chain7a. Reason: `get.spectra.automated()` and `get.spectral.variants()` take no gate list | user, 2026-11-10 |
| `opt188` gate routes (53) | two arms, because `get.spectra.automated()` takes no gate list: `opt188_pregate` = landmark polygon applied to cell controls and cell negative, gated copies, then `get.spectra.automated(10000, 500, 0.95)` and `get.spectral.variants()` unchanged; `opt188_legacyspec` = `define.flow.control(gate.list = list(lymphocyte = gate))` (density kept for beads and viability, per the bead decision) → `clean.controls(af.remove = FALSE)` → `get.fluorophore.spectra()`, variants from files without the landmark gate; both take the negative option and `n.passes` from the command line | user, 2026-11-10 |
| AF405 proof chain (60–65) | overlap plot; BV510-component knock-out and knock-in; dose α = 0, 0.25, 0.5, 1; spike-in BV510 at arcsinh 1, 2, 3 into raw CD4T/CD8T cells of all 5 donors (cells chosen with seed 42); per-cell anatomy of 10 cells (seed 42). Pre-registered: KO ≤ 0.15, KI ≥ 0.25, dose monotonic, spike-in recovery ≥ 0.9 under 1.8.8 | user, 2026-11-10 |
| Negatives for spread (34) | cells negative in every one of the 11 figure arms, each by its own cut; replaces "reference arm negatives" for spread, stain index and guardrails, because the reference cut trims only the reference's tail (fable review, 2026-11-10) | user, 2026-11-10 |
| G2 arm names | in the simulation, "as102_var" = 1.0.2 code + 1.8.8 spectra + variants (the s188_c102 recipe) and "as102_basic" = 1.0.2 code + 1.8.8 spectra; not the real-data arms of the same name | verified from args.txt, 2026-11-10 |
| Cleanup for plots | PeacoQC RemoveMargins, fluorescence min -21,761 / max 1e7, scatter at instrument range, `AF Index` excluded; shared across arms; plots only, not ranking | user, 2026-11-09/10 |
| `opt188` pre-gate scope (53, 49 chain7a) | the landmark polygon gates the stained cell controls only; the cell negative keeps its rectangle gate, because the polygon left it 578 of 6,374 events, below AutoSpectral's 1,000-event floor. Supersedes "cell controls and cell negative" in the gate-routes row | user, 2026-11-10 |
| Cell negative alternatives | `chain6b_neg` and opt188 `neg_h12` dropped: the H12 unstained file holds the same events as `Unstained (Cells).fcs` (identical `exprs`); only `neg_same_fsc4e5` is used | user, 2026-11-10 |
| `opt188` arms (53) | 6 arms: per route (pregate, legacyspec): joint n.passes/n.af.passes 1/1, joint 5/5, and legacy. Legacy takes no passes (AutoSpectralRcpp sends them only to the joint call), so legacy arm names carry none | user, 2026-11-10 (grill Q1-Q3, Q8) |
| `opt188` extraction (53) | extract spectra, AF and variants once per route; every arm of that route unmixes from the saved objects, so arm contrasts change only pipeline and passes | user, 2026-11-10 (grill Q9) |
| `opt188` choice | by design, not by outcome: opt188 = pregate route, joint, 1/1 (documented settings); the other 5 arms are sensitivity arms reported descriptively. Reason: no TNF readout separates signal from artefact without an unstimulated/FMO control; TNF+ among IFNg+ tracks the TNF+ level and the TNF-IFNg odds ratio reverses the ranking | user, 2026-11-10 (grill Q7; supersedes Q4/Q6 criterion) |
| Passes grid verdict (52) | descriptive: 1/1 kept by design; the grid shows how readouts shift with passes (its 1/1 arm is `as188_joint` itself) | user, 2026-11-10 (grill Q10) |
| `opt188` report readouts | CD8T, donor median with each donor as a point: TNF squeezed share, TNF+ above arcsinh 1, TNF+ among IFNg+, TNF+ among IFNg-, TNF-IFNg odds ratio, negative rSD ratio vs `as188_joint` for TNF, IFNg, Th2 and IL-10 | user, 2026-11-10 (grill Q12) |
| Coherence caveat (66) | new figure from existing scorer CSVs: TNF+ among IFNg+ and IFNg-, and the odds ratio, per arm; INDEX section 05 rewritten to match | user, 2026-11-10 (grill Q11) |
| Test 1, all-marker screen (75) | every fluorescence channel; parent CD3+ T cells for T-cell markers, live cells otherwise; arms `as188_joint`, `as102_var`, `s102_c188`, `s188_c102`, `swap_var102`, `swap_af102`; readouts squeezed share near 0 and negative rSD ratio vs `as188_joint` (39b method), on `as188_joint`'s gates | user, 2026-11-10 ("run all in parallel", accepting the proposed settings) |
| Test 2, per-dye × channel map (76) | the 32 existing per-dye arms; every channel scored as in test 1, on `as188_joint`'s gates | user, 2026-11-10 (same) |
| Test 3, variant code × controls (77) | AF405 only; 1.0.2 `get.spectral.variants` on the clean (zip rectangle) controls with the zip 1.0.2 settings (`som.dim` 10, `n.cells` 2000), seed 42; swapped into the 1.8.8 setup like `var102_only_Alexa_Fluor_405`; joint, speed slow; scored with 51 | user, 2026-11-10 (same; SBUV510/eFluor 450 extension left open) |
| Test 4, SBUV510 bead-gate arm (78) | SBUV510 reference and variants rebuilt from the widened bead gate (chain5 objects), swapped alone into the 1.8.8 setup; joint, slow; scored with 51 | user, 2026-11-10 (same) |
| Test 6, own SSM (82) | spillover spread matrix computed from the 1.8.8 single-stain controls with a stated direction (row = spreading dye, column = receiving detector/dye), compared with the zip raw SSM to fix its orientation | user, 2026-11-10 ("run in parallel", accepting must-do 1–4) |
| Test 7, AF origin (83) | 1.0.2 variant code on the AF405 control cells split into low-AF and high-AF halves by their unmixed AF (median split), seeds 42 and 43, otherwise test 3 settings; plus a figure of the variant spread shape vs the matched unstained AF spectra | user, 2026-11-10 (same) |
| Test 8, SBUV510 decomposition (84) | SBUV510 reference-only and variants-only swaps (1.0.2 zip objects) into the 1.8.8 setup, joint, slow, scored with 51 | user, 2026-11-10 (same) |
| Test 9, replicates (85) | test 3 rerun with seed 43; test 5 rerun with the cell AF step checked (log which AF step ran on cell controls; if none ran, report why) | user, 2026-11-10 (same) |
| Downstream arms (grill Q1) | carry both `as188_joint` and `opt188` downstream | user, 2026-11-10 (grill Q1) |
| Affected-marker flag (grill Q2) | a marker is flagged in a population if the donor-median squeezed share rises ≥0.10 vs `as188_joint` AND the negative rSD ratio ≤0.8, with the same sign in all 5 donors | user, 2026-11-10 (grill Q2) |
| Affected-marker flag, revised | the rSD criterion is dropped: flag if the donor-median squeezed share rises ≥0.10 vs `as188_joint`, same sign in all 5 donors. Reason: a pile at 0 beside a wider tail keeps rSD >0.8 (TNF in `s102_c188` +0.37 and the AF405-only swap +0.21 were missed). Supersedes the row above | user, 2026-11-10 |
| `af405_afremove` arm (script 87b) | `clean.controls(af.remove = TRUE)` with the chain7a landmark gates (`define.flow.control`, `gating.system = "landmarks"`); cleaned AF405 events written to FCS, then `get.spectra.automated` with the zip arguments and 1.8.8 `get.spectral.variants`; joint, slow; scored with 51. Inputs: controls `Reference Group Clean`; landmark gate recomputed on them (n.cells 2000, percentile 70); `clean.controls(downsample = FALSE)` with `positive.n` raised above the AF405 positive count (1e6), because `get.universal.negative` otherwise keeps only 500 positives; AF405 FCS holds the cleaned positive events only | user, 2026-11-10 |
| Th2 follow-up (grill Q3, script 86) | score Th2 (eFluor 450) in the existing AF405 proof arms; no new unmixing | user, 2026-11-10 (grill Q3) |
| Forced AF removal on AF405 (grill Q4, Q9; script 87) | two arms: `clean.controls(af.remove = TRUE)` and `get.af.spectra(refine = TRUE)`, each rebuilding AF405 with 1.8.8 code; all other arguments as in the zip; joint, slow; scored with 51 | user, 2026-11-10 (grill Q4, Q9) |
| Internal-negative TNF check (grill Q5) | not run | user, 2026-11-10 (grill Q5) |
| Spike-in (grill Q6, script 88) | Th2, IL-10 and TNF; +1 arcsinh (as the TNF spike); arms `as188_joint`, `as102_var`, `s188_c102`; 5 donors; figures per marker | user, 2026-11-10 (grill Q6) |
| Downstream comparison (grill Q8, Q11; script 89) | gating, Q2 flags and clustering (P10 settings, lineage markers) for `as188_joint` and `opt188`; descriptive, the user chooses | user, 2026-11-10 (grill Q8, Q11) |
| Commit (grill Q7, Q10) | experiments only (benchmarks, slurm wrappers, README), on `main`; package edits reviewed separately | user, 2026-11-10 (grill Q7, Q10) |
| Test 5, automated route on raw controls (79) | `get.spectra.automated` with the zip arguments on the raw 1.0.2 `Reference Group` controls; AF405 and SBUV510 reference offsets measured with the script 74 metric; no unmixing | user, 2026-11-10 (same) |
| UMAP kNN test (scripts 92, 92b) | UMAP only (`sw_dimred_run`, `n_neighbors` 15, `min_dist` 0.2, inherited defaults); arms `as188_joint`, `opt188_pregate_joint_p1af1`, `af405_afrefine`, `as102_var`, SpectroFlo (its rows pair exactly with ours); CD3+ T cells from `as188_joint`'s `CD3pos` rows, 10,000 per donor, seed 42; two marker sets (14 lineage; lineage + 17 functional: every other panel marker, viability and PE left out); z-scored and unscaled; per-cell kNN overlap at k = 30 vs `as188_joint` and vs SpectroFlo, median per donor; noise floor `as188_joint` re-embedded with seeds 1–3; descriptive, no threshold | user, 2026-11-10 (UMAP grill Q1–Q11) |
| Unmixing setting for this dataset (grill 2026-11-10, Q1) | `as188_joint` (1.8.8 zip), **final**. Evidence: single-dye spikes of 12 markers in CD4T/CD8T (88, 88b) and live CD3− cells (88e), TNF + GLUT1 co-spike (88d). Caveats: Th2 in CD8T (13.7% of cells lose >50% of a spike; Th2 analysed in CD4T), TNF in CD3− cells (6.8%), CS reads 3–4% high. Rejected: `af405_afrefine` (GLUT1 loss 15.3% in CD3− cells), `opt188` (Th2 CD8T 22.9%; CD38 shift merges clusters MC06/MC09) | user, 2026-11-11 (confirmed after co-spike and CD3− spike) |
| Th2 in CD8T (grill Q3) | reported with a caveat (11.8–16.1% of CD8T cells lose >50% of a Th2 spike under 1.8.8 zip); Th2 analysed in CD4T | user, 2026-11-11 |
| Co-spike (grill Q6, Q7; script 88d) | TNF +1 arcsinh with GLUT1 at +0, +1, +2 in the same CD4T and CD8T cells (seed 42, 1.8.8 dye spectra); setups `as188_joint`, `af405_afrefine`, `as102_var` (positive control); plus the existing TNF spike split by natural GLUT1 | user, 2026-11-11 |
| CD3− spike (grill Q8, Q9; script 88e) | the 12 spiked markers, +1 arcsinh, in live CD3− cells (`as188_joint` gates, seed 42); setups `as188_joint`, `af405_afrefine` | user, 2026-11-11 |
| Raw companions (grill Q4) | ggcyto hexbins for every summary figure: flagged markers vs GLUT1 per arm; lineage markers of the shifted clusters; spike before vs after per cell | user, 2026-11-11 |
| Spike passes (grill 2026-11-11; script 88h) | n.passes / n.af.passes 5/1, 1/5, 5/5 on `as188_joint` (joint, speed slow); 1/1 reference reuses the 88/88e c188 unmixes; Th2, IL-10 and TNF (control); CD4T, CD8T and CD3−; inherited levels 0–3, seed 42, lost = recovery < 0.5 at +1; descriptive, final setting stays 1/1 unless the user decides after | user, 2026-11-11 |

## Still to be set

- Nothing open. The 1.0.2 AF row may inflate its `max_sif` and condition number; the report notes it.

## Steps

1. `01_manifest-gate.R`: manifest, identity check, openCyto gates, shared event set, gate plots.
2. `02_metrics.R`: spread, matrix, artefacts, correlation, stain index.
3. `03_downstream.R`: clustering, ARI, cytokine frequencies.
4. `04_spectra.R`: spectra angles and `max_sif`.
5. `report.qmd`: figures and decision table.

Scripts use exported package functions only. Run with `experiments/slurm/run_script.sbatch <script.R> <run_dir>`.

## Outputs

`runs/2026-11-09_bmv-unmix-compare/`: `gate_counts.csv`, `shared_events.rds`, `quality.csv`, `negative_spread.csv`, `spread_matrix.csv`, `artefacts.csv`, `neg_correlation.csv`, `stain_index.csv`, `clusters_ari.csv`, `cytokine_freq.csv`, `spectra_compare.csv`, `figures/`, `report/provenance.csv`, `sessionInfo.txt`.
