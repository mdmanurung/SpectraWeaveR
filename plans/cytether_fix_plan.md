# Fix plan: HEB unpaired-batch coverage and CytofIn per-population arm (2026-11-10)

Context: stage B (ledger c020) left HEB cov90 at 0.70 (heb) and 0.80 (heb_pool_pop) on the shrink
stages and 0.68 / 0.82 on S17; dev stages pass (0.895). T6 EMD showed CytofIn meanshift worse than
no correction, explained by its one shift per marker and batch (diagnostic, no adapter bug).

## Approved rows

| ID | Decision | Source |
|---|---|---|
| F1 | Fix HEB unpaired-batch coverage first, then add the CytofIn arm | user 2026-11-10 |
| F2 | Add `cytofin_pop`: CytofIn meanshift within each FlowSOM population (the populations HEB uses), labelled adapted; the published global arm stays the primary competitor | user 2026-11-10 |
| F3 | HEB fix loop uses dev seeds 441-450 plus new dev seeds (count: open, F5); test seeds 431-440 and Z247 stay sealed | user 2026-11-10 |
| F4 | Diagnostic amendment to evaluate.py: coverage, RMSE and median z split by reference donors behind the effect (1 vs 2+); not scored | agent, non-scoring; benchmark card note |

| F5 | 10 extra dev seeds, 451-460 (20 dev seeds in total) | user 2026-11-10 |
| F6 | HEB fix: boundary-avoiding prior, penalised REML with a gamma(shape 2, rate -> 0) penalty on each variance-component SD (Chung et al., PMID 24092484; blme). F4 showed 1-donor effects with median |z| 1.6-99 (variance pinned at the floor); 2+ donor effects cover 0.885-0.904 | user 2026-11-10 |
| F8 | After F6 failed (ledger c022: prior lifts run variance, residual collapses instead), HEB residual = known sampling variance var(cells)/n per summary + one estimated extra variance (Fay-Herriot style); arms heb_pool_pop_kv and heb_pool_pop_kv_bap on dev seeds 441-460 | user 2026-11-10 |
| F9 | T6 adds `heb_kv_bap` (N4 model + F6 prior + F8 known variance) as an extra arm on EMD and OpenProblems; the approved `heb` arm stays | user 2026-11-10 |
| F10 | Next candidates, all approved: N1 heb_pool_all_kv_bap; N2 heb_afix_pop_kv_bap; N3 min_cells 20 and 50 as a labelled diagnostic (approved 10 stays); N4 T6 ablation heb_kv and heb_bap. Dev seeds 441-460 | user 2026-11-10 |
| F12 | Scale term (after T6 EMD decomposition: shape 51% of the gap to cycombine_ref; raw MAD ratios noisy, cross-fold r ~ 0.4): required argument scale_model = "none" or "log_sd". log_sd fits the lead's mixed model (donor fixed, study and run random, pooled per marker, bap, kv) to log SD per reference aliquot summary, known variance 1/(2(n-1)), no anchors; cells become c + exp(-delta_b) (x - shift - c), c = batch's reference median minus the shift. Arms heb_kv_bap_sc (T6) and heb_pool_pop_kv_bap_sc (synthetic). Pass: T6 pooled median-EMD gap to cycombine_ref shrinks by >= 0.06 and heb_kv_bap_sc <= heb_kv_bap in >= 3 of 4 folds; synthetic paired difference to the lead CI upper <= +0.01 on dev, shrink, S17, S18; S18 nd1 cov90 >= 0.85 | user 2026-11-10 (Q1-Q4) |
| F13 | Sub-cluster location (benchmark-side prototype, T6 only): heb_kv_bap_sc_cl = heb_kv_bap_sc on gated populations, then per population a k-unit SOM (seed 1, rlen 10) on stage-1 corrected calibrating cells, all cells mapped, k = min(4, floor(median calibrating cells per aliquot / 100)) >= 1; heb_kv_bap_sc refitted on raw values with population:subcluster labels; sub-clusters without an effect and unsplit populations keep stage 1. Pass: pooled median-EMD gap to cycombine_ref <= 0 and cl <= sc in >= 3 of 4 folds | user 2026-11-10 (Q1-Q4) |
| F11 | Fair cyCombine arm `cycombine_ref`: cyCombine package defaults (H1) fitted on calibrating reference cells only (per-batch scaling, SOM, ComBat), then its per batch x cluster x marker map clip(a + b * raw) applied to held-out cells. Matches batch_correct on 99.95% of cells when no cell is held out; held-out cells never move the fit. T6 EMD and OpenProblems | user 2026-11-10 ("fair cycombine") |
| F7 | Crash fix in the ported FlowSOM mapping metric: shared cell-type levels across splits; identical result whenever both splits hold every type | agent, no value change; deviation noted in code |

## Steps

1. Run the F4 diagnostic on the shrink and S17 dev sets (`scored_class`).
2. Ask F5 and F6 with the diagnostic in hand.
3. Run each HEB candidate on dev seeds; one ledger record per candidate.
4. Implement `cytofin_pop` in `adapter_cytofin.R`, add it to T6 (EMD 4 folds; OpenProblems 8 runs).
