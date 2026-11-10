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
| F7 | Crash fix in the ported FlowSOM mapping metric: shared cell-type levels across splits; identical result whenever both splits hold every type | agent, no value change; deviation noted in code |

## Steps

1. Run the F4 diagnostic on the shrink and S17 dev sets (`scored_class`).
2. Ask F5 and F6 with the diagnostic in hand.
3. Run each HEB candidate on dev seeds; one ledger record per candidate.
4. Implement `cytofin_pop` in `adapter_cytofin.R`, add it to T6 (EMD 4 folds; OpenProblems 8 runs).
