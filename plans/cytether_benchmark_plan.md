# cyTether benchmark plan (draft for approval, 2026-11-09)

Worktree: `.claude/worktrees/cytether-benchmark`, branch `cytether-benchmark`, base `6845d1d`.
Source bundle: `plans/cyTether_implementation_benchmark_bundle_v1.zip` (MIT).

## Question

Is cyTether (`paired_shift`, and joint HEB if built) not inferior to cyCombine, CytoNorm and CytofIn on real data with repeated reference donors?

## Approved choices

| # | Choice | Source |
|---|---|---|
| Q1 | Endpoint: median EMD on a held-out reference donor, with safety gates, plus OpenProblems metrics | user |
| Q2 | 10 seeds per simulated stage S0 to S12 | user (departs from 100 in ROADMAP C6) |
| Q3 | Primary datasets: FR-FCM-Z247 and Zenodo 5911417. Secondary: Z52G, ZYAJ, Z2QV, OpenProblems sets | user |
| Q4 | arcsinh cofactor by assay: 5 for CyTOF, 6000 for spectral | user |
| Q5 | Populations: FlowSOM fit on the anchor batch, then frozen | user |
| Q6 | Non-inferior only if the EMD criterion and the OpenProblems criteria both pass | user |
| Q7 | Margin: within 10% relative to the best competitor's median | user (no source for 10%) |
| Q8 | Anchor: batch sharing most reference donors, tie to earliest. Leave-one-donor-out | user |
| Q9 | HEB: 10 seeds everywhere, one 100-seed run for HEB vs `paired_shift` only | user |
| Q10 | Competitors in documented reference mode, pinned versions | user |
| Q11 | One family of 6 comparisons (2 datasets x 3 tools), Holm-adjusted, donor-level bootstrap | user |
| Q12 | 10% relative per OpenProblems metric, direction-aware | user |

## Safety gates (pass or fail, not scored)

Zero change to earlier corrected values; no use of `truth_EVALUATION_ONLY.npz` in fit code; abstain when no bridging reference; protected-contrast error reported.

## Open items

- Aliquot noise floor is logged as a diagnostic. It does not gate.
- `cyCombine`, `CytoNorm` and a CytofIn package are not in the pixi env. Pin each commit before any run.
- CyTOF cofactor 5 is not yet checked against the dataset papers.
- Z247 donor count must be at least 3 for leave-one-donor-out. Check after download.
- Each conclusion needs an inspected saved figure (plot-first rule). Figure list comes after the data is downloaded.

## Steps

1. Pin and install competitors in the pixi env.
2. Download and checksum Z247 and Zenodo 5911417 into `$SW_DATA_DIR`.
3. Run `sw_design_audit()` on each. Stop if no donor-matched reference.
4. Export bundle S0 to S12 (C5) and compare `paired_shift` with the bundle's `reference_shift.npz`.
5. Run 10 seeds. Run the 100-seed HEB gate only after HEB exists (C6).
6. Real-data runs, then figures, then the report.
