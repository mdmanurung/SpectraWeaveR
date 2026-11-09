# cyTether + HEB benchmark strategy (draft, 2026-11-09)

Supersedes the "Steps" and "Open items" of `cytether_benchmark_plan.md`. Approved choices Q1 to Q12 there still stand unless a row below says `needs decision`.

Sources opened (main text only, supplements not opened):
- cyCombine: Pedersen 2022, PMC8971492.
- CytoNorm 2020: Van Gassen, PMC7078957.
- CytoNorm 2.0: Quintelier 2025, Cytometry A 107:69, PMID 39871681.
- CytofIn: Lo 2022, PMC8854441.
- OpenProblems `task_cyto_batch_integration` (MIT, tag v0.0.1-rc3), read from GitHub.

## What the four sources do and do not give us

| Source | Design we can reuse | Gap we must fill |
|---|---|---|
| cyCombine | EMD reduction and MAD score; Van Gassen, Trussart, DFCIb3 and iMUBAC sets; tool run only on datasets it was designed for | No split, no seeds, no tests. SOM is fit on corrected data. Best setting picked on the evaluated data. |
| CytoNorm 2020 | Train on donor A's controls, test on blinded donor B, both directions, 4 setups; EMD per gated population and marker; cluster-percentage CV check | One assay, 2 donors, no seeds, no named test. No MAD. |
| CytoNorm 2.0 | Aggregate-as-control mode; goal-quantile mode; EMD and MAD per gated population | Not a benchmark: no competitors, same files train and test |
| CytofIn | RMSD, R2 and subpopulation frequency against bead normalization; kBET and LISI; adjusted Rand index | Anchors not held out; stable-channel count tuned on test; bead data is a technical reference, not biological truth |
| OpenProblems | Replicate-split design; 7 metric components, 9 scores; three control methods that anchor scaling; MIT | Needs S3 data; Docker or Apptainer for the full pipeline; 9 scores carry about 33 to 44% missing results on the website |

None of them tests: frozen-model invariance, abstention, protected biological contrasts, interval coverage, or arrival order. The bundle's S0 to S12 do. So the strategy is: use their metrics so numbers compare with the papers, then add what they lack.

## Tracks

| Track | Data | Reference design | Question answered | Methods |
|---|---|---|---|---|
| T1 Synthetic | bundle S0 to S12, 10 seeds | held-out reference donor, truth in `truth_EVALUATION_ONLY.npz` (evaluation only) | recovery, abstention, invariance, contrasts | none, naive, paired_shift, HEB, CytoNorm, cyCombine, CytofIn adapters |
| T2 CytoNorm replica | FR-FCM-Z246/Z247, 10 plates, 2 donors, stimulated and not | train on donor A, test blinded donor B, then swap (4 setups) | does cyTether match or beat CytoNorm on its own design | paired_shift, HEB, CytoNorm, cyCombine, CytofIn |
| T3 CytofIn replica | Zenodo 5911417, 989 samples, 50 plates | unknown until the zip is inspected. The paper has one anchor aliquot per plate, so leaving a plate's anchor out leaves cyTether nothing to estimate from and it abstains. That tests abstention, not held-out agreement | agreement with bead normalization (RMSD). This endpoint differs from Q1, which needs a held-out reference donor | same |
| T4 cyCombine sets | Van Gassen, Trussart, DFCIb3 | technical replicates | EMD reduction and MAD, comparable to the cyCombine paper | same, only where design fits |
| T5 OpenProblems | human CyTOF, mouse spleen, Lille spectral | `is_control` samples | all 9 scores scaled on the control methods | cyTether wrapper writing `layers["integrated"]` |

Order of build: T1, then T2, then T5, then T3, then T4. Reason: T1 needs no download; T2 reproduces a published number; T5 reuses the same Van Gassen data.

## Metric set

| Metric | Definition (source) | Direction |
|---|---|---|
| Held-out EMD | pairwise EMD per gated population and marker, 0.1 bins, maximum over batches (CytoNorm 2020) | lower |
| EMD reduction | share of original EMD removed; pairs with both EMDs under 2 dropped (CytoNorm 2020, cyCombine) | higher |
| MAD score | median absolute change in MAD per cluster, marker, batch (cyCombine, CytoNorm 2.0) | lower |
| RMSD to bead | mean-expression RMSD to bead-normalized, 36 markers (CytofIn) | lower |
| OpenProblems 9 scores | as in the task repo, scaled between no_integration and perfect_integration | as listed |
| Safety gates | zero change to earlier output; no truth leakage; abstain without bridging reference; protected-contrast error | pass or fail |
| Interval coverage | 90% and 95% posterior intervals on synthetic truth (HEB only) | near nominal |

## HEB-specific tasks (what the papers never test)

1. Nested recovery: study, run and aliquot effects recovered against known truth on S3, S9, S11, checked against an `lme4` oracle (ROADMAP C6).
2. Ablation against `paired_shift` on S4, S9, S10, S11 with the 100-seed run (Q9).
3. Interval coverage on S4 and S9.
4. Cluster-validity check: CV of per-cluster cell percentage across references, as in CytoNorm 2020. HEB and `paired_shift` fit per population, so a split population is a failure to flag.
5. Abstention on S5 and S12.

## How gates and scores apply

Gates are requirements on cyTether and HEB only. Every competitor is scored on every metric, whether or not it can freeze. A competitor that fails a gate (cyCombine refits, so it fails frozen invariance) has that reported as a property, never as a disqualification. Otherwise "not inferior to cyCombine" would be true by exclusion.

## Fairness rules taken from the gaps

1. Settings are fixed before any result, or tuned on training donors only. cyCombine's paper recorded "the setting with the best overall performance", which tunes on the evaluated data. We do not copy that.
2. Each competitor runs in its documented reference mode. `cyCombine` has no stored correction model and refits, so it fails the frozen-invariance gate visibly. `CytoNorm` stores `fsom` and `clusterRes` and applies to new files. `Cytofin` saves anchor statistics, but each plate to normalize appears to need its own anchor file (function body not read).
3. Pin commits. CytoNorm `278f352` (GPL-2+, wrap only, never copy), Cytofin `cb3e24b` (MIT), cyCombine v0.2.5 `5b5903c` (MIT), CytofRUV `mtrussart/CytofRUV` `7088af8` (GPL-3, wrap only).
4. Each tool runs only where its design fits. Report "not applicable" rather than a bad score.
5. Resample reference donors and plates, never cells. Holm across the Q11 family.

## Decisions settled (user, 2026-11-09)

| # | Decision |
|---|---|
| D1 | Primary EMD: CytoNorm 2020 definition, absolute EMD, no "both under 2" rule. The reproduction gate keeps the paper's recipe. EMD reduction with the exclusion is secondary. |
| D2 | User downloads FR-FCM-Z246/Z247 in a browser and copies them here. No certificate bypass. md5 recorded on arrival. |
| D3 | Z247 has 2 donors: 2-fold swap, 4 setups (stimulated and not), as in CytoNorm 2020. |
| D4 | Reproduction gate: our CytoNorm, paper recipe, 10 seeds. Pass if 0.62 lies within our seed range widened by 0.01. Unverified tolerance. |
| D5 | OpenProblems metrics run natively on T2 and T3 data converted to its schema. The OpenProblems datasets (bucket closed) and cyCombine sets run only if data becomes available. |
| D6 | T2 uses the 10 manual gates shipped with Z247. FlowSOM-frozen elsewhere. |
| D7 | Claim = non-inferior to tool X on dataset Y, needing EMD and all 9 OpenProblems scores within margin. Holm over the 6 claims. T1 is its own family, by scenario. |
| D8 | Intervals on Z247 resample the 10 plates; the 4 setups are shown unpooled. Unverified method. |

## Seeds

Real data has no simulation seed, but FlowSOM, cyCombine (`seed = 473` default) and CytoNorm (`seed` argument) are stochastic. Fix one seed per tool and record it in the run manifest.

## Settled by the sources

Arcsinh cofactor 5 for CyTOF: cyCombine ("for all tools"), CytoNorm 2020 (Fig 2 caption) and CytofIn (gating) agree. This closes the earlier open item. Spectral stays 6000.

## Open facts

- Z246/Z247 file counts, populations and gates: not opened (site unreachable).
- Zenodo 5911417 zip contents unknown. It is 7.2 GB, md5 `0b1f9cd31d28a510e104846c0a0c94b5`.
- OpenProblems data sits on S3 only. `aws s3 ls --no-sign-request` on `resources/`, `resources_raw/` and `resources_test/` all return `AllAccessDisabled`. This route is closed unless the bucket reopens; T5 needs another data source.
- cyCombine dataset accessions (Trussart, DFCIb3, Krieg, Ogishi) are not in the main text.

## Findings after download and install (2026-11-09)

- Competitors installed in `/exports/para-lipg-hpc/mdmanurung/spectraweaver_testdata/Rlib_competitors` (outside the repo): cyCombine 0.2.5 (`5b5903c`), CytoNorm 2.0.12 (`278f352`), cytofin 0.0.0.9000 (`cb3e24b`). The Cytofin package installs as lowercase `cytofin`. Not yet loaded or run.
- Zenodo 5911417 downloaded, md5 `0b1f9cd31d28a510e104846c0a0c94b5` verified. Contents: 51 healthy control FCS files in 7 variants (Original, Bead, Beadlike, Meanshift, Meanshift_bulk, Variance, Z) for 6 cohorts, plus 500 `raw_DDPR_files` leukemia files, plus Figure_5 lymphoma files.
- The 989-sample bead-normalized leukemia ground truth is not in the zip: only raw patient files, bead versions exist for the 51 controls only. The paper's RMSD-to-bead on 989 samples cannot be reproduced from this record.
- Each cohort has one healthy control donor across its plates (cohort names differ: Healthy1, Healthy02, Healthy03, healthy_bm). No donor is measured across cohorts. So T3 has no held-out reference donor: it cannot test the Q1 endpoint. Any T3 agreement-with-bead score would fit and evaluate on the same control files.
- FlowRepository Z247: download form needs a CAPTCHA. User downloads it in a browser (D2). Page lists 40 FCS files (2 donors, 10 plates) plus `manual.RDS` and `manual.wsp`, with per-file md5.

- D9 (user): T3 is secondary. Agreement with bead on the 51 controls only, fit and evaluated on the same files, outside the Q1 endpoint and the Q11 claims.
