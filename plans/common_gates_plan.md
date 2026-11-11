# Common SpectroFlo gates, propagated by cell (user, 2026-11-11)

Request (user, 2026-11-11): gate CD4T and CD8T, and the non-debris, singlet and live gates, on SpectroFlo, then propagate the labels by cell to every other analysis so the arms are comparable.

Fact: rows pair across arms. `93_flag_raw_scatter.R:115` stops unless FSC-A is identical per row in every arm, and SpectroFlo pairs the same way (INDEX §16).

## Decisions (grill, 2026-11-11)

| ID | Decision | Source |
|---|---|---|
| G1 | Gate SpectroFlo once with the 19a template: debris, singlets and live by mindensity; CD3+ by mindensity; CD4T and CD8T by quadGate.tmix K = 3, pooled over donors. Same parameters as `19a_gate_arm.R`. | user, Q1 |
| G2 | Scope: everything. §01, §03–§07, §09–§13, §15–§25: every analysis that uses a gate is rerun on the SpectroFlo labels. | user, Q2 (departs: pick was headline set) |
| G3 | Replace the per-arm-gate results. The old numbers stay in their run folders; INDEX and the dashboard show common-gate numbers only. | user, Q3 (departs: pick was keep as sensitivity) |
| G4 | The deck stays paused until the re-gate is done. Artifact shell: https://claude.ai/artifact/HkUGBqTSpb6HYee4HMuaZE | user, Q4 |
| G5 | Spikes rebuilt: spiked cells are re-picked from the SpectroFlo labels and every spiked set is re-unmixed on SLURM (speed "slow"). Inherited: +1 arcsinh, seed 42, spike fraction 0.25, loss cut recovery < 0.5, co-spike levels +0/+1/+2. | user, Q5 (departs: pick was re-score) |
| G6 | Checkpoint: stop after gating SpectroFlo. The user approves the gate figures per donor before propagation. | user, Q6 |
| G7 | Cytokine cut stays at arcsinh 1, cofactor 6000, in every arm. | inherited: `38_fixed_cut.R`; user, Q7 |
| G8 | Propagation is by row. Any arm whose rows do not pair with SpectroFlo (FSC-A and SSC-A per row) stops with an error; it is never matched approximately. | default |
| G6a | Checkpoint passed: the SpectroFlo gates are approved as is (user, 2026-11-11). The vertices come from `runs/2026-11-10_bmv-gated-pooled/gs/spectroflo`, in arcsinh. CD8T is CD8 > 0.758 & CD4 < 1.27. CD4T is CD8 < 0.758 & CD4 > 0.288. DN and DP take the rest. In A3, CD4T is 72.1%, CD8T 17.9%, DN 8.3% and DP 1.6% of CD3+. The labels are in `runs/2026-11-11_bmv-common-gates/labels.csv.gz`, and all 55 arm × donor pairs pair. | user, checkpoint |
| G1a | Every gate is pooled over all 5 donors (collapseDataForGating TRUE, groupBy 5). This covers nonDebris, singlets, live and CD3pos as well as the tmix. It supersedes G1/G6a, where only the tmix was pooled. The SpectroFlo gates are refitted in `runs/2026-11-11_bmv-common-gates-pooled/` (`97b_spectroflo_labels_pooled.R`), and the checkpoint repeats. | user, 2026-11-11: "make pooled gates on all samples", "do not gate per donor" |
| G7a | The cytokine cut is arcsinh 1.5 (cofactor 6000) for every cytokine, every arm and every section. This replaces G7's cut at 1. | user, 2026-11-11 (G12 answer, scope confirmed) |
| G10 | Re-clean from scratch. Fit the pooled SpectroFlo gates on the uncleaned files to get live. Drop margin- or negative-flagged events within that live set. The flag rules are inherited from `12_cleanup.R`: PeacoQC margins at −21,761 and 1e7 plus the scatter limits, or below −2 arcsinh in at least 3 fluorescence channels, in any of the 11 P arms. Write clean3 files for all 73 arms. Refit the pooled gates on clean3, with a checkpoint, then propagate those labels. | user, 2026-11-11 |
| G10a | The pooled template and methods are approved for the clean3 refit. The clean2 pooled run (`runs/2026-11-11_bmv-common-gates-pooled/`) is the reference for the cut shift. | user, 2026-11-11 |
| G11 | Figures for other arms colour each cell by its SpectroFlo label. Gate lines appear only on SpectroFlo panels. | user, 2026-11-11 |
| G10b | Clean3 pooled gates are approved as the common labels (user, 2026-11-11), and propagation and reruns are launched. The labels are in `runs/2026-11-11_bmv-common-gates-clean3/labels.csv.gz`. Cuts: live 1.841, CD3 0.928, CD4 0.293 and 1.277, CD8 0.761. There are 70 unique arms and 350 files, and all of them pair. | user, checkpoint |
| G13 | The negatives that define the negative rSD, and so the squeeze cut (0.5 × rSD), are cells below arcsinh 1.5. This applies in every batch and in the spikes. Batch 2 made this choice; the user confirmed it afterwards. | user, 2026-11-11 |
| G14 | §10's pre-registered knock-out prediction (CD8T squeezed ≤ 0.15) is reported as failed: 0.160 at cut 1.5. The baseline shift from 0.12 to 0.17 is stated beside it. The threshold is not changed after the fact. | user, 2026-11-11 |
| G15 | §03 uses 34c (`runs/2026-11-11_bmv-common-cut15-34c`) as the primary run, with per-arm, per-channel mindensity cuts refitted on the common gates. G13 applies to cytokine negatives, not to lineage channels that have their own fitted cuts. 34d, with 1.5 on all 32 channels, is kept as a note only. | user, 2026-11-11 |
| G16 | In §03 (34c), the cytokine channels BV510 and BV570 define negatives below arcsinh 1.5 (G13). Lineage channels keep their mindensity cuts (G15). Only the streak and tilt rows are rerun. | user, 2026-11-11 |
| G17 | The raw cluster companions (94) show clusters that shift the same way in all 5 donors and hold at least 1% of live cells in the 1.8.8 zip reference. The rule is stated in the caption. | user, 2026-11-11 |
| G9 | New scripts and new dated run folders only; no script a running SLURM job uses is edited. | standing rule |

## Steps

1. Gate SpectroFlo and export per-row labels (`sample, row, debris_ok, singlet, live, CD3pos, CD4T, CD8T, CD3neg_live`). Draw the gate figures per donor. **Checkpoint G6.**
2. In parallel with step 1: inventory every script and section that reads a gate, and name the change each needs.
3. After approval: propagate the labels to every arm's GatingSet or event tables, then rerun the inventory in parallel batches.
4. Rebuild the spiked raw sets from the SpectroFlo labels and re-unmix them on SLURM, then re-score.
5. Redraw the figures, rewrite INDEX and the dashboard, then build the deck.

## Inventory (step 2, 2026-11-11)

Every script gets membership with `gh_pop_get_indices()` from one of three gated dirs. All three use the clean2 rows (`cleaned/2026-11-09_bmv-unmix-compare-clean2`):

- **P** `runs/2026-11-10_bmv-gated-pooled/gs/<arm>` holds 11 arms with the pooled tmix.
- **A** `runs/2026-11-10_bmv-ablation-gated/gs/<arm>` holds 62 arms.
- **G** `runs/2026-11-10_bmv-gated/gs/as188_joint` is the older non-pooled gating. The 36, 37b, 39b, 43, 45b, 62 and 88* lineage uses it, not the pooled gates §02 describes.

Labels must carry both `row` (clean2) and `raw_row`, because the spike scripts and 34 work on raw rows. Step 3 builds `97_common_gs.R`: new GatingSets for every arm in P and A with the SpectroFlo labels as gates, keeping the node names. Each script copy then changes only its gated-dir knob.

Batches after 97 has run, each writing one new run dir per script:

1. **Fixed cut and cytokines:** 38, 40, 46, 47, 34, plus the 95 gate and viewer restyles.
2. **51 hub:** rescore every 51 call, then 66, 63, 67, 52b, 87 and 87b.
3. **Row-based scorers:** 36, 37b, 39b, 43, 45b, 69, 70, 75, 76 and 86, then 86a and 86b.
4. **Scatter and downstream:** 72 (4 sets), 89, 89b, 92 then 92b, and 94. 93 runs after 86b from batch 3.
5. **Spikes on SLURM:** the 88 build, the spiked-set unmix, then the score. Then 88b, 88e and 88d, then 88h, then 88c, 88f and 88g. Also 62, the spiked unmix, then 64. The `spk_<setup>_base` unmixes are reused.

Batches 1, 2, 3 and 5 are independent of each other. §04, §06, §08 and §14 have no gate step and stay as they are.

## Opus audit (2026-11-11): sound with fixes

The audit corrects the inventory in these places:

- **G equals P.** G has the same template as P and the same live, CD3pos, CD4T and CD8T counts per donor. It is not an older non-pooled gating.
- **Inputs the inventory missed.** 01 builds `shared_events.rds` from per-arm gates, and 02 fits `mindensity_cuts.csv` on those cells. Both feed 34, so both rerun before 34. 10, 19 and ten more 95 restyles were also missing; they go in a final batch 6.
- **Inline paths.** 75, 76, 86, 87 and 87b write their gate paths inline, so a knob does not reach them; each copy edits those lines.
- **Old run dirs.** 62, 63 and 64 write into the old `2026-11-10_bmv-af405-proof`. Every copy writes to a new dir with a must-not-exist guard.
- **The 51 race.** 51 writes the shared gated dir, and parallel calls race on it. It also gates any missing arm with its own gates without saying so. The 51 copy stops when an arm is missing and never calls 19a.
- **Stale readers.** 66, 67, 86a, 86b, 93, 88c, 88d, 88f, 88g, 88h and the 95 scripts read pinned old runs. Each copy repoints them to the new runs.
- **Batch order.** The batches are not independent. They run in this order: batch 0 (97_common_gs, full pairing check, then 01 and 02 rebuilt), then batches 1, 2, 3 and 5 in parallel, then batch 4, then batch 6.
- **Pairing check.** It must compare Time, FSC-A, FSC-H and SSC-A on all rows for every donor. So far only FSC-A on live rows was compared.
- **Consequences of the new gates.** The spike level assignment is reshuffled (already accepted under G5). The 92 UMAP sample and the 89 SOM sample change. The 51, 39b, 45b and 87 squeeze cuts are refitted on the new populations.

The audit also raised three open choices for the user. They are logged below as G10 to G12.
