# Spec: CD3 (reference-population) relativisation

Branch: `worktree-cytoscan-diagnostics`. Date: 2026-11-09.
Status: core functions and extensions E1 to E4 built from this spec.

## 1. Source

Glehr et al. 2026, EBioMedicine 131:106427, doi 10.1016/j.ebiom.2026.106427, open access (CC BY),
read as PMC13488037. Linked from https://www.colibri-cytometry.com/references. The formula is the
paper's Box 1. The authors' code (`ggrlab/cytobench`, `ggrlab/otcyto`) was not opened and its
licence is unknown, so our code is written from the paper only.

## 2. The method in the paper

For each sample s and each fluorescence channel x there is a separate control: the same specimen
stained with anti-CD3 (UCHT1) conjugated to the fluorochrome of channel x, acquired on the same
cytometer with the same settings as the full stain. The control is compensated and gated on
lymphocytes. Its events are split into CD3- and CD3+ by k-means with k = 2. Then

    m[s,x] = median of channel x over the CD3- events of the control
    M[s,x] = median of channel x over the CD3+ events of the control
    X~     = (X - m[s,x]) / (M[s,x] - m[s,x]) * 10^4

for every event of the full stain of sample s. It runs on compensated, linear intensities, before
any transform. It removes an additive offset (m) and a multiplicative gain (M - m) per channel.

## 3. Assumptions (from the paper's Box 1 and Discussion)

- A1. The reference marker's median on reference-positive cells is the same across donors.
- A2. All values sit in the cytometer's linear range.
- A3. The same clone in every fluorochrome; saturating antibody amounts.
- A4. The control is acquired like the sample (same cytometer, settings, time).
- A5. The reference population is present in every sample.

Our code cannot check A1 to A4 from the data. It reports the numbers that show when they break:
event counts, the gain, and the gain's spread across units, batches and time.

## 4. Formal workflow (implemented)

Notation: a *unit* is what one reference set applies to (a sample in the paper; a batch when the
controls are per batch). A *control map* is a table with one row per (unit, channel): which
control file calibrates that channel for that unit.

| Step | Function | Input | Output | Caller's choices (required, no default) |
|---|---|---|---|---|
| 0 | caller | raw FCS | compensated or unmixed, linear, gated controls and samples | gating, unmixing |
| 1 | `sw_relativise_split_kmeans()` | controls, control map | logical "reference-positive" per event per control | `split_transform`, `seed`, `nstart` |
| 1' | caller | controls | same list from a manual gate | the gate |
| 2 | `sw_relativise_fit()` | controls, control map, positive list | model: `m`, `M`, `gain`, event counts per unit and channel | `min_events` |
| 3 | `sw_relativise_diagnose()` | model, batch, order | long table: offset, gain, robust z of log gain, batch gain ratio, gain drift with order | `batch`, `order` (may be `NULL`) |
| 4 | `sw_relativise_apply()` | linear sample data, model, unit of each file | relativised data, same layout | `channels`, `low`, `high`, `input_scale` |
| 5 | caller | relativised data | transform, then `sw_batch_scores()` / `sw_diagnose()` before vs after | cofactor etc. |

Step 1 uses one channel of each control: the channel it calibrates. k-means runs on
`split_transform(values)` (for example `asinh(v / 6000)`; `identity` gives the linear scale). The
cluster with the larger centre is "positive". The paper does not say which scale its k-means used.

Step 2 stops when a unit-channel has fewer than `min_events` events on either side, or when
`M <= m`. Step 4 stops when a file's unit has no row for a requested channel, and when
`input_scale` is not `"linear"`, because the formula is only defined on linear values.

`low` and `high` generalise the paper's targets (0 and 10^4):
`X~ = low + (X - m) / (M - m) * (high - low)`.

## 5. Diagnostics (numbers only, no pass or fail)

`sw_relativise_diagnose()` gives, per channel:
- `offset` (m) and `gain` (M - m) per unit;
- `gain_z`: robust z of log(gain) across units (`(v - median) / (1.4826 * MAD)`);
- `gain_batch_ratio`: a batch's median gain over the median gain of all units (when `batch` given);
- `gain_order_rho`: Spearman correlation of log(gain) with acquisition order (when `order` given,
  4 or more units).
A batch with a different detector gain shows as `gain_batch_ratio` far from 1. A drifting
instrument shows as `gain_order_rho`. Under A1, a donor-level gain outlier points at the control,
not the biology.

## 6. Validation

- Unit tests: recovery of known m and M; exact mapping m to `low` and M to `high`; every stop
  condition; RNG state restored; diagnose finds a doubled gain.
- Benchmark `experiments/benchmarks/09_relativisation/`: linear simulation with per-batch offset and
  gain, per-sample CD3 controls, and donor-to-donor CD3 variation (A1 broken on purpose). Score:
  per-channel pairwise EMD between batches on an arcsinh scale, before and after. Figure saved and
  inspected.
- Result (run `experiments/runs/2026-11-09_relativisation-v2`, figure `fig_relativisation.png`,
  inspected; 4 batches x 8 donors, 6 channels, 5 seeds). Median largest batch-pair EMD per channel,
  arcsinh scale:

  | CD3 donor log SD | observed | per sample | per batch (E1) | no instrument |
  |---|---|---|---|---|
  | 0 | 0.179 | 0.037 | 0.038 | 0.034 |
  | 0.1 | 0.234 | 0.044 | 0.038 | 0.040 |
  | 0.3 | 0.164 | 0.066 | 0.034 | 0.035 |

  Median per-file OTD (blur 0.75) tells the same story: at log SD 0.3 it is 0.0322 observed,
  0.0236 per sample, 0.0103 per batch and 0.0097 with no instrument; every solve converged.
  Per-sample relativisation inherits donor CD3 variation; E1 does not, because its reference is one
  fixed material per batch. The simulation has no technical variation inside a batch, which E1
  cannot remove.
- Real data: not possible on the Aurora run folder. Its files are OMIQ-transformed (A2 unknown)
  and it has no CD3 controls.

## 7. Extensions (approved by the user 2026-11-09, built)

| # | Extension | Function | Why it matters | Risk |
|---|---|---|---|---|
| E1 | Per-batch controls: one single stain of a stable reference marker per channel and batch (the spectral unmixing single stains qualify when each channel uses the same marker in every batch) | `sw_relativise_batch_map()` | 40-colour panels cannot carry a CD3 control per channel and sample | Removes batch offset and gain only, not sample-level; needs the controls stained and run with the batch |
| E2 | Internal reference from the full stain, for the reference channel only | `sw_relativise_internal_map()` | No extra tubes | Only the reference channel is calibrated; its negative median depends on the sample's cell mix |
| E3 | `sw_diagnose()` modules `relativise` and `otd` | `sw_diagnose(..., relativise, otd)` | One table for all checks | `sw_diagnose()` gained two required arguments; callers updated |
| E4 | Optimal-transport distance per file to a pooled reference (paper's OTD) | `sw_sinkhorn_divergence()`, `sw_otd_scores()` | Paper's main score | Clean rewrite of the debiased Sinkhorn divergence (cost `|x-y|^2/2`, `epsilon = blur^2`, annealing factor 0.7). Checked: a translation by `d` gives `|d|^2/2` at any blur; symmetric; 0 on equal sets. Small blur converges slowly and is flagged |
