# Notes: technical checks from the Colibri Cytometry blog (O. Burton)

Date: 2026-11-09. Source: https://www.colibri-cytometry.com/blog and its post sitemap (105 posts).
Read in full by five parallel fetch agents: 47 posts. Not read: wet-lab posts on blocking, cytokines,
fixatives, GFP and conjugation, and a few announcements. Quotes stay in the agent reports; this file
keeps only what we can build. Every threshold below is the author's number, for one instrument or
dataset; none becomes a package default.

## Checks on unmixed per-file data (no controls needed)

| Check | What it catches | Author's number | Source post | In package? |
|---|---|---|---|---|
| Negative centring: median of a known-negative population per channel near 0 | unmixing skew, AF mismatch | none | what-s-an-unmixing-matrix; dealing-with-multiple-autofluorescences | no |
| Negative width: robust SD or CV of the negative per channel, vs the cohort | spread, AF over/under-correction | none | dealing-with-multiple-autofluorescences; unmixing-the-tortoise-beats-the-hare | partly (`sd_z` is whole-file) |
| Hyper-negative share: events far below 0 per channel | AF over-correction, wrong control | none | identifying-and-fixing-errors-in-flow-data | partly (`negative_share`) |
| Correlation of two channels within the negative population | angled negative, spillover error | none | panel-design part 1; tips-and-tricks-optimizing-unmixing | partly (`corr_shift` is whole-file) |
| Diagonal pattern between Spark UV 387 and RB670/BUV661/BV650 channels | dye-dye interaction | none | spark-uv-387 posts | no |
| Parent-channel signal in tandem-positive cells (PE, APC, PerCP) | tandem breakdown | none | blocking-tandem-breakdown; tips-tricks-dealing-with-tandems | no |
| Signal against Time; event rate; low-SSC burst | clog, sensitivity shift mid-run | flow rate < 1 µl/min was a clog sign | unmixing-should-you-use-a-universal-negative | yes (`rate_cv`, `max_gap`, `time_trend`) |
| Events at the scatter maximum | off-scale SSC clipping | none | setting-up-scatter | yes (`saturation_share`, needs `range_max`) |
| Debris share below a scatter threshold | wrong threshold, junk | junk was ~90% of events in one run | setting-up-scatter | no |
| Cross-entropy test between samples or batches | batch effect in embeddings | none (n = 3 trend) | doing-more-with-less-reducing-batch-effects | no (we have EMD) |
| Cluster stability over seeds and cell counts | unstable clusters | 40% of cells moved between two seeds (Rphenograph-style, k = 30) | the-peculiarities-of-phenograph | no |

## Checks that need the controls or the spectra

| Check | What it catches | Author's number | Source post | In package? |
|---|---|---|---|---|
| Mixing-matrix condition number vs fluorophore count | junk control, hard panel | warn when above the fluorophore count; "over 80" before cleaning | autospectral-in-honeychrome; autospectral-cleaning | yes, wrapped (`sw_unmix_condition_number()`) |
| Cosine similarity between fluorophore pairs | unresolvable pairs | 0.93 PE vs SYG570; 0.91 BV510 vs AF | panel-design part 1; spectral-mixing-matrix-normalization | via AutoSpectral |
| Hotspot matrix | unmixing-dependent spread | < 4 fine, 4 to 6 check, > 6 attention | autospectral-full-workflow-example | via AutoSpectral (`max_sif` 4/6) |
| Control spectrum vs reference spectrum (cosine) | wrong or dirty control | no cutoff given | autospectral-new-features; tips-tricks-benchmarking-controls | partly (`sw_unmix_check_controls()`) |
| AF signature too close to a fluorophore | contaminated unstained | cosine >= 0.995 excluded | autospectral-new-features | no |
| Normalised spectrum values below 0 | mismatched negative | any value < 0 | what-s-in-a-spectral-matrix | no |
| Control events and positive events | weak control | >= 5000 events warned; "a couple hundred" positives; 500 brightest used | autospectral-full-workflow-example; tips-and-tricks-optimizing-unmixing | no |
| Day-0 matrix vs fresh controls over days | unmixing drift | stored references held "out to a month" | stabilizing-your-unmixing(-drift) | no |
| Control variability projected through the unmixing matrix | where errors will appear | no formula; covariance preferred | predicting-unmixing-errors | no |
| Controls acquired fast vs slow | resolution loss | complexity 36.47 vs 35.09 | unmixing-the-tortoise-beats-the-hare | no |

## Reference-population relativisation

The blog's references page points to Glehr et al. 2026 (CD3 relativisation). Spec and status:
`plans/2026-11-09_reference-relativisation.md`. Built: `sw_relativise_*()`.

## Gaps the blog states

- No post gives a method to separate batch from biology beyond overnight staining and matched
  controls. Honeycluster lists batch correction as missing.
- The author's thresholds are visual or instrument-specific. Most checks above have no cutoff.
