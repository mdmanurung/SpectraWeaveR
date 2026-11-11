# Second multi-batch dataset: search result (G42, 2026-11-11)

Criteria (D5, G42): several batches; the same reference material in every batch; batch not
nested in donor, at least 4 donors per batch; per-batch single stains; raw linear FCS without an
application; spectral preferred, conventional accepted, mass cytometry excluded. Nothing is
downloaded before the user approves an accession.

No candidate meets every criterion on verified evidence. Every claim below cites the page it
came from; "unknown" means not read on any fetched page.

## Shortlist

| # | Dataset | Accession | Instrument | Batches | Bridge each batch | Donors per batch | Single stains | Access |
|---|---|---|---|---|---|---|---|---|
| 1 | den Braanker 2021 healthy-control spectral set | FR-FCM-Z4KT | Cytek Aurora, 5-laser, 31 colours | 6 | yes, 2 QC samples per run | unknown | one set for all batches, deposit unknown | FlowRepository, public; licence unknown |
| 2 | Public COVID-19 spectral set cited in the same paper | FR-FCM-Z3WR | unknown (36-colour panel) | unknown | unknown | samples randomised over batches | unknown | unknown |

Evidence for 1, from https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.768113/full:
"control samples in 6 batches measured at different dates"; "Every run, two of the same quality
control samples were acquired to correct for batch effects later on"; "measured on the 5-laser
Aurora spectral flow cytometer"; "Reference controls were acquired once and used for unmixing of
the multiple batches"; "FCS files are available from Flow Repository, FR-FCM-Z4KT". The paper
gives no donor or file count per batch. Risk: the deposit may hold only about 12 files (2 per
batch), which fails the 4-donors-per-batch rule.

Evidence for 2, same paper: "In the public dataset, FR-FCM-Z3WR, they randomly distributed their
samples over the batches." Nothing else verified.

## Record pages read over plain HTTP (user approved option b, 2026-11-11)

The HTTPS pages fail on an expired certificate; the same pages load over plain HTTP. The
file-list endpoint (`/list/<id>`) returns 403.

- FR-FCM-Z4KT (http://flowrepository.org/id/FR-FCM-Z4KT): "Spectral flow cytometry data from
  healthy controls with a 31-color panel"; experiment dates 2020-04-06 to 2021-04-06; "Cells
  were used as single stainings". The page lists no files and no batch or donor variables, so
  donors per batch stay unknown until the files are downloaded.
- FR-FCM-Z3WR (http://flowrepository.org/id/FR-FCM-Z3WR): "human PBMC in COVID-19", "36-color
  spectral flow cytometry", experiment dates 2020-05-01 to 2020-12-31, Duke University; "Single
  stained cells/beads". The page lists 84 sample files annotated by condition (healthy 10,
  exposed, infected, hospitalized), sex and time point, with repeat time points per donor (for
  example `5ABDCB_0`, `_14`, `_28`). It has no batch variable and no bridge sample; batches could
  only come from the `$DATE` keyword inside the files. Fails criterion 2 (bridge); criterion 3
  is plausible through repeat donors but unverified.

## Not verified

- The FlowRepository pages could not be fetched (expired TLS certificate), so file lists,
  donors per batch and licences of Z4KT and Z3WR are unknown. Check them in a browser.
- ImmPort SDY2540 (Aurora PBMC files) and SDY2268 / SDY2277 appeared only in search snippets.

## Rejected (pages fetched)

- Zenodo 17138653: synthetic batch effects, 15 files; not a real multi-batch set.
- Zenodo 15723074: Aurora, 7 donors, 56 files; no batches or bridge stated.
- Zenodo 23934346: files restricted.
- CytoNorm sets (FR-FCM-Z247 and others): mass cytometry, excluded.
- Frontiers 2026 17-colour PBMC study (bridge aliquot in each run): data only on request.
