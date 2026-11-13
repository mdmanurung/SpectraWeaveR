# Purpose: check every marker and dye named in a dataset config against that dataset's own control file.
#          The 3D, 24H and ex vivo panels differ, so a fluor is never assumed to carry the same marker in
#          another dataset (user 2026-11-12). A config entry that names a marker on a dye the panel gives to
#          another marker is a mismatch and stops the run.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/datasets/check_panel.R
# Inputs:  <data_dir>/autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv   fluorophore, marker [zip]
#          datasets/<ds>.json                                                       config values to check
# Outputs: a table on stdout, one row per config entry: item, config, panel, status
#          status: ok | mismatch (stops, exit 1) | slot_differs (dye kept as a slot, the panel says another marker)
#                  | absent (dye or marker not in the panel) | skipped (config value is null: not yet decided)
#          Nothing is written to disk.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# panel from the dataset's own control file: one row per dye, marker names normalised for matching
ctrl <- read_csv(
  data_path("autoSpectral_v1.8.8/Reference Group Clean/fcs_control_file.csv"),
  show_col_types = FALSE
)
norm <- \(x) str_to_lower(str_replace_all(x, "[^A-Za-z0-9]", ""))
panel <- ctrl |>
  filter(!fluorophore %in% c("AF", "Negative")) |>
  distinct(dye = fluorophore, marker) |>
  mutate(
    marker_key = map(marker, \(m) {
      unique(c(norm(m), norm(str_remove(m, "^.*_"))))
    }),
    channel = paste0(dye, "-A")
  )
stopifnot("one marker per dye" = !anyDuplicated(panel$dye))
dye_of_marker <- \(m) panel$dye[map_lgl(panel$marker_key, \(k) norm(m) %in% k)]
marker_of_dye <- \(d) panel$marker[panel$dye == d]
dye_of_channel <- \(ch) str_remove(ch, "-A$")

rows <- list()
add <- \(item, config, panel_says, status) {
  rows[[length(rows) + 1]] <<- tibble(
    item = item,
    config = config,
    panel = panel_says,
    status = status
  )
}
# a marker is expected on a dye: ok when the panel puts exactly that marker there
check_marker_on_dye <- \(item, marker, dye, kind = "marker") {
  if (!dye %in% panel$dye) {
    return(add(item, paste(marker, "on", dye), "dye not in panel", "absent"))
  }
  got <- marker_of_dye(dye)
  ok <- norm(marker) %in% unlist(panel$marker_key[panel$dye == dye])
  add(
    item,
    paste(marker, "on", dye),
    got,
    if (ok) {
      "ok"
    } else if (kind == "slot") {
      "slot_differs"
    } else {
      "mismatch"
    }
  )
}

# 1. lineage channels: marker name is the key, the channel must carry it
ch <- cfg("channels")
walk(c("CD3", "CD4", "CD8"), \(m) {
  check_marker_on_dye(paste0("channels$", m), m, dye_of_channel(ch[[m]]))
})

# 2. cytokine channels: a slot named after a marker. In a panel without it the slot is a dye kept on purpose
slots <- c(
  TNF_slot = "TNF",
  IFNg_slot = "IFNg",
  Th2_slot = "Th2",
  IL10_slot = "IL-10",
  IL2 = "IL-2",
  IL10 = "IL-10",
  Th2 = "Th2"
)
walk(names(slots), \(k) {
  if (!is.null(ch[[k]])) {
    check_marker_on_dye(
      paste0("channels$", k),
      slots[[k]],
      dye_of_channel(ch[[k]]),
      kind = "slot"
    )
  }
})

# 3. target labels: a label must be what the panel calls that dye
tg <- cfg("targets")
walk(
  list(
    c("tnf_dye", "tnf_label"),
    c("ifng_dye", "ifng_label"),
    c("th2_dye", "th2_label")
  ),
  \(p) {
    if (!is.null(tg[[p[[1]]]]) && !is.null(tg[[p[[2]]]])) {
      check_marker_on_dye(
        paste0("targets$", p[[2]]),
        tg[[p[[2]]]],
        tg[[p[[1]]]]
      )
    }
  }
)
if (!is.null(tg$streak_pair)) {
  walk(tg$streak_pair, \(ch1) {
    add(
      "targets$streak_pair",
      ch1,
      if (dye_of_channel(ch1) %in% panel$dye) {
        marker_of_dye(dye_of_channel(ch1))
      } else {
        "dye not in panel"
      },
      if (dye_of_channel(ch1) %in% panel$dye) "ok" else "absent"
    )
  })
}

# 4. marker -> dye / channel lists
dye_lists <- list(
  `spike_panels` = NULL,
  spike_marker_dyes = DS$spike_marker_dyes
)
for (nm in names(DS$spike_panels)) {
  dye_lists[[paste0("spike_panels$", nm)]] <- DS$spike_panels[[nm]]
}
walk(names(dye_lists), \(k) {
  v <- dye_lists[[k]]
  if (is.null(v)) {
    return(invisible())
  }
  iwalk(v, \(dye, marker) {
    check_marker_on_dye(paste0(k, "$", marker), marker, dye)
  })
})
for (k in c("cluster_markers_other", "lineage_channels_other")) {
  v <- DS[[k]]
  if (!is.null(v)) {
    iwalk(v, \(channel, marker) {
      check_marker_on_dye(
        paste0(k, "$", marker),
        marker,
        dye_of_channel(channel)
      )
    })
  }
}

# 5. marker names only: they must exist in the panel
umap <- DS$umap_markers
if (!is.null(umap)) {
  walk(unique(unlist(umap)), \(m) {
    add(
      "umap_markers",
      m,
      paste(marker_of_dye(dye_of_marker(m)), collapse = ","),
      if (length(dye_of_marker(m))) "ok" else "absent"
    )
  })
}
mp <- DS$marker_parent
if (!is.null(mp)) {
  mnames <- unique(c(names(mp$by_marker), unlist(mp$also_non_t)))
  walk(mnames, \(m) {
    add(
      "marker_parent",
      m,
      paste(marker_of_dye(dye_of_marker(m)), collapse = ","),
      if (length(dye_of_marker(m))) "ok" else "absent"
    )
  })
}

# 6. config keys that are still null: nothing to check
for (k in c(
  "spike_panels",
  "spike_marker_dyes",
  "cluster_markers_other",
  "lineage_channels_other",
  "umap_markers",
  "marker_parent"
)) {
  if (is.null(DS[[k]])) add(k, "null", "-", "skipped")
}

res <- bind_rows(rows) |>
  mutate(across(everything(), \(x) map_chr(x, \(y) paste(y, collapse = ","))))
cat(sprintf(
  "dataset %s: %d dyes in the control file\n",
  DS$dataset,
  nrow(panel)
))
print(count(res, status), n = Inf)
print(filter(res, status != "ok"), n = Inf, width = 200)
if (any(res$status %in% c("mismatch"))) {
  stop(
    "check_panel.R: config names a marker on a dye the panel gives to another marker (dataset ",
    DS$dataset,
    ")",
    call. = FALSE
  )
}
