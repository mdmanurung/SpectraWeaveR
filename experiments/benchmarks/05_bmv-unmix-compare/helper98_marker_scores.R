# Purpose: common-gates copy of helper75_marker-scores.R (plans/common_gates_plan.md, batch 3), shared by
#          75_all_marker_screen_common.R, 76_per_dye_channel_map_common.R and 86_th2_proof_arms_common.R. Rows come
#          from the clean3 manifest of the common gates; the negative cut and the above-cut readout are separate
#          (rSD of negatives, G13; share above, G7a). Cuts come from cut_of() per arm, sample and channel
#          (arcsinh, cofactor 6000); a caller that passes one number for every channel gets that number.
#          Original purpose: shared by 75 and 76. Loads an arm's events for one donor
#          in as188_joint's cleaned row order and scores every channel's negatives with the 39b / 51 metrics.
# Run:     source("experiments/benchmarks/05_bmv-unmix-compare/helper98_marker_scores.R") after the caller
#          has set COFACTOR and SQUEEZE_FRAC. SW_DATASET (3d, 24h or exv) is required; this file sources
#          datasets/dataset.R itself. score_channels() takes neg_cut, pos_cut and ref_neg_cut (one value or a
#          vector named by channel); left out, it reads the caller's globals NEG_CUT and POS_CUT.
# Inputs:  registry key common-gated: manifest.csv   clean3 FCS per arm [97c]
#          <data_dir>/autoSpectral_v1.8.8/AutoSpectral_unmixed   raw as188_joint
#          <data_dir>/unmix_cells/<arm>/   raw unmixed arms (24H: one subfolder per condition) [37, 45, 48 scripts]
#          datasets/<ds>_samples.csv   sample and file stem [00_samples.R]
# Outputs: none (functions only)

source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), cut_of(), samples_tbl()

POOLED_MANIFEST <- file.path(run_in("common-gated"), "manifest.csv") # batch 0: clean3 files of every arm
RAW_REF_DIR <- data_path("autoSpectral_v1.8.8/AutoSpectral_unmixed")
UNMIX_DIR <- data_path("unmix_cells")
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; maps raw rows to cleaned rows

read_x <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
}
key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))

# Raw unmixed file of one sample in one folder. Only files whose stem is a sample stem count, so stray .fcs files
# (e.g. controls/ inside an arm folder) are ignored; a missing or doubled file for the sample stops the run.
# recursive: 24H arm folders hold one subfolder per condition. AutoSpectral appends " AutoSpectral" or " AutoSpectral_v1-8-8" to the stem.
raw_file <- function(dir, sample_id) {
  all_stems <- samples_tbl()
  sample_row <- all_stems |>
    filter(sample == sample_id) |>
    select(sample, stem)
  arm_files <- tibble(
    file = list.files(
      dir,
      pattern = "[.]fcs$",
      recursive = TRUE,
      full.names = TRUE
    )
  ) |>
    mutate(
      stem = sub(
        " AutoSpectral.*$",
        "",
        tools::file_path_sans_ext(basename(file))
      )
    )
  # files with no sample stem (controls/, objects) are left out, and said so
  n_stray <- sum(!arm_files$stem %in% all_stems$stem)
  if (n_stray > 0) {
    message(sprintf(
      "raw_file(): %d FCS file(s) with no sample stem left out of %s",
      n_stray,
      dir
    ))
  }
  hit <- arm_files |>
    filter(stem %in% all_stems$stem) |>
    inner_join(sample_row, by = "stem", relationship = "many-to-one")
  if (nrow(hit) != 1) {
    stop(
      sprintf(
        "raw_file(): %d files for sample %s (stem %s) in %s, need exactly 1",
        nrow(hit),
        sample_id,
        paste(sample_row$stem, collapse = ", "),
        dir
      ),
      call. = FALSE
    )
  }
  hit$file
}

# Config marker_parent with its fields checked: a missing field would silently change which populations are scored.
# default = populations scored for every marker; by_marker = markers with their own parent; no_marker = parent of
# a channel without a marker; map_pop = the one parent of a multi-parent channel in 76; also_non_t = markers
# assigned to T cells but also expressed outside them (a note in marker_parent.csv).
marker_parent_cfg <- function() {
  mp <- cfg("marker_parent")
  absent <- setdiff(
    c("default", "no_marker", "map_pop", "also_non_t"),
    names(mp)
  )
  if (length(absent) > 0) {
    stop(
      "helper98_marker_scores.R: config marker_parent lacks ",
      paste(absent, collapse = ", "),
      " in dataset ",
      DS$dataset,
      call. = FALSE
    )
  }
  mp
}

# Cut per channel for one arm and sample, as a vector named by channel. kind: "pos" (share above) or "neg"
# (negative population). cut_of() stops when the dataset has no cut for a channel.
cuts_by_channel <- function(arm, sample, channels, kind) {
  map_dbl(set_names(channels), \(ch) cut_of(arm, sample, ch, kind))
}

# One arm and donor, all channels, linear, rows = as188_joint's cleaned rows. Pooled arms are read from
# their cleaned files; other arms are raw unmixed files masked to as188_joint's cleaned rows (as 51 does).
arm_matrix <- function(arm, sample, channels, pooled) {
  ref_clean <- pooled$file[
    pooled$arm == "as188_joint" & pooled$sample == sample
  ]
  ref_key <- key_of(read_x(ref_clean))
  if (arm %in% pooled$arm) {
    e <- read_x(pooled$file[pooled$arm == arm & pooled$sample == sample])
  } else {
    raw <- read_x(raw_file(file.path(UNMIX_DIR, arm), sample))
    raw_ref_key <- key_of(read_x(raw_file(RAW_REF_DIR, sample)))
    if (nrow(raw) != length(raw_ref_key)) {
      stop(
        sprintf("%s %s: raw rows differ from as188_joint raw", arm, sample),
        call. = FALSE
      )
    }
    e <- raw[raw_ref_key %in% ref_key, , drop = FALSE]
  }
  if (!identical(key_of(e), ref_key)) {
    stop(
      sprintf("%s %s: rows do not pair with as188_joint", arm, sample),
      call. = FALSE
    )
  }
  e[, channels, drop = FALSE]
}

# 39b / 51 readouts per channel on given rows: Nguyen rSD (q84 - q50) of events below arcsinh neg_cut, its ratio
# to as188_joint on the same rows (that arm's own cut, ref_neg_cut), share with |x| < SQUEEZE_FRAC x as188_joint's
# rSD, share above arcsinh pos_cut. A cut is one value or a vector named by channel.
nguyen <- \(v) unname(diff(quantile(v, c(0.5, 0.84))))
cut_for <- \(cuts, ch) if (length(cuts) == 1) cuts else cuts[[ch]]
score_channels <- function(
  x,
  x_ref,
  rows,
  neg_cut = NEG_CUT,
  pos_cut = POS_CUT,
  ref_neg_cut = neg_cut
) {
  colnames(x) |>
    set_names() |>
    map(\(ch) {
      lin <- x[rows, ch]
      lin_ref <- x_ref[rows, ch]
      ref_rsd <- nguyen(lin_ref[
        asinh(lin_ref / COFACTOR) < cut_for(ref_neg_cut, ch)
      ])
      tibble(
        n = length(rows),
        rsd_neg = nguyen(lin[asinh(lin / COFACTOR) < cut_for(neg_cut, ch)]),
        ref_rsd = ref_rsd,
        rsd_ratio = rsd_neg / ref_rsd,
        squeezed_share = mean(abs(lin) < SQUEEZE_FRAC * ref_rsd),
        pos_share = mean(asinh(lin / COFACTOR) > cut_for(pos_cut, ch))
      )
    }) |>
    list_rbind(names_to = "channel")
}
