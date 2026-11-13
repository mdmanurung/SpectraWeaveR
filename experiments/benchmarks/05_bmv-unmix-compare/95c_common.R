# Purpose: shared setup for the 95c_restyle_*.R scripts, which redraw the certified non-spike figures from
#          the common-gate (SpectroFlo labels, clean3) outputs of batches 0-4 (plans/common_gates_plan.md
#          G7a, G11, G13, G15). Loads the house theme and donor encoding, names the cut, and records each figure
#          in a per-script manifest part that 95c_manifest.R merges into MANIFEST_restyle_common.csv.
#          One restyle root per dataset: <SW_RUN_DATE>_bmv<run_tag>-restyle-common. The first script that
#          sources this file creates it; the other 95c scripts of the same run reuse it.
# Run:     SW_DATASET=3d SW_RUN_DATE=YYYY-MM-DD, then
#          source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R") from the repo root.
# Inputs:  datasets/dataset.R and datasets/<SW_DATASET>.json   DS, cfg(), run_in(), run_out(), samples_tbl()
#          experiments/benchmarks/05_bmv-unmix-compare/report_figures.R   donor palette, panel builders, save_panel
#          run folder key "restyle" (3d only, old certified restyle; read by 95c_manifest.R)
# Outputs: none; defines RUNS, OUT_ROOT, CUT_TXT, POS_LABEL, DONORS, DONOR_PAL, DONOR_SHAPES, fresh_dir(), add_manifest()
#          There is no scalar CUT: the cut differs by arm, sample and channel, so a script calls
#          cut_of(arm, sample, channel, "pos") for each cell.

source("experiments/benchmarks/05_bmv-unmix-compare/report_figures.R")
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")

# Parameters ---------------------------------------------------------------------------------------
RUNS <- "experiments/runs"
# The restyle root is shared by every 95c script, so it may already exist (run_out() stops on that).
.root <- file.path(
  RUNS,
  paste0(Sys.getenv("SW_RUN_DATE"), "_bmv", cfg("run_tag"), "-restyle-common")
)
OUT_ROOT <- if (dir.exists(.root)) .root else run_out("restyle-common")
# cut wording on every label and title, from the dataset's positive-cut rule
.pos <- cfg("pos_rule")
CUT_TXT <- switch(
  .pos$type,
  fixed = paste("arcsinh", .pos$value),
  control_quantile = paste0(
    "per-donor ",
    .pos$control_condition,
    " q",
    round(100 * .pos$prob, 2)
  ),
  mindensity = "mindensity",
  stop("95c_common.R: unknown pos_rule type '", .pos$type, "'", call. = FALSE)
)
# readout label of the positive cut in the 51c, 87 and 87b gap tables (same rule as POS_LABEL there); CUT_TXT
# is for titles and captions only
POS_LABEL <- if (identical(.pos$type, "fixed")) {
  sprintf("arcsinh %s", .pos$value)
} else {
  "the positive cut"
}

# Donor encoding: one colour and one shape per replicate, from the sample table. 3D and EXV key it on the
# sample (well code); with conditions (24H) the replicate is the donor.
# report_figures.R fixed the 3D five; these replace them in the global environment.
DONORS <- if (is.null(DS$conditions)) {
  samples_tbl()$sample
} else {
  unique(samples_tbl()$donor)
}
PAIRED <- sprintf(
  "Each point is one donor (colour and shape); the bar is the median of %d donors.",
  length(DONORS)
) # report_figures.R words this for 5 donors
DONOR_PAL <- categorical_palette(DONORS)
DONOR_SHAPES <- setNames(rep_len(c(16, 17, 15, 18, 25), length(DONORS)), DONORS)

# Conditions (24H only): a per-donor summary or facet adds "condition" to its keys; 3D and EXV have none.
COND <- if (is.null(DS$conditions)) NULL else "condition"
need_cond <- function(d, script) {
  if (!is.null(COND) && !"condition" %in% names(d)) {
    stop(
      script,
      ": table has no 'condition' column, needed for dataset ",
      DS$dataset,
      call. = FALSE
    )
  }
  invisible(d)
}

# Donor key and condition of a table with a `sample` column. 3D and EXV: the sample is the donor key and
# there is no condition. 24H: the sample joins the sample table (one row each, all must match) for the
# donor and, if the table has none, the condition (inner join; an unmatched sample stops).
sample_keys <- function(d) {
  if (is.null(COND)) {
    return(mutate(d, donor_key = .data$sample))
  }
  keys <- select(samples_tbl(), "sample", donor_key = "donor", "condition")
  if ("condition" %in% names(d)) {
    keys <- select(keys, -"condition")
  }
  inner_join(
    d,
    keys,
    by = "sample",
    relationship = "many-to-one",
    unmatched = c("error", "drop")
  )
}

# One panel per condition for a single-panel figure; NULL (adds nothing) without conditions.
cond_facet <- function() {
  if (is.null(COND)) {
    return(NULL)
  }
  facet_wrap(
    vars(.data$condition),
    nrow = 1,
    axes = "all",
    axis.labels = "margins"
  )
}

# A target the dataset lacks (for example Th2, IL-2 or IL-10 on EXV) stops the script; it is never skipped
# and never falls back to the 3D value. Returns cfg(group)[[key]] (a dye, label or channel).
need_target <- function(script, key, group = "targets") {
  v <- cfg(group)[[key]]
  if (is.null(v)) {
    stop(
      script,
      ": needs ",
      key,
      ", missing in dataset ",
      DS$dataset,
      call. = FALSE
    )
  }
  v
}

# Cut lines of a figure: one row per distinct (key, cut), where key is the facet column that holds the
# sample id and the cut comes from cut_of() for each arm in the table. Equal cuts across arms collapse.
cut_lines <- function(d, key, channel, kind = "pos") {
  d |>
    distinct(across(all_of(c("arm", key)))) |>
    mutate(
      cut = purrr::map2_dbl(.data$arm, .data[[key]], \(a, s) {
        cut_of(a, s, channel, kind)
      })
    ) |>
    distinct(across(all_of(c(key, "cut"))))
}

# Partial redraw: REDRAW_ONLY="key1,key2" reruns just those figures into the existing folder, overwriting
# their PNG, PDF and CSV (visual-pass fixes, 2026-11-11). Unset, a script writes a new folder as before.
REDRAW_ONLY <- Sys.getenv("REDRAW_ONLY")
ONLY <- if (nzchar(REDRAW_ONLY)) strsplit(REDRAW_ONLY, ",")[[1]] else NULL
only_key <- function(key) is.null(ONLY) | key %in% ONLY

# A section folder is new for each script run; a rerun deletes it by hand first.
fresh_dir <- function(section) {
  out <- file.path(OUT_ROOT, section)
  if (dir.exists(out) && !is.null(ONLY)) {
    return(out)
  }
  if (dir.exists(out)) {
    stop(
      "fresh_dir: ",
      out,
      " exists; remove it by hand to rerun.",
      call. = FALSE
    )
  }
  dir.create(out, recursive = TRUE)
  out
}

# One manifest line per figure: old figure path (as in audit_freshness.csv / the old restyle), the new
# figure path, the saved table or run dir it was drawn from, and the script. Parts merge in 95c_manifest.R.
add_manifest <- function(
  script,
  old_path,
  new_path,
  source,
  status = "redrawn",
  note = ""
) {
  dir.create(
    file.path(OUT_ROOT, "manifest_parts"),
    showWarnings = FALSE,
    recursive = TRUE
  )
  part <- file.path(OUT_ROOT, "manifest_parts", paste0(script, ".csv"))
  row <- tibble::tibble(
    old_figure = old_path,
    new_figure = new_path,
    source = source,
    script = script,
    status = status,
    note = note
  )
  # a partial redraw replaces the part's earlier row for the same figure (by new path, else old path)
  if (file.exists(part)) {
    old <- readr::read_csv(part, col_types = readr::cols(.default = "c"))
    same <- if (is.na(new_path)) {
      !is.na(old$old_figure) &
        old$old_figure == old_path &
        is.na(old$new_figure)
    } else {
      !is.na(old$new_figure) & old$new_figure == new_path
    }
    row <- dplyr::bind_rows(
      old[!same, ],
      dplyr::mutate(row, dplyr::across(dplyr::everything(), as.character))
    )
  }
  readr::write_csv(row, part)
}

# Arm display names -------------------------------------------------------------------------------------
# ARM_LABELS (arm_labels.R, 11 main arms) plus the crossed and swapped arms that report_figures.R leaves as
# raw codes. Every 95c figure labels arms through arm_label(), so one arm reads the same everywhere.
ARM_EXTRA <- c(
  swap_af102 = "1.8.8 + 1.0.2 AF",
  swap_var102 = "1.8.8 + 1.0.2 variants",
  swap_af102_var102 = "1.8.8 + 1.0.2 AF and variants",
  rev_af188 = "1.0.2 spectra, 1.8.8 AF",
  var102_only_SBUV510 = "1.8.8 + 1.0.2 SBUV510 variants",
  var102_only_Alexa_Fluor_405 = "1.8.8 + 1.0.2 AF405 variants",
  var102_only_SBUV510__Alexa_Fluor_405 = "1.8.8 + 1.0.2 SBUV510 and AF405 variants",
  opt188_pregate_joint_p1af1 = "opt188 (optimized)",
  af405proof_base188 = "1.8.8 (all 1.8.8 variants)",
  af405proof_dose_000 = "dose 0: AF405 = 1.0.2 reference only",
  af405proof_dose_025 = "dose 0.25",
  af405proof_dose_050 = "dose 0.5",
  af405proof_ki_bv510 = "knock-in: 1.8.8 AF405 + BV510 part",
  af405proof_ko_bv510 = "knock-out: 1.0.2 AF405 minus BV510 part",
  af405proof_recentre = "1.0.2 spread around 1.8.8 reference"
)
ARM_DISPLAY <- c(ARM_LABELS, ARM_EXTRA)
arm_label <- function(x) {
  if_else(x %in% names(ARM_DISPLAY), unname(ARM_DISPLAY[x]), x)
}

# Heatmaps with untested cells --------------------------------------------------------------------------
# A cell with no score (marker not scored in that population or arm) is added as NA, so it draws grey80
# instead of an empty gap. Factor levels keep the caller's order.
complete_grid <- function(d, row, col) {
  for (key in c(row, col)) {
    if (!is.factor(d[[key]])) d[[key]] <- factor(d[[key]], unique(d[[key]]))
  }
  tidyr::complete(d, .data[[row]], .data[[col]])
}
# Same diverging fill as plot_heatmap(), with na.value = "grey80"; added after it, it replaces its scale.
heat_fill <- function(limit, name) {
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-limit, limit),
    oob = scales::squish,
    na.value = "grey80",
    name = name
  )
}
NOT_TESTED <- "grey = not tested"
BOX_RULE <- "boxed: donor-median rise \u2265 0.10, rising in all 5 donors"
