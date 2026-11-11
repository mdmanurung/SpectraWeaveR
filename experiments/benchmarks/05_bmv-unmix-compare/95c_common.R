# Purpose: shared setup for the 95c_restyle_*.R scripts, which redraw the certified non-spike figures from
#          the common-gate (SpectroFlo labels, clean3, arcsinh 1.5) outputs of batches 0-4 (plans/common_gates_plan.md
#          G7a, G11, G13, G15). Loads the house theme and donor encoding, fixes the cut, and records each figure
#          in a per-script manifest part that 95c_manifest.R merges into MANIFEST_restyle_common.csv.
# Run:     source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R") from the repo root.
# Inputs:  experiments/benchmarks/05_bmv-unmix-compare/report_figures.R   donor palette, panel builders, save_panel
# Outputs: none; defines RUNS, OUT_ROOT, CUT, CUT_TXT, fresh_dir(), add_manifest()

source("experiments/benchmarks/05_bmv-unmix-compare/report_figures.R")

# Parameters ---------------------------------------------------------------------------------------
RUNS <- "experiments/runs"
OUT_ROOT <- "experiments/runs/2026-11-11_bmv-restyle-common"
OLD_ROOT <- "experiments/runs/2026-11-11_bmv-restyle"
CUT <- 1.5 # cytokine cut, arcsinh, cofactor 6000 (G7a); every label and line reads this
CUT_TXT <- "arcsinh 1.5"

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
