# Purpose: shared setup for the three 95c spike restyles (95c_restyle_spike.R, 95c_restyle_events_spike_raw.R,
#          95c_restyle_events_cospike_raw.R): the common-gate spike run paths, marker dyes, setup display names
#          through arm_label(), and the REDRAW_ONLY switch for spike figures (plans/common_gates_plan.md G5, G11).
#          The dataset comes from SW_DATASET (3d, 24h, exv) with SW_RUN_DATE; spike runs are registry keys.
# Run:     source("experiments/benchmarks/05_bmv-unmix-compare/95c_spike_common.R") from the repo root.
# Inputs:  experiments/benchmarks/05_bmv-unmix-compare/95c_common.R   RUNS, ONLY, fresh_dir(), add_manifest(), arm_label()
#          datasets/<SW_DATASET>.json   spike_marker_dyes (marker label -> dye, in figure order)
# Outputs: none; defines want(), SRC(), MARKERS, setup_label(), SETUP_ORDER, SETTING_ORDER, FIG_LEVEL, LOSS_FRAC, SEED

source("experiments/benchmarks/05_bmv-unmix-compare/95c_common.R")

# REDRAW_ONLY=spike redraws every spike figure; a list of figure stems redraws just those.
want <- \(stem) is.null(ONLY) || "spike" %in% ONLY || any(stem %in% ONLY)

# Parameters ---------------------------------------------------------------------------------------
FIG_LEVEL <- 1 # level the figures show (user, 2026-11-10: +1 arcsinh), inherited (G5)
LOSS_FRAC <- 0.5 # "losing most of the spike" = recovery below 0.5, inherited (G5)
SEED <- 42 # figure jitter, inherited (G5)
# run is a registry key (datasets/<SW_DATASET>_runs.csv), e.g. "spike-in-cg"
SRC <- \(run, file) file.path(run_in(run), file)
# marker = dye, from the control file (as 88, 88b, 88e); a dataset with a different panel lists its own
MARKERS <- unlist(cfg("spike_marker_dyes"))
# setup code -> arm key of ARM_DISPLAY; afrefine has no arm key and is named here
SETUP_ARM <- c(
  c188 = "as188_joint",
  afrefine = NA,
  opt188 = "opt188_pregate_joint_p1af1",
  s188_c102 = "s188_c102",
  s102_c102 = "as102_var"
)
SETUP_NAME <- c(afrefine = "AS 1.8.8 joint, refined AF")
setup_label <- \(x) {
  stopifnot("unknown setup code" = all(x %in% names(SETUP_ARM)))
  arm <- unname(SETUP_ARM[x])
  if_else(is.na(arm), unname(SETUP_NAME[x]), arm_label(arm))
}
SETUP_ORDER <- setup_label(names(SETUP_ARM)) # figure order, top to bottom
# n.passes / n.af.passes settings of 88h
SETTING_ORDER <- c(
  p1_af1 = "1 / 1 (reference)",
  p5_af1 = "5 / 1",
  p1_af5 = "1 / 5",
  p5_af5 = "5 / 5"
)
