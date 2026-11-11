# Purpose: score ablation arms on the common SpectroFlo gates (plans/common_gates_plan.md, batch 2). Copy of
#          51_score_ablation.R that only READS the common gated dir: it never gates (no 19a), never writes
#          gs/ or manifest, and stops when an arm is missing. Every arm carries the same SpectroFlo labels
#          by row (G8), so "own gate" and "as188_joint rows" are the same rows. Each call writes its own
#          run dir, so parallel calls cannot race.
#          Changes against 51 (all from the plan):
#            CUT 1 -> 1.5 (G7a): TNF+/IFNg+ cut, the negatives that define the rSD, and "above cut".
#            The squeeze cut (SQUEEZE_FRAC x as188_joint negative rSD) is refitted on the common pops;
#            squeeze_cut.csv gives it at cut 1.5 (used) and, for the old-vs-new comparison, at cut 1.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/51c_score_ablation_common.R <run_dir> <arms>
#            from the repo root. <arms> is a comma-separated list of arm names in the common gated dir.
#            <run_dir> must not exist. Example (Block A):
#            ... 51c_score_ablation_common.R experiments/runs/2026-11-11_bmv-common-ablation \
#                swap_spec102,swap_af102_var102,rev_var188,rev_af188,swap_thr102,regen188,bead_wide
# Inputs:  experiments/runs/2026-11-11_bmv-common-gated/README.md, manifest.csv, gs/<arm>   [batch 0]
#          nodes /nonDebris/singlets/live/CD3pos/{CD4T,CD8T,DP,DN}, children TNFpos and IFNgpos (not used:
#          the cytokine cut is the fixed 1.5, computed here)
# Outputs: <run_dir>/pct_pos.csv        TNF+ and IFNg+ % at cut 1.5 in CD4T and CD8T
#          <run_dir>/tnf_among_ifng.csv TNF+ % among IFNg+ cells
#          <run_dir>/scores.csv         squeezed share, share above 1.5, rSD of negatives; BV510 and BV570
#          <run_dir>/squeeze_cut.csv    squeeze cut per donor, pop, cytokine at cut 1.5 and at cut 1
#          <run_dir>/neg_spread.csv     Nguyen rSD ratio to as188_joint, TNF, IFNg, Th2, IL-10, shared negatives
#          <run_dir>/gap_closed.csv     per donor fraction of the as188_joint -> s102_c188 gap closed
#          <run_dir>/labels.csv         resembles / contributes / null per arm
#          <run_dir>/fig_gap_closed.{pdf,png,csv} gap closed per arm, donors as points
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))

# Parameters ------------------------------------------------------------------------------
COMMON <- "experiments/runs/2026-11-11_bmv-common-gated" # batch 0 output, read only
BASE_ARMS <- c("as188_joint", "s102_c188", "as102_var") # D1: start, chain target, final target
REF_ARM <- "as188_joint"
TARGET_ARM <- "s102_c188" # gap is as188_joint -> s102_c188 (D3)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
CUT <- 1.5 # fixed cytokine cut, arcsinh, every arm and donor (G7a, user, 2026-11-11)
CUT_OLD <- 1 # the pre-G7a cut, only to report the old squeeze cut on the new pops
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
CYTOKINES <- c(
  TNF = "BV510-A",
  IFNg = "BV570-A",
  Th2 = "eFluor 450-A",
  `IL-10` = "BB700-A"
)
STORED_ARCSINH <- c("BV510-A", "BV570-A") # 19a transforms these; the others are linear in the GatingSet
POPS <- c("CD4T", "CD8T")
POP_PATH <- "/nonDebris/singlets/live/CD3pos/" # full node path prefix, so CD4T/CD8T are unambiguous
PAIR_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # G8: rows must pair across arms
DONORS <- c("A3", "A6", "B11", "B2", "D1")
DONOR_SHAPES <- setNames(c(16, 17, 15, 18, 8), DONORS) # as the approved demo (2026-11-11)
# Pre-registered labels (user, 2026-11-10, D3)
RESEMBLE_MEDIAN <- 0.8
RESEMBLE_MIN <- 0.5
CONTRIB_FRAC <- 0.1
CONTRIB_DONORS <- 4

cli_args <- commandArgs(trailingOnly = TRUE)
run_dir <- cli_args[[1]]
new_arms <- strsplit(cli_args[[2]], ",")[[1]]
stopifnot("run dir must not exist" = !dir.exists(run_dir))
arms <- union(BASE_ARMS, new_arms)
stopifnot(
  "common dir has its README" = file.exists(file.path(COMMON, "README.md"))
)
common_manifest <- read_csv(
  file.path(COMMON, "manifest.csv"),
  show_col_types = FALSE
)
missing_arms <- setdiff(arms, common_manifest$arm)
if (length(missing_arms) > 0) {
  stop(
    "arms missing from the common manifest: ",
    toString(missing_arms),
    call. = FALSE
  )
}
no_gs <- arms[!dir.exists(file.path(COMMON, "gs", arms))]
if (length(no_gs) > 0) {
  stop("arms without gs/ in the common dir: ", toString(no_gs), call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# Per donor: every arm's events on one arcsinh scale, and population rows ------------------------
gating_sets <- arms |>
  set_names() |>
  map(\(a) load_gs(file.path(COMMON, "gs", a)))
events_of <- function(gs, s) {
  x <- exprs(gh_pop_get_data(gs[[s]], "root"))[, CYTOKINES]
  colnames(x) <- names(CYTOKINES)
  linear <- !CYTOKINES %in% STORED_ARCSINH
  x[, linear] <- asinh(x[, linear] / COFACTOR)
  x
}
by_donor <- DONORS |>
  set_names() |>
  map(\(s) {
    ev <- map(gating_sets, \(gs) events_of(gs, s))
    stopifnot("same events in every arm" = n_distinct(map_int(ev, nrow)) == 1)
    # G8: Time, FSC-A, FSC-H and SSC-A identical per row in every arm
    pair <- map(gating_sets, \(gs) {
      exprs(gh_pop_get_data(gs[[s]], "root"))[, PAIR_COLS]
    })
    stopifnot(
      "rows pair across arms" = all(map_lgl(pair, \(m) identical(m, pair[[1]])))
    )
    in_pop <- map(gating_sets, \(gs) {
      map(set_names(POPS), \(pop) {
        gh_pop_get_indices(gs[[s]], paste0(POP_PATH, pop))
      })
    })
    stopifnot(
      "labels are common to every arm" = all(map_lgl(in_pop, \(p) {
        identical(p, in_pop[[1]])
      }))
    )
    list(ev = ev, in_pop = in_pop)
  })

# Readouts at the fixed cut: % positive, TNF+ among IFNg+ -------------------------------------
own <- by_donor |>
  imap(\(d, s) {
    arms |>
      set_names() |>
      map(\(a) {
        POPS |>
          set_names() |>
          map(\(pop) {
            x <- d$ev[[a]][d$in_pop[[a]][[pop]], c("TNF", "IFNg")]
            tibble(
              n_parent = nrow(x),
              pct_tnf = 100 * mean(x[, "TNF"] > CUT),
              pct_ifng = 100 * mean(x[, "IFNg"] > CUT),
              n_ifng_pos = sum(x[, "IFNg"] > CUT),
              pct_tnf_among_ifng = 100 * mean(x[x[, "IFNg"] > CUT, "TNF"] > CUT)
            )
          }) |>
          list_rbind(names_to = "pop")
      }) |>
      list_rbind(names_to = "arm")
  }) |>
  list_rbind(names_to = "donor")
write_csv(
  select(own, donor, arm, pop, n_parent, pct_tnf, pct_ifng),
  file.path(run_dir, "pct_pos.csv")
)
write_csv(
  select(own, donor, arm, pop, n_ifng_pos, pct_tnf_among_ifng),
  file.path(run_dir, "tnf_among_ifng.csv")
)

# 39b readouts on the common population rows --------------------------------------------------
nguyen <- \(v) unname(diff(quantile(v, c(0.5, 0.84)))) # q84 - q50
scores <- by_donor |>
  imap(\(d, s) {
    arms |>
      set_names() |>
      map(\(a) {
        POPS |>
          set_names() |>
          map(\(pop) {
            rows <- d$in_pop[[REF_ARM]][[pop]]
            c("TNF", "IFNg") |>
              set_names() |>
              map(\(cy) {
                lin <- COFACTOR * sinh(d$ev[[a]][rows, cy])
                lin_ref <- COFACTOR * sinh(d$ev[[REF_ARM]][rows, cy])
                ref_rsd <- nguyen(lin_ref[asinh(lin_ref / COFACTOR) < CUT])
                tibble(
                  rsd_neg = nguyen(lin[asinh(lin / COFACTOR) < CUT]),
                  ref_rsd = ref_rsd,
                  ref_rsd_cut_old = nguyen(lin_ref[
                    asinh(lin_ref / COFACTOR) < CUT_OLD
                  ]),
                  pos_share = mean(asinh(lin / COFACTOR) > CUT),
                  squeezed_share = mean(abs(lin) < SQUEEZE_FRAC * ref_rsd)
                )
              }) |>
              list_rbind(names_to = "cytokine")
          }) |>
          list_rbind(names_to = "pop")
      }) |>
      list_rbind(names_to = "arm")
  }) |>
  list_rbind(names_to = "donor") |>
  mutate(rsd_ratio = rsd_neg / ref_rsd)
write_csv(select(scores, -ref_rsd_cut_old), file.path(run_dir, "scores.csv"))
# the squeeze cut itself: SQUEEZE_FRAC x as188_joint's negative rSD, linear units
scores |>
  dplyr::filter(arm == REF_ARM) |>
  transmute(
    donor,
    pop,
    cytokine,
    squeeze_cut = SQUEEZE_FRAC * ref_rsd,
    squeeze_cut_at_cut_old = SQUEEZE_FRAC * ref_rsd_cut_old
  ) |>
  write_csv(file.path(run_dir, "squeeze_cut.csv"))

# Negative spread on shared negatives --------------------------------------------------------
# Shared negatives per arm: in the parent and below the cut in that arm AND in each base arm, so an
# arm's ratio does not depend on which other ablation arms are scored in the same call (47 used all
# arms of its figure). Nguyen rSD on the linear scale, ratio to as188_joint on the same cells.
neg_spread <- by_donor |>
  imap(\(d, s) {
    arms |>
      set_names() |>
      map(\(a) {
        set_arms <- union(BASE_ARMS, a)
        POPS |>
          set_names() |>
          map(\(pop) {
            in_all <- reduce(map(d$in_pop[set_arms], pop), `&`)
            names(CYTOKINES) |>
              set_names() |>
              map(\(cy) {
                neg_all <- reduce(map(d$ev[set_arms], \(x) x[, cy] < CUT), `&`)
                keep <- in_all & neg_all
                tibble(
                  n_shared = sum(keep),
                  rsd = nguyen(COFACTOR * sinh(d$ev[[a]][keep, cy])),
                  rsd_ref = nguyen(COFACTOR * sinh(d$ev[[REF_ARM]][keep, cy]))
                )
              }) |>
              list_rbind(names_to = "cytokine")
          }) |>
          list_rbind(names_to = "pop")
      }) |>
      list_rbind(names_to = "arm")
  }) |>
  list_rbind(names_to = "donor") |>
  mutate(rsd_ratio = rsd / rsd_ref)
write_csv(neg_spread, file.path(run_dir, "neg_spread.csv"))

# Gap closed and labels (D3): CD8T BV510, squeezed share and share above the cut -----------------
gap <- scores |>
  dplyr::filter(pop == "CD8T", cytokine == "TNF") |>
  select(donor, arm, squeezed_share, pos_share) |>
  pivot_longer(
    c(squeezed_share, pos_share),
    names_to = "readout",
    values_to = "value"
  ) |>
  mutate(
    start = value[arm == REF_ARM],
    target = value[arm == TARGET_ARM],
    frac_closed = (value - start) / (target - start),
    .by = c(donor, readout)
  )
write_csv(gap, file.path(run_dir, "gap_closed.csv"))
labels <- gap |>
  summarise(
    median_frac = median(frac_closed),
    min_frac = min(frac_closed),
    max_frac = max(frac_closed),
    n_contrib = sum(frac_closed >= CONTRIB_FRAC),
    resembles_one = median_frac >= RESEMBLE_MEDIAN & min_frac >= RESEMBLE_MIN,
    .by = c(arm, readout)
  ) |>
  mutate(
    label = case_when(
      all(resembles_one) ~ "resembles 1.0.2",
      all(n_contrib >= CONTRIB_DONORS) ~ "contributes",
      .default = "null"
    ),
    .by = arm
  )
write_csv(labels, file.path(run_dir, "labels.csv"))

# Figure: gap closed per arm, donors as points ------------------------------------------------
READOUT_LABELS <- c(
  squeezed_share = "CD8T TNF squeezed near 0",
  pos_share = sprintf("CD8T TNF above arcsinh %s", CUT)
)
fig <- gap |>
  dplyr::filter(!arm %in% c(REF_ARM, TARGET_ARM)) |>
  mutate(
    arm = factor(
      arm,
      levels = rev(unique(c(
        setdiff(BASE_ARMS, c(REF_ARM, TARGET_ARM)),
        new_arms
      )))
    ),
    readout = factor(READOUT_LABELS[readout], levels = unname(READOUT_LABELS))
  )
fig_median <- fig |>
  summarise(median_frac = median(frac_closed), .by = c(arm, readout))
p <- ggplot(fig, aes(frac_closed, arm)) +
  geom_vline(
    xintercept = c(0, 1),
    colour = "grey50",
    linewidth = FIG_LINE / 2
  ) +
  geom_vline(
    xintercept = c(CONTRIB_FRAC, RESEMBLE_MIN, RESEMBLE_MEDIAN),
    linetype = "dashed",
    colour = "grey70",
    linewidth = FIG_LINE / 2
  ) +
  geom_point(aes(colour = donor, shape = donor), size = 1.3) +
  geom_point(data = fig_median, aes(median_frac, arm), shape = 18, size = 2.5) +
  scale_colour_manual(values = categorical_palette(DONORS), name = "Donor") +
  scale_shape_manual(values = DONOR_SHAPES, name = "Donor") +
  facet_wrap(~readout, axes = "all", axis.labels = "margins") +
  labs(x = "Fraction of the as188_joint to s102_c188 gap closed", y = NULL) +
  value_grid("x")
save_panel(
  p,
  fig,
  run_dir,
  "fig_gap_closed",
  width_mm = 183,
  height_mm = 25 + 5 * nlevels(fig$arm)
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
