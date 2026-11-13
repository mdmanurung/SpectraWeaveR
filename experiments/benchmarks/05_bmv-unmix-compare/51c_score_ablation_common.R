# Purpose: score ablation arms on the common SpectroFlo gates (plans/common_gates_plan.md, batch 2). Copy of
#          51_score_ablation.R that only READS the common gated dir: it never gates (no 19a), never writes
#          gs/ or manifest, and stops when an arm is missing. Every arm carries the same SpectroFlo labels
#          by row (G8), so "own gate" and "as188_joint rows" are the same rows. Each call writes its own
#          run dir, so parallel calls cannot race.
#          Changes against 51 (all from the plan):
#            CUT 1 -> 1.5 (G7a): TNF+/IFNg+ cut, the negatives that define the rSD, and "above cut".
#            The squeeze cut (SQUEEZE_FRAC x as188_joint negative rSD) is refitted on the common pops;
#            squeeze_cut.csv gives it at cut 1.5 (used) and, for the old-vs-new comparison, at cut 1.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/51c_score_ablation_common.R <run_dir> <arms>
#            from the repo root. <arms> is a comma-separated list of arm names in the common gated dir;
#            every arm must also be in the config list arms_cg. <run_dir> must not exist and is
#            run_out("<key>") of the caller (3D example below). Example (Block A):
#            ... 51c_score_ablation_common.R experiments/runs/2026-11-11_bmv-common-ablation \
#                swap_spec102,swap_af102_var102,rev_var188,rev_af188,swap_thr102,regen188,bead_wide
# Inputs:  run_in("common-gated")/README.md, manifest.csv, gs/<arm>   [batch 0, 97c to 97d]
#          nodes /nonDebris/singlets/live/CD3pos/{CD4T,CD8T,DP,DN}, children TNFpos and IFNgpos (not used:
#          the cytokine cuts come from cut_of(): 3D the fixed 1.5, computed here)
#          datasets/<SW_DATASET>.json (arms_cg, channels, targets, pos_rule, neg_rule), <ds>_samples.csv
#          run_in("common-cuts")/cuts.csv when the dataset's rule is not fixed [97e]
# Outputs: <run_dir>/pct_pos.csv        TNF+ and IFNg+ % at cut 1.5 in CD4T and CD8T
#          <run_dir>/tnf_among_ifng.csv TNF+ % among IFNg+ cells
#          <run_dir>/scores.csv         squeezed share, share above 1.5, rSD of negatives; BV510 and BV570
#          <run_dir>/squeeze_cut.csv    squeeze cut per donor, pop, cytokine at cut 1.5 and at cut 1
#          <run_dir>/neg_spread.csv     Nguyen rSD ratio to as188_joint, TNF, IFNg, Th2, IL-10, shared negatives
#          <run_dir>/gap_closed.csv     per donor fraction of the as188_joint -> s102_c188 gap closed (24H: + condition)
#          <run_dir>/labels.csv         resembles / contributes / null per arm (24H: per arm and condition)
#          <run_dir>/fig_gap_closed.{pdf,png,csv} gap closed per arm, donors as points
#          <run_dir>/sessionInfo.txt
#          The column `donor` holds the sample id (3D and EXV the well, 24H donor_condition) in every table.

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
LIB <- "experiments/benchmarks/05_bmv-unmix-compare/lib"
source(file.path(LIB, "theme.R"))
source(file.path(LIB, "palettes.R"))

# Config lookups ---------------------------------------------------------------------------
# a nested config value that the dataset may lack: stop naming this script and the value
need <- \(group, name) {
  v <- cfg(group)[[name]]
  if (is.null(v)) {
    stop(
      "51c_score_ablation_common.R: needs ",
      name,
      ", missing in dataset ",
      DS$dataset,
      call. = FALSE
    )
  }
  v
}

# Parameters ------------------------------------------------------------------------------
COMMON <- run_in("common-gated") # batch 0 output, read only
ARMS_CG <- cfg("arms_cg") # arms of the dataset's common gated manifest (97c)
BASE_ARMS <- c("as188_joint", "s102_c188", "as102_var") # D1: start, chain target, final target
REF_ARM <- "as188_joint"
TARGET_ARM <- "s102_c188" # gap is as188_joint -> s102_c188 (D3)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
# cytokine cuts: cut_of(arm, sample, channel, "pos" or "neg"); 3D is the fixed 1.5 in every arm (G7a, G13)
CUT_OLD <- 1 # the pre-G7a cut, only to report the old squeeze cut on the new pops
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
TNF_LABEL <- cfg("targets")$tnf_label # name of the TNF-slot readout in the tables
IFNG_LABEL <- cfg("targets")$ifng_label
CYTOKINES <- set_names(
  c(
    need("channels", "TNF_slot"),
    need("channels", "IFNg_slot"),
    need("channels", "Th2_slot"),
    need("channels", "IL10_slot")
  ),
  c(TNF_LABEL, IFNG_LABEL, "Th2", "IL-10")
)
STORED_ARCSINH <- unname(CYTOKINES[c(TNF_LABEL, IFNG_LABEL)]) # 19a transforms these; the others are linear in the GatingSet
# label of the positive cut in the figure: the fixed value, or the per-donor rule
POS_LABEL <- if (identical(cfg("pos_rule")$type, "fixed")) {
  sprintf("arcsinh %s", cfg("pos_rule")$value)
} else {
  "the positive cut"
}
POPS <- c("CD4T", "CD8T")
POP_PATH <- "/nonDebris/singlets/live/CD3pos/" # full node path prefix, so CD4T/CD8T are unambiguous
PAIR_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # G8: rows must pair across arms
# Samples: 3D and EXV the well, 24H donor_condition. Figures plot the sample (3D, EXV) or the donor over conditions (24H).
SAMPLES <- samples_tbl()
SAMPLE_IDS <- sort(SAMPLES$sample)
COND_KEY <- if (is.null(DS$conditions)) character() else "condition"
PLOT_DONOR <- set_names(
  if (is.null(DS$conditions)) SAMPLES$sample else SAMPLES$donor,
  SAMPLES$sample
)
DONORS <- sort(unique(PLOT_DONOR))
DONOR_SHAPES <- setNames(rep_len(c(16, 17, 15, 18, 8), length(DONORS)), DONORS) # as the approved demo (2026-11-11)
add_condition <- \(d) {
  if (is.null(DS$conditions)) {
    d
  } else {
    left_join(
      d,
      select(SAMPLES, donor = sample, condition),
      by = "donor",
      relationship = "many-to-one",
      unmatched = c("error", "drop")
    )
  }
}
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
off_config <- setdiff(arms, ARMS_CG)
if (length(off_config) > 0) {
  stop(
    "51c_score_ablation_common.R: arms not in config arms_cg: ",
    toString(off_config),
    call. = FALSE
  )
}
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
# Cuts per arm, sample and cytokine, read once (cut_of re-reads cuts.csv when the rule is not fixed)
cuts_of_kind <- \(kind) {
  arms |>
    set_names() |>
    map(\(a) {
      SAMPLE_IDS |>
        set_names() |>
        map(\(s) map_dbl(CYTOKINES, \(ch) cut_of(a, s, ch, kind)))
    })
}
POS_CUTS <- cuts_of_kind("pos")
NEG_CUTS <- cuts_of_kind("neg")

dir.create(run_dir, recursive = TRUE)

# Per donor: every arm's events on one arcsinh scale, and population rows ------------------------
gating_sets <- arms |>
  set_names() |>
  map(\(a) load_gs(file.path(COMMON, "gs", a)))
stopifnot(
  "gating sets hold the samples of samples_tbl()" = all(map_lgl(
    gating_sets,
    \(g) setequal(sampleNames(g), SAMPLE_IDS)
  ))
)
events_of <- function(gs, s) {
  x <- exprs(gh_pop_get_data(gs[[s]], "root"))[, CYTOKINES]
  colnames(x) <- names(CYTOKINES)
  linear <- !CYTOKINES %in% STORED_ARCSINH
  x[, linear] <- asinh(x[, linear] / COFACTOR)
  x
}
by_donor <- SAMPLE_IDS |>
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
        pos_tnf <- POS_CUTS[[a]][[s]][[TNF_LABEL]]
        pos_ifng <- POS_CUTS[[a]][[s]][[IFNG_LABEL]]
        POPS |>
          set_names() |>
          map(\(pop) {
            x <- d$ev[[a]][d$in_pop[[a]][[pop]], c(TNF_LABEL, IFNG_LABEL)]
            tibble(
              n_parent = nrow(x),
              pct_tnf = 100 * mean(x[, TNF_LABEL] > pos_tnf),
              pct_ifng = 100 * mean(x[, IFNG_LABEL] > pos_ifng),
              n_ifng_pos = sum(x[, IFNG_LABEL] > pos_ifng),
              pct_tnf_among_ifng = 100 *
                mean(x[x[, IFNG_LABEL] > pos_ifng, TNF_LABEL] > pos_tnf)
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
            c(TNF_LABEL, IFNG_LABEL) |>
              set_names() |>
              map(\(cy) {
                neg_ref <- NEG_CUTS[[REF_ARM]][[s]][[cy]]
                neg_own <- NEG_CUTS[[a]][[s]][[cy]]
                pos_own <- POS_CUTS[[a]][[s]][[cy]]
                lin <- COFACTOR * sinh(d$ev[[a]][rows, cy])
                lin_ref <- COFACTOR * sinh(d$ev[[REF_ARM]][rows, cy])
                ref_rsd <- nguyen(lin_ref[asinh(lin_ref / COFACTOR) < neg_ref])
                tibble(
                  rsd_neg = nguyen(lin[asinh(lin / COFACTOR) < neg_own]),
                  ref_rsd = ref_rsd,
                  ref_rsd_cut_old = nguyen(lin_ref[
                    asinh(lin_ref / COFACTOR) < CUT_OLD
                  ]),
                  pos_share = mean(asinh(lin / COFACTOR) > pos_own),
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
                neg_all <- reduce(
                  imap(d$ev[set_arms], \(x, k) {
                    x[, cy] < NEG_CUTS[[k]][[s]][[cy]]
                  }),
                  `&`
                )
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
  dplyr::filter(pop == "CD8T", cytokine == TNF_LABEL) |>
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
  ) |>
  add_condition()
write_csv(gap, file.path(run_dir, "gap_closed.csv"))
labels <- gap |>
  summarise(
    median_frac = median(frac_closed),
    min_frac = min(frac_closed),
    max_frac = max(frac_closed),
    n_contrib = sum(frac_closed >= CONTRIB_FRAC),
    resembles_one = median_frac >= RESEMBLE_MEDIAN & min_frac >= RESEMBLE_MIN,
    .by = c(arm, readout, any_of(COND_KEY))
  ) |>
  mutate(
    label = case_when(
      all(resembles_one) ~ "resembles 1.0.2",
      all(n_contrib >= CONTRIB_DONORS) ~ "contributes",
      .default = "null"
    ),
    .by = c(arm, any_of(COND_KEY))
  )
write_csv(labels, file.path(run_dir, "labels.csv"))

# Figure: gap closed per arm, donors as points ------------------------------------------------
READOUT_LABELS <- c(
  squeezed_share = sprintf("CD8T %s squeezed near 0", TNF_LABEL),
  pos_share = sprintf("CD8T %s above %s", TNF_LABEL, POS_LABEL)
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
    readout = factor(READOUT_LABELS[readout], levels = unname(READOUT_LABELS)),
    donor = unname(PLOT_DONOR[donor])
  )
fig_median <- fig |>
  summarise(
    median_frac = median(frac_closed),
    .by = c(arm, readout, any_of(COND_KEY))
  )
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
  facet_wrap(c("readout", COND_KEY), axes = "all", axis.labels = "margins") +
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
