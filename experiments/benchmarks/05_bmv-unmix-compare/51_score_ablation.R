# Purpose: score ablation arms with full gating (user, 2026-11-10, README D4) and label each by the
#          pre-registered gap criterion (D3). For each new arm: mask its raw unmixed files with the shared
#          cleanup rows of as188_joint, gate it with 19a_gate_arm.R into one shared gated dir, then score it
#          next to the three base arms (as188_joint, s102_c188, as102_var), whose GatingSets are copied
#          from 2026-11-10_bmv-gated-pooled.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/51_score_ablation.R <run_dir> <arms>
#            from the repo root. <arms> is a comma-separated list of arm folders under
#            experiments/bmv_test/compare_3D/data/unmix_cells/ (each with 5 raw unmixed files named as
#            as188_joint's). <run_dir> must not exist; one run dir per call. Arms already gated in the
#            shared gated dir are reused, not regated. Example (Block A):
#            ... 51_score_ablation.R experiments/runs/2026-11-10_bmv-ablation \
#                swap_spec102,swap_af102_var102,rev_var188,rev_af188,swap_thr102,regen188,bead_wide
# Inputs:  experiments/runs/2026-11-10_bmv-gated-pooled/manifest.csv, gs/<base arm>  [19a_gate_arm.R]
#          experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/AutoSpectral_unmixed  raw as188_joint
#          experiments/bmv_test/compare_3D/data/unmix_cells/<arm>/*.fcs  [48, 48b, 49 scripts]
# Outputs: experiments/bmv_test/compare_3D/data/cleaned/2026-11-10_bmv-ablation/<arm>/*.fcs  masked files
#          experiments/runs/2026-11-10_bmv-ablation-gated/  shared gated dir: manifest.csv, gs/<arm>,
#            gate_counts/, templates/ [19a]; manifest_rows/<arm>.csv
#          <run_dir>/pct_pos.csv        TNF+ and IFNg+ % at cut 1 in each arm's own CD4T and CD8T
#          <run_dir>/tnf_among_ifng.csv TNF+ % among IFNg+ cells, own gates
#          <run_dir>/scores.csv         squeezed share, share above 1, rSD of negatives; BV510 and BV570 on
#                                       as188_joint's CD4T and CD8T rows (39b method)
#          <run_dir>/neg_spread.csv     Nguyen rSD ratio to as188_joint, TNF, IFNg, Th2, IL-10, shared negatives
#          <run_dir>/gap_closed.csv     per donor fraction of the as188_joint -> s102_c188 gap closed
#          <run_dir>/labels.csv         resembles / contributes / null per arm
#          <run_dir>/fig_gap_closed.png gap closed per arm, donors as points (+ fig_gap_closed_source.csv)
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

# Parameters ------------------------------------------------------------------------------
DATA <- "experiments/bmv_test/compare_3D/data"
POOLED <- "experiments/runs/2026-11-10_bmv-gated-pooled" # base arms' cleaned files and GatingSets
GATED <- "experiments/runs/2026-11-10_bmv-ablation-gated" # shared gated dir for every ablation call
CLEANED <- file.path(DATA, "cleaned", "2026-11-10_bmv-ablation")
BASE_ARMS <- c("as188_joint", "s102_c188", "as102_var") # D1: start, chain target, final target
REF_ARM <- "as188_joint"
TARGET_ARM <- "s102_c188" # gap is as188_joint -> s102_c188 (D3)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
CUT <- 1 # fixed cytokine cut, arcsinh, every arm and donor (user, 2026-11-10)
SQUEEZE_FRAC <- 0.5 # squeezed = |x| below this fraction of as188_joint's rSD (39b; user, 2026-11-10)
CYTOKINES <- c(
  TNF = "BV510-A",
  IFNg = "BV570-A",
  Th2 = "eFluor 450-A",
  `IL-10` = "BB700-A"
)
STORED_ARCSINH <- c("BV510-A", "BV570-A") # 19a transforms these; the others are linear in the GatingSet
POPS <- c("CD4T", "CD8T")
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # identical across arms; maps raw rows to cleaned rows
VIABILITY <- "Zombie NIR-A" # every AutoSpectral arm (manifest of 2026-11-10_bmv-gated-pooled)
# Pre-registered labels (user, 2026-11-10, D3)
RESEMBLE_MEDIAN <- 0.8
RESEMBLE_MIN <- 0.5
CONTRIB_FRAC <- 0.1
CONTRIB_DONORS <- 4

cli_args <- commandArgs(trailingOnly = TRUE)
run_dir <- cli_args[[1]]
new_arms <- strsplit(cli_args[[2]], ",")[[1]]
stopifnot("run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)
walk(file.path(GATED, c("gs", "manifest_rows")), \(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
})
arms <- c(BASE_ARMS, new_arms)

# Base arms: copy GatingSets and manifest rows from the pooled run (once) ------------------------
pooled_manifest <- read_csv(
  file.path(POOLED, "manifest.csv"),
  show_col_types = FALSE
)
walk(BASE_ARMS, \(a) {
  if (!dir.exists(file.path(GATED, "gs", a))) {
    file.copy(
      file.path(POOLED, "gs", a),
      file.path(GATED, "gs"),
      recursive = TRUE
    )
    pooled_manifest |>
      dplyr::filter(arm == a) |>
      write_csv(file.path(GATED, "manifest_rows", paste0(a, ".csv")))
  }
})

# Shared cleanup rows: as188_joint raw rows kept in its cleaned file ------------------------------
read_x <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
}
key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))
ref_rows <- pooled_manifest |>
  dplyr::filter(arm == REF_ARM)
raw_ref_dir <- file.path(DATA, "autoSpectral_v1.8.8/AutoSpectral_unmixed")

# New arms: mask raw files, write manifest rows, gate with 19a -----------------------------------
rscript <- file.path(R.home("bin"), "Rscript")
for (a in setdiff(new_arms, BASE_ARMS)) {
  if (dir.exists(file.path(GATED, "gs", a))) {
    next
  }
  arm_clean_dir <- file.path(CLEANED, a)
  stopifnot(
    "cleaned dir for a new arm must not exist" = !dir.exists(arm_clean_dir)
  )
  dir.create(arm_clean_dir, recursive = TRUE)
  rows <- ref_rows |>
    mutate(
      arm = a,
      dir = file.path("unmix_cells", a),
      file = file.path(arm_clean_dir, basename(file))
    )
  pwalk(list(ref_rows$file, rows$file), \(ref_clean, out_file) {
    raw_ref <- read_x(file.path(raw_ref_dir, basename(ref_clean)))
    keep <- key_of(raw_ref) %in% key_of(read_x(ref_clean))
    ff <- read.FCS(
      file.path(DATA, "unmix_cells", a, basename(ref_clean)),
      transformation = FALSE,
      truncate_max_range = FALSE
    )
    stopifnot("arm rows equal as188_joint raw rows" = nrow(ff) == length(keep))
    write.FCS(ff[keep, ], out_file)
  })
  write_csv(rows, file.path(GATED, "manifest_rows", paste0(a, ".csv")))
  # 19a reads <gated>/manifest.csv; rebuild it from every arm's rows before each call
  list.files(file.path(GATED, "manifest_rows"), full.names = TRUE) |>
    map(\(f) read_csv(f, show_col_types = FALSE)) |>
    list_rbind() |>
    write_csv(file.path(GATED, "manifest.csv"))
  status <- system2(
    rscript,
    c("experiments/benchmarks/05_bmv-unmix-compare/19a_gate_arm.R", GATED, a)
  )
  if (status != 0) {
    stop("19a_gate_arm.R failed for ", a, " (exit ", status, ")", call. = FALSE)
  }
}

# Per donor: every arm's events on one arcsinh scale, and population rows ------------------------
gating_sets <- arms |>
  set_names() |>
  map(\(a) load_gs(file.path(GATED, "gs", a)))
donors <- sort(ref_rows$sample)
events_of <- function(gs, s) {
  x <- exprs(gh_pop_get_data(gs[[s]], "root"))[, CYTOKINES]
  colnames(x) <- names(CYTOKINES)
  linear <- !CYTOKINES %in% STORED_ARCSINH
  x[, linear] <- asinh(x[, linear] / COFACTOR)
  x
}
by_donor <- donors |>
  set_names() |>
  map(\(s) {
    ev <- map(gating_sets, \(gs) events_of(gs, s))
    stopifnot("same events in every arm" = n_distinct(map_int(ev, nrow)) == 1)
    in_pop <- map(gating_sets, \(gs) {
      map(set_names(POPS), \(pop) gh_pop_get_indices(gs[[s]], pop))
    })
    list(ev = ev, in_pop = in_pop)
  })

# Own-gate readouts: % positive at the cut, TNF+ among IFNg+ ------------------------------------
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

# 39b readouts on as188_joint's population rows (rows paired across arms) ------------------------
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
write_csv(scores, file.path(run_dir, "scores.csv"))

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

# Gap closed and labels (D3): CD8T BV510, squeezed share and share above 1 ------------------------
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

# Figure: gap closed per arm, donors as points --------------------------------------------------
fig <- gap |>
  dplyr::filter(!arm %in% c(REF_ARM, TARGET_ARM)) |>
  mutate(
    arm = factor(
      arm,
      levels = rev(c(setdiff(BASE_ARMS, c(REF_ARM, TARGET_ARM)), new_arms))
    ),
    readout = factor(
      readout,
      levels = c("squeezed_share", "pos_share"),
      labels = c("CD8T TNF squeezed near 0", "CD8T TNF above arcsinh 1")
    )
  )
write_csv(fig, file.path(run_dir, "fig_gap_closed_source.csv"))
fig_median <- fig |>
  summarise(median_frac = median(frac_closed), .by = c(arm, readout))
p <- ggplot(fig, aes(frac_closed, arm)) +
  geom_vline(xintercept = c(0, 1), colour = "grey50") +
  geom_vline(
    xintercept = c(CONTRIB_FRAC, RESEMBLE_MIN, RESEMBLE_MEDIAN),
    linetype = "dashed",
    colour = "grey70"
  ) +
  geom_point(aes(colour = donor), size = 2) +
  geom_point(data = fig_median, aes(median_frac, arm), shape = 18, size = 4) +
  scale_colour_brewer(palette = "Dark2", name = "Donor") +
  facet_wrap(~readout) +
  labs(x = "Fraction of the as188_joint to s102_c188 gap closed", y = NULL) +
  theme_bw(base_size = 11)
ggsave(
  file.path(run_dir, "fig_gap_closed.png"),
  p,
  width = 200,
  height = 30 + 9 * nlevels(fig$arm),
  units = "mm",
  dpi = 200,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
