#!/usr/bin/env Rscript
# Purpose: check whether the unmixing arms also differ for IL-2, IL-10 and Th2, as they do for TNF, or
#          agree, as they do for IFNg, on the common SpectroFlo gates (G1a, G7a). Per arm, donor and parent
#          (CD4T, CD8T): % positive at a fixed cut, density shape near 0, co-expression with TNF and IFNg,
#          Spearman rho for every cytokine pair and for three expected links (IFNg-Tbet, IL-10-FoxP3,
#          IL-2-CD25), and negative spread against as188_joint on cells negative in every arm.
#          Copy of 47_other_cytokines.R: new gated dir, cut 1.5, new run dir, Nature-style figures.
# Run:     SW_DATASET=<3d|24h> SW_RUN_DATE=<YYYY-MM-DD> \
#            pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/47c_other_cytokines_common.R [<run_dir>]
#          from the repo root, after batch 0 has written GATED_DIR/README.md. <run_dir> overrides the registry
#          output run_out("common-cut15-47c"); it must not exist. Stops on a dataset without TNF, IFNg,
#          IL-2, IL-10 and Th2 (EXV).
# Inputs:  run_in("common-gated")/manifest.csv, gs/<arm>/  common-gate GatingSets [97_common_gs.R]
#          cut_of(arm, sample, channel, "pos" and "neg")  positive and negative cuts [config pos_rule and
#            neg_rule; cuts.csv from 97e when not fixed]
#          experiments/benchmarks/05_bmv-unmix-compare/datasets/  config, run registry, sample table
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R, report_figures.R  labels, builders
# Outputs: <run_dir>/pct_pos.csv                % above the cut per arm, sample, parent, cytokine
#          <run_dir>/il2_coexpression.csv       % IL-2+ among TNF+ and among IFNg+
#          <run_dir>/rho_pairs.csv              Spearman rho per arm, donor, parent, pair
#          <run_dir>/neg_spread.csv             Nguyen rSD on shared negatives and ratio to as188_joint
#          <run_dir>/fig_a_pct_pos.{pdf,png,csv}            % positive, donors as points, median bar
#          <run_dir>/fig_b_density_<cytokine>.{pdf,png,csv} density near 0 per donor, parent and arm
#          <run_dir>/fig_c_il2_coexpression.{pdf,png,csv}   % IL-2+ among TNF+ and IFNg+
#          <run_dir>/fig_d_rho_cytokine_pairs.{pdf,png,csv} rho for the 10 cytokine pairs
#          <run_dir>/fig_e_rho_expected_links.{pdf,png,csv} rho for the 3 expected links
#          <run_dir>/fig_f_neg_spread.{pdf,png,csv}         rSD ratio to as188_joint
#          <run_dir>/fig_{a,b,c,d,e,f}_table.csv            tables handed to the builders
#          <run_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), run_in(), run_out(), cut_of()
source("experiments/benchmarks/05_bmv-unmix-compare/report_figures.R") # theme, save_panel, builders, ARM_LABELS

# Parameters ------------------------------------------------------------------------------
GATED_DIR <- run_in("common-gated") # batch 0 output (G9: read-only here)
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
# This script needs TNF, IFNg, IL-2, IL-10 and Th2; stop and name what the dataset lacks
EXTRA_CHANNEL_KEYS <- c(`IL-2` = "IL2", `IL-10` = "IL10", Th2 = "Th2") # names inside cfg("channels")
missing_targets <- c(
  if (!identical(cfg("targets")$tnf_label, "TNF")) "TNF",
  if (!identical(cfg("targets")$ifng_label, "IFNg")) "IFNg",
  names(EXTRA_CHANNEL_KEYS)[!EXTRA_CHANNEL_KEYS %in% names(cfg("channels"))]
)
if (length(missing_targets) > 0) {
  stop(
    "47c_other_cytokines_common.R: needs ",
    paste(missing_targets, collapse = ", "),
    ", missing in dataset ",
    DS$dataset,
    call. = FALSE
  )
}
# TNF and IFNg are stored as arcsinh in the GatingSets (19a); every other channel is stored linear.
# The cut comes from cut_of() per arm, sample and channel (3D: fixed arcsinh 1.5, G7a).
CYTOKINES <- c(
  TNF = cfg("channels")$TNF_slot,
  IFNg = cfg("channels")$IFNg_slot,
  unlist(cfg("channels")[EXTRA_CHANNEL_KEYS]) |>
    set_names(names(EXTRA_CHANNEL_KEYS))
)
STORED_ARCSINH <- c(cfg("channels")$TNF_slot, cfg("channels")$IFNg_slot)
# Cut text for titles: the fixed value when the rule is fixed (3D: 1.5), else the per-sample cut
CUT_VALUE_TXT <- if (identical(cfg("pos_rule")$type, "fixed")) {
  cfg("pos_rule")$value
} else {
  "per sample"
}
CUT_TXT <- if (identical(cfg("pos_rule")$type, "fixed")) {
  paste("arcsinh", CUT_VALUE_TXT)
} else {
  "the per-sample cut"
}
NEW_CYTOKINES <- c("IL-2", "IL-10", "Th2") # the three this script adds
# expected biological links, used as reference: rank correlation only, no cut
LINKS <- tribble(
  ~link           , ~from   , ~to     , ~parents          ,
  "IFNg - Tbet"   , "IFNg"  , "Tbet"  , c("CD4T", "CD8T") ,
  "IL-10 - FoxP3" , "IL-10" , "FoxP3" , "CD4T"            ,
  "IL-2 - CD25"   , "IL-2"  , "CD25"  , c("CD4T", "CD8T")
)
MARKERS <- c(Tbet = "PE-Cy5-A", FoxP3 = "RB780-A", CD25 = "SBUV605-A")
NODES <- c(
  CD4T = "/nonDebris/singlets/live/CD3pos/CD4T",
  CD8T = "/nonDebris/singlets/live/CD3pos/CD8T"
)
POPULATIONS <- names(NODES)
REF_ARM <- "as188_joint" # reference arm for the spread ratio (plan 2026-11-09)
ARM_ORDER <- unlist(VERSION_ROWS, use.names = FALSE)
VERSION_GROUP <- VERSION_ROWS |>
  imap(\(arms, group) set_names(rep(group, length(arms)), arms)) |>
  unname() |>
  unlist()
DENSITY_RANGE <- c(-1.5, 4) # arcsinh display range for densities; display only

run_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(run_dir)) {
  run_dir <- run_out("common-cut15-47c")
}
stopifnot("47c run dir must not exist" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)
manifest <- read_csv(
  file.path(GATED_DIR, "manifest.csv"),
  show_col_types = FALSE
)
stopifnot(
  "all 11 figure arms are in the gated manifest" = all(
    ARM_ORDER %in% manifest$arm
  ),
  "gated samples equal the sample table" = setequal(
    manifest$sample,
    samples_tbl()$sample
  )
)
donors <- sort(unique(manifest$sample))
gating_sets <- ARM_ORDER |>
  set_names() |>
  map(\(a) load_gs(file.path(GATED_DIR, "gs", a)))

# Per donor: every arm's events on one arcsinh scale -------------------------------------------
# All arms come from the same raw events, so row i is the same cell in every arm. Linear channels get
# arcsinh(x / 6000) in memory; TNF and IFNg are already arcsinh.
channels <- c(CYTOKINES, MARKERS)
events_of <- function(gs, s) {
  x <- exprs(gh_pop_get_data(gs[[s]], "root"))[, channels]
  colnames(x) <- names(channels)
  linear <- !channels %in% STORED_ARCSINH
  x[, linear] <- asinh(x[, linear] / COFACTOR)
  x
}
by_donor <- donors |>
  set_names() |>
  map(\(s) {
    ev <- map(gating_sets, \(gs) events_of(gs, s))
    stopifnot("same events in every arm" = n_distinct(map_int(ev, nrow)) == 1)
    in_pop <- map(gating_sets, \(gs) {
      map(set_names(POPULATIONS), \(pop) {
        gh_pop_get_indices(gs[[s]], NODES[[pop]])
      })
    })
    list(ev = ev, in_pop = in_pop)
  })

# Per arm, donor, parent: % positive, IL-2 co-expression, rank correlations ---------------------
pair_names <- combn(names(CYTOKINES), 2, simplify = FALSE)
per_cell_metrics <- by_donor |>
  imap(\(d, s) {
    ARM_ORDER |>
      set_names() |>
      map(\(a) {
        POPULATIONS |>
          set_names() |>
          map(\(pop) {
            x <- d$ev[[a]][d$in_pop[[a]][[pop]], ]
            cut_pos <- map_dbl(CYTOKINES, \(ch) cut_of(a, s, ch, "pos"))
            pos <- sweep(x[, names(CYTOKINES)], 2, cut_pos, ">")
            pairs <- pair_names |>
              map(\(p) {
                tibble(
                  pair = paste(p, collapse = " - "),
                  rho = cor(x[, p[[1]]], x[, p[[2]]], method = "spearman")
                )
              }) |>
              list_rbind()
            links <- LINKS |>
              dplyr::filter(map_lgl(parents, \(ps) pop %in% ps)) |>
              mutate(
                rho = map2_dbl(from, to, \(cx, cy) {
                  cor(x[, cx], x[, cy], method = "spearman")
                })
              ) |>
              select(pair = link, rho)
            list(
              pct = enframe(
                100 * colMeans(pos),
                name = "cytokine",
                value = "pct_pos"
              ),
              il2 = tibble(
                among = c("TNF+", "IFNg+"),
                n_among = c(sum(pos[, "TNF"]), sum(pos[, "IFNg"])),
                pct_il2 = 100 *
                  c(
                    mean(pos[pos[, "TNF"], "IL-2"]),
                    mean(pos[pos[, "IFNg"], "IL-2"])
                  )
              ),
              rho = bind_rows(
                mutate(pairs, kind = "cytokine pair"),
                mutate(links, kind = "expected link")
              )
            )
          })
      })
  })

# flatten donor -> arm -> parent into long tables
flatten_metric <- function(metric) {
  per_cell_metrics |>
    map(\(by_a) {
      map(by_a, \(by_p) list_rbind(map(by_p, metric), names_to = "parent")) |>
        list_rbind(names_to = "arm")
    }) |>
    list_rbind(names_to = "donor") |>
    mutate(arm = factor(arm, levels = ARM_ORDER))
}
pct_pos <- flatten_metric("pct")
il2_coexpression <- flatten_metric("il2")
rho_pairs <- flatten_metric("rho")
write_csv(pct_pos, file.path(run_dir, "pct_pos.csv"))
write_csv(il2_coexpression, file.path(run_dir, "il2_coexpression.csv"))
write_csv(rho_pairs, file.path(run_dir, "rho_pairs.csv"))

# Negative spread on cells negative in every arm ------------------------------------------------
# Shared negatives for a cytokine: in the parent in every arm and below the cut in every arm. Nguyen
# robust SD (84th - 50th percentile) on the linear scale, as in 34_rank_symmetric.R.
neg_spread <- by_donor |>
  imap(\(d, s) {
    POPULATIONS |>
      set_names() |>
      map(\(pop) {
        in_all <- reduce(map(d$in_pop, pop), `&`)
        names(CYTOKINES) |>
          set_names() |>
          map(\(cy) {
            neg_all <- reduce(
              imap(d$ev, \(x, a) {
                x[, cy] < cut_of(a, s, CYTOKINES[[cy]], "neg")
              }),
              `&`
            )
            keep <- in_all & neg_all
            map(d$ev, \(x) {
              lin <- COFACTOR * sinh(x[keep, cy])
              tibble(
                n_shared = sum(keep),
                rsd = unname(diff(quantile(lin, c(0.5, 0.84))))
              )
            }) |>
              list_rbind(names_to = "arm")
          }) |>
          list_rbind(names_to = "cytokine")
      }) |>
      list_rbind(names_to = "parent")
  }) |>
  list_rbind(names_to = "donor") |>
  mutate(
    rsd_ratio = rsd / rsd[arm == REF_ARM],
    .by = c(donor, parent, cytokine)
  ) |>
  mutate(arm = factor(arm, levels = ARM_ORDER))
write_csv(neg_spread, file.path(run_dir, "neg_spread.csv"))

# Figure a: % positive at the cut, new cytokines ----------------------------------------------
fig_a <- pct_pos |>
  dplyr::filter(cytokine %in% NEW_CYTOKINES) |>
  mutate(cytokine = factor(cytokine, levels = NEW_CYTOKINES))
write_csv(fig_a, file.path(run_dir, "fig_a_table.csv"))
redraw_panel(
  file.path(run_dir, "fig_a_table.csv"),
  "pct_pos",
  paste0("% positive at ", CUT_TXT),
  run_dir,
  "fig_a_pct_pos"
)

# Figure b: density near 0, one file per new cytokine -----------------------------------------
# lines per donor, panels per arm and parent; same events in every arm, so shape differences come
# from unmixing
fig_b <- by_donor |>
  imap(\(d, s) {
    ARM_ORDER |>
      set_names() |>
      map(\(a) {
        POPULATIONS |>
          set_names() |>
          map(\(pop) {
            x <- d$ev[[a]][d$in_pop[[a]][[pop]], NEW_CYTOKINES]
            NEW_CYTOKINES |>
              set_names() |>
              map(\(cy) {
                dens <- density(
                  x[, cy],
                  from = DENSITY_RANGE[[1]],
                  to = DENSITY_RANGE[[2]],
                  n = 256
                )
                tibble(value = dens$x, density = dens$y)
              }) |>
              list_rbind(names_to = "cytokine")
          }) |>
          list_rbind(names_to = "parent")
      }) |>
      list_rbind(names_to = "arm")
  }) |>
  list_rbind(names_to = "donor") |>
  mutate(version = VERSION_GROUP[arm])
write_csv(fig_b, file.path(run_dir, "fig_b_table.csv"))
for (cy in NEW_CYTOKINES) {
  redraw_panel(
    file.path(run_dir, "fig_b_table.csv"),
    "density_curves",
    paste0("fig_b_density_", cy),
    run_dir,
    paste0("fig_b_density_", cy)
  )
}

# Figure c: % IL-2+ among TNF+ and among IFNg+ ------------------------------------------------
write_csv(il2_coexpression, file.path(run_dir, "fig_c_table.csv"))
redraw_panel(
  file.path(run_dir, "fig_c_table.csv"),
  "il2",
  paste0("% IL-2+ among cytokine-positive cells, cut ", CUT_VALUE_TXT),
  run_dir,
  "fig_c_il2_coexpression"
)

# Figures d and e: Spearman rho, donors as points over a median bar ---------------------------
write_csv(
  dplyr::filter(rho_pairs, kind == "cytokine pair"),
  file.path(run_dir, "fig_d_table.csv")
)
write_csv(
  dplyr::filter(rho_pairs, kind == "expected link"),
  file.path(run_dir, "fig_e_table.csv")
)
redraw_panel(
  file.path(run_dir, "fig_d_table.csv"),
  "rho_pairs",
  "Rho for the 10 cytokine pairs",
  run_dir,
  "fig_d_rho_cytokine_pairs"
)
redraw_panel(
  file.path(run_dir, "fig_e_table.csv"),
  "rho_pairs",
  "Rho for the 3 expected links",
  run_dir,
  "fig_e_rho_expected_links"
)

# Figure f: negative spread ratio to as188_joint ----------------------------------------------
write_csv(neg_spread, file.path(run_dir, "fig_f_table.csv"))
redraw_panel(
  file.path(run_dir, "fig_f_table.csv"),
  "neg_spread",
  paste0("Negative rSD ratio to ", REF_ARM),
  run_dir,
  "fig_f_neg_spread"
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
