# Purpose: common gates, pooled version (plans/common_gates_plan.md step 1; user 2026-11-11: "make pooled gates
#          on all samples; do not gate per donor"). Gate SpectroFlo fresh with the 19a template and methods,
#          but with collapseDataForGating = TRUE and groupBy = 5 on EVERY row, so debris, singlets, live,
#          CD3+ and the CD4/CD8 quadrants each get one gate fitted on all 5 donors pooled and applied to
#          each donor. Exports per-event labels, checks that every GatingSet arm pairs row by row with
#          SpectroFlo (Time, FSC-A, FSC-H, SSC-A), compares the old per-donor gates to the pooled ones,
#          and draws the gate figures per donor. 97_spectroflo_labels.R (per-donor gates) is left as it was.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/97b_spectroflo_labels_pooled.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-10_bmv-gated-pooled/manifest.csv       FCS paths [01_manifest-gate.R]
#          experiments/runs/2026-11-10_bmv-gated-pooled/gs/<arm>/          GatingSets, pairing and old gates [19a]
#          experiments/runs/2026-11-10_bmv-ablation-gated/gs/<arm>/        GatingSets, pairing
#          experiments/runs/2026-11-11_bmv-common-gates/label_counts.csv   old per-donor % per gate [97]
#          experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R
# Outputs: <out_dir>/gs/spectroflo/, templates/spectroflo.csv   pooled GatingSet and the template used
#          <out_dir>/labels.csv.gz     one row per event and donor: sample, row (clean2), raw_row, FSC_A, SSC_A,
#                                      nondebris, singlet, live, CD3pos, CD4T, CD8T, CD3neg_live
#          <out_dir>/label_counts.csv  events per gate and donor, and % of parent
#          <out_dir>/quadrant_counts.csv  CD4T, CD8T, DP, DN events and % of CD3+ per donor
#          <out_dir>/cuts.csv          gate bounds per version (old per-donor, pooled), donor and node
#          <out_dir>/compare_old_pooled.csv  n and % of parent per gate and donor, old vs pooled
#          <out_dir>/pairing.csv       arm, sample, n_rows_equal, time_equal, fsca_equal, fsch_equal, ssca_equal
#          <out_dir>/figures/strategy_<donor>.png, .pdf, .csv       debris, singlets, live, CD3+ gates
#          <out_dir>/figures/gate_cd4_cd8_<donor>.png, .pdf, .csv   CD4 vs CD8, four quadrants shaded
#          <out_dir>/sessionInfo.txt

.libPaths(c(
  "experiments/bmv_test/compare_3D/data/lib_ggplot2_35", # ggplot2 3.5.2 for ggcyto
  "experiments/bmv_test/compare_3D/data/lib_cytoutils", # cytoUtils 0.1.1 as 19a
  .libPaths()
))
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(cytoUtils)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tibble)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE) # sw_gate_run
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters ---------------------------------------------------------------------------------------
GATED_POOLED <- "experiments/runs/2026-11-10_bmv-gated-pooled" # 19a gating of every arm
OUT_DIR <- "experiments/runs/2026-11-11_bmv-common-gates-pooled"
OLD_DIR <- "experiments/runs/2026-11-11_bmv-common-gates" # per-donor gating from 97
GATED_ABLATION <- "experiments/runs/2026-11-10_bmv-ablation-gated"
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
PAIR_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # all must be equal on every clean2 row
VIABILITY <- "LIVE DEAD NIR-A" # SpectroFlo viability channel (manifest.csv)
DEBRIS_RANGE <- "c(2e5, 1.5e6)" # as 19a_gate_arm.R
LINEAGE_RANGE <- "c(0, 3)" # as 19a_gate_arm.R
CD4_CD8_TMIX_ARGS <- "K=3, usePrior='no', quantile1=0.8, quantile3=0.8" # as 19a_gate_arm.R
TNF <- "BV510-A" # transformed with the others, as 19a_gate_arm.R
IFNG <- "BV570-A"
RAW_DIR <- "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Raw" # raw files, as 88_spike_in_markers.R
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # as 88_spike_in_markers.R:51; maps clean rows to raw rows
REF_ARM <- "spectroflo" # the arm whose gates all arms inherit (G1)
BINS <- 60 # hexagons per panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
CD4 <- "Spark Blue 574-A" # channels as in 19a_gate_arm.R
CD8 <- "SBUV445-A"
CD3 <- "SBUV510-A"
QUAD_FILL <- c(
  CD4T = "#2166ac",
  CD8T = "#b2182b",
  DP = "#762a83",
  DN = "#7f7f7f"
) # one colour per quadrant
QUADRANTS <- c(
  CD4T = paste0(CD4, "+", CD8, "-"),
  CD8T = paste0(CD4, "-", CD8, "+"),
  DP = paste0(CD4, "+", CD8, "+"),
  DN = paste0(CD4, "-", CD8, "-")
) # openCyto quadrant node names -> aliases, as 19a
AXIS_T <- c(-3, 6) # CD4 / CD8 display window, arcsinh; the % use every CD3+ cell
# gate node -> parent node, in hierarchy order; CD3neg_live is live & !CD3pos and is not a gate node
NODES <- c(
  nondebris = "/nonDebris",
  singlet = "/nonDebris/singlets",
  live = "/nonDebris/singlets/live",
  CD3pos = "/nonDebris/singlets/live/CD3pos",
  CD4T = "/nonDebris/singlets/live/CD3pos/CD4T",
  CD8T = "/nonDebris/singlets/live/CD3pos/CD8T"
)
PARENT <- c(
  nondebris = "all events",
  singlet = "nondebris",
  live = "singlet",
  CD3pos = "live",
  CD4T = "CD3pos",
  CD8T = "CD3pos",
  CD3neg_live = "live"
)

stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Gate SpectroFlo fresh, every row pooled over the 5 donors ------------------------------------------------
manifest <- read_csv(
  file.path(GATED_POOLED, "manifest.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(arm == REF_ARM) |>
  arrange(sample)
donors <- manifest$sample
stopifnot("5 donors" = length(donors) == 5)
# same methods and arguments as 19a_gate_arm.R; collapseDataForGating = TRUE and groupBy = 5 on every row
template <- tibble(
  alias = c("nonDebris", "singlets", "live", "CD3pos", "*"),
  pop = c("+", "+", "-", "+", "+/-+/-"),
  parent = c("root", "nonDebris", "singlets", "live", "CD3pos"),
  dims = c("FSC-A", "FSC-A,FSC-H", VIABILITY, CD3, paste(CD4, CD8, sep = ",")),
  gating_method = c(
    "mindensity",
    "singletGate",
    "mindensity",
    "mindensity",
    "quadGate.tmix"
  ),
  gating_args = c(
    paste0("gate_range=", DEBRIS_RANGE),
    "wider_gate=TRUE",
    NA,
    paste0("gate_range=", LINEAGE_RANGE),
    CD4_CD8_TMIX_ARGS
  ),
  collapseDataForGating = TRUE,
  groupBy = nrow(manifest),
  preprocessing_method = NA,
  preprocessing_args = NA
)
dir.create(file.path(OUT_DIR, "templates"))
template_file <- file.path(OUT_DIR, "templates", paste0(REF_ARM, ".csv"))
write_csv(template, template_file)
ref_gs <- sw_gate_run(
  fcs_files = manifest$file,
  gating_template = template_file,
  transform_channels = c(VIABILITY, CD3, CD4, CD8, TNF, IFNG),
  cofactor = COFACTOR
)
sampleNames(ref_gs) <- manifest$sample
iwalk(QUADRANTS, \(node, alias) {
  gs_pop_set_name(
    ref_gs,
    paste0("/nonDebris/singlets/live/CD3pos/", node),
    alias
  )
})
dir.create(file.path(OUT_DIR, "gs"))
save_gs(ref_gs, file.path(OUT_DIR, "gs", REF_ARM))

read_x <- \(path) {
  exprs(read.FCS(path, transformation = FALSE, truncate_max_range = FALSE))
}
key_of <- \(e) do.call(paste, c(as.data.frame(e[, KEY_COLS]), sep = "|"))

# Per-event labels, one table per donor ------------------------------------------------------------------
labels <- map(set_names(donors), \(d) {
  gh <- ref_gs[[d]]
  root <- exprs(gh_pop_get_data(gh, "root"))
  n_fcs <- as.integer(read.FCSheader(manifest$file[manifest$sample == d])[[1]][[
    "$TOT"
  ]])
  stopifnot("GatingSet rows are the FCS rows" = nrow(root) == n_fcs)
  # raw_row: row of the same event in the uncleaned file, matched on KEY_COLS as 88_spike_in_markers.R does.
  # The clean2 row set stays frozen; the check needs unique keys in both files, every clean key found,
  # and the same order, so clean row i maps to exactly one raw row.
  raw_key <- key_of(read_x(file.path(
    RAW_DIR,
    basename(manifest$file[manifest$sample == d])
  )))
  clean_key <- key_of(read_x(manifest$file[manifest$sample == d]))
  raw_row <- match(clean_key, raw_key)
  stopifnot(
    "clean-to-raw mapping is one-to-one and order-preserving" = !anyDuplicated(
      raw_key
    ) &&
      !anyDuplicated(clean_key) &&
      !anyNA(raw_row) &&
      !is.unsorted(raw_row, strictly = TRUE)
  )
  flags <- map(NODES, \(node) gh_pop_get_indices(gh, node)) |> as_tibble()
  tibble(
    sample = d,
    row = seq_len(nrow(root)),
    raw_row = raw_row,
    FSC_A = root[, "FSC-A"],
    SSC_A = root[, "SSC-A"]
  ) |>
    bind_cols(flags) |>
    mutate(CD3neg_live = live & !CD3pos)
}) |>
  list_rbind()
write_csv(labels, file.path(OUT_DIR, "labels.csv.gz"))

# Events per gate and donor, and % of parent -------------------------------------------------------------
n_events <- labels |>
  summarise(n = n(), .by = sample)
label_counts <- labels |>
  tidyr::pivot_longer(
    all_of(names(PARENT)),
    names_to = "gate",
    values_to = "flag"
  ) |>
  summarise(n = sum(flag), .by = c(sample, gate)) |>
  mutate(parent = PARENT[gate])
label_counts <- label_counts |>
  left_join(
    select(label_counts, sample, parent = gate, n_parent = n),
    by = join_by(sample, parent),
    relationship = "many-to-one"
  ) |>
  left_join(
    n_events,
    by = join_by(sample),
    suffix = c("", "_all"),
    relationship = "many-to-one"
  ) |>
  mutate(
    n_parent = if_else(parent == "all events", n_all, n_parent),
    pct_of_parent = 100 * n / n_parent,
    gate = factor(gate, levels = names(PARENT))
  ) |>
  select(-n_all) |>
  arrange(sample, gate)
write_csv(label_counts, file.path(OUT_DIR, "label_counts.csv"))

# Quadrant counts: DP and DN are not label columns, so they are counted here ---------------------------------
quad_counts <- map(set_names(donors), \(d) {
  gh <- ref_gs[[d]]
  n_cd3 <- sum(gh_pop_get_indices(gh, NODES[["CD3pos"]]))
  imap(QUADRANTS, \(node, alias) {
    n <- sum(gh_pop_get_indices(gh, paste0(NODES[["CD3pos"]], "/", alias)))
    tibble(quadrant = alias, n = n, n_cd3 = n_cd3, pct_of_cd3 = 100 * n / n_cd3)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "sample")
write_csv(quad_counts, file.path(OUT_DIR, "quadrant_counts.csv"))

# Gate bounds, old per-donor gates and pooled gates ------------------------------------------------------
# Bounds are in the GatingSet's coordinates: arcsinh (cofactor 6000) for the transformed channels, raw for scatter.
gate_nodes <- c(
  NODES,
  DP = paste0(NODES[["CD3pos"]], "/DP"),
  DN = paste0(NODES[["CD3pos"]], "/DN")
)
bounds_of <- \(gh, node) {
  g <- gh_pop_get_gate(gh, node)
  if (inherits(g, "polygonGate")) {
    b <- g@boundaries
    tibble(
      dim = colnames(b),
      lo = apply(b, 2, min),
      hi = apply(b, 2, max),
      vertices = paste(sprintf("%.6g,%.6g", b[, 1], b[, 2]), collapse = ";")
    )
  } else {
    tibble(
      dim = names(g@min),
      lo = unname(g@min),
      hi = unname(g@max),
      vertices = NA_character_
    )
  }
}
cuts_of <- \(gs) {
  map(set_names(sampleNames(gs)), \(d) {
    imap(gate_nodes, \(node, gate) bounds_of(gs[[d]], node)) |>
      list_rbind(names_to = "gate")
  }) |>
    list_rbind(names_to = "sample")
}
old_gs <- load_gs(file.path(GATED_POOLED, "gs", REF_ARM))
cuts <- list(old_per_donor = cuts_of(old_gs), pooled = cuts_of(ref_gs)) |>
  list_rbind(names_to = "version")
write_csv(cuts, file.path(OUT_DIR, "cuts.csv"))
# pooled gates are one gate per node: the cut must be identical on every donor
stopifnot(
  "pooled gates are identical on every donor" = cuts |>
    dplyr::filter(version == "pooled") |>
    summarise(n_cut = n_distinct(lo, hi, vertices), .by = c(gate, dim)) |>
    pull(n_cut) |>
    (\(x) all(x == 1))()
)

# Old vs pooled % of parent per gate and donor ------------------------------------------------------------
compare <- read_csv(
  file.path(OLD_DIR, "label_counts.csv"),
  show_col_types = FALSE
) |>
  select(sample, gate, n_old = n, pct_old = pct_of_parent) |>
  inner_join(
    select(label_counts, sample, gate, n_new = n, pct_new = pct_of_parent) |>
      mutate(gate = as.character(gate)),
    by = join_by(sample, gate),
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  mutate(delta_pct = pct_new - pct_old)
write_csv(compare, file.path(OUT_DIR, "compare_old_pooled.csv"))

# Row pairing of every arm with SpectroFlo (G8) ----------------------------------------------------------
# Reported, not enforced. Reference = the SpectroFlo clean2 FCS files; each arm's GatingSet root is compared
# on all rows. A GatingSet that cannot be read gets one row with its error in `note`.
ref_pair <- map(set_names(donors), \(d) {
  read_x(manifest$file[manifest$sample == d])[, PAIR_COLS, drop = FALSE]
})
arm_dirs <- c(
  list.files(file.path(GATED_POOLED, "gs"), full.names = TRUE),
  list.files(file.path(GATED_ABLATION, "gs"), full.names = TRUE)
)
pairing <- map(arm_dirs, \(dir) {
  arm <- paste(basename(dirname(dirname(dir))), basename(dir), sep = "/")
  tryCatch(
    {
      gs <- load_gs(dir)
      map(donors, \(d) {
        m <- exprs(gh_pop_get_data(gs[[d]], "root")[, PAIR_COLS])
        ref <- ref_pair[[d]]
        same_n <- nrow(m) == nrow(ref)
        same_col <- \(col) {
          same_n && identical(unname(m[, col]), unname(ref[, col]))
        }
        tibble(
          arm = arm,
          sample = d,
          n_rows_equal = same_n,
          time_equal = same_col("Time"),
          fsca_equal = same_col("FSC-A"),
          fsch_equal = same_col("FSC-H"),
          ssca_equal = same_col("SSC-A"),
          note = NA_character_
        )
      }) |>
        list_rbind()
    },
    error = \(e) {
      tibble(
        arm = arm,
        sample = NA_character_,
        n_rows_equal = NA,
        time_equal = NA,
        fsca_equal = NA,
        fsch_equal = NA,
        ssca_equal = NA,
        note = conditionMessage(e)
      )
    }
  )
}) |>
  list_rbind()
write_csv(pairing, file.path(OUT_DIR, "pairing.csv"))

# Gate figures per donor ---------------------------------------------------------------------------------
# Strategy: one ggcyto hexbin per gate, drawn on its parent population. The gates are pooled, so the same
# cut is drawn on every donor. The panel subtitle is the donor id only.
strategy_panels <- tribble(
  ~panel , ~node       , ~parent     , ~x      , ~y        , ~title     , ~alias      ,
       1 , "nonDebris" , "root"      , "FSC-A" , "SSC-A"   , "Debris"   , "nondebris" ,
       2 , "singlets"  , "nonDebris" , "FSC-A" , "FSC-H"   , "Singlets" , "singlet"   ,
       3 , "live"      , "singlets"  , "FSC-A" , VIABILITY , "Live"     , "live"      ,
       4 , "CD3pos"    , "live"      , CD3     , "SSC-A"   , "CD3+"     , "CD3pos"
)
hex_fill <- scale_fill_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
  trans = "sqrt",
  guide = "none"
)
panel_of <- \(gh, d, node, parent, x, y, title, pct) {
  ff <- cytoframe_to_flowFrame(gh_pop_get_data(gh, parent))
  identifier(ff) <- d
  p <- ggcyto(ff, aes(x = !!rlang::sym(x), y = !!rlang::sym(y))) +
    geom_hex(bins = BINS) +
    geom_gate(gh_pop_get_gate(gh, node), colour = "red", linewidth = FIG_LINE) +
    hex_fill +
    labs(title = sprintf("%s: %.1f%% of parent", title, pct)) +
    theme(aspect.ratio = 1)
  as.ggplot(p)
}
walk(donors, \(d) {
  panels <- strategy_panels |>
    mutate(
      pct = label_counts$pct_of_parent[match(
        paste(d, alias),
        paste(label_counts$sample, label_counts$gate)
      )]
    )
  plots <- pmap(panels, \(panel, node, parent, x, y, title, alias, pct) {
    panel_of(ref_gs[[d]], d, node, parent, x, y, title, pct)
  })
  cairo_pdf(tempfile(fileext = ".pdf"))
  g <- cowplot::plot_grid(plotlist = plots, nrow = 1)
  g <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label(
        sprintf("SpectroFlo pooled gating strategy, donor %s", d),
        fontface = "bold",
        size = 7,
        fontfamily = FIG_FONT,
        x = 0.01,
        hjust = 0
      ),
    g,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
  dev.off()
  save_panel(
    g,
    select(panels, panel, node, parent, x, y, pct),
    file.path(OUT_DIR, "figures"),
    paste0("strategy_", d),
    183,
    62
  )
})

# CD4 vs CD8 within CD3+: the four quadrant gates drawn as shaded rectangles from their own bounds
quad_rects <- cuts |>
  dplyr::filter(
    version == "pooled",
    sample == donors[[1]],
    gate %in% names(QUADRANTS)
  ) |>
  select(gate, dim, lo, hi) |>
  tidyr::pivot_wider(names_from = dim, values_from = c(lo, hi)) |>
  transmute(
    quadrant = gate,
    xmin = pmax(.data[[paste0("lo_", CD4)]], AXIS_T[[1]]),
    xmax = pmin(.data[[paste0("hi_", CD4)]], AXIS_T[[2]]),
    ymin = pmax(.data[[paste0("lo_", CD8)]], AXIS_T[[1]]),
    ymax = pmin(.data[[paste0("hi_", CD8)]], AXIS_T[[2]])
  )
# label anchor: the outer corner of each quadrant's window corner
quad_anchor <- tribble(
  ~quadrant , ~x   , ~y   , ~hjust , ~vjust ,
  "CD4T"    , Inf  , -Inf ,      1 ,      0 ,
  "CD8T"    , -Inf , Inf  ,      0 ,      1 ,
  "DP"      , Inf  , Inf  ,      1 ,      1 ,
  "DN"      , -Inf , -Inf ,      0 ,      0
)
# one filled rectangle per quadrant, constant fill so it does not touch the hexagon density scale
quad_layers <- pmap(quad_rects, \(quadrant, xmin, xmax, ymin, ymax) {
  geom_rect(
    data = tibble(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = QUAD_FILL[[quadrant]],
    alpha = 0.18,
    colour = QUAD_FILL[[quadrant]],
    linewidth = FIG_LINE
  )
})
walk(donors, \(d) {
  gh <- ref_gs[[d]]
  ff <- cytoframe_to_flowFrame(gh_pop_get_data(gh, "CD3pos"))
  identifier(ff) <- d
  pct <- quad_counts |>
    dplyr::filter(sample == d) |>
    left_join(
      quad_anchor,
      by = join_by(quadrant),
      relationship = "one-to-one",
      unmatched = "error"
    ) |>
    mutate(label = sprintf("%s %.1f%%", quadrant, pct_of_cd3))
  p <- ggcyto(ff, aes(x = !!rlang::sym(CD4), y = !!rlang::sym(CD8))) +
    geom_hex(bins = BINS) +
    hex_fill +
    quad_layers +
    geom_label(
      data = pct,
      aes(x = x, y = y, hjust = hjust, vjust = vjust, label = label),
      inherit.aes = FALSE,
      size = LABEL_PT,
      size.unit = "pt",
      family = FIG_FONT,
      fill = "white",
      alpha = 0.85,
      label.size = 0,
      label.padding = unit(0.1, "lines")
    ) +
    ggcyto_par_set(limits = list(x = AXIS_T, y = AXIS_T)) +
    labs(
      x = sprintf("CD4, %s (arcsinh)", CD4),
      y = sprintf("CD8, %s (arcsinh)", CD8),
      title = sprintf("SpectroFlo CD3+ T cells, donor %s", d)
    ) +
    theme(aspect.ratio = 1)
  cairo_pdf(tempfile(fileext = ".pdf"))
  g <- ggplotGrob(as.ggplot(p))
  dev.off()
  save_panel(
    g,
    select(pct, sample, quadrant, n, n_cd3, pct_of_cd3),
    file.path(OUT_DIR, "figures"),
    paste0("gate_cd4_cd8_", d),
    89,
    80
  )
})

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
