# Purpose: common gates (plans/common_gates_plan.md, step 1; G1, G6, G8, G9). Export the SpectroFlo gate
#          labels per event so later steps copy them to every other unmixing arm by row, and draw the gate
#          figures per donor for approval. The gates are the ones 19a_gate_arm.R fitted on SpectroFlo
#          (debris, singlets, live and CD3+ by mindensity; CD4T / CD8T by quadGate.tmix K = 3 pooled over the
#          5 donors). Nothing is refitted: the saved SpectroFlo GatingSet is read, after its template is
#          checked against the 19a template. Also checks, per arm and donor, that rows pair with SpectroFlo.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/97_spectroflo_labels.R
#            from the repo root.
# Inputs:  experiments/runs/2026-11-10_bmv-gated-pooled/gs/<arm>/             GatingSets [19a_gate_arm.R]
#          experiments/runs/2026-11-10_bmv-gated-pooled/templates/spectroflo.csv   template used [19a]
#          experiments/runs/2026-11-10_bmv-gated-pooled/manifest.csv          FCS paths [01_manifest-gate.R]
#          experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R
# Outputs: <out_dir>/labels.csv.gz     one row per event and donor: sample, row (clean2), raw_row, FSC_A, SSC_A, gate flags
#          <out_dir>/label_counts.csv  events per gate and donor, and % of parent
#          <out_dir>/pairing.csv       per arm and donor: n rows, FSC-A and SSC-A identical to SpectroFlo
#          <out_dir>/figures/strategy_<donor>.png, .pdf, .csv   debris, singlets, live, CD3+ gates
#          <out_dir>/figures/gate_cd4_cd8_<donor>.png, .pdf, .csv   CD4 vs CD8 with the quadrant gates
#          <out_dir>/sessionInfo.txt

.libPaths(c("experiments/bmv_test/compare_3D/data/lib_ggplot2_35", .libPaths())) # ggplot2 3.5.2 for ggcyto
suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tibble)
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters ---------------------------------------------------------------------------------------
GATED_POOLED <- "experiments/runs/2026-11-10_bmv-gated-pooled" # 19a gating of every arm
OUT_DIR <- "experiments/runs/2026-11-11_bmv-common-gates"
RAW_DIR <- "experiments/bmv_test/compare_3D/data/autoSpectral_v1.8.8/Raw" # raw files, as 88_spike_in_markers.R
KEY_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # as 88_spike_in_markers.R:51; maps clean rows to raw rows
REF_ARM <- "spectroflo" # the arm whose gates all arms inherit (G1)
BINS <- 60 # hexagons per panel; display only
LABEL_PT <- 6 # in-panel % text, inside Nature's 5-7 pt
CD4 <- "Spark Blue 574-A" # channels as in 19a_gate_arm.R
CD8 <- "SBUV445-A"
CD3 <- "SBUV510-A"
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

# Reuse the existing SpectroFlo gates only if the template is the 19a one --------------------------------
# 19a writes the template from its constants; here the expected cells are compared to the saved file.
template <- read_csv(
  file.path(GATED_POOLED, "templates", paste0(REF_ARM, ".csv")),
  show_col_types = FALSE
)
stopifnot(
  "template is the 19a template (G1)" = identical(
    template$gating_method,
    c("mindensity", "singletGate", "mindensity", "mindensity", "quadGate.tmix")
  ) &&
    identical(
      template$gating_args,
      c(
        "gate_range=c(2e5, 1.5e6)",
        "wider_gate=TRUE",
        NA,
        "gate_range=c(0, 3)",
        "K=3, usePrior='no', quantile1=0.8, quantile3=0.8"
      )
    ) &&
    identical(
      template$dims,
      c(
        "FSC-A",
        "FSC-A,FSC-H",
        "LIVE DEAD NIR-A",
        CD3,
        paste(CD4, CD8, sep = ",")
      )
    ) &&
    template$groupBy[[5]] == 5
)
ref_gs <- load_gs(file.path(GATED_POOLED, "gs", REF_ARM))
donors <- sampleNames(ref_gs)
manifest <- read_csv(
  file.path(GATED_POOLED, "manifest.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(arm == REF_ARM)
stopifnot(
  "manifest donors match the GatingSet" = setequal(manifest$sample, donors)
)

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

# Row pairing of every arm with SpectroFlo (G8) ----------------------------------------------------------
# Reported, not enforced: this step lists mismatches; propagation (step 3) stops on them.
arms <- list.files(file.path(GATED_POOLED, "gs"))
pairing <- map(arms, \(a) {
  gs <- load_gs(file.path(GATED_POOLED, "gs", a))
  stopifnot("arm has the SpectroFlo donors" = setequal(sampleNames(gs), donors))
  map(donors, \(d) {
    m <- exprs(gh_pop_get_data(gs[[d]], "root"))
    ref <- labels |>
      dplyr::filter(sample == d)
    tibble(
      arm = a,
      sample = d,
      n_rows = nrow(m),
      n_rows_ref = nrow(ref),
      same_n = nrow(m) == nrow(ref),
      same_fsc = same_n && identical(unname(m[, "FSC-A"]), ref$FSC_A),
      same_ssc = same_n && identical(unname(m[, "SSC-A"]), ref$SSC_A)
    )
  }) |>
    list_rbind()
}) |>
  list_rbind()
write_csv(pairing, file.path(OUT_DIR, "pairing.csv"))

# Gate figures per donor ---------------------------------------------------------------------------------
# Strategy: one ggcyto hexbin per gate, drawn on its parent population, gate from the SpectroFlo GatingSet.
strategy_panels <- tribble(
  ~panel , ~node       , ~parent     , ~x      , ~y                , ~title     ,
       1 , "nonDebris" , "root"      , "FSC-A" , "SSC-A"           , "Debris"   ,
       2 , "singlets"  , "nonDebris" , "FSC-A" , "FSC-H"           , "Singlets" ,
       3 , "live"      , "singlets"  , "FSC-A" , "LIVE DEAD NIR-A" , "Live"     ,
       4 , "CD3pos"    , "live"      , CD3     , "SSC-A"           , "CD3+"
)
panel_of <- \(gh, panel_row, pct) {
  ff <- cytoframe_to_flowFrame(gh_pop_get_data(gh, panel_row$parent))
  p <- ggcyto(
    ff,
    aes(x = !!rlang::sym(panel_row$x), y = !!rlang::sym(panel_row$y))
  ) +
    geom_hex(bins = BINS) +
    geom_gate(
      gh_pop_get_gate(gh, panel_row$node),
      colour = "red",
      linewidth = FIG_LINE
    ) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      trans = "sqrt",
      guide = "none"
    ) +
    labs(title = sprintf("%s: %.1f%% of parent", panel_row$title, pct)) +
    theme(aspect.ratio = 1)
  as.ggplot(p)
}
walk(donors, \(d) {
  gh <- ref_gs[[d]]
  pcts <- label_counts |>
    dplyr::filter(sample == d)
  panels <- strategy_panels |>
    mutate(
      alias = c("nondebris", "singlet", "live", "CD3pos"),
      pct = pcts$pct_of_parent[match(alias, pcts$gate)]
    )
  plots <- pmap(panels, \(panel, node, parent, x, y, title, alias, pct) {
    panel_of(
      gh,
      tibble(node = node, parent = parent, x = x, y = y, title = title),
      pct
    )
  })
  cairo_pdf(tempfile(fileext = ".pdf"))
  g <- cowplot::plot_grid(plotlist = plots, nrow = 1)
  g <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label(
        sprintf("SpectroFlo gating strategy, donor %s", d),
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

# CD4 vs CD8 within CD3+, with the two quadrant gates and their % of CD3+
walk(donors, \(d) {
  gh <- ref_gs[[d]]
  ff <- cytoframe_to_flowFrame(gh_pop_get_data(gh, "CD3pos"))
  pct <- label_counts |>
    dplyr::filter(sample == d, gate %in% c("CD4T", "CD8T")) |>
    mutate(
      label = sprintf("%s %.1f%%", gate, pct_of_parent),
      x = if_else(gate == "CD8T", -Inf, Inf),
      y = if_else(gate == "CD8T", Inf, -Inf),
      hjust = if_else(gate == "CD8T", 0, 1),
      vjust = if_else(gate == "CD8T", 1, 0)
    )
  p <- ggcyto(ff, aes(x = !!rlang::sym(CD4), y = !!rlang::sym(CD8))) +
    geom_hex(bins = BINS) +
    geom_gate(
      gh_pop_get_gate(gh, "CD4T"),
      colour = "red",
      linewidth = FIG_LINE
    ) +
    geom_gate(
      gh_pop_get_gate(gh, "CD8T"),
      colour = "red",
      linewidth = FIG_LINE
    ) +
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
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      trans = "sqrt",
      guide = "none"
    ) +
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
    select(pct, sample, gate, n, n_parent, pct_of_parent),
    file.path(OUT_DIR, "figures"),
    paste0("gate_cd4_cd8_", d),
    89,
    80
  )
})

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
