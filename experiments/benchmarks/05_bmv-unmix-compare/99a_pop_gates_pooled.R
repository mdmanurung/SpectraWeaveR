# Purpose: gate nonDebris, singlets, live and the major cell populations with openCyto on pooled samples of the
#          reference arm (AS 1.8.8), then pass the gate labels to the other arms by row pairing, so every arm is
#          scored on the same events. The populations come from datasets/<ds>.json (population_gates); the panels
#          differ per dataset, so no marker is assumed to sit on the same fluorophore as in 3D.
# Run:     SW_DATASET=<24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/99a_pop_gates_pooled.R [n_samples]
#          from the repo root; on SLURM: sbatch experiments/slurm/run_clean3.sbatch 99a_pop_gates_pooled.R
#          n_samples (optional) fits on the first n samples only, for a smoke test; it never registers a run.
# Inputs:  run_in("unmix-compare-final")/manifest.csv   arm, sample, file, viability of every arm [01]
#          datasets/<SW_DATASET>.json   gate_ref_arm, population_gates, viability_range, mindensity_method, channels
# Outputs: <out_dir> = run_out("pop-gates")
#          <out_dir>/gs/<ref arm>/             pooled GatingSet
#          <out_dir>/template.csv              openCyto template used
#          <out_dir>/labels.csv.gz             sample, row, FSC_A, SSC_A and one logical column per gate node
#          <out_dir>/label_counts.csv          events per node and sample, % of parent
#          <out_dir>/cuts.csv                  gate bounds per node (first sample; the fit is pooled); polygon gates give vertex bounds
#          <out_dir>/gate_flags.csv            gates whose median % of parent is under 0.05 or over 99.95
#          <out_dir>/pairing.csv               arm, sample, rows equal, Time/FSC-A/FSC-H/SSC-A equal (all arms vs ref)
#          <out_dir>/figures/gate_populations.png ggcyto hexbins: parent marker vs gated marker, fitted gates drawn
#          <out_dir>/figures/gate_strategy.png ggcyto hexbins of debris, singlets, live, CD3 and the CD4 x CD8 quadrants, pooled
#          <out_dir>/figures/pct_of_parent.png % of parent per node, one point per sample
#          <out_dir>/sessionInfo.txt

# dataset.R first: cfg("lib_dirs") must be known before the other libraries load
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
.libPaths(c(cfg("lib_dirs"), .libPaths()))
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
  library(tidyr)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE) # sw_gate_run
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters ---------------------------------------------------------------------------------------
REF_ARM <- cfg("gate_ref_arm") # arm whose pooled gates every arm inherits (user 2026-11-12)
POP_GATES <- cfg("population_gates") # one row per extra gate, from the config
stopifnot("population_gates is set for this dataset" = !is.null(POP_GATES))
COFACTOR <- 6000 # arcsinh cofactor for spectral flow (CLAUDE.md)
DEBRIS_RANGE <- "c(2e5, 1.5e6)" # as 19a_gate_arm.R and 3D
LINEAGE_RANGE <- "c(0, 3)" # CD3 search window, as 3D
CD4_CD8_TMIX_ARGS <- "K=3, usePrior='no', quantile1=0.8, quantile3=0.8" # as 3D
VIABILITY_ARGS <- paste0("gate_range=", cfg("viability_range")) # user 2026-11-12, page Q17
MIN_METHOD <- cfg("mindensity_method")
PAIR_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A")
N_PLOT <- 30000 # events per sample in the density figures; display only
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cli_n <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
OUT_DIR <- run_out("pop-gates")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Samples and files --------------------------------------------------------------------------------
manifest <- read_csv(
  file.path(run_in("unmix-compare-final"), "manifest.csv"),
  show_col_types = FALSE
)
smp <- samples_tbl()
stopifnot(
  "every arm has one file per sample" = all(
    count(manifest, arm)$n == nrow(smp)
  ),
  "reference arm is in the manifest" = REF_ARM %in% manifest$arm
)
src <- dplyr::filter(manifest, arm == REF_ARM) |> arrange(sample)
if (!is.na(cli_n)) {
  src <- slice_head(src, n = cli_n)
}
VIABILITY <- src$viability[[1]]
CD3 <- cfg("channels")$CD3
CD4 <- cfg("channels")$CD4
CD8 <- cfg("channels")$CD8

# Template: base hierarchy, then the config rows ----------------------------------------------------
path_of <- function(alias, rows) {
  # full gate path of an alias, built from the parents in the template
  parent <- rows$parent[rows$alias == alias]
  if (parent == "root") {
    paste0("/", alias)
  } else {
    paste0(path_of(parent, rows), "/", alias)
  }
}
base <- tibble(
  alias = c("nonDebris", "singlets", "live", "CD3pos", "*"),
  pop = c("+", "+", "-", "+", "+/-+/-"),
  parent = c("root", "nonDebris", "singlets", "live", "CD3pos"),
  dims = c("FSC-A", "FSC-A,FSC-H", VIABILITY, CD3, paste(CD4, CD8, sep = ",")),
  gating_method = c(
    MIN_METHOD,
    "singletGate",
    MIN_METHOD,
    MIN_METHOD,
    "quadGate.tmix"
  ),
  gating_args = c(
    paste0("gate_range=", DEBRIS_RANGE),
    "wider_gate=TRUE",
    VIABILITY_ARGS,
    paste0("gate_range=", LINEAGE_RANGE),
    CD4_CD8_TMIX_ARGS
  )
)
extra <- POP_GATES |>
  pmap(\(alias, parent, marker, channel, pop, method, gate_range, cut, ...) {
    is_ref <- identical(method, "refGate")
    tibble(
      alias = alias,
      pop = pop,
      parent = parent,
      dims = channel,
      gating_method = if (is_ref) {
        "refGate"
      } else if (identical(method, "fixed")) {
        "boundary" # openCyto rectangle gate [min, Inf]; the cut is the user's value
      } else if (identical(method, "tail")) {
        "tailgate" # cytoUtils gate_tail on the pooled parent events (user 2026-11-13)
      } else {
        MIN_METHOD
      },
      gating_args = if (is_ref) {
        NA_character_ # filled below once every path is known
      } else if (identical(method, "fixed")) {
        paste0("min=", cut)
      } else if (identical(method, "tail")) {
        "strict=TRUE" # strict TRUE and FALSE gave the same cut in 99b
      } else {
        paste0("gate_range=", gate_range)
      }
    )
  }) |>
  list_rbind()
template <- bind_rows(base, extra)
# a refGate row points at the positive gate of the same channel: CD3neg is the complement of CD3pos
ref_args <- \(i) {
  # path of the positive gate on the same channel as refGate row i (CD3neg points at CD3pos)
  pos <- template$alias[
    template$dims == template$dims[[i]] &
      template$pop == "+" &
      template$gating_method != "refGate"
  ]
  stopifnot("refGate has one positive gate on its channel" = length(pos) == 1)
  sub("^/", "", path_of(pos, template))
}
is_ref <- template$gating_method == "refGate"
template$gating_args[is_ref] <- map_chr(which(is_ref), ref_args)
template <- template |>
  mutate(
    collapseDataForGating = TRUE,
    groupBy = nrow(src), # openCyto pools only when groupBy is not empty
    preprocessing_method = NA,
    preprocessing_args = NA
  )
write_csv(template, file.path(OUT_DIR, "template.csv"))
# openCyto 2.18 does not register the tail gate; cytoUtils holds it as .tailgate, so register it under its template name
if (any(template$gating_method == "tailgate")) {
  openCyto::register_plugins(
    fun = getFromNamespace(".tailgate", "cytoUtils"),
    methodName = "tailgate",
    dep = NA_character_,
    "gating"
  )
}

# Fit on the pooled samples ---------------------------------------------------------------------------
gate_channels <- unique(c(
  VIABILITY,
  CD3,
  CD4,
  CD8,
  POP_GATES$channel
))
gs <- sw_gate_run(
  fcs_files = src$file,
  gating_template = file.path(OUT_DIR, "template.csv"),
  transform_channels = gate_channels,
  cofactor = COFACTOR
)
sampleNames(gs) <- src$sample
quads <- c(
  CD4T = paste0(CD4, "+", CD8, "-"),
  CD8T = paste0(CD4, "-", CD8, "+"),
  DP = paste0(CD4, "+", CD8, "+"),
  DN = paste0(CD4, "-", CD8, "-")
)
iwalk(quads, \(node, alias) {
  gs_pop_set_name(gs, paste0("/nonDebris/singlets/live/CD3pos/", node), alias)
})
dir.create(file.path(OUT_DIR, "gs"))
save_gs(gs, file.path(OUT_DIR, "gs", REF_ARM))
nodes <- setdiff(gs_get_pop_paths(gs[[1]], path = "auto"), "root")

# Labels, one table per sample -------------------------------------------------------------------------
labels <- map(set_names(src$sample), \(s) {
  gh <- gs[[s]]
  root <- exprs(gh_pop_get_data(gh, "root"))
  stopifnot(
    "GatingSet rows are the FCS rows" = nrow(root) ==
      as.integer(read.FCSheader(src$file[src$sample == s])[[1]][["$TOT"]])
  )
  tibble(
    sample = s,
    row = seq_len(nrow(root)),
    FSC_A = root[, "FSC-A"],
    SSC_A = root[, "SSC-A"]
  ) |>
    bind_cols(as_tibble(map(set_names(nodes), \(n) gh_pop_get_indices(gh, n))))
}) |>
  list_rbind() |>
  mutate(CD3neg_live = CD3neg) # the name 97e and the scoring scripts read
write_csv(labels, file.path(OUT_DIR, "labels.csv.gz"))

# parent of every node, for % of parent; the quadrant nodes sit under CD3pos
parent_of <- set_names(template$parent, template$alias)
parent_of <- parent_of[names(parent_of) != "*"]
parent_of[["nonDebris"]] <- "all events"
parent_of <- c(parent_of, set_names(rep("CD3pos", length(quads)), names(quads)))
label_counts <- labels |>
  pivot_longer(
    all_of(names(parent_of)),
    names_to = "gate",
    values_to = "flag"
  ) |>
  summarise(n = sum(flag), .by = c(sample, gate)) |>
  mutate(parent = parent_of[gate]) |>
  left_join(
    labels |> summarise(n_all = n(), .by = sample),
    by = join_by(sample),
    relationship = "many-to-one"
  )
label_counts <- label_counts |>
  left_join(
    select(label_counts, sample, parent = gate, n_parent = n),
    by = join_by(sample, parent),
    relationship = "many-to-one"
  ) |>
  mutate(
    n_parent = if_else(parent == "all events", n_all, n_parent),
    pct_of_parent = 100 * n / n_parent
  ) |>
  select(-n_all)
write_csv(label_counts, file.path(OUT_DIR, "label_counts.csv"))
# Gates that keep (almost) nothing or everything are reported, never silent: gate_flags.csv and a message
gate_flags <- label_counts |>
  summarise(median_pct = median(pct_of_parent), .by = gate) |>
  dplyr::filter(median_pct < 0.05 | median_pct > 99.95)
write_csv(gate_flags, file.path(OUT_DIR, "gate_flags.csv"))
if (nrow(gate_flags) > 0) {
  message(
    "99a_pop_gates_pooled.R: gate(s) keeping under 0.05% or over 99.95% of the parent: ",
    paste(gate_flags$gate, collapse = ", "),
    " (see gate_flags.csv)"
  )
}

# Gate bounds of the first sample (the fit is pooled, every sample has the same gate)
cuts <- map(set_names(nodes), \(n) {
  g <- gh_pop_get_gate(gs[[1]], n)
  if (is(g, "rectangleGate")) {
    tibble(
      type = "rectangle",
      dim = names(g@min),
      min = unname(g@min),
      max = unname(g@max)
    )
  } else if (is(g, "polygonGate")) {
    # vertex bounds per dimension, so the quadrant and singlet gates can be redrawn from cuts.csv
    tibble(
      type = "polygon",
      dim = colnames(g@boundaries),
      min = unname(apply(g@boundaries, 2, min)),
      max = unname(apply(g@boundaries, 2, max))
    )
  } else {
    tibble(
      type = class(g)[[1]],
      dim = NA_character_,
      min = NA_real_,
      max = NA_real_
    )
  }
}) |>
  list_rbind(names_to = "node")
write_csv(cuts, file.path(OUT_DIR, "cuts.csv"))

# Pairing: every arm has the reference arm's rows --------------------------------------------------------
read_key <- \(f) {
  exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))[,
    PAIR_COLS,
    drop = FALSE
  ]
}
ref_key <- parallel::mclapply(
  set_names(src$file, src$sample),
  read_key,
  mc.cores = N_CORES
)
pairing <- manifest |>
  dplyr::filter(sample %in% src$sample, arm != REF_ARM) |>
  mutate(
    result = parallel::mclapply(
      seq_len(n()),
      \(i) {
        m <- read_key(file[[i]])
        r <- ref_key[[sample[[i]]]]
        tibble(
          n_rows_equal = nrow(m) == nrow(r),
          time_equal = nrow(m) == nrow(r) &&
            identical(unname(m[, "Time"]), unname(r[, "Time"])),
          fsca_equal = nrow(m) == nrow(r) &&
            identical(unname(m[, "FSC-A"]), unname(r[, "FSC-A"])),
          fsch_equal = nrow(m) == nrow(r) &&
            identical(unname(m[, "FSC-H"]), unname(r[, "FSC-H"])),
          ssca_equal = nrow(m) == nrow(r) &&
            identical(unname(m[, "SSC-A"]), unname(r[, "SSC-A"]))
        )
      },
      mc.cores = N_CORES
    )
  ) |>
  select(arm, sample, result) |>
  unnest(result)
write_csv(pairing, file.path(OUT_DIR, "pairing.csv"))
stopifnot(
  "pairing covers every non-reference arm and sample" = nrow(pairing) ==
    (n_distinct(manifest$arm) - 1) * nrow(src),
  "pairing has the five check columns" = ncol(pairing) == 7
)
if (!all(unlist(select(pairing, -arm, -sample)))) {
  stop(
    "99a_pop_gates_pooled.R: an arm does not pair with the reference arm; see pairing.csv",
    call. = FALSE
  )
}

# marker of each gating channel, for axis labels
marker_of <- c(
  set_names(POP_GATES$marker, POP_GATES$channel),
  set_names(c("CD3", "CD4", "CD8", "viability"), c(CD3, CD4, CD8, VIABILITY))
)
# axis label "marker (fluorochrome channel)"; scatter channels keep their own name
lab_of <- \(ch) {
  if (ch %in% names(marker_of)) sprintf("%s (%s)", marker_of[[ch]], ch) else ch
}
# plain population names for titles
pop_label <- c(
  set_names(POP_GATES$label, POP_GATES$alias),
  nonDebris = "non-debris events",
  singlets = "singlets",
  live = "live cells",
  CD3pos = "CD3+ T cells",
  CD4T = "CD4+ T cells",
  CD8T = "CD8+ T cells",
  DP = "CD4+CD8+ T cells",
  DN = "CD4-CD8- T cells"
)
# Figures: ggcyto on the fitted gate objects ---------------------------------------------------------------------
# One pooled flowFrame per parent population (N_PLOT events per sample); the gate drawn is the openCyto gate.
set.seed(1)
pooled_ff <- function(parent, chans) {
  map(set_names(src$sample), \(s) {
    m <- exprs(gh_pop_get_data(gs[[s]], parent))[, chans, drop = FALSE]
    if (nrow(m) > N_PLOT) m[sample(nrow(m), N_PLOT), , drop = FALSE] else m
  }) |>
    reduce(rbind) |>
    flowFrame() |>
    (\(ff) {
      identifier(ff) <- "pooled"
      # marker names as channel descriptions, so ggcyto labels read "SBUV510-A CD3" instead of repeating the channel
      desc <- unname(marker_of[colnames(ff)])
      ff@parameters@data$desc <- desc
      ff
    })()
}
aes_of <- \(x, y = NULL) {
  if (is.null(y)) {
    aes(x = !!rlang::sym(x))
  } else {
    aes(x = !!rlang::sym(x), y = !!rlang::sym(y))
  }
}
ref_aliases <- template$alias[template$gating_method == "refGate"] # complements of a plotted gate
hex_fill <- scale_fill_viridis_c(trans = "log10", guide = "none")
# Population gates as ggcyto hexbins: the parent's gating marker (x) against the gated marker (y), drawn on the
# events of the parent's parent, with the parent gate and the node gate from the fitted GatingSet
pop_rows <- template |>
  dplyr::filter(
    alias %in% extra$alias,
    gating_method != "refGate"
  )
pop_panels <- pmap(pop_rows, \(alias, parent, dims, ...) {
  x_ch <- template$dims[template$alias == parent]
  grand <- template$parent[template$alias == parent]
  pct <- median(label_counts$pct_of_parent[label_counts$gate == alias])
  title <- paste(
    strwrap(
      sprintf(
        "Gate: %s | gated within: %s | %.1f%% of parent (median)",
        pop_label[[alias]],
        pop_label[[parent]],
        pct
      ),
      70
    ),
    collapse = "\n"
  )
  # view 1: the parent's gating marker against the gated marker, events of the parent's parent
  p1 <- ggcyto(pooled_ff(grand, c(x_ch, dims)), aes_of(x_ch, dims)) +
    geom_hex(bins = 60) +
    geom_gate(gh_pop_get_gate(gs[[1]], parent), colour = "red") +
    geom_gate(gh_pop_get_gate(gs[[1]], alias), colour = "red") +
    hex_fill +
    labs(
      x = lab_of(x_ch),
      y = lab_of(dims),
      title = title,
      subtitle = sprintf(
        "events shown: %s",
        if (grand %in% names(pop_label)) pop_label[[grand]] else grand
      )
    ) +
    theme(
      plot.title = element_text(size = 6),
      plot.subtitle = element_text(size = 6)
    )
  # view 2: SSC-A against the gated marker, events of the parent only, so the distribution is easy to see
  p2 <- ggcyto(pooled_ff(parent, c("SSC-A", dims)), aes_of("SSC-A", dims)) +
    geom_hex(bins = 60) +
    geom_gate(gh_pop_get_gate(gs[[1]], alias), colour = "red") +
    hex_fill +
    labs(
      x = "SSC-A",
      y = lab_of(dims),
      title = sprintf("Same gate | events shown: %s", pop_label[[parent]]),
      subtitle = " "
    ) +
    theme(
      plot.title = element_text(size = 6),
      plot.subtitle = element_text(size = 6)
    )
  list(as.ggplot(p1), as.ggplot(p2))
})
# one file per gate for the decision page, and the overview
walk2(pop_panels, pop_rows$alias, \(pr, a) {
  ggsave(
    file.path(OUT_DIR, "figures", paste0("gate_", a, ".png")),
    cowplot::plot_grid(pr[[1]], pr[[2]], ncol = 2),
    width = 183,
    height = 90,
    units = "mm",
    dpi = 150
  )
})
pop_panels <- unlist(pop_panels, recursive = FALSE)
ggsave(
  file.path(OUT_DIR, "figures", "gate_populations.png"),
  cowplot::plot_grid(plotlist = pop_panels, ncol = 2),
  width = 183,
  height = 30 + 62 * length(pop_panels) / 2,
  units = "mm",
  dpi = 150
)
# 2D strategy: debris, singlets, live, CD3 and the CD4 x CD8 quadrants, pooled, with the fitted gates
strategy <- tribble(
  ~node       , ~parent     , ~x      , ~y        ,
  "nonDebris" , "root"      , "FSC-A" , "SSC-A"   ,
  "singlets"  , "nonDebris" , "FSC-A" , "FSC-H"   ,
  "live"      , "singlets"  , "FSC-A" , VIABILITY ,
  "CD3pos"    , "live"      , CD3     , "SSC-A"
)
strategy_panels <- pmap(strategy, \(node, parent, x, y) {
  ff <- pooled_ff(parent, c(x, y))
  p <- ggcyto(ff, aes_of(x, y)) +
    geom_hex(bins = 60) +
    geom_gate(gh_pop_get_gate(gs[[1]], node), colour = "red") +
    hex_fill +
    labs(
      x = lab_of(x),
      y = lab_of(y),
      title = sprintf(
        "Gate: %s | %.1f%% of parent (median)",
        pop_label[[node]],
        median(label_counts$pct_of_parent[label_counts$gate == node])
      ),
      subtitle = sprintf(
        "events shown: %s",
        if (parent == "root") "all events" else pop_label[[parent]]
      )
    ) +
    theme(
      plot.title = element_text(size = 6),
      plot.subtitle = element_text(size = 6)
    )
  as.ggplot(p)
})
walk2(strategy_panels, strategy$node, \(p, n) {
  ggsave(
    file.path(OUT_DIR, "figures", paste0("gate_", n, ".png")),
    p,
    width = 95,
    height = 90,
    units = "mm",
    dpi = 150
  )
})
quad_colour <- c(
  CD4T = "#2166ac",
  CD8T = "#b2182b",
  DP = "#762a83",
  DN = "#7f7f7f"
)
AXIS_T <- c(-3, 6) # CD4 / CD8 display window, arcsinh; the gates and counts use every CD3+ cell
p_quad <- ggcyto(pooled_ff("CD3pos", c(CD4, CD8)), aes_of(CD4, CD8)) +
  geom_hex(bins = 60) +
  hex_fill
for (q in names(quads)) {
  p_quad <- p_quad +
    geom_gate(
      gh_pop_get_gate(gs[[1]], q),
      colour = quad_colour[[q]],
      fill = quad_colour[[q]],
      alpha = 0.12
    )
}
p_quad <- p_quad +
  ggcyto_par_set(limits = list(x = AXIS_T, y = AXIS_T)) +
  labs(
    x = lab_of(CD4),
    y = lab_of(CD8),
    title = "Gates: CD4+ / CD8+ quadrants | gated within: CD3+ T cells",
    subtitle = "events shown: CD3+ T cells"
  ) +
  theme(
    plot.title = element_text(size = 6),
    plot.subtitle = element_text(size = 6)
  )
p_quad <- as.ggplot(p_quad)
ggsave(
  file.path(OUT_DIR, "figures", "gate_CD4_CD8.png"),
  p_quad,
  width = 95,
  height = 90,
  units = "mm",
  dpi = 150
)
ggsave(
  file.path(OUT_DIR, "figures", "gate_strategy.png"),
  cowplot::plot_grid(plotlist = c(strategy_panels, list(p_quad)), ncol = 2),
  width = 183,
  height = 200,
  units = "mm",
  dpi = 150
)
# % of parent per node and sample, one point per sample; colour = condition when the dataset has one
pct_tbl <- label_counts |>
  left_join(
    smp |>
      mutate(
        group = if ("condition" %in% names(smp) && !all(is.na(smp$condition))) {
          as.character(condition)
        } else {
          "all"
        }
      ) |>
      select(sample, group),
    by = join_by(sample),
    relationship = "many-to-one"
  ) |>
  mutate(gate = factor(gate, levels = rev(unique(label_counts$gate))))
p_pct <- ggplot(pct_tbl, aes(pct_of_parent, gate, colour = group)) +
  geom_point(position = position_jitter(height = 0.15, seed = 1), size = 0.8) +
  labs(x = "% of parent", y = NULL, colour = NULL) +
  theme(
    legend.position = if (n_distinct(pct_tbl$group) == 1) "none" else "right"
  )
ggsave(
  file.path(OUT_DIR, "figures", "pct_of_parent.png"),
  p_pct,
  width = 120,
  height = 100,
  units = "mm",
  dpi = 150
)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
