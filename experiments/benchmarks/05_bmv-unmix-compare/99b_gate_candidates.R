# Purpose: show candidate gates for populations whose marker has no clear valley, so the user can pick from plots.
#          Candidates: openCyto mindensity2 (three search windows), gate_quantile, a negative-peak cut (my
#          construction: density mode + k SD of the values left of the mode), cytoUtils gate_tail (strict TRUE and
#          FALSE), and, for the ex vivo CD1c gate, the quadGate.tmix and transitional (flowClust.2d) 2D gates of
#          the openCyto HowToAutoGating vignette. No candidate is chosen here; every setting below is a proposal.
# Run:     SW_DATASET=<24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/99b_gate_candidates.R
#          from the repo root; on SLURM: sbatch experiments/slurm/run_common_gs.sbatch 99b_gate_candidates.R
# Inputs:  run_in("pop-gates")/gs/<ref arm>, template.csv   pooled GatingSet and template [99a]
#          datasets/<SW_DATASET>.json   gate_ref_arm, population_gates (channel, marker, label of each target)
# Outputs: <out_dir> = run_out("gate-candidates")
#          <out_dir>/candidate_cuts.csv        node, candidate, kind, cut (NA for 2D gates), % of parent pooled and per sample
#          <out_dir>/candidate_failures.csv    candidates that did not run, with the error text
#          <out_dir>/figures/context_<node>.png  distribution with the parent gate only (parent marker view and SSC-A view)
#          <out_dir>/figures/grid_<node>.png     one small hexbin per candidate (parent marker vs target), three per row
#          <out_dir>/figures/grid_ssc_<node>.png the same candidates against SSC-A
#          <out_dir>/sessionInfo.txt

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
})
source("experiments/benchmarks/05_bmv-unmix-compare/lib/theme.R")

# Parameters ---------------------------------------------------------------------------------------
TARGETS <- list(
  `24h` = "TCRgdpos", # user 2026-11-12: show plots; also try the tail gate on CD3+ cells
  exv = "CD1cpos" # user 2026-11-12: try the transitional gate or quadGate.tmix (vignette); fixed 1.0 is the config
)[[cfg("dataset")]]
WINDOWS <- list(`0 to 3` = c(0, 3), `0.3 to 2.5` = c(0.3, 2.5), `none` = NULL) # mindensity2 search windows
QUANTILES <- c(0.99, 0.999) # gate_quantile probs
PEAK_K <- c(3, 4) # negative-peak cut: mode + k * SD left of the mode
TAIL_STRICT <- c(TRUE, FALSE) # cytoUtils gate_tail strict argument; other arguments stay at their defaults
TMIX_K <- 3 # quadGate.tmix components, as the vignette and the CD4 x CD8 gate
TRANS_K <- 6 # flowClust.2d components, as the vignette
TRANS_TARGET <- c(3, 2) # transitional gate target (x, y) in arcsinh units: the CD11c-high, CD1c-high region in the hexbins
TRANS_QUANTILE <- 0.95 # as the vignette
TRANS_TRANSLATION <- 0.15 # as the vignette
FIT_N <- 2e5 # events used to fit the 2D gates (random, seed FIT_SEED)
FIT_SEED <- 2
MIN_METHOD <- cfg("mindensity_method")
N_PLOT <- 30000 # events per sample in the figure; display only
OUT_DIR <- run_out("gate-candidates")
stopifnot("run dir is new" = !dir.exists(OUT_DIR))
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE)

# Target, parent and the events around them -----------------------------------------------------------------
gs <- load_gs(file.path(run_in("pop-gates"), "gs", cfg("gate_ref_arm")))
rows <- cfg("population_gates")
node <- rows[rows$alias == TARGETS, ]
stopifnot("one target row" = nrow(node) == 1)
template <- read_csv(
  file.path(run_in("pop-gates"), "template.csv"),
  show_col_types = FALSE
)
parent_gate <- gh_pop_get_gate(gs[[1]], node$parent)
stopifnot("parent gate is a rectangle" = is(parent_gate, "rectangleGate"))
x_channel <- names(parent_gate@min)
x_cut <- unname(parent_gate@min)
y_channel <- node$channel
grand <- template$parent[template$alias == node$parent]
samples <- sampleNames(gs)
# events of the parent's parent: x (parent marker), y (target), SSC-A; the parent is x >= x_cut
gp <- map(set_names(samples), \(s) {
  exprs(gh_pop_get_data(gs[[s]], grand))[,
    c(x_channel, y_channel, "SSC-A"),
    drop = FALSE
  ]
})
vals <- map(gp, \(m) m[m[, x_channel] >= x_cut, y_channel]) |>
  purrr::keep(\(v) length(v) > 0) # a sample with no events in the parent has no %
pooled <- unlist(vals, use.names = FALSE)

# Labels for the figures --------------------------------------------------------------------------------------
marker_of <- c(
  set_names(rows$marker, rows$channel),
  set_names(
    c("CD3", "CD4", "CD8"),
    c(cfg("channels")$CD3, cfg("channels")$CD4, cfg("channels")$CD8)
  )
)
lab_of <- \(ch) {
  if (ch %in% names(marker_of)) sprintf("%s (%s)", marker_of[[ch]], ch) else ch
}
pop_label <- c(
  set_names(rows$label, rows$alias),
  live = "live cells",
  CD3pos = "CD3+ T cells",
  nonDebris = "non-debris events",
  singlets = "singlets"
)
name_of <- \(a) if (a %in% names(pop_label)) pop_label[[a]] else a

# Candidate cuts on the pooled parent values ------------------------------------------------------------------
ff1d <- \(x) flowFrame(matrix(x, ncol = 1, dimnames = list(NULL, y_channel)))
mi2 <- getExportedValue("openCyto", paste0("gate_", MIN_METHOD))
mi2_cut <- \(w) {
  args <- list(ff1d(pooled), channel = y_channel)
  if (!is.null(w)) {
    args$gate_range <- w
  }
  do.call(mi2, args)@min[[1]]
}
dens <- density(pooled)
mode_x <- dens$x[which.max(dens$y)]
sd_left <- sqrt(mean((pooled[pooled <= mode_x] - mode_x)^2))
failures <- list()
safely <- \(label, expr) {
  tryCatch(expr, error = \(e) {
    failures[[length(failures) + 1]] <<- tibble(
      candidate = label,
      error = conditionMessage(e)
    )
    message(
      "99b_gate_candidates.R: candidate '",
      label,
      "' failed: ",
      conditionMessage(e)
    )
    NULL
  })
}
cut_cands <- c(
  set_names(
    map_dbl(WINDOWS, mi2_cut),
    paste0(MIN_METHOD, ", window ", names(WINDOWS))
  ),
  set_names(
    map_dbl(QUANTILES, \(p) {
      openCyto::gate_quantile(
        ff1d(pooled),
        channel = y_channel,
        probs = p
      )@min[[1]]
    }),
    paste0("gate_quantile ", QUANTILES)
  ),
  set_names(
    mode_x + PEAK_K * sd_left,
    paste0("negative peak + ", PEAK_K, " SD")
  ),
  compact(set_names(
    map(TAIL_STRICT, \(st) {
      g <- safely(
        paste0("gate_tail strict=", st),
        cytoUtils::gate_tail(ff1d(pooled), channel = y_channel, strict = st)
      )
      if (is.null(g)) NULL else unname(g@min)
    }),
    paste0("gate_tail strict=", TAIL_STRICT)
  )) |>
    unlist()
)
if (!is.na(node$cut)) {
  cut_cands <- c(
    cut_cands,
    set_names(node$cut, paste0("fixed ", node$cut, " (config)"))
  )
}
rect_gate <- \(cut) {
  rectangleGate(set_names(
    list(c(x_cut, Inf), c(cut, Inf)),
    c(x_channel, y_channel)
  ))
}

# 2D candidates (ex vivo CD1c only): quadGate.tmix and the transitional gate, fitted on the parent's parent events
gates_2d <- list()
if (identical(TARGETS, "CD1cpos")) {
  all_m <- reduce(gp, rbind)[, c(x_channel, y_channel)]
  set.seed(FIT_SEED)
  fit_fr <- flowFrame(all_m[
    sample(nrow(all_m), min(FIT_N, nrow(all_m))),
    ,
    drop = FALSE
  ])
  q <- safely(
    "quadGate.tmix",
    openCyto:::.quadGate.tmix(
      fit_fr,
      channels = c(x_channel, y_channel),
      K = TMIX_K,
      usePrior = "no"
    )
  )
  if (!is.null(q)) {
    want <- paste0(paste0(c(x_channel, y_channel), "+"), collapse = "")
    hit <- Filter(
      \(g) identical(g@filterId, want),
      lapply(seq_along(q), \(i) q[[i]])
    )
    stopifnot("one ++ quadrant gate" = length(hit) == 1)
    gates_2d[[sprintf("quadGate.tmix K=%d (++ quadrant)", TMIX_K)]] <- hit[[1]]
  }
  tr <- safely(
    "transitional",
    openCyto:::.flowClust.2d(
      fit_fr,
      channels = c(x_channel, y_channel),
      K = TRANS_K,
      transitional = TRUE,
      target = TRANS_TARGET,
      quantile = TRANS_QUANTILE,
      translation = TRANS_TRANSLATION,
      pp_res = NULL
    )
  )
  if (!is.null(tr)) {
    gates_2d[[sprintf(
      "transitional K=%d, target (%.1f, %.1f)",
      TRANS_K,
      TRANS_TARGET[[1]],
      TRANS_TARGET[[2]]
    )]] <- tr
  }
}

# Share of the parent's events inside each candidate (pooled and per sample) ------------------------------------
in_gate <- \(g, m) {
  as.logical(
    flowCore::filter(
      flowFrame(m[, c(x_channel, y_channel), drop = FALSE]),
      g
    )@subSet
  )
}
share <- \(g) {
  per <- map_dbl(gp, \(m) {
    p <- m[m[, x_channel] >= x_cut, , drop = FALSE]
    if (nrow(p) == 0) NA_real_ else 100 * mean(in_gate(g, p))
  })
  pool <- reduce(gp, rbind)
  pool <- pool[pool[, x_channel] >= x_cut, , drop = FALSE]
  tibble(
    pct_pooled = 100 * mean(in_gate(g, pool)),
    pct_sample_min = min(per, na.rm = TRUE),
    pct_sample_max = max(per, na.rm = TRUE)
  )
}
res_cut <- imap(cut_cands, \(cut, name) {
  bind_cols(
    tibble(node = TARGETS, candidate = name, kind = "cut", cut = cut),
    share(rect_gate(cut))
  )
}) |>
  list_rbind()
res_2d <- imap(gates_2d, \(g, name) {
  bind_cols(
    tibble(node = TARGETS, candidate = name, kind = "2D gate", cut = NA_real_),
    share(g)
  )
}) |>
  list_rbind()
res <- bind_rows(res_cut, res_2d)
write_csv(res, file.path(OUT_DIR, "candidate_cuts.csv"))
write_csv(
  if (length(failures) > 0) {
    list_rbind(failures)
  } else {
    tibble(candidate = character(), error = character())
  },
  file.path(OUT_DIR, "candidate_failures.csv")
)

# Figures: ggcyto hexbins, one row per candidate ---------------------------------------------------------------------
set.seed(1)
sub_m <- \(m) {
  if (nrow(m) > N_PLOT) m[sample(nrow(m), N_PLOT), , drop = FALSE] else m
}
as_ff <- \(m) {
  ff <- flowFrame(m)
  identifier(ff) <- "pooled"
  ff@parameters@data$desc <- unname(marker_of[colnames(ff)])
  ff
}
m1 <- reduce(map(gp, sub_m), rbind) # grand events: x, y, SSC-A
ff1 <- as_ff(m1[, c(x_channel, y_channel)])
par_m <- m1[m1[, x_channel] >= x_cut, , drop = FALSE]
ff2 <- as_ff(par_m[, c("SSC-A", y_channel)])
y_window <- range(
  quantile(pooled, c(0.0005, 0.9995)),
  res$cut[!is.na(res$cut) & res$cut > quantile(pooled, 0.0005)]
) # keeps every candidate cut in view
x_window <- quantile(m1[, x_channel], c(0.0005, 0.9995))
ssc_window <- quantile(m1[, "SSC-A"], c(0.0005, 0.9995))
hex_fill <- scale_fill_viridis_c(trans = "log10", guide = "none")
view <- \(ff, x, x_lim, gate, title, subtitle) {
  p <- ggcyto(ff, aes(x = !!rlang::sym(x), y = !!rlang::sym(y_channel))) +
    geom_hex(bins = 70) +
    hex_fill +
    ggcyto_par_set(limits = list(x = x_lim, y = y_window)) +
    labs(
      x = lab_of(x),
      y = lab_of(y_channel),
      title = title,
      subtitle = subtitle
    ) +
    theme(
      plot.title = element_text(size = 6),
      plot.subtitle = element_text(size = 6)
    )
  if (!is.null(gate)) {
    p <- p + geom_gate(gate, colour = "#b2182b", fill = NA, linewidth = 0.5)
  }
  as.ggplot(p)
}
head_txt <- sprintf(
  "Gate: %s | gated within: %s",
  name_of(TARGETS),
  name_of(node$parent)
)
cand_panels <- \(
  candidate,
  kind,
  cut,
  pct_pooled,
  pct_sample_min,
  pct_sample_max,
  ...
) {
  title <- paste(
    c(
      strwrap(head_txt, 65),
      strwrap(
        sprintf(
          "%s: %.1f%% of parent (samples %.1f to %.1f)",
          candidate,
          pct_pooled,
          pct_sample_min,
          pct_sample_max
        ),
        65
      )
    ),
    collapse = "\n"
  )
  sub1 <- sprintf("events shown: %s", name_of(grand))
  if (kind == "cut") {
    g1 <- rectangleGate(set_names(
      list(c(x_cut, x_window[[2]]), c(cut, y_window[[2]])),
      c(x_channel, y_channel)
    ))
    g2 <- rectangleGate(set_names(
      list(ssc_window, c(cut, y_window[[2]])),
      c("SSC-A", y_channel)
    ))
    list(
      view(ff1, x_channel, x_window, g1, title, sub1),
      view(
        ff2,
        "SSC-A",
        ssc_window,
        g2,
        "same cut",
        sprintf("events shown: %s", name_of(node$parent))
      )
    )
  } else {
    g <- gates_2d[[candidate]]
    inside <- in_gate(g, m1[, c(x_channel, y_channel), drop = FALSE]) &
      m1[, x_channel] >= x_cut
    ff_in <- as_ff(m1[inside, c("SSC-A", y_channel), drop = FALSE])
    list(
      view(ff1, x_channel, x_window, g, title, sub1),
      view(
        ff_in,
        "SSC-A",
        ssc_window,
        NULL,
        "same gate",
        "events shown: inside the gate and the parent"
      )
    )
  }
}
# Compact layout for the decision page: one small panel per candidate, three per row ------------------------------
# context: the distribution with the parent gate only; grid: parent-marker view of every candidate;
# grid_ssc: the same candidates against SSC-A
context_panels <- list(
  view(
    ff1,
    x_channel,
    x_window,
    parent_gate,
    "parent gate only, no candidate",
    sprintf("events shown: %s", name_of(grand))
  ),
  view(
    ff2,
    "SSC-A",
    ssc_window,
    NULL,
    "no candidate",
    sprintf("events shown: %s", name_of(node$parent))
  )
)
head_title <- \(txt) {
  cowplot::ggdraw() +
    cowplot::draw_label(txt, size = 8, fontface = "bold", x = 0.01, hjust = 0)
}
short_panels <- pmap(
  res,
  \(candidate, kind, cut, pct_pooled, ...) {
    ttl <- paste0(
      paste(strwrap(candidate, 38), collapse = "\n"),
      "\n",
      if (kind == "cut") sprintf("cut %.2f", cut) else "2D gate",
      ": ",
      sprintf("%.1f%% of parent", pct_pooled)
    )
    if (kind == "cut") {
      g1 <- rectangleGate(set_names(
        list(c(x_cut, x_window[[2]]), c(cut, y_window[[2]])),
        c(x_channel, y_channel)
      ))
      g2 <- rectangleGate(set_names(
        list(ssc_window, c(cut, y_window[[2]])),
        c("SSC-A", y_channel)
      ))
      list(
        view(ff1, x_channel, x_window, g1, ttl, NULL),
        view(ff2, "SSC-A", ssc_window, g2, ttl, NULL)
      )
    } else {
      g <- gates_2d[[candidate]]
      inside <- in_gate(g, m1[, c(x_channel, y_channel), drop = FALSE]) &
        m1[, x_channel] >= x_cut
      ff_in <- as_ff(m1[inside, c("SSC-A", y_channel), drop = FALSE])
      list(
        view(ff1, x_channel, x_window, g, ttl, NULL),
        view(ff_in, "SSC-A", ssc_window, NULL, ttl, NULL)
      )
    }
  }
)
n_rows <- ceiling(length(short_panels) / 3)
grid_fig <- \(panels, txt, file) {
  ggsave(
    file.path(OUT_DIR, "figures", file),
    cowplot::plot_grid(
      head_title(txt),
      cowplot::plot_grid(plotlist = panels, ncol = 3),
      ncol = 1,
      rel_heights = c(0.05, 1)
    ),
    width = 240,
    height = 18 + 72 * n_rows,
    units = "mm",
    dpi = 150,
    bg = "white"
  )
}
head_txt2 <- sprintf(
  "Gate: %s | gated within: %s",
  name_of(TARGETS),
  name_of(node$parent)
)
grid_fig(
  map(short_panels, 1),
  sprintf(
    "%s | events shown: %s | red = candidate gate",
    head_txt2,
    name_of(grand)
  ),
  paste0("grid_", TARGETS, ".png")
)
grid_fig(
  map(short_panels, 2),
  sprintf(
    "%s | SSC-A view | events shown: %s (for 2D gates: events inside the gate)",
    head_txt2,
    name_of(node$parent)
  ),
  paste0("grid_ssc_", TARGETS, ".png")
)
ggsave(
  file.path(OUT_DIR, "figures", paste0("context_", TARGETS, ".png")),
  cowplot::plot_grid(
    head_title(sprintf("%s | distribution before any candidate", head_txt2)),
    cowplot::plot_grid(plotlist = context_panels, ncol = 2),
    ncol = 1,
    rel_heights = c(0.08, 1)
  ),
  width = 183,
  height = 95,
  units = "mm",
  dpi = 150,
  bg = "white"
)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
