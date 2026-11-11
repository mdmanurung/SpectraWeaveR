# Purpose: quantify how each unmixing arm changes the TNF and IFNg staining pattern in CD4 T and CD8 T
#          cells. Cells are paired across arms (same raw events, same row order), so population and
#          cytokine+ membership come from the reference arm's GatingSet only and every arm is scored on
#          the same cells. Descriptive only: nothing is graded or ranked.
#          Derived choice (not from a plan row): pile-up fraction = share of reference negatives with
#          |linear value| < PILEUP_K x the reference arm's Nguyen rSD for that donor, population, cytokine.
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/31_cytokine_benchmark.R
#          from the repo root. The run dir must not exist.
# Inputs:  experiments/runs/2026-11-10_bmv-gated/gs/<arm>/   saved GatingSets [19a_gate_arm.R]
#          experiments/runs/2026-11-10_bmv-gated/manifest.csv  arm and sample list [01_manifest-gate.R]
#          experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R  display labels
# Outputs: experiments/runs/2026-11-10_bmv-cytokine/cytokine_metrics.csv    one row per arm, donor,
#            population, cytokine: all metrics, pinned and own-gate flags
#          experiments/runs/2026-11-10_bmv-cytokine/biology_metrics.csv     one row per arm, donor: biology check
#          experiments/runs/2026-11-10_bmv-cytokine/fig_pct_pos_own.png     (+ _source.csv)
#          experiments/runs/2026-11-10_bmv-cytokine/fig_stain_index.png     (+ _source.csv)
#          experiments/runs/2026-11-10_bmv-cytokine/fig_neg_rsd.png         (+ _source.csv)
#          experiments/runs/2026-11-10_bmv-cytokine/fig_pileup.png          (+ _source.csv)
#          experiments/runs/2026-11-10_bmv-cytokine/fig_emd.png             (+ _source.csv)
#          experiments/runs/2026-11-10_bmv-cytokine/fig_biology.png         (+ _source.csv)
#          experiments/runs/2026-11-10_bmv-cytokine/sessionInfo_31.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})
source("experiments/benchmarks/05_bmv-unmix-compare/arm_labels.R")

# Parameters ------------------------------------------------------------------------------
GATED_DIR <- "experiments/runs/2026-11-10_bmv-gated" # input run (user)
RUN_DIR <- "experiments/runs/2026-11-10_bmv-cytokine" # new run dir (user, never overwrite)
REFERENCE_ARM <- "as188_joint" # negatives and positives from this arm only (user, 2026-11-10)
COFACTOR <- 6000 # arcsinh cofactor, spectral flow (CLAUDE.md)
MIN_N <- 100 # a metric is NA when its population has fewer events (user)
PILEUP_K <- 0.5 # pile-up: |linear value| < K x reference rSD (derived choice, see header)
HYPERNEG_K <- 3 # hypernegative: below -K x reference rSD (as in 28_rank_metrics.R)
PIN_TOL <- 0.02 # reference cut within this of a window edge counts as pinned (user)
CUT_WINDOW <- c(1, 4) # TNF mindensity gate_range in 19a_gate_arm.R, arcsinh (also applied to IFNg)
FAIL_PCT_HIGH <- 50 # own gate marks > this % of a population positive: flagged as failed (my choice)
TNF_CH <- "BV510-A" # TNF channel (19a_gate_arm.R)
IFNG_CH <- "BV570-A" # IFNg channel (19a_gate_arm.R)
CYTOKINES <- c(TNF = TNF_CH, IFNg = IFNG_CH)
POPS <- c(CD4T = "CD4T", CD8T = "CD8T")
POP_PATH <- "/nonDebris/singlets/live/CD3pos/"
DONORS <- c("A3", "A6", "B11", "B2", "D1")
ARM_GROUPS <- list(
  SpectroFlo = "spectroflo",
  `AS 1.0.2` = c("as102_basic", "as102_var", "best102_neg", "best102_af"),
  `AS 1.8.8` = c("as188_joint", "as188_fast_csv", "best188_p1", "best188_p5"),
  Crossed = c("s188_c102", "s102_c188")
) # four code-version groups (user); order on the x axis
DONOR_COLOURS <- c(
  A3 = "#0072B2",
  A6 = "#D55E00",
  B11 = "#009E73",
  B2 = "#CC79A7",
  D1 = "#E69F00"
) # Okabe-Ito, one colour per donor

stopifnot("run dir must not exist" = !dir.exists(RUN_DIR))
dir.create(RUN_DIR, recursive = TRUE)

arms <- unlist(ARM_GROUPS, use.names = FALSE)
manifest <- read_csv(
  file.path(GATED_DIR, "manifest.csv"),
  show_col_types = FALSE
)
stopifnot("arms match manifest" = setequal(arms, manifest$arm))
stopifnot("arms have labels" = all(arms %in% names(ARM_LABELS)))

# Helpers ---------------------------------------------------------------------------------
nguyen_rsd <- \(x) unname(quantile(x, 0.84) - quantile(x, 0.5)) # Nguyen robust SD, linear scale
lin <- \(x) sinh(x) * COFACTOR
# Pinned: reference cut within PIN_TOL of a window edge
is_pinned <- \(cut) {
  abs(cut - CUT_WINDOW[[1]]) <= PIN_TOL | abs(cut - CUT_WINDOW[[2]]) <= PIN_TOL
}

# Reference arm: membership, values, cuts, rSD per donor ---------------------------------
gs_ref <- load_gs(file.path(GATED_DIR, "gs", REFERENCE_ARM))
stopifnot("reference has all donors" = setequal(sampleNames(gs_ref), DONORS))

ref <- map(set_names(DONORS), \(d) {
  gh <- gs_ref[[d]]
  x <- exprs(gh_pop_get_data(gh, "root"))[, unname(CYTOKINES), drop = FALSE]
  cells <- expand_grid(pop = names(POPS), cyto = names(CYTOKINES)) |>
    pmap(\(pop, cyto) {
      pop_path <- paste0(POP_PATH, POPS[[pop]])
      in_pop <- gh_pop_get_indices(gh, pop_path)
      is_pos <- gh_pop_get_indices(gh, paste0(pop_path, "/", cyto, "pos"))
      v <- x[, CYTOKINES[[cyto]]]
      cut <- gh_pop_get_gate(gh, paste0(pop_path, "/", cyto, "pos"))@min[[1]]
      rsd <- nguyen_rsd(lin(v[in_pop & !is_pos]))
      list(
        in_pop = in_pop,
        pos = in_pop & is_pos,
        neg = in_pop & !is_pos,
        cut = cut,
        rsd = rsd
      )
    })
  names(cells) <- paste(
    rep(names(POPS), each = length(CYTOKINES)),
    names(CYTOKINES),
    sep = "|"
  )
  list(n_root = nrow(x), cells = cells, x = x)
})

# One arm: reference-cell metrics and own-gate metrics -----------------------------------
score_arm <- function(arm) {
  gs <- load_gs(file.path(GATED_DIR, "gs", arm))
  stopifnot("arm has all donors" = setequal(sampleNames(gs), DONORS))
  per_donor <- map(set_names(DONORS), \(d) {
    gh <- gs[[d]]
    x <- exprs(gh_pop_get_data(gh, "root"))[, unname(CYTOKINES), drop = FALSE]
    stopifnot("same root events as reference" = nrow(x) == ref[[d]]$n_root)
    own <- function(pop, cyto) {
      pop_path <- paste0(POP_PATH, POPS[[pop]])
      list(
        in_pop = gh_pop_get_indices(gh, pop_path),
        pos = gh_pop_get_indices(gh, pop_path) &
          gh_pop_get_indices(gh, paste0(pop_path, "/", cyto, "pos"))
      )
    }
    metrics <- expand_grid(pop = names(POPS), cyto = names(CYTOKINES)) |>
      pmap(\(pop, cyto) {
        r <- ref[[d]]$cells[[paste(pop, cyto, sep = "|")]]
        v <- x[, CYTOKINES[[cyto]]]
        neg_lin <- lin(v[r$neg])
        o <- own(pop, cyto)
        n_neg <- sum(r$neg)
        n_pos <- sum(r$pos)
        n_pop <- sum(r$in_pop)
        rsd <- if (n_neg >= MIN_N) nguyen_rsd(neg_lin) else NA_real_
        q <- if (n_neg >= MIN_N) {
          quantile(neg_lin, c(0.025, 0.975), names = FALSE)
        } else {
          c(NA_real_, NA_real_)
        }
        tibble(
          arm = arm,
          donor = d,
          pop = pop,
          cyto = cyto,
          cut_ref = r$cut,
          pinned = is_pinned(r$cut),
          n_pop_ref = n_pop,
          n_neg_ref = n_neg,
          n_pos_ref = n_pos,
          n_pop_own = sum(o$in_pop),
          n_pos_own = sum(o$pos),
          pct_pos_own = if (sum(o$in_pop) >= MIN_N) {
            100 * sum(o$pos) / sum(o$in_pop)
          } else {
            NA_real_
          },
          neg_rsd = rsd,
          neg_rsd_q = (q[[2]] - q[[1]]) / 3.92, # check: 2.5-97.5 / 3.92
          stain_index = if (n_neg >= MIN_N && n_pos >= MIN_N) {
            (median(lin(v[r$pos])) - median(neg_lin)) / (2 * rsd)
          } else {
            NA_real_
          },
          pileup_frac = if (n_neg >= MIN_N && !is.na(r$rsd)) {
            mean(abs(neg_lin) < PILEUP_K * r$rsd)
          } else {
            NA_real_
          },
          hyperneg_pct = if (n_neg >= MIN_N && !is.na(r$rsd)) {
            100 * mean(neg_lin < -HYPERNEG_K * r$rsd)
          } else {
            NA_real_
          },
          # 1D earth mover's distance on arcsinh values, same cells in both arms
          emd = if (n_pop >= MIN_N) {
            mean(abs(
              sort(v[r$in_pop]) - sort(ref[[d]]$x[r$in_pop, CYTOKINES[[cyto]]])
            ))
          } else {
            NA_real_
          }
        )
      }) |>
      list_rbind()
    # Biology check from the arm's own gates
    own_idx <- expand_grid(pop = names(POPS), cyto = names(CYTOKINES)) |>
      pmap(\(pop, cyto) {
        set_names(list(own(pop, cyto)), paste(pop, cyto, sep = "|"))
      }) |>
      list_flatten()
    pct_tnf_neg <- \(pop) {
      ifng <- own_idx[[paste(pop, "IFNg", sep = "|")]]$pos
      tnf <- own_idx[[paste(pop, "TNF", sep = "|")]]$pos
      if (sum(ifng) >= MIN_N) {
        100 * (1 - sum(ifng & tnf) / sum(ifng))
      } else {
        NA_real_
      }
    }
    pct_ifng <- \(pop) {
      o <- own_idx[[paste(pop, "IFNg", sep = "|")]]
      if (sum(o$in_pop) >= MIN_N) 100 * sum(o$pos) / sum(o$in_pop) else NA_real_
    }
    biology <- tibble(
      arm = arm,
      donor = d,
      pct_tnf_neg_in_ifng_cd4t = pct_tnf_neg("CD4T"),
      pct_tnf_neg_in_ifng_cd8t = pct_tnf_neg("CD8T"),
      ifng_ratio_cd8t_cd4t = pct_ifng("CD8T") / pct_ifng("CD4T")
    )
    list(metrics = metrics, biology = biology)
  })
  list(
    metrics = map(per_donor, "metrics") |> list_rbind(),
    biology = map(per_donor, "biology") |> list_rbind()
  )
}

scored <- map(set_names(arms), score_arm)
metrics <- map(scored, "metrics") |>
  list_rbind() |>
  mutate(
    own_failed = !is.na(pct_pos_own) &
      (pct_pos_own > FAIL_PCT_HIGH | n_pos_own == 0),
    arm = factor(arm, levels = arms),
    donor = factor(donor, levels = DONORS),
    pop = factor(pop, levels = names(POPS)),
    cyto = factor(cyto, levels = names(CYTOKINES))
  )
biology <- map(scored, "biology") |>
  list_rbind() |>
  mutate(
    arm = factor(arm, levels = arms),
    donor = factor(donor, levels = DONORS)
  )

metric_cols <- c(
  "pct_pos_own",
  "neg_rsd",
  "stain_index",
  "pileup_frac",
  "hyperneg_pct",
  "emd"
)
stopifnot(
  "no metric column is all NA" = !any(map_lgl(metrics[metric_cols], \(v) {
    all(is.na(v))
  }))
)
stopifnot(
  "biology not all NA" = !all(is.na(biology$ifng_ratio_cd8t_cd4t)),
  "metric rows: 11 arms x 5 donors x 2 pops x 2 cytokines" = nrow(metrics) ==
    220
)

write_csv(metrics, file.path(RUN_DIR, "cytokine_metrics.csv"))
write_csv(biology, file.path(RUN_DIR, "biology_metrics.csv"))

# Figures ---------------------------------------------------------------------------------
# Plot frame has columns arm, donor, value, open, row, col. Open symbol = pinned reference gate (or a
# failed own gate, for the metrics that depend on the arm's own gate).
arm_caption <- str_c(
  imap_chr(ARM_GROUPS, \(a, g) str_c(g, ": ", str_c(a, collapse = ", "))),
  collapse = "\n"
)
label_caption <- str_c(
  arms,
  " = ",
  str_replace_all(ARM_LABELS[arms], "\n", " | "),
  collapse = "\n"
)
group_breaks <- cumsum(lengths(ARM_GROUPS)) + 0.5
group_breaks <- group_breaks[-length(group_breaks)]
theme_set(theme_bw(base_size = 8))

plot_metric <- function(frame, ylab, file, extra_note) {
  p <- frame |>
    ggplot(aes(arm, value)) +
    geom_vline(
      xintercept = group_breaks,
      colour = "grey70",
      linetype = "dashed"
    ) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.7,
      linewidth = 0.3,
      colour = "grey20"
    ) +
    geom_point(
      aes(colour = donor, shape = open),
      position = position_jitter(width = 0.15, height = 0, seed = 1),
      size = 1.6,
      stroke = 0.6,
      na.rm = TRUE
    ) +
    (if (n_distinct(frame$row) == 1) {
      facet_wrap(vars(col), scales = "free_y")
    } else {
      facet_grid(row ~ col, scales = "free_y")
    }) +
    scale_colour_manual(values = DONOR_COLOURS, name = "Donor") +
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 1),
      name = "Open = flagged"
    ) +
    labs(
      x = NULL,
      y = ylab,
      caption = str_c(
        "Crossbar = median over donors; points = donors; groups: ",
        str_replace_all(arm_caption, "\n", "; "),
        ".\n",
        extra_note,
        "\n",
        label_caption
      )
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      plot.caption = element_text(hjust = 0, size = 5),
      plot.caption.position = "plot"
    )
  ggsave(
    file.path(RUN_DIR, file),
    p,
    width = 220,
    height = 200,
    units = "mm",
    dpi = 200
  )
  write_csv(
    frame,
    file.path(RUN_DIR, str_replace(file, "\\.png$", "_source.csv"))
  )
}

cyto_frame <- \(col, flag_own) {
  metrics |>
    transmute(
      arm,
      donor,
      value = .data[[col]],
      open = if (flag_own) pinned | own_failed else pinned,
      row = pop,
      col = cyto
    )
}
note_ref <- "Open = reference gate pinned (cut within 0.02 of 1 or 4 arcsinh). Scored on as188_joint cells."
note_own <- "Open = reference gate pinned, or own gate marks > 50% positive or 0 events positive."

plot_metric(
  cyto_frame("pct_pos_own", TRUE),
  "% positive, own gate",
  "fig_pct_pos_own.png",
  note_own
)
plot_metric(
  cyto_frame("stain_index", FALSE),
  "Stain index (linear)",
  "fig_stain_index.png",
  note_ref
)
plot_metric(
  cyto_frame("neg_rsd", FALSE),
  "Negative Nguyen rSD (linear, MFI)",
  "fig_neg_rsd.png",
  note_ref
)
plot_metric(
  cyto_frame("pileup_frac", FALSE),
  "Fraction of negatives with |x| < 0.5 x reference rSD",
  "fig_pileup.png",
  note_ref
)
plot_metric(
  cyto_frame("emd", FALSE),
  "EMD vs as188_joint (arcsinh units)",
  "fig_emd.png",
  note_ref
)

# Own-gate failure per arm, donor, population and cytokine (any failed gate behind a value flags it)
failed <- metrics |>
  summarise(failed = any(own_failed), .by = c(arm, donor, pop, cyto)) |>
  pivot_wider(names_from = c(pop, cyto), values_from = failed, names_sep = "_")
bio_flag <- biology |>
  left_join(
    failed,
    by = join_by(arm, donor),
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  transmute(
    arm,
    donor,
    pct_tnf_neg_in_ifng_cd4t = CD4T_TNF | CD4T_IFNg,
    pct_tnf_neg_in_ifng_cd8t = CD8T_TNF | CD8T_IFNg,
    ifng_ratio_cd8t_cd4t = CD4T_IFNg | CD8T_IFNg
  ) |>
  pivot_longer(-c(arm, donor), names_to = "panel", values_to = "open")
bio_frame <- biology |>
  pivot_longer(-c(arm, donor), names_to = "panel", values_to = "value") |>
  left_join(
    bio_flag,
    by = join_by(arm, donor, panel),
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  mutate(
    col = factor(
      panel,
      levels = c(
        "pct_tnf_neg_in_ifng_cd4t",
        "pct_tnf_neg_in_ifng_cd8t",
        "ifng_ratio_cd8t_cd4t"
      ),
      labels = c(
        "% TNF- among IFNg+, CD4T",
        "% TNF- among IFNg+, CD8T",
        "IFNg+ % ratio CD8T / CD4T"
      )
    ),
    row = "own gates"
  ) |>
  select(arm, donor, value, open, row, col)
plot_metric(
  bio_frame,
  "Value (percent, or ratio in right panel)",
  "fig_biology.png",
  "Each arm's own gates. Open = an own gate behind the value marks > 50% positive or 0 events positive."
)

writeLines(
  capture.output(sessionInfo()),
  file.path(RUN_DIR, "sessionInfo_31.txt")
)
