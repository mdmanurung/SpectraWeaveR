# Purpose: one GatingSet per arm on its clean3 files that carries the common SpectroFlo labels as gates
#          (plans/common_gates_plan.md: G7a, G8, G9, G10b, G11). Membership is by row: each node is a logical
#          vector from labels.csv.gz, so every arm is gated on exactly the SpectroFlo cells. Cytokine gates
#          TNFpos and IFNgpos are a fixed cut at arcsinh 1.5 (cofactor 6000) in the arm's own values.
#          Stops unless the arm's rows pair with SpectroFlo (Time, FSC-A, FSC-H, SSC-A on every row) and
#          every per-gate count equals the SpectroFlo label count. Replaces the pooled (P), ablation (A) and
#          older (G) gated dirs (runs/2026-11-10_bmv-gated-pooled, -ablation-gated, -gated).
#          One script for 3D, 24H and EXV. The cytokine cut of each arm, sample and channel comes from the
#          config rule (pos_rule: fixed, control_quantile or mindensity; neg_rule likewise) and is written to
#          cuts.csv, which cut_of() in datasets/dataset.R reads through the registry key "common-cuts".
# Run:     SW_DATASET=3d SW_RUN_DATE=2026-11-11 pixi run Rscript \
#            experiments/benchmarks/05_bmv-unmix-compare/97e_common_gs.R [out_dir] [arm1,arm2,...]
#          from the repo root; on SLURM: sbatch experiments/slurm/run_common_gs.sbatch 97e_common_gs.R
#          SW_DATASET is 3d, 24h or exv and SW_RUN_DATE is YYYY-MM-DD (neither has a default).
#          out_dir defaults to run_out("common-gated"); the arm list (default: all arms in the manifest) is for
#          smoke tests. Register the out dir under both keys common-gated and common-cuts after approval.
# Inputs:  registry keys of datasets/<SW_DATASET>_runs.csv:
#          run_in(cfg("manifest_key"))/manifest.csv   FCS of every arm, re-joined to the samples on the file stem
#                                  (3D: clean3 [97c]; 24H and ex vivo: unmix-compare-final [01], uncleaned files)
#          run_in(cfg("label_key"))/labels.csv.gz   gate labels (3D: SpectroFlo [97d], with DP and DN from gs/spectroflo
#                                  and quadrant_counts.csv; 24H and ex vivo: pooled AS 1.8.8 [99a])
#          run_in(cfg("label_key"))/label_counts.csv   expected counts
#          datasets/<SW_DATASET>_samples.csv, <SW_DATASET>.json   samples, channels, targets, pos_rule, neg_rule
# Outputs: <out_dir>/gs/<arm>/              GatingSet, one sample per donor, arcsinh on 6 channels
#          <out_dir>/gate_counts/<arm>.csv  events and parent events per node and donor
#          <out_dir>/manifest.csv           arm, dir, viability, file, sample (the clean3 manifest)
#          <out_dir>/pairing.csv            arm, sample, n_rows_equal, time/fsca/fsch/ssca_equal, counts_equal
#          <out_dir>/cuts.csv               arm, sample, channel, pos_cut, neg_cut; with a neg_rule of type "both"
#                                           one neg_cut_<variant> column per variant and neg_cut = NA; without a
#                                           neg_rule neg_cut = NA
#          <out_dir>/README.md              what this dir replaces and the node paths
#          <out_dir>/sessionInfo.txt

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(dplyr, exclude = "filter")
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R")
options(warn = 1) # warnings go to the log as they happen

# Parameters ---------------------------------------------------------------------------------------
OUT_KEY <- "common-gated" # registry key of the out dir; run_out() makes the new dated dir, which must not exist
MANIFEST_KEY <- cfg("manifest_key") # 3D: clean3 (97c); 24H and ex vivo: unmix-compare-final (step 01)
LABEL_KEY <- cfg("label_key") # 3D: common-gates-clean3 (97d); 24H and ex vivo: pop-gates (99a)
POOLED_LABELS <- identical(LABEL_KEY, "pop-gates") # 99a labels: uncleaned rows, no raw_row, DP and DN inside the table
CLEAN3_RUN <- run_in(MANIFEST_KEY)
LABEL_RUN <- run_in(LABEL_KEY)
REF_ARM <- "spectroflo" # the arm whose rows every arm must pair with
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
PAIR_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # must be equal on every row (audit, plan)
CD3 <- cfg("channels")$CD3 # channels as in 19a_gate_arm.R
CD4 <- cfg("channels")$CD4
CD8 <- cfg("channels")$CD8
# gate names come from the target labels: TNF and IFNg in 3D and 24H, CD1c and HLA-DR in EXV
CYTOKINES <- set_names(
  c(cfg("channels")$TNF_slot, cfg("channels")$IFNg_slot),
  c(cfg("targets")$tnf_label, cfg("targets")$ifng_label)
)
POPULATIONS <- c("CD4T", "CD8T")
POS_RULE <- cfg("pos_rule") # cut rule for the positive gates and shares (datasets/<ds>.json)
NEG_RULE <- DS$neg_rule # null while the user has not decided (EXV): neg_cut is NA, never a default
T_PARENT <- "/nonDebris/singlets/live/CD3pos"
# node table in hierarchy order: label column, alias and parent path
NODES <- tribble(
  ~label        , ~alias        , ~parent                    ,
  "nondebris"   , "nonDebris"   , "root"                     ,
  "singlet"     , "singlets"    , "/nonDebris"               ,
  "live"        , "live"        , "/nonDebris/singlets"      ,
  "CD3pos"      , "CD3pos"      , "/nonDebris/singlets/live" ,
  "CD4T"        , "CD4T"        , T_PARENT                   ,
  "CD8T"        , "CD8T"        , T_PARENT                   ,
  "DP"          , "DP"          , T_PARENT                   ,
  "DN"          , "DN"          , T_PARENT                   ,
  "CD3neg_live" , "CD3neg_live" , "/nonDebris/singlets/live"
) |>
  mutate(
    path = if_else(
      parent == "root",
      paste0("/", alias),
      paste0(parent, "/", alias)
    )
  )

cli_args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(cli_args) >= 1) cli_args[[1]] else run_out(OUT_KEY)
stopifnot("run dir is new" = !dir.exists(out_dir))
walk(file.path(out_dir, c("gs", "gate_counts")), \(d) {
  dir.create(d, recursive = TRUE)
})

# Labels, with DP and DN from the SpectroFlo GatingSet ---------------------------------------------
# the clean3 files are joined to the sample table on the file stem, one file per sample and arm
smp <- samples_tbl()
manifest <- read_csv(
  file.path(CLEAN3_RUN, "manifest.csv"),
  show_col_types = FALSE
) |>
  mutate(
    stem = basename(file) |>
      stringr::str_remove("\\.fcs$") |>
      stringr::str_remove("[ _]*AutoSpectral.*$")
  ) |>
  select(-sample) |>
  inner_join(
    select(smp, sample, stem),
    by = join_by(stem),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  select(-stem)
arms <- if (length(cli_args) >= 2) {
  strsplit(cli_args[[2]], ",")[[1]]
} else {
  unique(manifest$arm)
}
stopifnot(
  "arms are in the manifest" = all(arms %in% manifest$arm),
  REF_ARM %in% arms
)
donors <- sort(unique(manifest$sample), method = "radix") # C order, as arrange()
stopifnot("donors are the sample table" = setequal(donors, smp$sample))

labels <- read_csv(
  file.path(LABEL_RUN, "labels.csv.gz"),
  show_col_types = FALSE
)
if (POOLED_LABELS) {
  # 99a names: nonDebris and singlets columns; DP and DN are columns already
  labels <- labels |>
    rename(nondebris = nonDebris, singlet = singlets) |>
    split(~sample) |>
    map(\(d) {
      stopifnot("label rows are 1..n in order" = all(d$row == seq_len(nrow(d))))
      d
    })
} else {
  ref_gs <- load_gs(file.path(LABEL_RUN, "gs", REF_ARM))
  labels <- labels |>
    split(~sample) |>
    map(\(d) {
      gh <- ref_gs[[unique(d$sample)]]
      stopifnot(
        "label rows are 1..n in order" = all(d$row == seq_len(nrow(d)))
      )
      stopifnot(
        "label table and SpectroFlo gs have the same events" = nrow(d) ==
          nrow(gh_pop_get_data(gh, "root"))
      )
      d |>
        mutate(
          DP = gh_pop_get_indices(gh, paste0(T_PARENT, "/DP")),
          DN = gh_pop_get_indices(gh, paste0(T_PARENT, "/DN"))
        )
    })
}
stopifnot("labels cover the samples" = setequal(names(labels), donors))
# every child is inside its parent, so the logical gates form a valid hierarchy
walk(labels, \(d) {
  stopifnot(
    "labels are nested" = !any(d$singlet & !d$nondebris),
    !any(d$live & !d$singlet),
    !any(d$CD3pos & !d$live),
    !any((d$CD4T | d$CD8T | d$DP | d$DN) & !d$CD3pos),
    !any(d$CD3neg_live & !d$live)
  )
})

# Expected count per donor and node: label sums, checked against the 97d count tables
expected <- labels |>
  map(\(d) {
    d |>
      summarise(across(all_of(NODES$label), sum)) |>
      pivot_longer(everything(), names_to = "label", values_to = "n")
  }) |>
  list_rbind(names_to = "sample")
count_97d <- if (POOLED_LABELS) {
  # 99a label_counts.csv holds every node of the new hierarchy; keep the ones the common gates use
  read_csv(file.path(LABEL_RUN, "label_counts.csv"), show_col_types = FALSE) |>
    mutate(
      gate = case_when(
        gate == "nonDebris" ~ "nondebris",
        gate == "singlets" ~ "singlet",
        gate == "CD3neg" ~ "CD3neg_live",
        .default = gate
      )
    ) |>
    dplyr::filter(gate %in% NODES$label) |>
    select(sample, label = gate, n)
} else {
  bind_rows(
    read_csv(
      file.path(LABEL_RUN, "label_counts.csv"),
      show_col_types = FALSE
    ) |>
      select(sample, label = gate, n),
    read_csv(
      file.path(LABEL_RUN, "quadrant_counts.csv"),
      show_col_types = FALSE
    ) |>
      dplyr::filter(quadrant %in% c("DP", "DN")) |>
      select(sample, label = quadrant, n)
  )
}
stopifnot(
  "label sums equal the 97d count tables" = nrow(
    inner_join(
      expected,
      count_97d,
      by = join_by(sample, label, n),
      relationship = "one-to-one"
    )
  ) ==
    nrow(count_97d)
)

# Reference rows for the pairing check -------------------------------------------------------------
read_ff <- \(path) {
  read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
}
ref_files <- manifest |>
  dplyr::filter(arm == REF_ARM) |>
  arrange(sample)
ref_pair <- ref_files$file |>
  set_names(ref_files$sample) |>
  map(\(f) exprs(read_ff(f))[, PAIR_COLS, drop = FALSE])

# Cut rules ----------------------------------------------------------------------------------------
# A rule (cfg("pos_rule"), DS$neg_rule) is one of
#   fixed             value
#   control_quantile  per donor, quantile `prob` of the values of `control_arm` in sample(s) of condition
#                     `control_condition`, inside the label node `parent`; applied to every arm and to every
#                     sample (condition) of that donor
#   mindensity        openCyto mindensity valley per channel, pooled over all samples in the label node `parent`,
#                     in the values of the arm itself, or of `reference_arm` when the rule names one (then every
#                     arm gets that arm's cut); needs `gate_range` in the rule
#   by_channel        `rules`: one rule per channel (name = channel), each a fixed or mindensity rule; `parent`
#                     is the default parent node of the sub-rules
# neg_rule type "both" lists variants, each a rule or {same_as: "pos_rule"}; every variant gets its own column
# neg_cut_<name> in cuts.csv and neg_cut is NA.
resolve_rule <- \(r) if (is.null(r$same_as)) r else DS[[r$same_as]]
neg_rules <- if (is.null(NEG_RULE)) {
  list()
} else if (identical(NEG_RULE$type, "both")) {
  NEG_RULE$variants |>
    pmap(\(...) discard(list(...), \(x) length(x) == 1 && is.na(x))) |>
    map(resolve_rule) |>
    set_names(paste0("neg_cut_", NEG_RULE$variants$name))
} else {
  list(neg_cut = resolve_rule(NEG_RULE))
}
cut_rules <- c(list(pos_cut = POS_RULE), neg_rules)
CUT_CHANNELS <- unique(c(unname(CYTOKINES), unlist(map(cut_rules, "channels"))))
POS_CHANNELS <- if (is.null(POS_RULE$channels)) {
  CUT_CHANNELS
} else {
  POS_RULE$channels
}
# the positive gates sit under CD4T and CD8T, as every reader of the gate paths expects; a mindensity rule is
# fitted in its own parent node (EXV: CD3neg_live), so its gates are added there as well
GATE_PARENTS <- c(
  paste0(T_PARENT, "/", POPULATIONS),
  if (POS_RULE$type %in% c("mindensity", "by_channel")) {
    NODES$path[NODES$label == POS_RULE$parent]
  }
) |>
  unique()
stopifnot(
  "pos_rule parent is a label node" = is.null(POS_RULE$parent) ||
    POS_RULE$parent %in% NODES$label,
  "every cytokine channel has a pos_rule cut" = all(
    unname(CYTOKINES) %in% POS_CHANNELS
  )
)

# arcsinh values of `chans` in the events of label node `parent`, one matrix per sample
parent_values <- function(frame_of, rule_parent, chans, arm_name) {
  imap(frame_of, \(ff, s) {
    if (!all(chans %in% colnames(ff))) {
      stop(
        sprintf(
          "97e_common_gs.R: arm %s, sample %s lacks channel(s) %s",
          arm_name,
          s,
          paste(setdiff(chans, colnames(ff)), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    asinh(exprs(ff)[labels[[s]][[rule_parent]], chans, drop = FALSE] / COFACTOR)
  })
}

# control_quantile: one cut per donor and channel from the control arm; cached per rule, so the pos_rule and
# a negative variant that points at it read the control files once
ctrl_cache <- new.env()
control_cuts <- function(rule) {
  key <- paste(deparse(rule), collapse = "")
  if (is.null(ctrl_cache[[key]])) {
    if (!identical(rule$by, "donor")) {
      stop(
        "97e_common_gs.R: control_quantile supports by = donor only",
        call. = FALSE
      )
    }
    ctrl <- manifest |>
      dplyr::filter(arm == rule$control_arm) |>
      inner_join(
        dplyr::filter(smp, condition == rule$control_condition),
        by = join_by(sample),
        relationship = "one-to-one",
        unmatched = c("drop", "error") # other conditions drop out; every control sample needs a file
      )
    if (nrow(ctrl) != n_distinct(smp$donor)) {
      stop(
        sprintf(
          "97e_common_gs.R: control_quantile needs one %s sample of arm %s per donor, found %d for %d donors",
          rule$control_condition,
          rule$control_arm,
          nrow(ctrl),
          n_distinct(smp$donor)
        ),
        call. = FALSE
      )
    }
    ctrl_cache[[key]] <- pmap(ctrl, \(file, sample, donor, ...) {
      vals <- parent_values(
        list(read_ff(file)) |> set_names(sample),
        rule$parent,
        rule$channels,
        rule$control_arm
      )[[sample]]
      tibble(
        donor = donor,
        channel = rule$channels,
        cut = unname(apply(vals, 2, quantile, probs = rule$prob))
      )
    }) |>
      list_rbind()
  }
  ctrl_cache[[key]]
}

# frames of one named arm, all samples, cached: the reference for a mindensity rule with `reference_arm`
ref_cache <- new.env()
reference_frames <- function(arm_name) {
  if (is.null(ref_cache[[arm_name]])) {
    files <- manifest |>
      dplyr::filter(arm == arm_name) |>
      arrange(sample)
    stopifnot(
      "reference arm has one file per sample" = identical(files$sample, donors)
    )
    ref_cache[[arm_name]] <- set_names(map(files$file, read_ff), donors)
  }
  ref_cache[[arm_name]]
}

# cut of one rule: sample, channel, cut for this arm
rule_tbl <- function(rule, arm_name, frame_of) {
  chans <- if (is.null(rule$channels)) CUT_CHANNELS else rule$channels
  switch(
    rule$type,
    fixed = expand_grid(sample = donors, channel = chans, cut = rule$value),
    control_quantile = control_cuts(rule) |>
      inner_join(
        select(smp, sample, donor),
        by = join_by(donor),
        relationship = "many-to-many", # donor cuts (one per channel) x the donor's samples
        unmatched = "error"
      ) |>
      select(sample, channel, cut),
    mindensity = {
      if (is.null(rule$gate_range)) {
        stop(
          "97e_common_gs.R: mindensity rule needs gate_range in the pos_rule of dataset ",
          DS$dataset,
          " (the mindensity search window is the user's choice)",
          call. = FALSE
        )
      }
      fit_arm <- if (is.null(rule$reference_arm)) {
        arm_name
      } else {
        rule$reference_arm
      }
      fit_frames <- if (is.null(rule$reference_arm)) {
        frame_of
      } else {
        reference_frames(fit_arm)
      }
      pooled <- parent_values(fit_frames, rule$parent, chans, fit_arm) |>
        reduce(rbind)
      valley_cut <- map_dbl(set_names(chans), \(ch) {
        valley <- getExportedValue(
          "openCyto",
          paste0("gate_", cfg("mindensity_method"))
        )(
          flowFrame(pooled[, ch, drop = FALSE]),
          channel = ch,
          gate_range = rule$gate_range
        )@min[[1]]
        if (valley <= min(pooled[, ch])) NA_real_ else valley # a cut at the data minimum is no valley
      })
      expand_grid(sample = donors, channel = chans) |>
        mutate(cut = unname(valley_cut[channel]))
    },
    by_channel = imap(rule$rules, \(sub, ch) {
      sub$channels <- ch
      if (is.null(sub$parent)) {
        sub$parent <- rule$parent
      }
      rule_tbl(sub, arm_name, frame_of)
    }) |>
      list_rbind(),
    stop(
      "97e_common_gs.R: unknown cut rule type '",
      rule$type,
      "'",
      call. = FALSE
    )
  )
}

# cuts.csv rows of one arm: pos_cut, neg_cut and one neg_cut_<variant> per variant
arm_cut_table <- function(arm_name, frame_of) {
  tbl <- cut_rules |>
    imap(\(rule, col) {
      rule_tbl(rule, arm_name, frame_of) |>
        select(sample, channel, cut) |>
        set_names(c("sample", "channel", col))
    }) |>
    reduce(
      \(x, y) {
        left_join(
          x,
          y,
          by = join_by(sample, channel),
          relationship = "one-to-one"
        )
      },
      .init = expand_grid(sample = donors, channel = CUT_CHANNELS)
    )
  if (!"neg_cut" %in% names(tbl)) {
    tbl <- mutate(tbl, neg_cut = NA_real_, .after = pos_cut)
  }
  no_cut <- dplyr::filter(tbl, channel %in% POS_CHANNELS, is.na(pos_cut))
  if (nrow(no_cut) > 0) {
    stop(
      sprintf(
        "97e_common_gs.R: arm %s has no positive cut (no mindensity valley?) for %s",
        arm_name,
        paste(unique(no_cut$channel), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  mutate(tbl, arm = arm_name, .before = 1)
}

# One arm ------------------------------------------------------------------------------------------
gate_arm <- function(arm_name) {
  files <- manifest |>
    dplyr::filter(arm == arm_name) |>
    arrange(sample)
  stopifnot("arm has one file per sample" = identical(files$sample, donors))
  frames <- map(files$file, read_ff)

  # G8: rows pair with SpectroFlo on every row, or stop
  pairing <- map2(frames, donors, \(ff, d) {
    m <- exprs(ff)[, PAIR_COLS, drop = FALSE]
    ref <- ref_pair[[d]]
    same_n <- nrow(m) == nrow(ref)
    same_col <- \(col) same_n && identical(unname(m[, col]), unname(ref[, col]))
    tibble(
      arm = arm_name,
      sample = d,
      n_rows_equal = same_n,
      time_equal = same_col("Time"),
      fsca_equal = same_col("FSC-A"),
      fsch_equal = same_col("FSC-H"),
      ssca_equal = same_col("SSC-A")
    )
  }) |>
    list_rbind()
  if (!all(unlist(select(pairing, -arm, -sample)))) {
    stop(
      sprintf(
        "arm %s does not pair with SpectroFlo:\n%s",
        arm_name,
        paste(capture.output(print(pairing)), collapse = "\n")
      ),
      call. = FALSE
    )
  }
  stopifnot(
    "label length equals event count" = all(map2_lgl(frames, donors, \(ff, d) {
      nrow(ff) == nrow(labels[[d]])
    }))
  )

  # cut of each sample and channel from the config rules (fixed 1.5 in 3D, G7a); read before the frames
  # go into the GatingSet
  arm_cuts <- arm_cut_table(arm_name, set_names(frames, donors))

  gs <- GatingSet(flowSet(set_names(frames, donors)))
  sampleNames(gs) <- donors
  # same channels as 19a_gate_arm.R, so the arcsinh values match
  gs <- transform(
    gs,
    transformerList(
      c(files$viability[[1]], CD3, CD4, CD8, unname(CYTOKINES)),
      flow_trans("arcsinh", \(x) asinh(x / COFACTOR), \(x) sinh(x) * COFACTOR)
    )
  )
  # SpectroFlo labels as gates, in hierarchy order
  pwalk(NODES, \(label, alias, parent, path) {
    gate <- map(set_names(donors), \(d) labels[[d]][[label]])
    gs_pop_add(gs, gate, name = alias, parent = parent)
  })
  # cytokine cut in this arm's own values, one rectangle gate per sample (G7a: fixed 1.5 in 3D)
  for (parent_path in GATE_PARENTS) {
    iwalk(CYTOKINES, \(ch, cy) {
      cut_by_sample <- arm_cuts |>
        dplyr::filter(channel == ch) |>
        pull(pos_cut, name = sample)
      gates <- map(cut_by_sample[donors], \(cut) {
        rectangleGate(
          filterId = paste0(cy, "pos"),
          .gate = set_names(list(c(cut, Inf)), ch)
        )
      })
      gs_pop_add(gs, gates, name = paste0(cy, "pos"), parent = parent_path)
    })
  }
  recompute(gs)

  counts <- gs_pop_get_count_fast(gs) |>
    as_tibble() |>
    mutate(arm = arm_name, .before = 1)
  # per-gate counts equal the SpectroFlo label counts
  got <- counts |>
    inner_join(
      NODES,
      by = join_by(Population == path),
      relationship = "many-to-one"
    ) |>
    transmute(sample = name, label, n = Count)
  stopifnot(
    "every label node is in the GatingSet" = nrow(got) ==
      nrow(NODES) * length(donors),
    "gate counts equal the SpectroFlo label counts" = nrow(
      inner_join(
        got,
        expected,
        by = join_by(sample, label, n),
        relationship = "one-to-one"
      )
    ) ==
      nrow(expected)
  )
  write_csv(counts, file.path(out_dir, "gate_counts", paste0(arm_name, ".csv")))
  save_gs(gs, file.path(out_dir, "gs", arm_name))
  message(sprintf(
    "%s: %d nodes x %d donors, counts equal",
    arm_name,
    n_distinct(counts$Population),
    length(donors)
  ))
  list(pairing = mutate(pairing, counts_equal = TRUE), cuts = arm_cuts)
}

arm_results <- map(arms, \(a) {
  res <- gate_arm(a)
  gc() # release the arm's GatingSet and its temporary h5 files
  res
})
pairing <- list_rbind(map(arm_results, "pairing"))
write_csv(pairing, file.path(out_dir, "pairing.csv"))
write_csv(list_rbind(map(arm_results, "cuts")), file.path(out_dir, "cuts.csv"))
write_csv(
  dplyr::filter(manifest, arm %in% arms),
  file.path(out_dir, "manifest.csv")
)

# README -------------------------------------------------------------------------------------------
writeLines(
  c(
    sprintf("# Common-gate GatingSets (dataset %s)", DS$dataset),
    "",
    "This one dir replaces P, A and G (3D: runs/2026-11-10_bmv-gated-pooled, runs/2026-11-10_bmv-ablation-gated",
    "and runs/2026-11-10_bmv-gated). Analysis scripts load `gs/<arm>` from here. The old dirs stay as they are.",
    "",
    sprintf(
      "%d arms x %d samples, clean3 rows (%s). Made by 97e_common_gs.R.",
      length(arms),
      length(donors),
      CLEAN3_RUN
    ),
    sprintf(
      "Membership is by row from the SpectroFlo pooled labels (%s/labels.csv.gz),",
      LABEL_RUN
    ),
    "so every arm is gated on the same cells. Rows pair with SpectroFlo on Time, FSC-A, FSC-H and SSC-A (pairing.csv).",
    "Per-gate counts equal the SpectroFlo label counts exactly (gate_counts/<arm>.csv).",
    "",
    "Node paths:",
    "",
    paste0("- `", NODES$path, "`"),
    paste0(
      "- `",
      rep(GATE_PARENTS, each = length(CYTOKINES)),
      "/",
      rep(paste0(names(CYTOKINES), "pos"), times = length(GATE_PARENTS)),
      "`"
    ),
    "",
    if (identical(POS_RULE$type, "fixed")) {
      sprintf(
        "%s are a fixed cut at arcsinh %s (cofactor %d) on %s, in each arm's own values (G7a).",
        paste0(names(CYTOKINES), "pos", collapse = " and "),
        POS_RULE$value,
        COFACTOR,
        paste(CYTOKINES, collapse = " and ")
      )
    } else {
      sprintf(
        "%s use the cut of pos_rule type %s on %s, per arm, sample and channel (cuts.csv). The gates sit under %s.",
        paste0(names(CYTOKINES), "pos", collapse = " and "),
        POS_RULE$type,
        paste(CYTOKINES, collapse = " and "),
        paste(GATE_PARENTS, collapse = ", ")
      )
    },
    "They are rectangle gates (> cut). CD3neg_live is live and not CD3pos.",
    sprintf(
      "Transformed channels: viability, CD3, CD4, CD8, %s (arcsinh, cofactor %d); others are raw.",
      paste(CYTOKINES, collapse = ", "),
      COFACTOR
    )
  ),
  file.path(out_dir, "README.md")
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
