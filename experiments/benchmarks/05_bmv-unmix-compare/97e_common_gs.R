# Purpose: one GatingSet per arm on its clean3 files that carries the common SpectroFlo labels as gates
#          (plans/common_gates_plan.md: G7a, G8, G9, G10b, G11). Membership is by row: each node is a logical
#          vector from labels.csv.gz, so every arm is gated on exactly the SpectroFlo cells. Cytokine gates
#          TNFpos and IFNgpos are a fixed cut at arcsinh 1.5 (cofactor 6000) in the arm's own values.
#          Stops unless the arm's rows pair with SpectroFlo (Time, FSC-A, FSC-H, SSC-A on every row) and
#          every per-gate count equals the SpectroFlo label count. Replaces the pooled (P), ablation (A) and
#          older (G) gated dirs (runs/2026-11-10_bmv-gated-pooled, -ablation-gated, -gated).
# Run:     pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/97e_common_gs.R [out_dir] [arm1,arm2,...]
#          from the repo root; on SLURM: sbatch experiments/slurm/run_common_gs.sbatch 97e_common_gs.R
#          out_dir defaults to OUT_DIR; the arm list (default: all arms in the manifest) is for smoke tests.
# Inputs:  experiments/runs/2026-11-11_bmv-clean3/manifest.csv   clean3 FCS of every arm [97c]
#          experiments/runs/2026-11-11_bmv-common-gates-clean3/labels.csv.gz   SpectroFlo labels [97d]
#          experiments/runs/2026-11-11_bmv-common-gates-clean3/gs/spectroflo   DP and DN membership [97d]
#          experiments/runs/2026-11-11_bmv-common-gates-clean3/{label_counts,quadrant_counts}.csv   expected counts
# Outputs: <out_dir>/gs/<arm>/              GatingSet, one sample per donor, arcsinh on 6 channels
#          <out_dir>/gate_counts/<arm>.csv  events and parent events per node and donor
#          <out_dir>/manifest.csv           arm, dir, viability, file, sample (the clean3 manifest)
#          <out_dir>/pairing.csv            arm, sample, n_rows_equal, time/fsca/fsch/ssca_equal, counts_equal
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
options(warn = 1) # warnings go to the log as they happen

# Parameters ---------------------------------------------------------------------------------------
OUT_DIR <- "experiments/runs/2026-11-11_bmv-common-gated" # new dated dir; must not exist
CLEAN3_RUN <- "experiments/runs/2026-11-11_bmv-clean3"
LABEL_RUN <- "experiments/runs/2026-11-11_bmv-common-gates-clean3" # labels and SpectroFlo gs from 97d
REF_ARM <- "spectroflo" # the arm whose rows every arm must pair with
COFACTOR <- 6000 # arcsinh cofactor (CLAUDE.md, confirmed by user)
CYTOKINE_CUT <- 1.5 # arcsinh, every cytokine, arm and section (G7a, user 2026-11-11)
PAIR_COLS <- c("Time", "FSC-A", "FSC-H", "SSC-A") # must be equal on every row (audit, plan)
CD3 <- "SBUV510-A" # channels as in 19a_gate_arm.R
CD4 <- "Spark Blue 574-A"
CD8 <- "SBUV445-A"
CYTOKINES <- c(TNF = "BV510-A", IFNg = "BV570-A")
POPULATIONS <- c("CD4T", "CD8T")
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
out_dir <- if (length(cli_args) >= 1) cli_args[[1]] else OUT_DIR
stopifnot("run dir is new" = !dir.exists(out_dir))
walk(file.path(out_dir, c("gs", "gate_counts")), \(d) {
  dir.create(d, recursive = TRUE)
})

# Labels, with DP and DN from the SpectroFlo GatingSet ---------------------------------------------
manifest <- read_csv(
  file.path(CLEAN3_RUN, "manifest.csv"),
  show_col_types = FALSE
)
arms <- if (length(cli_args) >= 2) {
  strsplit(cli_args[[2]], ",")[[1]]
} else {
  unique(manifest$arm)
}
stopifnot(
  "arms are in the manifest" = all(arms %in% manifest$arm),
  REF_ARM %in% arms
)
donors <- sort(unique(manifest$sample))
stopifnot("5 donors" = length(donors) == 5)

ref_gs <- load_gs(file.path(LABEL_RUN, "gs", REF_ARM))
labels <- read_csv(
  file.path(LABEL_RUN, "labels.csv.gz"),
  show_col_types = FALSE
) |>
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
stopifnot("labels cover the 5 donors" = setequal(names(labels), donors))
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
count_97d <- bind_rows(
  read_csv(file.path(LABEL_RUN, "label_counts.csv"), show_col_types = FALSE) |>
    select(sample, label = gate, n),
  read_csv(
    file.path(LABEL_RUN, "quadrant_counts.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(quadrant %in% c("DP", "DN")) |>
    select(sample, label = quadrant, n)
)
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

# One arm ------------------------------------------------------------------------------------------
gate_arm <- function(arm_name) {
  files <- manifest |>
    dplyr::filter(arm == arm_name) |>
    arrange(sample)
  stopifnot("arm has the 5 donors" = identical(files$sample, donors))
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
  # fixed cytokine cut under CD4T and CD8T, in this arm's own values (G7a)
  for (pop in POPULATIONS) {
    iwalk(CYTOKINES, \(channel, cy) {
      gs_pop_add(
        gs,
        rectangleGate(
          filterId = paste0(cy, "pos"),
          .gate = set_names(list(c(CYTOKINE_CUT, Inf)), channel)
        ),
        parent = paste0(T_PARENT, "/", pop)
      )
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
  mutate(pairing, counts_equal = TRUE)
}

pairing <- map(arms, \(a) {
  res <- gate_arm(a)
  gc() # release the arm's GatingSet and its temporary h5 files
  res
}) |>
  list_rbind()
write_csv(pairing, file.path(out_dir, "pairing.csv"))
write_csv(
  dplyr::filter(manifest, arm %in% arms),
  file.path(out_dir, "manifest.csv")
)

# README -------------------------------------------------------------------------------------------
writeLines(
  c(
    "# Common-gate GatingSets (2026-11-11)",
    "",
    "This one dir replaces P, A and G: runs/2026-11-10_bmv-gated-pooled, runs/2026-11-10_bmv-ablation-gated",
    "and runs/2026-11-10_bmv-gated. Analysis scripts load `gs/<arm>` from here. The old dirs stay as they are.",
    "",
    sprintf(
      "%d arms x 5 donors, clean3 rows (runs/2026-11-11_bmv-clean3). Made by 97e_common_gs.R.",
      length(arms)
    ),
    "Membership is by row from the SpectroFlo pooled labels (runs/2026-11-11_bmv-common-gates-clean3/labels.csv.gz),",
    "so every arm is gated on the same cells. Rows pair with SpectroFlo on Time, FSC-A, FSC-H and SSC-A (pairing.csv).",
    "Per-gate counts equal the SpectroFlo label counts exactly (gate_counts/<arm>.csv).",
    "",
    "Node paths:",
    "",
    paste0("- `", NODES$path, "`"),
    paste0(
      "- `",
      rep(paste0(T_PARENT, "/", POPULATIONS), each = length(CYTOKINES)),
      "/",
      rep(paste0(names(CYTOKINES), "pos"), times = length(POPULATIONS)),
      "`"
    ),
    "",
    sprintf(
      "TNFpos and IFNgpos are a fixed cut at arcsinh %s (cofactor %d) on %s, in each arm's own values (G7a).",
      CYTOKINE_CUT,
      COFACTOR,
      paste(CYTOKINES, collapse = " and ")
    ),
    "They are rectangle gates (> cut), not data-driven. CD3neg_live is live and not CD3pos.",
    "Transformed channels: viability, CD3, CD4, CD8, BV510-A, BV570-A (arcsinh, cofactor 6000); others are raw."
  ),
  file.path(out_dir, "README.md")
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
