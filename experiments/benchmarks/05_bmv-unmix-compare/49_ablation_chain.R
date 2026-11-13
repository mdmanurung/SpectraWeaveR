# Purpose: ablation chain, block B (user plan, 2026-11-10). Start from the zip AutoSpectral 1.8.8 setup and
#          switch one more component to its 1.0.2 behaviour at each step, regenerating every object from that
#          step down, inside 1.8.8 code. Each arm then unmixes the raw files (joint, speed slow) so it can be scored
#          against as188_joint (start) and s102_c188 (target) by 51_score_ablation.R.
#            chain1_var     variants with 1.0.2 size: som.dim 10, n.cells 2000
#            chain2_filt    + variant event filter loosened: sim.threshold = sim.threshold.floor = 0.9
#                           (the kNN background matching has no off switch, so 1.0.2 variants are not reachable)
#            chain3_af      + AF regenerated without dedup / contaminant removal / unmixed AF:
#                           deduplicate = FALSE, remove.contaminants = FALSE, use.unmixed = FALSE
#            chain4_spec    + dye spectra by the 1.0.2 route inside 1.8.8: define.flow.control, clean.controls
#                           (universal.negative = TRUE, af.remove = FALSE), get.fluorophore.spectra
#            chain5_bead    + bead controls re-gated from the raw files with the widened rectangle
#            chain6a_neg    + cell controls and Unstained (Cells) re-gated from raw with a wider cell rectangle
#            chain6b_neg    dropped (user, 2026-11-10): H12 unstained holds the same events as Unstained (Cells)
#            chain7a_pregate chain6a + stained cell controls pre-gated with a landmark polygon (Unstained: rectangle only)
#                           (define.gate.landmarks), written as gated FCS copies used by every later step
#            chain7b_legacy chain6a + landmark gate passed to define.flow.control(gating.system = "landmarks")
#                           for the spectra only (get.spectral.variants takes no gate list)
#            chain8_legacy  chain7a + pipeline "legacy" (1.8.8's rewritten optimiser, not 1.0.2 code)
#          chain4 on use CONTROL_GATE_FIX, a control-file copy: the zip file shares one gate.name between
#          Unstained (Cells) and Negative_1 (Beads), which 1.8.8 define.flow.control rejects (user, 2026-11-10).
#            chain1_var_s43 chain1_var with seed 43, to bound regeneration noise
#          Step "check" runs after all arms: every arm's objects must differ from the previous step's.
# Run:     SW_DATASET=<3d|24h|exv> SW_RUN_DATE=<YYYY-MM-DD> pixi run Rscript experiments/benchmarks/05_bmv-unmix-compare/49_ablation_chain.R <arm | check> [condition]
#          from the repo root (SLURM: experiments/slurm/run_ablation_chain.sbatch). Output goes to
#          <data_dir>/unmix_cells/<arm>[/<condition>], which must not exist.
#          condition = last argument, required when the dataset has conditions (24H: one task per condition;
#          the arm's objects are regenerated in each condition folder, and check reads one condition's folders),
#          refused otherwise. SW_RUN_DATE is only read by scripts that write a run folder.
# Inputs:  <data_dir>/autoSpectral_v1.8.8/ [user zip]
#            Reference Group Clean/*.fcs + fcs_control_file.csv   controls after the zip rectangles
#            Reference Group/*.fcs                                 raw controls (chain5 on)
#          run key cell-gating: fcs_control_file_landmark.csv  chain7 on [50_cell_gate_plots.R, user-chosen gates]
#          run key ablation-chain: fcs_control_file_gatefix.csv  chain4 on; a hand-made control-file copy [user]
#            Reference Group/unstained_AF/*.fcs                    unstained sample for AF (the one file)
#            table_spectra/ config-stamp dye spectra and AF spectra
#            figure_spectral_variants/Spectral_variants.rds        zip variants (previous step of chain1)
#            Raw[/<condition>]/*.fcs, AutoSpectral_unmixed[/<condition>]/*.fcs  samples to unmix; as188_joint for verify
#          config control_gates: bead_wide, cell_wide rectangles and landmark n_cells, percentile
#          stamp file names (spectra, AF) from config stamps; samples from datasets/<ds>_samples.csv [00_samples.R]
# Outputs: <out_dir>/*.fcs                unmixed files
#          <out_dir>/objects/spectra.rds, af_spectra.rds, spectral_variants.rds   the arm's spectra set
#          <out_dir>/object_check.csv     variants per dye, BV510 variant angles, thresholds, AF rows, angles
#          <out_dir>/args.txt, verify.csv, sessionInfo.txt
#          <out_dir>/controls/, control_gate_log.csv   chain5 on: rebuilt controls and events kept per file
#          unmix_cells/chain_object_check[_<condition>].csv  step "check": differs-from-previous per arm

suppressPackageStartupMessages({
  library(AutoSpectral)
  library(flowCore)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Dataset config --------------------------------------------------------------------------------
source("experiments/benchmarks/05_bmv-unmix-compare/datasets/dataset.R") # DS, cfg(), data_path(), run_in(), run_out(), samples_tbl()

# Parameters ------------------------------------------------------------------------------------
DATA <- normalizePath(data_path()) # absolute: the script changes its working directory below
PIPELINE <- "joint" # zip 1.8.8 script (inherited: user zip)
SPEED <- "slow" # user rule 2026-11-10: always slow; under pipeline joint speed has no effect (repro188 fast = s188_c188_slow on A3)
N_CANDIDATES <- 10000L # zip get.spectra.automated call (not rerun here; documents the 143803 spectra)
VAR_SOM_DIM <- 10L # 1.0.2 get.spectral.variants value (factors.csv F5; user plan D12)
VAR_N_CELLS <- 2000L # same
SIM_THRESHOLD <- 0.9 # F7 loosened: sim.threshold = sim.threshold.floor (user plan D6)
# rectangles and landmark settings are dataset config (3D: bead_wide FSC 2e5-9e5, SSC 1.3e6-2.25e6, user 2026-11-10 O1;
# cell_wide FSC 4e5-4e6, SSC 1e5-1.6e6, user 2026-11-10 chain6; landmark n_cells 2000, percentile 70, user 2026-11-10 chain7)
BEAD_FSC <- cfg("control_gates")$bead_wide$fsc
BEAD_SSC <- cfg("control_gates")$bead_wide$ssc
CELL_FSC <- cfg("control_gates")$cell_wide$fsc
CELL_SSC <- cfg("control_gates")$cell_wide$ssc
LANDMARK_FILE <- normalizePath(
  file.path(run_in("cell-gating"), "fcs_control_file_landmark.csv")
) # gate.name lymphocyte (user); written by 50_cell_gate_plots.R
LANDMARK_N_CELLS <- cfg("control_gates")$landmark$n_cells # define.gate.landmarks n.cells
LANDMARK_PERCENTILE <- cfg("control_gates")$landmark$percentile # define.gate.landmarks percentile
# Control-file gate fix for define.flow.control (chain4 on): the zip file puts Unstained (Cells) and
# Negative_1 (Beads) in one gate.name, which 1.8.8 rejects. Fix in a copy, used from chain4 on (chain7 on: the
# landmark copy, which also groups Unstained (Cells) with the cells).
CONTROL_GATE_FIX <- normalizePath(
  file.path(run_in("ablation-chain"), "fcs_control_file_gatefix.csv")
) # hand-made copy: Unstained (Cells) moved to smallGate_2 (user, 2026-11-10)

# arm: cumulative switches (each step adds one to the previous), previous step and extraction seed
# (user plan D5: asp$bird.seed inside AutoSpectral, set.seed() around extraction)
s1 <- list(
  landmark_pregate = FALSE,
  filt = FALSE,
  af = FALSE,
  spec = FALSE,
  bead = FALSE,
  cell = FALSE,
  neg_h12 = FALSE,
  landmark = FALSE,
  pipeline = PIPELINE,
  seed = 42
)
s2 <- modifyList(s1, list(filt = TRUE))
s3 <- modifyList(s2, list(af = TRUE))
s4 <- modifyList(s3, list(spec = TRUE))
s5 <- modifyList(s4, list(bead = TRUE))
s6a <- modifyList(s5, list(cell = TRUE))
s6b <- modifyList(s6a, list(neg_h12 = TRUE)) # alternative to 6a: H12 as the cell universal negative
s7a <- modifyList(s6a, list(landmark_pregate = TRUE)) # user, 2026-11-10
s7b <- modifyList(s6a, list(landmark = TRUE)) # user, 2026-11-10
# 1.8.8 "legacy" is its rewritten per-cell optimiser, not 1.0.2 code (user, 2026-11-10)
s8 <- modifyList(s7a, list(pipeline = "legacy"))
ARMS <- list(
  chain1_var = c(prev = "zip", s1),
  chain2_filt = c(prev = "chain1_var", s2),
  chain3_af = c(prev = "chain2_filt", s3),
  chain4_spec = c(prev = "chain3_af", s4),
  chain5_bead = c(prev = "chain4_spec", s5),
  chain6a_neg = c(prev = "chain5_bead", s6a),
  chain7a_pregate = c(prev = "chain6a_neg", s7a),
  chain7b_legacy = c(prev = "chain6a_neg", s7b),
  chain8_legacy = c(prev = "chain7a_pregate", s8),
  chain1_var_s43 = c(prev = "chain1_var", modifyList(s1, list(seed = 43)))
)

# Inputs ---------------------------------------------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
arm <- cli_args[[1]]
stopifnot("arm is one of ARMS or check" = arm %in% c(names(ARMS), "check"))
# 24H: one task per condition subfolder (argument 2); 3D and EXV have no condition level
cond <- if (length(cli_args) >= 2) cli_args[[2]] else NULL
if (is.null(DS$conditions) && !is.null(cond)) {
  stop(
    "49_ablation_chain.R: dataset ",
    DS$dataset,
    " has no conditions; drop argument 2 ",
    cond,
    call. = FALSE
  )
}
if (!is.null(DS$conditions) && (is.null(cond) || !cond %in% DS$conditions)) {
  stop(
    "49_ablation_chain.R: needs argument 2, a condition of dataset ",
    DS$dataset,
    ": ",
    paste(DS$conditions, collapse = ", "),
    call. = FALSE
  )
}
samples <- samples_tbl()
if (!is.null(cond)) {
  samples <- filter(samples, .data$condition == .env$cond)
}
stopifnot(
  "1.8.8 code loaded" = as.character(packageVersion("AutoSpectral")) == "1.8.8"
)
# arm output folder; 24H adds the condition level
arm_dir <- \(a) {
  paste(c(file.path(DATA, "unmix_cells", a), cond), collapse = "/")
}
d188 <- file.path(DATA, "autoSpectral_v1.8.8")
clean_dir <- file.path(d188, "Reference Group Clean")
raw_ctrl_dir <- file.path(d188, "Reference Group")
ctrl_csv <- "fcs_control_file.csv"
h12_file <- list.files(
  file.path(raw_ctrl_dir, "unstained_AF"),
  pattern = "\\.fcs$",
  full.names = TRUE
)
stopifnot("one unstained AF file" = length(h12_file) == 1)
spectra_zip_file <- file.path(
  d188,
  "table_spectra",
  cfg("stamps")$as188_spectra
)
af_zip_file <- file.path(d188, "table_spectra", cfg("stamps")$as188_af)
# raw csv read: read.spectra() rescales rows (user, 2026-11-10)
read_raw <- \(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}
spectra_zip <- read_raw(spectra_zip_file)
# angle in degrees between each row of m and the same-named row of ref
angle_to <- \(m, ref) {
  common <- intersect(rownames(m), rownames(ref))
  map_dbl(set_names(common), \(f) {
    a <- m[f, colnames(ref)]
    b <- ref[f, ]
    acos(min(1, sum(a * b) / sqrt(sum(a^2) * sum(b^2)))) * 180 / pi
  })
}

# Step "check": each arm's objects differ from its previous step's -----------------------------------
if (arm == "check") {
  load_objects <- \(a) {
    if (a == "zip") {
      return(list(
        spectra = spectra_zip,
        af = read_raw(af_zip_file),
        variants = readRDS(file.path(
          d188,
          "figure_spectral_variants",
          "Spectral_variants.rds"
        ))
      ))
    }
    obj_dir <- file.path(arm_dir(a), "objects")
    list(
      spectra = readRDS(file.path(obj_dir, "spectra.rds")),
      af = readRDS(file.path(obj_dir, "af_spectra.rds")),
      variants = readRDS(file.path(obj_dir, "spectral_variants.rds"))
    )
  }
  # only arms that finished, with a finished previous step (chain4/5 can be blocked; see args)
  has_objects <- \(a) {
    a == "zip" ||
      file.exists(file.path(arm_dir(a), "objects", "spectral_variants.rds"))
  }
  done <- keep(ARMS, \(a) has_objects(a$prev))
  done <- done[map_lgl(names(done), has_objects)]
  check <- imap(done, \(a, name) {
    now <- load_objects(name)
    before <- load_objects(a$prev)
    tibble(
      arm = name,
      prev = a$prev,
      spectra_same = isTRUE(all.equal(
        now$spectra,
        before$spectra,
        check.attributes = FALSE
      )),
      af_same = isTRUE(all.equal(now$af, before$af, check.attributes = FALSE)),
      variants_same = isTRUE(all.equal(
        now$variants$variants,
        before$variants$variants
      )),
      thresholds_same = isTRUE(all.equal(
        now$variants$thresholds,
        before$variants$thresholds
      )),
      # chain8_legacy changes only the unmix pipeline, so its objects may equal chain7a's
      unmix_only = a$prev != "zip" && a$pipeline != ARMS[[a$prev]]$pipeline
    )
  }) |>
    list_rbind()
  write_csv(
    check,
    file.path(
      DATA,
      "unmix_cells",
      paste0(
        "chain_object_check",
        if (!is.null(cond)) paste0("_", cond),
        ".csv"
      )
    )
  )
  stopifnot(
    "every arm changes at least one object vs its previous step" = all(
      check$unmix_only |
        !(check$spectra_same &
          check$af_same &
          check$variants_same &
          check$thresholds_same)
    )
  )
  quit(save = "no")
}

a <- ARMS[[arm]]
out_dir <- arm_dir(arm)
stopifnot("out dir must not exist" = !dir.exists(out_dir))
obj_dir <- file.path(out_dir, "objects")
dir.create(obj_dir, recursive = TRUE)
setwd(obj_dir) # AutoSpectral writes tables and figures relative to the working directory

# Controls: zip Reference Group Clean, or a rebuilt folder (chain5 on) -------------------------------
# chain5_bead  beads re-gated from raw with the widened rectangle; cells from Reference Group Clean
# chain6a_neg  + cell controls and Unstained (Cells) re-gated from raw with the wider cell rectangle
# chain6b_neg  + H12 unstained, inside the same cell rectangle, as the cell universal negative
# chain7a/7b   + the landmark control-file copy (gate.name lymphocyte for cells and Unstained)
# Zombie NIR stays the 84,528-event Reference Group Clean file in every arm (user, 2026-11-10, D13).
control_dir <- clean_dir
if (a$bead) {
  control_dir <- file.path(out_dir, "controls")
  dir.create(control_dir)
  ctrl <- read.csv(
    if (a$landmark || a$landmark_pregate) LANDMARK_FILE else CONTROL_GATE_FIX,
    check.names = FALSE
  )
  bead_files <- ctrl$filename[ctrl$control.type == "beads"]
  regate_cells <- if (a$cell) {
    setdiff(
      ctrl$filename[ctrl$control.type == "cells"],
      "Zombie NIR (Cells).fcs"
    )
  }
  copy_files <- setdiff(ctrl$filename, c(bead_files, regate_cells))
  file.copy(
    file.path(clean_dir, copy_files),
    file.path(control_dir, copy_files)
  )
  # raw control inside a rectangle, written under the control-file name
  regate <- \(src, dest, fsc, ssc, kind) {
    ff <- read.FCS(src, transformation = FALSE, truncate_max_range = FALSE)
    kept <- Subset(
      ff,
      rectangleGate(filterId = kind, list("FSC-A" = fsc, "SSC-A" = ssc))
    )
    write.FCS(kept, file.path(control_dir, dest))
    tibble(file = dest, gate = kind, n_raw = nrow(ff), n_kept = nrow(kept))
  }
  gate_log <- bind_rows(
    map(bead_files, \(f) {
      regate(file.path(raw_ctrl_dir, f), f, BEAD_FSC, BEAD_SSC, "bead_wide")
    }),
    map(regate_cells, \(f) {
      regate(file.path(raw_ctrl_dir, f), f, CELL_FSC, CELL_SSC, "cell_wide")
    })
  )
  if (a$neg_h12) {
    # H12 joins the control file as the cell universal negative, a copy of the Unstained (Cells) row
    h12_name <- "H12 Unstained (Cells).fcs"
    gate_log <- bind_rows(
      gate_log,
      regate(h12_file, h12_name, CELL_FSC, CELL_SSC, "cell_wide")
    )
    h12_row <- ctrl[ctrl$filename == "Unstained (Cells).fcs", ]
    h12_row$filename <- h12_name
    ctrl <- rbind(ctrl, h12_row)
    ctrl$universal.negative[ctrl$control.type == "cells"] <- h12_name
  }
  write_csv(gate_log, file.path(out_dir, "control_gate_log.csv"))
  write.csv(ctrl, file.path(control_dir, ctrl_csv), row.names = FALSE)
}
control_file <- file.path(control_dir, ctrl_csv)
if (a$spec && !a$bead) {
  control_file <- CONTROL_GATE_FIX # chain4: zip controls, gate-fix control-file copy
}
asp <- get.autospectral.param(cytometer = "aurora", figures = FALSE)
t0 <- Sys.time()
set.seed(a$seed)

# Landmark gate for the lymphocyte group (chain7a, chain7b) ------------------------------------------
landmark_gate <- if (a$landmark || a$landmark_pregate) {
  define.gate.landmarks(
    control.file = control_file,
    control.dir = control_dir,
    asp = asp,
    gate.name = "lymphocyte",
    n.cells = LANDMARK_N_CELLS,
    percentile = LANDMARK_PERCENTILE,
    output.dir = file.path(obj_dir, "figure_gate")
  )
}
# chain7a: apply the polygon to the stained lymphocyte-group controls in place, so define.flow.control,
# clean.controls and get.spectral.variants all read the gated copies. The cell negative (fluorophore AF)
# keeps its rectangle gate only (user, 2026-11-10): the polygon left it 578 of 6,374 events, below
# AutoSpectral's 1,000-event floor.
if (a$landmark_pregate) {
  lymph_gate <- polygonGate(
    filterId = "landmark_lymphocyte",
    .gate = cbind("FSC-A" = landmark_gate$x, "SSC-A" = landmark_gate$y)
  )
  ctrl_now <- read.csv(control_file, check.names = FALSE)
  pregate_log <- map(
    ctrl_now$filename[
      ctrl_now$gate.name == "lymphocyte" & ctrl_now$fluorophore != "AF"
    ],
    \(f) {
      ff <- read.FCS(
        file.path(control_dir, f),
        transformation = FALSE,
        truncate_max_range = FALSE
      )
      kept <- Subset(ff, lymph_gate)
      write.FCS(kept, file.path(control_dir, f))
      tibble(file = f, n_before = nrow(ff), n_kept = nrow(kept))
    }
  ) |>
    list_rbind()
  write_csv(pregate_log, file.path(out_dir, "landmark_pregate_log.csv"))
  saveRDS(landmark_gate, file.path(obj_dir, "landmark_gate.rds"))
}

# Dye spectra: zip 143803, or the 1.0.2 route inside 1.8.8 (chain4 on) -------------------------------
spectra <- spectra_zip
if (a$spec) {
  # chain7b: the landmark gate enters through define.flow.control; every other arm uses its default
  # density gating (chain7a: on the landmark pre-gated files)
  flow_control_spec <- if (a$landmark) {
    define.flow.control(
      control_dir,
      control_file,
      asp,
      gating.system = "landmarks",
      gate.list = list(lymphocyte = landmark_gate)
    )
  } else {
    define.flow.control(control_dir, control_file, asp)
  }
  clean_control <- clean.controls(
    flow_control_spec,
    asp,
    universal.negative = TRUE,
    af.remove = FALSE,
    main.figures = FALSE
  )
  spectra <- get.fluorophore.spectra(
    clean_control,
    asp,
    use.clean.expr = TRUE,
    title = "Cells Matching Negative 188",
    figures = FALSE
  )
  spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE] # as in the 06/07/16/26 scripts
}

# AF spectra: zip 143933, or regenerated without dedup and contaminant removal (chain3 on) ---------------
af_spectra <- read_raw(af_zip_file)
if (a$af) {
  af_spectra <- get.af.spectra(
    unstained.sample = h12_file,
    asp = asp,
    spectra = spectra,
    refine = FALSE, # zip value
    deduplicate = FALSE,
    remove.contaminants = FALSE,
    use.unmixed = FALSE,
    figures = FALSE
  )
}

# Variants: 1.0.2 size in every arm; filter loosened from chain2 on ---------------------------------
variant_args <- list(
  control.dir = control_dir,
  control.def.file = control_file,
  asp = asp,
  spectra = spectra,
  som.dim = VAR_SOM_DIM,
  n.cells = VAR_N_CELLS,
  figures = FALSE
)
if (a$filt) {
  variant_args <- c(
    variant_args,
    list(sim.threshold = SIM_THRESHOLD, sim.threshold.floor = SIM_THRESHOLD)
  )
}
spectra_variants <- do.call(get.spectral.variants, variant_args)
extract_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
saveRDS(spectra, file.path(obj_dir, "spectra.rds"))
saveRDS(af_spectra, file.path(obj_dir, "af_spectra.rds"))
saveRDS(spectra_variants, file.path(obj_dir, "spectral_variants.rds"))

# Object check: what this arm's spectra set looks like -------------------------------------------------
# angle of each BV510 variant to the zip 143803 BV510 signature
bv510_var <- spectra_variants$variants[[cfg("targets")$tnf_dye]] # TNF dye, BV510 in 3D
bv510_ref <- spectra_zip[cfg("targets")$tnf_dye, colnames(bv510_var)]
bv510_cos <- drop(bv510_var %*% bv510_ref) /
  (sqrt(rowSums(bv510_var^2)) * sqrt(sum(bv510_ref^2)))
bv510_angles <- acos(pmin(1, bv510_cos)) * 180 / pi
object_check <- bind_rows(
  tibble(
    item = "variants_n",
    dye = names(spectra_variants$variants),
    value = map_dbl(spectra_variants$variants, nrow)
  ),
  tibble(
    item = "threshold",
    dye = names(spectra_variants$thresholds),
    value = unname(spectra_variants$thresholds)
  ),
  tibble(
    item = "bv510_variant_angle_deg",
    dye = paste0(cfg("targets")$tnf_dye, "_variant_", seq_along(bv510_angles)),
    value = unname(bv510_angles)
  ),
  tibble(
    item = "spectrum_angle_vs_143803_deg",
    dye = rownames(spectra),
    value = unname(angle_to(spectra, spectra_zip)[rownames(spectra)])
  ),
  tibble(item = "af_rows", dye = NA_character_, value = nrow(af_spectra))
) |>
  mutate(arm = arm, .before = 1)
write_csv(object_check, file.path(out_dir, "object_check.csv"))

# Unmix ---------------------------------------------------------------------------------------------
writeLines(
  c(
    paste("arm:", arm),
    paste("previous step:", a$prev),
    paste(
      "seed (set.seed around extraction):",
      a$seed,
      " asp$bird.seed:",
      asp$bird.seed
    ),
    paste("controls:", control_dir),
    paste("control file:", control_file),
    paste(
      "spectra:",
      if (a$spec) {
        "get.fluorophore.spectra after define.flow.control + clean.controls(universal.negative = TRUE, af.remove = FALSE)"
      } else {
        spectra_zip_file
      }
    ),
    paste(
      "af:",
      if (a$af) {
        "get.af.spectra(H12, refine = FALSE, deduplicate = FALSE, remove.contaminants = FALSE, use.unmixed = FALSE)"
      } else {
        af_zip_file
      }
    ),
    paste(
      "variants: get.spectral.variants(som.dim =",
      VAR_SOM_DIM,
      ", n.cells =",
      VAR_N_CELLS,
      if (a$filt) {
        paste(", sim.threshold = sim.threshold.floor =", SIM_THRESHOLD)
      } else {
        ", sim thresholds default"
      },
      ")"
    ),
    paste(
      "bead gate:",
      if (a$bead) {
        paste(
          "raw beads in FSC",
          paste(BEAD_FSC, collapse = "-"),
          "SSC",
          paste(BEAD_SSC, collapse = "-"),
          "; cells and Zombie NIR (84,528 events) from Reference Group Clean"
        )
      } else {
        "zip rectangle (Reference Group Clean)"
      }
    ),
    paste(
      "variants per dye:",
      paste(map_int(spectra_variants$variants, nrow), collapse = ",")
    ),
    paste("af rows:", nrow(af_spectra)),
    paste(
      "pipeline:",
      a$pipeline,
      if (a$pipeline == "legacy") {
        "(1.8.8 rewritten per-cell optimiser, not 1.0.2 code)"
      },
      " speed:",
      SPEED
    ),
    paste(
      "cell controls:",
      if (a$cell) {
        paste(
          "raw inside FSC",
          paste(CELL_FSC, collapse = "-"),
          "SSC",
          paste(CELL_SSC, collapse = "-")
        )
      } else {
        "zip rectangle"
      }
    ),
    paste(
      "cell universal negative:",
      if (a$neg_h12) {
        "H12 unstained inside the cell rectangle"
      } else {
        "Unstained (Cells)"
      }
    ),
    paste(
      "cell gating:",
      if (a$landmark || a$landmark_pregate) {
        paste(
          "define.gate.landmarks(n.cells =",
          LANDMARK_N_CELLS,
          ", percentile =",
          LANDMARK_PERCENTILE,
          ") with",
          LANDMARK_FILE,
          if (a$landmark_pregate) {
            "; polygon applied to the lymphocyte-group files (pre-gated copies), then define.flow.control density gating as chain6a"
          } else {
            "; gate passed to define.flow.control(gating.system = 'landmarks', gate.list) for spectra; variants read chain6a-rectangle files"
          }
        )
      } else {
        "zip control file, density gating in define.flow.control (chain4 on)"
      }
    ),
    "csv read: raw as.matrix(read.csv(row.names = 1, check.names = FALSE))",
    paste("extraction minutes:", extract_min)
  ),
  file.path(out_dir, "args.txt")
)
flow_control <- reload.flow.control(control_dir, control_file, asp)
t1 <- Sys.time()
unmix.folder(
  fcs.dir = paste(c(file.path(d188, "Raw"), cond), collapse = "/"),
  spectra = spectra,
  asp = asp,
  flow.control = flow_control,
  method = "AutoSpectral",
  af.spectra = af_spectra,
  spectra.variants = spectra_variants,
  pipeline = a$pipeline,
  speed = SPEED,
  output.dir = out_dir,
  parallel = TRUE,
  verbose = TRUE
)
cat(
  "unmix minutes:",
  round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1),
  "\n",
  file = file.path(out_dir, "args.txt"),
  append = TRUE
)

# Verify: same files, event counts and channel names as as188_joint ---------------------------------
ref_dir <- paste(
  c(file.path(d188, "AutoSpectral_unmixed"), cond),
  collapse = "/"
)
read_hdr <- \(path) {
  ff <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  list(n = nrow(ff), channels = colnames(ff))
}
# one unmixed file per sample: the stem is the delivered file name without its " AutoSpectral*.fcs" suffix (EXV: _v1-8-8);
# the join stops on any other file
verify <- tibble(
  ref_file = list.files(ref_dir, pattern = "\\.fcs$", recursive = TRUE)
) |>
  mutate(stem = sub(" AutoSpectral[^ ]*\\.fcs$", "", ref_file)) |>
  inner_join(
    select(samples, stem),
    by = "stem",
    relationship = "one-to-one",
    unmatched = "error"
  ) |>
  mutate(
    file = paste0(stem, " AutoSpectral.fcs"), # the arm's own file name
    ref = map(file.path(ref_dir, ref_file), read_hdr),
    new = map(file.path(out_dir, file), read_hdr),
    n_ref = map_int(ref, "n"),
    n_new = map_int(new, "n"),
    channels_equal = map2_lgl(ref, new, \(r, n) {
      identical(r$channels, n$channels)
    })
  ) |>
  select(file, n_ref, n_new, channels_equal)
write_csv(verify, file.path(out_dir, "verify.csv"))
stopifnot(
  "one unmixed file per sample" = nrow(verify) == nrow(samples),
  "event counts equal as188_joint" = all(verify$n_ref == verify$n_new),
  "channel names equal as188_joint" = all(verify$channels_equal)
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
