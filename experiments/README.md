# experiments/

Everything that is not the package. Reproducibility, benchmarks and real-data runs live here and nowhere else.
`R CMD build` ignores this folder (`.Rbuildignore`), and `experiments/runs/` is not tracked by git.

## Rules

1. **The package never reads this folder.** Nothing under `R/`, `tests/`, `inst/` or `vignettes/` refers to a path in `experiments/`.
2. **Experiments use the package only through its exported functions** (`library(SpectraWeaveR)` or `pkgload::load_all()`). They never `source()` files from `R/`.
3. **Unit tests stay in `tests/`.** They use mock data and must pass without `flowCore` or `AutoSpectral`. A check that needs real data or the full stack goes in `integration/`.
4. **Data stay outside the repository.** `datasets/registry.csv` lists URL, checksum and use. `datasets/download.sh` fetches and verifies them into `$SW_DATA_DIR`.
5. **Every run gets a new folder** `runs/YYYY-MM-DD_<slug>/` that must not exist beforehand. A finished run folder is never edited. Each run writes `report/provenance.csv` and `sessionInfo.txt`. Known gap: the git commit of the script is not recorded, and `provenance.csv` shows the SpectraWeaveR remote sha as `NA`.
6. **Scripts are `run.R` and `config*.R`** inside the benchmark folder they belong to (the `NN_verb-noun.R` naming was planned and not used; known gap, files not renamed). The shared runner is in `benchmarks/_common/`. Parameters sit at the top with a reason each. Impactful choices (thresholds, sample sizes, seeds, methods) are set in the benchmark's `README.md` plan, not hidden in code.
7. **Compute goes through SLURM** (`slurm/run_script.sbatch`), not the login node.

## Layout

```
experiments/
  README.md
  datasets/
    registry.csv            name, URL, checksum, size, use
    download.sh             fetch + verify into $SW_DATA_DIR
  benchmarks/
    01_mendeley-7color/     one folder per benchmark
      README.md             question, claim, plan, approved choices, status
      run.R                 driver; takes the run directory as its one argument
  integration/              real-data smoke tests of the package (need the pixi env)
  slurm/
    run_script.sbatch       sbatch <script.R> <run_dir>
  runs/                     outputs, one dated folder per run (git-ignored)
```

## Run one benchmark

```bash
pixi install                                   # once, from the repository root
bash experiments/datasets/download.sh autospectral_mendeley
sbatch experiments/slurm/run_script.sbatch \
  experiments/benchmarks/01_mendeley-7color/run.R \
  experiments/runs/$(date +%F)_mendeley-7color
```
