#!/usr/bin/env bash
# Purpose: install the GitHub-only and CRAN-only R packages that conda does not carry
#          into the pixi environment (library: .pixi/Rlib).
# Run:     pixi install && pixi run setup     (or: bash tools/setup_pixi_env.sh)
# Outputs: .pixi/Rlib/ (R packages), .pixi/github_installs.tsv (package, repo, commit)
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
export R_LIBS_USER="$PWD/.pixi/Rlib" R_ENVIRON_USER="$PWD/.pixi/Renviron.none" R_PROFILE_USER="$PWD/.pixi/Rprofile.none" R_MAKEVARS_USER="$PWD/.pixi/Makevars.none"
mkdir -p "$R_LIBS_USER"
LOG="$PWD/.pixi/github_installs.tsv"
printf 'package\trepo\tref\tstatus\n' > "$LOG"

# Pinned refs. Empty = default branch head; the commit actually installed is recorded in the log.
AUTOSPECTRAL_REF=3c504c4560cf1f0fd5dce298eaabaafaebfae9be   # version 1.8.8

# CRAN packages with no conda build.
Rscript -e 'for (p in c("mcptools","btw")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")' \
  || echo "CRAN extras failed (non-fatal)"

gh_install() { # package repo ref
  local pkg="$1" repo="$2" ref="${3:-}"
  local spec="$repo"; [ -n "$ref" ] && spec="$repo@$ref"
  if Rscript -e "remotes::install_github('$spec', upgrade = 'never', dependencies = NA)"; then
    local sha; sha=$(Rscript -e "d <- packageDescription('$pkg'); cat(if (is.null(d\$RemoteSha)) NA else d\$RemoteSha)")
    printf '%s\t%s\t%s\tok\n' "$pkg" "$repo" "$sha" >> "$LOG"
  else
    printf '%s\t%s\t%s\tFAILED\n' "$pkg" "$repo" "$ref" >> "$LOG"
  fi
}

gh_install cyCombine        biosurf/cyCombine
gh_install AutoSpectral     DrCytometer/AutoSpectral "$AUTOSPECTRAL_REF"
gh_install AutoSpectralRcpp DrCytometer/AutoSpectralRcpp
gh_install flowstate        nlaniewski/flowstate
gh_install spectracle       nlaniewski/spectracle
gh_install FastPG           sararselitsky/FastPG   # best effort
cat "$LOG"
