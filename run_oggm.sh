#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MINIFORGE="${KARAKORAM_MINIFORGE:-$PROJECT/miniforge3}"

# Activate the project's Conda environment
source "$MINIFORGE/etc/profile.d/conda.sh"
conda activate "$MINIFORGE/envs/oggm_env"

# Redirect project-specific configuration and caches
export HOME="$PROJECT/.runtime-home"
export XDG_CACHE_HOME="$PROJECT/runtime/cache"
export MPLCONFIGDIR="$PROJECT/runtime/matplotlib"
export TMPDIR="$PROJECT/runtime/tmp"
export OGGM_WORKDIR="$PROJECT/working"
export OGGM_DOWNLOAD_CACHE="$HOME/OGGM/download_cache"
export OGGM_DOWNLOAD_CACHE_RO=False
export OGGM_EXTRACT_DIR="$HOME/OGGM"
export OGGM_RGI_DIR="$HOME/OGGM/rgi"
export CONDA_PKGS_DIRS="$MINIFORGE/pkgs"

mkdir -p "$HOME" "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$TMPDIR"

cd "$PROJECT"

exec python "$@"
