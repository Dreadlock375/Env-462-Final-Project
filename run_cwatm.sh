#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MINIFORGE="${KARAKORAM_MINIFORGE:-$PROJECT/miniforge3}"

# Keep application-level configuration and caches local.
export HOME="$PROJECT/.runtime-home"

export XDG_CACHE_HOME="$PROJECT/runtime/cache"
export XDG_CONFIG_HOME="$PROJECT/runtime/config"
export XDG_DATA_HOME="$PROJECT/runtime/share"

export TMPDIR="$PROJECT/runtime/tmp"
export PYTHONPYCACHEPREFIX="$PROJECT/runtime/python-cache"

export CONDA_PKGS_DIRS="$MINIFORGE/pkgs"
export CONDA_ENVS_PATH="$MINIFORGE/envs"

mkdir -p \
  "$HOME" \
  "$XDG_CACHE_HOME" \
  "$XDG_CONFIG_HOME" \
  "$XDG_DATA_HOME" \
  "$TMPDIR" \
  "$PYTHONPYCACHEPREFIX"

# Load Miniforge and activate CWatM.
source "$MINIFORGE/etc/profile.d/conda.sh"
conda activate "$MINIFORGE/envs/cwatm_env"

# CWatM expects to run from its source directory.
cd "$PROJECT/cwatm/source"

exec python run_cwatm.py "$@"
