#!/usr/bin/env bash
# Historical workflow; run from any directory. No old scripts are needed.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
MINIFORGE="$PWD/miniforge3"
if [[ -n "${KARAKORAM_MINIFORGE:-}" && "$KARAKORAM_MINIFORGE" != "$MINIFORGE" ]]; then
  echo 'External Miniforge is not supported for the independent workflow.' >&2; exit 1
fi
unset PYTHONPATH PYTHONHOME
export HOME="$PWD/.runtime-home"
export XDG_CACHE_HOME="$PWD/runtime/cache" XDG_CONFIG_HOME="$PWD/runtime/config" XDG_DATA_HOME="$PWD/runtime/share"
export TMPDIR="$PWD/runtime/tmp" MPLCONFIGDIR="$PWD/runtime/matplotlib" PYTHONNOUSERSITE=1
export CONDA_PKGS_DIRS="$MINIFORGE/pkgs" CONDA_ENVS_PATH="$MINIFORGE/envs"
mkdir -p "$HOME" "$XDG_CACHE_HOME" "$XDG_CONFIG_HOME" "$XDG_DATA_HOME" "$TMPDIR" "$MPLCONFIGDIR"
case "${1:-help}" in
  inputs) shift; python3 -m karakoram.setup_inputs "$@" ;;
  rebuild-reference|rebuild-source)
    mode="$1"
    # Validate originals before expensive network or computation.
    python3 -m karakoram.setup_inputs verify-observations
    python3 -m karakoram.setup_inputs bootstrap
    python3 -m karakoram.setup_inputs acquire-static
    # An explicit request to replay the documented seven omissions and full calibration.
    bash workflow.sh basin --fetch-rgi
    bash workflow.sh observations
    if ! find data/climate/era5/raw/ecmwf_arco -name COMPLETE.json -print -quit 2>/dev/null | grep -q .; then
      bash workflow.sh acquire-climate
    fi
    bash workflow.sh prepare-climate
    bash workflow.sh glaciers --fetch-prepro --prune-gdirs --reproduce-omissions
    bash workflow.sh export --accept-omissions
    if [[ "$mode" == rebuild-reference ]]; then
      bash workflow.sh historical --repeat-calibration --require-reference-match
    else
      # Retain comparison report, but distinguish independently reacquired inputs
      # from an exact-byte reference replay. Calibration selection is still checked.
      bash workflow.sh historical --repeat-calibration
    fi
    ;;
  basin)
    shift
    # Acquire RGI inputs only if requested and missing; computation uses local files.
    if [[ "$#" -eq 1 && "$1" == --fetch-rgi ]]; then
      if [[ ! -d data/processed/basin ]]; then
        bash run_oggm.sh -m karakoram.setup_inputs prepare-oggm basin --download
      fi
      bash run_oggm.sh -m karakoram.basin
    else
      bash run_oggm.sh -m karakoram.basin "$@"
    fi
    ;;
  observations) shift; "$MINIFORGE/envs/cwatm_env/bin/python" -m karakoram.observations "$@" ;;
  acquire-climate) shift; "$MINIFORGE/envs/cwatm_env/bin/python" -m karakoram.acquire_climate "$@" ;;
  prepare-climate) shift; "$MINIFORGE/envs/cwatm_env/bin/python" -m karakoram.prepare_climate "$@" ;;
  glaciers)
    shift
    acquire=(); omit=(); self_test=0
    for arg in "$@"; do
      case "$arg" in
        --fetch-prepro) acquire=(--download) ;;
        --reproduce-omissions) omit=(--reproduce-omissions) ;;
        --self-test|-h|--help) self_test=1 ;;
      esac
    done
    if [[ "$self_test" -eq 0 && "${#acquire[@]}" -ne 0 ]]; then
      bash run_oggm.sh -m karakoram.setup_inputs prepare-oggm glaciers "${acquire[@]}" "${omit[@]}"
    fi
    bash run_oggm.sh -m karakoram.glaciers "$@"
    ;;
  export) shift; bash run_oggm.sh -m karakoram.coupling "$@" ;;
  historical) shift; "$MINIFORGE/envs/cwatm_env/bin/python" -m karakoram.hydrology "$@" ;;
  test) shift; bash run_oggm.sh -m unittest discover -s tests "$@" ;;
  *) printf '%s\n' 'Usage: bash workflow.sh {inputs|rebuild-reference|rebuild-source|basin|observations|acquire-climate|prepare-climate|glaciers|export|historical|test} [options]' ;;
esac
