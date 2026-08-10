#!/usr/bin/env bash
set -euo pipefail

# Shell-only entry point for the Ka/Ks reviewer sensitivity check.
# Optional: pass a conda environment name as the first argument.
# Example:
#   bash MS_REVISION/DNDS_SENSITIVITY/scripts/run_dnds_sensitivity.sh
#   bash MS_REVISION/DNDS_SENSITIVITY/scripts/run_dnds_sensitivity.sh dnds_sensitivity

CONDA_ENV="${1:-}"

if [[ -n "${CONDA_ENV}" ]]; then
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "${CONDA_ENV}"
  elif command -v mamba >/dev/null 2>&1; then
    eval "$(mamba shell.bash hook)"
    mamba activate "${CONDA_ENV}"
  else
    echo "ERROR: CONDA_ENV was supplied, but neither conda nor mamba is available on PATH." >&2
    exit 1
  fi
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REV_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_DIR="$(cd "${REV_DIR}/../.." && pwd)"
RESULTS_DIR="${REV_DIR}/results"
PLOTS_DIR="${REV_DIR}/plots"
FONT_CACHE_DIR="${TMPDIR:-/tmp}/dnds_sensitivity_fontcache_$(id -u)"

mkdir -p "${RESULTS_DIR}" "${PLOTS_DIR}" "${FONT_CACHE_DIR}"
export XDG_CACHE_HOME="${FONT_CACHE_DIR}"

python3 "${SCRIPT_DIR}/build_dnds_sensitivity.py" \
  --kaks-dir "${REPO_DIR}/Results/DNDS/PAL2NAL_DEDUP_KAKS" \
  --align-dir "${REPO_DIR}/Results/DNDS/PAL2NAL_DEDUP" \
  --classification "${REPO_DIR}/Manuscript_Figures/Figure_2/Data/twilight_classification_updated.tab" \
  --outdir "${RESULTS_DIR}"

Rscript "${SCRIPT_DIR}/plot_dnds_sensitivity.R" "${RESULTS_DIR}" "${PLOTS_DIR}"

echo "Sensitivity tables: ${RESULTS_DIR}"
echo "Sensitivity plots:  ${PLOTS_DIR}"
