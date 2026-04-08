#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

make -j

if command -v conda >/dev/null 2>&1 && conda run -n "${CONDA_ENV_NAME:-dem-solver}" python -c "import sys" >/dev/null 2>&1; then
	PYTHON_RUN=(conda run -n "${CONDA_ENV_NAME:-dem-solver}" python)
else
	PYTHON_RUN=("${PYTHON:-python3}")
fi

"${PYTHON_RUN[@]}" python/run_part19_2_study.py

echo "Saved: part19_2_timestep_sensitivity.csv"
echo "Saved: part19_2_trajectory_error_overlay.png"
echo "Saved: part19_2_error_vs_dt.png"
echo "Saved: part19_2_energy_evolution.png"
echo "Saved: part19_2_energy_drift_vs_dt.png"
echo "Saved: part19_2_richardson.csv"
echo "Saved: part19_2_critical_stability.csv"
echo "Saved: part19_2_critical_energy.png"
echo "Saved: part19_2_findings.tex"
