#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

echo "[1/7] Verification figures"
scripts/run_all_tests.sh

echo "[2/7] OpenMP performance sweep (Part 10-14)"
scripts/run_part10_14.sh

echo "[3/7] Serial hotspot profile (Part 12)"
scripts/run_part12_profile.sh

echo "[4/7] Strong/weak scaling"
scripts/run_scaling_modes.sh

echo "[5/7] Neighbor-search study (Part 18)"
scripts/run_part18_study.sh

echo "[6/7] Timestep and stability study (Part 19.2)"
scripts/run_part19_2_study.sh

echo "[7/7] LaTeX compile checks"
if [[ -x scripts/compile_latex_all.sh ]]; then
	ALLOW_MISSING_PDFLATEX=1 scripts/compile_latex_all.sh
else
	echo "Skipped LaTeX compile: scripts/compile_latex_all.sh not found."
fi

echo "All requested artifacts generated."
