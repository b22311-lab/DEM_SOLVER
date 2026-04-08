#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

cfg=${1:-inputs/part10_n200.in}
if [[ ! -f "${cfg}" && -f "inputs/${cfg}" ]]; then
  cfg="inputs/${cfg}"
fi
threads=${2:-2}
case_name=$(basename "${cfg}" .in)

serial_dir="parallel_validation_outputs/${case_name}/serial"
parallel_dir="parallel_validation_outputs/${case_name}/omp${threads}"
mkdir -p "${serial_dir}" "${parallel_dir}"

make -j

rm -f output_*.vtk
OMP_NUM_THREADS=1 ./demo.exe "${cfg}" > "${serial_dir}/run.log" 2>&1
mv output_*.vtk "${serial_dir}/"

rm -f output_*.vtk
OMP_NUM_THREADS=${threads} ./demo.exe "${cfg}" > "${parallel_dir}/run.log" 2>&1
mv output_*.vtk "${parallel_dir}/"

if command -v conda >/dev/null 2>&1 && conda run -n "${CONDA_ENV_NAME:-dem-solver}" python -c "import sys" >/dev/null 2>&1; then
  PYTHON_RUN=(conda run -n "${CONDA_ENV_NAME:-dem-solver}" python)
else
  PYTHON_RUN=("${PYTHON:-python3}")
fi

case_report="parallel_validation_report_${case_name}_t${threads}.csv"
case_log="parallel_validation_${case_name}_t${threads}.log"

"${PYTHON_RUN[@]}" python/compare_vtk_outputs.py \
  --serial-dir "${serial_dir}" \
  --parallel-dir "${parallel_dir}" \
  --case "${case_name}" \
  --threads "${threads}" \
  --output "${case_report}" \
  > "${case_log}"

if [[ ! -f parallel_validation_report.csv ]]; then
    cp "${case_report}" parallel_validation_report.csv
else
    tail -n +2 "${case_report}" >> parallel_validation_report.csv
fi

echo "Saved serial outputs: ${serial_dir}"
echo "Saved parallel outputs: ${parallel_dir}"
echo "Saved validation summary: ${case_report}"
echo "Saved validation log: ${case_log}"
echo "Updated aggregate report: parallel_validation_report.csv"
