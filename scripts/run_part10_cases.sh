#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

cases=(inputs/part10_n200.in inputs/part10_n1000.in inputs/part10_n5000.in)

make -j
mkdir -p part10_outputs

for cfg in "${cases[@]}"; do
    case_name=$(basename "${cfg}" .in)
    out_dir="part10_outputs/${case_name}"

    rm -f output_*.vtk
    mkdir -p "${out_dir}"

    OMP_NUM_THREADS=1 ./demo.exe "${cfg}" > "${out_dir}/run.log" 2>&1
    mv output_*.vtk "${out_dir}/"

    echo "Saved ${case_name} outputs to ${out_dir}"
    awk '$1=="PROFILE" {print $0}' "${out_dir}/run.log"
done
