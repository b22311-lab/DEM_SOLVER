#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

threads_list=${THREADS_LIST:-"1 6 8"}
repeats=${REPEATS:-5}
cases=(inputs/part10_n200.in inputs/part10_n1000.in inputs/part10_n5000.in)

mkdir -p perf_logs
make -j

export OMP_PROC_BIND="${OMP_PROC_BIND:-true}"
export OMP_PLACES="${OMP_PLACES:-cores}"

if command -v conda >/dev/null 2>&1 && conda run -n "${CONDA_ENV_NAME:-dem-solver}" python -c "import sys" >/dev/null 2>&1; then
    PYTHON_RUN=(conda run -n "${CONDA_ENV_NAME:-dem-solver}" python)
else
    PYTHON_RUN=("${PYTHON:-python3}")
fi

run_python_inline() {
    local tmp_py
    tmp_py=$(mktemp)
    cat > "${tmp_py}"
    "${PYTHON_RUN[@]}" "${tmp_py}"
    rm -f "${tmp_py}"
}

tmp_cfg_dir=$(mktemp -d)
cleanup() {
    rm -rf "${tmp_cfg_dir}"
}
trap cleanup EXIT

make_timing_cfg() {
    local src="$1"
    local dst="${tmp_cfg_dir}/$(basename "${src}")"
    awk '
        /^[[:space:]]*write_interval[[:space:]]*=/ {
            print "    write_interval=1000000000"
            next
        }
        { print }
    ' "${src}" > "${dst}"
    echo "${dst}"
}

echo "case,threads,repeat,total_runtime_s,initialize_s,gravity_s,particle_contacts_s,wall_contacts_s,integration_s,output_s" > performance_raw_runs.csv

for cfg in "${cases[@]}"; do
    case_name=$(basename "${cfg}" .in)
    cfg_timing=$(make_timing_cfg "${cfg}")

    for t in ${threads_list}; do
        export OMP_NUM_THREADS=${t}

        for r in $(seq 1 "${repeats}"); do
            rm -f output_*.vtk

            log_file="perf_logs/${case_name}_t${t}_r${r}.log"
            if ! ./demo.exe "${cfg_timing}" > "${log_file}" 2>&1; then
                echo "Run failed for ${case_name}, t=${t}, repeat=${r}. Check ${log_file}."
                echo "${case_name},${t},${r},nan,nan,nan,nan,nan,nan,nan" >> performance_raw_runs.csv
                continue
            fi

            total=$(awk '$1=="PROFILE" && $2=="total_runtime_s" {print $3}' "${log_file}" | tail -1)
            t_init=$(awk '$1=="PROFILE" && $2=="initialize_s" {print $3}' "${log_file}" | tail -1)
            t_grav=$(awk '$1=="PROFILE" && $2=="gravity_s" {print $3}' "${log_file}" | tail -1)
            t_pp=$(awk '$1=="PROFILE" && $2=="particle_contacts_s" {print $3}' "${log_file}" | tail -1)
            t_walls=$(awk '$1=="PROFILE" && $2=="wall_contacts_s" {print $3}' "${log_file}" | tail -1)
            t_int=$(awk '$1=="PROFILE" && $2=="integration_s" {print $3}' "${log_file}" | tail -1)
            t_out=$(awk '$1=="PROFILE" && $2=="output_s" {print $3}' "${log_file}" | tail -1)

            if [[ -z "${total}" ]]; then
                echo "Missing profile output for ${case_name}, t=${t}, repeat=${r}."
                echo "${case_name},${t},${r},nan,nan,nan,nan,nan,nan,nan" >> performance_raw_runs.csv
                continue
            fi

            echo "${case_name},${t},${r},${total},${t_init},${t_grav},${t_pp},${t_walls},${t_int},${t_out}" >> performance_raw_runs.csv
            echo "Finished ${case_name}, t=${t}, repeat=${r}: total_runtime=${total} s"
        done
    done
done

run_python_inline <<'PY'
import csv
import math
import statistics
from collections import defaultdict

metrics = [
    "total_runtime_s",
    "initialize_s",
    "gravity_s",
    "particle_contacts_s",
    "wall_contacts_s",
    "integration_s",
    "output_s",
]

bucket = defaultdict(lambda: {m: [] for m in metrics})
order = []

with open("performance_raw_runs.csv", "r", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = (row["case"], int(row["threads"]))
        if key not in bucket:
            order.append(key)
        for m in metrics:
            try:
                v = float(row[m])
            except (TypeError, ValueError):
                continue
            if math.isfinite(v):
                bucket[key][m].append(v)

with open("performance_raw.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["case", "threads", *metrics])
    writer.writeheader()
    for case, threads in order:
        out = {"case": case, "threads": threads}
        for m in metrics:
            vals = bucket[(case, threads)][m]
            out[m] = statistics.median(vals) if vals else float("nan")
        writer.writerow(out)
PY

"${PYTHON_RUN[@]}" python/analyze_speedup.py performance_raw.csv

echo "Benchmark settings: threads=${threads_list}, repeats=${repeats}, OMP_PROC_BIND=${OMP_PROC_BIND}, OMP_PLACES=${OMP_PLACES}"
echo "Timing runs use write_interval=1000000000 to minimize output overhead."
echo "Saved: performance_raw_runs.csv"
echo "Saved: performance_raw.csv"
echo "Saved: performance_speedup.csv"
echo "Saved: performance_amdahl_fit.csv"
echo "Saved: speedup_plot.png, efficiency_plot.png"
