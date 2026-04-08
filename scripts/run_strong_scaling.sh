#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

threads_list=${THREADS_LIST:-"1 6 8"}
repeats=${REPEATS:-5}
cfg=${STRONG_CONFIG:-"inputs/part10_n30000.in"}
if [[ ! -f "${cfg}" && -f "inputs/${cfg}" ]]; then
    cfg="inputs/${cfg}"
fi

mkdir -p strong_scaling_logs
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

cfg_timing="${tmp_cfg_dir}/$(basename "${cfg}")"
awk '
    /^[[:space:]]*write_interval[[:space:]]*=/ {
        print "    write_interval=1000000000"
        next
    }
    { print }
' "${cfg}" > "${cfg_timing}"

nparticles=$(awk -F= '/nparticles/ {gsub(/[ ,]/, "", $2); print $2; exit}' "${cfg_timing}")

run_csv="strong_scaling_raw_runs.csv"
raw_csv="strong_scaling_raw.csv"
echo "case,nparticles,threads,repeat,total_runtime_s,particle_contacts_s,contact_candidates_total,contacts_detected_total" > "${run_csv}"

for t in ${threads_list}; do
    export OMP_NUM_THREADS=${t}

    for r in $(seq 1 "${repeats}"); do
        rm -f output_*.vtk

        log_file="strong_scaling_logs/$(basename "${cfg}" .in)_t${t}_r${r}.log"
        if ! ./demo.exe "${cfg_timing}" > "${log_file}" 2>&1; then
            echo "Strong scaling run failed: t=${t}, repeat=${r}."
            echo "$(basename "${cfg}" .in),${nparticles},${t},${r},nan,nan,nan,nan" >> "${run_csv}"
            continue
        fi

        total=$(awk '$1=="PROFILE" && $2=="total_runtime_s" {print $3}' "${log_file}" | tail -1)
        t_pp=$(awk '$1=="PROFILE" && $2=="particle_contacts_s" {print $3}' "${log_file}" | tail -1)
        cand=$(awk '$1=="PROFILE" && $2=="contact_candidates_total" {print $3}' "${log_file}" | tail -1)
        cont=$(awk '$1=="PROFILE" && $2=="contacts_detected_total" {print $3}' "${log_file}" | tail -1)

        echo "$(basename "${cfg}" .in),${nparticles},${t},${r},${total},${t_pp},${cand},${cont}" >> "${run_csv}"
        echo "Strong scaling: threads=${t}, repeat=${r}, total=${total}s"
    done
done

run_python_inline <<'PY'
import csv
import math
import statistics
from collections import defaultdict

run_csv = "strong_scaling_raw_runs.csv"
raw_csv = "strong_scaling_raw.csv"

metrics = [
    "total_runtime_s",
    "particle_contacts_s",
    "contact_candidates_total",
    "contacts_detected_total",
]

b = defaultdict(lambda: {m: [] for m in metrics})
meta = {}
order = []

with open(run_csv, "r", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = int(row["threads"])
        if key not in b:
            order.append(key)
        meta[key] = {"case": row["case"], "nparticles": row["nparticles"]}
        for m in metrics:
            try:
                v = float(row[m])
            except (TypeError, ValueError):
                continue
            if math.isfinite(v):
                b[key][m].append(v)

with open(raw_csv, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["case", "nparticles", "threads", *metrics],
    )
    writer.writeheader()
    for t in order:
        out = {
            "case": meta[t]["case"],
            "nparticles": meta[t]["nparticles"],
            "threads": t,
        }
        for m in metrics:
            vals = b[t][m]
            out[m] = statistics.median(vals) if vals else float("nan")
        writer.writerow(out)
PY

"${PYTHON_RUN[@]}" python/analyze_strong_scaling.py "${raw_csv}"

echo "Benchmark settings: threads=${threads_list}, repeats=${repeats}, OMP_PROC_BIND=${OMP_PROC_BIND}, OMP_PLACES=${OMP_PLACES}"
echo "Timing runs use write_interval=1000000000 to minimize output overhead."
echo "Saved: ${run_csv}"
echo "Saved: ${raw_csv}"
echo "Saved: strong_scaling_summary.csv"
echo "Saved: strong_scaling_speedup.png, strong_scaling_efficiency.png"
echo "Saved: strong_scaling_findings.tex"
