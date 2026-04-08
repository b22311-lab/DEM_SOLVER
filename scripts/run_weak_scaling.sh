#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

threads_list=${THREADS_LIST:-"1 6 8"}
repeats=${REPEATS:-5}
base_n=${BASE_N_PER_THREAD:-500}
base_len=${BASE_DOMAIN_LENGTH:-21.2132034356}

mkdir -p weak_scaling_logs weak_scaling_configs
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

run_csv="weak_scaling_raw_runs.csv"
raw_csv="weak_scaling_raw.csv"

echo "case,threads,nparticles,domain_x,domain_y,repeat,total_runtime_s,particle_contacts_s,contact_candidates_total,contacts_detected_total" > "${run_csv}"

for t in ${threads_list}; do
    n=$((base_n * t))
    len=$(awk -v l="${base_len}" -v p="${t}" 'BEGIN{printf "%.8f", l*sqrt(p)}')

    cfg="weak_scaling_configs/weak_n${n}_t${t}.in"
    cat > "${cfg}" <<EOF
! Weak scaling auto-generated case
&domain
    dimension=2
    xr_in=${len}
    yr_in=${len}
/

&material_properties
    k=2000
    gamma=1.0
    density=1000.0
    rad_init=0.05
    nparticles=${n}
/

&particle_initialization
/

&simulation_control
    time_step=0.001
    end_time=0.02
    write_interval=1000000000
    verbose=.false.
    contact_search_method=0
    cell_size_factor=2.5
/

&gravity
    gravity_axis=2
    gravity_magnitude=9.81
/
/
EOF

    export OMP_NUM_THREADS=${t}

    for r in $(seq 1 "${repeats}"); do
        rm -f output_*.vtk

        log_file="weak_scaling_logs/weak_n${n}_t${t}_r${r}.log"
        if ! ./demo.exe "${cfg}" > "${log_file}" 2>&1; then
            echo "Weak scaling run failed: t=${t}, repeat=${r}."
            echo "weak_scaling,${t},${n},${len},${len},${r},nan,nan,nan,nan" >> "${run_csv}"
            continue
        fi

        total=$(awk '$1=="PROFILE" && $2=="total_runtime_s" {print $3}' "${log_file}" | tail -1)
        t_pp=$(awk '$1=="PROFILE" && $2=="particle_contacts_s" {print $3}' "${log_file}" | tail -1)
        cand=$(awk '$1=="PROFILE" && $2=="contact_candidates_total" {print $3}' "${log_file}" | tail -1)
        cont=$(awk '$1=="PROFILE" && $2=="contacts_detected_total" {print $3}' "${log_file}" | tail -1)

        echo "weak_scaling,${t},${n},${len},${len},${r},${total},${t_pp},${cand},${cont}" >> "${run_csv}"
        echo "Weak scaling: threads=${t}, repeat=${r}, N=${n}, runtime=${total}s"
    done
done

run_python_inline <<'PY'
import csv
import math
import statistics
from collections import defaultdict

run_csv = "weak_scaling_raw_runs.csv"
raw_csv = "weak_scaling_raw.csv"

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
        meta[key] = {
            "case": row["case"],
            "nparticles": row["nparticles"],
            "domain_x": row["domain_x"],
            "domain_y": row["domain_y"],
        }
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
        fieldnames=["case", "threads", "nparticles", "domain_x", "domain_y", *metrics],
    )
    writer.writeheader()
    for t in order:
        out = {
            "case": meta[t]["case"],
            "threads": t,
            "nparticles": meta[t]["nparticles"],
            "domain_x": meta[t]["domain_x"],
            "domain_y": meta[t]["domain_y"],
        }
        for m in metrics:
            vals = b[t][m]
            out[m] = statistics.median(vals) if vals else float("nan")
        writer.writerow(out)
PY

"${PYTHON_RUN[@]}" python/analyze_weak_scaling.py "${raw_csv}"

echo "Benchmark settings: threads=${threads_list}, repeats=${repeats}, OMP_PROC_BIND=${OMP_PROC_BIND}, OMP_PLACES=${OMP_PLACES}"
echo "Timing runs use write_interval=1000000000 to minimize output overhead."
echo "Saved: ${run_csv}"
echo "Saved: ${raw_csv}"
echo "Saved: weak_scaling_summary.csv"
echo "Saved: weak_scaling_runtime.png, weak_scaling_efficiency.png"
echo "Saved: weak_scaling_findings.tex"
