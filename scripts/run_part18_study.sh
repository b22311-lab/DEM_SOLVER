#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

configs=(
  inputs/part18_n200_allpairs.in
  inputs/part18_n200_cell.in
  inputs/part18_n1000_allpairs.in
  inputs/part18_n1000_cell.in
  inputs/part18_n5000_allpairs.in
  inputs/part18_n5000_cell.in
)

repeats=${REPEATS:-5}
part18_threads=${PART18_THREADS:-1}

make -j
mkdir -p part18_logs part18_outputs

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

run_csv="part18_neighbor_search_raw_runs.csv"
raw_csv="part18_neighbor_search_raw.csv"

echo "case,nparticles,method_id,method,repeat,total_runtime_s,particle_contacts_s,contact_candidates_total,contacts_detected_total,loop_count" > "${run_csv}"

for cfg in "${configs[@]}"; do
  case_name=$(basename "${cfg}" .in)
  cfg_timing=$(make_timing_cfg "${cfg}")

  nparticles=$(awk -F= '/nparticles/ {gsub(/[ ,]/, "", $2); print $2; exit}' "${cfg_timing}")
  loop_count=$(awk -F= '
/time_step/ {gsub(/[ ,]/, "", $2); dt=$2}
/end_time/ {gsub(/[ ,]/, "", $2); et=$2}
END {if (dt+0 > 0) printf "%d", et/dt; else print 0}
' "${cfg_timing}")

  for r in $(seq 1 "${repeats}"); do
    log_file="part18_logs/${case_name}_r${r}.log"

    rm -f output_*.vtk
    OMP_NUM_THREADS=${part18_threads} ./demo.exe "${cfg_timing}" > "${log_file}" 2>&1

    method_id=$(awk '$1=="PROFILE" && $2=="contact_search_method" {print int($3); exit}' "${log_file}")
    total=$(awk '$1=="PROFILE" && $2=="total_runtime_s" {print $3}' "${log_file}" | tail -1)
    t_pp=$(awk '$1=="PROFILE" && $2=="particle_contacts_s" {print $3}' "${log_file}" | tail -1)
    cand=$(awk '$1=="PROFILE" && $2=="contact_candidates_total" {print $3}' "${log_file}" | tail -1)
    cont=$(awk '$1=="PROFILE" && $2=="contacts_detected_total" {print $3}' "${log_file}" | tail -1)

    method="all_pairs"
    if [[ "${method_id}" == "1" ]]; then
      method="cell_linked"
    fi

    echo "${case_name},${nparticles},${method_id},${method},${r},${total},${t_pp},${cand},${cont},${loop_count}" >> "${run_csv}"
    echo "Finished ${case_name}: repeat=${r}, method=${method}, runtime=${total}s"
  done
done

run_python_inline <<'PY'
import csv
import math
import statistics
from collections import defaultdict

run_csv = "part18_neighbor_search_raw_runs.csv"
raw_csv = "part18_neighbor_search_raw.csv"

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
        key = row["case"]
        if key not in b:
            order.append(key)
        meta[key] = {
            "nparticles": row["nparticles"],
            "method_id": row["method_id"],
            "method": row["method"],
            "loop_count": row["loop_count"],
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
        fieldnames=[
            "case",
            "nparticles",
            "method_id",
            "method",
            "total_runtime_s",
            "particle_contacts_s",
            "contact_candidates_total",
            "contacts_detected_total",
            "loop_count",
        ],
    )
    writer.writeheader()
    for case in order:
        out = {
            "case": case,
            "nparticles": meta[case]["nparticles"],
            "method_id": meta[case]["method_id"],
            "method": meta[case]["method"],
            "loop_count": meta[case]["loop_count"],
        }
        for m in metrics:
            vals = b[case][m]
            out[m] = statistics.median(vals) if vals else float("nan")
        writer.writerow(out)
PY

"${PYTHON_RUN[@]}" python/analyze_part18.py "${raw_csv}"

echo "Benchmark settings: repeats=${repeats}, OMP_NUM_THREADS=${part18_threads}"
echo "Timing runs use write_interval=1000000000 to minimize output overhead."
echo "Saved: ${run_csv}"
echo "Saved: ${raw_csv}"
echo "Saved: part18_neighbor_search_summary.csv"
echo "Saved: part18_contact_runtime.png, part18_candidates.png, part18_speedup.png"
echo "Saved: part18_findings.tex"
