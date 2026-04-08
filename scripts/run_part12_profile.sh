#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

cfg=${1:-inputs/part10_n5000.in}
if [[ ! -f "${cfg}" && -f "inputs/${cfg}" ]]; then
    cfg="inputs/${cfg}"
fi
case_name=$(basename "${cfg}" .in)

mkdir -p part12_profile

make -j

export OMP_NUM_THREADS=1
rm -f output_*.vtk

log_file="part12_profile/${case_name}_serial.log"
./demo.exe "${cfg}" > "${log_file}" 2>&1

get_profile_value() {
    local key="$1"
    awk -v k="${key}" '$1=="PROFILE" && $2==k {print $3}' "${log_file}" | tail -1
}

total=$(get_profile_value total_runtime_s)
t_init=$(get_profile_value initialize_s)
t_grav=$(get_profile_value gravity_s)
t_pp=$(get_profile_value particle_contacts_s)
t_walls=$(get_profile_value wall_contacts_s)
t_int=$(get_profile_value integration_s)
t_out=$(get_profile_value output_s)

pct_init=$(get_profile_value pct_initialize)
pct_grav=$(get_profile_value pct_gravity)
pct_pp=$(get_profile_value pct_particle_contacts)
pct_walls=$(get_profile_value pct_wall_contacts)
pct_int=$(get_profile_value pct_integration)
pct_out=$(get_profile_value pct_output)

bottleneck=$(awk '
$1=="PROFILE" && $2 ~ /_s$/ && $2 != "total_runtime_s" {
    if ($3+0 > max) {max=$3; name=$2}
}
END {print name}
' "${log_file}")

cat > part12_serial_profile.csv <<EOF
case,total_runtime_s,initialize_s,gravity_s,particle_contacts_s,wall_contacts_s,integration_s,output_s,pct_initialize,pct_gravity,pct_particle_contacts,pct_wall_contacts,pct_integration,pct_output,bottleneck
${case_name},${total},${t_init},${t_grav},${t_pp},${t_walls},${t_int},${t_out},${pct_init},${pct_grav},${pct_pp},${pct_walls},${pct_int},${pct_out},${bottleneck}
EOF

cat > part12_profile_report.txt <<EOF
Part 12 Serial Profiling Report
Case: ${case_name}
Log file: ${log_file}

Total runtime [s]: ${total}

Runtime distribution [s]:
- initialize_s: ${t_init}
- gravity_s: ${t_grav}
- particle_contacts_s: ${t_pp}
- wall_contacts_s: ${t_walls}
- integration_s: ${t_int}
- output_s: ${t_out}

Runtime distribution [%]:
- pct_initialize: ${pct_init}
- pct_gravity: ${pct_grav}
- pct_particle_contacts: ${pct_pp}
- pct_wall_contacts: ${pct_walls}
- pct_integration: ${pct_int}
- pct_output: ${pct_out}

Identified bottleneck: ${bottleneck}
EOF

echo "Saved serial profile log: ${log_file}"
echo "Saved serial profile CSV: part12_serial_profile.csv"
echo "Saved serial profile report: part12_profile_report.txt"
