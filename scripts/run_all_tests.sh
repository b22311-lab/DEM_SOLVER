#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

run_case() {
	local cfg="$1"
	local test_name="$2"
	local figure_name="$3"

	rm -f output_*.vtk
	./demo.exe "${cfg}" >/dev/null 2>&1
	python3 python/test_verification.py --test "${test_name}" --config "${cfg}" >/dev/null 2>&1
	echo "Generated: ${figure_name}"
}

run_case "inputs/test1_freefall_config.in" "test1" "test1_freefall_verification.png"
run_case "inputs/test2_constant_velocity_config.in" "test2" "test2_constant_velocity_verification.png"
run_case "inputs/test3_particle_bounce_config_stable3.in" "test3" "test3_particle_bounce_verification.png"

echo "Verification image generation complete."