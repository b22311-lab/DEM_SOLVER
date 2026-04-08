#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

MAIN_TEX=${MAIN_TEX:-ieee_dem_whitepaper.tex}
BUILD_DIR=${LATEX_BUILD_DIR:-latex_build}
BASE_NAME="${MAIN_TEX%.tex}"

if [[ ! -f "${MAIN_TEX}" ]]; then
    echo "Missing LaTeX source: ${MAIN_TEX}"
    exit 1
fi

if ! command -v pdflatex >/dev/null 2>&1; then
    echo "pdflatex not found. Install a LaTeX distribution with IEEEtran support."
    exit 1
fi

mkdir -p "${BUILD_DIR}"

pdflatex -interaction=nonstopmode -halt-on-error -output-directory "${BUILD_DIR}" "${MAIN_TEX}"
pdflatex -interaction=nonstopmode -halt-on-error -output-directory "${BUILD_DIR}" "${MAIN_TEX}"

cp "${BUILD_DIR}/${BASE_NAME}.pdf" "${BASE_NAME}.pdf"
if [[ -f "${BUILD_DIR}/${BASE_NAME}.aux" ]]; then
    cp "${BUILD_DIR}/${BASE_NAME}.aux" "${BASE_NAME}.aux"
fi
if [[ -f "${BUILD_DIR}/${BASE_NAME}.log" ]]; then
    cp "${BUILD_DIR}/${BASE_NAME}.log" "${BASE_NAME}.log"
fi

echo "Saved: ${BASE_NAME}.pdf"
echo "Saved: ${BASE_NAME}.aux"
echo "Saved: ${BASE_NAME}.log"
