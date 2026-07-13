#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${ROOT_DIR}/two-steps-calibration"
NM_FIT="${OUT_DIR}/nelder_mead_fit.rds"

mkdir -p "${OUT_DIR}"
cd "${ROOT_DIR}"

echo "Step 1/2: Nelder-Mead MAP calibration"
Rscript nelder_mead.R "${OUT_DIR}"

echo "Step 2/2: HMC calibration initialized from Nelder-Mead estimates"
Rscript hmc.R "${OUT_DIR}" "${NM_FIT}"

echo "Two-step calibration complete. Outputs are in ${OUT_DIR}"
