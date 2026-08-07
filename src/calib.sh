#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${ROOT_DIR}/two-steps-calibration"

CHAINS="${CHAINS:-4}"
CORES="${CORES:-4}"
WARMUP="${WARMUP:-500}"
ITER="${ITER:-600}"
PPC="${PPC:-600}"

mkdir -p "${OUT_DIR}"
cd "${ROOT_DIR}"

echo "Step 1/2: Nelder-Mead MAP calibration"
Rscript nelder_mead.R --out-dir="${OUT_DIR}" 2>&1 | tee "${OUT_DIR}/nelder_mead.log"

echo "Step 2/2: HMC calibration using Nelder-Mead-derived priors"
Rscript hmc.R \
  --nm-fit="${OUT_DIR}/nelder_mead_fit.rds" \
  --out-dir="${OUT_DIR}" \
  --chains="${CHAINS}" \
  --cores="${CORES}" \
  --warmup="${WARMUP}" \
  --iter="${ITER}" \
  --ppc="${PPC}" \
  2>&1 | tee "${OUT_DIR}/hmc.log"

echo "Complete. Outputs are in ${OUT_DIR}"
