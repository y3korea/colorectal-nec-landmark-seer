#!/usr/bin/env bash
# One-command run of the full pipeline (see run_all.py --help).
#   ./run_all.sh            full run with the paper's settings (B = 2000)
#   ./run_all.sh --quick    smoke test (small bootstrap; e.g. on synthetic data from tools/make_synthetic_seer_export.py)
# Configure with NEC_DATA_DIR, NEC_WORK_DIR, NEC_RESULTS_DIR, NEC_FIG_DIR, NEC_BOOT (see README.md).
# Relative NEC_*_DIR values are resolved against the directory you start from (run_all.py passes absolute paths on).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-python3}"
exec "$PY" "$HERE/run_all.py" "$@"
