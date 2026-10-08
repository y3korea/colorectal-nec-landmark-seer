# -*- coding: utf-8 -*-
"""Shared path / run configuration for every Python step of the pipeline.

All locations are configured through environment variables; relative values are resolved against the
current working directory (run_all.py resolves them once, against the directory it is started from, and passes
absolute paths to every step); unset variables fall back to the defaults below (relative to the repository root).

    NEC_DATA_DIR     folder holding the four SEER*Stat case-listing exports            default ./data/raw
    NEC_WORK_DIR     derived cohort files (cohort_*.csv, flow_all.csv; patient-level!) default ./work
    NEC_RESULTS_DIR  aggregate results (JSON / Markdown / CSV)                        default ./results
    NEC_FIG_DIR      figures (analysis figures; manuscript figures in manuscript_final/) default ./figures/output
    NEC_BOOT         bootstrap replicates of the main analyses                        default 2000 (as in the paper);
                     secondary bootstraps (500 / 1000 in the paper) are scaled proportionally
    NEC_QUICK        1 = smoke-test mode (NEC_BOOT defaults to 40, anchors are informational)
    NEC_STRICT       1 = stop when a reproduction anchor of the paper is not met (default: 1, except in
                     QUICK mode or on synthetic data, where it defaults to 0)
    NEC_FILE_NEN / NEC_FILE_AC / NEC_FILE_LUNG / NEC_FILE_PANC
                     file names of the four exports inside NEC_DATA_DIR (defaults: the authors' names)
"""
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _truthy(v):
    return str(v).strip().lower() in ("1", "true", "yes", "y", "on")


def _dir(env, default):
    v = os.environ.get(env, "").strip()
    p = Path(v).expanduser() if v else REPO / default
    return str(p.resolve())


DATA_DIR = _dir("NEC_DATA_DIR", "data/raw")
WORK_DIR = _dir("NEC_WORK_DIR", "work")
RESULTS_DIR = _dir("NEC_RESULTS_DIR", "results")
FIG_DIR = _dir("NEC_FIG_DIR", "figures/output")
MANUSCRIPT_FIG_DIR = os.path.join(FIG_DIR, "manuscript_final")
PY_DIR = str(REPO / "python")
R_DIR = str(REPO / "R")

# trailing-separator variants used by the analysis scripts (string concatenation, as in the original code)
DATA = WORK_DIR + os.sep          # derived cohort files
RES = RESULTS_DIR + os.sep
FIGS = FIG_DIR + os.sep
RAW = DATA_DIR + os.sep

RAW_FILES = {
    "NEN": os.environ.get("NEC_FILE_NEN", "colorectal_NEN_SEER17_2000_2023_full_20260806.txt"),
    "AC": os.environ.get("NEC_FILE_AC", "colorectal_adenocarcinoma_SEER17_2000_2023_full_20260928.txt"),
    "LUNG": os.environ.get("NEC_FILE_LUNG", "SCLC_lung_SEER17_2000_2023_full_20260928.txt"),
    "PANC": os.environ.get("NEC_FILE_PANC", "pancreatic_adenocarcinoma_SEER17_2000_2023_full_20260928.txt"),
}

SYNTHETIC_MARKER = "SYNTHETIC_DATA_README.txt"
SYNTHETIC = os.path.exists(os.path.join(DATA_DIR, SYNTHETIC_MARKER))
QUICK = _truthy(os.environ.get("NEC_QUICK", "0"))
PAPER_B = 2000
BOOT = int(os.environ["NEC_BOOT"]) if os.environ.get("NEC_BOOT", "").strip() else (40 if QUICK else PAPER_B)
STRICT = _truthy(os.environ["NEC_STRICT"]) if os.environ.get("NEC_STRICT", "").strip() else not (QUICK or SYNTHETIC)


def B(paper_value=PAPER_B):
    """Number of bootstrap replicates for an analysis that used `paper_value` replicates in the paper.
    Identity at the default NEC_BOOT = 2000; otherwise scaled proportionally (minimum 4)."""
    if BOOT == PAPER_B:
        return int(paper_value)
    return max(4, int(round(paper_value * BOOT / PAPER_B)))


def anchor_fail(msg):
    """Reproduction anchor not met: stop in strict mode, warn otherwise (synthetic data / smoke tests)."""
    if STRICT:
        sys.exit(msg)
    print(f"[non-strict mode] ANCHOR NOT MET (continuing): {msg}", flush=True)


def relpath(p):
    """Path as stored in result JSON files: relative to the repository root when inside it."""
    try:
        return str(Path(p).resolve().relative_to(REPO))
    except ValueError:
        return str(p)


for _d in (RESULTS_DIR, FIG_DIR):
    os.makedirs(_d, exist_ok=True)

if __name__ == "__main__":
    for k in ("REPO", "DATA_DIR", "WORK_DIR", "RESULTS_DIR", "FIG_DIR", "BOOT", "QUICK", "STRICT", "SYNTHETIC"):
        print(f"{k:12s} {globals()[k]}")
