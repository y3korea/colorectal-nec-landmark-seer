#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""One-command pipeline: SEER*Stat exports -> cohorts -> analyses A01-A12 -> manuscript figures -> checks.

    python3 run_all.py                 # full run, paper settings (B = 2000), stops at the first failing step
    python3 run_all.py --quick         # smoke test: small bootstrap (NEC_BOOT=40), anchors informational, keep going
    python3 run_all.py --from A06      # resume from a step (cohorts in NEC_WORK_DIR are reused)
    python3 run_all.py --only A01,A03  # run selected steps
    python3 run_all.py --list          # show the steps

Paths and settings come from environment variables (see python/config.py and README):
NEC_DATA_DIR, NEC_WORK_DIR, NEC_RESULTS_DIR, NEC_FIG_DIR, NEC_BOOT, NEC_QUICK, NEC_STRICT, NEC_CORES.
Relative NEC_*_DIR values are resolved against the directory you start from and passed to every step as absolute paths.
Per-step logs: <NEC_RESULTS_DIR>/logs/.  Exit code 0 only if every step (and, in a full run, every check) passed.
"""
import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent

# (id, interpreter, script, description)
STEPS = [
    ("COHORTS", "R", "R/00_build_cohorts.R", "build cohorts from the four SEER exports (patient-level, NEC_WORK_DIR)"),
    ("FLOW", "py", "python/flow_counts.py", "aggregate cohort counts for the flow figure"),
    ("A01", "py", "python/A01_landmark_main.py", "landmark analysis, Table 2, plateau, sensitivity analyses"),
    ("A02", "py", "python/A02_window_sensitivity_model_comparison.py", "window-width sensitivity, plateau model comparison"),
    ("A03", "py", "python/A03_equivalence_tost.py", "equivalence (TOST) framing"),
    ("A04", "py", "python/A04_expanded_cohort_codes_subtypes.py", "case definitions, subtypes, era (+ R/A04 rate ratios)"),
    ("A05R", "R", "R/A05_heterogeneity_std.R", "standardized curves (multi-state Cox)"),
    ("A05", "py", "python/A05_heterogeneity_stratified.py", "population heterogeneity, strata, time-updated chemotherapy"),
    ("A06", "R", "R/A06_piecewise_poisson_suite.R", "piecewise-exponential (Poisson) suite"),
    ("A07", "py", "python/A07_split_sample_validation.py", "split-sample (calendar) validation (+ R/A07)"),
    ("A08", "R", "R/A08_standardization_extended.R", "standardization extended, PH check, Fine-Gray"),
    ("A09", "R", "R/A09_smoothed_hazards_muhaz.R", "kernel-smoothed hazards (muhaz)"),
    ("A09MD", "py", "python/A09_write_md.py", "A09 markdown summary"),
    ("A10", "py", "python/A10_data_quality_censoring.py", "data quality and censoring tables"),
    ("A11", "py", "python/A11_cross_tumour_comparison.py", "cross-tumour comparison (Table 3)"),
    ("A12", "R", "R/A12_lcnec_rate_ratios.R", "lung LCNEC rate ratios and shape tests"),
    ("FIG", "py", "python/FIG_manuscript_final.py", "manuscript Figures 1-3 and Supplementary Figures S1-S7"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="smoke test: NEC_BOOT=40 unless set, non-strict anchors, keep going")
    ap.add_argument("--boot", type=int, help="bootstrap replicates of the main analyses (sets NEC_BOOT)")
    ap.add_argument("--strict", action="store_true", help="stop when a reproduction anchor is not met (NEC_STRICT=1)")
    ap.add_argument("--keep-going", action="store_true", help="continue after a failing step")
    ap.add_argument("--from", dest="from_", metavar="STEP", help="start at this step")
    ap.add_argument("--only", metavar="STEPS", help="comma-separated steps to run")
    ap.add_argument("--no-tests", action="store_true", help="skip tests/check_anchors.py, tests/compare_with_reference.py and tests/check_small_cell_text.py")
    ap.add_argument("--list", action="store_true", help="list the steps and exit")
    a = ap.parse_args()
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    if a.list:
        for sid, kind, script, desc in STEPS:
            print(f"{sid:8s} {kind:3s} {script:52s} {desc}")
        return 0
    if a.quick:
        os.environ["NEC_QUICK"] = "1"
    if a.boot:
        os.environ["NEC_BOOT"] = str(a.boot)
    if a.strict:
        os.environ["NEC_STRICT"] = "1"
    sys.path.insert(0, str(REPO / "python"))
    import config  # noqa: E402  (reads the environment set above)

    keep_going = a.keep_going or config.QUICK
    ids = [s[0] for s in STEPS]
    sel = list(STEPS)
    if a.from_:
        if a.from_ not in ids: ap.error(f"unknown step {a.from_}; choose from {ids}")
        sel = STEPS[ids.index(a.from_):]
    if a.only:
        want = [x.strip() for x in a.only.split(",") if x.strip()]
        bad = [x for x in want if x not in ids]
        if bad: ap.error(f"unknown step(s) {bad}; choose from {ids}")
        sel = [s for s in STEPS if s[0] in want]

    print("colorectal-nec-landmark-seer pipeline")
    for k in ("DATA_DIR", "WORK_DIR", "RESULTS_DIR", "FIG_DIR", "BOOT", "QUICK", "STRICT", "SYNTHETIC"):
        print(f"  {k:12s} {getattr(config, k)}")
    if config.SYNTHETIC:
        print("  NOTE: synthetic data detected (SYNTHETIC_DATA_README.txt): results are NOT those of the paper.")
    if shutil.which("Rscript") is None:
        print("ERROR: Rscript not found on PATH (R >= 4.3 required; see R/install_packages.R)"); return 2
    if any(s[0] == "COHORTS" for s in sel):
        missing = [f for f in config.RAW_FILES.values() if not os.path.exists(os.path.join(config.DATA_DIR, f))]
        if missing:
            print("ERROR: SEER exports not found in NEC_DATA_DIR:\n  " + "\n  ".join(missing) +
                  "\nSee data_spec/SEER_export_specification.md, or set NEC_FILE_NEN/NEC_FILE_AC/NEC_FILE_LUNG/NEC_FILE_PANC."); return 2

    logdir = Path(config.RESULTS_DIR) / "logs"; logdir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ); env.setdefault("OMP_NUM_THREADS", "1"); env.setdefault("OPENBLAS_NUM_THREADS", "1"); env.setdefault("MKL_NUM_THREADS", "1")
    env["PYTHONUNBUFFERED"] = "1"
    # Pass the resolved absolute locations to every step, so that relative NEC_* values are interpreted once,
    # against the directory from which run_all.py / run_all.sh was started (steps run with cwd = repository root).
    env.update({"NEC_DATA_DIR": config.DATA_DIR, "NEC_WORK_DIR": config.WORK_DIR,
                "NEC_RESULTS_DIR": config.RESULTS_DIR, "NEC_FIG_DIR": config.FIG_DIR})
    manifest = {"started": time.strftime("%Y-%m-%d %H:%M:%S"), "python": platform.python_version(), "platform": platform.platform(),
                "boot": config.BOOT, "quick": config.QUICK, "strict": config.STRICT, "synthetic": config.SYNTHETIC, "steps": []}
    T0 = time.time(); failed = []
    for i, (sid, kind, script, desc) in enumerate(sel, 1):
        cmd = (["Rscript", str(REPO / script)] if kind == "R" else [sys.executable, str(REPO / script)])
        log = logdir / f"{ids.index(sid) + 1:02d}_{sid}.log"
        print(f"[{i:2d}/{len(sel)}] {sid:8s} {desc} ... ", end="", flush=True)
        t = time.time()
        with open(log, "w") as lf:
            rc = subprocess.run(cmd, cwd=str(REPO), env=env, stdout=lf, stderr=subprocess.STDOUT).returncode
        dt = time.time() - t
        ok = rc == 0
        print(f"{'ok' if ok else 'FAILED (exit ' + str(rc) + ')'}  {dt / 60:.1f} min  [log: {log.name}]")
        manifest["steps"].append({"step": sid, "script": script, "exit_code": rc, "minutes": round(dt / 60, 2)})
        if not ok:
            failed.append(sid)
            tail = log.read_text(errors="replace").splitlines()[-12:]
            print("        " + "\n        ".join(tail))
            if not keep_going:
                break
    manifest["pipeline_minutes"] = round((time.time() - T0) / 60, 2)
    print(f"\nPipeline: {len(sel) - len(failed)}/{len(sel)} steps ok in {manifest['pipeline_minutes']} min" + (f"; failed: {failed}" if failed else ""))

    tests_ok = True
    if not a.no_tests and (not failed or keep_going):
        soft = not config.STRICT
        print("\n== tests/check_anchors.py" + (" (informational: non-strict run)" if soft else ""))
        r1 = subprocess.run([sys.executable, str(REPO / "tests/check_anchors.py")] + (["--soft"] if soft else []), cwd=str(REPO), env=env).returncode
        manifest["check_anchors_exit"] = r1
        compare = (config.BOOT == config.PAPER_B) and not config.SYNTHETIC
        if compare:
            print("\n== tests/compare_with_reference.py")
            r2 = subprocess.run([sys.executable, str(REPO / "tests/compare_with_reference.py")], cwd=str(REPO), env=env).returncode
            manifest["compare_with_reference_exit"] = r2
        else:
            r2 = 0
            print("\n== tests/compare_with_reference.py skipped (bootstrap size differs from the paper or synthetic data)")
        print("\n== tests/check_small_cell_text.py (code, docs, notebooks and the results/*.md reports)")
        r3 = subprocess.run([sys.executable, str(REPO / "tests/check_small_cell_text.py"), "--also", str(config.RESULTS_DIR)],
                            cwd=str(REPO), env=env).returncode
        manifest["check_small_cell_text_exit"] = r3
        tests_ok = (r1 == 0 and r2 == 0 and r3 == 0) or soft
    manifest["total_minutes"] = round((time.time() - T0) / 60, 2)
    with open(logdir / "run_manifest.json", "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"\nTotal wall-clock {manifest['total_minutes']} min. Manifest: {logdir / 'run_manifest.json'}")
    print("Outputs in results/ are unsuppressed working files. Before publishing any of them run "
          "python3 tools/suppress_small_cells.py (SEER DUA small-cell rule); Graphical Abstract: python3 python/figures/graphical_abstract.py")
    return 0 if (not failed and tests_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
