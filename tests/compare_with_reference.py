# -*- coding: utf-8 -*-
"""Compare every produced results/*.json (and the manuscript figure_data.json) with reference_results/*.json.

    python3 tests/compare_with_reference.py [--tol 1e-6] [--soft] [--report path.json]

Rules
  * numbers (point estimates, CIs, P-values, counts): |produced - reference| <= tol (default 1e-6), or relative
    difference <= 1e-9 for large magnitudes; NaN/None must coincide
  * any key containing "seed": exact equality
  * strings: exact equality; file paths are compared by base name (absolute locations differ between machines)
  * dict keys must coincide; list lengths must coincide
  * not compared (machine / run specific): keys matching runtime*, date, software, R_version, muhaz_version,
    font, font_warnings, output_dir
  * suppressed cells: reference_results/ holds publication-safe copies made by tools/suppress_small_cells.py, in
    which cells based on 1-4 patients (SEER Research Data Use Agreement, item 14) are replaced by the marker
    "suppressed (SEER small-cell rule)"; such fields are skipped (counted and reported, not compared) and the metadata key
    "_suppression" is ignored. Every other value is compared as above; the produced output must still contain the
    suppressed field (same keys, same list lengths).
Every reference file must be produced; produced JSON files without a reference are listed but not compared.
Exit code 1 if any comparison fails (unless --soft).
"""
import argparse
import glob
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import config  # noqa: E402

REF_DIR = os.path.join(str(config.REPO), "reference_results")
IGNORE = re.compile(r"(^|_)runtime|^date$|^software$|^R_version$|^muhaz_version$|^font$|^font_warnings$|^output_dir$|^_suppression$")
SUPPRESSED = "suppressed (SEER small-cell rule)"   # marker written by tools/suppress_small_cells.py
PATHLIKE = re.compile(r"[\\/].*\.(png|pdf|json|csv|md|txt)$", re.I)


class Cmp:
    def __init__(self, tol):
        self.tol = tol; self.n = 0; self.fails = []; self.max_abs = 0.0; self.skipped = 0

    def fail(self, path, msg):
        self.fails.append(f"{path}: {msg}")

    def num(self, a, b, path, exact=False):
        self.n += 1
        fa, fb = float(a), float(b)
        if math.isnan(fa) or math.isnan(fb):
            if not (math.isnan(fa) and math.isnan(fb)):
                self.fail(path, f"{a!r} vs reference {b!r}")
            return
        if math.isinf(fa) or math.isinf(fb):
            if fa != fb:
                self.fail(path, f"{a!r} vs reference {b!r}")
            return
        d = abs(fa - fb)
        self.max_abs = max(self.max_abs, d) if not exact else self.max_abs
        if exact:
            if fa != fb:
                self.fail(path, f"seed {a!r} vs reference {b!r} (exact match required)")
        elif d > self.tol and d > 1e-9 * max(abs(fa), abs(fb)):
            self.fail(path, f"{a!r} vs reference {b!r} (|diff| = {d:.3g})")

    def walk(self, a, b, path="", seed=False):
        if isinstance(b, str) and b == SUPPRESSED:   # publication-safe reference: small cell, not compared
            self.skipped += 1
            return
        if isinstance(b, dict):
            if not isinstance(a, dict):
                return self.fail(path, f"type {type(a).__name__} vs reference dict")
            ka = {k for k in a if not IGNORE.search(k)}; kb = {k for k in b if not IGNORE.search(k)}
            for k in sorted(kb - ka): self.fail(f"{path}/{k}", "missing in produced output")
            for k in sorted(ka - kb): self.fail(f"{path}/{k}", "not in reference")
            for k in sorted(ka & kb):
                self.walk(a[k], b[k], f"{path}/{k}", seed or ("seed" in k.lower()))
            return
        if isinstance(b, list):
            if not isinstance(a, list):
                return self.fail(path, f"type {type(a).__name__} vs reference list")
            if len(a) != len(b):
                return self.fail(path, f"length {len(a)} vs reference {len(b)}")
            for i, (x, y) in enumerate(zip(a, b)):
                self.walk(x, y, f"{path}[{i}]", seed)
            return
        if isinstance(b, bool) or isinstance(a, bool):
            self.n += 1
            if a is not b and a != b:
                self.fail(path, f"{a!r} vs reference {b!r}")
            return
        if b is None or a is None:
            self.n += 1
            na = a is None or (isinstance(a, float) and math.isnan(a)); nb = b is None or (isinstance(b, float) and math.isnan(b))
            if not (na and nb):
                self.fail(path, f"{a!r} vs reference {b!r}")
            return
        if isinstance(b, (int, float)):
            if not isinstance(a, (int, float)):
                return self.fail(path, f"{a!r} vs reference {b!r}")
            return self.num(a, b, path, exact=seed)
        if isinstance(b, str):
            self.n += 1
            if a == b:
                return
            if isinstance(a, str) and (PATHLIKE.search(a) or PATHLIKE.search(b)) and os.path.basename(a.replace("\\", "/")) == os.path.basename(b.replace("\\", "/")):
                return
            return self.fail(path, f"{str(a)[:80]!r} vs reference {str(b)[:80]!r}")
        self.n += 1
        if a != b:
            self.fail(path, f"{a!r} vs reference {b!r}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--soft", action="store_true")
    ap.add_argument("--report", default=os.path.join(config.RESULTS_DIR, "logs", "compare_with_reference.json"))
    args = ap.parse_args()
    produced = {os.path.basename(p)[:-5]: p for p in glob.glob(os.path.join(config.RESULTS_DIR, "*.json"))}
    fd = os.path.join(config.MANUSCRIPT_FIG_DIR, "figure_data.json")
    if os.path.exists(fd):
        produced["figure_data"] = fd
    refs = {os.path.basename(p)[:-5]: p for p in sorted(glob.glob(os.path.join(REF_DIR, "*.json")))}
    print(f"Comparing {config.RESULTS_DIR} (+ figure_data.json) with {REF_DIR}; tolerance {args.tol:g}, seeds exact")
    report = {"tolerance": args.tol, "files": {}}
    total_fail = 0; total_skipped = 0
    for name, rp in refs.items():
        if name not in produced:
            print(f"  [FAIL] {name:16s} not produced"); total_fail += 1
            report["files"][name] = {"status": "missing"}; continue
        with open(produced[name], encoding="utf-8") as f: a = json.load(f)
        with open(rp, encoding="utf-8") as f: b = json.load(f)
        c = Cmp(args.tol); c.walk(a, b, name)
        total_fail += bool(c.fails)
        print(f"  [{'PASS' if not c.fails else 'FAIL'}] {name:16s} {c.n:6d} values compared, max |diff| {c.max_abs:.2e}, {len(c.fails)} mismatches"
              + (f"; {c.skipped} suppressed small-cell fields skipped" if c.skipped else ""))
        for m in c.fails[:15]: print("         ", m)
        if len(c.fails) > 15: print(f"          ... {len(c.fails) - 15} more")
        report["files"][name] = {"status": "pass" if not c.fails else "fail", "values_compared": c.n, "max_abs_diff": c.max_abs,
                                 "suppressed_skipped": c.skipped, "mismatches": c.fails[:200]}
        total_skipped += c.skipped
    extra = sorted(set(produced) - set(refs))
    if extra: print(f"  (produced without reference, not compared: {', '.join(extra)})")
    report["not_compared"] = extra; report["all_pass"] = total_fail == 0; report["suppressed_skipped"] = total_skipped
    os.makedirs(os.path.dirname(args.report), exist_ok=True)
    with open(args.report, "w") as f: json.dump(report, f, indent=1)
    if total_skipped:
        print(f"  ({total_skipped} small-cell fields are suppressed in reference_results/ and were not compared; see tools/suppress_small_cells.py)")
    print(f"\n{len(refs) - total_fail}/{len(refs)} reference files reproduced. Report: {args.report}")
    if total_fail and not args.soft:
        sys.exit(1)


if __name__ == "__main__":
    main()
