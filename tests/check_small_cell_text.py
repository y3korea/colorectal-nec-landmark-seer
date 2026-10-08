# -*- coding: utf-8 -*-
"""Scan text files for demographic or geographic statistics based on 1-4 patients (SEER DUA item 14).

tools/suppress_small_cells.py audits the result JSON files. This check covers free text: the code (which writes
sentences into the Markdown reports), the documentation, the notebooks and, optionally, the Markdown reports a run
wrote to results/. It flags a race/ethnicity, sex, age or geography label followed within the same clause by a count of
1-4 (forms "n = k", "(n=k)" and "<label> k (x%)" with k = 1-4) and JSON-style pairs "<label>": k. The scan is
heuristic; review every hit, reword without the count (or pool the category) and re-run.

    python3 tests/check_small_cell_text.py                    # files tracked by git (or all repository files)
    python3 tests/check_small_cell_text.py --also results     # plus the Markdown reports of a run
    python3 tests/check_small_cell_text.py --root DIR         # another tree (e.g. an extracted older commit)
Exit code 0 = no hits, 1 = hits.
"""
import argparse
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = (".py", ".R", ".r", ".md", ".ipynb", ".cff", ".txt", ".json", ".sh")
SKIP_DIRS = {".git", "work", "data", "results", "figures", "__pycache__", ".ipynb_checkpoints"}
LABEL = (r"(?:race|racial|ethnic\w*|hispanic|white|black|asian|pacific islander|american indian|alaska native|aian|"
         r"unknown race|males?|females?|men|women|sex|aged?|years old|registry|registries|county|counties|state)")
COUNT = r"(?:\(\s*)?n\s*=\s*[1-4](?![\d,.])"
PATTERNS = [
    # label ... n = 1-4 within one clause
    re.compile(r"\b" + LABEL + r"\b[^.;:\n]{0,60}?" + COUNT, re.I),
    # label followed directly by a bare count of 1-4 and a percentage in parentheses
    re.compile(r"\b" + LABEL + r"\b[^.;:\n\d]{0,40}?\s[1-4]\s*\(\s*\d+(?:\.\d+)?\s*%?\s*\)", re.I),
    # JSON / dict pairs: a quoted key containing a label, then ": k" with k = 1-4
    re.compile(r"[\"'][^\"'\n]{0,80}\b" + LABEL + r"\b[^\"'\n]{0,80}[\"']\s*:\s*[1-4](?![\d.])", re.I),
]
# Lines that state the rule itself are not statistics.
ALLOW = re.compile(r"1-4 patients|1–4 patients|counts 1 to 4|1 to 4|suppressed \(<5\)|SMALL = range|is_small|"
                   r"DEMOG_LABEL|PATTERNS|LABEL =", re.I)


def tracked_files(root):
    try:
        top = subprocess.run(["git", "-C", root, "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                             check=True).stdout.strip()
        if os.path.realpath(top) == os.path.realpath(root):
            out = subprocess.run(["git", "-C", root, "ls-files"], capture_output=True, text=True, check=True).stdout
            files = [os.path.join(root, f) for f in out.splitlines() if f]
            if files:
                return files
    except Exception:
        pass
    files = []
    for d, dirs, fs in os.walk(root):
        dirs[:] = [x for x in dirs if x not in SKIP_DIRS]
        files += [os.path.join(d, f) for f in fs]
    return files


def scan(path):
    hits = []
    try:
        with open(path, encoding="utf-8") as f:
            for i, line in enumerate(f, 1):
                if ALLOW.search(line):
                    continue
                for p in PATTERNS:
                    m = p.search(line)
                    if m:
                        hits.append((i, m.group(0)[:120]))
                        break
    except (UnicodeDecodeError, OSError):
        pass
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--also", nargs="*", default=[], help="extra folders to scan (e.g. results)")
    ap.add_argument("--root", default=ROOT, help="repository tree to scan (default: this repository)")
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    files = [f for f in tracked_files(root) if f.endswith(EXT) and os.path.isfile(f)
             and os.path.relpath(f, root) != os.path.join("tests", "check_small_cell_text.py")]
    for d in a.also:
        for dd, _, fs in os.walk(d):
            files += [os.path.join(dd, f) for f in fs if f.endswith((".md", ".txt"))]
    total = 0
    for f in sorted(set(files)):
        for ln, txt in scan(f):
            total += 1
            print(f"{os.path.relpath(f, root)}:{ln}: {txt}")
    if total:
        print(f"\n{total} possible demographic/geographic statistic(s) based on 1-4 patients; review and reword.")
        sys.exit(1)
    print(f"small-cell text check: {len(set(files))} files scanned, no demographic/geographic counts of 1-4 found")


if __name__ == "__main__":
    main()
