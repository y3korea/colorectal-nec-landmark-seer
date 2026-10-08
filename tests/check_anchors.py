# -*- coding: utf-8 -*-
"""Reproduction anchors: assert the key numbers of the manuscript in the produced results (NEC_RESULTS_DIR).

    python3 tests/check_anchors.py            # exit code 1 if any anchor is not met
    python3 tests/check_anchors.py --soft     # report only (exit code 0), used for synthetic / QUICK runs

Anchors (rounded as printed in the manuscript):
  distant colorectal NEC, 12-month cancer-death probability (Aalen-Johansen) at landmarks 0/3/6/12/24 months
  = 72.29 / 71.62 / 72.33 / 60.03 / 30.51 %; main NEC n = 1,417, distant n = 908; numbers at risk 908/553/396/169/48;
  broken-stick plateau end 6.75 months; adjusted RR months 6-9 vs 3-6 = 1.49 (1.17-1.90);
  distant SCLC 69.98 / 69.49 / 61.93, colorectal adenocarcinoma 40.59 / 33.92 / 34.41 and pancreatic adenocarcinoma
  77.12 / 69.72 / 65.38 % at landmarks 0/6/12; A12 (lung LCNEC) rate ratios and hazard-shape tests.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import config  # noqa: E402

RES = config.RES
_cache = {}


def J(name):
    if name not in _cache:
        with open(os.path.join(RES, name + ".json"), encoding="utf-8") as f:
            _cache[name] = json.load(f)
    return _cache[name]


def r2(x):
    return None if x is None else round(float(x), 2)


A = []  # (group, label, getter, expected, tolerance)


def add(group, label, getter, expected, tol=0.0051):
    A.append((group, label, getter, expected, tol))


# ---- distant colorectal NEC: landmark probabilities, cohort sizes, plateau (A01) ----
for j, (s, v) in enumerate(zip((0, 3, 6, 12, 24), (72.29, 71.62, 72.33, 60.03, 30.51))):
    add("NEC landmark (A01)", f"12-mo cancer-death % at landmark {s}", lambda s=s: J("A01")["monthly_distant"]["cancer_pct"][s], v)
for s, v in zip((0, 3, 6, 12, 24), (908, 553, 396, 169, 48)):
    add("NEC landmark (A01)", f"n at risk at landmark {s}", lambda s=s: J("A01")["monthly_distant"]["n_at_risk"][s], v, 0)
add("NEC cohort (A01)", "main NEC n", lambda: J("A01")["n_main_nec"], 1417, 0)
add("NEC cohort (A01)", "distant NEC n", lambda: J("A01")["n_distant"], 908, 0)
add("NEC plateau (A01)", "broken-stick plateau end, months", lambda: J("A01")["monthly_distant"]["plateau_cancer"]["c"], 6.75, 1e-9)
add("NEC plateau (A01)", "plateau end 95% CI lower", lambda: J("A01")["monthly_distant"]["plateau_cancer"]["c_ci"][0], 4.00, 1e-9)
add("NEC plateau (A01)", "plateau end 95% CI upper", lambda: J("A01")["monthly_distant"]["plateau_cancer"]["c_ci"][1], 9.75, 1e-9)
# ---- adjusted rate ratio 6-9 vs 3-6 (A06; same model in A04 R part) ----
for q, v in (("RR", 1.49), ("lo", 1.17), ("hi", 1.90)):
    add("NEC RR 6-9 vs 3-6 (A06)", q, lambda q=q: r2(J("A06")["reproduction"]["adj_RR"]["r69_36"][q]), v)
    add("NEC RR 6-9 vs 3-6 (A04)", q, lambda q=q: r2(J("A04_rr")["main"]["model_3mo"]["rr_6_9_vs_3_6"][q]), v)
# ---- comparators (A03 point estimates; cross-checked in A11) ----
for t, vals in (("SCLC", (69.98, 69.49, 61.93)), ("AC", (40.59, 33.92, 34.41)), ("PANC", (77.12, 69.72, 65.38))):
    for s, v in zip((0, 6, 12), vals):
        add(f"{t} landmark (A03)", f"12-mo cancer-death % at landmark {s}", lambda t=t, s=s: J("A03")["point_cif12"][t][str(s)], v)
        add(f"{t} landmark (A11)", f"12-mo cancer-death % at landmark {s}", lambda t=t, s=s: J("A11")["table3"][t]["cancer"]["landmark"][str(s)]["pct"], v)
# ---- A12: lung LCNEC as comparator ----
add("A12 lung LCNEC", "n distant LCNEC", lambda: J("A12")["n"]["LCNEC"], 3185, 0)
add("A12 lung LCNEC", "n colorectal small-cell / large-cell", lambda: [J("A12")["n"]["CRC_small"], J("A12")["n"]["CRC_large"]], [599, 309], 0)
for q, v in (("RR", 0.97), ("lo", 0.84), ("hi", 1.11)):
    add("A12 lung LCNEC", f"LCNEC RR 6-9 vs 3-6 {q}", lambda q=q: r2(J("A12")["lcnec"]["rr_6_9_vs_3_6"][q]), v)
for q, v in (("RR", 0.79), ("lo", 0.72), ("hi", 0.86)):
    add("A12 lung LCNEC", f"LCNEC RR 6-15 vs 0-6 {q}", lambda q=q: r2(J("A12")["lcnec"]["rr_6_15_vs_0_6"][q]), v)
add("A12 shape tests", "NEC vs SCLC W (6 df)", lambda: r2(J("A12")["shape_tests"]["NEC_vs_SCLC_all_years"]["W"]), 12.05)
add("A12 shape tests", "NEC vs SCLC P", lambda: round(J("A12")["shape_tests"]["NEC_vs_SCLC_all_years"]["p"], 3), 0.061, 0.00051)
add("A12 shape tests", "NEC vs LCNEC W, all years", lambda: r2(J("A12")["shape_tests"]["NEC_vs_LCNEC_all_years"]["W"]), 15.35)
add("A12 shape tests", "NEC vs LCNEC P, all years", lambda: round(J("A12")["shape_tests"]["NEC_vs_LCNEC_all_years"]["p"], 3), 0.018, 0.00051)
add("A12 shape tests", "NEC vs LCNEC W, 2016-2023", lambda: r2(J("A12")["shape_tests"]["NEC_vs_LCNEC_2016_2023"]["W"]), 19.77)
add("A12 shape tests", "NEC/LCNEC ratio of RR(6-9 vs 3-6)", lambda: r2(J("A12")["shape_tests"]["NEC_vs_LCNEC_all_years"]["ratio_r69_36_NEC_over_cmp"]["RR"]), 1.54)
for cell, v in (("CRC_small", 1.75), ("CRC_large", 1.18), ("lung_small", 1.38), ("lung_large", 0.97)):
    add("A12 morphology x organ", f"RR 6-9 vs 3-6, {cell}", lambda cell=cell: r2(J("A12")["morphology_by_organ"]["rr_6_9_vs_3_6"][cell]["RR"]), v)
add("A12 morphology x organ", "small/large ratio, lung", lambda: r2(J("A12")["morphology_by_organ"]["ratios"]["small_over_large_lung"]["RR"]), 1.43)
add("A12 morphology x organ", "morphology x organ interaction", lambda: r2(J("A12")["morphology_by_organ"]["ratios"]["morphology_x_organ_interaction"]["RR"]), 1.04)
add("A12 internal", "all A12 reproduction anchors passed", lambda: J("A12")["reproduction"]["all_ok"], True, 0)


def close(a, b, tol):
    if isinstance(b, bool) or isinstance(a, bool):
        return a is b or a == b
    if isinstance(b, list):
        return isinstance(a, list) and len(a) == len(b) and all(close(x, y, tol) for x, y in zip(a, b))
    try:
        return abs(float(a) - float(b)) <= tol + 1e-12
    except (TypeError, ValueError):
        return a == b


def main():
    soft = "--soft" in sys.argv
    print(f"Checking reproduction anchors in {RES}")
    n_fail = 0
    for group, label, getter, exp, tol in A:
        try:
            got = getter(); ok = close(got, exp, tol)
        except Exception as e:  # missing file / key: anchor not met
            got, ok = f"<{type(e).__name__}: {e}>", False
        n_fail += not ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {group:28s} {label:42s} expected {exp!s:>12}  got {got!s}")
    print(f"\n{len(A) - n_fail}/{len(A)} anchors met.")
    if n_fail and not soft:
        sys.exit(1)


if __name__ == "__main__":
    main()
