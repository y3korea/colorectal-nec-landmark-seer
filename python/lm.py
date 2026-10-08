# -*- coding: utf-8 -*-
"""Shared landmark / Aalen-Johansen toolkit (independent re-implementation of the authors' AJ-curve / plateau / bootstrap script).
Time unit: half-days (k = round(2*days)); months = days/(365.25/12). ev: 0 censored, 1 cancer death, 2 other death.
"""
import numpy as np, pandas as pd
import config
MD = 365.25 / 12
DATA = config.DATA          # derived cohort files (NEC_WORK_DIR)

def load(name):
    """name in {'NEN_all','AC','LUNG','PANC','PANC_cyt','LUNG_cyt','AC_cyt'}"""
    return pd.read_csv(DATA + f"cohort_{name}.csv", dtype={"pid": str, "hist": str, "s1": str, "s2": str, "seq": str, "grade17": str, "grade18": str,
                                                           "tt_days": str, "smonths": str, "ntum": str}, low_memory=False)

NEC_CODES = ["8013", "8041", "8042", "8043", "8044", "8045"]
SCLC_CODES = ["8041", "8042", "8043", "8044", "8045"]
ISLET = [str(c) for c in range(8150, 8158)]

def nec_main(nen):      return nen[nen["hist"].isin(NEC_CODES)]
def nec_expanded(nen):
    """pre-specified expanded definition: main codes + 8246/3 with registry grade 3 or 4"""
    g3 = nen.grade17.isin(["Poorly differentiated; Grade III", "Undifferentiated; anaplastic; Grade IV"]) | nen.grade18.isin(["Site-specific grade system category (3)", "Site-specific grade system category (4)"])
    return nen[nen["hist"].isin(NEC_CODES) | ((nen["hist"] == "8246") & g3)]
def sclc_main(lung):    return lung[lung["hist"].isin(SCLC_CODES)]
def panc_main(panc):    return panc[~panc["hist"].isin(ISLET)]

def prep(df, kmax_months=36 + 12):
    KMAX = int(np.floor(2 * MD * kmax_months)) + 1
    d = np.maximum(df["sdays"].to_numpy(float), 0.5)
    k = np.rint(2 * d).astype(np.int64)
    return np.minimum(k, KMAX + 1), df["ev"].to_numpy(int), KMAX

def curves(k, ev, KMAX, landmarks=np.arange(0, 25), window=12, w=None):
    """Return (CIF cancer death %, all-cause death %) within `window` months after each landmark (patients alive & under follow-up at the landmark)."""
    if w is None: w = np.ones(len(k))
    cnt = np.bincount(k, weights=w, minlength=KMAX + 2)
    d1 = np.bincount(k[ev == 1], weights=w[ev == 1], minlength=KMAX + 2)
    d2 = np.bincount(k[ev == 2], weights=w[ev == 2], minlength=KMAX + 2)
    n = np.cumsum(cnt[::-1])[::-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        h1 = np.where(n > 0, d1 / n, 0.0); ha = np.where(n > 0, (d1 + d2) / n, 0.0)
    cif = np.full(len(landmarks), np.nan); allc = np.full(len(landmarks), np.nan); nrisk = np.zeros(len(landmarks))
    for j, s in enumerate(landmarks):
        lo = int(np.floor(2 * MD * s)) + 1; hi = min(int(np.floor(2 * MD * (s + window))), KMAX + 1)
        nrisk[j] = n[lo] if lo < len(n) else 0
        if nrisk[j] == 0: continue
        surv_prev = np.concatenate(([1.0], np.cumprod(1 - ha[lo:hi + 1])[:-1]))
        cif[j] = np.sum(surv_prev * h1[lo:hi + 1]) * 100
        allc[j] = (1 - np.prod(1 - ha[lo:hi + 1])) * 100
    return cif, allc, nrisk

GRID = np.arange(0, 15.0001, 0.25)
def hinge_many(Y, m=np.arange(0, 19).astype(float)):
    """broken-stick y = b0 + b1*max(0, m-c); grid search c in [0,15] step 0.25; first argmin (ties -> smallest c)"""
    sse = np.empty((Y.shape[0], len(GRID)))
    for g, c0 in enumerate(GRID):
        X = np.column_stack([np.ones_like(m), np.maximum(0, m - c0)])
        P = X @ np.linalg.pinv(X)
        res = Y - Y @ P.T
        sse[:, g] = np.sum(res ** 2, axis=1)
    return GRID[np.argmin(sse, axis=1)]

def bootstrap(df, B=2000, seed=20260928, landmarks=np.arange(0, 25), window=12):
    k, ev, KMAX = prep(df); n = len(k)
    rng = np.random.default_rng(seed)
    bc = np.empty((B, len(landmarks))); ba = np.empty((B, len(landmarks)))
    for b in range(B):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        bc[b], ba[b], _ = curves(k, ev, KMAX, landmarks, window, w)
    return bc, ba

def crude_rates(df, cuts=(0, 3, 6, 9, 12, 15, 18, 24)):
    """cancer deaths per person-month in intervals (cause-specific, other deaths censor at death)"""
    tm = df.tm.to_numpy(float); ev = df.ev.to_numpy(int); out = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        at = tm > a
        expo = np.clip(np.minimum(tm, b) - a, 0, None)[at].sum()
        d = ((ev == 1) & (tm > a) & (tm <= b)).sum()
        out.append({"start": a, "end": b, "deaths": int(d), "person_months": float(expo), "rate": float(d / expo) if expo > 0 else np.nan})
    return pd.DataFrame(out)

def pct(v, q=(2.5, 97.5)): return [float(x) for x in np.nanpercentile(v, q)]

def surgery(df):
    """primary-site surgery (authors' rule): returns Series: 1 yes, 0 no, -1 unknown/NA"""
    s1 = df.s1.fillna("Blank(s)"); s2 = df.s2.fillna("Blank(s)")
    no = (s1 == "00") | (s2 == "A000")
    unk = s1.isin(["99", "98"]) | s2.isin(["A990", "A980"]) | (s1.isin(["Blank(s)", ""]) & s2.isin(["Blank(s)", ""]))
    out = pd.Series(1, index=df.index); out[no] = 0; out[unk & ~no] = -1
    return out
