#!/usr/bin/env python3
"""compare_super_bp_jaccard.py
============================
Compare two super-region sets at the BASE-PAIR level (robust to small boundary
shifts, unlike exact-coordinate Jaccard).

usage:
  python3 compare_super_bp_jaccard.py \
      --a rose_stitch_v3/super/H1_H3K27me3_t2500_signal.tsv \
      --b rose_out/H1_H3K27me3_t3000/H1_H3K27me3_peaks_SuperStitched.table.txt \
      --a-super-col isSuper
"""
import argparse, os
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--a", required=True)
    p.add_argument("--b", required=True)
    p.add_argument("--a-super-col", default="isSuper")
    p.add_argument("--b-filter", default="", help="optional: only rows where this col > 0")
    return p.parse_args()


def read_any(f):
    raw = open(f, encoding="utf-8", errors="ignore").read().splitlines()
    hi = 0
    for i, l in enumerate(raw[:20]):
        if "CHROM" in l.upper() or (chr(9) in l and ("start" in l.lower())):
            hi = i
            break
    sep = chr(9) if chr(9) in raw[hi] else r"\s+"
    t = pd.read_csv(f, sep=sep, skiprows=hi, engine="python")
    t.columns = [str(c).strip().strip(chr(34)) for c in t.columns]
    ren = {}
    for c in t.columns:
        cu = c.upper()
        if cu == "CHROM":
            ren[c] = "chr"
        elif cu == "START":
            ren[c] = "start"
        elif cu in ("STOP", "END"):
            ren[c] = "stop"
    return t.rename(columns=ren)


def merge_iv(df):
    """{chr: merged [(s,e), ...]}"""
    out = {}
    for c, sub in df.groupby("chr"):
        a = np.sort(sub[["start", "stop"]].values.astype(np.int64), axis=0)
        ms, me = [], []
        s, e = a[0]
        for s2, e2 in a[1:]:
            if s2 <= e:
                e = max(e, e2)
            else:
                ms.append(s); me.append(e); s, e = s2, e2
        ms.append(s); me.append(e)
        out[c] = (np.array(ms), np.array(me))
    return out


def bp(intervals):
    return {c: int((e - s).sum()) for c, (s, e) in intervals.items()}


def inter_bp(A, B):
    tot = {}
    for c in set(A) & set(B):
        a_s, a_e = A[c]; b_s, b_e = B[c]
        i = j = 0; L = 0
        while i < len(a_s) and j < len(b_s):
            lo = max(a_s[i], b_s[j]); hi = min(a_e[i], b_e[j])
            if hi > lo:
                L += int(hi - lo)
            if a_e[i] < b_e[j]:
                i += 1
            else:
                j += 1
        tot[c] = L
    return sum(tot.values())


def main():
    a = parse_args()
    A = read_any(a.a); B = read_any(a.b)
    if a.a_super_col in A.columns:
        A = A[A[a.a_super_col].astype(float) > 0]
    if a.b_filter and a.b_filter in B.columns:
        B = B[B[a.b_filter].astype(float) > 0]
    A = A[["chr", "start", "stop"]].astype({"start": int, "stop": int})
    B = B[["chr", "start", "stop"]].astype({"start": int, "stop": int})
    IA, IB = merge_iv(A), merge_iv(B)
    ba, bb = bp(IA), bp(IB)
    ta, tb = sum(ba.values()), sum(bb.values())
    ti = inter_bp(IA, IB)
    tu = ta + tb - ti
    print("A: %d regions  %.2f Mb" % (len(A), ta / 1e6))
    print("B: %d regions  %.2f Mb" % (len(B), tb / 1e6))
    print("intersection: %.2f Mb" % (ti / 1e6))
    print("bp Jaccard  : %.4f" % (ti / tu if tu else float("nan")))
    print("frac of A covered by B: %.4f" % (ti / ta if ta else float("nan")))
    print("frac of B covered by A: %.4f" % (ti / tb if tb else float("nan")))


if __name__ == "__main__":
    main()
