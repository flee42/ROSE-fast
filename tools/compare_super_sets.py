#!/usr/bin/env python3
"""compare_super_sets.py
======================
Compare the final SUPER region set of a ROSE run against a super set derived
from the same REGION_MAP with our own cutoff implementation.

Two modes
---------
  --mode cutoff   : take ROSE's *_AllStitched.table.txt, recompute the cutoff from
                    its own rank/control columns, and compare with ROSE's
                    *_SuperStitched.table.txt   -> tests the CUTOFF step only
  --mode sets     : compare two explicit region tables (ours vs ROSE) -> Jaccard

usage:
  python3 compare_super_sets.py --mode cutoff --rose-outdir rose_out/K562_H3K27me3_t3000
  python3 compare_super_sets.py --mode sets --a ours_SuperStitched.table.txt --b rose_SuperStitched.table.txt
"""
import argparse, glob, os, sys
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["cutoff", "sets"], default="cutoff")
    p.add_argument("--rose-outdir", default="")
    p.add_argument("--a", default="")
    p.add_argument("--b", default="")
    p.add_argument("--outdir", default="")
    return p.parse_args()


def read_rose_table(f):
    raw = open(f, encoding="utf-8", errors="ignore").read().splitlines()
    hi = 0
    for i, l in enumerate(raw[:20]):
        if "CHROM" in l.upper():
            hi = i
            break
    sep = chr(9) if chr(9) in raw[hi] else r"\s+"
    t = pd.read_csv(f, sep=sep, skiprows=hi, engine="python")
    t.columns = [str(c).strip().strip(chr(34)) for c in t.columns]
    ren = {}
    for c in t.columns:
        cu = c.upper()
        if cu in ("REGION_ID", "ID"):
            ren[c] = "id"
        elif cu == "CHROM":
            ren[c] = "chr"
        elif cu == "START":
            ren[c] = "start"
        elif cu in ("STOP", "END"):
            ren[c] = "stop"
    t = t.rename(columns=ren)
    t["start"] = t["start"].astype(int); t["stop"] = t["stop"].astype(int)
    return t


def key_set(t):
    return set(zip(t["chr"], t["start"], t["stop"]))


def jacc(a, b):
    u = len(a | b)
    return (len(a & b) / u) if u else float("nan")


def our_cutoff(v):
    """exhaustive version of ROSE_callSuper.R::calculate_cutoff"""
    v = np.sort(np.asarray(v, dtype=float).copy())
    v[v < 0] = 0.0
    n = len(v)
    slope = (v.max() - v.min()) / n
    xs = np.arange(1, n + 1, dtype=float)
    best = None
    for x in range(1, n + 1):
        b = v[x - 1] - slope * x
        obj = int(np.sum(v <= xs * slope + b))
        if best is None or obj < best[1]:
            best = (x, obj)
    return float(v[best[0] - 1]), best[0]


def main():
    a = parse_args()
    if a.mode == "sets":
        A = key_set(read_rose_table(a.a)); B = key_set(read_rose_table(a.b))
        print("ours  : %d" % len(A)); print("ROSE  : %d" % len(B))
        print("inter : %d" % len(A & B))
        print("Jaccard: %.4f" % jacc(A, B))
        print("only ours: %d   only ROSE: %d" % (len(A - B), len(B - A)))
        return

    d = a.rose_outdir
    allf = glob.glob(os.path.join(d, "*_AllStitched.table.txt"))
    supf = glob.glob(os.path.join(d, "*_SuperStitched.table.txt"))
    if not allf or not supf:
        raise SystemExit("need *_AllStitched.table.txt and *_SuperStitched.table.txt in %s" % d)
    T = read_rose_table(allf[0]); S = read_rose_table(supf[0])
    sigcols = [c for c in T.columns if c not in ("id", "chr", "start", "stop")
               and T[c].dtype.kind in "if"]
    # ---- pick the right columns by NAME ----
    ctl_cols = [c for c in T.columns if "control" in c.lower()]
    rank_cols = [c for c in T.columns if c.lower().endswith(".bam") and "control" not in c.lower()]
    print("rank cols   :", rank_cols)
    print("control cols:", ctl_cols)
    print("all cols    :", list(T.columns))
    if not rank_cols or not ctl_cols:
        raise SystemExit("cannot identify rank / control signal columns")
    sig = T[rank_cols[0]].astype(float) - T[ctl_cols[0]].astype(float)
    sig = sig.clip(lower=0)
    print("signal = %s - %s" % (rank_cols[0], ctl_cols[0]))

    cut, xstar = our_cutoff(sig.values)
    ours_flag = (sig >= cut).values
    print("cutoff = %.4f (x*=%d)  ours_frac = %.4f" % (cut, xstar, ours_flag.mean()))

    if "isSuper" in T.columns:
        rose_flag = T["isSuper"].astype(int).values.astype(bool)
        inter = int((ours_flag & rose_flag).sum())
        u = int((ours_flag | rose_flag).sum())
        print()
        print("ROW-WISE comparison (same region set, n=%d)" % len(T))
        print("  ROSE isSuper : %d (%.4f)" % (rose_flag.sum(), rose_flag.mean()))
        print("  ours         : %d (%.4f)" % (ours_flag.sum(), ours_flag.mean()))
        print("  intersection : %d" % inter)
        print("  Jaccard      : %.4f" % (inter / u if u else float("nan")))
        print("  precision    : %.4f" % (inter / max(ours_flag.sum(), 1)))
        print("  recall       : %.4f" % (inter / max(rose_flag.sum(), 1)))

    # ---- set comparison against the SuperStitched table ----
    ours = T.loc[ours_flag, ["chr", "start", "stop"]]
    A, B = key_set(ours), key_set(S)
    print()
    print("SET comparison")
    print("  ours super : %d" % len(A))
    print("  ROSE super : %d" % len(B))
    print("  inter      : %d" % len(A & B))
    print("  Jaccard    : %.4f" % jacc(A, B))
    print("  only ours  : %d   only ROSE: %d" % (len(A - B), len(B - A)))
    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)
        ours.to_csv(os.path.join(a.outdir, "ours_super_from_ROSE_signal.tsv"),
                    sep="\t", index=False)


if __name__ == "__main__":
    main()
