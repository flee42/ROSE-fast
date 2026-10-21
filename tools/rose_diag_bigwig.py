#!/usr/bin/env python3
"""rose_diag_bigwig.py -- diagnose why bw.stats() fails for the ROSE signal.

For one (cell, mark) it compares chromosome names between the peak BED and the
bigWig, then probes a few intervals with single-interval and list-mode calls.

usage:
  python3 rose_diag_bigwig.py --peak-manifest manifest.tsv \
                              --bigwig-manifest bigwig_manifest.tsv \
                              [--cell K562 --mark H3K27me3]
"""
import argparse
import sys
import numpy as np
import pandas as pd

TAB = chr(9)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peak-manifest", required=True,
                   help="TSV with cell, mark, bed")
    p.add_argument("--bigwig-manifest", required=True,
                   help="TSV with cell, mark, bigwig")
    p.add_argument("--cell", default="", help="inspect this cell (default: first match)")
    p.add_argument("--mark", default="", help="inspect this mark (default: first match)")
    p.add_argument("--n-peaks", type=int, default=5, help="peaks to use for the list-mode probe")
    return p.parse_args()


def main():
    a = parse_args()
    try:
        import pyBigWig
    except ImportError:
        sys.exit("pyBigWig is not installed:  pip install pyBigWig")

    peak_man = pd.read_csv(a.peak_manifest, sep=TAB, low_memory=False)
    bw_man = pd.read_csv(a.bigwig_manifest, sep=TAB, low_memory=False)
    print("peak manifest columns  :", list(peak_man.columns))
    print("bigwig manifest columns:", list(bw_man.columns))

    peak_man["_k"] = peak_man["cell"].astype(str) + "|" + peak_man["mark"].astype(str)
    bw_man["_k"] = bw_man["cell"].astype(str) + "|" + bw_man["mark"].astype(str)

    keys = [k for k in bw_man["_k"] if k in set(peak_man["_k"])]
    print("matched (cell|mark): %d  e.g. %s" % (len(keys), keys[:3]))
    if not keys:
        sys.exit("no matching (cell, mark) between the two manifests")

    if a.cell and a.mark:
        want = "%s|%s" % (a.cell, a.mark)
        if want not in keys:
            sys.exit("%s not found; available: %s" % (want, keys[:10]))
        k = want
    else:
        k = keys[0]
    cell, mark = k.split("|", 1)

    bed = str(peak_man[peak_man["_k"] == k].iloc[0]["bed"])
    bwp = str(bw_man[bw_man["_k"] == k].iloc[0]["bigwig"])
    print()
    print("=== %s / %s ===" % (cell, mark))
    print("peaks :", bed)
    print("bigwig:", bwp)

    peaks = pd.read_csv(bed, sep=TAB, header=None, usecols=[0, 1, 2],
                        names=["chr", "start", "end"])
    print()
    print("peak chroms (%d): %s" % (peaks["chr"].nunique(),
                                    list(pd.unique(peaks["chr"]))[:6]))
    print("peak start/end dtype:", peaks["start"].dtype)
    print(peaks.head(3).to_string(index=False))

    bw = pyBigWig.open(bwp)
    ch = bw.chroms()
    print()
    print("bw chroms (%d): %s" % (len(ch), list(ch)[:6]))
    print("bw header:", bw.header())

    r = peaks.iloc[0]
    c0, s0, e0 = str(r["chr"]), int(r["start"]), int(r["end"])
    print()
    print("--- probe single interval: %s:%d-%d ---" % (c0, s0, e0))
    cands = [c0, c0[3:] if c0.startswith("chr") else "chr" + c0]
    for cc in cands:
        if cc not in ch:
            print("  chrom %r NOT in bigWig" % cc)
            continue
        try:
            print("  bw.stats(%r,%d,%d) -> %r"
                  % (cc, s0, e0, bw.stats(cc, s0, e0, type="mean")))
        except Exception as ex:
            print("  bw.stats(%r) EXCEPTION: %r" % (cc, ex))

    print()
    print("--- probe LIST interval mode (first %d peaks) ---" % a.n_peaks)
    sub = peaks.head(a.n_peaks)
    try:
        print(bw.stats(c0, sub["start"].tolist(), sub["end"].tolist(), type="mean"))
    except Exception as ex:
        print("  LIST MODE EXCEPTION: %r" % (ex,))

    print()
    print("--- global non-empty check on the first chromosome ---")
    if ch:
        c1 = list(ch)[0]
        v = bw.stats(c1, 0, min(ch[c1], 10_000_000), nBins=10, type="mean")
        print("  %s 0-10Mb nBins=10 -> %r" % (c1, v))
    bw.close()


if __name__ == "__main__":
    main()
