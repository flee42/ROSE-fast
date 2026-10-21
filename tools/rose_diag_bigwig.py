#!/usr/bin/env python3
"""rose_diag_bigwig.py -- why does bw.stats() fail for the ROSE signal?

Tests one (cell, mark): compares chromosome names between the peak file and the
bigWig, and probes a few intervals with both single and list calls.
"""
import sys, pandas as pd, numpy as np, pyBigWig

MAN = sys.argv[1] if len(sys.argv) > 1 else "domwidth_tier12.tsv"
BWM = sys.argv[2] if len(sys.argv) > 2 else "cv_sensitivity_v3/manifest.tsv"

peak_man = pd.read_csv(MAN, sep="\t")
bw_man = pd.read_csv(BWM, sep="\t")
print("peak manifest columns:", list(peak_man.columns))
print("bigwig manifest columns:", list(bw_man.columns))

bw_man["k"] = bw_man["cell"].astype(str) + "|" + bw_man["mark"].astype(str)
peak_man["k"] = peak_man["cell"].astype(str) + "|" + peak_man["mark"].astype(str)
common = [k for k in bw_man["k"] if k in set(peak_man["k"])]
print("matched (cell|mark):", len(common), "e.g.", common[:3])
if not common:
    sys.exit("no matching (cell,mark) between the two manifests")

k = common[0]
cell, mark = k.split("|")
bed = str(peak_man[peak_man["k"] == k].iloc[0]["bed"])
bwp = str(bw_man[bw_man["k"] == k].iloc[0]["bigwig"])
print("\n=== %s / %s ===" % (cell, mark))
print("peaks :", bed)
print("bigwig:", bwp)

peaks = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chr", "start", "end"])
print("\npeak chroms (%d): %s" % (peaks["chr"].nunique(), list(pd.unique(peaks["chr"]))[:6]))
print("peak start/end dtype:", peaks["start"].dtype)
print(peaks.head(3).to_string(index=False))

bw = pyBigWig.open(bwp)
ch = bw.chroms()
print("\nbw chroms (%d): %s" % (len(ch), list(ch)[:6]))
print("bw header:", bw.header())

r = peaks.iloc[0]
c0, s0, e0 = str(r["chr"]), int(r["start"]), int(r["end"])
print("\n--- probe single interval: %s:%d-%d ---" % (c0, s0, e0))
for cc in [c0, c0[3:] if c0.startswith("chr") else "chr" + c0]:
    if cc not in ch:
        print("  chrom %r NOT in bigWig" % cc); continue
    try:
        print("  bw.stats(%r,%d,%d) -> %r" % (cc, s0, e0, bw.stats(cc, s0, e0, type="mean")))
    except Exception as ex:
        print("  bw.stats(%r) EXCEPTION: %r" % (cc, ex))

print("\n--- probe LIST interval mode (first 5 peaks) ---")
sub = peaks.head(5)
try:
    print(bw.stats(c0, sub["start"].tolist(), sub["end"].tolist(), type="mean"))
except Exception as ex:
    print("  LIST MODE EXCEPTION: %r" % (ex,))

print("\n--- global non-empty check on the first chromosome ---")
if ch:
    c1 = list(ch)[0]
    v = bw.stats(c1, 0, min(ch[c1], 10_000_000), nBins=10, type="mean")
    print("  %s 0-10Mb nBins=10 -> %r" % (c1, v))
bw.close()
