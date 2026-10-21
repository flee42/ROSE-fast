#!/usr/bin/env python3
"""verify_fast_gff.py -- diagnose the difference between ROSE_bamToGFF.py and
fast_bam_to_gff.py outputs.

usage: python3 verify_fast_gff.py /tmp/old.gff /tmp/new.gff
"""
import sys
import numpy as np
import pandas as pd


def load(f):
    d = pd.read_csv(f, sep="\t")
    d.columns = ["id", "locus", "val"]
    d["val"] = pd.to_numeric(d["val"], errors="coerce")
    d["chr"] = d["locus"].str.split("(").str[0]
    se = d["locus"].str.extract(r":(\d+)-(\d+)$").astype(np.int64)
    d["start"], d["end"] = se[0], se[1]
    d["span"] = d["end"] - d["start"]
    return d


a = load(sys.argv[1]); b = load(sys.argv[2])
print("rows: old=%d new=%d" % (len(a), len(b)))
m = a.merge(b, on=["id", "locus"], suffixes=("_old", "_new"))
print("matched: %d" % len(m))
d = m["val_new"] - m["val_old"]
print("\nold mean=%.4f  new mean=%.4f  ratio=%.4f" %
      (m.val_old.mean(), m.val_new.mean(), m.val_new.mean() / max(m.val_old.mean(), 1e-12)))
print("diff: max|.|=%.4f  median|.|=%.4f" % (d.abs().max(), d.abs().median()))
nz = m[(m.val_old > 0) & (m.val_new > 0)]
r = nz.val_new / nz.val_old
print("\nratio new/old (both>0, n=%d):" % len(nz))
print("  p10=%.4f  median=%.4f  p90=%.4f" % (r.quantile(.1), r.median(), r.quantile(.9)))
print("  pearson=%.4f  spearman=%.4f" %
      (np.corrcoef(nz.val_old, nz.val_new)[0, 1],
       np.corrcoef(nz.val_old.rank(), nz.val_new.rank())[0, 1]))
print("\nzero counts: old=%d new=%d both-zero=%d" %
      ((m.val_old == 0).sum(), (m.val_new == 0).sum(),
       ((m.val_old == 0) & (m.val_new == 0)).sum()))
print("\nfirst 6 rows:")
print(m[["id", "span_old", "val_old", "val_new"]].head(6).to_string(index=False))
print("\nlargest 6 differences:")
w = m.assign(ad=d.abs()).nlargest(6, "ad")
print(w[["id", "span_old", "val_old", "val_new", "ad"]].to_string(index=False))
print("\ndiff vs region length (spearman): %.4f" %
      np.corrcoef(m["span_old"].rank(), d.abs().rank())[0, 1])
