#!/usr/bin/env python3
"""verify_fast_gff.py -- compare two ROSE per-region GFFs (old vs fast mapper).

The GFF written by ROSE_bamToGFF.py / fast_bam_to_gff.py is 3 columns:
    <region_id>  <locus "chr:start-end,start-end,...">  <value>

Reports row counts, how many regions match exactly, the value ratio and
correlation, and where the largest differences sit.

usage:
  python3 verify_fast_gff.py old.gff new.gff
"""
import argparse
import sys
import numpy as np
import pandas as pd

TAB = chr(9)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("old", help="GFF produced by ROSE_bamToGFF.py")
    p.add_argument("new", help="GFF produced by fast_bam_to_gff.py")
    p.add_argument("--top", type=int, default=6, help="how many worst rows to print")
    return p.parse_args()


def load(f):
    d = pd.read_csv(f, sep=TAB, header=None, usecols=[0, 1, 2],
                    names=["id", "locus", "val"])
    d["val"] = pd.to_numeric(d["val"], errors="coerce")
    d["chr"] = d["locus"].str.split("(").str[0]
    se = d["locus"].str.extract(r":(\d+)-(\d+)$").astype(np.int64)
    d["start"], d["end"] = se[0], se[1]
    d["span"] = d["end"] - d["start"]
    return d


def main():
    a = parse_args()
    old, new = load(a.old), load(a.new)
    print("rows: old=%d new=%d" % (len(old), len(new)))

    m = old.merge(new, on=["id", "locus"], suffixes=("_old", "_new"))
    print("matched: %d" % len(m))
    if m.empty:
        sys.exit("no common (id, locus) rows -- are these the same input GFF?")

    d = m["val_new"] - m["val_old"]
    print()
    print("old mean=%.4f  new mean=%.4f  ratio=%.6f"
          % (m.val_old.mean(), m.val_new.mean(),
             m.val_new.mean() / max(m.val_old.mean(), 1e-12)))
    print("diff: max|.|=%.4f  median|.|=%.4f" % (d.abs().max(), d.abs().median()))

    nz = m[(m.val_old > 0) & (m.val_new > 0)]
    if len(nz):
        r = nz.val_new / nz.val_old
        print()
        print("ratio new/old (both>0, n=%d):" % len(nz))
        print("  p10=%.4f  median=%.4f  p90=%.4f"
              % (r.quantile(.1), r.median(), r.quantile(.9)))
        print("  pearson=%.6f  spearman=%.6f"
              % (np.corrcoef(nz.val_old, nz.val_new)[0, 1],
                 np.corrcoef(nz.val_old.rank(), nz.val_new.rank())[0, 1]))

    print()
    print("zero counts: old=%d new=%d both-zero=%d"
          % ((m.val_old == 0).sum(), (m.val_new == 0).sum(),
             ((m.val_old == 0) & (m.val_new == 0)).sum()))

    print()
    print("first %d rows:" % a.top)
    print(m[["id", "span_old", "val_old", "val_new"]].head(a.top).to_string(index=False))

    print()
    print("largest %d differences:" % a.top)
    w = m.assign(ad=d.abs()).nlargest(a.top, "ad")
    print(w[["id", "span_old", "val_old", "val_new", "ad"]].to_string(index=False))

    if len(m) > 2:
        print()
        print("diff vs region length (spearman): %.4f"
              % np.corrcoef(m["span_old"].rank(), d.abs().rank())[0, 1])


if __name__ == "__main__":
    main()
