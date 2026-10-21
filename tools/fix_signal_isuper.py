#!/usr/bin/env python3
"""fix_signal_isuper.py
Rebuild the isSuper column of *_signal.tsv by exact coordinate matching against
ROSE's own *_SuperStitched.table.txt.  Nothing else in the file is touched.

This is deliberately independent of REGION_ID parsing: R's AllStitched table is
re-sorted by stitchedPeakRank, and any positional/ID-based back-mapping is fragile.

usage: python3 fix_signal_isuper.py rose_bw_v3/super_R
"""
import glob, os, sys
import numpy as np
import pandas as pd

TAB = chr(9)


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "rose_bw_v3/super_R"
    n_ok = n_bad = n_skip = 0
    rows = []
    for si in sorted(glob.glob(os.path.join(root, "*_signal.tsv"))):
        base = si[: -len("_signal.tsv")]
        sp = base + "_SuperStitched.table.txt"
        if not os.path.exists(sp):
            n_skip += 1
            continue
        s = pd.read_csv(sp, sep=TAB, comment="#")
        s.columns = [str(c).strip() for c in s.columns]
        keys = set(zip(s["CHROM"].astype(str),
                       s["START"].astype(int),
                       s["STOP"].astype(int)))
        d = pd.read_csv(si, sep=TAB)
        new = np.array([1 if (str(c), int(a), int(b)) in keys else 0
                        for c, a, b in zip(d["chr"], d["start"], d["end"])], dtype=int)
        old = int(d["isSuper"].sum()) if "isSuper" in d.columns else -1
        hit = int((new == 1).sum())
        if hit != len(keys):
            print("  !! %s: matched %d of %d super regions -- NOT written"
                  % (os.path.basename(base), hit, len(keys)))
            n_bad += 1
            continue
        d["isSuper"] = new
        d.to_csv(si, sep=TAB, index=False)
        n_ok += 1
        rows.append(dict(sample=os.path.basename(base), n=len(d),
                         n_super_old=old, n_super_new=int(new.sum()),
                         overlap_with_old=int(((new == 1) & (d.get("isSuper", new) == new)).sum())
                         if old >= 0 else -1))
    T = pd.DataFrame(rows)
    if len(T):
        pd.set_option("display.width", 200); pd.set_option("display.max_rows", 100)
        print(T.to_string(index=False))
    print()
    print("rewritten=%d  unmatched=%d  no_super_table=%d" % (n_ok, n_bad, n_skip))


if __name__ == "__main__":
    main()
