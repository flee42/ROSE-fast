#!/usr/bin/env python3
"""rose_call_super_R_native.py
Use ROSE's ORIGINAL R script (ROSE_callSuper.R) for the cutoff and super calling,
with a fold-change bigWig supplying the signal instead of read densities.

Build a REGION_MAP with 7 columns
    REGION_ID  CHROM  START  STOP  NUM_LOCI  CONSTITUENT_SIZE  <bigwig signal>
and call:  Rscript ROSE_callSuper.R <outdir> <region_map> <enhancerName> NONE
(wceName='NONE' -> R uses column 7 directly, no control subtraction).

If an existing *_signal.tsv is present its `signal` column is reused (--reuse-signal),
so the bigWig does not have to be rescanned.

usage:
  python3 rose_call_super_R_native.py \
      --stitched-dir rose_bw_v3/stitch/stitched \
      --manifest rose_bam/peaks_matched/manifest.tsv \
      --bigwig-manifest rose_bam/peaks_matched/bigwig_manifest.tsv \
      --marks H3K27me3 --outdir rose_bw_v3/super_R \
      --signal-from rose_bw_v3/super --reuse-signal
"""
import argparse, glob, gzip, os, re, subprocess, sys
import numpy as np
import pandas as pd

TAB = chr(9)
NL = chr(10)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--rose-dir", default="ROSE")
    p.add_argument("--stitched-dir", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--bigwig-manifest", default="")
    p.add_argument("--marks", default="H3K27me3,H3K27ac,H3K4me1")
    p.add_argument("--outdir", default="rose_super_R")
    p.add_argument("--reuse-signal", action="store_true")
    p.add_argument("--signal-from", default="")
    p.add_argument("--no-plot", action="store_true")
    return p.parse_args()


def read_bed3(path):
    rows = []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip(NL).split(TAB)
            if len(p) < 3:
                continue
            try:
                rows.append((p[0], int(p[1]), int(p[2])))
            except ValueError:
                continue
    return rows


def score_with_bigwig(regions, peaks, bw_path):
    """sum over constituent peaks of (mean fold-change * peak length)"""
    import pyBigWig
    bw = pyBigWig.open(bw_path)
    bwc = set(bw.chroms())

    def one(c, s, e):
        cn = c if c in bwc else (c[3:] if str(c).lower().startswith("chr") else "chr" + str(c))
        if cn not in bwc:
            return np.nan
        v = bw.stats(cn, int(s), int(e), type="mean")
        return float(v[0]) if v and v[0] is not None else np.nan

    out = np.full(len(regions), np.nan)
    idx_of = {k: i for i, k in enumerate(regions.index)}
    for c, rsub in regions.groupby("chr"):
        psub = peaks[peaks["chr"] == c]
        if not len(psub):
            continue
        ps = np.sort(psub["start"].values.astype(np.int64))
        pe = psub["end"].values.astype(np.int64)[np.argsort(psub["start"].values)]
        for idx, r in rsub.iterrows():
            i0 = np.searchsorted(ps, r["start"], side="left")
            i1 = np.searchsorted(ps, r["end"], side="right")
            if i1 <= i0:
                out[idx_of[idx]] = 0.0
                continue
            sig, tot = 0.0, 0
            for k in range(i0, i1):
                v = one(c, ps[k], pe[k])
                if np.isfinite(v):
                    sig += v * (pe[k] - ps[k])
                    tot += 1
            out[idx_of[idx]] = sig if tot else np.nan
    bw.close()
    return out


def main():
    a = parse_args()
    rose = os.path.abspath(a.rose_dir)
    out = os.path.abspath(a.outdir)
    os.makedirs(out, exist_ok=True)
    rscript = os.path.join(rose, "bin", "ROSE_callSuper.R")
    if not os.path.exists(rscript):
        raise SystemExit("ROSE_callSuper.R not found at " + rscript)

    man = pd.read_csv(a.manifest, sep=TAB)
    man["_k"] = man["cell"].astype(str).str.lower() + "|" + man["mark"].astype(str)
    bwmap = {}
    if a.bigwig_manifest and os.path.exists(a.bigwig_manifest):
        b = pd.read_csv(a.bigwig_manifest, sep=TAB)
        b["_k"] = b["cell"].astype(str).str.lower() + "|" + b["mark"].astype(str)
        bwmap = dict(zip(b["_k"], b["bigwig"]))
    sigdir = a.signal_from or out
    marks = [m.strip() for m in a.marks.split(",") if m.strip()]

    print("%-38s %-10s %8s %10s %8s %s"
          % ("sample", "mark", "regions", "cutoff", "n_super", "src"), flush=True)
    rows = []
    for f in sorted(glob.glob(os.path.join(a.stitched_dir, "*.bed.gz"))):
        m = re.match(r"(.+)_(H3K\w+)_t(\d+)_s(\d+)\.bed\.gz$", os.path.basename(f))
        if not m:
            continue
        cell, mark, tb = m.group(1), m.group(2), int(m.group(3))
        if mark not in marks:
            continue
        regions = pd.read_csv(f, sep=TAB)
        base = "%s_%s_t%d" % (cell, mark, tb)
        sig, src = None, ""
        if a.reuse_signal:
            cand = os.path.join(sigdir, base + "_signal.tsv")
            if os.path.exists(cand):
                prev = pd.read_csv(cand, sep=TAB)
                if "signal" in prev.columns and len(prev) == len(regions):
                    sig = prev["signal"].values
                    src = "reused"
        if sig is None:
            key = cell.lower() + "|" + mark
            bw = bwmap.get(key)
            if not bw or not os.path.exists(bw):
                print("  !! no bigWig for %s %s" % (cell, mark), flush=True)
                continue
            phit = man[man["_k"] == key]
            if not phit.empty:
                peaks = pd.DataFrame(read_bed3(str(phit.iloc[0]["bed"])),
                                     columns=["chr", "start", "end"])
            else:
                peaks = regions[["chr", "start", "end"]]
            sig = score_with_bigwig(regions, peaks, bw)
            src = "bigWig"
        if np.isfinite(sig).sum() < 10:
            print("  !! %s %s: too few finite signals" % (cell, mark), flush=True)
            continue

        sigcol = os.path.basename(bwmap.get(cell.lower() + "|" + mark, "bigwig"))
        rm = pd.DataFrame({
            "REGION_ID": ["%s_%d" % (base, i) for i in range(len(regions))],
            "CHROM": regions["chr"].values,
            "START": regions["start"].values,
            "STOP": regions["end"].values,
            "NUM_LOCI": regions["n_peaks"].values if "n_peaks" in regions.columns else 1,
            "CONSTITUENT_SIZE": (regions["end"] - regions["start"]).values,
            sigcol: [("NA" if not np.isfinite(v) else v) for v in sig],
        })
        rm_path = os.path.join(out, base + "_REGION_MAP.txt")
        rm.to_csv(rm_path, sep=TAB, index=False)

        cmd = ["Rscript", rscript, out + os.sep, rm_path, base, "NONE"]
        rc = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if rc != 0:
            print("  !! R failed rc=%d for %s" % (rc, base), flush=True)
            continue

        al = os.path.join(out, base + "_AllStitched.table.txt")
        isSuper = np.zeros(len(regions), dtype=int)
        cut = np.nan
        if os.path.exists(al):
            with open(al, encoding="utf-8", errors="ignore") as fh:
                head = "".join(fh.readline() for _ in range(10))
            mm = re.search(r"cutoff of\s+([0-9.eE+-]+)", head)
            if mm:
                cut = float(mm.group(1))
            t = pd.read_csv(al, sep=TAB, comment="#")
            t.columns = [str(c).strip() for c in t.columns]
            # ROSE_callSuper.R sorts AllStitched by stitchedPeakRank -> we MUST map
            # isSuper back by REGION_ID, not by row position.
        # Exact coordinate match against ROSE's own SuperStitched table -- immune to
        # the stitchedPeakRank re-sort inside writeSuperEnhancer_table.
        _sp = os.path.join(out, base + "_SuperStitched.table.txt")
        if os.path.exists(_sp):
            _s = pd.read_csv(_sp, sep=TAB, comment="#")
            _s.columns = [str(c).strip() for c in _s.columns]
            _keys = set(zip(_s["CHROM"].astype(str), _s["START"].astype(int),
                            _s["STOP"].astype(int)))
            isSuper = np.array([1 if (str(c), int(a), int(b)) in _keys else 0
                                for c, a, b in zip(regions["chr"], regions["start"],
                                                   regions["end"])], dtype=int)

        order = np.argsort(-np.nan_to_num(sig, nan=-1))
        rank = np.empty(len(sig), dtype=int)
        rank[order] = np.arange(1, len(sig) + 1)
        sdf = regions.copy()
        sdf["signal"] = sig
        sdf["rank"] = rank
        sdf["isSuper"] = isSuper
        sdf.to_csv(os.path.join(out, base + "_signal.tsv"), sep=TAB, index=False)

        sp = os.path.join(out, base + "_SuperStitched.table.txt")
        ns = int(isSuper.sum())
        if os.path.exists(sp):
            with open(sp, encoding="utf-8", errors="ignore") as fh:
                nrow_sp = sum(1 for ln in fh if ln.strip() and not ln.startswith("#")) - 1
            if nrow_sp >= 0 and nrow_sp != ns:
                print("  !! %s: isSuper count %d != SuperStitched rows %d -- SKIP write"
                      % (base, ns, nrow_sp), flush=True)
                continue
        print("%-38s %-10s %8d %10.1f %8d %s"
              % (cell, mark, len(regions), cut, ns, src), flush=True)
        rows.append(dict(cell=cell, mark=mark, n_regions=len(regions), cutoff=cut,
                         n_super=ns, frac_super=ns / len(regions), src=src))

    if rows:
        T = pd.DataFrame(rows)
        T.to_csv(os.path.join(out, "call_super_R_summary.tsv"), sep=TAB, index=False)
        print()
        print(T.groupby("mark")[["n_regions", "n_super", "frac_super"]]
              .median().round(4).to_string(), flush=True)
    print()
    print("written", out, flush=True)


if __name__ == "__main__":
    main()
