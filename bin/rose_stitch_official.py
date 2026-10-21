#!/usr/bin/env python3
"""rose_stitch_official.py
Stitch constituent peaks with ROSE's OWN code (ROSE_utils.bedToGFF +
ROSE_main.regionStitching) and emit a BED for rose_call_super_v3.py.

Handles gzipped input peaks (ROSE_utils.parseTable cannot read gzip) by
materialising a plain BED once per sample.

usage:
  python3 rose_stitch_official.py --manifest rose_bw_all/manifest.tsv \
      --marks H3K27me3 --annot ROSE/annotation/hg38_refseq.ucsc \
      --tss 3000 --stitch 12500 --outdir rose_bw_all/stitch
"""
import argparse, gzip, os, sys
import numpy as np
import pandas as pd

TAB = chr(9)
NL = chr(10)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--rose-dir", default="ROSE")
    p.add_argument("--manifest", required=True)
    p.add_argument("--marks", default="H3K27me3")
    p.add_argument("--annot", default="ROSE/annotation/hg38_refseq.ucsc")
    p.add_argument("--tss", type=int, default=3000)
    p.add_argument("--stitch", type=int, default=12500)
    p.add_argument("--outdir", default="rose_bw_all/stitch")
    p.add_argument("--resume", action="store_true",
                   help="skip a sample whose stitched BED already exists")
    return p.parse_args()


def materialise_plain(bed, plain_dir):
    """decompress a .gz peak file to a plain 6-column BED (cached)"""
    if not str(bed).endswith(".gz"):
        return bed
    os.makedirs(plain_dir, exist_ok=True)
    base = os.path.basename(bed)
    for ext in (".bed.gz", ".narrowPeak.gz", ".broadPeak.gz"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    dst = os.path.join(plain_dir, base + ".bed")
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        return dst
    with gzip.open(bed, "rt", errors="ignore") as fh, open(dst, "w", newline=NL) as o:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip(NL).split(TAB)
            if len(f) < 3:
                continue
            o.write(TAB.join(f[:6]) + NL)
    return dst


def main():
    a = parse_args()
    rose = os.path.abspath(a.rose_dir)
    sys.path.insert(0, os.path.join(rose, "lib"))
    sys.path.insert(0, os.path.join(rose, "bin"))
    import ROSE_utils
    import ROSE_main

    out = os.path.abspath(a.outdir)
    gdir = os.path.join(out, "gff")
    bdir = os.path.join(out, "stitched")
    pdir = os.path.join(out, "peaks_plain")
    for d in (gdir, bdir, pdir):
        os.makedirs(d, exist_ok=True)

    man = pd.read_csv(a.manifest, sep=TAB)
    marks = [m.strip() for m in a.marks.split(",") if m.strip()]
    print("%-38s %-10s %8s %8s %8s" % ("cell", "mark", "peaks", "kept", "regions"), flush=True)
    rows = []
    for r in man.itertuples():
        if r.mark not in marks:
            continue
        bed = str(r.bed)
        if not os.path.exists(bed):
            print("  !! missing %s" % bed, flush=True)
            continue
        dst = os.path.join(bdir, "%s_%s_t%d_s%d.bed.gz" % (r.cell, r.mark, a.tss, a.stitch))
        if a.resume and os.path.exists(dst) and os.path.getsize(dst) > 0:
            print("%-38s %-10s  skip (exists)" % (r.cell, r.mark), flush=True)
            continue
        bed = materialise_plain(bed, pdir)
        gff = os.path.join(gdir, "%s_%s.gff" % (r.cell, r.mark))
        if not os.path.exists(gff):
            ROSE_utils.bedToGFF(bed, gff)
        n_peak = sum(1 for _ in open(gff))
        stitched, _ = ROSE_main.regionStitching(gff, a.stitch, a.tss, a.annot,
                                                removeTSS=bool(a.tss))
        g = ROSE_utils.locusCollectionToGFF(stitched)
        d = pd.DataFrame({"chr": [x[0] for x in g],
                          "start": [int(x[3]) for x in g],
                          "end": [int(x[4]) for x in g]})
        pk = pd.read_csv(gff, sep=TAB, header=None, usecols=[0, 3, 4],
                         names=["chr", "start", "end"])
        d = d.sort_values(["chr", "start"]).reset_index(drop=True)
        parts = []
        for c, sub in pk.groupby("chr"):
            s = np.sort(sub["start"].values.astype(np.int64))
            e = sub["end"].values.astype(np.int64)[np.argsort(sub["start"].values)]
            rr = d[d.chr == c]
            i0 = np.searchsorted(s, rr["start"].values, side="left")
            i1 = np.searchsorted(s, rr["end"].values, side="right")
            parts.append(pd.Series(i1 - i0, index=rr.index))
        d["n_peaks"] = pd.concat(parts).sort_index() if parts else 0
        d["span_bp"] = d["end"] - d["start"]
        d.to_csv(dst, sep=TAB, index=False, compression="gzip")
        print("%-38s %-10s %8d %8d %8d" % (r.cell, r.mark, n_peak, len(pk), len(d)), flush=True)
        rows.append(dict(cell=r.cell, mark=r.mark, n_peaks_init=n_peak,
                         n_peaks=len(pk), n_regions=len(d),
                         median_span_kb=d.span_bp.median() / 1000,
                         total_mb=d.span_bp.sum() / 1e6))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(out, "stitch_summary.tsv"), sep=TAB, index=False)
    if len(T):
        print(flush=True)
        print("=== median per mark ===", flush=True)
        print(T.groupby("mark")[["n_regions", "median_span_kb", "total_mb"]].median().round(2).to_string(), flush=True)
        print(flush=True)
    print("written", bdir, flush=True)


if __name__ == "__main__":
    main()
