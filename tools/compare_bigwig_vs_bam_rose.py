#!/usr/bin/env python3
"""compare_bigwig_vs_bam_rose.py  (v2: full ROSE regionStitching incl. the
>2-TSS revert step)

usage:
  python3 compare_bigwig_vs_bam_rose.py \
      --rose-outdir rose_out/K562_H3K27me3_t3000 \
      --peaks-bed rose_bam/peaks/K562_H3K27me3_peaks.bed \
      --gtf gencode.v38.annotation.gtf --tss 3000 --stitch 12500 \
      --outdir rose_out/compare_K562
"""
import argparse, glob, os, sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--rose-outdir", required=True)
    p.add_argument("--peaks-bed", required=True)
    p.add_argument("--gtf", default="gencode.v38.annotation.gtf")
    p.add_argument("--tss", type=int, default=3000)
    p.add_argument("--stitch", type=int, default=12500)
    p.add_argument("--protein-only", action="store_true", default=False,
                   help="ROSE uses ALL refseq genes, so default is False")
    p.add_argument("--no-revert", action="store_true", help="skip ROSE step 3")
    p.add_argument("--outdir", default="rose_out/compare")
    return p.parse_args()



def read_bed3(path):
    """minimal BED reader (used to be imported from the deprecated stitching port)"""
    rows = []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip(chr(10)).split(chr(9))
            if len(p) < 3:
                continue
            try:
                rows.append((p[0], int(p[1]), int(p[2])))
            except ValueError:
                continue
    return rows


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from rose_region_stitching import (load_genes, _per_chr, drop_contained,
                                       stitch, revert_multi_tss)

    g = glob.glob(os.path.join(a.rose_outdir, "gff", "*_STITCHED*.gff"))
    if not g:
        raise SystemExit("no *_STITCHED*.gff in %s/gff" % a.rose_outdir)
    rose = pd.read_csv(g[0], sep="\t", header=None)
    rose = rose.rename(columns={0: "chr", 3: "start", 4: "end"})[["chr", "start", "end"]]
    rose = rose.astype({"start": int, "end": int}).sort_values(["chr", "start"]).reset_index(drop=True)

    d0 = read_bed3(a.peaks_bed)
    print("input peaks        :", len(d0))
    genes = load_genes(a.gtf, a.protein_only)
    by = _per_chr(genes)
    print("genes loaded       : %d  (protein_only=%s)" % (len(genes), a.protein_only))

    d = drop_contained(d0, by, a.tss) if a.tss > 0 else d0
    print("after TSS excl     :", len(d))
    st = stitch(d, a.stitch)
    print("stitched (plain)   :", len(st))
    if a.tss > 0 and not a.no_revert:
        st = revert_multi_tss(st, d, by)
        print("after >2-TSS revert:", len(st))
    ours = st.sort_values(["chr", "start"]).reset_index(drop=True)

    print("\nstitched  ours=%d   ROSE=%d" % (len(ours), len(rose)))
    n = min(len(ours), len(rose))
    same = ((ours["chr"].values[:n] == rose["chr"].values[:n]) &
            (ours["start"].values[:n] == rose["start"].values[:n]) &
            (ours["end"].values[:n] == rose["end"].values[:n]))
    # ---- set-based comparison (positions shift when one side has an extra region) ----
    key = lambda d: set(zip(d["chr"], d["start"], d["end"]))
    A, B = key(ours), key(rose)
    inter = len(A & B)
    union = len(A | B)
    print("exact-region intersection : %d" % inter)
    print("Jaccard                   : %.4f  (%d / %d)" % (inter / union if union else 0, inter, union))
    print("only in ours              : %d" % len(A - B))
    print("only in ROSE              : %d" % len(B - A))

    # tolerant match: any ROSE region within +/-tol bp with same chrom
    tol = 1
    Bd = {}
    for c, s, e in B:
        Bd.setdefault(c, []).append((s, e))
    for c in Bd:
        Bd[c].sort()
    hit = 0
    for c, s, e in A:
        ok = False
        for s2, e2 in Bd.get(c, []):
            if abs(s2 - s) <= tol and abs(e2 - e) <= tol:
                ok = True
                break
        hit += ok
    print("tolerant match (+/-%d bp) : %d / %d  (%.2f%%)" % (tol, hit, len(A), 100 * hit / max(len(A), 1)))

    if A - B:
        ex = sorted(A - B)[:3]
        print("  e.g. only in ours :", ex)
    if B - A:
        ex = sorted(B - A)[:3]
        print("  e.g. only in ROSE :", ex)

    pd.DataFrame(dict(chr=ours["chr"].values[:n], start=ours["start"].values[:n],
                      end=ours["end"].values[:n], same=same)).to_csv(
        os.path.join(a.outdir, "stitched_identity.tsv"), sep="\t", index=False)
    print("\nwritten", os.path.join(a.outdir, "stitched_identity.tsv"))


if __name__ == "__main__":
    main()
