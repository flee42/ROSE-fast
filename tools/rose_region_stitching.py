#!/usr/bin/env python3
"""rose_region_stitching.py -- exact port of ROSE_main.regionStitching()

ROSE does THREE things, not one:
  1. drop constituent loci **contained** within a TSS +/- tssWindow locus
  2. stitchCollection(stitchWindow)
  3. for each stitched locus overlapping >2 DISTINCT gene TSS (TSS +/- 50 bp):
     REVERT it to its original (unstitched) constituent loci

Step 3 is what makes ROSE's region count larger than a plain stitch.
"""
import numpy as np
import pandas as pd


def _load_refseq(path):
    """ROSE annotation/hg*.refseq.ucsc -- one row per transcript.
    cols: name(2) chrom(3) strand(4) txStart(5) txEnd(6) ... name2(13)=gene symbol"""
    out = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            p = line.rstrip(chr(10)).split(chr(9))
            if len(p) < 13:
                continue
            chrom, strand = p[2], p[3]
            try:
                tss = int(p[4]) if strand == "+" else int(p[5])
            except ValueError:
                continue
            out.append((chrom, tss, p[12]))
    return out


def load_genes(path, protein_only=False):
    """dispatch: ROSE refseq .ucsc  |  GENCODE .gtf"""
    if str(path).lower().endswith(".ucsc") or "refseq" in str(path).lower():
        g = _load_refseq(path)
        print("  [annot] refseq: %d transcripts, %d unique gene symbols"
              % (len(g), len({x[2] for x in g})))
        return g
    out = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip(chr(10)).split(chr(9))
            if len(p) < 9 or p[2] != "gene":
                continue
            attr = p[8]
            if protein_only and "protein_coding" not in attr:
                continue
            name = ""
            for kv in attr.split(";"):
                kv = kv.strip()
                if kv.startswith("gene_name "):
                    name = kv.split(" ", 1)[1].strip(chr(34))
                    break
            if not name:
                for kv in attr.split(";"):
                    kv = kv.strip()
                    if kv.startswith("gene_id "):
                        name = kv.split(" ", 1)[1].strip(chr(34))
                        break
            tss = int(p[3]) if p[6] == "+" else int(p[4])
            out.append((p[0], tss, name))
    print("  [annot] gencode: %d genes, %d unique names" % (len(out), len({x[2] for x in out})))
    return out


def _per_chr(genes):
    d = {}
    for c, t, n in genes:
        d.setdefault(c, []).append((t, n))
    for c in d:
        d[c].sort()
    return d


def drop_contained(df, by_chr, tss_window):
    """drop constituents fully inside a TSS +/- window (ROSE step 1).
    pos arrays are precomputed ONCE per chromosome (was O(Nreg x Ntss))."""
    pos_by = {c: np.array([x[0] for x in v], dtype=np.int64) for c, v in by_chr.items()}
    keep = np.ones(len(df), dtype=bool)
    starts = df["start"].values.astype(np.int64)
    ends = df["end"].values.astype(np.int64)
    chroms = df["chr"].values
    for c in np.unique(chroms):
        pos = pos_by.get(c)
        if pos is None or pos.size == 0:
            continue
        idx = np.flatnonzero(chroms == c)
        lo = np.searchsorted(pos, starts[idx] - tss_window, side="left")
        hi = np.searchsorted(pos, ends[idx] + tss_window, side="right")
        for k, a, b in zip(idx, lo, hi):
            for t in pos[a:b]:
                if (t - tss_window) <= starts[k] and ends[k] <= (t + tss_window):
                    keep[k] = False
                    break
    return df[keep].reset_index(drop=True)


def stitch(df, stitch_bp):
    out = []
    for c, sub in df.groupby("chr", sort=False):
        a = sub["start"].values.astype(np.int64); b = sub["end"].values.astype(np.int64)
        o = np.argsort(a); a, b = a[o], b[o]
        s, e = a[0], b[0]
        for i in range(1, len(a)):
            if a[i] - e <= stitch_bp:
                e = max(e, b[i])
            else:
                out.append([c, int(s), int(e)]); s, e = a[i], b[i]
        out.append([c, int(s), int(e)])
    return pd.DataFrame(out, columns=["chr", "start", "end"])


def revert_multi_tss(stitched, original, by_chr, max_genes=2, tss_half=50):
    """ROSE step 3: a stitched locus overlapping >max_genes distinct TSS is
    replaced by the original (unstitched) constituent loci it covers.
    Returns chr/start/end/n_peaks (same schema as stitch())."""
    pos_by = {c: np.array([x[0] for x in v], dtype=np.int64) for c, v in by_chr.items()}
    name_by = {c: [x[1] for x in v] for c, v in by_chr.items()}
    orig_by = {}
    for c, g in original.groupby("chr"):
        g = g.sort_values("start")
        orig_by[c] = (g["start"].values.astype(np.int64),
                      g["end"].values.astype(np.int64),
                      g["start"].values.astype(np.int64))
    fixed = []
    for r in stitched.itertuples():
        pos = pos_by.get(r.chr)
        names = set()
        if pos is not None and pos.size:
            lo = np.searchsorted(pos, r.start - tss_half, side="left")
            hi = np.searchsorted(pos, r.end + tss_half, side="right")
            nm = name_by[r.chr]
            for k in range(lo, hi):
                t = pos[k]
                if (t + tss_half) >= r.start and (t - tss_half) <= r.end:
                    names.add(nm[k])
        og = orig_by.get(r.chr)
        if og is None:
            fixed.append([r.chr, int(r.start), int(r.end), 0]); continue
        o_start, o_end, _ = og
        hi = int(np.searchsorted(o_start, r.end, side="left"))
        if hi == 0:
            n_cov = 0
        else:
            m = o_end[:hi] > r.start
            n_cov = int(m.sum())
        if hi == 0:
            m = np.zeros(0, dtype=bool)
        if len(names) > max_genes:
            for k in np.flatnonzero(m if hi else np.zeros(0, dtype=bool)):
                fixed.append([r.chr, int(o_start[k]), int(o_end[k]), 1])
        else:
            fixed.append([r.chr, int(r.start), int(r.end), n_cov])
    out = pd.DataFrame(fixed, columns=["chr", "start", "end", "n_peaks"])
    return out.drop_duplicates().sort_values(["chr", "start"]).reset_index(drop=True)
