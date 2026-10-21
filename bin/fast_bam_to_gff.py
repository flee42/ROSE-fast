#!/usr/bin/env python3
"""fast_bam_to_gff.py -- drop-in replacement for ROSE_bamToGFF.py (same CLI,
same output), but O(genome) instead of O(regions x bases), with chunk-level
multiprocessing.

Why ROSE's version is slow: for EVERY input region it does a pysam fetch, builds
a per-base python dict (`for x in range(start, end+1): h[x]+=1`) and then scans
all keys to sum inside the region -> cost ~ NREG x (fetch + per-base dict ops).

This version: read each genomic chunk once into numpy arrays, build a
difference array, cumsum -> exact per-base coverage, then O(1) per region.

usage:
  fast_bam_to_gff.py -f 1 -e 200 -r -m 1 -w 8 -b x.bam -i in.gff -o out.gff
"""
import argparse, sys
from collections import defaultdict
import numpy as np

try:
    import pysam
except ImportError:
    pysam = None


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-b", "--bam", required=True)
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("-f", "--floor", type=int, default=0)
    p.add_argument("-e", "--extension", type=int, default=200)
    p.add_argument("-r", "--rpm", action="store_true")
    p.add_argument("-m", "--matrix", type=int, default=1)
    p.add_argument("-s", "--sense", default="both")
    p.add_argument("-w", "--workers", type=int, default=4)
    p.add_argument("--chunk-mb", type=float, default=5.0)
    return p.parse_args()


def read_gff(path):
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0].startswith("#") or p[0] == "GENE_ID":
                continue
            try:
                rows.append((p[1], p[0], int(p[3]), int(p[4]), p[6]))
            except ValueError:
                continue
    return rows


def _chunk_worker(task):
    """coverage for one (chrom, window) chunk; returns {region_index: density}"""
    (bam_path, chrom, wstart, wend, regs, ext, sense, floor, mmr) = task
    import pysam as _ps
    bam = _ps.AlignmentFile(bam_path, "rb")
    wp = sense in ("+", "both", ".")
    wm = sense in ("-", "both", ".")
    st, en = [], []
    for r in bam.fetch(chrom, max(0, wstart - ext), wend + ext):   # ROSE searchLocus = region +/- ext
        if r.is_unmapped:
            continue
        if r.is_reverse:
            if not wm:
                continue
            s = max(wstart, r.reference_start - ext); e = min(wend, r.reference_end)
        else:
            if not wp:
                continue
            s = max(wstart, r.reference_start); e = min(wend, r.reference_end + ext)
        if e > s:
            st.append(s); en.append(e)
    bam.close()
    if not st:
        return {idx: 0.0 for idx, _, _ in regs}
    st = np.asarray(st, dtype=np.int64)
    en = np.asarray(en, dtype=np.int64)
    st = np.clip(st, wstart, wend - 1)
    en = np.clip(en, wstart + 1, wend)
    n = wend - wstart + 3
    d = np.zeros(n, dtype=np.int32)
    np.add.at(d, st - wstart, 1)
    np.add.at(d, en - wstart + 1, -1)   # ROSE uses range(start, end+1): inclusive end
    cov = np.cumsum(d, dtype=np.int32)
    if floor > 0:
        cov = np.where(cov > floor, cov, 0)
    cum = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(cov, out=cum[1:])
    out = {}
    for idx, s, e in regs:
        ln = e - s
        if ln <= 1:
            out[idx] = 0.0; continue
        # ROSE: keys are s < x < e-1 , denominator = len-1
        a = max(s + 1, wstart) - wstart
        b = min(e - 1, wend) - wstart
        if b <= a:
            out[idx] = 0.0; continue
        out[idx] = ((int(cum[b] - cum[a]) / (ln - 1)) / mmr)

    return out


def main():
    a = parse_args()
    if pysam is None:
        sys.exit("pysam is required")
    if a.matrix != 1:
        sys.exit("only -m 1 is supported (that is what ROSE_main.py uses)")

    regions = read_gff(a.input)
    by_chr = defaultdict(list)
    for k, (rid, chrom, s, e, strand) in enumerate(regions):
        by_chr[chrom].append((k, s, e))

    bam = pysam.AlignmentFile(a.bam, "rb")
    refs = set(bam.references)
    has_chr = any(r.startswith("chr") for r in refs)

    def canon(c):
        if has_chr:
            return c if c.startswith("chr") else "chr" + c
        return c.replace("chr", "", 1) if c.startswith("chr") else c

    mmr = 1.0
    if a.rpm:
        tot = 0
        for line in pysam.idxstats(a.bam).splitlines():
            p = line.split("\t")
            if len(p) >= 3:
                try:
                    tot += int(p[2])
                except ValueError:
                    pass
        mmr = max(tot, 1) / 1e6

    chunks = int(max(1.0, a.chunk_mb) * 1e6)
    tasks = []
    for chrom, regs in by_chr.items():
        ch = canon(chrom)
        if ch not in refs:
            continue
        clen = bam.get_reference_length(ch)
        buckets = {}
        for k, s, e in regs:
            buckets.setdefault(s // chunks, []).append((k, s, e))
        for bidx, lst in buckets.items():
            wstart = bidx * chunks
            wend = min(clen, max(e for _, _, e in lst) + a.extension + 1)
            tasks.append((a.bam, ch, wstart, wend, lst, a.extension, a.sense, a.floor, mmr))
    bam.close()
    print("chunks: %d  regions: %d  workers: %d" % (len(tasks), len(regions), a.workers),
          file=sys.stderr)

    sig = {}
    def _run_serial():
        s = {}
        for k, t in enumerate(tasks):
            s.update(_chunk_worker(t))
            if (k + 1) % 100 == 0:
                print("  %d/%d chunks" % (k + 1, len(tasks)), file=sys.stderr)
        return s

    if a.workers and a.workers > 1:
        try:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(max_workers=a.workers) as pool:
                for k, res in enumerate(pool.map(_chunk_worker, tasks)):
                    sig.update(res)
                    if (k + 1) % 100 == 0:
                        print("  %d/%d chunks" % (k + 1, len(tasks)), file=sys.stderr)
        except Exception as e:
            print("!! process pool failed (%s: %s) -- falling back to SERIAL"
                  % (type(e).__name__, e), file=sys.stderr)
            sig = _run_serial()
    else:
        for k, t in enumerate(tasks):
            sig.update(_chunk_worker(t))
            if (k + 1) % 100 == 0:
                print("  %d/%d chunks" % (k + 1, len(tasks)), file=sys.stderr)

    bamfn = a.bam.replace("\\", "/").split("/")[-1]
    with open(a.output, "w", newline="\n") as fh:
        fh.write("GENE_ID\tlocusLine\tbin_1_%s\n" % bamfn)
        for k, (rid, chrom, s, e, strand) in enumerate(regions):
            fh.write("%s\t%s(%s):%d-%d\t%.4f\n" % (rid, chrom, strand, s, e, sig.get(k, 0.0)))
    print("written %s (%d regions)" % (a.output, len(regions)), file=sys.stderr)


if __name__ == "__main__":
    main()
