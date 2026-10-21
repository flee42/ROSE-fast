#!/usr/bin/env python3
"""rose_call_super_bigwig.py
======================
Rank ROSE-stitched regions and apply ROSE's exact cutoff, using a
fold-change-over-control bigWig in place of read densities.

ROSE signal per stitched region (ROSE_main.py):
    signal = SUM over its constituent peaks of (density_peak * length_peak)
Here `density_peak` is replaced by the mean fold-change over the peak, i.e.
    signal = SUM over constituent peaks of (mean_fold_change * peak_length)
which is proportional to ROSE's read-based signal when the control is locally
uniform and the bigWig is not clipped -> the RANK ORDER, and therefore the
position of the cutoff on the rank curve, is preserved.

Cutoff = ROSE's calculate_cutoff (ROSE_callSuper.R), ported verbatim:
    v <- sort(v); v[v<0] <- 0
    slope <- (max(v) - min(v)) / length(v)
    x* <- argmin_x  #{ i : v[i] <= slope*i + (v[x] - slope*x) }
    cutoff <- v[x*]
("slope = 1" is NOT what the code does.)

Outputs
-------
  {cell}_{mark}_t{tss}_signal.tsv   region, signal, rank, isSuper
  {cell}_{mark}_t{tss}_hockeystich.png
  call_super_summary.tsv            n_regions, n_super, frac_super, cutoff,
                                    cutoff_over_median, x_star/N, n_ties
  call_super.log

usage
-----
  python3 rose_call_super_bigwig.py \
      --stitched-dir rose_stitch_v3/stitched \
      --manifest domwidth_tier12.tsv \
      --bigwig-manifest bigwig_manifest.tsv \
      --marks H3K27ac,H3K27me3 --outdir rose_stitch_v3/super
"""
from __future__ import annotations
import argparse, glob, gzip, os, re
import numpy as np
import pandas as pd

TIER12 = ["CD14_positive_monocyte", "GM12878", "H1", "HeLa-S3", "HepG2",
          "IMR-90", "K562", "MCF-7", "SK_N_SH", "endothelial_cell_of_umbilical_vein"]


def norm(s):
    return "".join(ch for ch in str(s).lower() if ch.isalnum())


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_bed3(path):
    rows = []
    with _open(path) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                rows.append((p[0], int(p[1]), int(p[2])))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["chr", "start", "end"]).drop_duplicates()


def rose_cutoff(v):
    """verbatim port of ROSE_callSuper.R::calculate_cutoff (exhaustive scan)."""
    v = np.sort(np.asarray(v, dtype=float).copy())
    v[v < 0] = 0.0
    n = len(v)
    if n < 3:
        return np.nan, np.nan, np.nan
    slope = (v.max() - v.min()) / n
    xs = np.arange(1, n + 1, dtype=float)
    best_x, best_obj, n_ties = 1, None, 0
    for x in range(1, n + 1):
        b = v[x - 1] - slope * x
        obj = int(np.sum(v <= xs * slope + b))
        if best_obj is None or obj < best_obj:
            best_obj, best_x, n_ties = obj, x, 1
        elif obj == best_obj:
            n_ties += 1
    return slope, best_x, v[best_x - 1]


def peak_signals(regions, peaks, bw, log=None):
    """SUM over constituent peaks of (mean_fold_change * peak_length).

    Handles chromosome-name mismatches (chr1 vs 1) and logs the first failure
    instead of silently returning NaN.
    """
    bwchroms = bw.chroms()
    def bwname(c):
        if c in bwchroms:
            return c
        alt = c[3:] if str(c).lower().startswith("chr") else "chr" + str(c)
        return alt if alt in bwchroms else None

    out = np.full(len(regions), np.nan)
    warned = [False]

    def one(c, s, e):
        """mean fold-change over [s,e) -> float or nan"""
        cn = bwname(c)
        if cn is None:
            if not warned[0] and log:
                log("  !! chrom %r not in bigWig (bw has e.g. %s)" %
                    (c, list(bwchroms)[:4])); warned[0] = True
            return np.nan
        try:
            v = bw.stats(cn, int(s), int(e), type="mean")
        except Exception as ex:
            if not warned[0] and log:
                log("  !! bw.stats failed: %s %s-%s : %r" % (cn, s, e, ex)); warned[0] = True
            return np.nan
        if not v or v[0] is None:
            return np.nan
        return float(v[0])

    for c, rsub in regions.groupby("chr"):
        psub = peaks[peaks["chr"] == c]
        if not len(psub):
            continue
        ps = np.sort(psub["start"].values.astype(np.int64))
        pe = psub["end"].values.astype(np.int64)[np.argsort(psub["start"].values)]
        for idx, r in rsub.iterrows():
            i0 = np.searchsorted(ps, r["start"], side="left")
            i1 = np.searchsorted(ps, r["end"], side="right")
            pos = regions.index.get_loc(idx)
            if i1 <= i0:
                out[pos] = 0.0
                continue
            sig, tot = 0.0, 0
            for k in range(i0, i1):
                v = one(c, ps[k], pe[k])
                if np.isfinite(v):
                    sig += v * (pe[k] - ps[k]); tot += 1
            out[pos] = sig if tot else np.nan
    return out


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--stitched-dir", default="rose_stitch_v3/stitched")
    p.add_argument("--manifest", default="domwidth_tier12.tsv", help="cell, mark, bed (constituent peaks)")
    p.add_argument("--bigwig-manifest", required=True, help="cell, mark, bigwig (fold change over control)")
    p.add_argument("--marks", default="H3K27ac,H3K27me3")
    p.add_argument("--outdir", default="rose_stitch_v3/super")
    p.add_argument("--no-plot", action="store_true")
    return p.parse_args()


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    logf = open(os.path.join(a.outdir, "call_super.log"), "a")
    def log(*x):
        s = " ".join(str(i) for i in x)
        print(s); logf.write(s + "\n"); logf.flush()

    import pyBigWig
    peaks_man = pd.read_csv(a.manifest, sep="\t")
    bw_man = pd.read_csv(a.bigwig_manifest, sep="\t")
    bw_man["k"] = bw_man["cell"].astype(str).map(norm) + "|" + bw_man["mark"].astype(str)
    bwmap = dict(zip(bw_man["k"], bw_man["bigwig"]))
    marks = [m for m in a.marks.split(",") if m.strip()]
    files = sorted(glob.glob(os.path.join(a.stitched_dir, "*.bed.gz")))
    if not files:
        raise SystemExit("no stitched bed under %s" % a.stitched_dir)

    rows = []
    for f in files:
        m = re.match(r"(.+)_(H3K\w+)_t(\d+)_s(\d+)\.bed\.gz$", os.path.basename(f))
        if not m:
            continue
        cell, mark, tb, sb = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        if mark not in marks:
            continue
        bwpath = bwmap.get(norm(cell) + "|" + mark)
        if bwpath is None or not os.path.exists(bwpath):
            log("  !! no bigWig for %s %s" % (cell, mark)); continue
        phit = peaks_man[(peaks_man["mark"] == mark) &
                         (peaks_man["cell"].astype(str).map(norm) == norm(cell))]
        if not len(phit):
            log("  !! no peak file for %s %s" % (cell, mark)); continue
        peaks = read_bed3(str(phit.iloc[0]["bed"]))
        regions = pd.read_csv(f, sep="\t")
        bw = pyBigWig.open(bwpath)
        sig = peak_signals(regions, peaks, bw, log=log)
        bw.close()
        ok = np.isfinite(sig)
        if ok.sum() < 10:
            log("  !! %s %s: too few finite signals (%d/%d). DIAGNOSTICS:" % (cell, mark, int(ok.sum()), len(sig)))
            log("     bigWig: %s" % bwpath)
            _bw = pyBigWig.open(bwpath)
            log("     bw chroms (%d): %s" % (len(_bw.chroms()), list(_bw.chroms())[:5]))
            log("     peak chroms (%d): %s" % (peaks["chr"].nunique(), list(pd.unique(peaks["chr"]))[:5]))
            log("     peak count: %d | regions: %d" % (len(peaks), len(regions)))
            r0 = regions.iloc[0]
            log("     first region: %s:%s-%s | bw.stats ->" % (r0["chr"], r0["start"], r0["end"]))
            try:
                log("       %r" % (_bw.stats(str(r0["chr"]), int(r0["start"]), int(r0["end"]), type="mean"),))
            except Exception as ex:
                log("       EXCEPTION: %r" % (ex,))
            _bw.close()
            continue
        slope, xstar, cut = rose_cutoff(sig[ok])
        order = np.argsort(-np.nan_to_num(sig, nan=-1))
        rank = np.empty(len(sig), dtype=int); rank[order] = np.arange(1, len(sig) + 1)
        out = regions.copy()
        out["signal"] = sig
        out["rank"] = rank
        out["isSuper"] = (sig > cut).astype(int)
        stem = "%s_%s_t%d" % (cell, mark, tb)
        out.to_csv(os.path.join(a.outdir, stem + "_signal.tsv"), sep="\t", index=False)
        med = np.median(sig[ok])
        rows.append(dict(cell=cell, mark=mark, tss_bp=tb, stitch_bp=sb,
                         n_regions=int(ok.sum()), n_super=int((sig[ok] > cut).sum()),
                         frac_super=float((sig[ok] > cut).mean()),
                         cutoff=cut, cutoff_over_median=cut / med if med else np.nan,
                         x_star_frac=xstar / ok.sum(), n_ties=int(0),
                         n_constituent_peaks=len(peaks)))
        log("%-38s %-9s t=%4d  regions %5d  super %4d (%.1f%%)  cutoff %.3g (%.1fx median)"
            % (cell, mark, tb, ok.sum(), int((sig[ok] > cut).sum()),
               100 * (sig[ok] > cut).mean(), cut, cut / med if med else np.nan))
        if not a.no_plot:
            try:
                import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(5, 4.2))
                v = np.sort(sig[ok])[::-1]
                ax.set_xlabel("rank"); ax.set_ylabel("signal (mean fold-change x peak bp)")
                ax.set_title("%s %s\ncutoff %.3g (%.1fx median), %d/%d super"
                             % (cell, mark, cut, cut / med if med else np.nan,
                                int((sig[ok] > cut).sum()), ok.sum()), fontsize=8)
                fig.tight_layout()
                fig.savefig(os.path.join(a.outdir, stem + "_hockeystick.png"), dpi=140)
                plt.close(fig)
            except Exception as e:
                log("  !! plot failed: %s" % e)
    if not rows:
        raise SystemExit("nothing produced")
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(a.outdir, "call_super_summary.tsv"), sep="\t", index=False)
    log("\n== median over cells ==")
    log(R.groupby(["mark", "tss_bp"]).agg(
        n_cells=("cell", "nunique"), n_regions=("n_regions", "median"),
        n_super=("n_super", "median"), frac_super=("frac_super", "median"),
        cutoff_over_median=("cutoff_over_median", "median")).round(4).to_string())
    logf.close()


if __name__ == "__main__":
    main()
