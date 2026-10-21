#!/usr/bin/env python3
"""Unified runner for ROSE fast paths. Original ROSE_main.py is unchanged."""
import argparse
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent
ROSE = BIN.parent
DEFAULT_ANNOT = ROSE / "annotation" / "hg38_refseq.ucsc"


def call(cmd):
    """Run a command and exit on failure."""
    print("+", " ".join(str(x) for x in cmd), flush=True)
    rc = subprocess.call([str(x) for x in cmd])
    if rc:
        sys.exit(rc)


def add_rose_dir(p):
    p.add_argument("--rose-dir", default=str(ROSE))


def add_stitch_args(p):
    p.add_argument("--tss", type=int, default=3000, help="TSS window; 0 disables")
    p.add_argument("--stitch", type=int, default=12500, help="stitch window (bp)")


def cmd_bam(a):
    cmd = [
        sys.executable, BIN / "rose_run_one.py",
        "--rose-dir", a.rose_dir,
        "--bed", a.bed,
        "--rank-bam", a.rank_bam,
        "--outdir", a.outdir,
        "--stitch", a.stitch,
        "--tss", a.tss,
        "--gff-workers", a.gff_workers,
    ]
    if a.control_bam:
        cmd += ["--control-bam", a.control_bam]
    if a.annot:
        cmd += ["--annot", a.annot]
    if a.slow_gff:
        cmd.append("--slow-gff")
    if a.resume:
        cmd.append("--resume")
    call(cmd)


def cmd_stitch(a):
    annot = Path(a.annot).resolve() if a.annot else Path(a.rose_dir) / "annotation" / "hg38_refseq.ucsc"
    call([
        sys.executable, BIN / "rose_stitch_official.py",
        "--rose-dir", a.rose_dir,
        "--manifest", a.manifest,
        "--marks", a.marks,
        "--tss", a.tss,
        "--stitch", a.stitch,
        "--outdir", a.outdir,
        "--annot", annot,
    ])


def cmd_score(a):
    cmd = [
        sys.executable, BIN / "rose_call_super_R_native.py",
        "--stitched-dir", a.stitched_dir,
        "--manifest", a.manifest,
        "--bigwig-manifest", a.bigwig_manifest,
        "--marks", a.marks,
        "--outdir", a.outdir,
    ]
    if a.no_plot:
        cmd.append("--no-plot")
    call(cmd)


def cmd_bigwig(a):
    bed = Path(a.bed).resolve()
    bw = Path(a.bigwig).resolve()
    annot = Path(a.annot).resolve() if a.annot else DEFAULT_ANNOT
    if not bed.is_file():
        sys.exit("missing peak BED: %s" % bed)
    if not bw.is_file():
        sys.exit("missing bigWig: %s" % bw)
    if not annot.is_file():
        sys.exit("missing annotation: %s" % annot)

    out = Path(a.outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    peak_man = out / "peaks.tsv"
    bw_man = out / "bigwigs.tsv"
    with open(peak_man, "w", newline="\n") as fh:
        fh.write("cell\tmark\tbed\n%s\t%s\t%s\n" % (a.cell, a.mark, bed))
    with open(bw_man, "w", newline="\n") as fh:
        fh.write("cell\tmark\tbigwig\n%s\t%s\t%s\n" % (a.cell, a.mark, bw))

    stitch_dir = out / "stitch"
    super_dir = out / "super"
    call([
        sys.executable, BIN / "rose_stitch_official.py",
        "--rose-dir", a.rose_dir,
        "--manifest", peak_man,
        "--marks", a.mark,
        "--annot", annot,
        "--tss", a.tss,
        "--stitch", a.stitch,
        "--outdir", stitch_dir,
    ])
    # cutoff + super calling go through ROSE's ORIGINAL ROSE_callSuper.R
    # (7-column REGION_MAP, wceName='NONE') -> same code path as the BAM route.
    cmd = [
        sys.executable, BIN / "rose_call_super_R_native.py",
        "--rose-dir", a.rose_dir,
        "--stitched-dir", stitch_dir / "stitched",
        "--manifest", peak_man,
        "--bigwig-manifest", bw_man,
        "--marks", a.mark,
        "--outdir", super_dir,
    ]
    if a.no_plot:
        cmd.append("--no-plot")
    call(cmd)
    print("done:", super_dir)


def main():
    p = argparse.ArgumentParser(description="ROSE fast runner")
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("bam", help="fast BAM ROSE for one sample")
    add_rose_dir(b)
    b.add_argument("--bed", required=True)
    b.add_argument("--rank-bam", required=True)
    b.add_argument("--control-bam", default="")
    b.add_argument("--outdir", required=True)
    b.add_argument("--annot", default="")
    b.add_argument("--gff-workers", type=int, default=4)
    b.add_argument("--slow-gff", action="store_true", help="use original mapper")
    b.add_argument("--resume", action="store_true")
    add_stitch_args(b)
    b.set_defaults(func=cmd_bam)

    s = sub.add_parser("stitch", help="batch stitch with official ROSE code")
    add_rose_dir(s)
    s.add_argument("--manifest", required=True)
    s.add_argument("--marks", default="H3K27me3")
    s.add_argument("--annot", default="")
    s.add_argument("--outdir", required=True)
    add_stitch_args(s)
    s.set_defaults(func=cmd_stitch)

    c = sub.add_parser("score", help="batch score stitched BED with bigWig")
    c.add_argument("--stitched-dir", required=True)
    c.add_argument("--manifest", required=True)
    c.add_argument("--bigwig-manifest", required=True)
    c.add_argument("--marks", default="H3K27me3")
    c.add_argument("--outdir", required=True)
    c.add_argument("--no-plot", action="store_true")
    c.set_defaults(func=cmd_score)

    w = sub.add_parser("bigwig", help="stitch + score one sample from bigWig")
    add_rose_dir(w)
    w.add_argument("--cell", required=True)
    w.add_argument("--mark", required=True)
    w.add_argument("--bed", required=True)
    w.add_argument("--bigwig", required=True)
    w.add_argument("--outdir", required=True)
    w.add_argument("--annot", default="")
    w.add_argument("--no-plot", action="store_true")
    add_stitch_args(w)
    w.set_defaults(func=cmd_bigwig)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
