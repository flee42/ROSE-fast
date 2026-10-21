#!/usr/bin/env python3
"""rose_run_one.py -- run the whole ROSE pipeline for ONE (cell, mark), robustly.

Steps (all done here, no ROSE_main.py polling):
  1. BED  -> GFF                     (ROSE_utils.bedToGFF)
  2. regionStitching  -> stitchedGFF  (ROSE_main.regionStitching)
  3. ROSE_bamToGFF, 4 jobs in parallel (stitched+original) x (rank+control)
  4. mapCollection -> *_REGION_MAP.txt
  5. Rscript ROSE_callSuper.R -> AllStitched / SuperStitched

usage:
  python3 rose_run_one.py \
      --rose-dir ROSE \
      --bed rose_bam/peaks/K562_H3K27me3_peaks.bed \
      --rank-bam rose_bam/merged/K562/H3K27me3/rank.bam \
      --control-bam rose_bam/merged/K562/control_ENCSR000AKY.bam \
      --outdir rose_out/K562_H3K27me3 \
      --stitch 12500 --tss 3000
"""
import argparse, os, subprocess, sys, time


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--rose-dir", default="ROSE")
    p.add_argument("--bed", required=True)
    p.add_argument("--rank-bam", required=True)
    p.add_argument("--control-bam", default="")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stitch", type=int, default=12500)
    p.add_argument("--tss", type=int, default=0)
    p.add_argument("--annot", default="", help="default: <rose-dir>/annotation/hg38_refseq.ucsc")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--resume", action="store_true", help="reuse existing mappedGFF, skip mapping")
    p.add_argument("--fast-gff", action="store_true", default=True,
                   help="use fast_bam_to_gff.py instead of ROSE_bamToGFF.py (default)")
    p.add_argument("--slow-gff", dest="fast_gff", action="store_false")
    p.add_argument("--gff-workers", type=int, default=4,
                   help="fast_bam_to_gff workers; 2-4 is safe on memory-limited nodes")
    return p.parse_args()


def main():
    a = parse_args()
    rose = os.path.abspath(a.rose_dir)
    annot = a.annot or os.path.join(rose, "annotation", "hg38_refseq.ucsc")
    sys.path.insert(0, os.path.join(rose, "lib"))
    sys.path.insert(0, os.path.join(rose, "bin"))
    import ROSE_utils
    import ROSE_main

    out = os.path.abspath(a.outdir) + "/"
    gff_dir, mapped = out + "gff/", out + "mappedGFF/"
    os.makedirs(gff_dir, exist_ok=True); os.makedirs(mapped, exist_ok=True)
    base = os.path.basename(a.bed)
    for _ext in (".narrowPeak.gz", ".broadPeak.gz", ".bed.gz",
                 ".narrowPeak", ".broadPeak", ".bed"):
        if base.endswith(_ext):
            base = base[: -len(_ext)]
            break
    suffix = "_%gKB_STITCHED" % (a.stitch / 1000.0) + ("_TSS_DISTAL" if a.tss else "")
    stitchedName = base + suffix
    inputGFF = gff_dir + base + ".gff"
    stitchedGFF = gff_dir + stitchedName + ".gff"

    print("=" * 70); print("ROSE:", base); print("out :", out)

    # 1) BED -> GFF
    if not os.path.exists(inputGFF):
        ROSE_utils.bedToGFF(a.bed, inputGFF)
    print("[1] input GFF :", inputGFF, os.path.getsize(inputGFF), "bytes")

    # 2) stitching
    refCol = ROSE_utils.gffToLocusCollection(inputGFF)
    stitched, dbg = ROSE_main.regionStitching(inputGFF, a.stitch, a.tss, annot,
                                              removeTSS=bool(a.tss))
    if not os.path.exists(stitchedGFF):
        ROSE_utils.unParseTable(ROSE_utils.locusCollectionToGFF(stitched), stitchedGFF, "\t")
    print("[2] stitched  :", stitchedGFF, len(stitched.getLoci()), "loci")

    bams = [a.rank_bam] + ([a.control_bam] if a.control_bam else [])
    # ROSE's mapCollection only ever reads <stitchedName>_<bam>_MAPPED.gff, and
    # ROSE_geneMapper never touches the _MAPPED files at all -> the mapping of the
    # ORIGINAL gff is dead work (2 of the 4 runs).  Skip it.
    jobs = []
    for bam in bams:
        bn = os.path.basename(bam)
        jobs.append((bam, stitchedGFF, mapped + stitchedName + "_" + bn + "_MAPPED.gff"))

    # 3) bamToGFF
    todo = [j for j in jobs if not (a.resume and os.path.exists(j[2]))]
    if todo:
        print("[3] running %d ROSE_bamToGFF jobs ..." % len(todo))
        procs = []
        for bam, gff, dst in todo:
            if a.fast_gff:
                tool = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fast_bam_to_gff.py")
                cmd = [sys.executable, tool, "-f", "1", "-e", "200", "-r", "-m", "1",
                       "-w", str(a.gff_workers), "-b", bam, "-i", gff, "-o", dst]
            else:
                cmd = [sys.executable, os.path.join(rose, "bin", "ROSE_bamToGFF.py"),
                       "-f", "1", "-e", "200", "-r", "-m", "1", "-b", bam, "-i", gff, "-o", dst]
            log = open(dst + ".log", "w")
            procs.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT), dst, log))
        while procs:
            time.sleep(20)
            still = []
            for p, dst, log in procs:
                if p.poll() is None:
                    still.append((p, dst, log))
                else:
                    log.close()
                    ok = os.path.exists(dst)
                    print("    done rc=%s %s  %s bytes" % (p.returncode, os.path.basename(dst),
                          os.path.getsize(dst) if ok else "MISSING"))
            procs = still
    else:
        print("[3] all mappedGFF present -- skipped")

    # 4) mapCollection
    regionMap = out + stitchedName + "_REGION_MAP.txt"
    print("[4] mapCollection ->", regionMap)
    ROSE_main.mapCollection(stitched, refCol, bams, mapped, regionMap, refName=stitchedName)

    # 5) ROSE_callSuper.R
    ctl = os.path.basename(a.control_bam) if a.control_bam else "NONE"
    cmd = ["Rscript", os.path.join(rose, "bin", "ROSE_callSuper.R"), out, regionMap, base, ctl]
    print("[5]", " ".join(cmd))
    rc = subprocess.call(cmd)
    print("[5] exit", rc)
    print("\n--- outputs ---")
    for f in sorted(os.listdir(out)):
        p = os.path.join(out, f)
        print("   %-70s %s" % (f, os.path.getsize(p) if os.path.isfile(p) else "<dir>"))


if __name__ == "__main__":
    main()
