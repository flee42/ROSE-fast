#!/usr/bin/env python3
"""rose_super_io.py -- single safe reader for ROSE super-region tables.

Handles all three shapes we produce:
  a) ROSE native       <x>_SuperStitched.table.txt   (has '#' comment header lines)
  b) ROSE native       <x>_AllStitched.table.txt     (same, filtered by isSuper)
  c) bigWig proxy      <x>_signal.tsv                (plain header, has isSuper)

load_super(path) accepts a directory or a file and always returns
    (DataFrame[chr, start, end], cell_name)
containing ONLY the super regions.
"""
import glob, os
import pandas as pd

TAB = chr(9)


def read_rose_table(path):
    """ROSE writes several '#' description lines before the real header."""
    raw = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    hi = 0
    for i, l in enumerate(raw[:30]):
        u = l.upper()
        if ("CHROM" in u) or ("REGION_ID" in u):
            hi = i
            break
    sep = TAB if TAB in raw[hi] else r"\s+"
    t = pd.read_csv(path, sep=sep, skiprows=hi, engine="python")
    t.columns = [str(c).strip().strip('"') for c in t.columns]
    ren = {}
    for c in t.columns:
        cu = c.upper()
        if cu == "CHROM":
            ren[c] = "chr"
        elif cu == "START":
            ren[c] = "start"
        elif cu in ("STOP", "END"):
            ren[c] = "end"
        elif cu == "ISSUPER":
            ren[c] = "isSuper"
    t = t.rename(columns=ren)
    miss = {"chr", "start", "end"} - set(t.columns)
    if miss:
        raise ValueError("cannot find %s in %s (cols=%s)" % (miss, path, list(t.columns)[:12]))
    return t


def _cell_of(path):
    b = os.path.basename(path)
    for suf in ("_SuperStitched.table.txt", "_AllStitched.table.txt",
                "_Stitched_withSuper.bed", "_signal.tsv"):
        if b.endswith(suf):
            b = b[: -len(suf)]
            break
    for tok in ("_H3K27me3", "_H3K27ac", "_H3K4me1", "_H3K4me3", "_H3K9me3",
                "_H3K36me3", "_H2AK119ub"):
        if tok in b:
            return b.split(tok)[0]
    return b


def load_super(path):
    """dir or file -> (DataFrame[chr,start,end] of SUPER regions, cell name)"""
    cand = []
    if os.path.isdir(path):
        cand = (sorted(glob.glob(os.path.join(path, "*_SuperStitched.table.txt")))
                or sorted(glob.glob(os.path.join(path, "*_AllStitched.table.txt")))
                or sorted(glob.glob(os.path.join(path, "*_signal.tsv"))))
    elif os.path.isfile(path):
        if not path.endswith(("_SuperStitched.table.txt", "_AllStitched.table.txt", "_signal.tsv")):
            return None, None          # e.g. *_Plot_points.png from glob
        cand = [path]
    if not cand:
        return None, None
    f = cand[0]
    if f.endswith(".tsv"):
        t = pd.read_csv(f, sep=TAB)
        if "isSuper" in t.columns:
            t = t[t["isSuper"].astype(float) > 0]
        if not {"chr", "start", "end"}.issubset(t.columns):
            return None, None
        return t[["chr", "start", "end"]].astype({"start": int, "end": int}), _cell_of(f)
    t = read_rose_table(f)
    if "isSuper" in t.columns:
        t = t[t["isSuper"].astype(float) > 0]
    if t.empty:
        return None, None
    return t[["chr", "start", "end"]].astype({"start": int, "end": int}), _cell_of(f)


SUPER_SUFFIX_RANK = {"_SuperStitched.table.txt": 0,
                     "_AllStitched.table.txt": 1,
                     "_signal.tsv": 2}


def pick_one_per_cell(paths):
    """Several files can describe the SAME super set (SuperStitched / AllStitched
    filtered / signal.tsv). Keep exactly one per (cell, mark), preferring the
    ROSE-native SuperStitched table."""
    best = {}
    for p in paths:
        c = _cell_of(p)
        r = next((v for k, v in SUPER_SUFFIX_RANK.items() if p.endswith(k)), 9)
        if c not in best or r < best[c][0]:
            best[c] = (r, p)
    return sorted(v[1] for v in best.values())
