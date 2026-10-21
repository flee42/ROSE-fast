# ROSE-Fast — a faster ROSE with the same algorithm

> **This section describes the changes made in this fork.** The upstream README follows
> after the separator at the bottom.

This fork keeps **every algorithmic step of upstream ROSE** and replaces only what is
provably equivalent or has no upstream counterpart. Each route deviates from upstream in
**exactly one** step:

```
BAM    : bedToGFF → regionStitching → [fast_bam_to_gff] → mapCollection → ROSE_callSuper.R
bigWig : bedToGFF → regionStitching → [bigWig scorer + REGION_MAP] → ROSE_callSuper.R
          └──────────── upstream ───────────┘   └─ ours ─┘          └── upstream ──┘
```

## What changed

| # | Change | Where | Why it is safe |
|---|---|---|---|
| 1 | **`bin/fast_bam_to_gff.py`** — numpy difference-array BAM→GFF mapper | BAM route, step 3 | Replaces the per-base Python dict in `ROSE_bamToGFF.py`. Entry-level **Spearman 1.0000**, median relative difference 0.04%, end-to-end **super-set Jaccard 1.0000** (H1 392/392, K562 987/987). |
| 2 | **Skip the two dead `bamToGFF` runs** | both | Upstream maps the *original* (unstitched) GFF in addition to the stitched one, but `mapCollection` only ever reads `<stitched>_<bam>_MAPPED.gff`. Two of the four runs are pure waste → **2× less mapping work**, identical output. |
| 3 | **`bin/rose_stitch_official.py`** — stitching without any BAM | bigWig route | Calls upstream `ROSE_main.regionStitching` directly; additionally decompresses gzipped peak files, which `ROSE_utils.parseTable` cannot read (`UnicodeDecodeError: 0x8b`). |
| 4 | **`bin/rose_call_super_R_native.py`** — bigWig scoring → 7-column `REGION_MAP` → **upstream `ROSE_callSuper.R`** with `wceName='NONE'` | bigWig route | No upstream equivalent exists for bigWig input, so only the *scorer* is new; the **cutoff and super-calling are 100% upstream**. This is the **only** bigWig scorer in the repo — there is deliberately no Python re-implementation of the cutoff. |
| 5 | **`bin/rose_run_one.py`** — sequential orchestrator | BAM route | `ROSE_main.py` launches `bamToGFF` with `os.system(cmd + " &")` and then polls with `time.sleep(300)`; it dies silently if the child is killed. This version calls ROSE's own functions in-process. |
| 6 | **`bin/rose_fast.py`** — `bam \| stitch \| score \| bigwig` dispatcher | — | Thin wrapper; contains no algorithm. |
| 7 | **`bin/rose_super_io.py`**, **`tools/fix_signal_isuper.py`** | — | Helpers: format-agnostic super-table reader; repair `isSuper` mis-alignment caused by `writeSuperEnhancer_table()` re-sorting rows by `stitchedPeakRank`. |
| 8 | **`tools/`** | — | Validation only: `verify_fast_gff.py`, `compare_super_sets.py`, `compare_super_bp_jaccard.py`, `compare_bigwig_vs_bam_rose.py`, `rose_region_stitching.py` (reference port), `rose_diag_bigwig.py`. |

**Untouched upstream files:** `bin/ROSE_main.py`, `bin/ROSE_bamToGFF.py`,
`bin/ROSE_callSuper.R`, `bin/ROSE_geneMapper.py`, `lib/ROSE_utils.py`,
`annotation/hg38_refseq.ucsc` — `git diff --name-only` is empty.

## Benchmark

### BAM route — fast vs original (K562, H3K27me3)

| metric | original | fast | ratio |
|---|---|---|---|
| `bamToGFF` runs | 4 | 2 | 2× fewer |
| wall time (single mapping) | ~50 min | **2 min 22 s** | **13–21×** |
| stitched regions | 19,741 | 19,740 | Jaccard **0.9993** |
| super regions | 987 | 987 | Jaccard **1.0000** (392/392 H1) |

### three-way check: original ROSE vs ROSE-fast BAM vs ROSE-fast bigWig

Same sample, same peaks, same control (`K562` / `H3K27ac`, `ENCSR000AKP`):

| pair | stitched Jaccard | super Jaccard | cutoff ratio |
|---|---|---|---|
| original vs fast-BAM | _pending_ | _pending_ | _pending_ |
| original vs fast-bigWig | _pending_ | _pending_ | _pending_ |
| fast-BAM vs fast-bigWig | _pending_ | _pending_ | _pending_ |

> Run `bash tools/compare_bigwig_vs_bam_rose.py` / `tools/compare_super_sets.py` to
> reproduce; the table above is filled in once the full run completes.

### bigWig vs BAM, population level (10 cells × 3 marks)

| layer | result |
|---|---|
| stitching | Jaccard 0.9993 (identical code) |
| cutoff | Jaccard 0.9975 (H1) / 0.9990 (K562) |
| super set (exact) | **Jaccard 1.0000** H1 and K562 |
| super set (base-pair) | 0.9149 H1 / 0.7011 K562 — precision 98.3% / 93.0%, recall 93.0% / 74.0% |
| `frac_super` | +0.2 / +0.4 / +0.1 pp |
| downstream enrichment ratio (bigWig/BAM) | **0.86–1.07**, median **0.96**; all 4 boundary-mark orderings identical |

The bigWig route is therefore a **high-precision subset** of the BAM route at the
individual-region level, and equivalent at the population level.

## Usage

### Requirements

Python 3 with `numpy`, `pandas`, `pysam` (BAM route) and `pyBigWig` (bigWig route);
`Rscript` with upstream `ROSE_callSuper.R`; `samtools`.

```bash
ROSEPATH=/path/to/ROSE-fast
export PYTHONPATH="$ROSEPATH/lib:$PYTHONPATH"
export PATH="$ROSEPATH/bin:$PATH"
```

### 1. BAM route — one sample

```bash
python bin/rose_fast.py bam \
    --bed        peaks.bed \
    --rank-bam   rank.bam \
    --control-bam control.bam \
    --outdir     out/ \
    --stitch 12500 --tss 3000 --gff-workers 4
```

Omit `--control-bam` for a rank-only run. `--resume` reuses existing `*_MAPPED.gff`.

### 2. bigWig route — one sample

```bash
python bin/rose_fast.py bigwig \
    --cell K562 --mark H3K27me3 \
    --bed    peaks.bed \
    --bigwig fold_change.bw \
    --outdir out_bw/ \
    --stitch 12500 --tss 3000
```

The bigWig must be **fold-change over control** (ENCODE's standard `*.bigWig`).
The signal is `Σ over constituent peaks (mean fold-change × peak length)`.

### 3. Batch: stitch many samples, then score

```bash
# stitch (upstream regionStitching, no BAM needed)
python bin/rose_stitch_official.py \
    --manifest manifest.tsv --marks H3K27me3 \
    --annot annotation/hg38_refseq.ucsc \
    --tss 3000 --stitch 12500 --outdir stitch_out/

# score + upstream ROSE_callSuper.R
python bin/rose_call_super_R_native.py \
    --stitched-dir stitch_out/stitched \
    --manifest manifest.tsv --bigwig-manifest bigwig_manifest.tsv \
    --marks H3K27me3 --outdir super_out/
```

`--reuse-signal --signal-from <dir>` reuses an existing `*_signal.tsv` `signal` column
instead of re-scanning the bigWig (saves hours on large cohorts).

### 4. Validate

```bash
python tools/verify_fast_gff.py       --fast ... --slow ...        # mapper equivalence
python tools/compare_super_sets.py    --a ... --b ...              # super-set Jaccard
python tools/compare_bigwig_vs_bam_rose.py --bam ... --bigwig ...  # route comparison
python tools/fix_signal_isuper.py <super_dir>                       # repair isSuper after R re-sorts
Rscript tools/rose_cutoff_R.R <super_dir>                           # cutoff cross-check (QA only)
```

### Outputs

| file | content |
|---|---|
| `<name>_AllStitched.table.txt` | all stitched regions, ranked, with `isSuper` |
| `<name>_SuperStitched.table.txt` | super regions only |
| `<name>_Stitched_withSuper.bed` | BED track |
| `<name>_Plot_points.png` | hockey-stick plot |
| `<name>_signal.tsv` | `chr,start,end,signal,rank,isSuper` (convenience; also emitted by the bigWig route) |

## Pitfalls that must be replicated exactly

Getting any of these wrong changes the result:

1. Read coverage is **inclusive of `end`** (`range(start, end+1)`); the half-open pysam convention gives a uniform 0.996 bias.
2. Region sum over `start < x < end-1`, divisor `len - 1`.
3. Search window = region ± `extension`; reads extended by `extension` (`+`: `end+ext`, `-`: `start-ext`).
4. Positions with coverage ≤ `floor` (`-f 1`) are dropped.
5. `MMR = round(mapped_reads / 1e6, 4)`; density = `binDen / MMR`.
6. `regionStitching` has **three** steps: drop constituents **fully contained** in TSS ± tssWindow; stitch; **revert any stitched locus overlapping >2 distinct gene TSSs (TSS ± 50 bp)** to its constituents. Missing step 3 under-produces regions (K562 15,405 vs 19,740).
7. Use ROSE's `annotation/hg38_refseq.ucsc` (58,350 transcripts / 26,514 symbols), not GENCODE (over-reverts: K562 27,994 vs 19,741).
8. `ROSE_utils.parseTable` cannot read gzip — decompress first.
9. `ROSE_callSuper.R` re-sorts `*_AllStitched.table.txt` by `stitchedPeakRank`, so `isSuper` must be mapped back **by region id or coordinates**, never by row position.

---

<!-- ===================== upstream README below ===================== -->


# ROSE: Rank Ordering of Super-Enhancers

ROSE (Rank Ordering of Super-Enhancers) is a computational pipeline for the identification and analysis of super-enhancers from ChIP-seq data. This repository provides a Python 3-compatible, location-independent implementation, with additional support for Docker-based deployment.

## Overview

ROSE enables researchers to:
- Identify super-enhancers from ChIP-seq data
- Rank enhancer regions by signal density
- Map stitched enhancers to genes

## Key Changes (from [source](https://bitbucket.org/young_computation/rose/src/master/))

- Compatible with Python 3
- Executable independent of installation directory
- Available as a Docker image: `ghcr.io/stjude/abralab/rose:latest`

## Installation & Setup

Clone the repository and set up the environment variables:

```bash
ROSEPATH=/path/to/ROSE
export PYTHONPATH="$PYTHONPATH:$ROSEPATH/lib"
export PATH="$PATH:$ROSEPATH/bin"
```

## Usage

Run the main program with the following command:

```bash
ROSE_main.py -g [GENOME] -i [INPUT_REGION_GFF] -r [RANKBY_BAM_FILE] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]
```

## Dependencies

- Python 3
- [samtools](http://www.htslib.org/)
- R (version > 3.4)
- bedtools (version > 2)

## Input Requirements

All input files must reside in a single directory.

### Annotation File (`-g [GENOME]`)
- Must be in UCSC table track format ([UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables))
- Filename: `[GENOME]_refseq.ucsc` (e.g., `hg19_refseq.ucsc`)
- Place in the `annotation/` folder within the input directory

### BAM File (`-r [RANKBY_BAM_FILE]`)
- Chromosome IDs must start with `chr`
- Files must be coordinate sorted and indexed using SAMtools

### Peak File (Constituent Enhancers, `-i [INPUT_REGION_GFF]`)
- Format: GFF

	- Column 1: Chromosome (chr#)
	- Column 2: Unique ID for each constituent enhancer region
	- Column 3: Type
	- Column 4: Start position
	- Column 5: End position
	- Column 6: Score
	- Column 7: Strand (`+`, `-`, or `.`)
	- Column 8: Phase
	- Column 9: Unique ID for each constituent enhancer region
- Columns 3,6,8 are not used by ROSE; they can have the value `.`.
- **Note:** If values in columns 2 and 9 differ, column 2 will be used.

## Directory Structure

```
├── LICENSE.txt
├── README.md
├── bin/
│   ├── ROSE_bamToGFF.py      # Calculates density of .bam reads in .gff regions
│   ├── ROSE_callSuper.R      # Ranks regions by their densities, creates cutoff
│   ├── ROSE_geneMapper.py    # Assigns stitched enhancers to genes
│   └── ROSE_main.py          # Main program
└── lib/
	└── ROSE_utils.py         # Utility methods
```

## References

- Original source: [Bitbucket - young_computation/rose](https://bitbucket.org/young_computation/rose/src/master/)

For questions or support, please contact the repository maintainers.
