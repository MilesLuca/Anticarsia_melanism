# Repeat enrichment in matched genomic windows

**Paper:** Figure S12; percentile labels in Figure 3A  

Compare each GWA interval’s repeat-covered fraction with the genome-wide distribution of full-sized windows of the same length, advanced in 50-kb steps.

## Inputs

- The three representative genome FASTAs and their .fai indexes.
- The three filtered repeat GFFs from analysis 05.

## Outputs

- EarlGrey/TE_class_window_histograms/03_tables/te_class_observed_vs_windows_summary.tsv
- EarlGrey/TE_class_window_histograms/03_tables/te_class_window_values_long.tsv
- EarlGrey/TE_class_window_histograms/04_plots/*_te_class_histograms_full.{png,pdf}

## Software

Python 3, numpy, pandas, matplotlib and SAMtools (for absent FASTA indexes).

## Execution and interpretation

1. The percentile is 100 × fraction of null-window values less than or equal to the observed value. Union coverage prevents repeated counting of overlapping annotations within a class; Total is a separate union.
2. Final full-distribution output: A 87.3425th percentile, B 88.1995th, C 83.0855th. These correspond to top 12.7%, 11.8% and 16.9%.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Calculate enrichment and plot full distributions

**Source:** `EarlGrey/TE_class_window_histograms/05_scripts/run_te_class_window_histograms.py`  
**Save as:** `run_te_class_window_histograms.py`  
**Archive treatment:** Omit only source lines 381, 391, 392 (trimmed-figure call and its two output messages); calculations unchanged. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples.

```python
#!/usr/bin/env python3
import os
import csv
import math
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = "/path/to/project/EarlGrey/TE_class_window_histograms"
STEP = 50000

GENOMES = [
    {
        "genome_id": "dark_ref",
        "label": "Dark Reference",
        "fasta": "/path/to/project/Agem_genome/NCBI/GCF_050436995.1_ilAntGemm2_primary_genomic.fna",
        "gff": "/path/to/project/EarlGrey/Final_Dark_Genome_Run/o/ADGO_EarlGrey/ADGO_summaryFiles/ADGO.filteredRepeats.gff",
        "gwa_seqid": "NC_134752.1",
        "gwa_start": 2012000,
        "gwa_end": 2613000,
    },
    {
        "genome_id": "plain_ref",
        "label": "Plain Reference",
        "fasta": "/path/to/project/Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold/agem_light_06.ragtag.final.fasta",
        "gff": "/path/to/project/EarlGrey/Final_Light_Genome_Run/o/ALGO_EarlGrey/ALGO_summaryFiles/ALGO.filteredRepeats.gff",
        "gwa_seqid": "chr8",
        "gwa_start": 1942695,
        "gwa_end": 2491580,
    },
    {
        "genome_id": "dark01_h2",
        "label": "Dark01 h2",
        "fasta": "/path/to/project/Pacbio_Reseq_Hilo/Agem_Dark_01/hifiasm_from_hfaf_t1/haplotype_fastas/Agem_Dark_01.bp.hap2.p_ctg.fa",
        "gff": "/path/to/project/EarlGrey/Final_Dark_01_h2_haplogenome_run/o/AD1H2_EarlGrey/AD1H2_summaryFiles/AD1H2.filteredRepeats.gff",
        "gwa_seqid": "h2tg000012l",
        "gwa_start": 1913107,
        "gwa_end": 2487925,
    },
]

CLASS_ORDER = [
    "DNA",
    "LINE",
    "LTR",
    "Rolling_Circle",
    "SINE",
    "Simple_Repeat",
    "Unclassified",
    "Total",
]

CLASS_LABELS = {
    "DNA": "DNA",
    "LINE": "LINE",
    "LTR": "LTR",
    "Rolling_Circle": "Rolling Circle",
    "SINE": "SINE",
    "Simple_Repeat": "Simple Repeat",
    "Unclassified": "Unclassified",
    "Total": "Total",
}

COLORS = {
    "DNA": "#ff0000",
    "LINE": "#24b7df",
    "LTR": "#00a000",
    "Rolling_Circle": "#f4a300",
    "SINE": "#c0008f",
    "Simple_Repeat": "#12cbe8",
    "Unclassified": "#b3b3b3",
    "Total": "#3a3a3a",
}

os.makedirs(os.path.join(ROOT, "02_class_beds"), exist_ok=True)
os.makedirs(os.path.join(ROOT, "03_tables"), exist_ok=True)
os.makedirs(os.path.join(ROOT, "04_plots"), exist_ok=True)

def run(cmd):
    subprocess.run(cmd, check=True)

def ensure_fai(fasta):
    if not os.path.exists(fasta + ".fai"):
        run(["samtools", "faidx", fasta])

def read_fai_lengths(fasta):
    ensure_fai(fasta)
    lengths = {}
    with open(fasta + ".fai") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            if e > cur_e:
                cur_e = e
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged

def write_bed(path, interval_dict):
    with open(path, "w") as out:
        for seqid in sorted(interval_dict.keys()):
            for s, e in interval_dict[seqid]:
                out.write(f"{seqid}\t{s}\t{e}\n")

def map_class(feature):
    f = feature.strip()
    if f == "DNA" or f.startswith("DNA/"):
        return "DNA"
    if f.startswith("LINE"):
        return "LINE"
    if f.startswith("LTR"):
        return "LTR"
    if f == "RC" or f.startswith("RC/") or "Helitron" in f:
        return "Rolling_Circle"
    if f.startswith("SINE"):
        return "SINE"
    if f in {
        "Simple_repeat", "Satellite", "Low_complexity",
        "RNA", "rRNA", "tRNA", "snRNA", "scRNA", "srpRNA",
        "Satellite/centr", "Satellite/telo"
    }:
        return "Simple_Repeat"
    if f == "Unknown" or f.startswith("Unknown"):
        return "Unclassified"
    if f == "Penelope" or f.startswith("Penelope/"):
        return None
    return None

def parse_gff_to_class_intervals(gff):
    raw = {c: defaultdict(list) for c in CLASS_ORDER}
    with open(gff) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            seqid = parts[0]
            feature = parts[2]
            try:
                start0 = int(parts[3]) - 1
                end = int(parts[4])
            except ValueError:
                continue
            if start0 < 0:
                start0 = 0

            raw["Total"][seqid].append((start0, end))
            cls = map_class(feature)
            if cls is not None:
                raw[cls][seqid].append((start0, end))

    merged = {}
    for cls in CLASS_ORDER:
        merged[cls] = {}
        for seqid, ints in raw[cls].items():
            merged[cls][seqid] = merge_intervals(ints)
    return merged

def make_windows(lengths, window_size, step):
    windows = {}
    for seqid, L in lengths.items():
        seq_windows = []
        if L >= window_size:
            s = 0
            while s + window_size <= L:
                seq_windows.append((s, s + window_size))
                s += step
        windows[seqid] = seq_windows
    return windows

def coverage_percent_for_windows(windows_by_seqid, intervals_by_seqid, window_size):
    out = []
    for seqid, windows in windows_by_seqid.items():
        intervals = intervals_by_seqid.get(seqid, [])
        j = 0
        n = len(intervals)
        for ws, we in windows:
            while j < n and intervals[j][1] <= ws:
                j += 1
            k = j
            covered = 0
            while k < n and intervals[k][0] < we:
                s = max(ws, intervals[k][0])
                e = min(we, intervals[k][1])
                if e > s:
                    covered += (e - s)
                k += 1
            out.append(100.0 * covered / window_size)
    return np.array(out, dtype=float)

def coverage_percent_for_single_interval(seqid, start0, end, intervals_by_seqid):
    window_size = end - start0
    intervals = intervals_by_seqid.get(seqid, [])
    covered = 0
    for s0, e0 in intervals:
        if e0 <= start0:
            continue
        if s0 >= end:
            break
        s = max(start0, s0)
        e = min(end, e0)
        if e > s:
            covered += (e - s)
    return 100.0 * covered / window_size

def percentile_of_observed(values, obs):
    if len(values) == 0:
        return float("nan")
    return 100.0 * np.mean(values <= obs)

def plot_one_genome(summary_df, windows_df, genome_id, label, trim=False):
    ssub = summary_df[summary_df["genome_id"] == genome_id].copy()
    wsub = windows_df[windows_df["genome_id"] == genome_id].copy()

    fig, axes = plt.subplots(len(CLASS_ORDER), 1, figsize=(8.7, 12.8))
    fig.subplots_adjust(hspace=0.35, left=0.12, right=0.88, top=0.93, bottom=0.08)

    mode = "99.5th-percentile trimmed facets" if trim else "full-null histograms"
    fig.suptitle(f"{label}: Full GWA interval ({mode})", fontsize=14)

    for ax, cls in zip(axes, CLASS_ORDER):
        d = wsub[wsub["class_id"] == cls]["window_pct"].to_numpy(dtype=float)
        row = ssub[ssub["class_id"] == cls].iloc[0]
        obs = float(row["observed_pct"])
        pctile = float(row["percentile"])
        trim_x = float(row["trim_99_5_pct"])

        ax.set_facecolor("#EBEBEB")
        ax.grid(True, color="white", linewidth=1.0)
        ax.set_axisbelow(True)

        if trim:
            plot_vals = d[d <= trim_x]
            xmax = max(trim_x, obs) * 1.03 if max(trim_x, obs) > 0 else 1.0
        else:
            plot_vals = d
            xmax = max(float(np.max(d)) if len(d) else 0.0, obs) * 1.03
            if xmax <= 0:
                xmax = 1.0

        if len(plot_vals) == 0:
            plot_vals = np.array([0.0])

        ax.hist(plot_vals, bins=30, color=COLORS[cls], edgecolor="black", linewidth=0.8)
        ax.axvline(obs, color="red", linewidth=1.6)
        ax.set_xlim(0, xmax)

        ann = f"obs = {obs:.2f}%\npercentile = {pctile:.1f}%"
        ax.text(
            0.98, 0.80, ann,
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="square,pad=0.25")
        )

        ax.text(
            1.01, 0.5, CLASS_LABELS[cls],
            transform=ax.transAxes,
            rotation=-90,
            va="center", ha="left",
            fontsize=10,
            bbox=dict(facecolor="#D9D9D9", edgecolor="black", boxstyle="square,pad=0.25")
        )

        if ax is axes[-1]:
            ax.set_xlabel("Repeat-covered % of interval-sized window", fontsize=11)
        else:
            ax.set_xlabel("")

    fig.text(0.02, 0.5, "Number of genome windows", rotation=90, va="center", fontsize=11)

    suffix = "trim99_5" if trim else "full"
    png = os.path.join(ROOT, "04_plots", f"{genome_id}_te_class_histograms_{suffix}.png")
    pdf = os.path.join(ROOT, "04_plots", f"{genome_id}_te_class_histograms_{suffix}.pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

summary_rows = []
window_rows = []

meta_rows = []

for g in GENOMES:
    genome_id = g["genome_id"]
    label = g["label"]
    fasta = g["fasta"]
    gff = g["gff"]
    gwa_seqid = g["gwa_seqid"]
    gwa_start = int(g["gwa_start"])
    gwa_end = int(g["gwa_end"])
    gwa_size = gwa_end - gwa_start + 1
    gwa_start0 = gwa_start - 1

    outdir = os.path.join(ROOT, "02_class_beds", genome_id)
    os.makedirs(outdir, exist_ok=True)

    lengths = read_fai_lengths(fasta)
    if gwa_seqid not in lengths:
        raise ValueError(f"{gwa_seqid} not found in {fasta}")

    windows = make_windows(lengths, gwa_size, STEP)
    n_windows_total = sum(len(v) for v in windows.values())

    merged_intervals = parse_gff_to_class_intervals(gff)

    for cls in CLASS_ORDER:
        write_bed(os.path.join(outdir, f"{cls}.bed"), merged_intervals[cls])

    meta_rows.append({
        "genome_id": genome_id,
        "label": label,
        "gwa_seqid": gwa_seqid,
        "gwa_start": gwa_start,
        "gwa_end": gwa_end,
        "gwa_size": gwa_size,
        "step_size": STEP,
        "n_windows_total": n_windows_total,
    })

    for cls in CLASS_ORDER:
        values = coverage_percent_for_windows(windows, merged_intervals[cls], gwa_size)
        obs = coverage_percent_for_single_interval(gwa_seqid, gwa_start0, gwa_end, merged_intervals[cls])
        pctile = percentile_of_observed(values, obs)
        trim_995 = float(np.quantile(values, 0.995)) if len(values) else 0.0
        maxv = float(np.max(values)) if len(values) else 0.0

        for v in values:
            window_rows.append({
                "genome_id": genome_id,
                "label": label,
                "class_id": cls,
                "class_label": CLASS_LABELS[cls],
                "window_pct": round(float(v), 6),
            })

        summary_rows.append({
            "genome_id": genome_id,
            "label": label,
            "gwa_seqid": gwa_seqid,
            "gwa_start": gwa_start,
            "gwa_end": gwa_end,
            "gwa_size": gwa_size,
            "step_size": STEP,
            "n_windows": int(len(values)),
            "class_id": cls,
            "class_label": CLASS_LABELS[cls],
            "observed_pct": round(float(obs), 6),
            "percentile": round(float(pctile), 4),
            "trim_99_5_pct": round(trim_995, 6),
            "max_window_pct": round(maxv, 6),
        })

meta_df = pd.DataFrame(meta_rows)
summary_df = pd.DataFrame(summary_rows)
windows_df = pd.DataFrame(window_rows)

meta_df.to_csv(os.path.join(ROOT, "03_tables", "genome_metadata.tsv"), sep="\t", index=False)
summary_df.to_csv(os.path.join(ROOT, "03_tables", "te_class_observed_vs_windows_summary.tsv"), sep="\t", index=False)
windows_df.to_csv(os.path.join(ROOT, "03_tables", "te_class_window_values_long.tsv"), sep="\t", index=False)

for g in GENOMES:
    plot_one_genome(summary_df, windows_df, g["genome_id"], g["label"], trim=False)

print("WROTE:", os.path.join(ROOT, "03_tables", "genome_metadata.tsv"))
print("WROTE:", os.path.join(ROOT, "03_tables", "te_class_observed_vs_windows_summary.tsv"))
print("WROTE:", os.path.join(ROOT, "03_tables", "te_class_window_values_long.tsv"))

for g in GENOMES:
    gid = g["genome_id"]
    print("WROTE:", os.path.join(ROOT, "04_plots", f"{gid}_te_class_histograms_full.png"))
    print("WROTE:", os.path.join(ROOT, "04_plots", f"{gid}_te_class_histograms_full.pdf"))
```
