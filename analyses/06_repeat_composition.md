# Repeat composition and donut plots

**Paper:** Figure 3A  

Compute repeat-class base coverage in the whole Dark reference and each GWA interval, then display class proportions as donuts.

## Inputs

- The three filtered repeat GFFs from analysis 05.
- The Dark reference FASTA for genome-wide length.

## Outputs

- EarlGrey/TE_piecharts/01_tables/gwa_te_broadclass_bp.tsv
- EarlGrey/TE_piecharts/01_tables/dark_reference_genomewide_broadclass_bp.tsv
- EarlGrey/TE_piecharts/01_tables/gwa_plus_genomewide_broadclass_bp.tsv
- EarlGrey/TE_piecharts/05_donuts/gwa_te_donut_combined.png and .pdf

## Software

Python 3, pandas and matplotlib; standard-library csv and collections.

## Execution and interpretation

1. Run the two builders, then the table-combination excerpt, then the donut plotter. The excerpt comes from a mixed table-building/pie-plotting script; obsolete pie rendering is omitted.
2. The whole-genome denominator is 391,691,109 bp. The original plot order is whole genome, A/Dark reference, C/Plain reference, B/Dark01 h2.
3. A: 601,001 bp and 34.1347% repeats. B: 574,819 bp and 31.0927%. C: 548,886 bp and 31.0547%..
4. Donut broad classes include an Other category; S12 separates Simple Repeat.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Build interval composition tables

Uses union coverage within each repeat class and records input intervals and diagnostics.

**Source:** `EarlGrey/TE_piecharts/03_scripts/build_gwa_te_tables.py`  
**Save as:** `build_gwa_te_tables.py`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples.

```python
#!/usr/bin/env python3
import csv
import os
from collections import defaultdict, Counter

genomes = [
    {
        "genome": "Dark Reference",
        "seqid": "NC_134752.1",
        "gwa_start": 2012000,
        "gwa_end": 2613000,
        "gff": "/path/to/project/EarlGrey/Final_Dark_Genome_Run/o/ADGO_EarlGrey/ADGO_summaryFiles/ADGO.filteredRepeats.gff",
    },
    {
        "genome": "Plain Reference",
        "seqid": "chr8",
        "gwa_start": 1942695,
        "gwa_end": 2491580,
        "gff": "/path/to/project/EarlGrey/Final_Light_Genome_Run/o/ALGO_EarlGrey/ADGO_summaryFiles/../ALGO_summaryFiles/ALGO.filteredRepeats.gff",
    },
    {
        "genome": "Dark01 h2",
        "seqid": "h2tg000012l",
        "gwa_start": 1913107,
        "gwa_end": 2487925,
        "gff": "/path/to/project/EarlGrey/Final_Dark_01_h2_haplogenome_run/o/AD1H2_EarlGrey/AD1H2_summaryFiles/AD1H2.filteredRepeats.gff",
    },
]

for g in genomes:
    g["gff"] = os.path.normpath(g["gff"])

broad_order = [
    "DNA",
    "Rolling Circle",
    "Penelope",
    "LINE",
    "SINE",
    "LTR",
    "Other",
    "Unclassified",
    "Non-Repeat",
]

def map_feature(feature: str) -> str:
    f = feature.strip()

    if f == "DNA" or f.startswith("DNA/"):
        return "DNA"
    if f == "RC" or f.startswith("RC/") or "Helitron" in f:
        return "Rolling Circle"
    if f == "Penelope" or f.startswith("Penelope/"):
        return "Penelope"
    if f == "LINE" or f.startswith("LINE"):
        return "LINE"
    if f == "SINE" or f.startswith("SINE"):
        return "SINE"
    if f == "LTR" or f.startswith("LTR"):
        return "LTR"
    if f == "Unknown" or f.startswith("Unknown"):
        return "Unclassified"

    other_types = {
        "Simple_repeat",
        "Satellite",
        "Low_complexity",
        "RNA",
        "rRNA",
        "tRNA",
        "snRNA",
        "scRNA",
        "srpRNA",
        "Satellite/centr",
        "Satellite/telo",
    }
    if f in other_types:
        return "Other"

    return "Other"

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e + 1:
            if e > cur_e:
                cur_e = e
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged

def union_bp(intervals):
    merged = merge_intervals(intervals)
    return sum(e - s + 1 for s, e in merged)

outdir = "/path/to/project/EarlGrey/TE_piecharts/01_tables"
os.makedirs(outdir, exist_ok=True)

metadata_rows = []
size_rows = []
diag_rows = []
broad_rows = []
feature_rows = []
mapping_rows = []
clipped_rows = []

for g in genomes:
    genome = g["genome"]
    seqid = g["seqid"]
    gwa_start = g["gwa_start"]
    gwa_end = g["gwa_end"]
    gff = g["gff"]
    interval_bp = gwa_end - gwa_start + 1

    metadata_rows.append({
        "genome": genome,
        "seqid": seqid,
        "gwa_start": gwa_start,
        "gwa_end": gwa_end,
        "interval_bp": interval_bp,
        "gff": gff,
    })

    size_rows.append({
        "genome": genome,
        "seqid": seqid,
        "gwa_start": gwa_start,
        "gwa_end": gwa_end,
        "interval_bp": interval_bp,
    })

    class_intervals = defaultdict(list)
    feature_intervals = defaultdict(list)
    all_intervals = []
    feature_counts = Counter()
    mapping_counts = Counter()

    n_records_total = 0
    n_records_overlapping_interval = 0

    with open(gff) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue

            rec_seqid = parts[0]
            feature = parts[2]

            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue

            n_records_total += 1

            if rec_seqid != seqid:
                continue
            if start > gwa_end or end < gwa_start:
                continue

            n_records_overlapping_interval += 1
            cs = max(start, gwa_start)
            ce = min(end, gwa_end)
            if cs > ce:
                continue

            broad = map_feature(feature)
            bp = ce - cs + 1

            class_intervals[broad].append((cs, ce))
            feature_intervals[(feature, broad)].append((cs, ce))
            all_intervals.append((cs, ce))

            feature_counts[(feature, broad)] += 1
            mapping_counts[(feature, broad)] += 1

            clipped_rows.append({
                "genome": genome,
                "seqid": seqid,
                "gwa_start": gwa_start,
                "gwa_end": gwa_end,
                "feature": feature,
                "broad_class": broad,
                "clip_start": cs,
                "clip_end": ce,
                "clip_bp": bp,
            })

    total_repeat_union_bp = union_bp(all_intervals)
    nonrepeat_bp = interval_bp - total_repeat_union_bp

    merged_class_bp_sum = 0
    raw_class_bp_sum = 0
    per_class_union = {}

    for broad in broad_order:
        if broad == "Non-Repeat":
            continue
        ints = class_intervals.get(broad, [])
        bp_union = union_bp(ints)
        bp_raw = sum(e - s + 1 for s, e in ints)
        per_class_union[broad] = bp_union
        merged_class_bp_sum += bp_union
        raw_class_bp_sum += bp_raw

    cross_class_overlap_bp = merged_class_bp_sum - total_repeat_union_bp

    diag_rows.append({
        "genome": genome,
        "seqid": seqid,
        "gwa_start": gwa_start,
        "gwa_end": gwa_end,
        "interval_bp": interval_bp,
        "gff_records_total": n_records_total,
        "gff_records_overlapping_interval": n_records_overlapping_interval,
        "raw_repeat_bp_sum_clipped": raw_class_bp_sum,
        "merged_class_bp_sum": merged_class_bp_sum,
        "total_repeat_union_bp": total_repeat_union_bp,
        "cross_class_overlap_bp": cross_class_overlap_bp,
        "nonrepeat_bp": nonrepeat_bp,
        "repeat_pct_of_interval": round(100 * total_repeat_union_bp / interval_bp, 4),
        "nonrepeat_pct_of_interval": round(100 * nonrepeat_bp / interval_bp, 4),
    })

    for broad in broad_order:
        if broad == "Non-Repeat":
            bp = nonrepeat_bp
        else:
            bp = per_class_union.get(broad, 0)
        pct = 100 * bp / interval_bp if interval_bp > 0 else 0.0
        broad_rows.append({
            "genome": genome,
            "seqid": seqid,
            "gwa_start": gwa_start,
            "gwa_end": gwa_end,
            "interval_bp": interval_bp,
            "broad_class": broad,
            "bp_union": bp,
            "pct_interval": round(pct, 4),
        })

    for (feature, broad), ints in sorted(feature_intervals.items()):
        bp_union = union_bp(ints)
        bp_raw = sum(e - s + 1 for s, e in ints)
        pct = 100 * bp_union / interval_bp if interval_bp > 0 else 0.0
        feature_rows.append({
            "genome": genome,
            "seqid": seqid,
            "gwa_start": gwa_start,
            "gwa_end": gwa_end,
            "interval_bp": interval_bp,
            "original_feature": feature,
            "broad_class": broad,
            "n_records_clipped": feature_counts[(feature, broad)],
            "bp_raw_clipped": bp_raw,
            "bp_union": bp_union,
            "pct_interval": round(pct, 4),
        })

    for (feature, broad), n in sorted(mapping_counts.items()):
        mapping_rows.append({
            "genome": genome,
            "original_feature": feature,
            "broad_class": broad,
            "n_records_clipped": n,
        })

def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

write_tsv(
    os.path.join(outdir, "gwa_interval_metadata.tsv"),
    metadata_rows,
    ["genome", "seqid", "gwa_start", "gwa_end", "interval_bp", "gff"],
)

write_tsv(
    os.path.join(outdir, "gwa_interval_sizes.tsv"),
    size_rows,
    ["genome", "seqid", "gwa_start", "gwa_end", "interval_bp"],
)

write_tsv(
    os.path.join(outdir, "gwa_te_interval_diagnostics.tsv"),
    diag_rows,
    [
        "genome", "seqid", "gwa_start", "gwa_end", "interval_bp",
        "gff_records_total", "gff_records_overlapping_interval",
        "raw_repeat_bp_sum_clipped", "merged_class_bp_sum",
        "total_repeat_union_bp", "cross_class_overlap_bp",
        "nonrepeat_bp", "repeat_pct_of_interval", "nonrepeat_pct_of_interval"
    ],
)

write_tsv(
    os.path.join(outdir, "gwa_te_broadclass_bp.tsv"),
    broad_rows,
    ["genome", "seqid", "gwa_start", "gwa_end", "interval_bp", "broad_class", "bp_union", "pct_interval"],
)

write_tsv(
    os.path.join(outdir, "gwa_te_featuretype_bp.tsv"),
    feature_rows,
    [
        "genome", "seqid", "gwa_start", "gwa_end", "interval_bp",
        "original_feature", "broad_class", "n_records_clipped",
        "bp_raw_clipped", "bp_union", "pct_interval"
    ],
)

write_tsv(
    os.path.join(outdir, "gwa_te_featuretype_mapping.tsv"),
    mapping_rows,
    ["genome", "original_feature", "broad_class", "n_records_clipped"],
)

write_tsv(
    os.path.join(outdir, "gwa_te_clipped_records.tsv"),
    clipped_rows,
    ["genome", "seqid", "gwa_start", "gwa_end", "feature", "broad_class", "clip_start", "clip_end", "clip_bp"],
)

print("WROTE corrected GWA tables:")
for fn in [
    "gwa_interval_metadata.tsv",
    "gwa_interval_sizes.tsv",
    "gwa_te_interval_diagnostics.tsv",
    "gwa_te_broadclass_bp.tsv",
    "gwa_te_featuretype_bp.tsv",
    "gwa_te_featuretype_mapping.tsv",
    "gwa_te_clipped_records.tsv",
]:
    print(os.path.join(outdir, fn))
```

<a id="step-02"></a>

### 2. Build whole-genome composition table

Uses the full Dark reference, including sex chromosomes and unplaced scaffolds.

**Source:** `EarlGrey/TE_piecharts/03_scripts/build_dark_reference_genomewide_table.py`  
**Save as:** `build_dark_reference_genomewide_table.py`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples.

```python
#!/usr/bin/env python3
import csv
import os
from collections import defaultdict

dark_fasta = "/path/to/project/Agem_genome/NCBI/GCF_050436995.1_ilAntGemm2_primary_genomic.fna"
dark_gff = "/path/to/project/EarlGrey/Final_Dark_Genome_Run/o/ADGO_EarlGrey/ADGO_summaryFiles/ADGO.filteredRepeats.gff"

outdir = "/path/to/project/EarlGrey/TE_piecharts/01_tables"
os.makedirs(outdir, exist_ok=True)

broad_order = [
    "DNA",
    "Rolling Circle",
    "Penelope",
    "LINE",
    "SINE",
    "LTR",
    "Other",
    "Unclassified",
    "Non-Repeat",
]

def map_feature(feature: str) -> str:
    f = feature.strip()

    if f == "DNA" or f.startswith("DNA/"):
        return "DNA"
    if f == "RC" or f.startswith("RC/") or "Helitron" in f:
        return "Rolling Circle"
    if f == "Penelope" or f.startswith("Penelope/"):
        return "Penelope"
    if f == "LINE" or f.startswith("LINE"):
        return "LINE"
    if f == "SINE" or f.startswith("SINE"):
        return "SINE"
    if f == "LTR" or f.startswith("LTR"):
        return "LTR"
    if f == "Unknown" or f.startswith("Unknown"):
        return "Unclassified"

    other_types = {
        "Simple_repeat",
        "Satellite",
        "Low_complexity",
        "RNA",
        "rRNA",
        "tRNA",
        "snRNA",
        "scRNA",
        "srpRNA",
        "Satellite/centr",
        "Satellite/telo",
    }
    if f in other_types:
        return "Other"

    return "Other"

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e + 1:
            if e > cur_e:
                cur_e = e
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged

def union_bp_by_seqid(interval_dict):
    total = 0
    for seqid, intervals in interval_dict.items():
        merged = merge_intervals(intervals)
        total += sum(e - s + 1 for s, e in merged)
    return total

def fasta_lengths(path):
    lengths = {}
    name = None
    seq_len = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].strip().split()[0]
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if name is not None:
            lengths[name] = seq_len
    return lengths

genome = "Dark Reference Genome-wide"
seqid = "GENOMEWIDE"

lengths = fasta_lengths(dark_fasta)
interval_bp = sum(lengths.values())

class_intervals = defaultdict(lambda: defaultdict(list))
all_intervals = defaultdict(list)

n_records_total = 0
raw_repeat_bp_sum = 0

with open(dark_gff) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue

        rec_seqid = parts[0]
        feature = parts[2]

        try:
            start = int(parts[3])
            end = int(parts[4])
        except ValueError:
            continue

        if rec_seqid not in lengths:
            continue

        n_records_total += 1
        broad = map_feature(feature)

        class_intervals[broad][rec_seqid].append((start, end))
        all_intervals[rec_seqid].append((start, end))
        raw_repeat_bp_sum += (end - start + 1)

total_repeat_union_bp = union_bp_by_seqid(all_intervals)
nonrepeat_bp = interval_bp - total_repeat_union_bp

merged_class_bp_sum = 0
per_class_union = {}

for broad in broad_order:
    if broad == "Non-Repeat":
        continue
    bp = union_bp_by_seqid(class_intervals.get(broad, {}))
    per_class_union[broad] = bp
    merged_class_bp_sum += bp

cross_class_overlap_bp = merged_class_bp_sum - total_repeat_union_bp

diag_path = os.path.join(outdir, "dark_reference_genomewide_diagnostics.tsv")
broad_path = os.path.join(outdir, "dark_reference_genomewide_broadclass_bp.tsv")

with open(diag_path, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow([
        "genome", "seqid", "interval_bp", "gff_records_total",
        "raw_repeat_bp_sum", "merged_class_bp_sum",
        "total_repeat_union_bp", "cross_class_overlap_bp",
        "nonrepeat_bp", "repeat_pct_of_interval", "nonrepeat_pct_of_interval"
    ])
    w.writerow([
        genome,
        seqid,
        interval_bp,
        n_records_total,
        raw_repeat_bp_sum,
        merged_class_bp_sum,
        total_repeat_union_bp,
        cross_class_overlap_bp,
        nonrepeat_bp,
        round(100 * total_repeat_union_bp / interval_bp, 4),
        round(100 * nonrepeat_bp / interval_bp, 4),
    ])

with open(broad_path, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["genome", "seqid", "gwa_start", "gwa_end", "interval_bp", "broad_class", "bp_union", "pct_interval"])
    for broad in broad_order:
        if broad == "Non-Repeat":
            bp = nonrepeat_bp
        else:
            bp = per_class_union.get(broad, 0)
        pct = 100 * bp / interval_bp
        w.writerow([genome, seqid, 1, interval_bp, interval_bp, broad, bp, round(pct, 4)])

print("WROTE:", diag_path)
print("WROTE:", broad_path)
```

<a id="step-03"></a>

### 3. Combine the tables for plotting

Only source lines 1–76 are needed to create the combined input; pie-drawing functions and calls are excluded.

**Source:** `EarlGrey/TE_piecharts/03_scripts/plot_gwa_te_pies_with_genomewide.py`  
**Save as:** `plot_gwa_te_pies_with_genomewide.py`  
**Archive treatment:** Verbatim source lines 1–76; pie plotting omitted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples.

```python
#!/usr/bin/env python3
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

gwa_file = "/path/to/project/EarlGrey/TE_piecharts/01_tables/gwa_te_broadclass_bp.tsv"
genomewide_file = "/path/to/project/EarlGrey/TE_piecharts/01_tables/dark_reference_genomewide_broadclass_bp.tsv"
outdir = "/path/to/project/EarlGrey/TE_piecharts/02_plots"
tabledir = "/path/to/project/EarlGrey/TE_piecharts/01_tables"

os.makedirs(outdir, exist_ok=True)
os.makedirs(tabledir, exist_ok=True)

category_order = [
    "DNA",
    "Rolling Circle",
    "Penelope",
    "LINE",
    "SINE",
    "LTR",
    "Other",
    "Unclassified",
    "Non-Repeat",
]

plot_order = [
    "Dark Reference Genome-wide",
    "Dark Reference",
    "Plain Reference",
    "Dark01 h2",
]

colors = {
    "DNA": "#ff1f1f",
    "Rolling Circle": "#f28500",
    "Penelope": "#7b61b5",
    "LINE": "#1496d4",
    "SINE": "#b3006e",
    "LTR": "#008f2d",
    "Other": "#e3a1b2",
    "Unclassified": "#a6aaae",
    "Non-Repeat": "#000000",
}

def safe_name(x):
    return re.sub(r"[^A-Za-z0-9]+", "_", x).strip("_")

df_gwa = pd.read_csv(gwa_file, sep="\t")
df_gw = pd.read_csv(genomewide_file, sep="\t")

# Corrected wide tables for the three GWA intervals
index_cols = ["genome", "seqid", "gwa_start", "gwa_end", "interval_bp"]

wide_bp = (
    df_gwa.pivot_table(index=index_cols, columns="broad_class", values="bp_union", aggfunc="first")
    .reset_index()
)
wide_pct = (
    df_gwa.pivot_table(index=index_cols, columns="broad_class", values="pct_interval", aggfunc="first")
    .reset_index()
)

wide_bp = wide_bp[index_cols + category_order]
wide_pct = wide_pct[index_cols + category_order]

wide_bp.to_csv(os.path.join(tabledir, "gwa_te_broadclass_bp_wide.tsv"), sep="\t", index=False)
wide_pct.to_csv(os.path.join(tabledir, "gwa_te_broadclass_pct_wide.tsv"), sep="\t", index=False)

# Merge GWA + genome-wide for plotting
df = pd.concat([df_gw, df_gwa], ignore_index=True)
df["broad_class"] = pd.Categorical(df["broad_class"], categories=category_order, ordered=True)

merged_plot_table = os.path.join(tabledir, "gwa_plus_genomewide_broadclass_bp.tsv")
df.to_csv(merged_plot_table, sep="\t", index=False)
```

<a id="step-04"></a>

### 4. Draw the final donut graphics

**Source:** `EarlGrey/TE_piecharts/03_scripts/plot_gwa_te_donuts.py`  
**Save as:** `plot_gwa_te_donuts.py`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples.

```python
#!/usr/bin/env python3
import os
import re
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

infile = "/path/to/project/EarlGrey/TE_piecharts/01_tables/gwa_plus_genomewide_broadclass_bp.tsv"
outdir = "/path/to/project/EarlGrey/TE_piecharts/05_donuts"

os.makedirs(outdir, exist_ok=True)

category_order = [
    "DNA",
    "Rolling Circle",
    "Penelope",
    "LINE",
    "SINE",
    "LTR",
    "Other",
    "Unclassified",
    "Non-Repeat",
]

plot_order = [
    "Dark Reference Genome-wide",
    "Dark Reference",
    "Plain Reference",
    "Dark01 h2",
]

colors = {
    "DNA": "#ff1f1f",
    "Rolling Circle": "#f28500",
    "Penelope": "#7b61b5",
    "LINE": "#1496d4",
    "SINE": "#b3006e",
    "LTR": "#008f2d",
    "Other": "#e3a1b2",
    "Unclassified": "#a6aaae",
    "Non-Repeat": "#000000",
}

inside_label_min_pct = 4.0
outside_label_min_pct = 0.01
donut_width = 0.42

def safe_name(x):
    return re.sub(r"[^A-Za-z0-9]+", "_", x).strip("_")

df = pd.read_csv(infile, sep="\t")
df["broad_class"] = pd.Categorical(df["broad_class"], categories=category_order, ordered=True)
df["genome"] = pd.Categorical(df["genome"], categories=plot_order, ordered=True)
df = df.sort_values(["genome", "broad_class"]).copy()

def add_labels(ax, wedges, sub):
    for wedge, (_, row) in zip(wedges, sub.iterrows()):
        pct = float(row["pct_interval"])
        label = f"{pct:.1f}%"

        if pct < outside_label_min_pct:
            continue

        theta = (wedge.theta1 + wedge.theta2) / 2.0
        rad = math.radians(theta)
        x = math.cos(rad)
        y = math.sin(rad)

        if pct >= inside_label_min_pct:
            ax.text(
                0.78 * x,
                0.78 * y,
                label,
                ha="center",
                va="center",
                fontsize=10,
                color="white" if row["broad_class"] == "Non-Repeat" else "black",
                fontweight="bold",
            )
        else:
            ha = "left" if x >= 0 else "right"
            ax.annotate(
                label,
                xy=(0.98 * x, 0.98 * y),
                xytext=(1.23 * x, 1.23 * y),
                ha=ha,
                va="center",
                fontsize=9,
                arrowprops=dict(
                    arrowstyle="-",
                    lw=0.8,
                    color="black",
                    shrinkA=0,
                    shrinkB=0,
                    connectionstyle="arc3,rad=0.15",
                ),
            )

def title_long(sub):
    genome = str(sub["genome"].iloc[0])
    seqid = str(sub["seqid"].iloc[0])
    start = int(sub["gwa_start"].iloc[0])
    end = int(sub["gwa_end"].iloc[0])
    interval_bp = int(sub["interval_bp"].iloc[0])
    repeat_pct = 100.0 - float(sub.loc[sub["broad_class"] == "Non-Repeat", "pct_interval"].iloc[0])

    if genome == "Dark Reference Genome-wide":
        return f"{genome}\nwhole genome  |  total = {interval_bp:,} bp"
    return f"{genome}\n{seqid}:{start}-{end}  |  interval = {interval_bp:,} bp"

def title_short(sub):
    genome = str(sub["genome"].iloc[0])
    seqid = str(sub["seqid"].iloc[0])
    start = int(sub["gwa_start"].iloc[0])
    end = int(sub["gwa_end"].iloc[0])

    if genome == "Dark Reference Genome-wide":
        return f"{genome}\nwhole genome"
    return f"{genome}\n{seqid}:{start}-{end}"

def center_text(sub):
    repeat_pct = 100.0 - float(sub.loc[sub["broad_class"] == "Non-Repeat", "pct_interval"].iloc[0])
    return f"{repeat_pct:.1f}%\nTEs"

def legend_handles():
    return [Patch(facecolor=colors[c], edgecolor="none") for c in category_order]

def legend_labels(sub):
    labels = []
    for cat in category_order:
        row = sub[sub["broad_class"] == cat].iloc[0]
        labels.append(f"{cat} ({float(row['pct_interval']):.2f}%)")
    return labels

def plot_single(sub, outfile_base):
    sub = sub.sort_values("broad_class").copy()
    values = sub["bp_union"].tolist()
    pie_colors = [colors[c] for c in sub["broad_class"]]

    fig, ax = plt.subplots(figsize=(11.5, 7.2))
    wedges, _ = ax.pie(
        values,
        colors=pie_colors,
        startangle=90,
        counterclock=True,
        wedgeprops=dict(width=donut_width, edgecolor="white", linewidth=0.7),
    )
    add_labels(ax, wedges, sub)

    ax.text(
        0, 0, center_text(sub),
        ha="center", va="center",
        fontsize=16, fontweight="bold"
    )
    ax.set_aspect("equal")
    ax.set_title(title_long(sub), fontsize=13)

    ax.legend(
        legend_handles(),
        legend_labels(sub),
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=10.5,
        handlelength=1.2,
        handletextpad=0.5,
    )

    plt.tight_layout()
    fig.savefig(outfile_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(outfile_base + ".pdf", bbox_inches="tight")
    plt.close(fig)

def plot_combined(df_all, outfile_base):
    fig, axes = plt.subplots(1, len(plot_order), figsize=(24, 7.5))

    for ax, genome in zip(axes, plot_order):
        sub = df_all[df_all["genome"] == genome].sort_values("broad_class").copy()
        values = sub["bp_union"].tolist()
        pie_colors = [colors[c] for c in sub["broad_class"]]

        wedges, _ = ax.pie(
            values,
            colors=pie_colors,
            startangle=90,
            counterclock=True,
            wedgeprops=dict(width=donut_width, edgecolor="white", linewidth=0.7),
        )
        add_labels(ax, wedges, sub)

        ax.text(
            0, 0, center_text(sub),
            ha="center", va="center",
            fontsize=13, fontweight="bold"
        )
        ax.set_aspect("equal")
        ax.set_title(title_short(sub), fontsize=11)

    fig.legend(
        legend_handles(),
        category_order,
        loc="center right",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=11,
    )
    plt.tight_layout(rect=[0, 0, 0.90, 1])
    fig.savefig(outfile_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(outfile_base + ".pdf", bbox_inches="tight")
    plt.close(fig)

print("== Donut plotting input summary ==")
for genome in plot_order:
    sub = df[df["genome"] == genome].sort_values("broad_class")
    print(f"\n[{genome}]")
    repeat_pct = 100.0 - float(sub.loc[sub["broad_class"] == "Non-Repeat", "pct_interval"].iloc[0])
    print(f"overall_repeat_pct = {repeat_pct:.4f}")
    for _, r in sub.iterrows():
        print(f"{str(r['broad_class']):15s}  bp={int(r['bp_union']):10d}  pct={float(r['pct_interval']):7.4f}")

for genome in plot_order:
    sub = df[df["genome"] == genome].copy()
    outbase = os.path.join(outdir, f"gwa_te_donut_{safe_name(genome)}")
    plot_single(sub, outbase)
    print(f"\nWROTE: {outbase}.png")
    print(f"WROTE: {outbase}.pdf")

combined_base = os.path.join(outdir, "gwa_te_donut_combined")
plot_combined(df, combined_base)
print(f"\nWROTE: {combined_base}.png")
print(f"WROTE: {combined_base}.pdf")
```
