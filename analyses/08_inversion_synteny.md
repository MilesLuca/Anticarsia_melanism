# Inversion-region LASTZ synteny ribbons

**Paper:** Figure 3E  

Align anchored haplotype A/B and B/C intervals with LASTZ and draw curved ribbons with gene and repeat tracks.

## Inputs

The new preparation step uses the three full-genome FASTAs, RefSeq/Liftoff gene GFFs and final EarlGrey repeat GFFs listed in section 0. It creates the following inputs for the recovered alignment and plotting scripts:

- EarlGrey/Synteny_with_TEs/NC_134752.1_anchored.fasta
- EarlGrey/Synteny_with_TEs/Chr8_Dark01_anchored.fasta
- EarlGrey/Synteny_with_TEs/chr8_anchored.fasta
- The matching composite GFFs named in lastz_ribbons_curved.R.

## Outputs

- Extracted FASTAs and .fai indexes, combined plotting GFFs, preparation_summary.tsv and samtools_version.txt.
- EarlGrey/Synteny_with_TEs/05_lastz_ribbons/01_alignments/NC134752_vs_Dark01.lastz.axt
- EarlGrey/Synteny_with_TEs/05_lastz_ribbons/01_alignments/Dark01_vs_chr8.lastz.axt
- EarlGrey/Synteny_with_TEs/05_lastz_ribbons/02_plots/threeway_lastz_ribbons_abscoords_curved.pdf

## Software

Python 3 and SAMtools for input preparation; LASTZ; Bash; R and the libraries loaded by the plotting script.

## Execution and interpretation

1. LASTZ settings are --strand=both --gapped --chain with AXT output, matching Methods.
2. Recorded absolute intervals: A 2,429,727–2,564,732; B 2,339,720–2,450,507; C 2,368,367–2,446,666.
3. Use lastz_ribbons_curved.R.
4. Section 0 extracts all three intervals in the forward orientation and converts them to uppercase, matching the saved FASTAs. The Dark reference input is soft-masked; the saved plotting extracts are uppercase. It clips overlapping annotations to the interval and converts them to local, one-based inclusive coordinates.
5. The preparation code reproduces the three plotted gene backbones (GeneIDs 142974873, 142974719 and 142974627), their exons and all overlapping EarlGrey repeats. Nine curated exon records are transcribed from the saved Geneious-source annotations. The C cortex backbone is extended to span its curated exons, matching the saved track. Other gene models are outside the saved plot selection.
6. These are flat plotting GFFs containing gene, exon and repeat features. They preserve plotted coordinates, strands and repeat classes; they do not reconstruct ignored CDS/transcript records or the editor metadata in the original exports.
7. Run section 0, then section 1, then section 2. Point section 1's ROOT at the preparation output directory. Save and run the R script from its 05_lastz_ribbons subdirectory, where section 1 writes the alignment folder.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-00"></a>

### 0. Prepare sequences and combined annotation tracks

**Save as:** `prepare_fig3e_inputs.py`

Provide the nine input files below in one directory, using these filenames (copies or links are suitable). Use the same assembly and annotation versions as analyses [02](02_light_reference_and_gene_annotation.md) and [05](05_repeat_annotation.md). The FASTAs are uncompressed and must retain the sequence identifiers in the next table.

| Haplotype | FASTA input | Gene annotation input | Repeat annotation input |
|---|---|---|---|
| A | `A.fasta`: Dark RefSeq genome, GCF_050436995.1 | `A.genes.gff`: RefSeq genomic.gff | `A.repeats.gff`: ADGO.filteredRepeats.gff |
| B | `B.fasta`: Agem_Dark_01.bp.hap2.p_ctg.fa | `B.genes.gff`: Agem_Dark_01_hap2.liftoff.gff3 | `B.repeats.gff`: AD1H2.filteredRepeats.gff |
| C | `C.fasta`: agem_light_06.ragtag.final.fasta | `C.genes.gff`: agem_light_06.liftoff.gff3 | `C.repeats.gff`: ALGO.filteredRepeats.gff |

| Haplotype | Source sequence | Inclusive interval | Extracted length | Output prefix |
|---|---|---|---|---|
| A | NC_134752.1 | 2,429,727–2,564,732 | 135,006 bp | NC_134752.1_anchored |
| B | h2tg000012l | 2,339,720–2,450,507 | 110,788 bp | Chr8_Dark01_anchored |
| C | chr8 | 2,368,367–2,446,666 | 78,300 bp | chr8_anchored |

Run `python3 prepare_fig3e_inputs.py /path/to/inputs /path/to/new_output_directory`. The output directory must be new. The script writes one FASTA, FASTA index and plotting GFF per haplotype, plus a summary with coordinates, feature counts and sequence hashes. Input indexes are reused when present; otherwise source indexes are created inside the output directory. Sequence extraction follows the [SAMtools faidx interface](https://www.htslib.org/doc/samtools-faidx.html). Annotation positions use the one-based inclusive convention described in the [GFF specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).

The constants `GENE_IDS` and `CURATED_EXONS` explicitly encode the saved plot's annotation selection. The curated exon positions are local to these exact extracted intervals; update them if adapting the code to different intervals or assemblies. The historical Geneious exports providing these records are `EarlGrey/Synteny_with_TEs/NC_134752.1_anchored.gff`, `Chr8_Dark01_anchored.gff` and `chr8_anchored.gff`. Their hashes are recorded in the source manifest.


```python
#!/usr/bin/env python3
"""New Figure 3E input preparation, based on retained sequences/annotations.

Usage: python3 prepare_fig3e_inputs.py INPUT_DIRECTORY NEW_OUTPUT_DIRECTORY
Inputs: A.fasta, A.genes.gff, A.repeats.gff, and corresponding B/C files.
Requires Python 3 and samtools on PATH. Input FASTAs are uncompressed.
"""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path
import re
import subprocess
from urllib.parse import quote, unquote

# One-based inclusive assembly coordinates; all three extracts are forward.
REGIONS = {
    "A": ("NC_134752.1", 2429727, 2564732, "NC_134752.1_anchored"),
    "B": ("h2tg000012l", 2339720, 2450507, "Chr8_Dark01_anchored"),
    "C": ("chr8", 2368367, 2446666, "chr8_anchored"),
}
# These are the three gene backbones present in the saved plotting GFFs.
GENE_IDS = {"142974873", "142974719", "142974627"}
# Curated exons transcribed from the saved Geneious-source GFF records.
# Coordinates are LOCAL to the extracted intervals, one-based inclusive.
# None denotes the standalone ivory exon, which has no plotted gene backbone.
CURATED_EXONS = {
    "A": [(None, 127613, 128122, "+", "Ivory_Exon1")],
    "B": [(None, 98547, 99056, "+", "Ivory_Exon1")],
    "C": [
        ("142974627", 49267, 49406, "-", "Cort_E9"),
        ("142974627", 49962, 50105, "-", "Cort_E8"),
        ("142974627", 51319, 51526, "-", "Cort_E7"),
        ("142974627", 52012, 52222, "-", "Cort_E6"),
        ("142974627", 52836, 52984, "-", "Cort_E5"),
        ("142974627", 58381, 58428, "-", "Cort_E1"),
        (None, 71373, 71879, "+", "Ivory_Exon1"),
    ],
}


def attribute(text, key):
    for item in text.split(";"):
        name, sep, value = item.partition("=")
        if sep and name == key:
            return unquote(value)
    return ""


def clipped_features(path, seqid, start, end):
    """Include overlaps, clip at interval edges, and shift to local bases."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            row = line.rstrip("\r\n").split("\t")
            if len(row) != 9:
                raise ValueError(f"{path}:{number}: expected nine GFF columns")
            if row[0] != seqid:
                continue
            left, right = int(row[3]), int(row[4])
            if left < 1 or right < left:
                raise ValueError(f"{path}:{number}: invalid feature coordinates")
            if left > end or right < start:
                continue
            yield {"source": row[1], "type": row[2],
                   "start": max(left, start) - start + 1,
                   "end": min(right, end) - start + 1,
                   "strand": row[6], "attributes": row[8]}


def prepare_annotations(haplotype, genes_path, repeats_path):
    seqid, start, end, prefix = REGIONS[haplotype]
    genes, exons = {}, []
    for feature in clipped_features(genes_path, seqid, start, end):
        if feature["type"] not in {"gene", "exon"}:
            continue
        matches = re.findall(r"GeneID:(\d+)", attribute(feature["attributes"], "Dbxref"))
        selected = GENE_IDS.intersection(matches)
        if not selected:
            continue
        if len(selected) != 1:
            raise ValueError("Ambiguous gene assignment in input annotation")
        gene_id = selected.pop()
        feature["gene_id"] = gene_id
        feature["name"] = attribute(feature["attributes"], "Name") or f"LOC{gene_id}"
        if feature["type"] == "gene":
            if gene_id in genes:
                raise ValueError(f"Multiple models for GeneID:{gene_id} in interval")
            genes[gene_id] = feature
        else:
            exons.append(feature)
    if set(genes) != GENE_IDS:
        raise ValueError(f"{haplotype}: missing expected gene backbones")
    for gene_id, left, right, strand, name in CURATED_EXONS[haplotype]:
        # Avoid adding a second copy if the input already contains this exon.
        if any((f["start"], f["end"], f["strand"]) == (left, right, strand) for f in exons):
            continue
        exons.append({"source": "Geneious", "type": "exon", "start": left,
                      "end": right, "strand": strand, "gene_id": gene_id, "name": name})
    # The saved C cortex backbone spans the curated exons as well as Liftoff exons.
    for gene_id, gene in genes.items():
        children = [f for f in exons if f["gene_id"] == gene_id]
        if any(f["strand"] != gene["strand"] for f in children):
            raise ValueError(f"{haplotype}: inconsistent exon/gene strands")
        gene["start"] = min([gene["start"]] + [f["start"] for f in children])
        gene["end"] = max([gene["end"]] + [f["end"] for f in children])
    repeats = list(clipped_features(repeats_path, seqid, start, end))
    if not repeats or any(f["source"] != "Earl_Grey" for f in repeats):
        raise ValueError(f"{haplotype}: expected overlapping Earl_Grey repeat records")
    for feature in repeats:
        feature["name"] = attribute(feature["attributes"], "Name") or feature["type"]
    features = list(genes.values()) + exons + repeats
    if any(not 1 <= f["start"] <= f["end"] <= end-start+1 for f in features):
        raise ValueError(f"{haplotype}: annotation outside extracted sequence")
    return features, (len(genes), len(exons), len(repeats))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_directory", type=Path)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    input_dir = args.input_directory.resolve()
    out = args.output_directory.resolve()
    for hap in REGIONS:
        for suffix in ("fasta", "genes.gff", "repeats.gff"):
            if not (input_dir / f"{hap}.{suffix}").is_file():
                parser.error(f"Missing input: {input_dir / f'{hap}.{suffix}'}")
    out.mkdir(parents=True, exist_ok=False)
    summary = []
    for hap, (seqid, start, end, prefix) in REGIONS.items():
        fasta = (input_dir / f"{hap}.fasta").resolve()
        # Reuse an existing index; otherwise create one inside the output folder.
        index = Path(str(fasta) + ".fai")
        if not index.is_file():
            index = out / f"{hap}.source.fai"
            subprocess.run(["samtools", "faidx", "--fai-idx", str(index), str(fasta)], check=True)
        result = subprocess.run(["samtools", "faidx", "--fai-idx", str(index),
                                 str(fasta), f"{seqid}:{start}-{end}"],
                                check=True, capture_output=True, text=True)
        # Saved plotting FASTAs are uppercase, including the soft-masked A source.
        sequence = "".join(result.stdout.splitlines()[1:]).upper()
        if len(sequence) != end-start+1:
            raise ValueError(f"{hap}: extracted sequence length disagrees with interval")
        output_fasta = out / f"{prefix}.fasta"
        with output_fasta.open("w") as handle:
            handle.write(f">{prefix}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i:i+80] + "\n")
        subprocess.run(["samtools", "faidx", str(output_fasta)], check=True)
        features, counts = prepare_annotations(hap, input_dir / f"{hap}.genes.gff",
                                               input_dir / f"{hap}.repeats.gff")
        # Flat plotting GFF: only gene/exon/repeat features used by the R plotter.
        # CDS phases and transcript hierarchies are deliberately not reconstructed.
        with (out / f"{prefix}.gff").open("w") as handle:
            handle.write(f"##sequence-region {prefix} 1 {len(sequence)}\n")
            for number, feature in enumerate(features, 1):
                name = quote(feature["name"], safe=" :._-")
                attr = f"ID={hap}_plot_{number};Name={name}"
                handle.write("\t".join(map(str, [prefix, feature["source"], feature["type"],
                    feature["start"], feature["end"], ".", feature["strand"], ".", attr])) + "\n")
        summary.append([hap, seqid, start, end, "+", len(sequence), *counts,
                        hashlib.sha256(sequence.upper().encode()).hexdigest()])
    with (out / "preparation_summary.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["haplotype", "source_seqid", "start", "end", "orientation",
                         "length", "genes", "exons", "repeats", "sequence_sha256"])
        writer.writerows(summary)
    version = subprocess.run(["samtools", "--version"], check=True, capture_output=True, text=True)
    (out / "samtools_version.txt").write_text(version.stdout)
    print(f"Prepared Figure 3E FASTA/index/GFF files in {out}")


if __name__ == "__main__":
    main()
```

<a id="step-01"></a>

### 1. Generate pairwise chained alignments

Produces both AXT files required by the plotter.

**Source:** `EarlGrey/Synteny_with_TEs/run_lastz_ribbons.sh`  
**Save as:** `run_lastz_ribbons.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults.

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT=/path/to/project/EarlGrey/Synteny_with_TEs
OUTROOT="$ROOT/05_lastz_ribbons"
ALIGNDIR="$OUTROOT/01_alignments"
PLOTDIR="$OUTROOT/02_plots"

REF_TOP="$ROOT/NC_134752.1_anchored.fasta"
REF_MID="$ROOT/Chr8_Dark01_anchored.fasta"
REF_BOT="$ROOT/chr8_anchored.fasta"

mkdir -p "$ALIGNDIR" "$PLOTDIR"
cd "$ALIGNDIR"

{
  echo "=================================================="
  echo "LASTZ ribbon alignment run"
  echo "PWD: $(pwd)"
  echo "TOP: $REF_TOP"
  echo "MID: $REF_MID"
  echo "BOT: $REF_BOT"
  echo "=================================================="
  echo
} > diagnostics.txt

echo "[1/8] Checking FASTA inputs"
ls -lh "$REF_TOP" "$REF_MID" "$REF_BOT"

# --------------------------------------------------
# 1. NC_134752.1_anchored vs Chr8_Dark01_anchored
# --------------------------------------------------
echo "[2/8] Running lastz: NC_134752.1_anchored vs Chr8_Dark01_anchored"

lastz \
  "$REF_TOP" \
  "$REF_MID" \
  --strand=both \
  --gapped \
  --chain \
  --axt=NC134752_vs_Dark01.lastz.axt

echo "[3/8] Summarizing NC134752_vs_Dark01.lastz.axt"

{
  echo "NC134752_vs_Dark01.lastz.axt"
  awk '
    NF>0 && $1 ~ /^[0-9]+$/ {
      n++
      len = $4 - $3 + 1
      if ($8 == "+") f++
      else if ($8 == "-") r++
    }
    END {
      print "total_blocks\t" (n+0)
      print "forward\t" (f+0)
      print "reverse\t" (r+0)
    }
  ' NC134752_vs_Dark01.lastz.axt

  echo "-- top 20 longest --"
  awk '
    NF>0 && $1 ~ /^[0-9]+$/ {
      len = $4 - $3 + 1
      print len "\t" $0
    }
  ' NC134752_vs_Dark01.lastz.axt | sort -k1,1nr | awk 'NR<=20'
  echo
} >> diagnostics.txt

# --------------------------------------------------
# 2. Chr8_Dark01_anchored vs chr8_anchored
# --------------------------------------------------
echo "[4/8] Running lastz: Chr8_Dark01_anchored vs chr8_anchored"

lastz \
  "$REF_MID" \
  "$REF_BOT" \
  --strand=both \
  --gapped \
  --chain \
  --axt=Dark01_vs_chr8.lastz.axt

echo "[5/8] Summarizing Dark01_vs_chr8.lastz.axt"

{
  echo "Dark01_vs_chr8.lastz.axt"
  awk '
    NF>0 && $1 ~ /^[0-9]+$/ {
      n++
      len = $4 - $3 + 1
      if ($8 == "+") f++
      else if ($8 == "-") r++
    }
    END {
      print "total_blocks\t" (n+0)
      print "forward\t" (f+0)
      print "reverse\t" (r+0)
    }
  ' Dark01_vs_chr8.lastz.axt

  echo "-- top 20 longest --"
  awk '
    NF>0 && $1 ~ /^[0-9]+$/ {
      len = $4 - $3 + 1
      print len "\t" $0
    }
  ' Dark01_vs_chr8.lastz.axt | sort -k1,1nr | awk 'NR<=20'
  echo
} >> diagnostics.txt

echo "[6/8] Final diagnostics"
cat diagnostics.txt

echo "[7/8] Output files"
ls -lh *.axt diagnostics.txt

echo "[8/8] Done"
```

<a id="step-02"></a>

### 2. Render curved ribbons and annotations

Draws curved ribbons using the Figure 3E coordinate system.

**Source:** `EarlGrey/Synteny_with_TEs/05_lastz_ribbons/lastz_ribbons_curved.R`  
**Save as:** `lastz_ribbons_curved.R`  
**Archive treatment:** Complete source, verbatim.

```r
library(ggplot2)

# =========================================================
# inputs
# =========================================================
fai_top <- "../NC_134752.1_anchored.fasta.fai"
fai_mid <- "../Chr8_Dark01_anchored.fasta.fai"
fai_bot <- "../chr8_anchored.fasta.fai"

gff_top <- "../NC_134752.1_anchored.gff"
gff_mid <- "../Chr8_Dark01_anchored.gff"
gff_bot <- "../chr8_anchored.gff"

axt_top_mid <- "01_alignments/NC134752_vs_Dark01.lastz.axt"
axt_mid_bot <- "01_alignments/Dark01_vs_chr8.lastz.axt"

dir.create("02_plots", showWarnings = FALSE)

# =========================================================
# absolute coordinate ranges
# =========================================================
top_abs_start <- 2429727
top_abs_end   <- 2564732

mid_abs_start <- 2339720
mid_abs_end   <- 2450507

bot_abs_start <- 2368367
bot_abs_end   <- 2446666

# =========================================================
# helpers
# =========================================================
read_fai_len <- function(fai_file) {
  x <- read.table(fai_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  as.numeric(x[1, 2])
}

read_gff <- function(gff_file) {
  g <- read.table(
    gff_file,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  colnames(g) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
  g$start <- as.numeric(g$start)
  g$end   <- as.numeric(g$end)
  g
}

# AXT headers:
# id ref r_start r_end qry q_start q_end strand score
# For strand == "-", q_start/q_end are on the reverse-complement axis of query.
# Convert them back to forward genomic coordinates using query_size.
read_axt_headers <- function(axt_file, query_size) {
  lines <- trimws(readLines(axt_file))
  hdr <- grep(
    "^[0-9]+\\s+\\S+\\s+[0-9]+\\s+[0-9]+\\s+\\S+\\s+[0-9]+\\s+[0-9]+\\s+[+-]\\s+[0-9]+$",
    lines,
    value = TRUE
  )

  parts <- strsplit(hdr, "\\s+")
  df <- do.call(rbind, lapply(parts, function(x) {
    data.frame(
      block_id = as.integer(x[1]),
      ref      = x[2],
      r_start  = as.numeric(x[3]),
      r_end    = as.numeric(x[4]),
      qry      = x[5],
      q_start  = as.numeric(x[6]),
      q_end    = as.numeric(x[7]),
      strand   = x[8],
      score    = as.numeric(x[9]),
      stringsAsFactors = FALSE
    )
  }))

  df$orientation <- ifelse(df$strand == "+", "forward", "reverse")
  df$r_len <- df$r_end - df$r_start + 1
  df$q_len <- abs(df$q_end - df$q_start) + 1

  df$q_start_fwd <- ifelse(
    df$strand == "+",
    df$q_start,
    query_size - df$q_end + 1
  )

  df$q_end_fwd <- ifelse(
    df$strand == "+",
    df$q_end,
    query_size - df$q_start + 1
  )

  df
}

norm_start <- function(pos, len) (pos - 1) / len
norm_end   <- function(pos, len) pos / len

abs_to_norm <- function(abs_pos, abs_start, abs_end) {
  (abs_pos - abs_start) / (abs_end - abs_start)
}

pretty_ticks_with_ends <- function(abs_start, abs_end, n = 5) {
  ticks <- pretty(c(abs_start, abs_end), n = n)
  ticks <- ticks[ticks >= abs_start & ticks <= abs_end]
  ticks <- unique(c(abs_start, ticks, abs_end))
  sort(ticks)
}

build_axis_df <- function(abs_start, abs_end, n = 5) {
  ticks <- pretty_ticks_with_ends(abs_start, abs_end, n = n)
  data.frame(
    abs = ticks,
    x   = abs_to_norm(ticks, abs_start, abs_end),
    label = format(ticks, big.mark = ",", scientific = FALSE, trim = TRUE),
    stringsAsFactors = FALSE
  )
}

inherit_exon_strand <- function(exons, genes) {
  bad <- !(exons$strand %in% c("+", "-"))
  if (!any(bad)) return(exons)

  for (i in which(bad)) {
    ov <- genes[genes$start <= exons$end[i] & genes$end >= exons$start[i], ]
    if (nrow(ov) == 0) next
    overlaps <- pmin(exons$end[i], ov$end) - pmax(exons$start[i], ov$start) + 1
    overlaps[overlaps < 0] <- 0
    if (all(overlaps == 0)) next
    exons$strand[i] <- ov$strand[which.max(overlaps)]
  }
  exons
}

make_feature_rects <- function(df, seq_len) {
  out <- df
  out$xmin <- norm_start(out$start, seq_len)
  out$xmax <- norm_end(out$end, seq_len)
  out
}

make_gene_backbones <- function(genes, seq_len) {
  out <- genes
  out$x    <- norm_start(out$start, seq_len)
  out$xend <- norm_end(out$end, seq_len)
  out
}

add_const_col <- function(df, colname, value) {
  df[[colname]] <- rep(value, nrow(df))
  df
}

# =========================================================
# TE classification
# =========================================================
classify_te <- function(type_vec) {
  out <- rep("Other", length(type_vec))

  out[grepl("^RC/Helitron", type_vec)] <- "Helitron"

  out[grepl("^SINE/5S$", type_vec)]    <- "5S"
  out[grepl("^SINE/tRNA$", type_vec)]  <- "tRNA"
  out[grepl("^SINE", type_vec) & out == "Other"] <- "SINE"

  out[grepl("^DNA/TcMar", type_vec)]    <- "TcMar"
  out[grepl("^DNA/CMC", type_vec)]      <- "CMC"
  out[grepl("^DNA/hAT", type_vec)]      <- "hAT"
  out[grepl("^DNA/MULE", type_vec)]     <- "MULE"
  out[grepl("^DNA/PIF", type_vec)]      <- "PIF"
  out[grepl("^DNA/PiggyBac", type_vec)] <- "PiggyBac"
  out[grepl("^DNA/Maverick", type_vec)] <- "Maverick"
  out[grepl("^DNA/Academ", type_vec)]   <- "Academ"
  out[grepl("^DNA/Sola", type_vec)]     <- "Sola"
  out[grepl("^DNA/Zator", type_vec)]    <- "Zator"
  out[grepl("^DNA", type_vec) & out == "Other"] <- "DNA"

  out[grepl("^LINE/CR1", type_vec)]    <- "CR1"
  out[grepl("^LINE/Dong", type_vec)]   <- "Dong"
  out[grepl("^LINE/I", type_vec)]      <- "I"
  out[grepl("^LINE/L1", type_vec)]     <- "L1"
  out[grepl("^LINE/L2", type_vec)]     <- "L2"
  out[grepl("^LINE/Proto2", type_vec)] <- "Proto2"
  out[grepl("^LINE/R1", type_vec)]     <- "R1"
  out[grepl("^LINE/R2", type_vec)]     <- "R2"
  out[grepl("^LINE/RTE", type_vec)]    <- "RTE"
  out[grepl("^LINE", type_vec) & out == "Other"] <- "LINE"

  out[grepl("^LTR/Copia", type_vec)]   <- "Copia"
  out[grepl("^LTR/Gypsy", type_vec)]   <- "Gypsy"
  out[grepl("^LTR/Pao", type_vec)]     <- "Pao"
  out[grepl("^LTR", type_vec) & out == "Other"] <- "LTR"

  out[grepl("^PLE", type_vec)]         <- "PLE"

  out[type_vec == "Unknown"]           <- "Unknown"
  out[type_vec == "Simple_repeat"]     <- "Simple_repeat"
  out[type_vec == "Satellite"]         <- "Satellite"
  out[type_vec == "Low_complexity"]    <- "Low_complexity"

  out
}

te_palette <- c(
  "Academ"        = "#F8766D",
  "CMC"           = "#DB8E00",
  "DNA"           = "#AEA200",
  "hAT"           = "#64B200",
  "Maverick"      = "#00BD5C",
  "MULE"          = "#00C1A7",
  "PIF"           = "#00A6FF",
  "PiggyBac"      = "#00BADE",
  "Sola"          = "#B385FF",
  "TcMar"         = "#EF67EB",
  "Zator"         = "#FF63B6",

  "CR1"           = "#EFF3FF",
  "Dong"          = "#08306B",
  "I"             = "#2171B5",
  "L1"            = "#4292C6",
  "L2"            = "#6BAED6",
  "LINE"          = "#9ECAE1",
  "Proto2"        = "#C6DBEF",
  "R1"            = "#DEEBF7",
  "R2"            = "#F7FBFF",
  "RTE"           = "#DEEBF7",

  "Copia"         = "#238B45",
  "Gypsy"         = "#74C476",
  "LTR"           = "#BAE4B3",
  "Pao"           = "#EDF8E9",

  "5S"            = "#F03B20",
  "SINE"          = "#FEB24C",
  "tRNA"          = "#FFEDA0",

  "PLE"           = "#EFEDF5",
  "Helitron"      = "#FEE6CE",

  "Unknown"       = "#BDBDBD",
  "Satellite"     = "#969696",
  "Simple_repeat" = "#E0E0E0",
  "Low_complexity"= "#F2F2F2",
  "Other"         = "#CCCCCC"
)

legend_order <- c(
  "TcMar", "Helitron", "5S", "SINE", "tRNA",
  "CR1", "Dong", "I", "L1", "L2", "LINE", "Proto2", "R1", "R2", "RTE",
  "Copia", "Gypsy", "LTR", "Pao",
  "Academ", "CMC", "DNA", "hAT", "Maverick", "MULE", "PIF", "PiggyBac", "Sola", "Zator",
  "PLE",
  "Unknown", "Satellite", "Simple_repeat", "Low_complexity", "Other"
)

# =========================================================
# curved ribbon helpers
# =========================================================
bezier_curve <- function(p0, p1, p2, p3, n = 50) {
  t <- seq(0, 1, length.out = n)
  x <- (1 - t)^3 * p0[1] +
       3 * (1 - t)^2 * t * p1[1] +
       3 * (1 - t) * t^2 * p2[1] +
       t^3 * p3[1]
  y <- (1 - t)^3 * p0[2] +
       3 * (1 - t)^2 * t * p1[2] +
       3 * (1 - t) * t^2 * p2[2] +
       t^3 * p3[2]
  data.frame(x = x, y = y)
}

# curve_strength:
#   larger values = more bowing / pinching
# n:
#   number of interpolated points per side
make_curved_ribbon_polygons <- function(coords, ref_len, qry_len, y_ref, y_qry,
                                        pair_name, curve_strength = 0.22, n = 60) {
  polys <- vector("list", nrow(coords))

  for (i in seq_len(nrow(coords))) {
    x1 <- norm_start(coords$r_start[i], ref_len)
    x2 <- norm_end(coords$r_end[i], ref_len)

    q1 <- norm_start(coords$q_start_fwd[i], qry_len)
    q2 <- norm_end(coords$q_end_fwd[i], qry_len)

    # Which bottom x corresponds to ref start and ref end?
    if (coords$orientation[i] == "forward") {
      bottom_for_ref_start <- q1
      bottom_for_ref_end   <- q2
    } else {
      bottom_for_ref_start <- q2
      bottom_for_ref_end   <- q1
    }

    dy <- abs(y_ref - y_qry)
    bow <- dy * curve_strength

    # Side from top-left to mapped bottom position for ref start
    side1 <- bezier_curve(
      p0 = c(x1, y_ref),
      p1 = c(x1, y_ref - bow),
      p2 = c(bottom_for_ref_start, y_qry + bow),
      p3 = c(bottom_for_ref_start, y_qry),
      n = n
    )

    # Straight bottom edge from mapped start to mapped end
    bottom_edge <- data.frame(
      x = seq(bottom_for_ref_start, bottom_for_ref_end, length.out = 8),
      y = rep(y_qry, 8)
    )

    # Side from mapped bottom position for ref end back to top-right
    side2 <- bezier_curve(
      p0 = c(bottom_for_ref_end, y_qry),
      p1 = c(bottom_for_ref_end, y_qry + bow),
      p2 = c(x2, y_ref - bow),
      p3 = c(x2, y_ref),
      n = n
    )

    # Straight top edge back to top-left
    top_edge <- data.frame(
      x = seq(x2, x1, length.out = 8),
      y = rep(y_ref, 8)
    )

    poly <- rbind(side1, bottom_edge, side2, top_edge)

    poly$block_id <- coords$block_id[i]
    poly$pair <- pair_name
    poly$orientation <- coords$orientation[i]
    poly$poly_id <- paste0(pair_name, "_", coords$block_id[i])

    polys[[i]] <- poly
  }

  do.call(rbind, polys)
}

prepare_track <- function(gff, seq_len, y_base) {
  genes <- subset(gff, type == "gene")
  exons <- subset(gff, type == "exon")
  exons <- inherit_exon_strand(exons, genes)

  tes <- subset(gff, source == "Earl_Grey")
  tes$te_class <- classify_te(tes$type)

  genes_plus  <- subset(genes, strand == "+")
  genes_minus <- subset(genes, strand == "-")
  exons_plus  <- subset(exons, strand == "+")
  exons_minus <- subset(exons, strand == "-")

  gene_plus_df   <- add_const_col(make_gene_backbones(genes_plus, seq_len), "y", y_base + 0.20)
  gene_minus_df  <- add_const_col(make_gene_backbones(genes_minus, seq_len), "y", y_base + 0.00)

  exon_plus_df   <- add_const_col(add_const_col(make_feature_rects(exons_plus, seq_len), "ymin", y_base + 0.165), "ymax", y_base + 0.235)
  exon_minus_df  <- add_const_col(add_const_col(make_feature_rects(exons_minus, seq_len), "ymin", y_base - 0.035), "ymax", y_base + 0.035)

  te_df          <- add_const_col(add_const_col(make_feature_rects(tes, seq_len), "ymin", y_base + 0.075), "ymax", y_base + 0.135)

  list(
    te_df = te_df,
    gene_plus_df = gene_plus_df,
    gene_minus_df = gene_minus_df,
    exon_plus_df = exon_plus_df,
    exon_minus_df = exon_minus_df
  )
}

# =========================================================
# load data
# =========================================================
len_top <- read_fai_len(fai_top)
len_mid <- read_fai_len(fai_mid)
len_bot <- read_fai_len(fai_bot)

coords_tm <- read_axt_headers(axt_top_mid, query_size = len_mid)
coords_mb <- read_axt_headers(axt_mid_bot, query_size = len_bot)

track_top <- prepare_track(read_gff(gff_top), len_top, 2.80)
track_mid <- prepare_track(read_gff(gff_mid), len_mid, 1.50)
track_bot <- prepare_track(read_gff(gff_bot), len_bot, 0.20)

top_axis_df <- build_axis_df(top_abs_start, top_abs_end, n = 6)
mid_axis_df <- build_axis_df(mid_abs_start, mid_abs_end, n = 6)
bot_axis_df <- build_axis_df(bot_abs_start, bot_abs_end, n = 6)

cat("Length checks\n")
cat("-------------\n")
cat("Top abs span: ", top_abs_end - top_abs_start + 1, " | fasta length: ", len_top, "\n")
cat("Mid abs span: ", mid_abs_end - mid_abs_start + 1, " | fasta length: ", len_mid, "\n")
cat("Bot abs span: ", bot_abs_end - bot_abs_start + 1, " | fasta length: ", len_bot, "\n\n")

cat("LASTZ block counts\n")
cat("------------------\n")
cat("Top-mid blocks:", nrow(coords_tm), " forward:", sum(coords_tm$orientation == "forward"),
    " reverse:", sum(coords_tm$orientation == "reverse"), "\n")
cat("Mid-bot blocks:", nrow(coords_mb), " forward:", sum(coords_mb$orientation == "forward"),
    " reverse:", sum(coords_mb$orientation == "reverse"), "\n")

# =========================================================
# y layout
# =========================================================
y_top <- 2.80
y_mid <- 1.50
y_bot <- 0.20

y_tm_top <- 2.55
y_tm_bot <- 1.75
y_mb_top <- 1.25
y_mb_bot <- 0.45

# Build curved ribbons
ribbons_tm <- make_curved_ribbon_polygons(
  coords_tm, len_top, len_mid, y_tm_top, y_tm_bot,
  pair_name = "top_mid",
  curve_strength = 0.22,
  n = 60
)

ribbons_mb <- make_curved_ribbon_polygons(
  coords_mb, len_mid, len_bot, y_mb_top, y_mb_bot,
  pair_name = "mid_bot",
  curve_strength = 0.22,
  n = 60
)

ribbons_all <- rbind(ribbons_tm, ribbons_mb)

present_te_classes <- unique(c(track_top$te_df$te_class, track_mid$te_df$te_class, track_bot$te_df$te_class))
present_te_classes <- legend_order[legend_order %in% present_te_classes]

# =========================================================
# plot
# =========================================================
p <- ggplot() +
  # forward ribbons
  geom_polygon(
    data = subset(ribbons_all, orientation == "forward"),
    aes(x = x, y = y, group = poly_id),
    fill = "grey55",
    color = NA,
    alpha = 0.22
  ) +
  # reverse ribbons
  geom_polygon(
    data = subset(ribbons_all, orientation == "reverse"),
    aes(x = x, y = y, group = poly_id),
    fill = "#F28E2B",
    color = NA,
    alpha = 0.35
  ) +

  # track backgrounds
  annotate("rect", xmin = 0, xmax = 1, ymin = y_top - 0.08, ymax = y_top + 0.28,
           fill = "grey92", color = "grey75") +
  annotate("rect", xmin = 0, xmax = 1, ymin = y_mid - 0.08, ymax = y_mid + 0.28,
           fill = "grey92", color = "grey75") +
  annotate("rect", xmin = 0, xmax = 1, ymin = y_bot - 0.08, ymax = y_bot + 0.28,
           fill = "grey92", color = "grey75") +

  # top absolute axis
  geom_segment(
    data = top_axis_df,
    aes(x = x, xend = x, y = 3.08, yend = 3.04),
    inherit.aes = FALSE,
    linewidth = 0.28, color = "black"
  ) +
  geom_text(
    data = top_axis_df,
    aes(x = x, y = 3.12, label = label),
    inherit.aes = FALSE,
    size = 3, vjust = 0
  ) +

  # middle absolute axis
  geom_segment(
    data = mid_axis_df,
    aes(x = x, xend = x, y = 1.38, yend = 1.42),
    inherit.aes = FALSE,
    linewidth = 0.28, color = "black"
  ) +
  geom_text(
    data = mid_axis_df,
    aes(x = x, y = 1.34, label = label),
    inherit.aes = FALSE,
    size = 3, vjust = 1
  ) +

  # bottom absolute axis
  geom_segment(
    data = bot_axis_df,
    aes(x = x, xend = x, y = -0.085, yend = -0.045),
    inherit.aes = FALSE,
    linewidth = 0.28, color = "black"
  ) +
  geom_text(
    data = bot_axis_df,
    aes(x = x, y = -0.13, label = label),
    inherit.aes = FALSE,
    size = 3, vjust = 1
  ) +

  # gene backbones
  geom_segment(data = track_top$gene_plus_df,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.35, color = "black") +
  geom_segment(data = track_top$gene_minus_df,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.35, color = "black") +
  geom_segment(data = track_mid$gene_plus_df,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.35, color = "black") +
  geom_segment(data = track_mid$gene_minus_df,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.35, color = "black") +
  geom_segment(data = track_bot$gene_plus_df,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.35, color = "black") +
  geom_segment(data = track_bot$gene_minus_df,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.35, color = "black") +

  # TE lanes
  geom_rect(data = track_top$te_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = te_class),
            color = NA) +
  geom_rect(data = track_mid$te_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = te_class),
            color = NA) +
  geom_rect(data = track_bot$te_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = te_class),
            color = NA) +

  # exon lanes
  geom_rect(data = track_top$exon_plus_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "black", color = NA) +
  geom_rect(data = track_top$exon_minus_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "black", color = NA) +
  geom_rect(data = track_mid$exon_plus_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "black", color = NA) +
  geom_rect(data = track_mid$exon_minus_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "black", color = NA) +
  geom_rect(data = track_bot$exon_plus_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "black", color = NA) +
  geom_rect(data = track_bot$exon_minus_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "black", color = NA) +

  # labels
  annotate("text", x = 0.00, y = y_top + 0.255, label = "5'", hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = 1.00, y = y_top + 0.255, label = "3'", hjust = 1, vjust = 1, size = 4) +
  annotate("text", x = 0.00, y = y_mid + 0.255, label = "5'", hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = 1.00, y = y_mid + 0.255, label = "3'", hjust = 1, vjust = 1, size = 4) +
  annotate("text", x = 0.00, y = y_bot + 0.255, label = "5'", hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = 1.00, y = y_bot + 0.255, label = "3'", hjust = 1, vjust = 1, size = 4) +

  annotate("text", x = -0.045, y = y_top + 0.10, label = "NC_134752.1_anchored", hjust = 0, size = 4) +
  annotate("text", x = -0.045, y = y_mid + 0.10, label = "Chr8_Dark01_anchored", hjust = 0, size = 4) +
  annotate("text", x = -0.045, y = y_bot + 0.10, label = "chr8_anchored", hjust = 0, size = 4) +

  annotate("text", x = 0.00, y = 3.22, label = "Three-track synteny ribbons from LASTZ alignments", hjust = 0, vjust = 1, size = 4.4, fontface = "bold") +

  scale_fill_manual(
    values = te_palette,
    breaks = present_te_classes,
    drop = FALSE,
    name = "TE class"
  ) +

  coord_cartesian(xlim = c(-0.06, 1), ylim = c(-0.16, 3.24), expand = FALSE, clip = "off") +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(18, 20, 32, 40),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  guides(
    fill = guide_legend(
      nrow = 3,
      byrow = TRUE,
      override.aes = list(color = NA)
    )
  )

ggsave("02_plots/threeway_lastz_ribbons_abscoords_curved.pdf", p, width = 14, height = 8.8)
```
