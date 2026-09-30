# Light reference scaffolding and gene annotation

**Paper:** Light-reference Pool-GWAS (S7); annotated haplogenomes used in Figures 2–3 and S12–S14  

Order the Light assembly against the Dark reference and transfer reference gene annotations to the Light assembly and the Dark01 hap2 representative of haplotype B.

## Inputs

- Agem_genome/NCBI/GCF_050436995.1_ilAntGemm2_primary_genomic.fna
- Agem_genome/NCBI/genomic.gff
- Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/agem_light_06.fasta
- Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold/agem_light_06.ragtag.final.fasta
- Pacbio_Reseq_Hilo/Agem_Dark_01/hifiasm_from_hfaf_t1/haplotype_fastas/Agem_Dark_01.bp.hap2.p_ctg.fa

## Outputs

- RagTag scaffold FASTA, AGP and placement/confidence files; renamed final FASTA for annotation.
- Light and Dark01-hap2 Liftoff GFF3 files, polished annotations and unmapped-feature reports.

## Software

Bash, Python 3 (standard library), RagTag, minimap2, Liftoff and SAMtools. Make these executables available on PATH.

## Execution and interpretation

1. Source assembly `agem_light_06.fasta` available through dataset, which belongs to the paper's Light 4 (Light04). This workflow starts from that supplied FASTA. The original filename is retained; see the [sample-name mapping](../README.md#sample-names).
2. RagTag produces `ragtag.scaffold.fasta`. Sections 1a–1b supply the saved name map and a new helper that generates `agem_light_06.ragtag.final.fasta` for Liftoff. All 33 retained sequences and their order are unchanged; the helper reproduces the final file byte for byte.
3. The chromosome mapping and unplaced-reference files are preserved below as small input records. Do not substitute an inferred map.
4. This gene-transfer workflow does not replace the NCBI Eukaryotic Genome Annotation Pipeline used for the original reference.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Scaffold the Light contigs

Run against the Dark reference and retain the AGP and confidence outputs.

**Source:** `Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold_agem.sh`  
**Save as:** `ragtag_scaffold_agem.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
set -euo pipefail


REF=/path/to/project/Agem_genome/NCBI/GCF_050436995.1_ilAntGemm2_primary_genomic.fna
QRY=/path/to/project/Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/agem_light_06.fasta
OUTDIR=/path/to/project/Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold

mkdir -p "$OUTDIR"

echo "=== RagTag scaffold run ==="
echo "REF: $REF"
echo "QRY: $QRY"
echo "OUTDIR: $OUTDIR"
date

ragtag.py scaffold \
  -t 40 \
  -o "$OUTDIR" \
  "$REF" \
  "$QRY"

echo
echo "=== RagTag finished ==="
date
echo
echo "Output FASTA:"
ls -lh "$OUTDIR"/ragtag.scaffold.fasta
echo
echo "AGP:"
ls -lh "$OUTDIR"/ragtag.scaffold.agp
echo
echo "Confidence / placement files:"
ls -lh "$OUTDIR" | sed -n '1,20p'
```

<a id="step-01a"></a>

### 1a. Final chromosome-name map

**Source:** `Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold/ragtag_to_final_names.tsv`  
**Save as:** `ragtag_to_final_names.tsv`  
**Archive treatment:** Recovered input table, verbatim; two tab-separated columns without a header.

This saved map gives the final chromosome names. Its 41 rows include nine possible unplaced scaffolds mapped to `chrU`; only `NW_027466813.1_RagTag` occurs in the retained RagTag FASTA. The renaming script ignores unused map rows and rejects duplicate output names among records actually present. For a different assembly with several unplaced records, assign those records distinct names in the map.

```text
NC_134745.1_RagTag	chr1
NC_134746.1_RagTag	chr2
NC_134747.1_RagTag	chr3
NC_134748.1_RagTag	chr4
NC_134749.1_RagTag	chr5
NC_134750.1_RagTag	chr6
NC_134751.1_RagTag	chr7
NC_134752.1_RagTag	chr8
NC_134753.1_RagTag	chr9
NC_134754.1_RagTag	chr10
NC_134755.1_RagTag	chr11
NC_134756.1_RagTag	chr12
NC_134757.1_RagTag	chr13
NC_134758.1_RagTag	chr14
NC_134759.1_RagTag	chr15
NC_134760.1_RagTag	chr16
NC_134761.1_RagTag	chr17
NC_134762.1_RagTag	chr18
NC_134763.1_RagTag	chr19
NC_134764.1_RagTag	chr20
NC_134765.1_RagTag	chr21
NC_134766.1_RagTag	chr22
NC_134767.1_RagTag	chr23
NC_134768.1_RagTag	chr24
NC_134769.1_RagTag	chr25
NC_134770.1_RagTag	chr26
NC_134771.1_RagTag	chr27
NC_134772.1_RagTag	chr28
NC_134773.1_RagTag	chr29
NC_134774.1_RagTag	chr30
NC_134775.1_RagTag	chrW
NC_134776.1_RagTag	chrZ
NW_027466806.1_RagTag	chrU
NW_027466807.1_RagTag	chrU
NW_027466808.1_RagTag	chrU
NW_027466809.1_RagTag	chrU
NW_027466810.1_RagTag	chrU
NW_027466811.1_RagTag	chrU
NW_027466812.1_RagTag	chrU
NW_027466813.1_RagTag	chrU
NW_027466814.1_RagTag	chrU
```

<a id="step-01b"></a>

### 1b. Rename FASTA headers before annotation

**Save as:** `rename_fasta_headers.py`  

Run `python3 rename_fasta_headers.py ragtag.scaffold.fasta ragtag_to_final_names.tsv agem_light_06.ragtag.final.fasta`, replacing the three paths as needed. The output must be a new file. Use that output as the Light target FASTA for Liftoff in section 4. This step changes FASTA identifiers only; it does not rewrite AGP or annotation files.

The helper preserves sequence bytes, record order and header descriptions. It rejects missing mappings and duplicate output identifiers before creating the output. Running it on the retained 33-record RagTag FASTA produced a byte-for-byte match to the retained final FASTA, including `chr1`–`chr30`, `chrW`, `chrZ` and `chrU`. The saved map and the observed FASTA headers agree; no sequence editing or reorientation is needed.

```python
#!/usr/bin/env python3
"""Rename FASTA sequence identifiers using a two-column, headerless name map."""
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("fasta", type=Path)
parser.add_argument("name_map", type=Path)
parser.add_argument("output", type=Path, help="New output FASTA; must not exist")
args = parser.parse_args()

names = {}
with args.name_map.open("rb") as handle:
    for number, line in enumerate(handle, 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 2 or fields[0] in names:
            raise SystemExit(f"Invalid or repeated source name on map line {number}")
        names[fields[0]] = fields[1]


def rename_header(line):
    fields = line[1:].split(maxsplit=1)
    if not fields:
        raise SystemExit("Empty FASTA header")
    old = fields[0]
    if old not in names:
        raise SystemExit(f"No mapping for {old.decode()}")
    # Keep the description, whitespace and line ending after the identifier.
    return names[old], b">" + names[old] + line[1 + len(old):]


# Validate first, so missing mappings or duplicate names produce no output file.
seen = set()
with args.fasta.open("rb") as handle:
    for line in handle:
        if line.startswith(b">"):
            new, _ = rename_header(line)
            if new in seen:
                raise SystemExit(f"Duplicate output identifier: {new.decode()}")
            seen.add(new)
if not seen:
    raise SystemExit("No FASTA records found")

# Exclusive creation protects existing files. Sequence lines are copied verbatim.
with args.fasta.open("rb") as source, args.output.open("xb") as output:
    for line in source:
        output.write(rename_header(line)[1] if line.startswith(b">") else line)
print(f"Renamed {len(seen)} FASTA records: {args.output}")
```

<a id="step-02"></a>

### 2. Chromosome correspondence input

Save as liftoff.chroms.csv; original reference names are mapped to the Light assembly names.

**Source:** `Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold/liftoff.chroms.csv`  
**Save as:** `liftoff.chroms.csv`  
**Archive treatment:** Complete source, verbatim.

```text
NC_134745.1,chr1
NC_134746.1,chr2
NC_134747.1,chr3
NC_134748.1,chr4
NC_134749.1,chr5
NC_134750.1,chr6
NC_134751.1,chr7
NC_134752.1,chr8
NC_134753.1,chr9
NC_134754.1,chr10
NC_134755.1,chr11
NC_134756.1,chr12
NC_134757.1,chr13
NC_134758.1,chr14
NC_134759.1,chr15
NC_134760.1,chr16
NC_134761.1,chr17
NC_134762.1,chr18
NC_134763.1,chr19
NC_134764.1,chr20
NC_134765.1,chr21
NC_134766.1,chr22
NC_134767.1,chr23
NC_134768.1,chr24
NC_134769.1,chr25
NC_134770.1,chr26
NC_134771.1,chr27
NC_134772.1,chr28
NC_134773.1,chr29
NC_134774.1,chr30
```

<a id="step-03"></a>

### 3. Unplaced-reference input

Save as liftoff.unplaced.txt; this is the original supplied list.

**Source:** `Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold/liftoff.unplaced.txt`  
**Save as:** `liftoff.unplaced.txt`  
**Archive treatment:** Complete source, verbatim.

```text
NC_134775.1
NC_134776.1
NW_027466806.1
NW_027466807.1
NW_027466808.1
NW_027466809.1
NW_027466810.1
NW_027466811.1
NW_027466812.1
NW_027466813.1
NW_027466814.1
```

<a id="step-04"></a>

### 4. Transfer annotations to the finalized Light assembly

Requires the finalized/renamed FASTA and both mapping inputs.

**Source:** `Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold/liftoff_full.slurm`  
**Save as:** `liftoff_full.slurm`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-32}"
set -euo pipefail

REFDIR=/path/to/project/Agem_genome/NCBI
REFFA=$REFDIR/GCF_050436995.1_ilAntGemm2_primary_genomic.fna
REFGFF=$REFDIR/genomic.gff
TGTFA=/path/to/project/Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold/agem_light_06.ragtag.final.fasta

BASEDIR=/path/to/project/Pacbio_Reseq_Hilo/Sheina_Assemblies/scaffolding/ragtag_scaffold
CHROMS=$BASEDIR/liftoff.chroms.csv
UNPLACED=$BASEDIR/liftoff.unplaced.txt

OUTDIR=$BASEDIR/liftoff_full_run
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== runtime environment ==="
date
hostname
pwd
which liftoff
which minimap2 || true

echo
echo "=== input files ==="
ls -lh "$REFFA" "$REFGFF" "$TGTFA" "$CHROMS" "$UNPLACED"

echo
echo "=== chroms file ==="
cat "$CHROMS"

echo
echo "=== unplaced file ==="
cat "$UNPLACED"

echo
echo "=== input annotation counts ==="
grep -v '^#' "$REFGFF" | awk 'BEGIN{FS=OFS="\t"} $3=="gene"{n++} END{print "gene_features", n+0}'
grep -v '^#' "$REFGFF" | awk 'BEGIN{FS=OFS="\t"} $3=="mRNA" || $3=="transcript"{n++} END{print "transcript_features", n+0}'
grep -v '^#' "$REFGFF" | awk 'BEGIN{FS=OFS="\t"} $3=="CDS"{n++} END{print "CDS_features", n+0}'

echo
echo "=== starting liftoff ==="
date

liftoff \
  -g "$REFGFF" \
  -o agem_light_06.liftoff.gff3 \
  -u agem_light_06.unmapped.txt \
  -dir intermediate_files \
  -chroms "$CHROMS" \
  -unplaced "$UNPLACED" \
  -polish \
  -p "${THREADS}" \
  "$TGTFA" "$REFFA"

echo
echo "=== liftoff finished ==="
date

echo
echo "=== output files ==="
find . -maxdepth 2 -type f | sort

echo
echo "=== output annotation counts ==="
for f in agem_light_06.liftoff.gff3 agem_light_06.liftoff.gff3_polished agem_light_06.liftoff_polished.gff3; do
  if [ -s "$f" ]; then
    echo "--- $f ---"
    grep -v '^#' "$f" | awk 'BEGIN{FS=OFS="\t"} $3=="gene"{n++} END{print "gene_features", n+0}'
    grep -v '^#' "$f" | awk 'BEGIN{FS=OFS="\t"} $3=="mRNA" || $3=="transcript"{n++} END{print "transcript_features", n+0}'
    grep -c 'partial_mapping=True' "$f" || true
    grep -c 'low_identity=True' "$f" || true
    echo "top target seqids:"
    grep -v '^#' "$f" | cut -f1 | sort | uniq -c | sort -k1,1nr | head -20
  fi
done

echo
echo "=== unmapped summary ==="
if [ -s agem_light_06.unmapped.txt ]; then
  echo "nonblank_unmapped_lines"
  grep -vc '^[[:space:]]*$' agem_light_06.unmapped.txt
  echo
  echo "first_30_unmapped"
  head -30 agem_light_06.unmapped.txt
else
  echo "0"
fi
```

<a id="step-05"></a>

### 5. Transfer annotations to Dark01 hap2

Uses the full hap2 assembly as the target; the cortex-region contig is h2tg000012l.

**Source:** `Pacbio_Reseq_Hilo/Agem_Dark_01/hifiasm_from_hfaf_t1/haplotype_fastas/liftoff_Agem_Dark_01_hap2_full.slurm`  
**Save as:** `liftoff_Agem_Dark_01_hap2_full.slurm`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-32}"
set -euo pipefail

LIFTOFF="$(command -v liftoff)"
MINIMAP2="$(command -v minimap2)"
SAMTOOLS="$(command -v samtools)"

REFDIR=/path/to/project/Agem_genome/NCBI
REFFA=$REFDIR/GCF_050436995.1_ilAntGemm2_primary_genomic.fna
REFGFF=$REFDIR/genomic.gff
TGTFA=/path/to/project/Pacbio_Reseq_Hilo/Agem_Dark_01/hifiasm_from_hfaf_t1/haplotype_fastas/Agem_Dark_01.bp.hap2.p_ctg.fa

BASEDIR=/path/to/project/Pacbio_Reseq_Hilo/Agem_Dark_01/hifiasm_from_hfaf_t1/haplotype_fastas
OUTDIR=$BASEDIR/liftoff_Agem_Dark_01_hap2_run
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== runtime environment ==="
date
hostname
pwd
echo "LIFTOFF=$LIFTOFF"
echo "MINIMAP2=$MINIMAP2"
echo "SAMTOOLS=$SAMTOOLS"
"$LIFTOFF" -h | head -5 || true
"$MINIMAP2" --version || true
"$SAMTOOLS" --version | head -3 || true

echo
echo "=== input files ==="
ls -lh "$REFFA" "$REFGFF" "$TGTFA"

echo
echo "=== ensure fasta indexes exist ==="
"$SAMTOOLS" faidx "$REFFA"
"$SAMTOOLS" faidx "$TGTFA"

echo
echo "=== target fasta summary ==="
echo -n "target_contig_count "
wc -l < "$TGTFA.fai"
head -10 "$TGTFA.fai"

echo
echo "=== input annotation counts ==="
grep -v '^#' "$REFGFF" | awk 'BEGIN{FS=OFS="\t"} $3=="gene"{n++} END{print "gene_features", n+0}'
grep -v '^#' "$REFGFF" | awk 'BEGIN{FS=OFS="\t"} $3=="mRNA" || $3=="transcript"{n++} END{print "transcript_features", n+0}'
grep -v '^#' "$REFGFF" | awk 'BEGIN{FS=OFS="\t"} $3=="CDS"{n++} END{print "CDS_features", n+0}'

echo
echo "=== starting liftoff ==="
date

"$LIFTOFF" \
  -g "$REFGFF" \
  -o Agem_Dark_01_hap2.liftoff.gff3 \
  -u Agem_Dark_01_hap2.unmapped.txt \
  -dir intermediate_files \
  -m "$MINIMAP2" \
  -polish \
  -p "${THREADS}" \
  "$TGTFA" "$REFFA"

echo
echo "=== liftoff finished ==="
date

echo
echo "=== output files ==="
find . -maxdepth 2 -type f | sort

echo
echo "=== output annotation counts ==="
for f in Agem_Dark_01_hap2.liftoff.gff3 Agem_Dark_01_hap2.liftoff_polished.gff3 Agem_Dark_01_hap2.liftoff.gff3_polished; do
  if [ -s "$f" ]; then
    echo "--- $f ---"
    grep -v '^#' "$f" | awk 'BEGIN{FS=OFS="\t"} $3=="gene"{n++} END{print "gene_features", n+0}'
    grep -v '^#' "$f" | awk 'BEGIN{FS=OFS="\t"} $3=="mRNA" || $3=="transcript"{n++} END{print "transcript_features", n+0}'
    echo "partial_mapping_count"
    grep -c 'partial_mapping=True' "$f" || true
    echo "low_identity_count"
    grep -c 'low_identity=True' "$f" || true
    echo "top_target_seqids"
    grep -v '^#' "$f" | cut -f1 | sort | uniq -c | sort -k1,1nr | head -20
  fi
done

echo
echo "=== unmapped summary ==="
if [ -s Agem_Dark_01_hap2.unmapped.txt ]; then
  echo "nonblank_unmapped_lines"
  grep -vc '^[[:space:]]*$' Agem_Dark_01_hap2.unmapped.txt
  echo
  echo "first_30_unmapped"
  head -30 Agem_Dark_01_hap2.unmapped.txt
else
  echo 0
fi
```
