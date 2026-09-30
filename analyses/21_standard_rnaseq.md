# Standard RNA-seq alignment and ivory-region extraction


**Paper:** Figure S5A, RNA-seq alignment input  

Align the standard RNA-seq libraries together with STAR, retain mapped reads flagged as properly paired with mapping quality at least 20, and extract the cortex/ivory region.

## Inputs

- Reference FASTA: `Anticarsia_gemmatalis_primary_pseudohap_scaffold.fasta`, containing `Chromosome8`.
- Paired gzip-compressed FASTQs with suffixes `_R1.fastq.gz` and `_R2.fastq.gz` for AG09, AG10, AG11, AG12, AG13, AG16, AG17, AG18, AG19, AG20, Ag04, Ag05 and Ag06. Preserve capitalization and the matching order of the two mate lists.

## Outputs

- STAR genome index.
- `AgemnoindexAligned.sortedByCoord.out.bam` and STAR alignment logs.
- `Agem.filtered.bam` and its index.
- `Agem.filtered.chr8.ivory.bam` and its index, containing alignments overlapping `Chromosome8:2410000-2850000` (one-based, inclusive).

## Software

STAR `2.7.10a_alpha_220207` and SAMtools `1.15.1` are recorded in the selected BAM header; the STAR version also appears in the completed alignment log. These versions concern standard RNA-seq. The separate small-RNA workflow in analysis 18 records STAR 2.7.11b. Bash and gzip/zcat are required, with the tools available on PATH.

## Execution and interpretation

1. Set the same `ROOT` in all three scripts. The example layout places the reference in `Agem_RNA-Seq_Dec2024`, the FASTQs in its `Agem/fastqs` subdirectory, and outputs in `Agem/results_combined`. Edit these paths for your own files. Save each block under its stated filename and run sections 1–3 in order with Bash. Use a new output directory for a fresh run.
2. Section 1 recovers the index-generation command recorded inside the STAR alignment log: `genomeSAindexNbases 12`, with no annotation argument. Section 2 preserves the combined-library alignment parameters and both explicit FASTQ lists from `STAR_map.sh`.
3. Section 3 restores the filtering command that is commented out in the saved `samtools.sh` but recorded as executed in the BAM header: `-F 0x04 -f 0x2 -q 20`. It also restores indexing of the filtered BAM before region extraction. The filter operates on alignment records; it does not separately remove secondary, supplementary or duplicate-flagged alignments.
4. The default thread counts remain 40 for STAR and 32 for SAMtools. Set the `THREADS` environment variable to use a different count. The STAR sorting-memory limit remains the recorded 23,188,627,292 bytes.
5. The retained `Chromosome8` sequence is 14,130,858 bp and matches `NC_134752.1` from the paper's RefSeq assembly after conversion to uppercase. The uppercase sequence SHA-256 is `cfe67f13ab533eb851f9e0e1ec766f631d8cd18f24174343b46be67e13ec7d27`. This verifies chromosome-8 coordinates and orientation; it is not a comparison of the entire assemblies. Use a reference, STAR index and region identifier that agree with one another.

## Code

<a id="step-01"></a>

### 1. Build the STAR reference index

**Source:** `Agem_RNA-Seq_Dec2024/Agem/results_combined/AgemnoindexLog.out`, recorded `genomeGenerate` command  
**Save as:** `build_standard_rnaseq_star_index.sh`  
**Archive treatment:** Extract the executed command from the log; add editable paths, output-directory creation and an explicit thread setting. The index parameters are retained.

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/path/to/project"
THREADS="${THREADS:-40}"
GENOME="${ROOT}/Agem_RNA-Seq_Dec2024/Anticarsia_gemmatalis_primary_pseudohap_scaffold.fasta"
INDEX="${ROOT}/Agem_RNA-Seq_Dec2024/Agem/star_index"

mkdir -p "$INDEX"
STAR --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX" \
    --genomeFastaFiles "$GENOME" \
    --genomeSAindexNbases 12
```

<a id="step-02"></a>

### 2. Align the combined standard RNA-seq libraries

**Source:** `Agem_RNA-Seq_Dec2024/Agem/STAR_map.sh`  
**Save as:** `align_standard_rnaseq_combined.sh`  
**Archive treatment:** Remove scheduler and executable-installation settings. Use editable, quoted input/output paths and explicit threads. Retain both mate lists and all alignment settings. Remove unused variables from the original script.

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/path/to/project"
THREADS="${THREADS:-40}"
INDEX="${ROOT}/Agem_RNA-Seq_Dec2024/Agem/star_index"
READ_DIR="${ROOT}/Agem_RNA-Seq_Dec2024/Agem/fastqs"
OUTDIR="${ROOT}/Agem_RNA-Seq_Dec2024/Agem/results_combined"

READS1="AG09_R1.fastq.gz,AG10_R1.fastq.gz,AG11_R1.fastq.gz,AG12_R1.fastq.gz,AG13_R1.fastq.gz,AG16_R1.fastq.gz,AG17_R1.fastq.gz,AG18_R1.fastq.gz,AG19_R1.fastq.gz,AG20_R1.fastq.gz,Ag04_R1.fastq.gz,Ag05_R1.fastq.gz,Ag06_R1.fastq.gz"
READS2="AG09_R2.fastq.gz,AG10_R2.fastq.gz,AG11_R2.fastq.gz,AG12_R2.fastq.gz,AG13_R2.fastq.gz,AG16_R2.fastq.gz,AG17_R2.fastq.gz,AG18_R2.fastq.gz,AG19_R2.fastq.gz,AG20_R2.fastq.gz,Ag04_R2.fastq.gz,Ag05_R2.fastq.gz,Ag06_R2.fastq.gz"

mkdir -p "$OUTDIR"
cd "$READ_DIR"
STAR --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX" \
    --outFileNamePrefix "${OUTDIR}/Agemnoindex" \
    --readFilesCommand zcat \
    --readFilesIn "$READS1" "$READS2" \
    --outSAMtype BAM SortedByCoordinate \
    --alignSoftClipAtReferenceEnds No \
    --alignIntronMax 80351 \
    --alignMatesGapMax 8035 \
    --limitBAMsortRAM 23188627292
```

<a id="step-03"></a>

### 3. Filter alignments and extract the cortex/ivory region

**Save as:** `filter_extract_standard_rnaseq.sh`  
**Archive treatment:** Restore the commented filtering and filtered-BAM indexing lines. The BAM header verifies the filter and region-extraction commands. Remove scheduler/module settings, use an editable output path and expose the recorded thread count. Indexing of the unfiltered input is unnecessary for the initial whole-file filter and is omitted.

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/path/to/project"
THREADS="${THREADS:-32}"
OUTDIR="${ROOT}/Agem_RNA-Seq_Dec2024/Agem/results_combined"
REGION="Chromosome8:2410000-2850000"

cd "$OUTDIR"
samtools view -@ "$THREADS" -F 0x04 -f 0x2 -q 20 -b \
    AgemnoindexAligned.sortedByCoord.out.bam > Agem.filtered.bam
samtools index -@ "$THREADS" Agem.filtered.bam

samtools view -@ "$THREADS" -b Agem.filtered.bam "$REGION" \
    > Agem.filtered.chr8.ivory.bam
samtools index -@ "$THREADS" Agem.filtered.chr8.ivory.bam
```
