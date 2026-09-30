# HiFi haplotype assembly and locus identification

**Paper:** Figure 2D–E; Figures S10–S11; phased assembly resource  

Assemble the eight F2 females and reassemble the Dark reference individual while retaining haplotype-specific outputs. Locate cortex-containing contigs for downstream comparison.

## Inputs

- Pacbio_Reseq_Hilo/Agem_{Dark,Light}_{01,02,05,06}/*.bam (original multiplexed HiFi reads)
- Pacbio_Reseq_Hilo/Reference_genome_reassembly/sra_cache/SRR32404155/SRR32404155.sra
- Pacbio_Reseq_Hilo/cortex.fa

## Outputs

- Filtered HiFi FASTQs; hifiasm hap1/hap2 contig GFAs and unitig graphs.
- Haplotype FASTAs and cortex BLAST tables; section 5a supplies a new GFA-to-FASTA export using the documented hifiasm command. Locus FASTAs for Figure 2D/S10 were prepared by extracting and orienting the regions in Geneious.

## Software

Bash, AWK, HiFiAdapterFilt, BAMtools, pigz, BLAST+, SRA Toolkit, hifiasm. Successful assembly logs report hifiasm 0.25.0-r726; the source scripts specify BLAST+ 2.17.0+.

## Execution and interpretation

1. Run the Dark02 adapter-filtering step before its separate assembly step. Run the seven-sample script once for each index from 0 through 6.
2. The SRA export requires the SRR32404155 archive to have been downloaded into the stated cache. The reference reassembly command is included below.
3. The reassembly of the reference individual supplies phased locus comparisons. It is distinct from the original chromosome-level ilAntGemm2 assembly.
4. Before cortex BLAST, run section 5a to export the four FASTAs expected per individual (p_utg, r_utg, hap1 and hap2) from the sequence-bearing hifiasm GFA files.
5. Identifiers 01/02/05/06 correspond to publication individuals 1/2/3/4, respectively, for both Dark and Light morphs. The scripts retain the original input-file identifiers

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Filter Dark02 HiFi adapters

Run once for Dark02; input BAM and output prefix are explicit.

**Source:** `Pacbio_Reseq_Hilo/run_hfaf_1t_Agem_Dark_02_gpu.sh`  
**Save as:** `run_hfaf_1t_Agem_Dark_02_gpu.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
set -euo pipefail


BASE="/path/to/project/Pacbio_Reseq_Hilo"
SAMPLE="Agem_Dark_02"
SAMPLEDIR="${BASE}/${SAMPLE}"
BAM="m84125_260225_034835_s4.hifi_reads.bc2042.bam"
PREFIX="${BAM%.bam}"

OUTDIR="${SAMPLEDIR}/hfaf_test_gpu_t1"
THREADS=1

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

ln -sf "${SAMPLEDIR}/${BAM}" "${BAM}"

echo "Running HiFiAdapterFilt with THREADS=${THREADS}"
hifiadapterfilt.sh -p "${PREFIX}" -t "${THREADS}" -o "${OUTDIR}"

echo
echo "Done. Final files:"
ls -lh
```

<a id="step-02"></a>

### 2. Assemble Dark02

Consumes the filtered FASTQ from the preceding step.

**Source:** `Pacbio_Reseq_Hilo/run_hifiasm_Agem_Dark_02_from_hfaf.sh`  
**Save as:** `run_hifiasm_Agem_Dark_02_from_hfaf.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE="/path/to/project/Pacbio_Reseq_Hilo"
SAMPLE="Agem_Dark_02"

READS="${BASE}/${SAMPLE}/hfaf_test_gpu_t1/m84125_260225_034835_s4.hifi_reads.bc2042.filt.fastq.gz"
HIFIASM="$(command -v hifiasm)"

OUTDIR="${BASE}/${SAMPLE}/hifiasm_from_hfaf_t1"
PREFIX="${OUTDIR}/${SAMPLE}"

THREADS=${THREADS:-40}

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "======================================="
echo "Job started on: $(date)"
echo "PWD: $(pwd)"
echo "Sample: ${SAMPLE}"
echo "Reads: ${READS}"
echo "Threads: ${THREADS}"
echo "======================================="

if [[ ! -s "${READS}" ]]; then
    echo "ERROR: reads file not found: ${READS}" >&2
    exit 1
fi

if [[ ! -x "${HIFIASM}" ]]; then
    echo "ERROR: hifiasm binary not executable: ${HIFIASM}" >&2
    exit 1
fi

echo
echo "Running hifiasm..."
echo "${HIFIASM} -t ${THREADS} -o ${PREFIX} ${READS}"

"${HIFIASM}" \
    -t "${THREADS}" \
    -o "${PREFIX}" \
    "${READS}"

echo
echo "Assembly directory contents:"
ls -lh

echo
echo "Key GFA outputs:"
ls -lh "${PREFIX}"*.gfa || true

echo
echo "Job finished on: $(date)"
```

<a id="step-03"></a>

### 3. Filter and assemble the other seven females

Run with sample indices 0 through 6; each invocation filters and assembles one sample.

**Source:** `Pacbio_Reseq_Hilo/run_hifiasm_array_remaining_gpu.sh`  
**Save as:** `run_hifiasm_array_remaining_gpu.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Replace the job-array index with a required, range-checked sample-index argument. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-40}"
# Run once per sample: bash run_hifiasm_array_remaining_gpu.sh INDEX (INDEX 0 through 6).
TASK_ID="${1:?Supply a zero-based sample index (0-6)}"
if [[ ! "$TASK_ID" =~ ^[0-6]$ ]]; then
  echo "Sample index must be 0-6." >&2
  exit 2
fi
set -euo pipefail


BASE="/path/to/project/Pacbio_Reseq_Hilo"
HIFIASM="$(command -v hifiasm)"

SAMPLES=(
  "Agem_Dark_01"
  "Agem_Dark_05"
  "Agem_Dark_06"
  "Agem_Light_01"
  "Agem_Light_02"
  "Agem_Light_05"
  "Agem_Light_06"
)

SAMPLE="${SAMPLES[$TASK_ID]}"
SAMPLEDIR="${BASE}/${SAMPLE}"

THREADS_HIFIASM=${THREADS:-40}
THREADS_HFAF=1

echo "======================================="
echo "Job started on: $(date)"
echo "Sample: ${SAMPLE}"
echo "Sample dir: ${SAMPLEDIR}"
echo "Array task ID: ${TASK_ID}"
echo "======================================="

if [[ ! -d "${SAMPLEDIR}" ]]; then
    echo "ERROR: sample directory not found: ${SAMPLEDIR}" >&2
    exit 1
fi

BAM=$(find "${SAMPLEDIR}" -maxdepth 1 -type f -name "*.bam" | head -n 1 || true)

if [[ -z "${BAM}" || ! -s "${BAM}" ]]; then
    echo "ERROR: no BAM found in ${SAMPLEDIR}" >&2
    exit 1
fi

BAM_BASENAME=$(basename "${BAM}")
PREFIX="${BAM_BASENAME%.bam}"

if [[ ! -x "${HIFIASM}" ]]; then
    echo "ERROR: hifiasm binary not executable: ${HIFIASM}" >&2
    exit 1
fi

########################################
# Step 1: HiFiAdapterFilt
########################################
HFAF_DIR="${SAMPLEDIR}/hfaf_gpu_t1"
mkdir -p "${HFAF_DIR}"
cd "${HFAF_DIR}"

ln -sf "${BAM}" "${BAM_BASENAME}"

echo
echo "Step 1: running HiFiAdapterFilt"
echo "BAM: ${BAM_BASENAME}"
echo "Prefix: ${PREFIX}"
echo "HiFiAdapterFilt threads: ${THREADS_HFAF}"

hifiadapterfilt.sh -p "${PREFIX}" -t "${THREADS_HFAF}" -o "${HFAF_DIR}"

FILTERED_FASTQ="${HFAF_DIR}/${PREFIX}.filt.fastq.gz"
STATS_FILE="${HFAF_DIR}/${PREFIX}.stats"
BLOCKLIST_FILE="${HFAF_DIR}/${PREFIX}.blocklist"

if [[ ! -s "${FILTERED_FASTQ}" ]]; then
    echo "ERROR: filtered FASTQ not produced: ${FILTERED_FASTQ}" >&2
    exit 1
fi

echo
echo "HiFiAdapterFilt outputs:"
ls -lh "${FILTERED_FASTQ}" "${STATS_FILE}" "${BLOCKLIST_FILE}" || true

########################################
# Step 2: hifiasm
########################################
ASM_DIR="${SAMPLEDIR}/hifiasm_from_hfaf_t1"
mkdir -p "${ASM_DIR}"
cd "${ASM_DIR}"

OUT_PREFIX="${ASM_DIR}/${SAMPLE}"

echo
echo "Step 2: running hifiasm"
echo "Reads: ${FILTERED_FASTQ}"
echo "hifiasm threads: ${THREADS_HIFIASM}"
echo "Command: ${HIFIASM} -t ${THREADS_HIFIASM} -o ${OUT_PREFIX} ${FILTERED_FASTQ}"

"${HIFIASM}" \
    -t "${THREADS_HIFIASM}" \
    -o "${OUT_PREFIX}" \
    "${FILTERED_FASTQ}"

echo
echo "Final assembly directory contents:"
ls -lh

echo
echo "Key GFA outputs:"
ls -lh "${OUT_PREFIX}"*.gfa || true

echo
echo "Job finished on: $(date)"
echo "======================================="
```

<a id="step-04"></a>

### 4. Export reference HiFi reads

Exports the existing SRA cache to SRR32404155.fastq.

**Source:** `Pacbio_Reseq_Hilo/Reference_genome_reassembly/run_fasterqdump_hifi_only.sh`  
**Save as:** `run_fasterqdump_hifi_only.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE="/path/to/project/Pacbio_Reseq_Hilo/Reference_genome_reassembly"
THREADS=${THREADS:-8}

SRA="${BASE}/sra_cache/SRR32404155/SRR32404155.sra"
OUTDIR="${BASE}/HiFi_reads"

mkdir -p "${OUTDIR}"

echo "======================================="
echo "Job started on: $(date)"
echo "SRA: ${SRA}"
echo "OUTDIR: ${OUTDIR}"
echo "Threads: ${THREADS}"
echo "======================================="

if [[ ! -s "${SRA}" ]]; then
    echo "ERROR: SRA file not found: ${SRA}" >&2
    exit 1
fi

fasterq-dump \
  --threads "${THREADS}" \
  --outdir "${OUTDIR}" \
  "${SRA}"

echo
echo "Final contents:"
ls -lh "${OUTDIR}"

echo
echo "Line counts:"
wc -l "${OUTDIR}"/* || true

echo
echo "Job finished on: $(date)"
```

<a id="step-05"></a>

### 5. Reassemble the Dark reference individual

Reassemble the reference individual from its HiFi reads.

**Source:** `Pacbio_Reseq_Hilo/Reference_genome_reassembly/hifiasm_ref.sh`  
**Save as:** `hifiasm_ref.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE="/path/to/project/Pacbio_Reseq_Hilo/Reference_genome_reassembly"
READS="${BASE}/HiFi_reads/SRR32404155.fastq"
HIFIASM="$(command -v hifiasm)"

OUTDIR="${BASE}/hifiasm_reference"
PREFIX="${OUTDIR}/Reference_reassembly"

THREADS=${THREADS:-40}

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "======================================="
echo "Job started on: $(date)"
echo "PWD: $(pwd)"
echo "Reads: ${READS}"
echo "Threads: ${THREADS}"
echo "======================================="

if [[ ! -s "${READS}" ]]; then
    echo "ERROR: reads file not found: ${READS}" >&2
    exit 1
fi

if [[ ! -x "${HIFIASM}" ]]; then
    echo "ERROR: hifiasm binary not executable: ${HIFIASM}" >&2
    exit 1
fi

echo
echo "Running hifiasm..."
echo "${HIFIASM} -t ${THREADS} -o ${PREFIX} ${READS}"

"${HIFIASM}" \
    -t "${THREADS}" \
    -o "${PREFIX}" \
    "${READS}"

echo
echo "Assembly directory contents:"
ls -lh

echo
echo "Key GFA outputs:"
ls -lh "${PREFIX}"*.gfa || true

echo
echo "Job finished on: $(date)"
```

<a id="step-05a"></a>

### 5a. Export hifiasm GFA sequences to FASTA

**Save as:** `export_hifiasm_fastas.sh`  

Run `bash export_hifiasm_fastas.sh /path/to/assembly/SAMPLE /path/to/exports`, where `SAMPLE` is the output prefix passed to hifiasm. The script reads the sequence-bearing `.bp.hap1.p_ctg.gfa`, `.bp.hap2.p_ctg.gfa`, `.bp.p_utg.gfa` and `.bp.r_utg.gfa` files. It creates the matching `haplotype_fastas` and `untig_fastas` subdirectories. Use those directories for the BLAST inputs in section 6, adjusting the example paths as needed. Existing FASTAs are protected from overwriting.

The conversion exports GFA segment names and sequences in their original order. It does not rerun assembly or orient the locus sequences. Do not use `.noseq.gfa` files, which lack sequence. For the checked Dark01 assembly, all four exported FASTAs match the retained files byte for byte.

```bash
#!/usr/bin/env bash
set -euo pipefail
set -C  # Refuse to overwrite existing FASTAs.
# Usage: bash export_hifiasm_fastas.sh /path/to/SAMPLE /path/to/exports
PREFIX=${1:?Supply the hifiasm output prefix without .bp or .gfa}
OUTDIR=${2:?Supply an output directory for the exported FASTAs}
SAMPLE=$(basename "$PREFIX")
for kind in hap1.p_ctg hap2.p_ctg p_utg r_utg; do
    [[ -s "$PREFIX.bp.$kind.gfa" ]] || { echo "Missing $PREFIX.bp.$kind.gfa" >&2; exit 1; }
done
mkdir -p "$OUTDIR/haplotype_fastas" "$OUTDIR/untig_fastas"
for kind in hap1.p_ctg hap2.p_ctg p_utg r_utg; do
    case "$kind" in
        hap*) subdir=haplotype_fastas ;;
        *) subdir=untig_fastas ;;
    esac
    # Export sequence-bearing GFA segment records; retain their identifiers/order.
    awk -F '\t' '
        $1 == "S" {
            if (NF < 3 || $3 == "*") {
                print "Use a sequence-bearing GFA, not a .noseq.gfa" > "/dev/stderr"
                exit 1
            }
            print ">" $2
            print $3
            count++
        }
        END { if (!count) exit 1 }
    ' "$PREFIX.bp.$kind.gfa" > "$OUTDIR/$subdir/$SAMPLE.bp.$kind.fa"
done
```

<a id="step-06"></a>

### 6. Identify cortex-containing haplotype contigs

Uses the GFA-to-FASTA exports from section 5a. BLAST tables identify candidate contigs. Locus extraction and orientation for Figure 2D and Figure S10 were performed in Geneious, followed by FASTA export.

**Source:** `Pacbio_Reseq_Hilo/blast_cort_all.sh`  
**Save as:** `blast_cort_all.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
set -euo pipefail


BASE="/path/to/project/Pacbio_Reseq_Hilo"
QUERY="${BASE}/cortex.fa"
THREADS=${THREADS:-8}

SAMPLES=(
  Agem_Dark_01
  Agem_Dark_02
  Agem_Dark_05
  Agem_Dark_06
  Agem_Light_01
  Agem_Light_02
  Agem_Light_05
  Agem_Light_06
)

echo "======================================="
echo "Job started on: $(date)"
echo "Base: ${BASE}"
echo "Query: ${QUERY}"
echo "Threads: ${THREADS}"
echo "======================================="

if [[ ! -s "${QUERY}" ]]; then
    echo "ERROR: query file not found: ${QUERY}" >&2
    exit 1
fi

for SAMPLE in "${SAMPLES[@]}"; do
    ASM_DIR="${BASE}/${SAMPLE}/hifiasm_from_hfaf_t1"
    UNTIG_DIR="${ASM_DIR}/untig_fastas"
    HAP_DIR="${ASM_DIR}/haplotype_fastas"
    OUTDIR="${ASM_DIR}/cortex_blasts"

    P_UTG="${UNTIG_DIR}/${SAMPLE}.bp.p_utg.fa"
    R_UTG="${UNTIG_DIR}/${SAMPLE}.bp.r_utg.fa"
    HAP1="${HAP_DIR}/${SAMPLE}.bp.hap1.p_ctg.fa"
    HAP2="${HAP_DIR}/${SAMPLE}.bp.hap2.p_ctg.fa"

    mkdir -p "${OUTDIR}"

    echo
    echo "Processing ${SAMPLE} ..."
    echo "ASM_DIR=${ASM_DIR}"

    for fasta in "${P_UTG}" "${R_UTG}" "${HAP1}" "${HAP2}"; do
        if [[ ! -s "${fasta}" ]]; then
            echo "ERROR: missing FASTA: ${fasta}" >&2
            exit 1
        fi
    done

    makeblastdb -in "${P_UTG}" -dbtype nucl >/dev/null
    makeblastdb -in "${R_UTG}" -dbtype nucl >/dev/null
    makeblastdb -in "${HAP1}" -dbtype nucl >/dev/null
    makeblastdb -in "${HAP2}" -dbtype nucl >/dev/null

    blastn \
      -query "${QUERY}" \
      -db "${P_UTG}" \
      -num_threads "${THREADS}" \
      -out "${OUTDIR}/cortex_vs_${SAMPLE}.bp.p_utg.tsv" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
      -max_target_seqs 20

    blastn \
      -query "${QUERY}" \
      -db "${R_UTG}" \
      -num_threads "${THREADS}" \
      -out "${OUTDIR}/cortex_vs_${SAMPLE}.bp.r_utg.tsv" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
      -max_target_seqs 20

    blastn \
      -query "${QUERY}" \
      -db "${HAP1}" \
      -num_threads "${THREADS}" \
      -out "${OUTDIR}/cortex_vs_${SAMPLE}.bp.hap1.p_ctg.tsv" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
      -max_target_seqs 20

    blastn \
      -query "${QUERY}" \
      -db "${HAP2}" \
      -num_threads "${THREADS}" \
      -out "${OUTDIR}/cortex_vs_${SAMPLE}.bp.hap2.p_ctg.tsv" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
      -max_target_seqs 20

    echo "Top hits for ${SAMPLE}:"
    echo "=== p_utg ==="
    head -5 "${OUTDIR}/cortex_vs_${SAMPLE}.bp.p_utg.tsv" || true
    echo "=== r_utg ==="
    head -5 "${OUTDIR}/cortex_vs_${SAMPLE}.bp.r_utg.tsv" || true
    echo "=== hap1 ==="
    head -5 "${OUTDIR}/cortex_vs_${SAMPLE}.bp.hap1.p_ctg.tsv" || true
    echo "=== hap2 ==="
    head -5 "${OUTDIR}/cortex_vs_${SAMPLE}.bp.hap2.p_ctg.tsv" || true
done

echo
echo "Job finished on: $(date)"
echo "======================================="
```
