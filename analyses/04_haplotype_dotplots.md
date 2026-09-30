# Haplotype sequence dotplots


**Paper:** Figure 2D, Figure 3B–C and Figure S10  

Compare homologous locus sequences using BLASTN and visualize alignment identity and structural rearrangements.

## Inputs

- Pacbio_Reseq_Hilo/Dotplots/Reference_GWA_region.fasta
- Pacbio_Reseq_Hilo/Dotplots/*.fasta (prepared haplotype intervals)
- Pacbio_Reseq_Hilo/Dotplots/Zoomed_for_figure_inv_haplos/* (the three FASTAs named in the second script)

For Figure 2D and Figure S10, the locus sequences were extracted from the full assemblies, oriented as needed in Geneious, and exported as FASTA files. The plotting workflows start from these prepared FASTAs.

## Outputs

- Figure 2D broad-region BLAST alignments and haplotype dotplots in blastn2dotplots_fig2d_500bp/.
- S10-like 400-kb dotplots in blastn2dotplot_zoomed_pdf/.
- Absolute-coordinate inversion-region dotplots for Dark01-hap2/B and Light06/C.

## Software

Bash, Python 3, BLAST+ and blastn2dotplots.

## Execution and interpretation

1. All three plotting commands use minimum identity 70% and minimum alignment length 500 bp.
2. The S10 workflow extracts local bases 350,000–750,000 from interval FASTAs prepared in Geneious. This additional local crop is performed by the script below; the preceding locus extraction and orientation were performed in Geneious.
3. The inversion comparison uses an absolute-coordinate script.
4. Figure 2D uses the full prepared locus intervals, without the S10 crop. Section 3 adapts the recovered broad-region script to a 500-bp minimum alignment length.
5. The code below generates the alignment dotplots in Figure 3B–C.

The three sections are independent workflows for their respective panels. Save the relevant code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the selected script. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Generate the S10 region dotplots

Runs BLASTN and blastn2dotplots on the same local window for each prepared haplotype.

**Source:** `Pacbio_Reseq_Hilo/Dotplots/blastn2dotplots_zoomed.sh`  
**Save as:** `blastn2dotplots_zoomed.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-8}"
set -euo pipefail


WORKDIR="/path/to/project/Pacbio_Reseq_Hilo/Dotplots"
OUTROOT="${WORKDIR}/blastn2dotplos_zoomed"
PDFDIR="${WORKDIR}/blastn2dotplot_zoomed_pdf"

REF="Reference_GWA_region.fasta"

# Square zoom window
WIN_START=350000
WIN_END=750000

mkdir -p "${OUTROOT}" "${PDFDIR}"

cd "${WORKDIR}"

echo "Working directory: ${WORKDIR}"
echo "Reference: ${REF}"
echo "Window: ${WIN_START}-${WIN_END}"
echo "Per-sample output root: ${OUTROOT}"
echo "Collected PDF directory: ${PDFDIR}"

echo
echo "### Reference header and length"
grep "^>" "${REF}"
echo -n "${REF}  "
grep -v "^>" "${REF}" | tr -d '\n' | wc -c

# ------------------------------------------------------------------
# Make the reference window once
# ------------------------------------------------------------------
REF_WINDOW_FA="${OUTROOT}/Reference_GWA_region_${WIN_START}_${WIN_END}.fasta"

python3 - <<PY
from pathlib import Path

ref_fa = Path("${REF}")
out_fa = Path("${REF_WINDOW_FA}")
start = ${WIN_START}
end = ${WIN_END}

header = None
seq_parts = []

with open(ref_fa) as fh:
    for line in fh:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith(">"):
            if header is None:
                header = line[1:].split()[0]
            else:
                raise SystemExit("ERROR: Reference fasta contains more than one sequence; expected one.")
        else:
            seq_parts.append(line)

seq = "".join(seq_parts)
seqlen = len(seq)

if start < 1 or end > seqlen or start > end:
    raise SystemExit(f"ERROR: Requested reference window {start}-{end} is invalid for reference length {seqlen}")

subseq = seq[start-1:end]

with open(out_fa, "w") as out:
    out.write(f">Reference_GWA_region_{start}_{end}\n")
    for i in range(0, len(subseq), 80):
        out.write(subseq[i:i+80] + "\n")

print(f"Reference window written to: {out_fa}")
print(f"Reference window length: {len(subseq)} bp")
PY

echo
echo "### Reference window header and length"
grep "^>" "${REF_WINDOW_FA}"
grep -v "^>" "${REF_WINDOW_FA}" | tr -d '\n' | wc -c

# ------------------------------------------------------------------
# Process all sample FASTAs except the reference window fasta itself
# ------------------------------------------------------------------
shopt -s nullglob

for HAP in *.fasta; do
    if [[ "${HAP}" == "${REF}" ]]; then
        continue
    fi

    SAMPLE="${HAP%.fasta}"
    INDDIR="${OUTROOT}/${SAMPLE}"
    DBDIR="${INDDIR}/blast_db"
    PLOTDIR="${INDDIR}/plots"
    TMPDIR="${INDDIR}/tmp"

    mkdir -p "${DBDIR}" "${PLOTDIR}" "${TMPDIR}"

    echo
    echo "============================================================"
    echo "Processing sample: ${SAMPLE}"
    echo "FASTA: ${HAP}"
    echo "Output directory: ${INDDIR}"
    echo "============================================================"

    echo "Header:"
    grep "^>" "${HAP}"

    HAP_LEN=$(grep -v "^>" "${HAP}" | tr -d '\n' | wc -c)
    echo "Length: ${HAP_LEN}"

    if (( HAP_LEN < WIN_END )); then
        echo "Skipping ${SAMPLE}: length ${HAP_LEN} is shorter than WIN_END ${WIN_END}"
        continue
    fi

    HAP_WINDOW_FA="${TMPDIR}/${SAMPLE}_${WIN_START}_${WIN_END}.fasta"

    python3 - <<PY
from pathlib import Path

in_fa = Path("${HAP}")
out_fa = Path("${HAP_WINDOW_FA}")
start = ${WIN_START}
end = ${WIN_END}
out_header = "${SAMPLE}_${WIN_START}_${WIN_END}"

header = None
seq_parts = []

with open(in_fa) as fh:
    for line in fh:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith(">"):
            if header is None:
                header = line[1:].split()[0]
            else:
                raise SystemExit(f"ERROR: {in_fa} contains more than one sequence; expected one.")
        else:
            seq_parts.append(line)

seq = "".join(seq_parts)
seqlen = len(seq)

if start < 1 or end > seqlen or start > end:
    raise SystemExit(f"ERROR: Requested window {start}-{end} is invalid for {in_fa} length {seqlen}")

subseq = seq[start-1:end]

with open(out_fa, "w") as out:
    out.write(f">{out_header}\n")
    for i in range(0, len(subseq), 80):
        out.write(subseq[i:i+80] + "\n")

print(f"Wrote {out_fa} ({len(subseq)} bp)")
PY

    echo "Extracted haplotype window:"
    grep "^>" "${HAP_WINDOW_FA}"
    grep -v "^>" "${HAP_WINDOW_FA}" | tr -d '\n' | wc -c

    DBPREFIX="${DBDIR}/${SAMPLE}_${WIN_START}_${WIN_END}"

    makeblastdb \
        -in "${HAP_WINDOW_FA}" \
        -dbtype nucl \
        -parse_seqids \
        -out "${DBPREFIX}"

    BLAST_OUT="${INDDIR}/${SAMPLE}_vs_ReferenceWindow_${WIN_START}_${WIN_END}.blastn.tsv"

    blastn \
        -query "${REF_WINDOW_FA}" \
        -db "${DBPREFIX}" \
        -out "${BLAST_OUT}" \
        -outfmt "6 std qlen slen" \
        -num_threads "${THREADS}"

    echo "BLAST line count:"
    wc -l "${BLAST_OUT}"

    DB_TXT="${INDDIR}/db.txt"
    QUERY_TXT="${INDDIR}/query.txt"

    cat > "${DB_TXT}" <<EOF
${SAMPLE}_${WIN_START}_${WIN_END}	${SAMPLE}	${WIN_START}
EOF

    cat > "${QUERY_TXT}" <<EOF
Reference_GWA_region_${WIN_START}_${WIN_END}	Reference_GWA_region	${WIN_START}
EOF

    OUTPREFIX="${PLOTDIR}/${SAMPLE}_vs_ReferenceWindow_${WIN_START}_${WIN_END}"

    blastn2dotplots \
        -i1 "${DB_TXT}" \
        -i2 "${QUERY_TXT}" \
        --blastn "${BLAST_OUT}" \
        --out "${OUTPREFIX}" \
        --min_identity 70 \
        --min_alignlen 500 \
        --line_width 14.0 \
        --figure_size 6 6 \
        --font_size 10 \
        --tick_label_size 8 \
        --colormap 4

    PDF="${OUTPREFIX}.pdf"

    if [[ -f "${PDF}" ]]; then
        cp -f "${PDF}" "${PDFDIR}/"
        echo "Copied PDF to ${PDFDIR}/$(basename "${PDF}")"
    else
        echo "WARNING: PDF not found for ${SAMPLE}"
    fi
done

echo
echo "### Final collected PDFs"
ls -lh "${PDFDIR}"

echo
echo "Done."
```

<a id="step-02"></a>

### 2. Generate inversion-region dotplots in absolute coordinates

Draws the alignment dotplots using the Figure 3B–C coordinate ranges.

**Source:** `Pacbio_Reseq_Hilo/Dotplots/Zoomed_for_figure_inv_haplos/run_blastn2dotplots_zoomed_inv_absXY.sh`  
**Save as:** `run_blastn2dotplots_zoomed_inv_absXY.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults.

```bash
#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/path/to/project/Pacbio_Reseq_Hilo/Dotplots/Zoomed_for_figure_inv_haplos"
cd "${WORKDIR}"

mkdir -p 00_logs 01_seqtxt 02_blastn 03_dotplots

THREADS="${THREADS:-40}"

REF_FASTA="Dark_Reference_2.4to2.8Mb.fasta"
REF_ID="Dark_Reference_2.4to2.8Mb"
REF_START=2400000
REF_END=2800000

# sample_fasta|sample_id|sample_start|sample_end|label
PAIRS=(
  "Dark_01_INV_equivalent_to_2.4to2.8Mb.fasta|Dark_01_INV_equivalent_to_2.4to2.8Mb|2309343|2674611|DarkRef_vs_Dark01_INV"
  "Light_06_INVT_equivalent_to_2.4to2.8Mb.fasta|Light_06_INVT_equivalent_to_2.4to2.8Mb|2338747|2677702|DarkRef_vs_Light06_INVT"
)

check_cmd() {
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: ${cmd} not found in PATH on compute node." >&2
    exit 1
  fi
}

fasta_len() {
  awk '
    /^>/ {next}
    {
      gsub(/[[:space:]]/, "", $0)
      n += length($0)
    }
    END {print n}
  ' "$1"
}

echo "=== Preflight ==="
echo "Working directory: ${WORKDIR}"
echo "Threads: ${THREADS}"
date

check_cmd blastn
check_cmd makeblastdb
check_cmd blastn2dotplots

echo "blastn         : $(command -v blastn)"
echo "makeblastdb    : $(command -v makeblastdb)"
echo "blastn2dotplots: $(command -v blastn2dotplots)"
echo

if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR: Missing reference FASTA: ${REF_FASTA}" >&2
  exit 1
fi

REF_LEN=$(fasta_len "${REF_FASTA}")
REF_EXPECTED_LEN=$((REF_END - REF_START + 1))

echo "Reference FASTA        : ${REF_FASTA}"
echo "Reference ID           : ${REF_ID}"
echo "Reference abs interval : ${REF_START}-${REF_END}"
echo "Reference FASTA length : ${REF_LEN}"
echo "Reference expected len : ${REF_EXPECTED_LEN}"

if [[ "${REF_LEN}" -ne "${REF_EXPECTED_LEN}" ]]; then
  echo "ERROR: Reference FASTA length does not match supplied coordinates." >&2
  exit 1
fi
echo

for entry in "${PAIRS[@]}"; do
  IFS='|' read -r SAMPLE_FASTA SAMPLE_ID SAMPLE_START SAMPLE_END LABEL <<< "${entry}"

  if [[ ! -s "${SAMPLE_FASTA}" ]]; then
    echo "ERROR: Missing sample FASTA: ${SAMPLE_FASTA}" >&2
    exit 1
  fi

  SAMPLE_LEN=$(fasta_len "${SAMPLE_FASTA}")
  SAMPLE_EXPECTED_LEN=$((SAMPLE_END - SAMPLE_START + 1))

  echo "=== Processing ${LABEL} ==="
  echo "Sample FASTA        : ${SAMPLE_FASTA}"
  echo "Sample ID           : ${SAMPLE_ID}"
  echo "Sample abs interval : ${SAMPLE_START}-${SAMPLE_END}"
  echo "Sample FASTA length : ${SAMPLE_LEN}"
  echo "Sample expected len : ${SAMPLE_EXPECTED_LEN}"

  if [[ "${SAMPLE_LEN}" -ne "${SAMPLE_EXPECTED_LEN}" ]]; then
    echo "ERROR: Sample FASTA length does not match supplied coordinates for ${SAMPLE_ID}." >&2
    exit 1
  fi

  DB_PREFIX="02_blastn/${LABEL}"
  BLAST_OUT="02_blastn/${LABEL}.blastn.tsv"
  OUTPREFIX="03_dotplots/${LABEL}"

  printf "%s\t%s\t%s\n" \
    "${SAMPLE_ID}" "${SAMPLE_ID}" "${SAMPLE_START}" \
    > "01_seqtxt/${LABEL}.input1_db_rows.tsv"

  printf "%s\t%s\t%s\n" \
    "${REF_ID}" "${REF_ID}" "${REF_START}" \
    > "01_seqtxt/${LABEL}.input2_query_cols.tsv"

  echo "Making BLAST database for ${SAMPLE_FASTA}..."
  makeblastdb \
    -in "${SAMPLE_FASTA}" \
    -dbtype nucl \
    -out "${DB_PREFIX}" \
    > "02_blastn/${LABEL}.makeblastdb.log" 2>&1

  echo "Running blastn..."
  blastn \
    -query "${REF_FASTA}" \
    -db "${DB_PREFIX}" \
    -outfmt "6 std qlen slen" \
    -num_threads "${THREADS}" \
    -out "${BLAST_OUT}"

  echo "BLAST alignments: $(wc -l < "${BLAST_OUT}")"
  head -5 "${BLAST_OUT}" || true

  echo "Running blastn2dotplots..."
  blastn2dotplots \
    -i1 "01_seqtxt/${LABEL}.input1_db_rows.tsv" \
    -i2 "01_seqtxt/${LABEL}.input2_query_cols.tsv" \
    --blastn "${BLAST_OUT}" \
    --out "${OUTPREFIX}" \
    --min_identity 70 \
    --min_alignlen 500 \
    --line_width 7.0 \
    --figure_size 6 6 \
    --font_size 10 \
    --tick_label_size 8 \
    --colormap 4

  echo "Finished ${LABEL}"
  ls -lh "${OUTPREFIX}.pdf"
  echo
done

echo "=== All done ==="
date
ls -lh 03_dotplots/*.pdf
```

<a id="step-03"></a>

### 3. Generate the Figure 2D broad-region dotplots

Align both prepared haplotype intervals from each of the eight individuals against `Reference_GWA_region.fasta`, then plot alignments using 70% minimum identity and a 500-bp minimum alignment length. This workflow uses the full supplied intervals; it does not extract the narrower S10 window. Each FASTA must contain one prepared, oriented locus sequence. The recovered file identifiers are retained. The 01/02/05/06 samples correspond to paper individuals 1/2/3/4 for both morphs.

Set `WORKDIR` to the folder containing the reference FASTA and the 16 haplotype FASTAs listed below. Save this block as `run_blastn2dotplots_fig2d_500bp.sh` and run it with Bash. Outputs include BLAST tables and one PDF per haplotype, grouped by individual.

**Source:** `Pacbio_Reseq_Hilo/Dotplots/run_blastn2dotplots_all.sh`  
**Save as:** `run_blastn2dotplots_fig2d_500bp.sh`  
**Archive treatment:** set --min_alignlen to 500 bp; retain 70% minimum identity and the recovered BLAST/plotting options. Limit the sample list to the eight Figure 2D individuals, omitting the additional reference-reassembly comparison. Replace local paths and scheduler/environment setup with an editable working directory, PATH tools and a THREADS default of 8. Write outputs under blastn2dotplots_fig2d_500bp.

```bash
#!/usr/bin/env bash

set -euo pipefail

########################################
# USER SETTINGS
########################################
WORKDIR="/path/to/project/Pacbio_Reseq_Hilo/Dotplots"
OUTROOT="${WORKDIR}/blastn2dotplots_fig2d_500bp"
REF="Reference_GWA_region.fasta"

# If blastn2dotplots is not in PATH, replace this with the full path
B2D="blastn2dotplots"

THREADS="${THREADS:-8}"

########################################
# WORKING DIRECTORY
########################################
cd "${WORKDIR}"

########################################
# BASIC CHECKS
########################################
echo "Working directory: ${WORKDIR}"
echo "Output root: ${OUTROOT}"
echo "Reference: ${REF}"
echo "Threads: ${THREADS}"

command -v makeblastdb >/dev/null 2>&1 || { echo "ERROR: makeblastdb not found in PATH"; exit 1; }
command -v blastn >/dev/null 2>&1 || { echo "ERROR: blastn not found in PATH"; exit 1; }
command -v "${B2D}" >/dev/null 2>&1 || { echo "ERROR: ${B2D} not found in PATH"; exit 1; }

[[ -f "${REF}" ]] || { echo "ERROR: Reference FASTA ${REF} not found"; exit 1; }

mkdir -p "${OUTROOT}"

########################################
# FIGURE 2D: EIGHT INDIVIDUALS, TWO HAPLOTYPES EACH
########################################
declare -A SAMPLE_HAPS
SAMPLE_HAPS["Dark_01"]="Dark_01_Hap1_h1tg000014l.fasta Dark_01_Hap2_h2tg000012l.fasta"
SAMPLE_HAPS["Dark_02"]="Dark_02_Hap1_h1tg000022l.fasta Dark_02_Hap2_h2tg000023l.fasta"
SAMPLE_HAPS["Dark_05"]="Dark_05_Hap1_h1tg000028l.fasta Dark_05_Hap2_h2tg000028l.fasta"
SAMPLE_HAPS["Dark_06"]="Dark_06_Hap1_h1tg000022l.fasta Dark_06_Hap2_h2tg000022l.fasta"
SAMPLE_HAPS["Light_01"]="Light_01_Hap1_h1tg000024l.fasta Light_01_Hap2_h2tg000026l.fasta"
SAMPLE_HAPS["Light_02"]="Light_02_Hap1_h1tg000029l.fasta Light_02_Hap2_h2tg000023l.fasta"
SAMPLE_HAPS["Light_05"]="Light_05_Hap1_h1tg000007l.fasta Light_05_Hap2_h2tg000007l.fasta"
SAMPLE_HAPS["Light_06"]="Light_06_Hap1_h1tg000007l.fasta Light_06_Hap2_h2tg000007l.fasta"

########################################
# HELPER
########################################
get_fasta_header() {
    local fasta="$1"
    awk '/^>/{sub(/^>/,"",$1); print $1; exit}' "${fasta}"
}

########################################
# MAIN LOOP
########################################
for sample in "${!SAMPLE_HAPS[@]}"; do
    echo
    echo "=================================================="
    echo "Processing sample: ${sample}"
    echo "=================================================="

    sample_dir="${OUTROOT}/${sample}"
    mkdir -p "${sample_dir}"

    cp -f "${REF}" "${sample_dir}/"
    ref_copy="${sample_dir}/$(basename "${REF}")"
    ref_id="$(get_fasta_header "${ref_copy}")"

    echo "Reference header: ${ref_id}"

    makeblastdb \
        -in "${ref_copy}" \
        -dbtype nucl \
        -out "${sample_dir}/refdb"

    printf "%s\n" "${ref_id}" > "${sample_dir}/db.txt"

    for hap in ${SAMPLE_HAPS[$sample]}; do
        echo
        echo "  --- Haplotype: ${hap}"

        [[ -f "${hap}" ]] || { echo "ERROR: Missing haplotype FASTA ${hap}"; exit 1; }

        hap_base="${hap%.fasta}"
        hap_id="$(get_fasta_header "${hap}")"

        echo "  Query header: ${hap_id}"

        cp -f "${hap}" "${sample_dir}/"

        printf "%s\n" "${hap_id}" > "${sample_dir}/${hap_base}.query.txt"

        blastn \
            -query "${hap}" \
            -db "${sample_dir}/refdb" \
            -outfmt "6 std qlen slen" \
            -num_threads "${THREADS}" \
            > "${sample_dir}/${hap_base}.blastn.tsv"

        # Author-confirmed Figure 2D threshold: minimum alignment length 500 bp.
        "${B2D}" \
            -i1 "${sample_dir}/db.txt" \
            -i2 "${sample_dir}/${hap_base}.query.txt" \
            --blastn "${sample_dir}/${hap_base}.blastn.tsv" \
            --out "${sample_dir}/${hap_base}" \
            --min_identity 70 \
            --min_alignlen 500 \
            --line_width 13.0 \
            --figure_size 6 6 \
            --font_size 10 \
            --tick_label_size 8 \
            --colormap 4

        echo "  Created: ${sample_dir}/${hap_base}.pdf"
    done
done

echo
echo "All done."
echo "Outputs are in: ${OUTROOT}"
```
