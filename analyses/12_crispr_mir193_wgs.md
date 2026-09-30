# mir-193 CRISPR whole-genome genotyping coverage

**Paper:** Figure S15 coverage  

Trim and map CRISPR short reads, quantify locus coverage and normalize by the chromosome mean, with the Light pool providing a wild-type population baseline.

## Inputs

- CRISPR_genotyping/wgs/fastqs/*_R{1,2}_001.fastq.gz and jobs
- Agem_genome/NCBI/GCF_050436995.1_ilAntGemm2_primary_genomic.fna with BWA index.
- CRISPR_genotyping/wgs/bams_NC_134752.1/plain.NC_134752.1.bam (prepared wild-type pool subset)

## Outputs

- Trimmed FASTQs, sorted/indexed chromosome BAMs, 1-bp mosdepth outputs.
- Separate normalized locus-coverage PDFs for Ag-miR-R-9, Ag-miR-R-F2 and the Light pool.

## Software

fastp, BWA-MEM, SAMtools, mosdepth; R/dplyr/ggplot2/tibble. Inspected mutant BAM headers record BWA 0.7.17 and SAMtools 1.22. The WT-pool header separately records upstream BWA 0.7.19-r1273 and SAMtools 1.22.1, followed by chromosome extraction with SAMtools 1.22.

## Execution and interpretation

1. Run fastp/mapping first. The mapping-only script is the later rerun using those trimmed FASTQs; it replaces the mapping stage rather than being an additional biological analysis.
2. Section 4a supplies the WT-pool chromosome extraction recovered from the retained BAM header, plus index creation. Start from an indexed, coordinate-sorted Light-pool whole-genome BAM and write the `plain.NC_134752.1.bam` input expected by section 6.
3. Depth is computed in 1-bp bins and divided by the chromosome mean, matching the CRISPR WGS Methods. This differs from the Pool-GWAS median-depth definition.
4. The final mutant-allele sequences shown in Figure S15 were prepared manually in Geneious. The code here covers read mapping and coverage analysis.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Sample list for the per-sample scripts

Save as fastqs/jobs; one sample per line.

**Source:** `CRISPR_genotyping/wgs/fastqs/jobs`  
**Save as:** `jobs`  
**Archive treatment:** Complete source, verbatim.

```text
Ag-miR-R-9
Ag-miR-R-F2
```

<a id="step-02"></a>

### 2. Trim paired reads and map

Provides the trimming stage and initial mapping.

**Source:** `CRISPR_genotyping/wgs/array_fastp_bwa_agem_trimmed.sh`  
**Save as:** `array_fastp_bwa_agem_trimmed.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Replace the job-array index with a required, range-checked sample-index argument. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

# Run once per sample: bash array_fastp_bwa_agem_trimmed.sh INDEX (INDEX 0 through 1).
TASK_ID="${1:?Supply a zero-based sample index (0-1)}"
if [[ ! "$TASK_ID" =~ ^[0-1]$ ]]; then
  echo "Sample index must be 0-1." >&2
  exit 2
fi
set -euo pipefail


REF="/path/to/project/Agem_genome/NCBI/GCF_050436995.1_ilAntGemm2_primary_genomic.fna"
FASTQDIR="/path/to/project/CRISPR_genotyping/wgs/fastqs"
JOBS="${FASTQDIR}/jobs"

TRIMDIR="/path/to/project/CRISPR_genotyping/wgs/trimmed_fastqs"
BAMDIR="/path/to/project/CRISPR_genotyping/wgs/bams_trimmed"
STATSDIR="${BAMDIR}/alnstats"
QCDIR="/path/to/project/CRISPR_genotyping/wgs/fastp_reports"

mkdir -p "${TRIMDIR}" "${BAMDIR}" "${STATSDIR}" "${QCDIR}"

# ---- sanity: jobs file exists
if [[ ! -s "${JOBS}" ]]; then
  echo "ERROR: jobs file missing/empty: ${JOBS}" >&2
  echo "Create it with:" >&2
  echo "  ls -1 ${FASTQDIR}/*_R1_001.fastq.gz | sed 's/_R1_001\\.fastq\\.gz$//' | xargs -n1 basename | sort -u > ${FASTQDIR}/jobs" >&2
  exit 1
fi

# Get sample name for this task
SAMPLE=$(sed -n "$((TASK_ID+1))p" "${JOBS}")
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: SAMPLE is empty for TASK_ID=${TASK_ID} (check ${JOBS})" >&2
  exit 1
fi

echo "[$(date)] SAMPLE=${SAMPLE}"

R1="${FASTQDIR}/${SAMPLE}_R1_001.fastq.gz"
R2="${FASTQDIR}/${SAMPLE}_R2_001.fastq.gz"

# Output trimmed fastqs
TRIM_R1="${TRIMDIR}/${SAMPLE}.trim.R1.fastq.gz"
TRIM_R2="${TRIMDIR}/${SAMPLE}.trim.R2.fastq.gz"

# fastp reports
FP_HTML="${QCDIR}/${SAMPLE}.fastp.html"
FP_JSON="${QCDIR}/${SAMPLE}.fastp.json"

# Sanity checks
if [[ ! -f "${R1}" ]]; then echo "Missing R1: ${R1}" >&2; exit 1; fi
if [[ ! -f "${R2}" ]]; then echo "Missing R2: ${R2}" >&2; exit 1; fi

# Step 1: Trim (paired-end)
fastp \
  --thread 40 \
  --detect_adapter_for_pe \
  -i "${R1}" -I "${R2}" \
  -o "${TRIM_R1}" -O "${TRIM_R2}" \
  --html "${FP_HTML}" \
  --json "${FP_JSON}"

# Step 2: Map trimmed reads -> sorted BAM
OUTBAM="${BAMDIR}/${SAMPLE}.trimmed.sorted.bam"

bwa mem -t 40 "${REF}" "${TRIM_R1}" "${TRIM_R2}" \
  | samtools sort -@ 40 -o "${OUTBAM}" -

# Step 3: Stats + index
samtools flagstat -@ 40 "${OUTBAM}" > "${STATSDIR}/${SAMPLE}.stats.out"
samtools index -@ 40 "${OUTBAM}"

echo "[$(date)] Done: ${OUTBAM}"

```

<a id="step-03"></a>

### 3. Final mapping rerun from trimmed reads

Uses the same sample list and produces the final named BAMs.

**Source:** `CRISPR_genotyping/wgs/array_bwa_agem_trimmed.sh`  
**Save as:** `array_bwa_agem_trimmed.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Replace the job-array index with a required, range-checked sample-index argument. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

# Run once per sample: bash array_bwa_agem_trimmed.sh INDEX (INDEX 0 through 1).
TASK_ID="${1:?Supply a zero-based sample index (0-1)}"
if [[ ! "$TASK_ID" =~ ^[0-1]$ ]]; then
  echo "Sample index must be 0-1." >&2
  exit 2
fi
set -euo pipefail


REF="/path/to/project/Agem_genome/NCBI/GCF_050436995.1_ilAntGemm2_primary_genomic.fna"

FASTQDIR="/path/to/project/CRISPR_genotyping/wgs/fastqs"
JOBS="${FASTQDIR}/jobs"

TRIMDIR="/path/to/project/CRISPR_genotyping/wgs/trimmed_fastqs"
BAMDIR="/path/to/project/CRISPR_genotyping/wgs/bams_trimmed"
STATSDIR="${BAMDIR}/alnstats"

mkdir -p "${BAMDIR}" "${STATSDIR}"

# ---- sanity: jobs file exists
if [[ ! -s "${JOBS}" ]]; then
  echo "ERROR: jobs file missing/empty: ${JOBS}" >&2
  exit 1
fi

# Get sample name for this task (jobs is 1 sample per line)
SAMPLE=$(sed -n "$((TASK_ID+1))p" "${JOBS}")
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: SAMPLE is empty for TASK_ID=${TASK_ID} (check ${JOBS})" >&2
  exit 1
fi

echo "[$(date)] SAMPLE=${SAMPLE}"

# Trimmed FASTQs produced by your fastp script
TRIM_R1="${TRIMDIR}/${SAMPLE}.trim.R1.fastq.gz"
TRIM_R2="${TRIMDIR}/${SAMPLE}.trim.R2.fastq.gz"

# Input checks
if [[ ! -f "${TRIM_R1}" ]]; then echo "Missing TRIM_R1: ${TRIM_R1}" >&2; exit 1; fi
if [[ ! -f "${TRIM_R2}" ]]; then echo "Missing TRIM_R2: ${TRIM_R2}" >&2; exit 1; fi

OUTBAM="${BAMDIR}/${SAMPLE}.trimmed.sorted.bam"

echo "[$(date)] Mapping -> ${OUTBAM}"
bwa mem -t 40 "${REF}" "${TRIM_R1}" "${TRIM_R2}" \
  | samtools sort -@ 40 -o "${OUTBAM}" -

echo "[$(date)] flagstat"
samtools flagstat -@ 40 "${OUTBAM}" > "${STATSDIR}/${SAMPLE}.stats.out"

echo "[$(date)] indexing"
samtools index -@ 40 "${OUTBAM}"

echo "[$(date)] Done"

```

<a id="step-04"></a>

### 4. Extract chromosome 8

The output directory is corrected to `bams_NC_134752.1`, matching the downstream depth scripts; the original source file is unchanged.

**Source:** `CRISPR_genotyping/wgs/bams_trimmed/array_subset_bams_chrom.sh`  
**Save as:** `array_subset_bams_chrom.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Replace the job-array index with a required, range-checked sample-index argument. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

# Run once per sample: bash array_subset_bams_chrom.sh INDEX (INDEX 0 through 1).
TASK_ID="${1:?Supply a zero-based sample index (0-1)}"
if [[ ! "$TASK_ID" =~ ^[0-1]$ ]]; then
  echo "Sample index must be 0-1." >&2
  exit 2
fi
set -euo pipefail


# Inputs
BAMDIR="/path/to/project/CRISPR_genotyping/wgs/bams_trimmed"
FASTQDIR="/path/to/project/CRISPR_genotyping/wgs/fastqs"
JOBS="${FASTQDIR}/jobs"

# Output
OUTDIR="/path/to/project/CRISPR_genotyping/wgs/bams_NC_134752.1"
mkdir -p "${OUTDIR}"

REGION="NC_134752.1"

# Get sample
if [[ ! -s "${JOBS}" ]]; then
  echo "ERROR: jobs file missing/empty: ${JOBS}" >&2
  exit 1
fi

SAMPLE=$(sed -n "$((TASK_ID+1))p" "${JOBS}")
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: SAMPLE is empty for TASK_ID=${TASK_ID}" >&2
  exit 1
fi

INBAM="${BAMDIR}/${SAMPLE}.trimmed.sorted.bam"
if [[ ! -f "${INBAM}" ]]; then
  echo "ERROR: missing input BAM: ${INBAM}" >&2
  exit 1
fi

OUTBAM="${OUTDIR}/${SAMPLE}.NC_134752.1.sorted.bam"

echo "[$(date)] SAMPLE=${SAMPLE}"
echo "[$(date)] INBAM=${INBAM}"
echo "[$(date)] REGION=${REGION}"
echo "[$(date)] OUTBAM=${OUTBAM}"

# Ensure input BAM is indexed (samtools view uses index for random access)
if [[ ! -f "${INBAM}.bai" ]]; then
  echo "[$(date)] Index missing; indexing input BAM..."
  samtools index -@ 8 "${INBAM}"
fi

# Subset, then sort (view output is not guaranteed coordinate-sorted)
samtools view -@ 8 -b "${INBAM}" "${REGION}" \
  | samtools sort -@ 8 -o "${OUTBAM}" -

# Index for IGV
samtools index -@ 8 "${OUTBAM}"

# Quick sanity stats
samtools idxstats "${OUTBAM}" > "${OUTBAM}.idxstats.txt"
samtools flagstat -@ 8 "${OUTBAM}" > "${OUTBAM}.flagstat.txt"

echo "[$(date)] Done"
ls -lh "${OUTBAM}" "${OUTBAM}.bai"

```

<a id="step-04a"></a>

### 4a. Prepare the wild-type Light-pool chromosome BAM

**Save as:** `subset_wt_chromosome8.sh`  
**Source evidence:** The `@PG` header in the retained `plain.NC_134752.1.bam` records `samtools view -@40 -o plain.NC_134752.1.bam plain.bam NC_134752.1`, executed with SAMtools 1.22. The original standalone script was not found.  
**Archive treatment:** Recover that recorded extraction command into a portable wrapper with explicit input/output paths. Index creation and output-existence checks are archive additions.

Run `bash subset_wt_chromosome8.sh /path/to/plain.bam /path/to/bams_NC_134752.1/plain.NC_134752.1.bam`, using an existing output directory and a new output filename. The input must be a coordinate-sorted, indexed Light-pool whole-genome BAM. Section 6 consumes this output.

```bash
#!/usr/bin/env bash
set -euo pipefail
# Usage: bash subset_wt_chromosome8.sh indexed_plain.bam NEW_OUTPUT.bam
INBAM=${1:?Supply the indexed, coordinate-sorted Light-pool BAM}
OUTBAM=${2:?Supply a new output BAM path}
THREADS=${THREADS:-40}
[[ "$OUTBAM" == *.bam ]] || { echo "Output must end in .bam" >&2; exit 2; }
[[ ! -e "$OUTBAM" && ! -e "$OUTBAM.bai" ]] || { echo "Output already exists" >&2; exit 1; }
# Extraction parameters recovered from the retained WT BAM's @PG record.
samtools view -@ "$THREADS" -o "$OUTBAM" "$INBAM" NC_134752.1
# Index the new output for downstream use; indexing is an added archive command.
samtools index -@ "$THREADS" "$OUTBAM"
```

<a id="step-05"></a>

### 5. Measure mutant coverage in 1-bp bins

Run with sample index 0 or 1, using the explicit mutant BAM order.

**Source:** `CRISPR_genotyping/wgs/bams_NC_134752.1/mosdepth.sh`  
**Save as:** `mosdepth.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Replace the job-array index with a required, range-checked sample-index argument. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

# Run once per sample: bash mosdepth.sh INDEX (INDEX 0 through 1).
TASK_ID="${1:?Supply a zero-based sample index (0-1)}"
if [[ ! "$TASK_ID" =~ ^[0-1]$ ]]; then
  echo "Sample index must be 0-1." >&2
  exit 2
fi
BAMDIR="/path/to/project/CRISPR_genotyping/wgs/bams_NC_134752.1"
OUTDIR="${BAMDIR}/mosdepth_1bp"
CHROM="NC_134752.1"
THREADS=40

mkdir -p "${OUTDIR}"

# Explicitly list just the two BAMs (stable ordering, no glob surprises)
bam_files=(
  "${BAMDIR}/Ag-miR-R-9.NC_134752.1.sorted.bam"
  "${BAMDIR}/Ag-miR-R-F2.NC_134752.1.sorted.bam"
)

bam_file="${bam_files[${TASK_ID}]}"
base_name="$(basename "${bam_file}" .bam)"

# Safety checks
if [[ ! -s "${bam_file}" ]]; then
  echo "ERROR: BAM not found or empty: ${bam_file}" >&2
  exit 1
fi

# Confirm index exists (no re-indexing; fail fast if missing)
if [[ ! -s "${bam_file}.bai" && ! -s "${bam_file%.bam}.bai" ]]; then
  echo "ERROR: BAM index (.bai) not found for: ${bam_file}" >&2
  echo "Expected: ${bam_file}.bai or ${bam_file%.bam}.bai" >&2
  exit 1
fi

echo "Running mosdepth on: ${bam_file}"
echo "Chrom: ${CHROM} | Bin size: 1 bp | Threads: ${THREADS}"
echo "Output prefix: ${OUTDIR}/1bp_${base_name}"

mosdepth \
  --threads "${THREADS}" \
  --chrom "${CHROM}" \
  --by 1 \
  "${OUTDIR}/1bp_${base_name}" \
  "${bam_file}"

```

<a id="step-06"></a>

### 6. Measure wild-type pool coverage

Uses the WT-pool chromosome BAM prepared in section 4a.

**Source:** `CRISPR_genotyping/wgs/bams_NC_134752.1/mosdepth_poolWT.sh`  
**Save as:** `mosdepth_poolWT.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
BAMDIR="/path/to/project/CRISPR_genotyping/wgs/bams_NC_134752.1"
BAM="${BAMDIR}/plain.NC_134752.1.bam"

OUTDIR="${BAMDIR}/mosdepth_1bp"
CHROM="NC_134752.1"
THREADS=40

mkdir -p "${OUTDIR}"

# fail fast if missing
[[ -s "${BAM}" ]] || { echo "ERROR: missing BAM: ${BAM}" >&2; exit 1; }

# confirm index exists (no re-indexing)
if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then
  echo "ERROR: BAM index (.bai) not found for: ${BAM}" >&2
  echo "Expected: ${BAM}.bai or ${BAM%.bam}.bai" >&2
  exit 1
fi

base_name="$(basename "${BAM}" .bam)"
echo "Running mosdepth 1bp for pooled control: ${base_name}"

mosdepth \
  --threads "${THREADS}" \
  --chrom "${CHROM}" \
  --by 1 \
  "${OUTDIR}/1bp_${base_name}" \
  "${BAM}"

```

<a id="step-07"></a>

### 7. Normalize and plot the mir-193 region

Plots the miRNA arms and the recorded sgRNA cut site at 2,764,784.

**Source:** `CRISPR_genotyping/wgs/bams_NC_134752.1/mosdepth.R`  
**Save as:** `mosdepth.R`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples.

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

# --------------------
# User settings
# --------------------
chrom <- "NC_134752.1"

region_start <- 2764700
region_end   <- 2764840

# Feature annotations
ann <- tribble(
  ~feature,     ~start,   ~end,
  "miR193_5p",  2764726,  2764747,
  "miR193_3p",  2764762,  2764783
)

cut_site <- 2764784

# Mosdepth output directory
cov_dir <- "/path/to/project/CRISPR_genotyping/wgs/bams_NC_134752.1/mosdepth_1bp"

# Prefixes (everything before .regions.bed.gz / .mosdepth.summary.txt)
# Now includes pooled WT control: plain.NC_134752.1
prefixes <- c(
  file.path(cov_dir, "1bp_Ag-miR-R-9.NC_134752.1.sorted"),
  file.path(cov_dir, "1bp_Ag-miR-R-F2.NC_134752.1.sorted"),
  file.path(cov_dir, "1bp_plain.NC_134752.1")
)

# --------------------
# Helpers
# --------------------
read_chr_mean_from_summary <- function(summary_file, chrom_name) {
  # Your mosdepth summary has a header and "mean" is a named column
  s <- read.table(summary_file, header = TRUE, sep = "", stringsAsFactors = FALSE)
  row <- s[s$chrom == chrom_name, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Could not uniquely find chrom '", chrom_name, "' in: ", summary_file)
  }
  as.numeric(row$mean)
}

read_regions_window <- function(regions_file, chrom_name, start_bp, end_bp) {
  # mosdepth regions: chr, start, end, depth (no header)
  df <- read.table(gzfile(regions_file), header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "depth")

  df %>%
    filter(chr == chrom_name,
           start >= start_bp,
           end   <= end_bp) %>%
    mutate(pos = (start + end) / 2)
}

plot_one <- function(prefix) {
  regions_file <- paste0(prefix, ".regions.bed.gz")
  summary_file <- paste0(prefix, ".mosdepth.summary.txt")

  if (!file.exists(regions_file)) stop("Missing regions file: ", regions_file)
  if (!file.exists(summary_file)) stop("Missing summary file: ", summary_file)

  chr_mean <- read_chr_mean_from_summary(summary_file, chrom)

  df <- read_regions_window(regions_file, chrom, region_start, region_end) %>%
    mutate(depth_norm = depth / chr_mean)

  sample_name <- sub("^1bp_", "", basename(prefix))

  ymax <- max(df$depth_norm, na.rm = TRUE)
  y_bar <- ymax * 1.12  # space above the coverage for annotations

  ggplot(df, aes(x = pos, y = depth_norm)) +
    # Filled area under curve (light blue)
    geom_area(fill = "lightblue", alpha = 0.6) +
    # Thin outline for clarity
    geom_line(linewidth = 0.2) +

    # miR annotations as horizontal bars above the coverage
    geom_segment(
      data = ann,
      aes(x = start, xend = end, y = y_bar, yend = y_bar),
      inherit.aes = FALSE,
      linewidth = 2
    ) +
    geom_text(
      data = ann,
      aes(x = (start + end)/2, y = y_bar * 1.02, label = feature),
      inherit.aes = FALSE,
      size = 3,
      vjust = 0
    ) +

    # sgRNA cut site
    geom_vline(xintercept = cut_site, linetype = "dashed", linewidth = 0.4) +
    annotate("text",
             x = cut_site, y = y_bar * 0.92,
             label = "sgRNA cut", angle = 90, vjust = -0.3, size = 3) +

    coord_cartesian(
      xlim = c(region_start, region_end),
      ylim = c(0, y_bar * 1.08)
    ) +
    labs(
      title = sample_name,
      x = paste0(chrom, " position (bp)"),
      y = "Depth / chromosome mean depth"
    ) +
    theme_classic(base_size = 11)
}

# --------------------
# Make plots (separate PDFs)
# --------------------
for (pfx in prefixes) {
  p <- plot_one(pfx)
  out_pdf <- paste0(
    "coverage_chrnorm_",
    basename(pfx), "_",
    region_start, "_", region_end,
    ".pdf"
  )
  ggsave(out_pdf, p, width = 7, height = 4)
  message("Wrote: ", out_pdf)
}

```
