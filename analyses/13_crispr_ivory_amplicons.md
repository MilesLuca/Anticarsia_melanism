# ivory CRISPR Nanopore amplicon coverage

**Paper:** Figure S16B  

Trim Nanopore amplicon reads, map them to the ivory amplicon reference and display read depth around both CRISPR targets.

## Inputs

- CRISPR_genotyping/Nanopore_Amplicons/new_runs/raw_fastq[s]/*.fastq.gz
- CRISPR_genotyping/Nanopore_Amplicons/new_runs/amplicon_ref/ivory_amplicon_reference.fasta

## Outputs

- Porechop-trimmed FASTQs, sorted/indexed minimap2 BAMs, 5-bp mosdepth outputs and per-sample coverage PDFs.

## Software

Porechop, minimap2, SAMtools, mosdepth; R with dplyr, ggplot2, tibble and readr. The inspected new_runs BAM records minimap2 2.24-r1122 and SAMtools 1.22.

## Execution and interpretation

1. Use new_runs. The older directory is explicitly named old_runs_potential_wrong_ids and has not been used to supply final allele calls.
2. The retained mapping builds a default minimap2 index and uses -ax map-ont. The historical parameters are preserved.
3. mosdepth uses default mean depth in 5-bp bins. Plot positions are amplicon-local BED coordinates; cut boundaries are 129.5 and 248.5.
4. The final S16B comparison uses Plain, G2_R16, G2_R15 and G2_R11. The source plotter also lists Dark and G1_R24; these two unshown plot calls are omitted from the archive copy.
5. Retained BAM provenance records minimap2 2.24-r1122 for this mapping step.
6. The final mutant-allele sequences shown in Figure S16 were prepared manually in Geneious. The code here covers amplicon mapping and coverage analysis.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Trim and map amplicons

The filename’s original spelling is retained.

**Source:** `CRISPR_genotyping/Nanopore_Amplicons/new_runs/porechop_miminap_amplicons.sh`  
**Save as:** `porechop_miminap_amplicons.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-8}"
set -euo pipefail

WD="/path/to/project/CRISPR_genotyping/Nanopore_Amplicons/new_runs"
REF_DIR="${WD}/amplicon_ref"

# Autodetect raw fastq directory name (raw_fastq vs raw_fastqs)
RAW_DIR=""
for d in "${WD}/raw_fastq" "${WD}/raw_fastqs"; do
  if [[ -d "$d" ]]; then
    RAW_DIR="$d"
    break
  fi
done
[[ -n "${RAW_DIR}" ]] || { echo "ERROR: cannot find raw_fastq or raw_fastqs directory in ${WD}"; exit 1; }
echo "[INFO] Using RAW_DIR=${RAW_DIR}"

REF="${REF_DIR}/ivory_amplicon_reference.fasta"

TRIM_DIR="${WD}/trimmed_fastqs"
BAM_DIR="${WD}/bams"
QC_DIR="${WD}/qc"

mkdir -p "${TRIM_DIR}" "${BAM_DIR}" "${QC_DIR}"
cd "${WD}"

# --- Required tools must be on PATH ---
command -v porechop >/dev/null 2>&1 || { echo "ERROR: porechop not found in PATH"; exit 1; }
command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 not found in PATH"; exit 1; }
command -v samtools  >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH"; exit 1; }

[[ -s "${REF}" ]] || { echo "ERROR: Reference fasta not found: ${REF}"; exit 1; }

# Build minimap2 index once (default params)
MMI="${REF}.mmi"
if [[ ! -s "${MMI}" ]]; then
  echo "[INFO] Building minimap2 index (default params): ${MMI}"
  minimap2 -d "${MMI}" "${REF}"
fi

# Ensure fasta index exists
if [[ ! -s "${REF}.fai" ]]; then
  echo "[INFO] Building samtools faidx: ${REF}.fai"
  samtools faidx "${REF}"
fi

echo "[INFO] Starting trimming + mapping..."
echo "[INFO] Raw fastqs: ${RAW_DIR}"

shopt -s nullglob
FASTQS=("${RAW_DIR}"/*.fastq.gz)
if [[ ${#FASTQS[@]} -eq 0 ]]; then
  echo "ERROR: No .fastq.gz files found in ${RAW_DIR}"
  echo "[DEBUG] Directory listing:"
  ls -lh "${RAW_DIR}" || true
  exit 1
fi

for FQ in "${FASTQS[@]}"; do
  BN="$(basename "${FQ}" .fastq.gz)"
  TRIM_FQ="${TRIM_DIR}/${BN}.trimmed.fastq.gz"
  BAM="${BAM_DIR}/${BN}.sorted.bam"

  echo "------------------------------------------------------------"
  echo "[INFO] Sample: ${BN}"
  echo "[INFO] Input : ${FQ}"
  echo "[INFO] Trim  : ${TRIM_FQ}"
  echo "[INFO] BAM   : ${BAM}"

  # 1) Trim adapters (porechop)
  porechop \
    -i "${FQ}" \
    -o "${TRIM_FQ}" \
    --threads "${THREADS}" \
    > "${QC_DIR}/${BN}.porechop.log" 2>&1

  # 2) Map (original params) + sort + index
  minimap2 -t "${THREADS}" -ax map-ont "${MMI}" "${TRIM_FQ}" \
    | samtools sort -@ "${THREADS}" -o "${BAM}" -

  samtools index -@ "${THREADS}" "${BAM}"

  # 3) Quick QC
  samtools flagstat "${BAM}" > "${QC_DIR}/${BN}.flagstat.txt"
  samtools idxstats "${BAM}" > "${QC_DIR}/${BN}.idxstats.txt"
  samtools view "${BAM}" | awk '$6 ~ /D/ {c++} END{print c+0}' > "${QC_DIR}/${BN}.cigar_has_D.count.txt"

  echo "[INFO] Done: ${BN}"
done

echo "[INFO] All samples complete."
echo "[INFO] Trimmed fastqs in: ${TRIM_DIR}"
echo "[INFO] BAMs in         : ${BAM_DIR}"
echo "[INFO] QC in           : ${QC_DIR}"

```

<a id="step-02"></a>

### 2. Calculate 5-bp depth

Uses sorted BAMs and writes full amplicon depth bins.

**Source:** `CRISPR_genotyping/Nanopore_Amplicons/new_runs/mosdepth.sh`  
**Save as:** `mosdepth.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-4}"
set -euo pipefail

WD="/path/to/project/CRISPR_genotyping/Nanopore_Amplicons/new_runs"
BAM_DIR="${WD}/bams"
OUT_DIR="${WD}/mosdepth_5bp"
REF="${WD}/amplicon_ref/ivory_amplicon_reference.fasta"

mkdir -p "${OUT_DIR}"
cd "${WD}"

command -v mosdepth >/dev/null 2>&1 || { echo "ERROR: mosdepth not found in PATH"; exit 1; }
command -v samtools  >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH"; exit 1; }

[[ -s "${REF}" ]] || { echo "ERROR: missing reference: ${REF}"; exit 1; }

# Ensure reference is indexed (useful for coordinate sanity; not strictly required for mosdepth)
if [[ ! -s "${REF}.fai" ]]; then
  samtools faidx "${REF}"
fi

shopt -s nullglob
BAMS=("${BAM_DIR}"/*.sorted.bam)
if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "ERROR: no *.sorted.bam files found in ${BAM_DIR}"
  ls -lh "${BAM_DIR}" || true
  exit 1
fi

echo "[INFO] Running mosdepth --by 5 on ${#BAMS[@]} BAM(s)"
echo "[INFO] Output dir: ${OUT_DIR}"

for BAM in "${BAMS[@]}"; do
  BN="$(basename "${BAM}" .sorted.bam)"
  PREFIX="${OUT_DIR}/${BN}.5bp"

  [[ -s "${BAM}.bai" ]] || samtools index -@ "${THREADS}" "${BAM}"

  echo "------------------------------------------------------------"
  echo "[INFO] Sample: ${BN}"
  echo "[INFO] BAM   : ${BAM}"
  echo "[INFO] Prefix: ${PREFIX}"

  # --by 5 => 5bp bins in *.regions.bed.gz
  # --fast-mode is fine for amplicons; remove if you prefer full mode
  mosdepth \
    --threads "${THREADS}" \
    --by 5 \
    "${PREFIX}" \
    "${BAM}"

done

echo "[INFO] Done. Regions files are in: ${OUT_DIR}"
echo "[INFO] Example output: *.regions.bed.gz  *.mosdepth.summary.txt"

```

<a id="step-03"></a>

### 3. Plot the four published coverage profiles

Unshown Dark/G1_R24 plot calls are omitted; calculations and plotting functions are unchanged.

**Source:** `CRISPR_genotyping/Nanopore_Amplicons/new_runs/mosdepth.R`  
**Save as:** `mosdepth.R`  
**Archive treatment:** Omit Dark and G1_R24 from the plotting prefix list; remove the preceding trailing comma. Preserve all per-sample computations. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples.

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(readr)
})

# --------------------
# User settings
# --------------------

# FASTA record name (matches amplicon_ref/ivory_amplicon_reference.fasta.fai)
chrom <- "ivory_Ref"

# Path to the FASTA index so we can pull the exact amplicon length
fai_file <- "/path/to/project/CRISPR_genotyping/Nanopore_Amplicons/new_runs/amplicon_ref/ivory_amplicon_reference.fasta.fai"

# Between-base cut boundaries (0-based axis): between 129/130 and between 248/249
cut_site1_boundary <- 129.5
cut_site2_boundary <- 248.5

# Optional feature annotations (0-based coordinates; leave empty if none)
ann <- tribble(
  ~feature, ~start, ~end
  # "sgRNA1_target", 120, 140,
  # "sgRNA2_target", 239, 259
)

# Mosdepth output directory and sample prefixes
cov_dir <- "/path/to/project/CRISPR_genotyping/Nanopore_Amplicons/new_runs/mosdepth_5bp"

# Prefixes are everything before .regions.bed.gz
prefixes <- c(
  file.path(cov_dir, "Plain.5bp"),
  file.path(cov_dir, "G2_R11.5bp"),
  file.path(cov_dir, "G2_R15.5bp"),
  file.path(cov_dir, "G2_R16.5bp")
)

# --------------------
# Helpers
# --------------------

get_amplicon_length_from_fai <- function(fai_path, chrom_name) {
  if (!file.exists(fai_path)) stop("Missing .fai file: ", fai_path)
  fai <- read_tsv(
    fai_path,
    col_names = c("seqname", "length", "offset", "linebases", "linewidth"),
    col_types = cols(
      seqname = col_character(),
      length = col_double(),
      offset = col_double(),
      linebases = col_double(),
      linewidth = col_double()
    ),
    progress = FALSE
  )
  row <- fai %>% filter(seqname == chrom_name)
  if (nrow(row) != 1) stop("Could not uniquely find '", chrom_name, "' in: ", fai_path)
  as.numeric(row$length[1])
}

read_regions_full <- function(regions_file, chrom_name) {
  # mosdepth regions: chr, start, end, depth (0-based BED start; end exclusive)
  df <- read_tsv(regions_file,
                 col_names = c("chr", "start", "end", "depth"),
                 col_types = cols(
                   chr = col_character(),
                   start = col_double(),
                   end = col_double(),
                   depth = col_double()
                 ),
                 progress = FALSE)

  df %>%
    filter(chr == chrom_name) %>%
    mutate(pos = (start + end) / 2)
}

plot_one <- function(prefix, amplicon_len) {
  regions_file <- paste0(prefix, ".regions.bed.gz")
  if (!file.exists(regions_file)) stop("Missing regions file: ", regions_file)

  df <- read_regions_full(regions_file, chrom)

  sample_name <- sub("\\.5bp$", "", basename(prefix))

  ymax <- max(df$depth, na.rm = TRUE)
  if (!is.finite(ymax)) ymax <- 0
  y_bar <- if (ymax > 0) ymax * 1.12 else 1

  ggplot(df, aes(x = pos, y = depth)) +
    geom_area(fill = "lightblue", alpha = 0.6) +
    geom_line(linewidth = 0.25) +

    # Feature bars (optional)
    {if (nrow(ann) > 0) geom_segment(
      data = ann,
      aes(x = start, xend = end, y = y_bar, yend = y_bar),
      inherit.aes = FALSE,
      linewidth = 2
    )} +
    {if (nrow(ann) > 0) geom_text(
      data = ann,
      aes(x = (start + end)/2, y = y_bar * 1.02, label = feature),
      inherit.aes = FALSE,
      size = 3,
      vjust = 0
    )} +

    # Exact between-base cut boundaries
    geom_vline(xintercept = cut_site1_boundary, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = cut_site2_boundary, linetype = "dashed", linewidth = 0.5) +
    annotate("text",
             x = cut_site1_boundary, y = y_bar * 0.92,
             label = "cut 1 (129|130)", angle = 90, vjust = -0.3, size = 3) +
    annotate("text",
             x = cut_site2_boundary, y = y_bar * 0.92,
             label = "cut 2 (248|249)", angle = 90, vjust = -0.3, size = 3) +

    coord_cartesian(
      xlim = c(0, amplicon_len),
      ylim = c(0, y_bar * 1.08)
    ) +
    labs(
      title = sample_name,
      x = paste0(chrom, " position (0-based; length=", amplicon_len, " bp)"),
      y = "Depth (mosdepth, 5 bp bins)"
    ) +
    theme_classic(base_size = 11)
}

# --------------------
# Run
# --------------------

amplicon_len <- get_amplicon_length_from_fai(fai_file, chrom)
message("[INFO] Amplicon length from FAI: ", amplicon_len, " bp")

for (pfx in prefixes) {
  p <- plot_one(pfx, amplicon_len)
  out_pdf <- paste0("coverage_5bp_", basename(pfx), "_full_amplicon.pdf")
  ggsave(out_pdf, p, width = 7, height = 4)
  message("Wrote: ", out_pdf)
}

```
