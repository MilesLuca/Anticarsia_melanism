# Haplotype variant calling and genotype visualization


**Paper:** Figure 2E; Figure S11B  

Align each haploid assembly to chromosome 8, retain calls with exactly one aligned haplotype sequence (DP=1), merge samples and visualize informative biallelic variants.

## Inputs

- Pacbio_Reseq_Hilo/ivory_haplotypes/NC_134752.1.fa
- Pacbio_Reseq_Hilo/ivory_haplotypes/Agem*.fa (prepared cortex-containing haplotype FASTAs)
- Pacbio_Reseq_Hilo/ivory_haplotypes/04_genotype_plot/Agem_haplotypes.popmap.tsv

Input identifiers 01/02/05/06 correspond to paper individuals 1/2/3/4, respectively, for both morphs. VCF sample names, filenames and the popmap below retain the original identifiers so they remain consistent. This renumbering does not change hap1/hap2 identity or specify the final A/B/C row order.

## Outputs

- Pacbio_Reseq_Hilo/ivory_haplotypes/01_bams/*.mm2asm20nlj.sort.bam
- Pacbio_Reseq_Hilo/ivory_haplotypes/03_merged_DP1_vcf/Agem_haplotypes.DP1.merged.varsites.biallelic.vcf.gz
- Agem_haplotypes_ROI_genotype_plot.nomissing.GWA_marked.pdf

## Software

minimap2, SAMtools, bcftools, R, GenotypePlot and ggplot2. The final VCF header records bcftools/HTSlib 1.15.1.

## Execution and interpretation

1. Run scripts 01 → 02 → 03 → 05. Script 05 reads the merged VCF directly; script 04 is not a dependency of this saved result.
2. The final VCF header confirms mpileup --skip-indels --min-MQ 60, haploid calling, DP=1 filtering, merge and biallelic filtering.
3. Individual haplotype VCFs proceed from DP=1 filtering to merging.
4. The latest boundary-marked plotting script retains the core genotype settings used by the earlier plots.
5. The plotting boundaries are 2,016,199–2,611,793; the TE/population-statistic interval is 2,012,000–2,613,000.
6. The pooled haplotype-frequency contribution and separately reconstructed pie helper are included in [analysis 15](15_poolseq_haplotype_frequencies.md).

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Align and index haploid assemblies

Uses minimap2 asm20 with --no-long-join and per-haplotype read groups.

**Source:** `Pacbio_Reseq_Hilo/ivory_haplotypes/scripts/01_align_sort_index.sh`  
**Save as:** `01_align_sort_index.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-40}"
set -euo pipefail

cd /path/to/project/Pacbio_Reseq_Hilo/ivory_haplotypes


REF=NC_134752.1.fa
OUTDIR=01_bams
TMPDIR=${PWD}/tmp_sort

mkdir -p "$OUTDIR" "$TMPDIR"

echo "Working directory: $(pwd)"
echo "Reference: $REF"
echo "Start time: $(date)"

for FASTA in Agem_*.fa; do
    SAMPLE=$(basename "$FASTA" .fa)

    BAM="$OUTDIR/${SAMPLE}.mm2asm20nlj.bam"
    SORTBAM="$OUTDIR/${SAMPLE}.mm2asm20nlj.sort.bam"

    echo
    echo "=============================="
    echo "Processing $SAMPLE"
    echo "Input fasta: $FASTA"
    echo "Time: $(date)"

    minimap2 \
        -ax asm20 \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
        --no-long-join \
        -t ${THREADS} \
        "$REF" "$FASTA" \
    | samtools view -Sb - > "$BAM"

    samtools sort \
        -@ ${THREADS} \
        -T "$TMPDIR/${SAMPLE}" \
        "$BAM" > "$SORTBAM"

    samtools index "$SORTBAM"

    rm -f "$BAM"

    samtools flagstat "$SORTBAM" > "$OUTDIR/${SAMPLE}.flagstat.txt"
done

echo
echo "Sorted BAM count:"
ls -1 "$OUTDIR"/*.mm2asm20nlj.sort.bam | wc -l

echo
echo "BAM index count:"
ls -1 "$OUTDIR"/*.mm2asm20nlj.sort.bam.bai | wc -l

echo
echo "Done: $(date)"
```

<a id="step-02"></a>

### 2. Call haploid variants and filter DP=1

Produces one VCF per haplotype.

**Source:** `Pacbio_Reseq_Hilo/ivory_haplotypes/scripts/02_call_DP1_vcfs.sh`  
**Save as:** `02_call_DP1_vcfs.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-40}"
set -euo pipefail

cd /path/to/project/Pacbio_Reseq_Hilo/ivory_haplotypes


REF=NC_134752.1.fa
INDIR=01_bams
OUTDIR=02_DP1_vcfs

mkdir -p "$OUTDIR"

echo "Working directory: $(pwd)"
echo "Reference: $REF"
echo "Start time: $(date)"

for BAM in "$INDIR"/*.mm2asm20nlj.sort.bam; do
    SAMPLE=$(basename "$BAM" .mm2asm20nlj.sort.bam)
    OUTVCF="$OUTDIR/${SAMPLE}.DP1.vcf.gz"

    echo
    echo "=============================="
    echo "Processing $SAMPLE"
    echo "Input bam: $BAM"
    echo "Time: $(date)"

    bcftools mpileup \
        --threads ${THREADS} \
        --skip-indels \
        -f "$REF" \
        --min-MQ 60 \
        --annotate FORMAT/DP \
        "$BAM" \
    | bcftools call \
        --threads ${THREADS} \
        -m \
        --ploidy 1 \
    | bcftools filter \
        --threads ${THREADS} \
        -i 'FORMAT/DP == 1' \
    | bgzip -c > "$OUTVCF"

    tabix -f -p vcf "$OUTVCF"

    bcftools stats "$OUTVCF" > "$OUTDIR/${SAMPLE}.DP1.stats.txt"
done

echo
echo "VCF count:"
ls -1 "$OUTDIR"/*.DP1.vcf.gz | wc -l

echo
echo "Tabix index count:"
ls -1 "$OUTDIR"/*.DP1.vcf.gz.tbi | wc -l

echo
echo "Done: $(date)"
```

<a id="step-03"></a>

### 3. Merge the per-haplotype VCFs

Retains reference reassembly haplotypes alongside the sixteen F2 haplotypes.

**Source:** `Pacbio_Reseq_Hilo/ivory_haplotypes/scripts/03_merge_DP1_vcfs.sh`  
**Save as:** `03_merge_DP1_vcfs.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-32}"
set -euo pipefail

cd /path/to/project/Pacbio_Reseq_Hilo/ivory_haplotypes


INDIR=02_DP1_vcfs
OUTDIR=03_merged_DP1_vcf

mkdir -p "$OUTDIR"

ls -1 "$INDIR"/*.DP1.vcf.gz | sort > "$OUTDIR/dp1_vcf_list.txt"

bcftools merge \
    --threads ${THREADS} \
    -l "$OUTDIR/dp1_vcf_list.txt" \
    -Oz \
    -o "$OUTDIR/Agem_haplotypes.DP1.merged.vcf.gz"

tabix -f -p vcf "$OUTDIR/Agem_haplotypes.DP1.merged.vcf.gz"

echo
echo "Merged samples:"
bcftools query -l "$OUTDIR/Agem_haplotypes.DP1.merged.vcf.gz"

echo
echo "Done: $(date)"
```

<a id="step-04"></a>

### 4. Keep informative biallelic variants

The final input/output filenames agree with the retained VCF header.

**Source:** `Pacbio_Reseq_Hilo/ivory_haplotypes/scripts/05_filter_biallelic_variant_sites.sh`  
**Save as:** `05_filter_biallelic_variant_sites.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash

THREADS="${THREADS:-8}"
set -euo pipefail

cd /path/to/project/Pacbio_Reseq_Hilo/ivory_haplotypes


INDIR=03_merged_DP1_vcf
INVCF="$INDIR/Agem_haplotypes.DP1.merged.vcf.gz"
OUTVCF="$INDIR/Agem_haplotypes.DP1.merged.varsites.biallelic.vcf.gz"

echo "Working directory: $(pwd)"
echo "Input VCF: $INVCF"
echo "Output VCF: $OUTVCF"
echo "Start time: $(date)"

if [[ ! -f "$INVCF" ]]; then
    echo "ERROR: Input VCF not found: $INVCF"
    exit 1
fi

bcftools view \
    --threads ${THREADS} \
    -m2 -M2 \
    -v snps,indels \
    -i 'ALT!="." && ALT!="*"' \
    -Oz \
    -o "$OUTVCF" \
    "$INVCF"

tabix -f -p vcf "$OUTVCF"

echo
echo "Input site count:"
bcftools view -H "$INVCF" | wc -l

echo
echo "Biallelic variant-only site count:"
bcftools view -H "$OUTVCF" | wc -l

echo
echo "Done: $(date)"
```

<a id="step-05"></a>

### 5. GenotypePlot sample groups

Save as the original popmap filename. This groups Dark reference, Dark and Light samples.

**Source:** `Pacbio_Reseq_Hilo/ivory_haplotypes/04_genotype_plot/Agem_haplotypes.popmap.tsv`  
**Save as:** `Agem_haplotypes.popmap.tsv`  
**Archive treatment:** Complete source, verbatim.

```text
Agem_Reference_reassembly.hap1_cortex_contig	Dark_Ref
Agem_Reference_reassembly.hap2_cortex_contig	Dark_Ref
Agem_Dark_01.hap1_cortex_contig	Dark
Agem_Dark_01.hap2_cortex_contig	Dark
Agem_Dark_02.hap1_cortex_contig	Dark
Agem_Dark_02.hap2_cortex_contig	Dark
Agem_Dark_05.hap1_cortex_contig	Dark
Agem_Dark_05.hap2_cortex_contig	Dark
Agem_Dark_06.hap1_cortex_contig	Dark
Agem_Dark_06.hap2_cortex_contig	Dark
Agem_Light_01.hap1_cortex_contig	Light
Agem_Light_01.hap2_cortex_contig	Light
Agem_Light_02.hap1_cortex_contig	Light
Agem_Light_02.hap2_cortex_contig	Light
Agem_Light_05.hap1_cortex_contig	Light
Agem_Light_05.hap2_cortex_contig	Light
Agem_Light_06.hap1_cortex_contig	Light
Agem_Light_06.hap2_cortex_contig	Light
```

<a id="step-06"></a>

### 6. Plot genotypes and GWA boundaries

Latest recovered version with GWA boundaries marked.

**Source:** `Pacbio_Reseq_Hilo/ivory_haplotypes/04_genotype_plot/run_genotype_plot_coordinates.R`  
**Save as:** `run_genotype_plot_coordinates.R`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples.

```r
library(GenotypePlot)
library(ggplot2)

setwd("/path/to/project/Pacbio_Reseq_Hilo/ivory_haplotypes/04_genotype_plot")

vcf_file <- "/path/to/project/Pacbio_Reseq_Hilo/ivory_haplotypes/03_merged_DP1_vcf/Agem_haplotypes.DP1.merged.varsites.biallelic.vcf.gz"
popmap_file <- "Agem_haplotypes.popmap.tsv"

popmap <- read.table(
  popmap_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(popmap) <- c("ind", "pop")

gp <- genotype_plot(
  vcf = vcf_file,
  popmap = popmap,
  chr = "NC_134752.1",
  start = 1900000,
  end = 2810000,
  is.haploid = TRUE,
  is.multi_allelic = FALSE,
  cluster = FALSE,
  plot_phased = FALSE,
  missingness = 0.01,
  invariant_filter = TRUE,
  colour_scheme = c("#134c84", "#b0207c", "#f89c30")
)

## ----------------------------------------------------------------------------
## GWA boundaries from Pool-seq
## ----------------------------------------------------------------------------

gwa_start <- 2016199L
gwa_end   <- 2611793L

## ----------------------------------------------------------------------------
## Find nearest plotted SNP positions from the genotype panel itself
## ----------------------------------------------------------------------------

geno_x <- sort(unique(ggplot_build(gp$genotypes)$data[[1]]$x))

if (length(geno_x) == 0) {
  stop("Could not recover plotted SNP x positions from gp$genotypes.")
}

gwa_start_plot <- geno_x[which.min(abs(geno_x - gwa_start))]
gwa_end_plot   <- geno_x[which.min(abs(geno_x - gwa_end))]

cat("Requested GWA start:", gwa_start, "\n")
cat("Nearest plotted position:", gwa_start_plot,
    " | delta:", gwa_start_plot - gwa_start, "\n")

cat("Requested GWA end:", gwa_end, "\n")
cat("Nearest plotted position:", gwa_end_plot,
    " | delta:", gwa_end_plot - gwa_end, "\n")

## ----------------------------------------------------------------------------
## Add boundary lines directly to the positions and genotype panels
## ----------------------------------------------------------------------------

line_colour <- "black"
line_type   <- "dashed"
line_width  <- 0.7

gp$positions <- gp$positions +
  geom_vline(
    xintercept = c(gwa_start_plot, gwa_end_plot),
    linetype = line_type,
    linewidth = line_width,
    colour = line_colour
  )

gp$genotypes <- gp$genotypes +
  geom_vline(
    xintercept = c(gwa_start_plot, gwa_end_plot),
    linetype = line_type,
    linewidth = line_width,
    colour = line_colour
  )

## ----------------------------------------------------------------------------
## Save combined plot
## ----------------------------------------------------------------------------

pdf("Agem_haplotypes_ROI_genotype_plot.nomissing.GWA_marked.pdf",
    width = 14, height = 8)
print(combine_genotype_plot(gp))
dev.off()
```
