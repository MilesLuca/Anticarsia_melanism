# Light-reference Pool-seq depth and chromosome normalization

**Paper:** Figure S7C  

Calculate median read depth in non-overlapping 100-bp windows across Light-reference chromosome 8. Divide each window by its pool's median non-zero window depth across the entire chromosome, then plot the locus profiles. Figure S7C used the same depth calculation and normalization as [Figure 2A](19_poolseq_median_depth.md).

These two blocks were written for the archive on 14 September 2026. The original S7C calculation, normalization and plotting scripts were not recovered. Retained whole-genome files named `median100bp_dark.newgenome.regions.bed` and `median100bp_plain.newgenome.regions.bed` supply the saved window depths; `plain` denotes the Light pool. The new normalizer was run directly on those files.

## Inputs and software

- For calculation from reads: coordinate-sorted Dark and Light pool BAMs aligned to the same Light reference, after the intended upstream read filtering. Their chromosome names and lengths must agree.
- Alternatively, start with existing four-column mosdepth BED files containing chromosome, zero-based start, exclusive end and median depth. Each must include the complete chromosome-wide 100-bp grid. Files may contain additional chromosomes and may be plain text or gzip-compressed.
- The retained depth files use `ptg000006l_1` for chromosome 8, length **13,938,792 bp**. The final assembly uses `chr8`; the correspondence is verified below. Both scripts default to the retained identifier. Set the chromosome argument explicitly when using `chr8`-named inputs.
- Bash, SAMtools, mosdepth and base R, with command-line tools available on `PATH`.

The calculation uses mosdepth `--by 100 --use-median`, exclusion flag 1796, mapping-quality threshold 0 and normal CIGAR/paired-read overlap handling, matching the explicit implementation defaults in analysis 19. The supplied BWA-MEM/MAPQ-20 preparation is now included in [analysis 14](14_poolseq_gwas.md), with `agem_light_06.fasta` as its contributed reference. Identity of those alignment outputs with the historical depth-input BAMs and the exact historical depth command remain unverified. Supply BAMs with the intended filtering already applied.

## Chromosome correspondence

| Reference record | Retained source | Length and orientation |
|---|---|---|
| `ptg000006l_1` | `agem_light_06.fasta` and the depth BEDs | 13,938,792 bp; forward |
| `chr8` | `agem_light_06.ragtag.final.fasta` | Identical sequence; forward, zero coordinate offset |

The original contig and chromosome 8 in both retained copies of the final FASTA are identical, including case. Their sequence SHA-256 is `aac948518dfe3bfcf6c686c591ec78832ff54c5d81f74df7309057305b88ac8b`. The retained provenance table independently records `chr8 / ptg000006l_1 / +`. Thus these depth coordinates can be displayed as chromosome 8 without reversal or an offset. This check covers chromosome 8 only; it does not establish whole-reference identity or equivalence of fresh whole-genome alignments. The Light_06 assembly belongs to the paper's Light 4 individual; see the [sample mapping](../README.md#sample-names).

## Calculation and retained-data check

For each pool, `baseline = median(window_median_depth[window_median_depth > 0])`, using the complete chromosome. Every window is then divided by that pool's baseline. Zero-depth windows remain in the normalized output; only the baseline calculation excludes them. The final 92-bp window is retained with the same weight as other windows when calculating the baseline.

| Pool | All chromosome-8 windows | Non-zero windows | Baseline from retained depths |
|---|---|---|---|
| Dark | 139,388 | 137,530 | 87 |
| Light (`plain`) | 139,388 | 137,568 | 116 |

These are newly calculated checks of the saved depths under the confirmed method, rather than recovered historical normalization outputs. The script computes baselines from its inputs; the values above are not hard-coded.

The plot defaults to 1.75–2.85 Mb, around the interval displayed in S7C. Plot bounds affect display only. The new two-panel scatter plot shows every selected normalized window and a dashed baseline at 1. Colours, axis scaling and display bounds are implementation choices. The original smoothing settings were not recovered, so this implementation does not add a fitted curve or claim to reproduce the exact published artwork.

## Run order and outputs

1. Save section 1 as `calculate_light_reference_depth.sh`. Run `bash calculate_light_reference_depth.sh /path/to/dark.bam /path/to/light.bam /path/to/new_depth_output`. An optional fourth argument sets the chromosome name. The script creates missing BAM indexes and writes chromosome metadata, tool versions and `Dark.median100bp.regions.bed.gz` / `Light.median100bp.regions.bed.gz`. Set `THREADS` to change the default of 4.
2. Save section 2 as `normalize_plot_light_reference_depth.R`. Run `Rscript normalize_plot_light_reference_depth.R /path/to/new_depth_output/Dark.median100bp.regions.bed.gz /path/to/new_depth_output/Light.median100bp.regions.bed.gz /path/to/new_plot_output`. To use the retained window files directly, replace the first two arguments with `median100bp_dark.newgenome.regions.bed` and `median100bp_plain.newgenome.regions.bed` in their location on your computer; skip section 1.
3. For final-assembly chromosome names, use `chr8` as the calculation script's fourth argument and append `chr8 13938792` to the R command. To change plot bounds as well, append both bounds, for example `chr8 13938792 1800000 2800000`. Use the chromosome length from the matching reference index or BAM header. Both scripts require a new output directory.

The normalizer writes `pool_chromosome_baselines.tsv`, `pool_normalized_depth.tsv` for the complete chromosome, `plot_settings.tsv`, `R_sessionInfo.txt` and `FigureS7C_normalized_median_depth.pdf`. Its input checks reject missing/duplicate windows, an incorrect chromosome length, negative/non-finite depth and a chromosome with no non-zero windows.

## Code

<a id="step-01"></a>

### 1. Calculate median depth across Light-reference chromosome 8

**Save as:** `calculate_light_reference_depth.sh`  

```bash
#!/usr/bin/env bash
set -euo pipefail

if (( $# < 3 || $# > 4 )); then
    echo "Usage: bash $0 DARK.bam LIGHT.bam OUTPUT_DIR [CHROMOSOME]" >&2
    exit 2
fi
DARK_BAM=$1
LIGHT_BAM=$2
OUTDIR=$3
CHROM=${4:-ptg000006l_1}
THREADS=${THREADS:-4}

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive." >&2; exit 2; }
for tool in samtools mosdepth; do
    command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done

chromosome_length() {
    samtools view -H "$1" | awk -F '\t' -v chrom="$CHROM" '
        $1 == "@SQ" {
            name = ""; len = ""
            for (i = 2; i <= NF; i++) {
                if ($i ~ /^SN:/) name = substr($i, 4)
                if ($i ~ /^LN:/) len = substr($i, 4)
            }
            if (name == chrom) print len
        }'
}

samtools quickcheck -v "$DARK_BAM" "$LIGHT_BAM"
DARK_LENGTH=$(chromosome_length "$DARK_BAM")
LIGHT_LENGTH=$(chromosome_length "$LIGHT_BAM")
if [[ ! "$DARK_LENGTH" =~ ^[1-9][0-9]*$ || "$DARK_LENGTH" != "$LIGHT_LENGTH" ]]; then
    echo "Both BAMs must contain $CHROM with the same positive chromosome length." >&2
    exit 1
fi

[[ ! -e "$OUTDIR" ]] || { echo "Use a new output directory: $OUTDIR" >&2; exit 1; }
mkdir -p "$OUTDIR"
printf 'chromosome\tlength_bp\n%s\t%s\n' "$CHROM" "$DARK_LENGTH" > "$OUTDIR/chromosome.tsv"
{
    mosdepth --version
    samtools --version
} > "$OUTDIR/tool_versions.txt" 2>&1

for pool in Dark Light; do
    if [[ "$pool" == Dark ]]; then bam=$DARK_BAM; else bam=$LIGHT_BAM; fi
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" && ! -f "${bam}.csi" ]]; then
        samtools index -@ "$THREADS" "$bam"
    fi
    mosdepth --threads "$THREADS" --chrom "$CHROM" \
        --by 100 --use-median --no-per-base --flag 1796 --mapq 0 \
        "$OUTDIR/${pool}.median100bp" "$bam"
done
```

<a id="step-02"></a>

### 2. Normalize complete chromosome depths and plot the locus

**Save as:** `normalize_plot_light_reference_depth.R`  
**Provenance:** Newly written base-R implementation, validated on the retained depth files and synthetic inputs. Whole-genome files are read in chunks and only the requested chromosome is retained.

```r
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 5L, 7L)) {
  stop(paste("Usage: Rscript normalize_plot_light_reference_depth.R",
             "DARK.regions.bed[.gz] LIGHT.regions.bed[.gz] OUTPUT_DIR",
             "[CHROMOSOME LENGTH_BP [START_BP END_BP]]"))
}
input_paths <- setNames(args[1:2], c("Dark", "Light"))
output_dir <- args[3]
chromosome <- if (length(args) >= 5L) args[4] else "ptg000006l_1"
chromosome_length <- if (length(args) >= 5L) as.numeric(args[5]) else 13938792
plot_start <- if (length(args) == 7L) as.numeric(args[6]) else 1750000
plot_end <- if (length(args) == 7L) as.numeric(args[7]) else 2850000
if (!nzchar(chromosome) || grepl("[[:space:]]", chromosome) ||
    !is.finite(chromosome_length) || chromosome_length <= 0 ||
    chromosome_length != floor(chromosome_length)) stop("Invalid chromosome or length.")
if (any(!is.finite(c(plot_start, plot_end))) || plot_start < 0 ||
    plot_end <= plot_start || plot_end > chromosome_length) {
  stop("Plot limits must be increasing and within the chromosome.")
}
if (file.exists(output_dir)) stop("Use a new output directory: ", output_dir)

read_pool <- function(pool) {
  path <- input_paths[[pool]]
  connection <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  # Stream whole-genome files, retaining only the requested chromosome.
  selected <- list()
  repeat {
    chunk <- readLines(connection, n = 100000L, warn = FALSE)
    if (!length(chunk)) break
    keep <- startsWith(chunk, paste0(chromosome, "\t"))
    if (any(keep)) selected[[length(selected) + 1L]] <- chunk[keep]
  }
  if (!length(selected)) stop(pool, ": chromosome not found: ", chromosome)
  windows <- read.table(text = unlist(selected, use.names = FALSE), sep = "\t",
                        header = FALSE, comment.char = "", quote = "",
                        colClasses = "character", fill = FALSE)
  if (ncol(windows) != 4L) stop(pool, ": expected four-column mosdepth window BED.")
  names(windows) <- c("chromosome", "start", "end", "median_depth")
  for (column in c("start", "end", "median_depth")) {
    windows[[column]] <- as.numeric(windows[[column]])
    if (any(!is.finite(windows[[column]]))) stop(pool, ": invalid ", column)
  }
  if (any(windows$median_depth < 0)) stop(pool, ": negative depth.")
  windows <- windows[order(windows$start), ]
  expected_start <- seq(0, chromosome_length - 1, by = 100)
  if (nrow(windows) != length(expected_start) ||
      any(windows$start != expected_start) ||
      any(windows$end != pmin(expected_start + 100, chromosome_length))) {
    stop(pool, ": supply the complete chromosome-wide 100-bp window grid.")
  }

  nonzero <- windows$median_depth > 0
  if (!any(nonzero)) stop(pool, ": no non-zero windows for normalization.")
  baseline <- median(windows$median_depth[nonzero])
  windows$pool <- pool
  windows$baseline_depth <- baseline
  windows$normalized_depth <- windows$median_depth / baseline
  windows$midpoint_bp <- (windows$start + windows$end) / 2
  summary <- data.frame(pool = pool, chromosome = chromosome,
                        chromosome_length_bp = chromosome_length,
                        windows_total = nrow(windows), windows_nonzero = sum(nonzero),
                        baseline_depth = baseline)
  list(windows = windows, summary = summary)
}

# Compute baselines before selecting a region for display.
pools <- lapply(names(input_paths), read_pool)
normalized <- do.call(rbind, lapply(pools, `[[`, "windows"))
baselines <- do.call(rbind, lapply(pools, `[[`, "summary"))
plot_data <- normalized[normalized$end > plot_start & normalized$start < plot_end, ]
if (!nrow(plot_data)) stop("No windows overlap the plotting interval.")
if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create output directory.")
write.table(baselines, file.path(output_dir, "pool_chromosome_baselines.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(normalized, file.path(output_dir, "pool_normalized_depth.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(chromosome = chromosome, display_chromosome = "chr8",
                       plot_start_bp = plot_start, plot_end_bp = plot_end),
            file.path(output_dir, "plot_settings.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

# A new scatter display of the confirmed calculation; all selected depths are shown.
colours <- c(Dark = "#675542", Light = "#9e8872")
ymax <- max(1.05, plot_data$normalized_depth) * 1.05
pdf(file.path(output_dir, "FigureS7C_normalized_median_depth.pdf"), width = 10, height = 5)
par(mfrow = c(2, 1), mar = c(3.5, 4.5, 2, 1), las = 1)
for (pool in names(input_paths)) {
  x <- plot_data[plot_data$pool == pool, ]
  plot(x$midpoint_bp / 1e6, x$normalized_depth, pch = 16, cex = 0.35,
       col = adjustcolor(colours[[pool]], alpha.f = 0.35),
       xlim = c(plot_start, plot_end) / 1e6, ylim = c(0, ymax), xaxs = "i", yaxs = "i",
       xlab = "Chromosome 8, Light reference (Mb)", ylab = "Normalized median depth",
       main = paste(pool, "pool"), bty = "l")
  abline(h = 1, lty = 2, lwd = 1)
}
dev.off()
```
