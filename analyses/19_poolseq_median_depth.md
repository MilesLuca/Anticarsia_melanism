# Pool-seq median depth and chromosome normalization

**Paper:** Figure 2A  

Calculate median read depth in non-overlapping 100-bp windows across chromosome 8. Normalize each pool independently by the median of its non-zero window depths across the whole chromosome, then plot the locus profiles.

This code was written for the archive on 14 September 2026, following their confirmation of median depth and whole-chromosome normalization and the definition in STAR Methods. It is a new implementation; the original Figure 2A scripts were not recovered.

## Inputs and software

- Final, coordinate-sorted Dark and Light pool BAMs aligned to the same Dark reference. Supply the BAMs after the intended upstream read filtering.
- Chromosome name `NC_134752.1` for Dark-reference chromosome 8. Both BAM headers must contain this chromosome with the same length.
- Bash, SAMtools, mosdepth and base R. Tools must be available on `PATH`; no cluster configuration is required.

The supplied BWA-MEM/MAPQ-20 preparation in [analysis 14](14_poolseq_gwas.md) can generate compatible BAM inputs using the Dark reference. These settings for both reference runs; identity with the historical depth-input BAMs has not been independently verified.

The depth command uses mosdepth's `--use-median --by 100`. The chromosome restriction covers the full chromosome. It retains standard mosdepth read handling: exclusion flag 1796, mapping-quality threshold 0, and normal CIGAR/paired-read overlap handling. These defaults do not replace the upstream filtering of the supplied BAMs. See the [mosdepth documentation](https://github.com/brentp/mosdepth#usage).

## Calculation and outputs

For pool *p*, let *dᵢ* be the median per-base depth in chromosome-8 window *i*. The baseline is `median(dᵢ where dᵢ > 0)`, using all chromosome-8 windows. The normalized depth is `dᵢ / baseline`.

- All windows, including zero-depth windows and a shorter terminal window if present, are retained in the normalized output.
- Zero-depth windows are excluded only from the baseline calculation.
- The whole-chromosome baseline is calculated before selecting the plotting interval. Each pool has its own baseline.
- The scripts write raw window medians, chromosome metadata, tool versions, a baseline table, a table of normalized depths for the whole chromosome, and a two-panel PDF.
- The default plot spans 1.5–3.5 Mb and shows window points with a dashed baseline at 1. Plotting limits are configurable and do not affect normalization. Display settings are new choices for this implementation.

## Run order

1. Save the first block as `calculate_pool_median_depth.sh`. Run `bash calculate_pool_median_depth.sh /path/to/dark.bam /path/to/light.bam /path/to/depth_output`. An optional fourth argument changes the chromosome; `THREADS=4` is the default. The script creates missing BAM indexes.
2. Save the second block as `normalize_plot_pool_depth.R`. Run `Rscript normalize_plot_pool_depth.R /path/to/depth_output /path/to/plot_output`. Optional third and fourth arguments set the plot start and end in bp, for example `1500000 3500000`.

Use whole-chromosome window files from the first step. The second step checks the complete 100-bp window grid against the chromosome length, so cropped inputs cannot silently change the baseline.

## Code

<a id="step-01"></a>

### 1. Calculate median depth across chromosome 8

**Save as:** `calculate_pool_median_depth.sh`  

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
CHROM=${4:-NC_134752.1}
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

### 2. Normalize by chromosome depth and plot

**Save as:** `normalize_plot_pool_depth.R`  

```r
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 4L)) {
  stop("Usage: Rscript normalize_plot_pool_depth.R DEPTH_DIR OUTPUT_DIR [START_BP END_BP]")
}
depth_dir <- args[1]
output_dir <- args[2]
plot_start <- if (length(args) == 4L) as.numeric(args[3]) else 1500000
plot_end <- if (length(args) == 4L) as.numeric(args[4]) else 3500000
if (any(!is.finite(c(plot_start, plot_end))) || plot_start < 0 || plot_end <= plot_start) {
  stop("Plot limits must be finite, non-negative and increasing.")
}

metadata <- read.delim(file.path(depth_dir, "chromosome.tsv"), stringsAsFactors = FALSE)
if (!identical(names(metadata), c("chromosome", "length_bp")) || nrow(metadata) != 1L) {
  stop("chromosome.tsv must contain one chromosome and its length_bp.")
}
chromosome <- as.character(metadata$chromosome[1])
chromosome_length <- as.numeric(metadata$length_bp[1])
if (is.na(chromosome) || !nzchar(chromosome) || !is.finite(chromosome_length) ||
    chromosome_length <= 0 || chromosome_length != floor(chromosome_length)) {
  stop("Invalid chromosome name or length.")
}
if (plot_end > chromosome_length) stop("Plot limits exceed the chromosome length.")

read_pool <- function(pool) {
  path <- file.path(depth_dir, paste0(pool, ".median100bp.regions.bed.gz"))
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  windows <- read.table(connection, sep = "\t", header = FALSE,
                        comment.char = "", quote = "", stringsAsFactors = FALSE)
  if (ncol(windows) != 4L || nrow(windows) == 0L) {
    stop(pool, ": expected a four-column mosdepth window BED.")
  }
  names(windows) <- c("chromosome", "start", "end", "median_depth")
  for (column in c("start", "end", "median_depth")) {
    windows[[column]] <- as.numeric(windows[[column]])
    if (any(!is.finite(windows[[column]]))) stop(pool, ": invalid ", column)
  }
  if (anyNA(windows$chromosome) || any(windows$chromosome != chromosome) ||
      any(windows$median_depth < 0)) stop(pool, ": wrong chromosome or negative depth.")
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

# Calculate each baseline from the whole chromosome before restricting the plot.
pools <- lapply(c("Dark", "Light"), read_pool)
normalized <- do.call(rbind, lapply(pools, `[[`, "windows"))
baselines <- do.call(rbind, lapply(pools, `[[`, "summary"))
plot_data <- normalized[normalized$end > plot_start & normalized$start < plot_end, ]
if (nrow(plot_data) == 0L) stop("No windows overlap the plotting interval.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(baselines, file.path(output_dir, "pool_chromosome_baselines.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(normalized, file.path(output_dir, "pool_normalized_depth.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

# Use the same y-axis scale for both pools; show every selected window.
colours <- c(Dark = "#675542", Light = "#9e8872")
ymax <- max(1.05, plot_data$normalized_depth) * 1.05
pdf(file.path(output_dir, "Figure2A_normalized_median_depth.pdf"), width = 10, height = 5)
par(mfrow = c(2, 1), mar = c(3.5, 4.5, 2, 1), las = 1)
for (pool in c("Dark", "Light")) {
  x <- plot_data[plot_data$pool == pool, ]
  plot(x$midpoint_bp / 1e6, x$normalized_depth, pch = 16, cex = 0.35,
       col = adjustcolor(colours[[pool]], alpha.f = 0.35),
       xlim = c(plot_start, plot_end) / 1e6, ylim = c(0, ymax), xaxs = "i", yaxs = "i",
       xlab = paste0(chromosome, " position (Mb)"), ylab = "Normalized median depth",
       main = paste(pool, "pool"), bty = "l")
  abline(h = 1, lty = 2, lwd = 1)
}
dev.off()
```
