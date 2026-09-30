# Pool-seq GWAS mapping, statistics and plotting

**Paper:** Figure 1E; GWAS components of Figures 2C, 3B–C and S7–S9; upstream pool fractions for Figure 2F and S11  

## Inputs and software

- Paired pool reads: `Ag-PLAIN_R1_001.fastq.gz`, `Ag-PLAIN_R2_001.fastq.gz`, `Ag-DARK_R1_001.fastq.gz` and `Ag-DARK_R2_001.fastq.gz`. Plain means Light. The Light R1 may instead be uncompressed; section 1 compresses a copy when needed. The pools contain 23 Light and 25 Dark females, as reported in Methods; these scripts do not take pool-size parameters.
- The matching indexed reference FASTA. Build its BWA index with `bwa index /path/to/reference.fasta` if needed. For the Light run, the supplied source specifically uses `agem_light_06.fasta`; retain its contig names. Its chromosome-8 correspondence with the final Light assembly is described in [analysis 22](22_light_reference_depth.md). The Dark-reference resource is identified in the [README](../README.md); its original local FASTA filename was not supplied.
- Bash, BWA, SAMtools, pigz, Java and Perl, available on `PATH`. Use a complete [PoPoolation2 installation](https://github.com/popgenvienna/popoolation2), including `mpileup2sync.jar`, both Perl scripts and their `Modules` directory. The Fisher implementation also requires [Text::NSP](https://metacpan.org/pod/distribution/Text-NSP/lib/Text/NSP/Measures/2D/Fisher.pm).

## Settings and output interpretation

| Stage | Supplied settings and interpretation |
|---|---|
| Alignment | BWA-MEM, 40 threads by default; retain BAM alignments with MAPQ >=20 and coordinate-sort. The BAM filter has no additional duplicate, proper-pair or secondary/supplementary exclusion. Downstream tools retain their own read-handling defaults. |
| Pileup | `samtools mpileup -B`, Light BAM first and Dark BAM second. `-B` disables BAQ. |
| Conversion | Java with assertions and 7-GB heap; Sanger qualities, minimum base quality 20, 40 threads by default. The SYNC sample order is Light then Dark, with counts in A:T:C:G:N:deletion order. |
| SNP frequencies and Fisher tests | `--min-count 6 --min-coverage 50 --max-coverage 200`. The minimum minor-allele count is assessed across the pools together; coverage limits apply to each pool. The two commands read the same SYNC independently. They do not rewrite it as a filtered input. |
| Fisher windows | No window flags are supplied. The checked upstream defaults are **1-bp window and 1-bp step**, so the calculation is per-site. |


## Run order

Save the five blocks under their stated filenames. Use one new output directory for the Light-reference run and another for the Dark-reference run. Within each run, both pool scripts must receive the same reference and output directory. `THREADS` defaults to 40 for alignment/conversion and can be reduced for the available hardware.

1. Run `bash bwa.plain.newgenome.sh /path/to/reference.fasta /path/to/fastq_files /path/to/run_outputs`.
2. Run `bash bwa.dark.newgenome.sh /path/to/reference.fasta /path/to/fastq_files /path/to/run_outputs`.
3. Run `bash sync.sh /path/to/run_outputs /path/to/popoolation2/mpileup2sync.jar`.
4. After verifying conversion completed, run `bash snpFreqDif.sh /path/to/run_outputs /path/to/popoolation2`.
5. Run `bash FisherExact.sh /path/to/run_outputs /path/to/popoolation2`.

## Code

<a id="step-01"></a>

### 1. Align the Light pool

**Save as:** `bwa.plain.newgenome.sh`  

```bash
#!/usr/bin/env bash
set -euo pipefail
if (( $# != 3 )); then
    echo "Usage: bash $0 REFERENCE.fasta FASTQ_DIR OUTPUT_DIR" >&2
    exit 2
fi
REFERENCE=$1
FASTQ_DIR=$2
OUTDIR=$3
THREADS=${THREADS:-40}
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive." >&2; exit 2; }
mkdir -p "$OUTDIR"
R1="$FASTQ_DIR/Ag-PLAIN_R1_001.fastq.gz"
R2="$FASTQ_DIR/Ag-PLAIN_R2_001.fastq.gz"
# Use compressed input if available; otherwise preserve the original FASTQ.
if [[ ! -s "$R1" ]]; then
    R1="$OUTDIR/Ag-PLAIN_R1_001.fastq.gz"
    [[ ! -e "$R1" ]] || { echo "Compressed output already exists: $R1" >&2; exit 1; }
    pigz -p 8 -c "$FASTQ_DIR/Ag-PLAIN_R1_001.fastq" > "$R1"
fi
[[ -s "$REFERENCE" && -s "$R1" && -s "$R2" ]] || { echo "Missing reference or FASTQ input." >&2; exit 1; }
[[ ! -e "$OUTDIR/plain.newgenome.sam" && ! -e "$OUTDIR/plain.newgenome.bam" ]] || { echo "Use a new alignment output location." >&2; exit 1; }
# The matching BWA index must already exist beside REFERENCE.
bwa mem -t "$THREADS" "$REFERENCE" "$R1" "$R2" > "$OUTDIR/plain.newgenome.sam"
samtools view -q 20 -bS "$OUTDIR/plain.newgenome.sam" | \
    samtools sort - -o "$OUTDIR/plain.newgenome.bam"
```

<a id="step-02"></a>

### 2. Align the Dark pool

**Save as:** `bwa.dark.newgenome.sh`  

```bash
#!/usr/bin/env bash
set -euo pipefail
if (( $# != 3 )); then
    echo "Usage: bash $0 REFERENCE.fasta FASTQ_DIR OUTPUT_DIR" >&2
    exit 2
fi
REFERENCE=$1
FASTQ_DIR=$2
OUTDIR=$3
THREADS=${THREADS:-40}
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive." >&2; exit 2; }
mkdir -p "$OUTDIR"
R1="$FASTQ_DIR/Ag-DARK_R1_001.fastq.gz"
R2="$FASTQ_DIR/Ag-DARK_R2_001.fastq.gz"
[[ -s "$REFERENCE" && -s "$R1" && -s "$R2" ]] || { echo "Missing reference or FASTQ input." >&2; exit 1; }
[[ ! -e "$OUTDIR/dark.newgenome.sam" && ! -e "$OUTDIR/dark.newgenome.bam" ]] || { echo "Use a new alignment output location." >&2; exit 1; }
# The matching BWA index must already exist beside REFERENCE.
bwa mem -t "$THREADS" "$REFERENCE" "$R1" "$R2" > "$OUTDIR/dark.newgenome.sam"
samtools view -q 20 -bS "$OUTDIR/dark.newgenome.sam" | \
    samtools sort - -o "$OUTDIR/dark.newgenome.bam"
samtools index "$OUTDIR/dark.newgenome.bam"
```

<a id="step-03"></a>

### 3. Generate the two-pool SYNC file

**Save as:** `sync.sh`  

```bash
#!/usr/bin/env bash
set -euo pipefail
if (( $# != 2 )); then
    echo "Usage: bash $0 WORK_DIR /path/to/popoolation2/mpileup2sync.jar" >&2
    exit 2
fi
WORK_DIR=$1
# Resolve the jar before changing directory; relative and absolute paths work.
SYNC_JAR=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
THREADS=${THREADS:-40}
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive." >&2; exit 2; }
[[ -s "$SYNC_JAR" ]] || { echo "Missing mpileup2sync.jar." >&2; exit 1; }
cd "$WORK_DIR"
[[ ! -e PlainDark.mpileup && ! -e PlainDark_java.sync ]] || { echo "Output already exists; use a fresh work directory." >&2; exit 1; }
# Restored from the supplied commented command. Keep Light first, Dark second.
# -B disables BAQ; the supplied command has no reference (-f) argument.
samtools mpileup -B plain.newgenome.bam dark.newgenome.bam > PlainDark.mpileup
java -ea -Xmx7g -jar "$SYNC_JAR" \
    --input PlainDark.mpileup --output PlainDark_java.sync \
    --fastq-type sanger --min-qual 20 --threads "$THREADS"
```

<a id="step-04"></a>

### 4. Calculate SNP allele fractions

**Save as:** `snpFreqDif.sh`  

```bash
#!/usr/bin/env bash
set -euo pipefail
if (( $# != 2 )); then
    echo "Usage: bash $0 WORK_DIR /path/to/popoolation2" >&2
    exit 2
fi
WORK_DIR=$1
POPOOLATION2_DIR=$(cd "$2" && pwd)
cd "$WORK_DIR"
[[ -s PlainDark_java.sync ]] || { echo "Missing PlainDark_java.sync." >&2; exit 1; }
for output in PlainDark_rc PlainDark_pwc PlainDark.params; do
    [[ ! -e "$output" ]] || { echo "Output already exists: $output" >&2; exit 1; }
done
perl "$POPOOLATION2_DIR/snp-frequency-diff.pl" --input PlainDark_java.sync --output-prefix PlainDark \
    --min-count 6 --min-coverage 50 --max-coverage 200
```

<a id="step-05"></a>

### 5. Calculate per-site Fisher tests

**Save as:** `FisherExact.sh`  

```bash
#!/usr/bin/env bash
set -euo pipefail
if (( $# != 2 )); then
    echo "Usage: bash $0 WORK_DIR /path/to/popoolation2" >&2
    exit 2
fi
WORK_DIR=$1
POPOOLATION2_DIR=$(cd "$2" && pwd)
cd "$WORK_DIR"
[[ -s PlainDark_java.sync ]] || { echo "Missing PlainDark_java.sync." >&2; exit 1; }
for output in PlainDark_10000_1000.fet PlainDark_10000_1000.fet.params; do
    [[ ! -e "$output" ]] || { echo "Output already exists: $output" >&2; exit 1; }
done
# No window arguments were supplied: the checked defaults are 1 bp / 1 bp.
# The retained output filename does not specify the analysis window settings.
perl "$POPOOLATION2_DIR/fisher-test.pl" --input PlainDark_java.sync --output PlainDark_10000_1000.fet \
    --min-count 6 --min-coverage 50 --max-coverage 200 --suppress-noninformative
```

<a id="step-06"></a>

### 6. Plot the Dark-reference per-site GWAS and chromosome-8 detail views

**Save as:** `plot_dark_persite_gwas.R`  

Run `Rscript plot_dark_persite_gwas.R /path/to/Dark_reference_run/PlainDark_10000_1000.fet /path/to/sizes.genome /path/to/NC_134752.1.gff /path/to/new_plot_outputs`. The first argument can instead be the prepared `PlainDark_perSite.fet`. The reference-specific output directory from section 5 supplies the raw Fisher file; its retained filename does not change the per-site settings. This plotting example uses Dark-reference sequence names and coordinates.

Use R with ggplot2, dplyr, readr, cowplot, ggrastr and scales, including ggrastr's dependencies. R must support the Cairo PDF device. Main point layers are rasterized at 600 dpi; the three further detail exports use vector points. The plot definitions, dimensions, colors, site values, genomic intervals and annotation calculations are retained.

| Input | Contents and interpretation |
|---|---|
| Per-site FET | Headerless six-column table: scaffold, one-based position, retained SNP count, covered fraction, average coverage and `-log10(P)`. |
| `sizes.genome` | Two columns, scaffold and length, for the same Dark reference. Preserve chromosome order: `NC_134752.1` must be the eighth `NC_` record. It can be prepared from the first two columns of the matching FASTA index. |
| `NC_134752.1.gff` | Nine-column annotation with gene, exon and/or CDS records in the same coordinates. The named files retained with analysis 10 are candidates. |


| Export (all filenames begin `PlainDark_perSite.`) | Scope and size |
|---|---|
| `fullGenome_chr8Genes_chr8Zoom_shadedGWA.pdf` | Genome, annotation and chr8:1,500,000–3,500,000, with GWA shading at 2,012,000–2,613,000; 14 × 12 inches |
| `genomewide.pdf` | Genome-wide view; 14 × 3 inches |
| `chr8_zoom_with_annotation_shadedGWA.pdf` | Annotation and 1.5–3.5 Mb profile with GWA shading; 14 × 3 inches |
| `chr8_further_zoom_nonrasterized.pdf` | 1,750,707–2,893,709; 10 × 5 inches |
| `chr8_2.4Mb_2.7Mb_nonrasterized.pdf` | 2,400,000–2,700,000; 15 × 3.5 inches |
| `chr8_2.4Mb_2.7Mb_with_annotation_nonrasterized.pdf` | **2,400,000–2,800,000**, as specified by the actual code, despite the retained filename/comment; 654 × 265 PDF points |

The script also writes `PlainDark_perSite.parsed.tsv`, preserving the six named fields with numeric statistics. The main views share the genome-wide y maximum; detail views use their own regional maxima. Inclusive zoom boundaries, GWA shading, exon/CDS handling and all six PDF names are retained. The final detail view deliberately retains the supplied 2.8-Mb endpoint; no boundary correction is inferred from its filename.


```r
## ============================================================================
## Per-site Pool-seq GWAS plot for Dark assembly
## File: PlainDark_perSite.fet
##
## Plots:
##   1. Genome-wide Manhattan plot
##   2. Chr8 annotation panel
##   3. Chr8 zoom plot with shaded GWA interval
##
## Notes:
##   - Column 6 is already raw -log10(P)
##   - No significance cutoff is plotted
##   - No points are highlighted in red
##   - Points are rasterized with ggrastr for Illustrator-friendly PDFs
## ============================================================================

# Archive setup: explicit inputs and a new output directory; no scheduler needed.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript plot_dark_persite_gwas.R FET_FILE SIZES_FILE GFF_FILE NEW_OUTPUT_DIR")
}
input_paths <- vapply(args[1:3], normalizePath, character(1), mustWork = TRUE)
if (file.exists(args[4])) stop("Use a new output directory.")
if (!dir.create(args[4], recursive = TRUE)) stop("Cannot create output directory.")
setwd(normalizePath(args[4], mustWork = TRUE))

library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)
library(ggrastr)
library(scales)

## ----------------------------------------------------------------------------
## 1. Input files
## ----------------------------------------------------------------------------

fet_file   <- input_paths[1]
sizes_file <- input_paths[2]
gff_file   <- input_paths[3]

## ----------------------------------------------------------------------------
## 2. Read per-site FET data
##    Last column is already raw -log10(P)
## ----------------------------------------------------------------------------

fet_data <- read_tsv(
  fet_file,
  col_names = FALSE,
  col_types = cols(
    X1 = col_character(),
    X2 = col_double(),
    X3 = col_double(),
    X4 = col_double(),
    X5 = col_double(),
    X6 = col_character()  # accept numeric or raw PoPoolation2 1:2=<value>
  )
)

# Archive input checks: six columns and a nonempty per-site table.
readr::stop_for_problems(fet_data)
if (ncol(fet_data) != 6L || nrow(fet_data) == 0L) {
  stop("FET input must contain six columns and at least one site.")
}
colnames(fet_data) <- c(
  "scaffold",
  "pos",
  "n_sites",
  "freq_diff",
  "avg_cov",
  "logp"
)

# New format adapter: remove only the known two-pool tag; do not log-transform.
parsed_logp <- readr::parse_double(sub("^1:2=", "", fet_data$logp))
readr::stop_for_problems(parsed_logp)
fet_data$logp <- as.numeric(parsed_logp)
if (any(!is.finite(fet_data$logp) | fet_data$logp < 0) ||
    any(!is.finite(fet_data$pos) | fet_data$pos < 1 |
          fet_data$pos != floor(fet_data$pos)) ||
    any(!is.finite(fet_data$n_sites) | fet_data$n_sites != 1)) {
  stop("Require finite non-negative -log10(P), valid positions and one SNP per row.")
}

cat("\n==================== Diagnostics ====================\n")
cat("Input file:", fet_file, "\n")
cat("Number of sites:", nrow(fet_data), "\n")
cat("Minimum -log10(P):", min(fet_data$logp, na.rm = TRUE), "\n")
cat("Maximum -log10(P):", max(fet_data$logp, na.rm = TRUE), "\n")
cat("====================================================\n\n")

## ----------------------------------------------------------------------------
## 3. Read genome sizes
## ----------------------------------------------------------------------------

scaffold_lengths <- read_tsv(
  sizes_file,
  col_names = FALSE,
  col_types = cols(
    X1 = col_character(),
    X2 = col_double()
  )
)

readr::stop_for_problems(scaffold_lengths)
if (ncol(scaffold_lengths) != 2L) stop("Sizes input must have two columns.")
colnames(scaffold_lengths) <- c("scaffold", "length")
if (anyNA(scaffold_lengths$scaffold) || anyDuplicated(scaffold_lengths$scaffold) ||
    any(!is.finite(scaffold_lengths$length) | scaffold_lengths$length <= 0)) {
  stop("Sizes input needs unique scaffold names and positive finite lengths.")
}
nc_order <- scaffold_lengths$scaffold[grepl("^NC_", scaffold_lengths$scaffold)]
if (length(nc_order) < 8L || nc_order[8] != "NC_134752.1") {
  stop("Preserve Dark-reference chromosome order: chromosome 8 must be NC_134752.1.")
}
site_lengths <- scaffold_lengths$length[match(fet_data$scaffold, scaffold_lengths$scaffold)]
if (anyNA(site_lengths) || any(fet_data$pos > site_lengths)) {
  stop("FET scaffolds/positions must match the supplied genome sizes.")
}

# NC_ = chromosomes, NW_ = unplaced scaffolds
scaffold_lengths <- scaffold_lengths %>%
  mutate(type = ifelse(grepl("^NC_", scaffold), "chrom", "unplaced"))

## ----------------------------------------------------------------------------
## 4. Build chromosome/scaffold map
##    - each NC_ scaffold becomes chromosome 1..N in file order
##    - all NW_ scaffolds concatenated into pseudochromosome "Unplaced"
## ----------------------------------------------------------------------------

chrom_scaffolds <- scaffold_lengths %>%
  filter(type == "chrom") %>%
  mutate(
    chrom = as.character(row_number()),
    offset_within_chr = 0
  )

unplaced_scaffolds <- scaffold_lengths %>%
  filter(type == "unplaced") %>%
  arrange(scaffold) %>%
  mutate(
    chrom = "Unplaced",
    offset_within_chr = lag(cumsum(length), default = 0)
  )

scaffold_map <- bind_rows(chrom_scaffolds, unplaced_scaffolds)

## ----------------------------------------------------------------------------
## 5. Attach chromosome placement to FET data
## ----------------------------------------------------------------------------

fet_data <- fet_data %>%
  left_join(
    scaffold_map %>% select(scaffold, chrom, offset_within_chr),
    by = "scaffold"
  ) %>%
  mutate(
    pos_chr = pos + offset_within_chr
  )

## ----------------------------------------------------------------------------
## 6. Build chromosome lengths and cumulative positions
## ----------------------------------------------------------------------------

chrom_lengths <- scaffold_map %>%
  group_by(chrom) %>%
  summarise(chr_len = sum(length), .groups = "drop")

n_chrom <- sum(scaffold_lengths$type == "chrom")
chrom_levels <- c(as.character(1:n_chrom), "Unplaced")

chrom_lengths <- chrom_lengths %>%
  mutate(chrom = factor(chrom, levels = chrom_levels)) %>%
  arrange(chrom) %>%
  mutate(
    cum_start   = lag(cumsum(chr_len), default = 0),
    chrom_index = row_number()
  )

fet_data <- fet_data %>%
  mutate(chrom = factor(chrom, levels = chrom_levels)) %>%
  left_join(
    chrom_lengths %>% select(chrom, cum_start, chrom_index),
    by = "chrom"
  ) %>%
  mutate(
    pos_cum = pos_chr + cum_start,
    shade   = factor(chrom_index %% 2)
  )

## ----------------------------------------------------------------------------
## 7. Save parsed output table
## ----------------------------------------------------------------------------

write_tsv(
  fet_data %>%
    select(scaffold, pos, n_sites, freq_diff, avg_cov, logp),
  "PlainDark_perSite.parsed.tsv"
)

## ----------------------------------------------------------------------------
## 8. Axis/background helpers
## ----------------------------------------------------------------------------

axis_df <- chrom_lengths %>%
  mutate(center = cum_start + chr_len / 2)

background_rects <- chrom_lengths %>%
  mutate(
    shade = factor(chrom_index %% 2),
    xmin  = cum_start,
    xmax  = cum_start + chr_len
  )

max_logp <- max(fet_data$logp, na.rm = TRUE)

## ----------------------------------------------------------------------------
## 9. Genome-wide Manhattan plot
##    No cutoff line, no red highlighted points
## ----------------------------------------------------------------------------

p_manhattan <- ggplot() +
  geom_rect(
    data = background_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = shade),
    alpha = 0.18
  ) +
  ggrastr::geom_point_rast(
    data = fet_data,
    aes(x = pos_cum, y = logp, color = shade),
    size = 0.55,
    alpha = 0.75,
    raster.dpi = 600
  ) +
  scale_fill_manual(values = c("0" = "grey92", "1" = "grey82")) +
  scale_color_manual(values = c("0" = "grey45", "1" = "black")) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$chrom,
    expand = c(0.002, 0)
  ) +
  scale_y_continuous(
    limits = c(0, max_logp * 1.03),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line.x = element_line(linewidth = 0.4),
    axis.line.y = element_line(linewidth = 0.4)
  )

## ----------------------------------------------------------------------------
## 10. Read chr8 GFF annotation
## ----------------------------------------------------------------------------

start_region <- 1.5e6
end_region   <- 3.5e6
top <- 0
bot <- -30

ANNOT <- read.table(
  gff_file,
  sep = "\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)

names(ANNOT) <- c(
  "contig", "source", "type",
  "con_start", "con_end",
  "score", "strand", "phase", "descr"
)

ANNOT <- subset(
  ANNOT,
  contig == "NC_134752.1" &
    con_end   >= start_region &
    con_start <= end_region
)

## ----------------------------------------------------------------------------
## 11. Build gene and exon/CDS plotting data
## ----------------------------------------------------------------------------

ANNOTgene <- data.frame(
  xmin = numeric(),
  xmax = numeric(),
  ymin = numeric(),
  ymax = numeric()
)

for (g in seq_len(nrow(ANNOT))) {
  if (ANNOT$type[g] == "gene" && ANNOT$strand[g] == "-") {
    ANNOTgene[nrow(ANNOTgene) + 1, ] <- c(
      ANNOT$con_start[g],
      ANNOT$con_end[g],
      bot + (abs(top + bot)) * 0.25,
      bot + (abs(top + bot)) * 0.25
    )
  }
  if (ANNOT$type[g] == "gene" && ANNOT$strand[g] == "+") {
    ANNOTgene[nrow(ANNOTgene) + 1, ] <- c(
      ANNOT$con_start[g],
      ANNOT$con_end[g],
      bot + (abs(top + bot)) * 0.75,
      bot + (abs(top + bot)) * 0.75
    )
  }
}

colnames(ANNOTgene) <- c("xmin", "xmax", "ymin", "ymax")

ANNOTexon <- data.frame(
  xmin = numeric(),
  xmax = numeric(),
  ymin = numeric(),
  ymax = numeric()
)

for (g in seq_len(nrow(ANNOT))) {
  if (ANNOT$type[g] %in% c("exon", "CDS") && ANNOT$strand[g] == "-") {
    ANNOTexon[nrow(ANNOTexon) + 1, ] <- c(
      ANNOT$con_start[g],
      ANNOT$con_end[g],
      bot + (abs(top + bot)) * 0.10,
      bot + (abs(top + bot)) * 0.40
    )
  }
  if (ANNOT$type[g] %in% c("exon", "CDS") && ANNOT$strand[g] == "+") {
    ANNOTexon[nrow(ANNOTexon) + 1, ] <- c(
      ANNOT$con_start[g],
      ANNOT$con_end[g],
      bot + (abs(top + bot)) * 0.60,
      bot + (abs(top + bot)) * 0.90
    )
  }
}

colnames(ANNOTexon) <- c("xmin", "xmax", "ymin", "ymax")

## ----------------------------------------------------------------------------
## 12. Annotation panel
## ----------------------------------------------------------------------------

p_annot <- ggplot() +
  xlim(start_region, end_region) +
  annotate(
    "rect",
    xmin = ANNOTexon$xmin,
    xmax = ANNOTexon$xmax,
    ymin = ANNOTexon$ymin,
    ymax = ANNOTexon$ymax
  ) +
  geom_segment(
    data = ANNOTgene,
    aes(x = xmin, xend = xmax, y = ymin, yend = ymax)
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

## ----------------------------------------------------------------------------
## 13. Chr8 zoom plot with shaded GWA interval
## ----------------------------------------------------------------------------

x_min <- start_region
x_max <- end_region

gwa_start <- 2012000
gwa_end   <- 2613000

chr8_zoom <- fet_data %>%
  filter(
    chrom == "8",
    pos_chr >= x_min,
    pos_chr <= x_max
  )

p_chr8 <- ggplot(chr8_zoom, aes(x = pos_chr, y = logp)) +
  annotate(
    "rect",
    xmin = gwa_start,
    xmax = gwa_end,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey70",
    alpha = 0.25
  ) +
  ggrastr::geom_point_rast(
    color = "grey35",
    size = 0.65,
    alpha = 0.80,
    raster.dpi = 600
  ) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    labels = comma
  ) +
  scale_y_continuous(
    limits = c(0, max_logp * 1.03),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Chromosome 8 position (bp)",
    y = expression(-log[10](P))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.4),
    axis.line.y = element_line(linewidth = 0.4)
  )

## ----------------------------------------------------------------------------
## 14. Combine full figure
## ----------------------------------------------------------------------------

combined_all <- plot_grid(
  p_manhattan,
  p_annot,
  p_chr8,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(2, 1, 2)
)

## ----------------------------------------------------------------------------
## 15. Save figures
## ----------------------------------------------------------------------------

ggsave(
  filename = "PlainDark_perSite.fullGenome_chr8Genes_chr8Zoom_shadedGWA.pdf",
  plot = combined_all,
  width = 14,
  height = 12,
  units = "in",
  device = cairo_pdf
)

zoom_with_annot <- plot_grid(
  p_annot,
  p_chr8,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 3)
)

ggsave(
  filename = "PlainDark_perSite.genomewide.pdf",
  plot = p_manhattan,
  width = 14,
  height = 3,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "PlainDark_perSite.chr8_zoom_with_annotation_shadedGWA.pdf",
  plot = zoom_with_annot,
  width = 14,
  height = 3,
  units = "in",
  device = cairo_pdf
)

## ----------------------------------------------------------------------------
## Further zoomed-in DARK reference chr8 plot (non-rasterized PDF)
## Assumes fet_data already exists from the earlier dark-reference script
## ----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(scales)

dark_start <- 1750707
dark_end   <- 2893709

max_logp_zoom <- max(
  fet_data$logp[fet_data$chrom == "8" &
                  fet_data$pos_chr >= dark_start &
                  fet_data$pos_chr <= dark_end],
  na.rm = TRUE
)

dark_zoom <- fet_data %>%
  filter(
    chrom == "8",
    pos_chr >= dark_start,
    pos_chr <= dark_end
  )

p_dark_zoom <- ggplot(dark_zoom, aes(x = pos_chr, y = logp)) +
  geom_point(
    color = "grey35",
    size = 0.4,
    alpha = 0.8
  ) +
  scale_x_continuous(
    limits = c(dark_start, dark_end),
    labels = comma
  ) +
  scale_y_continuous(
    limits = c(0, max_logp_zoom * 1.03),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Chromosome 8 position (bp)",
    y = expression(-log[10](P))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.4),
    axis.line.y = element_line(linewidth = 0.4)
  )

ggsave(
  filename = "PlainDark_perSite.chr8_further_zoom_nonrasterized.pdf",
  plot = p_dark_zoom,
  width = 10,
  height = 5,
  units = "in",
  device = cairo_pdf
)

## ----------------------------------------------------------------------------
## Further zoomed-in DARK reference chr8 plot (non-rasterized PDF)
## Zoom = 2.4 Mb to 2.7 Mb
## Assumes fet_data already exists from the earlier dark-reference script
## ----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(scales)

dark_start <- 2400000
dark_end   <- 2700000

max_logp_zoom <- max(
  fet_data$logp[fet_data$chrom == "8" &
                  fet_data$pos_chr >= dark_start &
                  fet_data$pos_chr <= dark_end],
  na.rm = TRUE
)

dark_zoom <- fet_data %>%
  filter(
    chrom == "8",
    pos_chr >= dark_start,
    pos_chr <= dark_end
  )

p_dark_zoom <- ggplot(dark_zoom, aes(x = pos_chr, y = logp)) +
  geom_point(
    color = "grey35",
    size = 0.4,
    alpha = 0.8
  ) +
  scale_x_continuous(
    limits = c(dark_start, dark_end),
    labels = comma
  ) +
  scale_y_continuous(
    limits = c(0, max_logp_zoom * 1.03),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Chromosome 8 position (bp)",
    y = expression(-log[10](P))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.4),
    axis.line.y = element_line(linewidth = 0.4)
  )

ggsave(
  filename = "PlainDark_perSite.chr8_2.4Mb_2.7Mb_nonrasterized.pdf",
  plot = p_dark_zoom,
  width = 15,
  height = 3.5,
  units = "in",
  device = cairo_pdf
)

## ----------------------------------------------------------------------------
## DARK reference chr8 zoom: 2.4–2.7 Mb
## With annotation panel above
## Non-rasterized PDF
## Export dimensions matched exactly to Zoomed_in_GWAS_for_dotplot_supp.pdf
## Assumes fet_data already exists from the earlier dark-reference script
## ----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)

dark_start <- 2400000
dark_end   <- 2800000

## exact dimensions of the attached PDF
pdf_width_in  <- 654 / 72
pdf_height_in <- 265 / 72

## ----------------------------------------------------------------------------
## 1. Read / rebuild annotation for this region
## ----------------------------------------------------------------------------

# Reuse the supplied GFF path, including for the final zoom.
gff_file <- input_paths[3]

ANNOT <- read.table(
  gff_file,
  sep = "\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)

names(ANNOT) <- c(
  "contig", "source", "type",
  "con_start", "con_end",
  "score", "strand", "phase", "descr"
)

ANNOT <- subset(
  ANNOT,
  contig == "NC_134752.1" &
    con_end   >= dark_start &
    con_start <= dark_end
)

top <- 0
bot <- -30

ANNOTgene <- data.frame(
  xmin = numeric(),
  xmax = numeric(),
  ymin = numeric(),
  ymax = numeric()
)

for (g in seq_len(nrow(ANNOT))) {
  if (ANNOT$type[g] == "gene" && ANNOT$strand[g] == "-") {
    ANNOTgene[nrow(ANNOTgene) + 1, ] <- c(
      ANNOT$con_start[g],
      ANNOT$con_end[g],
      bot + (abs(top + bot)) * 0.25,
      bot + (abs(top + bot)) * 0.25
    )
  }
  if (ANNOT$type[g] == "gene" && ANNOT$strand[g] == "+") {
    ANNOTgene[nrow(ANNOTgene) + 1, ] <- c(
      ANNOT$con_start[g],
      ANNOT$con_end[g],
      bot + (abs(top + bot)) * 0.75,
      bot + (abs(top + bot)) * 0.75
    )
  }
}

colnames(ANNOTgene) <- c("xmin", "xmax", "ymin", "ymax")

ANNOTexon <- data.frame(
  xmin = numeric(),
  xmax = numeric(),
  ymin = numeric(),
  ymax = numeric()
)

for (g in seq_len(nrow(ANNOT))) {
  if (ANNOT$type[g] %in% c("exon", "CDS") && ANNOT$strand[g] == "-") {
    ANNOTexon[nrow(ANNOTexon) + 1, ] <- c(
      ANNOT$con_start[g],
      ANNOT$con_end[g],
      bot + (abs(top + bot)) * 0.10,
      bot + (abs(top + bot)) * 0.40
    )
  }
  if (ANNOT$type[g] %in% c("exon", "CDS") && ANNOT$strand[g] == "+") {
    ANNOTexon[nrow(ANNOTexon) + 1, ] <- c(
      ANNOT$con_start[g],
      ANNOT$con_end[g],
      bot + (abs(top + bot)) * 0.60,
      bot + (abs(top + bot)) * 0.90
    )
  }
}

colnames(ANNOTexon) <- c("xmin", "xmax", "ymin", "ymax")

p_annot_zoom <- ggplot() +
  xlim(dark_start, dark_end) +
  annotate(
    "rect",
    xmin = ANNOTexon$xmin,
    xmax = ANNOTexon$xmax,
    ymin = ANNOTexon$ymin,
    ymax = ANNOTexon$ymax
  ) +
  geom_segment(
    data = ANNOTgene,
    aes(x = xmin, xend = xmax, y = ymin, yend = ymax)
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

## ----------------------------------------------------------------------------
## 2. GWAS panel (non-rasterized)
## ----------------------------------------------------------------------------

dark_zoom <- fet_data %>%
  filter(
    as.character(chrom) == "8",
    pos_chr >= dark_start,
    pos_chr <= dark_end
  )

max_logp_zoom <- max(dark_zoom$logp, na.rm = TRUE)

p_dark_zoom <- ggplot(dark_zoom, aes(x = pos_chr, y = logp)) +
  geom_point(
    color = "grey35",
    size = 0.4,
    alpha = 0.8
  ) +
  scale_x_continuous(
    limits = c(dark_start, dark_end),
    labels = comma
  ) +
  scale_y_continuous(
    limits = c(0, max_logp_zoom * 1.03),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Chromosome 8 position (bp)",
    y = expression(-log[10](P))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.4),
    axis.line.y = element_line(linewidth = 0.4)
  )

## ----------------------------------------------------------------------------
## 3. Combine and save at exact same size as attached PDF
## ----------------------------------------------------------------------------

combined_zoom <- plot_grid(
  p_annot_zoom,
  p_dark_zoom,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 3)
)

ggsave(
  filename = "PlainDark_perSite.chr8_2.4Mb_2.7Mb_with_annotation_nonrasterized.pdf",
  plot = combined_zoom,
  width = pdf_width_in,
  height = pdf_height_in,
  units = "in",
  device = cairo_pdf
)
```
