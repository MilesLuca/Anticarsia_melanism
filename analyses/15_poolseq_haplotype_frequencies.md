# Pooled haplotype frequencies and matching genotype tracks

**Paper:** Figure 2F; Figure S11A–C  

## Inputs and software

If starting from raw PoPoolation2 counts, run section 0 first to create a new input directory containing `PlainDark_rc_NC_134752.1_withRef`, then add the other listed inputs there. Save both R code blocks under their stated filenames and the example input in section 3 as `genoplot/popmap.txt` inside that directory. The main script takes an input directory, a new output directory and the path to the pie helper. Run `Rscript pool_haplotype_frequencies.R /path/to/inputs /path/to/new_outputs /path/to/make_haplotype_pies.R`.

Place these files in the input directory, with the `genoplot` subdirectory shown below. All tabular inputs have a header, as in the supplied script.

| Input | Required contents |
|---|---|
| `PlainDark_rc_NC_134752.1_withRef` | Prepared from Dark-reference counts and the matching FASTA by [section 0](#step-00). Pool counts are expressed as numerator/denominator strings. Required R column names after `read.table` name conversion: `chr`, `pos`, `rc`, `major_alleles.maa.`, `minor_alleles.mia.`, `maa_1`, `maa_2`, `mia_1`, `mia_2`. Major/minor identity strings contain one base per pool, in order. |
| `Agem_haplotypes.DP1.merged.varsites.biallelic_positions.txt` | Haplotype-VCF positions on `NC_134752.1`, including numeric `pos`. Supply the same site list used by the collaborator. |
| `genoplot/Agem_haplotypes.DP1.merged.varsites.biallelic.vcf.gz` and its index | Biallelic haploid VCF. The identically named merged VCF is retained locally and its generation is documented in [analysis 03](03_haplotype_variants_and_genotypes.md); identity with the collaborator's copy is not yet verified. |
| `genoplot/popmap.txt` | Headered two-column sample/population map. Adapt these labels for another dataset. |
| `genoplot/Agem_popB_private_alt_sites.tsv` | B-informative site list, including `POS` (renamed to `pos` by the code). |
| `genoplot/Agem_popC_private_alt_sites.OnlyFixed.tsv` | C-informative fixed-site list, including `POS`. |

The delivered code assumes that these position tables are chromosome-8 tables and merges them by position alone. Their rows and column names must not create duplicate/suffixed frequency columns. The raw pool-count generation code is now supplied in [analysis 14](14_poolseq_gwas.md), with Light first and Dark second and the same mapping settings for either reference. This frequency workflow requires the Dark-reference run. The raw output has `N` reference bases and a `##chr` header. Section 0 now supplies the collaborator-derived reference-base preparation with corrected header handling, connecting that output to `PlainDark_rc_NC_134752.1_withRef`. The B/C-list generation and exact pool/diagnostic files were not found in the project. An existing popmap is now supplied as an example; its sample names match the retained VCF. An older single-individual private-site table was not substituted for the required group-level B/C inputs.

Software: R with dplyr, tidyr, ggplot2, scales, vcfR, GenotypePlot, cowplot, patchwork and grid; bcftools must be available on `PATH` for the filename-based GenotypePlot call. Historical package versions were not supplied. GenotypePlot's [upstream implementation](https://github.com/JimWhiting91/genotype_plot/blob/master/R/genotype_plot.R) was checked to interpret its missingness, polymorphism and position handling; that current source is not evidence of the collaborator's installed version.

## Calculations retained from the contribution

- Pool 1 is `Plain`/Light; pool 2 is Dark. Each numerator/denominator string is converted to a fraction. The code matches the reference base to each pool's major or minor base. If it matches neither, the code assigns reference fraction 0 and alternate fraction 1 for that pool.
- The alternate-frequency threshold is 0. Sites with alternate fraction 1 in both pools are removed before the position-table left join. The reference-ACGT-only filter is commented out and remains disabled.
- The supplied VCF reconstruction selects chromosome 8:1,900,000–2,810,000 inclusively, requires per-site missingness **strictly below 0.01**, then removes invariant sites. It creates `SNP_pos` before that object is used. The frequency position table uses strict outer boundaries (`pos > 1900000 & pos < 2810000`) and is restricted to the reconstructed positions.
- A left join preserves positions without usable pool fractions as `NA`; they become blank in the frequency track. Blanks indicate positions without a retained Pool-Seq allele-frequency estimate after filtering and matching to the haplotype SNP set. The supplied code does not independently test the number of pool alleles.
- B and C sites are joined to their supplied diagnostic lists. A sites are the remaining matched positions outside the union of those lists. All three groups use **`2350000 < pos < 2600000`**, a 250-kb span with strict boundaries.
- B/C distributions use alternate fractions; the A distribution uses **reference** fractions at A-group sites. Each pool/group summary is the median, excluding `NA`. Thus the A boxplot is not calculated as `1 - median(B) - median(C)`. The archived y-axis label is "Allele frequency", covering reference fractions for A and alternate fractions for B/C.
- Boxplots use the supplied `geom_boxplot` defaults; no confidence-interval calculation is supplied. The jitter has no fixed random seed, so point placement may vary between runs. GWA and ivory decoration code remains commented out.

The SNP reconstruction uses all VCF samples; the genotype plotting call uses the popmap. These must select the same samples for the reconstructed site set to match the plot filters. The section 3 example contains all 18 samples in the retained VCF, satisfying that sample-set requirement. The supplied layout creates `geno_aln` but ultimately passes the original `a` genotype plot to `plot_grid`. Those calls are retained. Check panel alignment when adapting the inputs.

## Added pie-data construction

STAR Methods specify the median frequency across private fixed diagnostic SNPs for B and C. The helper uses the B/C medians already calculated by the contributed script. For each pool it sets `B = median(B-site alternate fractions)` and `C = median(C-site alternate fractions)`.

It then sets **`A = 1 - B - C`**. 

## Archive adaptations and outputs

The source transcription is retained separately in the working review folder. Its checksum identifies that saved transcription, with whitespace normalized; it is not a hash of an original collaborator file. The main script moves the supplied `SNP_pos` reconstruction before its use, removes its redundant first VCF read, adds missing package imports, comments out two bare filename expressions, and provides explicit input/output paths. It calls the separately identified new pie helper. Scientific boundaries, remapping, missing-value handling, grouping and medians remain unchanged. Plot definitions are retained apart from the documented y-axis wording correction to "Allele frequency".

Outputs include the supplied `final_plot.pdf` and `final_plot_ribbon.pdf`, plus added exports of retained VCF sites, matched pool frequencies, A/B/C site frequencies, medians, pie proportions, diagnostic-site tracks and R session information. Separate frequency-track and boxplot PDFs are also exported for inspection.

Synthetic fixtures show identical intermediate calculation tables from the source and adapted code. They cover allele remapping, exact boundaries, missing pool rows, B/C membership, A reference fractions and medians. The new helper and the contributed box/pie/frequency components were executed. The full VCF-to-combined-figure workflow was not run: the original pool/diagnostic inputs are unavailable and the full plotting dependency set was not installed for those checks. The subsequently added example map was checked separately against the retained VCF header and the script's R reader. The reported 2,312 B sites, 2,305 C sites and biological frequencies remain unverified against the supplied code. Diagnostic-site generation remains outstanding.

## Code

<a id="step-00"></a>

### 0. Add reference bases to the pool-count table

**Save as:** `prepare_pool_reference_bases.sh`  

Run `bash prepare_pool_reference_bases.sh /path/to/Dark_reference_run/PlainDark_rc /path/to/GCF_050436995.1_ilAntGemm2_primary_genomic.fna /path/to/new_allele_frequency_inputs`. An optional fourth argument changes the chromosome; the default `NC_134752.1` matches the downstream figure script. Use the raw count table from [analysis 14, section 4](14_poolseq_gwas.md#step-04) and the same Dark-reference FASTA used for that alignment. Bash, AWK and bedtools must be available on `PATH`.

The contributed method converts one-based SNP positions to one-base BED intervals `[pos - 1, pos)`, extracts forward reference sequence with `bedtools getfasta -tab`, and places each base in column 3 (`rc`). The archive preserves all selected chromosome-8 rows, their order, allele identities, count fractions and other fields. It does not add an allele-frequency or ACGT-only filter; lowercase and ambiguous reference bases remain available for the existing R code to handle.

The archived preparation corrects three input-handling problems in the supplied shell fragment:

- Select the chromosome by an exact first-column match and preserve the table header separately. As written, the supplied `grep` removes the PoPoolation2 header, while its later `NR==1` branch and `read.table(..., h=TRUE)` still assume a header. The correction annotates and retains the first SNP as data.
- Normalize `##chr` to `chr` so the existing R reader can read the header. Preserve the header's field count: the checked two-pool raw output has 13 columns. Forcing `NF=14` would retain an added coordinate field; the archive carries any genuinely named extra input columns without adding one.
- Verify every bedtools coordinate identifier before attaching its base. Missing, skipped or out-of-range FASTA records cause an error instead of shifting subsequent annotations. The final output filename is created only after this check succeeds.

These are archive corrections to the supplied preparation fragment. They do not establish whether the collaborator performed additional manual header edits in the historical analysis. The original biological `withRef` table was not available for comparison.

The new directory contains `PlainDark_rc_NC_134752.1_withRef`, the selected chromosome count rows, one-base BED coordinates, the bedtools coordinate/base table and the normalized header. Add the position table and `genoplot` inputs listed above, then run section 1 using this directory as `INPUT_DIR`. An existing prepared `withRef` file can be used directly instead.

The preparation passed a synthetic bedtools run and the downstream R reader, covering first/last reference bases, the first SNP, exact chromosome selection, preserved field counts/order and missing-reference errors. See [validation](../VALIDATION.md).

```bash
#!/usr/bin/env bash
set -euo pipefail

if (( $# < 3 || $# > 4 )); then
    echo "Usage: bash $0 PlainDark_rc REFERENCE.fna NEW_OUTPUT_DIR [CHROMOSOME]" >&2
    exit 2
fi
COUNTS=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
REFERENCE=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
OUTDIR=$3
CHROM=${4:-NC_134752.1}
[[ "$CHROM" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "Unsupported chromosome name." >&2; exit 2; }
[[ -s "$COUNTS" && -s "$REFERENCE" ]] || { echo "Missing count table or reference FASTA." >&2; exit 1; }
command -v bedtools >/dev/null || { echo "bedtools must be on PATH." >&2; exit 1; }
[[ ! -e "$OUTDIR" ]] || { echo "Use a new output directory." >&2; exit 1; }
mkdir -p "$OUTDIR"
cd "$OUTDIR"

HEADER=pool_counts.header.tsv
SELECTED="PlainDark_${CHROM}_rc"
BED="coords.${CHROM}.bed"
BASES="PlainDark_rc.ref.${CHROM}.txt"
RESULT="PlainDark_rc_${CHROM}_withRef"

# Preserve a real header and select an exact chromosome match.
# Keep all existing fields; the checked two-pool output has 13, not 14.
awk -v chrom="$CHROM" -v header="$HEADER" -v selected="$SELECTED" -v bed="$BED" '
BEGIN { OFS="\t" }
function fail(message) { print message > "/dev/stderr"; failed=1; exit 1 }
NR==1 {
    sub(/^#+/, "", $1)
    required="chr pos rc allele_count allele_states deletion_sum snp_type major_alleles(maa) minor_alleles(mia) maa_1 maa_2 mia_1 mia_2"
    split(required, names, " ")
    if (NF < 13) fail("Expected a headered two-pool PoPoolation2 count table.")
    for (i=1; i<=13; i++) {
        if ($i != names[i]) fail("Unexpected count-table header field: " $i)
    }
    columns=NF
    print > header
    next
}
$1==chrom {
    if (NF != columns) fail("Count-table row has an inconsistent field count.")
    if ($2 !~ /^[1-9][0-9]*$/) fail("SNP positions must be positive one-based integers.")
    $1=$1
    print > selected
    printf "%s\t%d\t%d\n", $1, $2-1, $2 > bed
    kept++
}
END {
    if (!failed && !kept) {
        print "No sites found for the requested chromosome." > "/dev/stderr"
        exit 1
    }
}' "$COUNTS"

# The one-base BED interval [pos-1,pos) extracts the forward reference base.
bedtools getfasta -fi "$REFERENCE" -bed "$BED" -tab > "$BASES"

NCOLS=$(awk 'NR==1 { print NF; exit }' "$HEADER")
cat "$HEADER" > "${RESULT}.tmp"
# Check coordinate IDs before assigning bases, so skipped/out-of-range FASTA
# records cannot silently shift the remaining annotations by one row.
paste "$SELECTED" "$BASES" | awk -F '\t' -v ncols="$NCOLS" '
BEGIN { OFS="\t" }
function fail(message) { print message > "/dev/stderr"; exit 1 }
{
    if (NF != ncols+2) fail("Count rows and reference-base rows do not match.")
    expected=sprintf("%s:%d-%d", $1, $2-1, $2)
    if ($(ncols+1) != expected) fail("Reference-base coordinate mismatch: " expected)
    base=$(ncols+2)
    if (length(base) != 1 || base !~ /^[A-Za-z]$/) fail("Expected exactly one reference base.")
    $3=base
    NF=ncols
    print
}' >> "${RESULT}.tmp"
mv "${RESULT}.tmp" "$RESULT"
echo "Prepared $PWD/$RESULT"
```

<a id="step-01"></a>

### 1. Contributed frequency calculations and plotting sequence

**Save as:** `pool_haplotype_frequencies.R`  
**Dependency:** Save the helper in section 2 before running this script.

```r
# Adaptations: ordered dependencies, explicit libraries/paths, output exports,
# and a separately documented new pie-data helper. Scientific filters are retained.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript pool_haplotype_frequencies.R INPUT_DIR NEW_OUTPUT_DIR PIE_HELPER_R")
}
input_dir <- normalizePath(args[1], mustWork = TRUE)
pie_helper <- normalizePath(args[3], mustWork = TRUE)
if (file.exists(args[2])) stop("Use a new output directory.")
if (!dir.create(args[2], recursive = TRUE)) stop("Cannot create output directory.")
output_dir <- normalizePath(args[2], mustWork = TRUE)
source(pie_helper)
setwd(input_dir)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(vcfR)
library(GenotypePlot)
library(cowplot)
library(patchwork)
library(grid)

# Input paths below are relative to INPUT_DIR; preserve the collaborator's names.
# Reconstruct SNP_pos before using it to select positions for the frequency plots.
vcf_in <- read.vcfR(
  "genoplot/Agem_haplotypes.DP1.merged.varsites.biallelic.vcf.gz",
  verbose = FALSE
)

# Match genotype_plot args
target_chr <- "NC_134752.1"
start_pos <- 1900000
end_pos <- 2810000
missingness_cutoff <- 0.01

# 0. Subset to chr/start/end
vcf_in <- vcf_in[
  vcf_in@fix[, "CHROM"] == target_chr &
    as.numeric(vcf_in@fix[, "POS"]) >= start_pos &
    as.numeric(vcf_in@fix[, "POS"]) <= end_pos,
]

# 1. Apply genotype_plot missingness filter: missingness = 0.01
per_site_missing <- apply(extract.gt(vcf_in), 1, function(x) sum(is.na(x)))
per_site_missing <- per_site_missing / (ncol(vcf_in@gt) - 1)

keep_missing <- which(per_site_missing < missingness_cutoff)
vcf_in2 <- vcf_in[keep_missing, ]

# 2. Apply genotype_plot invariant_filter = TRUE
vcf_in3 <- vcf_in2[is.polymorphic(vcf_in2, na.omit = TRUE), ]

# 3. Recreate internal SNP_pos table
SNP_pos <- data.frame(
  chrom = unique(vcf_in3@fix[, "CHROM"]),
  BP = as.integer(vcf_in3@fix[, "POS"])
)

SNP_pos$GEN_pos <- seq(
  min(SNP_pos$BP),
  max(SNP_pos$BP),
  by = (max(SNP_pos$BP) - min(SNP_pos$BP)) / nrow(SNP_pos)
)[1:nrow(SNP_pos)]

SNP_pos <- SNP_pos %>%
  mutate(geno_order = row_number())

## code for making the frequency plot to match the genoplot


chr_use    <- "NC_134752.1"
start_pos  <- 1900000
end_pos    <- 2810000

gwas_start <- 2012000
gwas_end   <- 2613000

ivory_start <- 2557339
ivory_end   <- 2766471

alt_threshold <- 0

tmp <- read.table('PlainDark_rc_NC_134752.1_withRef',h=T)
tmp$thingy <- NULL

convert_frac <- function(x) {
  sapply(strsplit(x, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
}

tmp$maa_1 <- convert_frac(tmp$maa_1)
tmp$maa_2 <- convert_frac(tmp$maa_2)
tmp$mia_1 <- convert_frac(tmp$mia_1)
tmp$mia_2 <- convert_frac(tmp$mia_2)

plot_dat <- tmp %>%
  mutate(
    rc = toupper(rc),
    # population-specific major/minor allele identities
    maj1 = substr(major_alleles.maa., 1, 1),
    maj2 = substr(major_alleles.maa., 2, 2),
    min1 = substr(minor_alleles.mia., 1, 1),
    min2 = substr(minor_alleles.mia., 2, 2)
  ) %>%
  filter(
    chr == chr_use,
    pos >= start_pos,
    pos <= end_pos
#    rc %in% c("A", "C", "G", "T")   # drop rows where reference is N
  ) %>%
  mutate(
    # if rc matches neither observed allele, treat as alt = 1
    ref_1 = case_when(
      rc == maj1 ~ maa_1,
      rc == min1 ~ mia_1,
      TRUE ~ 0
    ),
    alt_1 = case_when(
      rc == maj1 ~ mia_1,
      rc == min1 ~ maa_1,
      TRUE ~ 1
    ),
    ref_2 = case_when(
      rc == maj2 ~ maa_2,
      rc == min2 ~ mia_2,
      TRUE ~ 0
    ),
    alt_2 = case_when(
      rc == maj2 ~ mia_2,
      rc == min2 ~ maa_2,
      TRUE ~ 1
    )
  ) %>%
  # drop low-alt sites after remapping, then remove sites where both pops are alt == 1
  filter(
    (alt_1 >= alt_threshold | alt_2 >= alt_threshold) &
    !(alt_1 == 1 & alt_2 == 1)
  ) %>%
  arrange(pos) %>%
  mutate(
    variant_index = row_number(),
    in_gwas = pos >= gwas_start & pos <= gwas_end,
    in_ivory = pos >= ivory_start & pos <= ivory_end
  )

##first import the list of positions luca used in the genoplot for the main figure
vcf_positions <- read.table('Agem_haplotypes.DP1.merged.varsites.biallelic_positions.txt',h=T)
retained_pos <- SNP_pos$BP

vcf_positions <- subset(vcf_positions, pos > start_pos & pos < end_pos & pos %in% retained_pos)
### merge it with plot_dat
freq_sites_matching_vcf <- merge(vcf_positions, plot_dat, by='pos',all.x=TRUE)
### crop to interval used in that plot
#freq_sites_matching_vcf <- freq_sites_matching_vcf[freq_sites_matching_vcf$pos > 1900000 & freq_sites_matching_vcf$pos < 2810000,]

#plot_dat <- subset(freq_sites_matching_vcf, pos %in% vcf_positions$pos)

freq_sites_matching_vcf <- freq_sites_matching_vcf %>%
  arrange(pos) %>%
  mutate(
    variant_index = row_number(),
    in_gwas = pos >= gwas_start & pos <= gwas_end,
    in_ivory = pos >= ivory_start & pos <= ivory_end
  )
## switch NAs for Ns

# Re-indexed GWAS interval
gwas_idx <- freq_sites_matching_vcf %>%
  filter(in_gwas) %>%
  summarise(
    xmin = min(variant_index),
    xmax = max(variant_index)
  )

# Re-indexed ivory interval
ivory_idx <- freq_sites_matching_vcf %>%
  filter(in_ivory) %>%
  summarise(
    xmin = min(variant_index),
    xmax = max(variant_index)
  )

# then resume plot
plot_long <- freq_sites_matching_vcf %>%
  transmute(
    variant_index,
    genomic_pos = pos,
    pop1_ref = ref_1,
    pop1_alt = alt_1,
    pop2_ref = ref_2,
    pop2_alt = alt_2
  ) %>%
  pivot_longer(
    cols = c(pop1_ref, pop1_alt, pop2_ref, pop2_alt),
    names_to = c("population", "allele_class"),
    names_sep = "_",
    values_to = "fraction"
  ) %>%
  mutate(
    population = recode(population,
      pop1 = "Plain",
      pop2 = "Dark"
    ),
    allele_class = recode(allele_class,
      ref = "Reference",
      alt = "Alternate"
    )
  )

p <- ggplot(plot_long, aes(x = variant_index, y = fraction, fill = allele_class)) +
  geom_col(width = 1) +
  facet_wrap(~population, ncol = 1) +
  scale_fill_manual(
    values = c(
      Reference = "#134c84",
      Alternate = "#f89c30"
    )
  ) +
  scale_y_continuous(
    limits = c(-0.12, 1.05),
    expand = c(0, 0)
  ) +
  labs(
    x = "Variant site index",
    y = "Allele fraction",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# Add GWAS interval marking if any retained sites remain in that interval
#if (nrow(gwas_idx) > 0 && is.finite(gwas_idx$xmin) && is.finite(gwas_idx$xmax)) {
#  p <- p +
#    annotate(
#      "rect",
#      xmin = gwas_idx$xmin - 0.5,
#      xmax = gwas_idx$xmax + 0.5,
#      ymin = 0,
#      ymax = 1.05,
#      fill = "grey70",
#      alpha = 0.18
#    ) +
#    geom_vline(
#      xintercept = c(gwas_idx$xmin - 0.5, gwas_idx$xmax + 0.5),
#      linetype = "dashed",
#      linewidth = 0.5,
#      color = "black"
#    ) +
#    annotate(
#      "text",
#      x = mean(c(gwas_idx$xmin, gwas_idx$xmax)),
#      y = 1.03,
#      label = "GWAS interval",
#      size = 4,
#      fontface = "bold"
#    )
#}

# Add ivory annotation track below each panel if any retained sites remain in that interval
#if (nrow(ivory_idx) > 0 && is.finite(ivory_idx$xmin) && is.finite(ivory_idx$xmax)) {
#  p <- p +
#    annotate(
#      "rect",
#      xmin = ivory_idx$xmin - 0.5,
#      xmax = ivory_idx$xmax + 0.5,
#      ymin = -0.09,
#      ymax = -0.03,
#      fill = "#E8E1C8",
#      color = NA
#    ) +
#    annotate(
#      "text",
#      x = mean(c(ivory_idx$xmin, ivory_idx$xmax)),
#      y = -0.105,
#      label = "ivory",
#      size = 3.5,
#      fontface = "bold"
#    )
#}

p

#### the genoplot code
#### using the old genome
popmap <- read.table('genoplot/popmap.txt', h=T)
new_plot <- genotype_plot(
  vcf = "genoplot/Agem_haplotypes.DP1.merged.varsites.biallelic.vcf.gz",
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

combine_genotype_plot(new_plot)
a <- new_plot$genotypes

### makin it line up

# Make both plots use identical x limits and no x expansion
xlim_sites <- range(freq_sites_matching_vcf$variant_index, na.rm = TRUE)

p_aln <- p +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    # remove axes
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    # remove facet labels Plain / Dark
    strip.text = element_blank(),
    strip.background = element_blank(),
    # bring facet panels close together
    panel.spacing = unit(0, "lines"),
    # keep left/right margins, shrink top/bottom
    plot.margin = margin(0, 5.5, 0, 5.5),
    legend.position = "none"
    )
p_aln

geno_aln <- a +
  coord_cartesian(xlim = xlim_sites, clip = "off") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    plot.margin = margin(2, 5.5, 2, 5.5)
  )

aligned <- align_plots(
  p_aln,
  geno_aln,
  align = "v",
  axis = "lr"
)

plot_grid(
  aligned[[1]],
  a,
  ncol = 1,
  rel_heights = c(0.5, 2),  # make top plot shorter → less visual gap
  align = "v",
  axis = "lr"
)

# further frequency calculatons

# Agem_popB_private_alt_sites.tsv
# Agem_popC_private_alt_sites.OnlyFixed.tsv

bAlleles <- read.table('genoplot/Agem_popB_private_alt_sites.tsv',h=T)
cAlleles <- read.table('genoplot/Agem_popC_private_alt_sites.OnlyFixed.tsv', h=T)

# standardise column names
names(bAlleles)[names(bAlleles) == "POS"] <- "pos"
names(cAlleles)[names(cAlleles) == "POS"] <- "pos"

# merge + filter
bAlleles <- bAlleles %>%
  merge(freq_sites_matching_vcf, by = "pos") %>%
  filter(pos < 2600000 & pos > 2350000)

cAlleles <- cAlleles %>%
  merge(freq_sites_matching_vcf, by = "pos") %>%
  filter(pos < 2600000 & pos > 2350000)

## get A positions
all_sites <- freq_sites_matching_vcf %>%
  filter(pos < 2600000 & pos > 2350000)
b_pos <- bAlleles %>% distinct(pos)
c_pos <- cAlleles %>% distinct(pos)
# define A as everything not in B or C
aAlleles <- all_sites %>%
  anti_join(bind_rows(b_pos, c_pos), by = "pos")


abc_cols <- c(
  A = "#514235",
  B = "#3191bc",
  C = "#9e8872"
)

plot_df <- bind_rows(
  aAlleles %>% transmute(group = "A", pos, alt_1 = ref_1, alt_2 = ref_2),
  bAlleles %>% transmute(group = "B", pos, alt_1, alt_2),
  cAlleles %>% transmute(group = "C", pos, alt_1, alt_2)
) %>%
  pivot_longer(
    cols = c(alt_1, alt_2),
    names_to = "population",
    values_to = "alt_freq"
  ) %>%
  mutate(
    population = recode(population,
      alt_1 = "Light",
      alt_2 = "Dark"
    ),
    population = factor(population, levels = c("Dark", "Light")),
    group = factor(group, levels = c("A", "B", "C")),
    x_group = interaction(population, group, sep = "_"),
    x_group = factor(x_group, levels = c("Dark_A","Dark_B", "Dark_C", "Light_A","Light_B", "Light_C"))
  )

med_df <- plot_df %>%
  group_by(population, group, x_group) %>%
  summarise(
    median_alt = median(alt_freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = factor(group, levels = names(abc_cols)),
    label_hjust = case_when(
      group == "A" ~ 1.15,
      group == "B" ~ 1.15,
      TRUE ~ -0.15
    )
  )

box_p <- ggplot(plot_df, aes(x = x_group, y = alt_freq, fill = group)) +
  geom_jitter(
    aes(color = group),
    width = 0.08,
    size = 1.2,
    alpha = 0.08,
    stroke = 0
  ) +
  geom_boxplot(
    width = 0.22,
    linewidth = 0.6,
    outlier.shape = NA,
    color = "black",
    alpha = 0.95
  ) +
  geom_text(
    data = med_df,
    aes(
      x = x_group,
      y = median_alt,
      label = sprintf("%.2f", median_alt),
      hjust = label_hjust
    ),
    size = 3.4,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = abc_cols,
    breaks = c("A", "B", "C"),
    limits = c("A", "B", "C")
  ) +
  scale_color_manual(
    values = abc_cols,
    breaks = c("A", "B", "C"),
    limits = c("A", "B", "C")
  ) +
  scale_x_discrete(
    labels = c(
      Dark_A  = "A\nDark",
      Dark_B  = "B",
      Dark_C  = "C",
      Light_A = "A\nLight",
      Light_B = "B",
      Light_C = "C"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Allele frequency",
    fill = NULL,
    color = NULL
  ) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold", lineheight = 0.9),
    axis.line = element_line(linewidth = 0.55),
    axis.ticks = element_line(linewidth = 0.55),
    plot.margin = margin(3, 4, 3, 3)
  )

pie_df <- make_haplotype_pie_data(med_df)

pie_p <- ggplot(pie_df, aes(x = "", y = fraction, fill = component)) +
         geom_col(width = 1, color = "white", linewidth = 0.35) +
         geom_text( aes(y = label_pos, label = label), x = 1, fontface = "bold", size = 3.6, color = "white" ) +
         coord_polar(theta = "y") +
         facet_wrap(~population, nrow = 1) +
         scale_fill_manual( values = abc_cols, breaks = names(abc_cols), limits = names(abc_cols), drop = FALSE ) +
         theme_void(base_size = 13) + theme( legend.position = "none", strip.text = element_text(face = "bold", margin = margin(0, 0, 1, 0)), panel.spacing = unit(0.2, "lines"), plot.margin = margin(2, 2, 2, 2) )

final_p <- box_p / pie_p + plot_layout(heights = c(2.1, 0.9)) & theme(plot.margin = margin(2, 2, 2, 2))
final_p

ggsave(file.path(output_dir, "final_plot.pdf"), plot = final_p, width = 6, height = 4, units = "in")
### now make ribbon plots of used sites.

# build plotting data
ribbon_df <- bind_rows(
  aAlleles %>% transmute(variant_index, track = "A"),
  bAlleles %>% transmute(variant_index, track = "B"),
  cAlleles %>% transmute(variant_index, track = "C")
)

ribbon_p <- ggplot(ribbon_df) +
  geom_rect(
    aes(
      xmin = variant_index - 0.4,
      xmax = variant_index + 0.4,
      ymin = 0,
      ymax = 1,
      fill = track
    ),
    color = NA
  ) +
  facet_wrap(~track, ncol = 1) +
  scale_fill_manual(values = c(A = "#514235", B = "#3191bc", C = "#9e8872")) +
  scale_x_continuous(
    limits = xlim_sites,
    expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    strip.text = element_blank(),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 5.5, 0, 5.5)
  )

final_ribbon <- plot_grid(
  aligned[[1]],
  a,
  ribbon_p,
  ncol = 1,
  rel_heights = c(0.5, 2, 0.5),  # make top plot shorter → less visual gap
  align = "v",
  axis = "lr")

ggsave(file.path(output_dir, "final_plot_ribbon.pdf"), plot = final_ribbon, width = 6, height = 4, units = "in")


# Archive additions: export intermediate calculations and separate plot components.
write.table(SNP_pos, file.path(output_dir, "retained_vcf_sites.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(freq_sites_matching_vcf, file.path(output_dir, "matched_pool_frequencies.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(plot_df, file.path(output_dir, "ABC_site_frequencies.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(med_df, file.path(output_dir, "ABC_site_medians.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pie_df, file.path(output_dir, "ABC_pie_proportions.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)
write.table(ribbon_df, file.path(output_dir, "ABC_site_tracks.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
ggsave(file.path(output_dir, "allele_frequency_tracks.pdf"), p, width = 10, height = 3)
ggsave(file.path(output_dir, "ABC_frequency_boxplots.pdf"), box_p, width = 6, height = 3)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
```

<a id="step-02"></a>

### 2. New pie-data helper

**Save as:** `make_haplotype_pies.R`  

```r
# New archive implementation: B/C medians follow STAR Methods.
# This differs from the A-site reference-frequency distribution in the boxplot.
make_haplotype_pie_data <- function(med_df) {
  required <- c("population", "group", "median_alt")
  if (!all(required %in% names(med_df))) stop("Missing median-summary columns.")
  rows <- lapply(c("Dark", "Light"), function(pop) {
    x <- med_df[as.character(med_df$population) == pop &
                  as.character(med_df$group) %in% c("B", "C"), ]
    if (nrow(x) != 2L || anyDuplicated(as.character(x$group))) {
      stop(pop, ": require one median for each of B and C.")
    }
    bc <- x$median_alt[match(c("B", "C"), as.character(x$group))]
    if (!is.numeric(bc) || any(!is.finite(bc)) || any(bc < 0 | bc > 1)) {
      stop(pop, ": B/C medians must be finite fractions between 0 and 1.")
    }
    if (sum(bc) > 1) stop(pop, ": B + C exceeds 1; review site selection and inputs.")
    fractions <- c(1 - sum(bc), bc)
    data.frame(population = pop, component = c("A", "B", "C"),
               fraction = fractions,
               # ggplot's default stacking puts the last factor level at the bottom.
               label_pos = 1 - cumsum(fractions) + fractions / 2,
               label = ifelse(fractions == 0, "",
                              paste0(c("A", "B", "C"), "\n",
                                     sprintf("%.0f%%", 100 * fractions))))
  })
  result <- do.call(rbind, rows)
  result$population <- factor(result$population, levels = c("Dark", "Light"))
  result$component <- factor(result$component, levels = c("A", "B", "C"))
  result
}
```

<a id="step-03"></a>

### 3. Example genotype-plot sample map

**Save as:** `genoplot/popmap.txt` inside `INPUT_DIR`  
**Source:** `Pacbio_Reseq_Hilo/ivory_haplotypes/04_genotype_plot/Agem_haplotypes.popmap.byhaplo.tsv`  

Create the `genoplot` input subdirectory if needed and save this tab-separated block there. It preserves all 18 source records and their order, adding only the `ind`/`pop` header required by the contributed `read.table(..., h=TRUE)` call. Both columns are read as sample identifiers and plotting-group labels by the existing workflow. The names match all 18 samples in the retained merged VCF; the R reader was checked to retain the first sample and every subsequent record.

This example labels individual haplotypes and retains the original 01/02/05/06 identifiers. See the [publication-name mapping](../README.md#sample-names). It is sufficient for sharing the workflow, and is not presented as the final S11B A/B/C grouping or order. When adapting to another VCF, replace the sample names/labels consistently and ensure that the map and SNP reconstruction use the same sample set.

```text
ind	pop
Agem_Reference_reassembly.hap1_cortex_contig	Ref_h1
Agem_Reference_reassembly.hap2_cortex_contig	Ref_h2
Agem_Dark_06.hap1_cortex_contig	Dark_06_h1
Agem_Dark_05.hap1_cortex_contig	Dark_05_h1
Agem_Dark_02.hap1_cortex_contig	Dark_02_h1
Agem_Dark_01.hap1_cortex_contig	Dark_01_h1
Agem_Dark_01.hap2_cortex_contig	Dark_01_h2
Agem_Dark_05.hap2_cortex_contig	Dark_05_h2
Agem_Dark_06.hap2_cortex_contig	Dark_06_h2
Agem_Dark_02.hap2_cortex_contig	Dark_02_h2
Agem_Light_02.hap1_cortex_contig	Light_02_h1
Agem_Light_05.hap1_cortex_contig	Light_05_h1
Agem_Light_06.hap1_cortex_contig	Light_06_h1
Agem_Light_05.hap2_cortex_contig	Light_05_h2
Agem_Light_06.hap2_cortex_contig	Light_06_h2
Agem_Light_01.hap1_cortex_contig	Light_01_h1
Agem_Light_01.hap2_cortex_contig	Light_01_h2
Agem_Light_02.hap2_cortex_contig	Light_02_h2
```
