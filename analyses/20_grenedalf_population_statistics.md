# Calculate Pool-seq population statistics with grenedalf

**Paper:** Figure 2B; population-statistic panels of S8–S9  

Calculate genome-wide FST, between-pool divergence (DXY) and total nucleotide diversity from the filtered Dark/Light allele-count dataset, in 10-kb windows with a 1-kb slide. Outputs use the filenames expected by analysis 10.


## Method and supporting evidence

The retained FST and pi tables support the **unbiased Hudson** implementation. Across 389,709 windows with finite comparable values, saved FST agrees with `1 - pi-within / pi-between` to within 0.0001 in every window (median absolute difference approximately 0.00000090). The alternative Nei identity has a median discrepancy of approximately 0.00734. This is evidence from the numerical outputs, not a recovered command log.

The command uses `fst --write-pi-tables`: `pi-between` supplies the DXY panel and `pi-total` supplies the total-diversity panel, matching the recovered plotting inputs. `pi-within` is an intermediate output used in calculating FST. See the [official FST documentation](https://github.com/lczech/grenedalf/wiki/Subcommand:-fst).

## Inputs

- A coordinate-sorted PoPoolation2 `.sync` or `.sync.gz` file containing the final filtered allele counts used for Pool-seq analysis, aligned to the Dark reference. Retain callable invariant positions as well as SNPs: they contribute to the per-base diversity denominator. The archived result tables contain invariant-position counts.
- The matching reference FASTA index (`.fai`), specifying chromosome names, lengths and ordering.
- The exact names identifying the Light and Dark samples in that input. A grenedalf sync header can provide these names. Without a header, names are derived from the filename and sample-column number, for example `pools.1` and `pools.2` for `pools.sync`. Pass the correct Light and Dark names explicitly. The contributed SYNC command in [analysis 14](14_poolseq_gwas.md) puts Light first and Dark second; identity with the original grenedalf input is unconfirmed. See [input formats and sample naming](https://github.com/lczech/grenedalf/wiki/Input).


## Explicit settings

| Setting | Value in this implementation | Basis |
|---|---|---|
| Estimator | `unbiased-hudson` | Supported by the saved FST/pi identities. |
| Window width / stride | 10,000 / 1,000 bp | STAR Methods and saved window coordinates. |
| Pool sizes | Light 46; Dark 50 haploid copies | Derived from 23 and 25 diploid individuals. These are diploid-locus values; the original parameter file was not recovered. |
| Window denominator | `valid-loci` | include callable invariant positions in the per-base denominator. |
| Minimum per-pool read depth | 2 | Documented filtering setting; permits the read-sampling correction. |
| Maximum per-pool read depth | 0, meaning no additional maximum | Prepared input is expected to carry the original filtering. |
| Minimum total minor-allele count | 0, meaning no additional cutoff | Prepared input is expected to carry the original filtering. |
| GWA window assignment | Fully inside chr8:2,012,000–2,613,000 | Applied by analysis 10. |
| Negative FST | Removed before plotting and summaries | Applied by analysis 10; raw grenedalf output is retained. |

 Pool sizes follow the stated numbers of diploid individuals. This remains a new implementation; the original command and parameter file were not recovered. Window averaging affects the absolute DXY and pi scale; see [windowing and averaging](https://github.com/lczech/grenedalf/wiki/Windowing). The haploid pool-size variables can be changed for loci with different ploidy.

## Software and execution

Use Bash and grenedalf. The script was checked with [grenedalf 0.6.3](https://github.com/lczech/grenedalf/releases/tag/v0.6.3); this is the validation version, not a claim about the original run.

Tool citation: [Czech, Spence and Expósito-Alonso (2024), *grenedalf: population genetic statistics for the next generation of pool sequencing*](https://doi.org/10.1093/bioinformatics/btae508), reference 107 in the manuscript.

Save the block as `run_grenedalf_statistics.sh`. Run it with five arguments:

`bash run_grenedalf_statistics.sh /path/to/pools.sync.gz /path/to/reference.fna.fai LIGHT_INPUT_NAME DARK_INPUT_NAME /path/to/new_results`

For a headerless `pools.sync.gz` with Light first and Dark second, use `pools.1 pools.2` for the two names. For a named input, use the names in its header. Do not infer sample order from the filename alone.

The output directory must be new. Set `WINDOW_AVERAGE_POLICY`, `MIN_READ_DEPTH`, `MAX_READ_DEPTH`, `SNP_MIN_COUNT`, `LIGHT_POOL_SIZE`, `DARK_POOL_SIZE` or `THREADS` before running to override the documented defaults. For example, `MIN_READ_DEPTH=10 bash run_grenedalf_statistics.sh ...` applies an additional per-pool depth threshold of 10.

The output includes `PlDkfst.csv`, `PlDkpi-between.csv`, `PlDkpi-total.csv`, the auxiliary `PlDkpi-within.csv`, `sizes.genome`, sample mappings, settings, the executed command, tool-version information and a calculation log. Coordinates are one-based and inclusive, matching analysis 10 and the [documented output format](https://github.com/lczech/grenedalf/wiki/Output).

For plotting, point analysis 10 at the new output directory and supply its separately listed GFF and windowed-Fisher input. Its existing negative-FST filter and fully-contained-window classification then apply to these results.

## Code

<a id="step-01"></a>

### 1. Calculate the two-pool window statistics

**Save as:** `run_grenedalf_statistics.sh`  

```bash
#!/usr/bin/env bash
set -euo pipefail

if (( $# != 5 )); then
    echo "Usage: bash $0 INPUT.sync[.gz] REFERENCE.fai LIGHT_NAME DARK_NAME NEW_OUTPUT_DIR" >&2
    exit 2
fi
SYNC_INPUT=$1
REFERENCE_FAI=$2
LIGHT_NAME=$3
DARK_NAME=$4
OUTDIR=$5

THREADS=${THREADS:-4}
LIGHT_POOL_SIZE=${LIGHT_POOL_SIZE:-46}
DARK_POOL_SIZE=${DARK_POOL_SIZE:-50}
WINDOW_AVERAGE_POLICY=${WINDOW_AVERAGE_POLICY:-valid-loci}
MIN_READ_DEPTH=${MIN_READ_DEPTH:-2}
MAX_READ_DEPTH=${MAX_READ_DEPTH:-0}
SNP_MIN_COUNT=${SNP_MIN_COUNT:-0}

command -v grenedalf >/dev/null || { echo "grenedalf must be on PATH." >&2; exit 1; }
[[ -s "$SYNC_INPUT" && -s "$REFERENCE_FAI" ]] || { echo "Missing sync input or reference index." >&2; exit 1; }
if [[ -z "$LIGHT_NAME" || -z "$DARK_NAME" || "$LIGHT_NAME" == "$DARK_NAME" ]]; then
    echo "Supply distinct Light and Dark sample names from the input." >&2
    exit 2
fi
for value in "$THREADS" "$LIGHT_POOL_SIZE" "$DARK_POOL_SIZE" "$MIN_READ_DEPTH" "$MAX_READ_DEPTH" "$SNP_MIN_COUNT"; do
    [[ "$value" =~ ^(0|[1-9][0-9]*)$ ]] || { echo "Counts and threads must be non-negative integers." >&2; exit 2; }
done
if (( THREADS < 1 || LIGHT_POOL_SIZE < 2 || DARK_POOL_SIZE < 2 || MIN_READ_DEPTH < 2 )); then
    echo "Use at least one thread and at least two haploid copies/reads for the sampling corrections." >&2
    exit 2
fi
if (( MAX_READ_DEPTH != 0 && MAX_READ_DEPTH < MIN_READ_DEPTH )); then
    echo "MAX_READ_DEPTH must be zero (disabled) or at least MIN_READ_DEPTH." >&2
    exit 2
fi
case "$WINDOW_AVERAGE_POLICY" in
    valid-loci|available-loci|window-length) ;;
    *) echo "Choose valid-loci, available-loci or window-length for per-base averages." >&2; exit 2 ;;
esac
[[ ! -e "$OUTDIR" ]] || { echo "Use a new output directory: $OUTDIR" >&2; exit 2; }
mkdir -p "$OUTDIR"

# Assign pool sizes by explicit sample identity, independently of input-column order.
printf '%s\tLight\n%s\tDark\n' "$LIGHT_NAME" "$DARK_NAME" > "$OUTDIR/rename_samples.tsv"
printf 'Light\t%s\nDark\t%s\n' "$LIGHT_POOL_SIZE" "$DARK_POOL_SIZE" > "$OUTDIR/pool_sizes.tsv"
cut -f 1,2 "$REFERENCE_FAI" > "$OUTDIR/sizes.genome"
grenedalf version > "$OUTDIR/grenedalf_version.txt" 2>&1

{
    printf 'setting\tvalue\n'
    printf 'input\t%s\nreference_fai\t%s\n' "$SYNC_INPUT" "$REFERENCE_FAI"
    printf 'light_input_name\t%s\ndark_input_name\t%s\n' "$LIGHT_NAME" "$DARK_NAME"
    printf 'light_haploid_copies\t%s\ndark_haploid_copies\t%s\n' "$LIGHT_POOL_SIZE" "$DARK_POOL_SIZE"
    printf 'method\tunbiased-hudson\nwindow_width_bp\t10000\nwindow_stride_bp\t1000\n'
    printf 'window_average_policy\t%s\n' "$WINDOW_AVERAGE_POLICY"
    printf 'minimum_read_depth\t%s\nmaximum_read_depth\t%s\n' "$MIN_READ_DEPTH" "$MAX_READ_DEPTH"
    printf 'total_snp_min_count\t%s\nthreads\t%s\n' "$SNP_MIN_COUNT" "$THREADS"
} > "$OUTDIR/run_settings.tsv"

cmd=(
    grenedalf fst
    --sync-path "$SYNC_INPUT"
    --reference-genome-fai "$REFERENCE_FAI"
    --rename-samples-list "$OUTDIR/rename_samples.tsv"
    --filter-samples-include Light,Dark
    --comparand Light
    --second-comparand Dark
    --pool-sizes "$OUTDIR/pool_sizes.tsv"
    --method unbiased-hudson
    --window-type interval
    --window-interval-width 10000
    --window-interval-stride 1000
    --window-average-policy "$WINDOW_AVERAGE_POLICY"
    --filter-sample-min-read-depth "$MIN_READ_DEPTH"
    --filter-sample-max-read-depth "$MAX_READ_DEPTH"
    --filter-total-snp-min-count "$SNP_MIN_COUNT"
    --write-pi-tables
    --separator-char comma
    --file-prefix PlDk
    --out-dir "$OUTDIR"
    --threads "$THREADS"
    --log-file "$OUTDIR/grenedalf_fst.log"
)
printf '%q ' "${cmd[@]}" > "$OUTDIR/command.txt"
printf '\n' >> "$OUTDIR/command.txt"
"${cmd[@]}"
```
