# Repeat library generation and annotation

**Paper:** Figure 3A/E; Figures S12–S14; deposited repeat annotations  

Build an A. gemmatalis repeat library from the Dark reference, then annotate the three representative genomes with the same library and the Dfam Obtectomera collection.

## Inputs

- Dark reference FASTA; finalized Light FASTA; full Dark01-hap2 FASTA (see analysis 02).
- EarlGrey/09_dark_genome_run/Agem_Dark_Genome_EarlGrey/Agem_Dark_Genome_strainer/Agem_Dark_Genome-families.fa.strained.FROZEN
- Dfam Obtectomera repeat collection, selected with `-r Obtectomera` in the final annotation scripts. The inspected RepeatMasker installation contains Dfam 3.9.

## Outputs

- EarlGrey/Final_Dark_Genome_Run/o/ADGO_EarlGrey/ADGO_summaryFiles/ADGO.filteredRepeats.gff
- EarlGrey/Final_Light_Genome_Run/o/ALGO_EarlGrey/ALGO_summaryFiles/ALGO.filteredRepeats.gff
- EarlGrey/Final_Dark_01_h2_haplogenome_run/o/AD1H2_EarlGrey/AD1H2_summaryFiles/AD1H2.filteredRepeats.gff

## Software

EarlGrey 7.2.1 as specified in Methods; earlGreyAnnotationOnly, RepeatModeler2, RepeatMasker, TRF and their configured dependencies. The installed FamDB metadata identifies Dfam 3.9 (10 March 2025). The final scripts and retained libraries establish selection of Obtectomera.

## Execution and interpretation

1. Use one working directory per representative genome, with G.fa as the genome FASTA and L.fa as the frozen repeat library.
2. G.fa must point to the correct representative genome and L.fa to the same frozen species library in all three runs. Copies or links to the supplied inputs can be used.
3. Representative mapping: ADGO = A/Dark reference; AD1H2 = B/Dark01 hap2; ALGO = C/Light reference. The example paths use the Final_* output folders read by the downstream code.
4. These annotation runs are the source for the downstream GFF files; earlier region tests and alternative repeat-age analyses are excluded.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Generate the species repeat library

Produces the Dark-reference EarlGrey run from which the frozen library was retained.

**Source:** `EarlGrey/run_earlgrey_dark_genome.sh`  
**Save as:** `run_earlgrey_dark_genome.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults.

```bash
#!/usr/bin/env bash
set -euo pipefail

cd /path/to/project/EarlGrey

mkdir -p 09_dark_genome_run

export EARLGREY_MP_START=fork
export TMPDIR=/tmp/${USER}/eg
export TEMP=$TMPDIR
export TMP=$TMPDIR
mkdir -p "$TMPDIR"

echo "Start: $(date)"
echo "PWD: $(pwd)"
echo "EARLGREY_MP_START=$EARLGREY_MP_START"
echo "TMPDIR=$TMPDIR"

which earlGrey
which python
which R
which jq

earlGrey \
  -g /path/to/project/Agem_genome/NCBI/GCF_050436995.1_ilAntGemm2_primary_genomic.fna \
  -s Agem_Dark_Genome \
  -o /path/to/project/EarlGrey/09_dark_genome_run \
  -t 40

echo "End: $(date)"
```

<a id="step-02"></a>

### 2. Annotate haplotype A / Dark reference

Short working directory egf; species/output prefix ADGO.

**Source:** `EarlGrey/Final_Dark_Genome_Run/run_eg_fix_obtect_short.sh`  
**Save as:** `run_eg_fix_obtect_short.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults.

```bash
#!/usr/bin/env bash
set -euo pipefail

cd /path/to/project/EarlGrey/Final_Dark_Genome_Run

export EARLGREY_MP_START=fork
export TMPDIR=/tmp/${USER}/eg
export TEMP=$TMPDIR
export TMP=$TMPDIR
mkdir -p "$TMPDIR"

GENOME=/path/to/project/EarlGrey/Final_Dark_Genome_Run/G.fa
LIB=/path/to/project/EarlGrey/Final_Dark_Genome_Run/L.fa
OUTDIR=/path/to/project/EarlGrey/Final_Dark_Genome_Run/o
SPECIES=ADGO

echo "Start: $(date)"
echo "PWD: $(pwd)"
echo "TMPDIR=$TMPDIR"
echo "EARLGREY_MP_START=$EARLGREY_MP_START"

which earlGreyAnnotationOnly
which RepeatMasker
which trf

echo
echo "=== inputs ==="
ls -lh "$GENOME" "$LIB"

echo
echo "=== path budget sanity check ==="
python - <<'PY'
p = "/path/to/project/EarlGrey/Final_Dark_Genome_Run/o/ADGO_EarlGrey/ADGO_RepeatMasker_Against_Custom_Library/RM_123456/G.fa.prep_batch-6.masked"
print("example TRF input path length:", len(p))
print(p)
PY

earlGreyAnnotationOnly \
  -g "$GENOME" \
  -s "$SPECIES" \
  -o "$OUTDIR" \
  -l "$LIB" \
  -r Obtectomera \
  -t 40 \
  -m no \
  -d no \
  -e no

echo "End: $(date)"
```

<a id="step-03"></a>

### 3. Annotate haplotype B / Dark01 hap2

Short working directory egd01h2; prefix AD1H2.

**Source:** `EarlGrey/Final_Dark_01_h2_haplogenome_run/run_dark01_h2_short.sh`  
**Save as:** `run_dark01_h2_short.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults.

```bash
#!/usr/bin/env bash
set -euo pipefail

cd /path/to/project/EarlGrey/Final_Dark_01_h2_haplogenome_run

export EARLGREY_MP_START=fork
export TMPDIR=/tmp/${USER}/eg
export TEMP=$TMPDIR
export TMP=$TMPDIR
mkdir -p "$TMPDIR"

GENOME=/path/to/project/EarlGrey/Final_Dark_01_h2_haplogenome_run/G.fa
LIB=/path/to/project/EarlGrey/Final_Dark_01_h2_haplogenome_run/L.fa
OUTDIR=/path/to/project/EarlGrey/Final_Dark_01_h2_haplogenome_run/o
SPECIES=AD1H2

echo "Start: $(date)"
echo "PWD: $(pwd)"
echo "TMPDIR=$TMPDIR"
echo "EARLGREY_MP_START=$EARLGREY_MP_START"

which earlGreyAnnotationOnly
which RepeatMasker
which trf
which python

echo
echo "=== inputs ==="
ls -lh "$GENOME" "$LIB"

echo
echo "=== path budget sanity check ==="
python - <<'PY'
p = "/path/to/project/EarlGrey/Final_Dark_01_h2_haplogenome_run/o/AD1H2_EarlGrey/AD1H2_RepeatMasker_Against_Custom_Library/RM_123456/G.fa.prep_batch-6.masked"
print("example TRF input path length:", len(p))
print(p)
PY

earlGreyAnnotationOnly \
  -g "$GENOME" \
  -s "$SPECIES" \
  -o "$OUTDIR" \
  -l "$LIB" \
  -r Obtectomera \
  -t 40 \
  -m no \
  -d no \
  -e no

echo "End: $(date)"
```

<a id="step-04"></a>

### 4. Annotate haplotype C / Light reference

Short working directory egl; prefix ALGO.

**Source:** `EarlGrey/Final_Light_Genome_Run/run_light_obtect_short.sh`  
**Save as:** `run_light_obtect_short.sh`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults.

```bash
#!/usr/bin/env bash
set -euo pipefail

cd /path/to/project/EarlGrey/Final_Light_Genome_Run

export EARLGREY_MP_START=fork
export TMPDIR=/tmp/${USER}/eg
export TEMP=$TMPDIR
export TMP=$TMPDIR
mkdir -p "$TMPDIR"

GENOME=/path/to/project/EarlGrey/Final_Light_Genome_Run/G.fa
LIB=/path/to/project/EarlGrey/Final_Light_Genome_Run/L.fa
OUTDIR=/path/to/project/EarlGrey/Final_Light_Genome_Run/o
SPECIES=ALGO

echo "Start: $(date)"
echo "PWD: $(pwd)"
echo "TMPDIR=$TMPDIR"
echo "EARLGREY_MP_START=$EARLGREY_MP_START"

which earlGreyAnnotationOnly
which RepeatMasker
which trf
which python

echo
echo "=== inputs ==="
ls -lh "$GENOME" "$LIB"

echo
echo "=== path budget sanity check ==="
python - <<'PY'
p = "/path/to/project/EarlGrey/Final_Light_Genome_Run/o/ALGO_EarlGrey/ALGO_RepeatMasker_Against_Custom_Library/RM_123456/G.fa.prep_batch-6.masked"
print("example TRF input path length:", len(p))
print(p)
PY

earlGreyAnnotationOnly \
  -g "$GENOME" \
  -s "$SPECIES" \
  -o "$OUTDIR" \
  -l "$LIB" \
  -r Obtectomera \
  -t 40 \
  -m no \
  -d no \
  -e no

echo "End: $(date)"
```
