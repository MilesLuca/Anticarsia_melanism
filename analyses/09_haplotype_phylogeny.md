# Haplotype phylogeny and occupancy sensitivity

**Paper:** Figure S14C–D  

Align homologous cortex-region sequences from A/B/C and three erebid outgroups, test alignment occupancy filters and infer maximum-likelihood trees.

## Inputs

- Haplotype_trees/Cortex_region_fastas.fasta (six sequences: A/B/C, Catocala fraxini, Dysgonia rogenhoferi, Eilema depressum)

The six cortex-region sequences were extracted and combined in Geneious, then exported as `Cortex_region_fastas.fasta`. The code starts from this prepared FASTA. Haplotype B is subsequently reverse-complemented by the included script before MAFFT alignment.

## Outputs

- B-reverse-complemented MAFFT alignment; five occupancy-filtered alignments.
- Haplotype_trees/04_trees/trim_*.treefile and associated logs
- Haplotype_trees/05_plots/trim_sensitivity_trees.pdf

## Software

MAFFT, Python 3, IQ-TREE and R/ape. Matching saved final and sensitivity runs report IQ-TREE 2.3.6.

## Execution and interpretation

1. Run from your Haplotype_trees working directory; create 01_alignment before the header-cleaning step and 05_plots before plotting.
2. Only the header-cleaning portion of the initial low-memory alignment script is needed as an input-preparation step. The final alignment is made by 03_revcomp_B_and_mafft.slurm.
3. The final alignment uses MAFFT --retree 2 --maxiterate 0 --thread 1 with B reverse-complemented. Occupancy filters treat -, N and ? as absent.
4. The strict all-six-taxa alignment has 5,459 columns. Other recorded lengths: A/B/C+1 outgroup 8,429; A/B/C+2 outgroups 7,833; at least 5/6 11,038; at least 4/6 18,021.
5. IQ-TREE settings are -m MFP -B 1000 -alrt 1000 -T 1. The strict tree’s reported 44.7/52 support is present in the saved result. Model/seed provenance is retained in the original logs; a fresh stochastic run need not give identical support.
6. The +2-outgroup run and the sensitivity loop generate complementary outputs; run both before plotting. The strict tree is the published S14C panel, despite earlier comments calling +2 outgroups the main tree.
7. Geneious input preparation is required.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Prepare cleaned input headers

**Source:** `Haplotype_trees/scripts/01_mafft_align_lowmem.slurm`  
**Save as:** `prepare_cortex_headers.sh`  
**Archive treatment:** Verbatim source lines 32–40: header cleaning only. Run from Haplotype_trees with 01_alignment already created.

```bash
awk '
  /^>/ {
    gsub(/[[:space:]]+$/, "", $0)
    gsub(/[[:space:]]+/, "_", $0)
    print
    next
  }
  { print }
' Cortex_region_fastas.fasta > 01_alignment/Cortex_region_fastas.clean.fasta
```

<a id="step-02"></a>

### 2. Orient B and run the final MAFFT alignment

Consumes the cleaned six-sequence FASTA.

**Source:** `Haplotype_trees/scripts/03_revcomp_B_and_mafft.slurm`  
**Save as:** `03_revcomp_B_and_mafft.slurm`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
set -euo pipefail

cd /path/to/project/Haplotype_trees

mkdir -p 01_alignment 00_logs


echo "Running on: $(hostname)"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo "MAFFT path: $(command -v mafft || echo 'MAFFT_NOT_FOUND')"
echo "MAFFT version:"
mafft --version || true

# Make a new FASTA where only Haplotype_B is reverse-complemented.
python - <<'PY'
from pathlib import Path

inp = Path("01_alignment/Cortex_region_fastas.clean.fasta")
out = Path("01_alignment/Cortex_region_fastas.B_revcomp.clean.fasta")

comp = str.maketrans(
    "ACGTNacgtn",
    "TGCANtgcan"
)

records = []
name = None
seq = []

with inp.open() as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name = line[1:]
            seq = []
        else:
            seq.append(line)
    if name is not None:
        records.append((name, "".join(seq)))

with out.open("w") as oh:
    for name, s in records:
        if name == "Haplotype_B_cortex_region":
            s = s.translate(comp)[::-1]
            name = name + "_revcomp"
        oh.write(f">{name}\n")
        s = s.upper()
        for i in range(0, len(s), 80):
            oh.write(s[i:i+80] + "\n")

print(f"Wrote {out}")
for name, s in records:
    print(f"{name}\t{len(s)}")
PY

rm -f 01_alignment/Cortex_region_fastas.B_revcomp.mafft.lowmem.fa
rm -f 01_alignment/Cortex_region_fastas.B_revcomp.mafft.lowmem.stderr
rm -f 01_alignment/alignment_sequence_stats.B_revcomp.lowmem.tsv
rm -f 01_alignment/aligned_headers.B_revcomp.lowmem.txt

# Low-memory MAFFT alignment, single-threaded.
mafft \
  --thread 1 \
  --retree 2 \
  --maxiterate 0 \
  01_alignment/Cortex_region_fastas.B_revcomp.clean.fasta \
  > 01_alignment/Cortex_region_fastas.B_revcomp.mafft.lowmem.fa \
  2> 01_alignment/Cortex_region_fastas.B_revcomp.mafft.lowmem.stderr

# Alignment diagnostics.
python - <<'PY'
from pathlib import Path

fasta = Path("01_alignment/Cortex_region_fastas.B_revcomp.mafft.lowmem.fa")
out = Path("01_alignment/alignment_sequence_stats.B_revcomp.lowmem.tsv")

records = []
name = None
seq = []

with fasta.open() as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name = line[1:]
            seq = []
        else:
            seq.append(line)
    if name is not None:
        records.append((name, "".join(seq)))

with out.open("w") as oh:
    oh.write("name\taligned_length\tgap_count\tgap_fraction\tN_count\n")
    for name, s in records:
        gaps = s.count("-")
        n_count = s.upper().count("N")
        gap_fraction = gaps / len(s) if s else 0
        oh.write(f"{name}\t{len(s)}\t{gaps}\t{gap_fraction:.5f}\t{n_count}\n")

print(f"Aligned records: {len(records)}")
if records:
    lengths = sorted(set(len(s) for _, s in records))
    print(f"Unique aligned lengths: {lengths}")

for name, s in records:
    gaps = s.count("-")
    print(f"{name}\taligned={len(s)}\tgaps={gaps}\tgap_fraction={gaps/len(s):.5f}")
PY

grep "^>" 01_alignment/Cortex_region_fastas.B_revcomp.mafft.lowmem.fa \
  > 01_alignment/aligned_headers.B_revcomp.lowmem.txt

echo "Finished: $(date)"
echo "Output alignment: 01_alignment/Cortex_region_fastas.B_revcomp.mafft.lowmem.fa"
```

<a id="step-03"></a>

### 3. Create occupancy-filtered alignments

Produces all five sensitivity alignments described above.

**Source:** `Haplotype_trees/scripts/06_make_candidate_trimmed_alignments.py`  
**Save as:** `06_make_candidate_trimmed_alignments.py`  
**Archive treatment:** Complete source, verbatim.

```python
from pathlib import Path

aln = Path("01_alignment/Cortex_region_fastas.B_revcomp.mafft.lowmem.fa")
outdir = Path("03_trimmed")
outdir.mkdir(exist_ok=True)

records = []
name = None
seq = []

with aln.open() as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq).upper()))
            name = line[1:]
            seq = []
        else:
            seq.append(line)
    if name is not None:
        records.append((name, "".join(seq).upper()))

names = [n for n, s in records]
seqs = dict(records)
L = len(next(iter(seqs.values())))

haplotypes = [
    "Haplotype_A_cortex_region",
    "Haplotype_B_cortex_region_revcomp",
    "Haplotype_C_cortex_region",
]

outgroups = [
    "C_fra_cortex_region",
    "D_rog_cortex_region",
    "E_dep_cortex_region",
]

def present(c):
    return c not in "-N?"

def write_alignment(label, keep_cols):
    out = outdir / f"{label}.fa"

    with out.open("w") as oh:
        for name in names:
            s = seqs[name]
            trimmed = "".join(s[i] for i in keep_cols)
            oh.write(f">{name}\n")
            for j in range(0, len(trimmed), 80):
                oh.write(trimmed[j:j+80] + "\n")

def summarize(label, keep_cols):
    ncols = len(keep_cols)
    variable = 0
    informative = 0

    for i in keep_cols:
        chars = [seqs[n][i] for n in names if present(seqs[n][i])]
        alleles = set(chars)

        if len(alleles) > 1:
            variable += 1

        counts = {a: chars.count(a) for a in alleles}
        if sum(v >= 2 for v in counts.values()) >= 2:
            informative += 1

    print(f"{label}\tcolumns={ncols}\tvariable={variable}\tparsimony_informative={informative}")

    for name in names:
        s = seqs[name]
        chars = [s[i] for i in keep_cols]
        gaps = chars.count("-")
        called = sum(present(c) for c in chars)
        gap_fraction = gaps / ncols if ncols else 0
        print(f"  {name}\tcalled={called}\tgaps={gaps}\tgap_fraction={gap_fraction:.4f}")

# Candidate 1: all six taxa present
all6 = [
    i for i in range(L)
    if all(present(seqs[n][i]) for n in names)
]

# Candidate 2: all three haplotypes present, and at least one outgroup present
haps_plus_1out = [
    i for i in range(L)
    if all(present(seqs[n][i]) for n in haplotypes)
    and any(present(seqs[n][i]) for n in outgroups)
]

# Candidate 3: all three haplotypes present, and at least two outgroups present
haps_plus_2out = [
    i for i in range(L)
    if all(present(seqs[n][i]) for n in haplotypes)
    and sum(present(seqs[n][i]) for n in outgroups) >= 2
]

# Candidate 4: at least five of six taxa present
atleast5of6 = [
    i for i in range(L)
    if sum(present(seqs[n][i]) for n in names) >= 5
]

# Candidate 5: at least four of six taxa present
atleast4of6 = [
    i for i in range(L)
    if sum(present(seqs[n][i]) for n in names) >= 4
]

candidates = {
    "trim_all6_present": all6,
    "trim_haplotypes_plus_1outgroup": haps_plus_1out,
    "trim_haplotypes_plus_2outgroups": haps_plus_2out,
    "trim_atleast5of6": atleast5of6,
    "trim_atleast4of6": atleast4of6,
}

print("source_alignment", aln)
print("source_length", L)
print("\nCandidate trimmed alignments:")

for label, cols in candidates.items():
    write_alignment(label, cols)
    summarize(label, cols)
    print()
```

<a id="step-04"></a>

### 4. Infer the A/B/C plus two-outgroup tree

Complements the four trees generated in the next step.

**Source:** `Haplotype_trees/scripts/07_iqtree_haplotypes_plus_2outgroups.slurm`  
**Save as:** `07_iqtree_haplotypes_plus_2outgroups.slurm`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
set -euo pipefail

cd /path/to/project/Haplotype_trees

mkdir -p 04_trees 00_logs


echo "Running on: $(hostname)"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo "iqtree2 path: $(command -v iqtree2 || echo 'IQTREE2_NOT_FOUND')"
echo "iqtree path: $(command -v iqtree || echo 'IQTREE_NOT_FOUND')"

ALIGNMENT="03_trimmed/trim_haplotypes_plus_2outgroups.fa"
PREFIX="04_trees/trim_haplotypes_plus_2outgroups"

if command -v iqtree2 >/dev/null 2>&1; then
    IQTREE=iqtree2
elif command -v iqtree >/dev/null 2>&1; then
    IQTREE=iqtree
else
    echo "ERROR: neither iqtree2 nor iqtree found in PATH" >&2
    exit 1
fi

echo "Using: $IQTREE"
$IQTREE --version || true

# Remove old outputs with this prefix, if present
rm -f ${PREFIX}.*

# First tree: model selection + ultrafast bootstrap + SH-aLRT
# Single-threaded by request.
$IQTREE \
  -s "$ALIGNMENT" \
  -st DNA \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -T 1 \
  --prefix "$PREFIX"

echo "Finished: $(date)"
echo "Treefile: ${PREFIX}.treefile"
echo "Logfile: ${PREFIX}.log"
```

<a id="step-05"></a>

### 5. Infer the remaining sensitivity trees

Includes the strict all-six-taxa tree used in S14C.

**Source:** `Haplotype_trees/scripts/08_iqtree_trim_sensitivity.slurm`  
**Save as:** `08_iqtree_trim_sensitivity.slurm`  
**Archive treatment:** Complete analysis source with environment settings adapted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Remove scheduler directives and job logging; use ordinary Bash and explicit thread defaults. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash
#!/usr/bin/env bash
set -euo pipefail

cd /path/to/project/Haplotype_trees

mkdir -p 04_trees 00_logs


if command -v iqtree2 >/dev/null 2>&1; then
    IQTREE=iqtree2
elif command -v iqtree >/dev/null 2>&1; then
    IQTREE=iqtree
else
    echo "ERROR: neither iqtree2 nor iqtree found in PATH" >&2
    exit 1
fi

echo "Running on: $(hostname)"
echo "Started: $(date)"
echo "Using: $IQTREE"
$IQTREE --version || true

for ALIGNMENT in \
    03_trimmed/trim_all6_present.fa \
    03_trimmed/trim_haplotypes_plus_1outgroup.fa \
    03_trimmed/trim_atleast5of6.fa \
    03_trimmed/trim_atleast4of6.fa
do
    BASENAME=$(basename "$ALIGNMENT" .fa)
    PREFIX="04_trees/${BASENAME}"

    echo
    echo "========================================"
    echo "Running IQ-TREE on $ALIGNMENT"
    echo "Prefix: $PREFIX"
    echo "========================================"

    rm -f ${PREFIX}.*

    $IQTREE \
      -s "$ALIGNMENT" \
      -st DNA \
      -m MFP \
      -B 1000 \
      -alrt 1000 \
      -T 1 \
      --prefix "$PREFIX"
done

echo
echo "Finished: $(date)"
echo "Treefiles:"
ls -lh 04_trees/trim_*.treefile
```

<a id="step-06"></a>

### 6. Plot rooted trees with support labels

Requires 05_plots to exist; root display uses Catocala fraxini. All five results are plotted for comparison.

**Source:** `Haplotype_trees/scripts/09_plot_trim_sensitivity_trees.R`  
**Save as:** `09_plot_trim_sensitivity_trees.R`  
**Archive treatment:** Complete source, verbatim.

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

tree_files <- c(
  all6 = "04_trees/trim_all6_present.treefile",
  hap1out = "04_trees/trim_haplotypes_plus_1outgroup.treefile",
  hap2out = "04_trees/trim_haplotypes_plus_2outgroups.treefile",
  atleast5 = "04_trees/trim_atleast5of6.treefile",
  atleast4 = "04_trees/trim_atleast4of6.treefile"
)

titles <- c(
  all6 = "All six taxa present",
  hap1out = "A/B/C present + at least one outgroup",
  hap2out = "A/B/C present + at least two outgroups",
  atleast5 = "At least five of six taxa present",
  atleast4 = "At least four of six taxa present"
)

descriptions <- c(
  all6 = "Strictest alignment: only columns where all haplotypes and all outgroups have bases.",
  hap1out = "Focal-haplotype alignment: A, B, and C retained at every site; at least one outgroup anchors rooting.",
  hap2out = "Balanced main tree: A, B, and C retained at every site; at least two outgroups present.",
  atleast5 = "Less strict occupancy filter: columns retained if five or more taxa have bases.",
  atleast4 = "Most permissive occupancy filter: columns retained if four or more taxa have bases."
)

pretty_labels <- c(
  C_fra_cortex_region = "Catocala fraxini",
  D_rog_cortex_region = "Dysgonia rogenhoferi",
  E_dep_cortex_region = "Eilema depressum",
  Haplotype_A_cortex_region = "Haplotype A",
  Haplotype_B_cortex_region_revcomp = "Haplotype B (revcomp)",
  Haplotype_C_cortex_region = "Haplotype C"
)

for (f in tree_files) {
  if (!file.exists(f)) {
    stop("Missing tree file: ", f)
  }
}

trees <- lapply(tree_files, read.tree)

# Root each tree on Catocala for visual consistency.
# This does not change the underlying inference; it just gives a consistent display.
trees <- lapply(trees, function(tr) {
  root(tr, outgroup = "C_fra_cortex_region", resolve.root = TRUE)
})

for (i in seq_along(trees)) {
  trees[[i]]$tip.label <- ifelse(
    trees[[i]]$tip.label %in% names(pretty_labels),
    pretty_labels[trees[[i]]$tip.label],
    trees[[i]]$tip.label
  )
}

pdf("05_plots/trim_sensitivity_trees.pdf", width = 11, height = 14)
par(mfrow = c(5, 1), mar = c(1.2, 1.2, 4.8, 1.2), xpd = NA)

for (nm in names(trees)) {
  tr <- trees[[nm]]

  plot(
    tr,
    type = "phylogram",
    use.edge.length = TRUE,
    cex = 0.85,
    label.offset = 0.015,
    no.margin = TRUE
  )

  title(main = titles[[nm]], line = 3.1, cex.main = 1.05)
  mtext(descriptions[[nm]], side = 3, line = 1.8, cex = 0.78)

  # Show branch support labels where IQ-TREE stored them as node labels.
  nodelabels(
    text = tr$node.label,
    frame = "none",
    cex = 0.65,
    adj = c(1.1, -0.2)
  )

  add.scale.bar(length = NULL, cex = 0.7)
}

dev.off()

cat("Wrote: 05_plots/trim_sensitivity_trees.pdf\n")
```
