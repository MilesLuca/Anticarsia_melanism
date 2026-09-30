# Code for *Hyperdivergent haplotypes control melanic camouflage in a polymorphic moth*

Companion code for the forthcoming *Current Biology* paper on melanic camouflage in *Anticarsia gemmatalis*.

**Start with the analysis index below.** This repository contains the 21 annotated workflows from [Zenodo version v1](https://zenodo.org/records/22755120), together with the original code and dataset READMEs.

- **Archived code and data:** [DOI 10.5281/zenodo.22755120](https://doi.org/10.5281/zenodo.22755120).
- **Large datasets:** download [Archive.zip (2.2 GB)](https://zenodo.org/records/22755120/files/Archive.zip?download=1) from Zenodo. See [README_DATASET.TXT](README_DATASET.TXT) for its contents.
- **Original documentation:** [README_code.md](README_code.md), preserved as deposited.
- **License:** [CC BY-SA 4.0](LICENSE.md), as listed on the Zenodo record.
- **Provenance:** the workflows and original READMEs match the published Zenodo checksums; see [transfer verification](VALIDATION.md).

The code is distributed as annotated Markdown documents containing Bash, Python and R code blocks. Save each block under its stated filename and follow that workflow's input, software and execution instructions. These workflows require separately obtained datasets and appropriate computing resources.

## Citation

Please cite the companion paper and the archived code/data when using this work. The paper's final publication details will be added when available.

> Livraghi, L. (2026). *Hyperdivergent haplotypes control melanic camouflage in a polymorphic moth* (v1) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.22755120

GitHub-readable citation metadata is provided in [CITATION.cff](CITATION.cff). The original code README retains its historical OSF reference; this GitHub copy was obtained from the Zenodo deposit above.

## Analysis index

| Analysis | Paper connection |
|---|---|
| [01 — HiFi haplotype assembly and locus identification](analyses/01_hifi_haplotype_assembly.md) | Figure 2D–E; Figures S10–S11; phased assembly resource |
| [02 — Light reference scaffolding and gene annotation](analyses/02_light_reference_and_gene_annotation.md) | Light-reference Pool-GWAS (S7); annotated haplogenomes used in Figures 2–3 and S12–S14 |
| [03 — Haplotype variant calling and genotype visualization](analyses/03_haplotype_variants_and_genotypes.md) | Figure 2E; Figure S11B |
| [04 — Haplotype sequence dotplots](analyses/04_haplotype_dotplots.md) | Figure 2D, Figure 3B–C and Figure S10 |
| [05 — Repeat library generation and annotation](analyses/05_repeat_annotation.md) | Figure 3A/E; Figures S12–S14; deposited repeat annotations |
| [06 — Repeat composition and donut plots](analyses/06_repeat_composition.md) | Figure 3A |
| [07 — Repeat enrichment in matched genomic windows](analyses/07_repeat_enrichment.md) | Figure S12; percentile labels in Figure 3A |
| [08 — Inversion-region LASTZ synteny ribbons](analyses/08_inversion_synteny.md) | Figure 3E |
| [09 — Haplotype phylogeny and occupancy sensitivity](analyses/09_haplotype_phylogeny.md) | Figure S14C–D |
| [10 — Pool-seq population-statistic plots](analyses/10_poolseq_population_statistics.md) | Figure 2B; population-statistic panels of S8–S9 |
| [11 — Pool-seq coverage at cortex and in exon/non-exon windows](analyses/11_poolseq_coverage.md) | Figure S6 |
| [12 — mir-193 CRISPR whole-genome genotyping coverage](analyses/12_crispr_mir193_wgs.md) | Figure S15 coverage |
| [13 — ivory CRISPR Nanopore amplicon coverage](analyses/13_crispr_ivory_amplicons.md) | Figure S16B |
| [14 — Pool-seq GWAS mapping, statistics and plotting](analyses/14_poolseq_gwas.md) | Figure 1E; GWAS components of 2C, 3B–C and S7–S9; upstream pool fractions |
| [15 — Pooled haplotype frequencies and matching genotype tracks](analyses/15_poolseq_haplotype_frequencies.md) | Figure 2F; Figure S11A–C |
| [17 — Merian-element chromosome painting](analyses/17_merian_element_painting.md) | Figure S4A–C |
| [18 — Small-RNA sequencing and miRNA annotation](analyses/18_small_rnaseq.md) | Figure S5B, E–H |
| [19 — Pool-seq median depth and chromosome normalization](analyses/19_poolseq_median_depth.md) | Figure 2A |
| [20 — Calculate Pool-seq population statistics with grenedalf](analyses/20_grenedalf_population_statistics.md) | Figure 2B; population-statistic panels of S8–S9 |
| [21 — Standard RNA-seq alignment and ivory-region extraction](analyses/21_standard_rnaseq.md) | Figure S5A alignment input |
| [22 — Light-reference Pool-seq depth and chromosome normalization](analyses/22_light_reference_depth.md) | Figure S7C |


## Using the code

1. Choose a workflow from the index. Read its input and software requirements.
2. Save each fenced code block under its stated filename. Blocks marked as input records contain the supplied sample lists, chromosome mappings or reference sequences and should be saved as data files.
3. Replace `/path/to/project` and `/path/to/tools` with your own directories. The folder names illustrate how inputs and outputs connect; adapt the file locations consistently across dependent steps. Create the working and output directories used by each script.
4. Install the listed command-line tools and R/Python packages. Command-line tools must be available on `PATH`. Where an exact executed version is known, it is recorded in the workflow notes.
5. Run Bash blocks with `bash filename`, Python blocks with `python3 filename`, and R blocks with `Rscript filename`. Some Bash source filenames retain a `.slurm` extension; run them with Bash in the same way. No scheduler is required by the distributed scripts.
6. For scripts that process one sample at a time, provide the zero-based sample index as the first argument, for example `bash script.sh 0`. The permitted indices are stated and checked at the top of those scripts. Run each listed index before continuing to the next stage. Adjust thread counts to your available resources.
7. Follow the execution order within each workflow and supply its named input files. The package includes code and small input records; sequencing reads, assemblies, BAMs and large annotation files must be obtained separately.

## Data and naming conventions

| Resource | Manuscript accession or input identity |
|---|---|
| Dark reference primary assembly | `GCA_050436995.1`; annotated RefSeq `GCF_050436995.1` |
| Reference alternate haplotype | `GCA_050436975.1` |
| Reference HiFi reads used for reassembly | `SRR32404155` |
| Eight F2 HiFi individuals | `PRJNA1450809`; phased assemblies and annotations in the Zenodo dataset |
| Dark/Light pooled sequencing | `PRJNA1338508` |
| Standard RNA-seq | Original FASTQ identifiers and reference FASTA listed in analysis 21 |
| Small-RNA sequencing | `PRJNA1165401`; supplied adapter and reference sequences in analysis 18 |
| Light pseudoassembly | `agem_light_06.ragtag.final.fasta`, from the paper's Light 4 (Light04); available in the Zenodo dataset |
| Repeat library | Frozen Dark-reference EarlGrey library named in analysis 05 |
| Merian assignments | Pinned reference table and contributed R plotter in analysis 17 |
| CRISPR reads, six-taxon locus FASTA and prepared BED/GFF inputs | Required file identities are listed in the corresponding workflows; reads available from NCBI |

Accessions are those recorded in the manuscript. `Plain` denotes the Light morph. Haplotype A is the Dark reference (`ADGO`); B is Dark01 hap2 (`AD1H2`, contig `h2tg000012l`); C is the Light reference (`ALGO`, `chr8`). Dark-reference chromosome 8 is `NC_134752.1`.

<a id="sample-names"></a>

The eight HiFi individuals were renumbered consecutively in the paper for clarity. The following mapping applies to both haplotypes of each individual:

| Identifier suffix in code/input files | Dark individual in paper | Light individual in paper |
|---|---|---|
| 01 | Dark 1 | Light 1 |
| 02 | Dark 2 | Light 2 |
| 05 | Dark 3 | Light 3 |
| 06 | Dark 4 | Light 4 |

The scripts retain the original identifiers to match their input filenames and VCF sample names; use this table to interpret their labels. For example, `Agem_Light_06`, `Light_06` and the assembly prefix `agem_light_06` refer to the paper's Light 4 (Light04). Haplotype suffixes `hap1` and `hap2` keep their original meanings. This mapping concerns the HiFi individuals.

GFF/VCF positions and quoted genomic intervals are generally one-based; BED starts are zero-based and ends are exclusive. Some plots use SNP index or amplicon-local coordinates. Preserve each script's explicit conversion and reference sequence when adapting paths or inputs.