README for 10.5281/zenodo.22755120 
Last updated: Apr 28, 2026

Data repository accompanying:


Livraghi L., Cassily J.L., Brayer M., Carter A.T., Alqassar J.D., Sim S.B., Geib S.M., Hanly J.J., Martin A. Hyperdivergent haplotypes control melanic camouflage in a polymorphic moth.

Please consider citing this manuscript accordingly if you use the dataset.


====================================================================
OVERVIEW
====================================================================


This repository contains genome assemblies, phased haplotype sequences,
and annotation files used to analyze structural and sequence variation
at the ivory:mir-193 melanism locus in the moth Anticarsia gemmatalis.

Files are organized into folders corresponding to whole-genome reference
assemblies and locus-specific datasets used in Figures 2 and 3 of the manuscript.


Haplotype nomenclature:


A = Dark reference-like haplotype
B = alternate Dark-associated haplotype
C = Light-associated haplotype

The following folders include gene annotations (from NCBI RefSeq, or their LiftOff version), and TE annotations generated with EarlGrey:
* Genome_HaplotypeA_Dark_reference/
* Genome_HaplotypeB_Dark/
* Genome_HaplotypeC_Light_reference/


====================================================================
RELATION TO MANUSCRIPT FIGURES
====================================================================


Figure 2: Phased haplotypes and hyperdivergent locus comparisons
-> Fig2_16_phased_haplotypes/ 


Figure 3: Structural rearrangements and inversion analyses
-> Fig3_ivory_scaffolds/




====================================================================
FOLDER CONTENTS
====================================================================


1_Genome_HaplotypeA_Dark_reference/


Whole-genome reference primary assembly representing haplotype A.
This individual was used in the genome assembly and annotation process consisted of female pupa that resulted from 3 generations of Dark x Dark crossing, and is homozygous for the Dark Haplotype A.

This genome corresponds to the NCBI RefSeq genome for the species (primary assembly), and the sequence and RefSeq annotation files from NCBI are shared here for convenience. 
The TE annotation file generated with Earl Grey is added.

Files:
Dark_Reference_TEannotation.gff
GCF_050436995.1_ilAntGemm2_primary_genomic.fna.gz
GCF_050436995.1_ilAntGemm2_primary_genomic.gff.gz


--------------------------------------------------------------------


2_Genome_HaplotypeB_Dark/


Whole-genome assembly carrying haplotype B at the melanism locus (hap2 of Dark_01 individual), LiftOff annotation transfer from the RefSeq annotation, and TE annotation

Files:
Agem_Dark_01_hap2.p_ctg.fa
Agem_Dark_01_hap2.liftoff.gff3_polished
Agem_Dark_01_hap2_TEannotations.gff


--------------------------------------------------------------------


3_Genome_HaplotypeC_Light_reference/


Whole-genome Light reference assembly carrying haplotype C, with its LiftOff annotation transfer from the RefSeq annotation, and TE annotation. The FASTA file corresponds to the HifiASM primary assembly (pseudo-assembly) of individual Light_04. 



Files:


Light_04_liftoff_annotation.gff3
Light_04_ragtag_pseudoassembly.fasta
Light_04_reference_HaplotypeC_TEannotations.gff




--------------------------------------------------------------------


4_Fig2_16_phased_haplotypes/


Contains 16 phased haplotype contigs derived from 8 PacBio HiFi sequenced
individuals (4 Dark, 4 Light), each with hap1 and hap2 assemblies.


Files:


Agem_Dark_01.bp.hap1.p_ctg.fa
Agem_Dark_01.bp.hap2.p_ctg.fa
Agem_Dark_02.bp.hap1.p_ctg.fa
Agem_Dark_02.bp.hap2.p_ctg.fa
Agem_Dark_03.bp.hap1.p_ctg.fa
Agem_Dark_03.bp.hap2.p_ctg.fa
Agem_Dark_04.bp.hap1.p_ctg.fa
Agem_Dark_04.bp.hap2.p_ctg.fa
Agem_Light_01.bp.hap1.p_ctg.fa
Agem_Light_01.bp.hap2.p_ctg.fa
Agem_Light_02.bp.hap1.p_ctg.fa
Agem_Light_02.bp.hap2.p_ctg.fa
Agem_Light_03.bp.hap1.p_ctg.fa
Agem_Light_03.bp.hap2.p_ctg.fa
Agem_Light_04.bp.hap1.p_ctg.fa
Agem_Light_04.bp.hap2.p_ctg.fa


--------------------------------------------------------------------


5_Fig3_ivory_scaffolds/


Contains extracted scaffold/contig sequences centered on the cortex /
ivory:mir-193 candidate interval for all phased haplotypes and references.


Files:
Agem_Dark_01.hap1_cortex_contig.fa
Agem_Dark_01.hap2_cortex_contig.fa
Agem_Dark_02.hap1_cortex_contig.fa
Agem_Dark_02.hap2_cortex_contig.fa
Agem_Dark_03.hap1_cortex_contig.fa
Agem_Dark_03.hap2_cortex_contig.fa
Agem_Dark_04.hap1_cortex_contig.fa
Agem_Dark_04.hap2_cortex_contig.fa
Agem_Light_01.hap1_cortex_contig.fa
Agem_Light_01.hap2_cortex_contig.fa
Agem_Light_02.hap1_cortex_contig.fa
Agem_Light_02.hap2_cortex_contig.fa
Agem_Light_03.hap1_cortex_contig.fa
Agem_Light_03.hap2_cortex_contig.fa
Agem_Light_04.hap1_cortex_contig.fa
Agem_Light_04.hap2_cortex_contig.fa
Agem_Reference_reassembly.hap1_cortex_contig.fa
Agem_Reference_reassembly.hap2_cortex_contig.fa


====================================================================
FILE FORMAT NOTES
====================================================================


.fa / .fasta / .fna     Genome or contig nucleotide FASTA files
.gff / .gff3           Gene or repeat annotations in GFF format
.gz                    gzip-compressed files




====================================================================
END
====================================================================