# WES-Preprocessing

The exome sequencing data were aligned to the human GRCh38 reference genome using Bowtie 2 (v2.5.2)6). Low-quality bases and PCR duplicates were removed. Realignment and recalibration were performed using GATK4 (v4.4.0.0)7). SNV candidates were identified in both tumor and control samples. The MutSigCV (v1.41)8) tool was utilized to discover significantly mutated genes among somatic SNVs in the 76 microsatellite-stable EOGC patients.


# RNA-seq Preprocessing

The .sra files were subjected to fasterq-dump using SRA Toolkit (v3.0.7). Trimming of adapters and removal of low-quality sequences were executed with Trimmomatic (v0.39)9). Transcriptome alignment was performed utilizing STAR (v2.7.11a)10), employing a genome index generated from the GRCh38 (Human) sequence. Finally, quantification of mRNA reads was achieved using RSEM (v1.3.3)11).


# Proteomics-preprocessing

All data were processed by PE-MMR for precursor mass correction and refinement. Resultant MS/MS data were identified as peptide-spectrum matches (PSMs) by MS-GF+ search engine (v.9387) at composite database (DB) of UniProt DB and specific modifications depending on the type of proteome.


### Reference

6) Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012;9:357-359.
7) McKenna A, Hanna M, Banks E et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20:1297-303.
8) Lawrence, M., Stojanov, P., Polak, P. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature. 2013;499:14–218.
9) Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30(15):2114-2120.
10) Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.
11) Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics. 2011;12:323.


