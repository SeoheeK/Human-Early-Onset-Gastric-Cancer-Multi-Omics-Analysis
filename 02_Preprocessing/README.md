# Preprocessing
This portion of the code has been implemented by citing the work of another individual, and as such, we will keep it private.

## WES-Preprocessing

The exome sequencing data were aligned to the human GRCh38 reference genome using Bowtie 2 (v2.5.2)(1). Low-quality bases and PCR duplicates were removed. Realignment and recalibration were performed using GATK4 (v4.4.0.0)(2). SNV candidates were identified in both tumor and control samples. The MutSigCV (v1.41)(3) tool was utilized to discover significantly mutated genes among somatic SNVs in the 76 microsatellite-stable EOGC patients. 

## RNA-seq Preprocessing

The .sra files were subjected to fasterq-dump using SRA Toolkit (v3.0.7). Trimming of adapters and removal of low-quality sequences were executed with Trimmomatic (v0.39)(4). Transcriptome alignment was performed utilizing STAR (v2.7.11a)(5), employing a genome index generated from the GRCh38 (Human) sequence. Finally, quantification of mRNA reads was achieved using RSEM (v1.3.3)(6). 

## Proteomics-preprocessing

 All LC-MS/MS (Liquid Chromatography Tandem Mass Spectrometry) data were processed with Post-Experiment Monoisotopic Mass Refinement (PE-MMR) for precursor mass correction and refinement. Tandem Mass Spectrometry (MS/MS) data for tissue pairs were analyzed for peptide identification using the MS-GF+ search engine (v.9387) with composite database (DB) of UniProt DB. Subsequently, ResultMerger (v5.4.16) was employed to consolidate multiple .mzid files into a single result file. PIPRegister (v0.6) was utilized to calculate Precursor Ion Purity (PIP) from mzXML files, measuring the confidence of scans. PSMs (Peptide Spectrum Matches) were then validated by filtering based on a PIP threshold (PIP > 70). Quantile normalization was applied, followed by the computation of Normalized Fold Changes. For additional analysis, PSMs with a False Discovery Rate (FDR) less than 0.01 were used. 

## Citation
If you want more detailed information about the preprocessing steps, please refer to the parent project 'Proteogenomic Characterization of Human Early-Onset Gastric Cancer' which serves as the foundation for this project(https://www.cell.com/cancer-cell/fulltext/S1535-6108(18)30574-9#secsectitle0080).


### Reference

1) Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012;9:357-359.
2) McKenna A, Hanna M, Ba5nks E et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation D3NA sequencing data. Genome Res. 2010;20:1297-303.
3) Lawrence, M., Stojanov, P., Polak, P. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature. 2013;499:14–218.
4)  Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30(15):2114-2120.
5) Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.
6) Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics. 2011;12:323.



