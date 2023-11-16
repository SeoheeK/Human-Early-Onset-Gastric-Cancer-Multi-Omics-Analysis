# Guideline for Users

## 01 Data Acquistion
Exome and mRNA sequencing data were obtained from the NCBI SRA (PRJNA505380 and PRJNA508414) database(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122401)
Global, phospho-, and N-glycoproteomics data were retrieved from the CPTAC data portal (PDC000214~216)(https://pdc.cancer.gov/pdc/browse/filters/pdc_study_id:PDC000214%7CPDC000215%7CPDC000216).

## 02 Data Pre-Processing
### Pre-processing

#### WES-Preprocessing

The exome sequencing data were aligned to the human GRCh38 reference genome using Bowtie 2 (v2.5.2)(1). Low-quality bases and PCR duplicates were removed by using Trimmomatic(v0.39)(2). Realignment and recalibration were performed using GATK4 (v4.4.0.0)(3). SNV candidates were identified in both tumor and control samples. To identify somatic mutations, we used Mutect (v1.1.7)(4) and Strelka (v1.0.7)(5). The MutSigCV (v1.41)(6) tool was utilized to discover significantly mutated genes among somatic SNVs in the 76 microsatellite-stable EOGC patients. 

#### RNA-seq Preprocessing

The .sra files were subjected to fasterq-dump using SRA Toolkit (v3.0.7). Trimming of adapters and removal of low-quality sequences were executed with Trimmomatic (v0.39)(4). Transcriptome alignment was performed utilizing STAR (v2.7.11a)(5), employing a genome index generated from the GRCh38 (Human) sequence. Finally, quantification of mRNA reads was achieved using RSEM (v1.3.3)(6). 

#### Proteomics-preprocessing

All LC-MS/MS (Liquid Chromatography Tandem Mass Spectrometry) data were processed with Post-Experiment Monoisotopic Mass Refinement (PE-MMR) for precursor mass correction and refinement. Tandem Mass Spectrometry (MS/MS) data for tissue pairs were analyzed for peptide identification using the MS-GF+ search engine (v.9387) with composite database (DB) of UniProt DB. Subsequently, ResultMerger (v5.4.16) was employed to consolidate multiple .mzid files into a single result file. PIPRegister (v0.6) was utilized to calculate Precursor Ion Purity (PIP) from mzXML files, measuring the confidence of scans. PSMs (Peptide Spectrum Matches) were then validated by filtering based on a PIP threshold (PIP > 70). Quantile normalization was applied, followed by the computation of Normalized Fold Changes. For additional analysis, PSMs with a False Discovery Rate (FDR) less than 0.01 were used. 

#### Citation
If you want more detailed information about the preprocessing steps, please refer to the parent project 'Proteogenomic Characterization of Human Early-Onset Gastric Cancer' which serves as the foundation for this project (https://www.cell.com/cancer-cell/fulltext/S1535-6108(18)30574-9#secsectitle0080).


## 03 Data Analysis
### Correlation profiles and pathway analysis 

#### (1) Mutations & modified peptides 
Phosphopeptides and N-glycopeptides detected in more than 50% out of 80 patients were used for this analysis. Log2 scale fold-changes between tumor and matched normal samples were normalized using quantile normalization. Significantly changed peptides between mutation and wildtype samples with p<0.05 were selected (Wilcoxon rank-sum test). We further filtered the data by choosing the peptides with median in the mutation group >0, median in the non-mutation group <0.
The selected significant peptides mapped with HGNC symbols were subjected to over-representation analysis using ConsensusPathDB Pathway-based sets from Wikipathways, KEGG, PID, and Biocarta were utilized. 

#### (2) mRNAs and proteins

Spearman's correlation coefficient was computed using mRNA and protein T/N log2(fold-change) data for which protein and mRNA abundances were available in at least 30% of patients. Patients or pairs of mRNAs and proteins with FDR less than 0.01 (by Benjamini-Hochberg method) were regarded as having significant correlations.
Enrichment analysis encompassed KEGG pathways and was executed on two separate sets of 500 genes: one set characterized by the strongest correlations, and the other marked by the weakest correlations. These enrichment analyses were conducted through DAVID(7).


### Subtyping EOGC patients

#### (1) Clustering by each dataset

Clustering analyses were executed on the filtered samples, as outlined in the figure. Following quantile normalization, only molecules falling within the top 10-20% of MADs were then subjected to CNMF clustering, utilizing CancerSubtypes (v1.26.0)(7). The ultimate number of clusters was established based on the value of 'k' that yielded the most stable k-cluster decomposition.

#### (2) Differentially expressed molecules

Molecular signatures defining the subtypes through clustering were identified. Initially, we compared the log2(fold-change) in the subtype to those in other subtypes using a t-test, and we selected molecules with a p-value or adjusted p-value <0.05. We applied additional filtering criteria, retaining molecules that met the following conditions: median within the subtype > 0, median in the remaining patients < 0.

#### (3) Integrative clustering

Integrative clustering was conducted based on the results obtained from molecular subtyping, employing ConsensusClusterPlus (v.1.64.0)(8). Initially, each subtype was transformed into an indicator vector, where '1' represented samples belonging to the subtype, and '0' denoted samples that did not belong to it. Four types of data were consolidated into a single indicator matrix. The clustering was performed using sample resampling at an 80% rate, with 1000 iterations of hierarchical clustering and the Pearson correlation as the dissimilarity measure.


### Reference

1) Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012;9:357-359.
2) McKenna A, Hanna M, Ba5nks E et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation D3NA sequencing data. Genome Res. 2010;20:1297-303.
3) Lawrence, M., Stojanov, P., Polak, P. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature. 2013;499:14–218.
4)  Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30(15):2114-2120.
5) Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.
6) Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics. 2011;12:323.
7) Huang DW, Sherman BT, Tan Q, et al. The DAVID Gene Functional Classification Tool: a novel biological module-centric algorithm to functionally analyze large gene lists. Genome Biol. 2007;8(9):R183. 
8) Wilkerson MD, Hayes DN. ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. Bioinformatics. 2010;26(12):1572-1573.
