# Correlation profiles and pathway analysis 

## Mutations & modified peptides 
Phosphopeptides and N-glycopeptides detected in more than 50% out of 80 patients were used for this analysis. Log2 scale fold-changes between tumor and matched normal samples were normalized using quantile normalization. Significantly changed peptides between mutation and wildtype samples with p<0.05 were selected (Wilcoxon rank-sum test). We further filtered the data by choosing the peptides with 1) median in the mutation group >0, 2) median in the non-mutation group <0.
The selected significant peptides mapped with HGNC symbols were subjected to over-representation analysis using ConsensusPathDB Pathway-based sets from Wikipathways, KEGG, PID, and Biocarta were utilized. 

## mRNAs and proteins

Spearman's correlation coefficient was computed using mRNA and protein T/N log2(fold-change) data for which protein and mRNA abundances were available in at least 30% of patients. Patients or pairs of mRNAs and proteins with FDR less than 0.01 (by Benjamini-Hochberg method) were regarded as having significant correlations.
Enrichment analysis encompassed KEGG pathways and was executed on two separate sets of 500 genes: one set characterized by the strongest correlations, and the other marked by the weakest correlations. These enrichment analyses were conducted through DAVID12).


# Subtyping EOGC patients

## Clustering by each dataset

Clustering analyses were executed on the filtered samples, as outlined in the figure. Following quantile normalization, only molecules falling within the top 10-20% of MADs were then subjected to CNMF clustering, utilizing CancerSubtypes (v1.26.0)13). The ultimate number of clusters was established based on the value of 'k' that yielded the most stable k-cluster decomposition.

## Differentially expressed molecules

Molecular signatures defining the subtypes through clustering were identified. Initially, we compared the log2(fold-change) in the subtype to those in other subtypes using a t-test, and we selected molecules with a p-value or adjusted p-value <0.05. We applied additional filtering criteria, retaining molecules that met the following conditions: 1) median within the subtype > 0, 2) median in the remaining patients < 0.

## Integrative clustering

Integrative clustering was conducted based on the results obtained from molecular subtyping, employing ConsensusClusterPlus (v.1.64.0)14). Initially, each subtype was transformed into an indicator vector, where '1' represented samples belonging to the subtype, and '0' denoted samples that did not belong to it. Four types of data were consolidated into a single indicator matrix. The clustering was performed using sample resampling at an 80% rate, with 1000 iterations of hierarchical clustering and the Pearson correlation as the dissimilarity measure. 


### References
1) Merchant SJ, Kim J, Choi AH, Sun V, Chao J, Nelson R. A rising trend in the incidence of advanced gastric cancer in young Hispanic men. Gastric Cancer. 2017;20(2):226-234.
2) De B, Rhome R, Jairam V, et al. Gastric adenocarcinoma in young adult patients: patterns of care and survival in the United States. Gastric Cancer. 2018;21(6):889-899.
3) Kulig J, Popiela T, Kolodziejczyk P, et al. Clinicopathological profile and long-term outcome in young adults with gastric cancer: multicenter evaluation of 214 patients. Langenbecks Arch Surg. 2008;393(1):37-43.
4) Isobe T, Hashimoto K, Kizaki J, et al. Characteristics and prognosis of gastric cancer in young patients. Oncol Rep. 2013;30(1):43-49. 
5) Mun DG, Bhin J, Kim S, et al. Proteogenomic Characterization of Human Early-Onset Gastric Cancer. Cancer Cell. 2019;35(1):111-124.e10.
6) Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012;9:357-359.
7) McKenna A, Hanna M, Banks E et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20:1297-303.
8) Lawrence, M., Stojanov, P., Polak, P. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature. 2013;499:14–218.
9) Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30(15):2114-2120.
10) Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.
11) Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics. 2011;12:323.
12) Xu T, Le TD, Liu L, et al. CancerSubtypes: an R/Bioconductor package for molecular cancer subtype identification, validation and visualization. Bioinformatics. 2017;33(19):3131-3133.
13) Huang DW, Sherman BT, Tan Q, et al. The DAVID Gene Functional Classification Tool: a novel biological module-centric algorithm to functionally analyze large gene lists. Genome Biol. 2007;8(9):R183. 
14) Wilkerson MD, Hayes DN. ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. Bioinformatics. 2010;26(12):1572-1573.


