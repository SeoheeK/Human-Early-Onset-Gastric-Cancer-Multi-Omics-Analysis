# Data Analysis
## Correlation profiles and pathway analysis 

### (1) Mutations & modified peptides 
Phosphopeptides and N-glycopeptides detected in more than 50% out of 80 patients were used for this analysis. Log2 scale fold-changes between tumor and matched normal samples were normalized using quantile normalization. Significantly changed peptides between mutation and wildtype samples with p<0.05 were selected (Wilcoxon rank-sum test). We further filtered the data by choosing the peptides with median in the mutation group >0, median in the non-mutation group <0.
The selected significant peptides mapped with HGNC symbols were subjected to over-representation analysis using ConsensusPathDB Pathway-based sets from Wikipathways, KEGG, PID, and Biocarta were utilized. 

If you want to conduct this analysis, utilize the code provided in Section 01 & 02 & 03 and the preprocessed data provided in 00_DB folder.

### (2) mRNAs and proteins

Spearman's correlation coefficient was computed using mRNA and protein T/N log2(fold-change) data for which protein and mRNA abundances were available in at least 30% of patients. Patients or pairs of mRNAs and proteins with FDR less than 0.01 (by Benjamini-Hochberg method) were regarded as having significant correlations.
Enrichment analysis encompassed KEGG pathways and was executed on two separate sets of 500 genes: one set characterized by the strongest correlations, and the other marked by the weakest correlations. These enrichment analyses were conducted through DAVID.

If you want to conduct this analysis, utilize the code provided in Section 04 and the preprocessed data provided in 00_DB folder.


## Subtyping EOGC patients

### (1) Clustering by each dataset

Clustering analyses were executed on the filtered samples, as outlined in the figure. Following quantile normalization, only molecules falling within the top 10-20% of MADs were then subjected to CNMF clustering, utilizing CancerSubtypes (v1.26.0). The ultimate number of clusters was established based on the value of 'k' that yielded the most stable k-cluster decomposition.

### (2) Differentially expressed molecules

Molecular signatures defining the subtypes through clustering were identified. Initially, we compared the log2(fold-change) in the subtype to those in other subtypes using a t-test, and we selected molecules with a p-value or adjusted p-value <0.05. We applied additional filtering criteria, retaining molecules that met the following conditions: median within the subtype > 0, median in the remaining patients < 0.

### (3) Integrative clustering

Integrative clustering was conducted based on the results obtained from molecular subtyping, employing ConsensusClusterPlus (v.1.64.0). Initially, each subtype was transformed into an indicator vector, where '1' represented samples belonging to the subtype, and '0' denoted samples that did not belong to it. Four types of data were consolidated into a single indicator matrix. The clustering was performed using sample resampling at an 80% rate, with 1000 iterations of hierarchical clustering and the Pearson correlation as the dissimilarity measure.

If you want to conduct this analysis, utilize the code provided in Section 05 and the preprocessed data provided in 00_DB folder.
