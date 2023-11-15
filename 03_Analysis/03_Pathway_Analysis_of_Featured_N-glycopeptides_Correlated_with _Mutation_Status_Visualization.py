#-------------------------------------------------------------------------------------------------
# 03 Correlation between Somatic Mutations and N-glycopeptides
#-------------------------------------------------------------------------------------------------


"""
purpose:
Cellular processes represented by proteins whose N-glycosylation levels 
correlated with nonsynonymous somatic SNVs in ARID1A. 

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#consensuspathdb 결과값을 input으로 받아오기
df = pd.read_csv("./ORA_results_glyco_arid1a.tab",sep="\t")

df = df[["pathway","p-value"]]

#log10취하기
df["log10_pvalue"]=np.log10(df["p-value"])

#abs로 x축 양수로 만들기
df["-log10(pvalue)"] = abs(df["log10_pvalue"])

#barplot 그리기
fig, ax = plt.subplots()
name = df["pathway"][::-1]
log = df["-log10(pvalue)"][::-1]
ax.barh(name,log)
plt.xlabel('-log10(p value)')
plt.savefig("glyco-log10_enrichment.png")