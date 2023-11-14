"""
purpose: 
Cellular processes represented by proteins whose phosphorylation levels 
correlated with nonsynonymous somatic SNVs in ARID1A, CDH1, and RHOA.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#consensuspathdb 결과값을 input으로 받아오기
arid1a_df = pd.read_csv("./ORA_results_arid1a.tab",sep="\t")
rhoa_df = pd.read_csv("./ORA_results_rhoa.tab",sep="\t")
cdh1_df = pd.read_csv("./ORA_results_cdh1.tab",sep="\t")

#필요한 컬럼만 가져오기
arid1a_df = arid1a_df[["pathway","q-value"]]
cdh1_df = cdh1_df[["pathway","q-value"]]
rhoa_df = rhoa_df[["pathway","q-value"]]

#log10취하기
arid1a_df["log10_qvalue"]=np.log10(arid1a_df["q-value"])
rhoa_df["log10_qvalue"]=np.log10(rhoa_df["q-value"])
cdh1_df["log10_qvalue"]=np.log10(cdh1_df["q-value"])

#abs로 x축 양수로 만들기
arid1a_df["-log10(qvalue)"] = abs(arid1a_df["log10_qvalue"])

rhoa_df["-log10(qvalue)"] = abs(rhoa_df["log10_qvalue"])

cdh1_df["-log10(qvalue)"] = abs(cdh1_df["log10_qvalue"])

#input data from consensuspathdb
arid1a_df = pd.read_csv("arid1a_phospho_005_pqvalue.csv")
cdh1_df = pd.read_csv("cdh1_phospho_005_pqvalue.csv")
rhoa_df = pd.read_csv("rhoa_phospho_005_pqvalue.csv")

#pathway name과 qvalue 컬럼만 가져오기
arid1a_q = arid1a_df[["pathway","-log10(qvalue)"]]
cdh1_q = cdh1_df[["pathway","-log10(qvalue)"]]
rhoa_q = rhoa_df[["pathway","-log10(qvalue)"]]

arid1a_q_sort = arid1a_q.sort_values(by="-log10(qvalue)", ascending=False)
cdh1_q_sort = cdh1_q.sort_values(by="-log10(qvalue)", ascending=False)
rhoa_q_sort = rhoa_q.sort_values(by="-log10(qvalue)", ascending=False)

#top10 가져오기
cdh1_top10 = cdh1_q_sort.iloc[:10,:]
arid1a_top10 = arid1a_q_sort.iloc[:10,:]
rhoa_top10 = rhoa_q_sort.iloc[:10,:]

#세가지 gene에 대한 pathway의 합집합 생성
all_df = pd.DataFrame(list(arid1a_top10["pathway"].to_list()+cdh1_top10["pathway"].to_list()+rhoa_top10["pathway"].to_list()),columns=["pathway"])

#세가지 gene에 대해 all_df 결과에 있는 값들 가져오기
cdh1_select = cdh1_q[cdh1_q["pathway"].isin(all_df["pathway"])]
arid1a_select = arid1a_q[arid1a_q["pathway"].isin(all_df["pathway"])]
rhoa_select = rhoa_q[rhoa_q["pathway"].isin(all_df["pathway"])]

#3개 유전자 결과 병합
all_df_cdh1 =pd.merge(cdh1_select, all_df, on="pathway", how = "outer").rename(columns={"-log10(qvalue)":"CDH1"})
all_df_rhoa =pd.merge(rhoa_select, all_df_cdh1, on="pathway", how = "outer").rename(columns={"-log10(qvalue)":"RHOA"})
all_df_arid1a =pd.merge(arid1a_select, all_df_rhoa, on="pathway", how = "outer").rename(columns={"-log10(qvalue)":"ARID1A"})

#결과 저장
all_df_arid1a.drop_duplicates("pathway").to_csv("-log10pvalue_3genes_pathway.csv",index=False)