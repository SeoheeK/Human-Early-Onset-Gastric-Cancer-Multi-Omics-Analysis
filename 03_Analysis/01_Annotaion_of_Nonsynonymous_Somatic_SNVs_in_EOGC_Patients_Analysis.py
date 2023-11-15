"""
purpose:  
1. mutation types of significantly mutated genes in each patient
2. clinical and histological parameters for each patient 

"""

import pandas as pd
import numpy as np

#1. mutation type
#input data 가져오기: 각 mutation 별로 mutation type 명시
df = pd.read_excel("WES_somaticVariants.xlsx")

#관심있는 유전자에 해당하는 값만 가져오기
df_selected = df[df.Gene.isin(["CDH1","TP53","MUC5B","ARID1A","BANP","RHOA"])]

#Patient 컬럼 값으로 다른 데이터와 동일한 형식의 patient id로 만들어주기
df_selected["patientnum"] = df_selected["Patient"].str[:-1]
df_selected["patientnum"] = str("T")+df_selected["patientnum"]
df_selected["patientnum"] = df_selected['patientnum'].apply(lambda x: f'N{int(x[1:])-1}T{x[1:]}')

#dataframe을 구성할 세가지 컬럼만 가져오기
df_3col = df_selected[["patientnum","Gene","MutationType"]]

#위 데이터 프레임에 patient id에 해당하는 mutation type 정보를 gene name 컬럼으로 새로 생성해 가져오기
df_3col["MUC5B"] = df_3col.loc[df_3col['Gene'] == 'MUC5B', 'MutationType']
df_3col["ARID1A"] = df_3col.loc[df_3col['Gene'] == 'ARID1A', 'MutationType']
df_3col["BANP"] = df_3col.loc[df_3col['Gene'] == 'BANP', 'MutationType']
df_3col["RHOA"] = df_3col.loc[df_3col['Gene'] == 'RHOA', 'MutationType']
df_3col["CDH1"] = df_3col.loc[df_3col['Gene'] == 'CDH1', 'MutationType']
df_3col["TP53"] = df_3col.loc[df_3col['Gene'] == 'TP53', 'MutationType']

#원하는 정보만 추출
new = df_3col[["patientnum","MUC5B","ARID1A","BANP","RHOA","CDH1","TP53"]]

#중복값 제거
new_drop= new.drop_duplicates()

#patient id 하나당 환자 mutation type 정보로 정리하기
new_group = new.groupby(['patientnum']).agg(lambda x: x.dropna().unique())

no_mutation = list(set(patient_meta.sampleID.unique()) - set(new_group.index.unique()))

no_mut = pd.DataFrame(no_mutation,columns=["patientnum"])

#NaN값 정리
df_fillna = new_group.applymap(lambda x: np.nan if len(x) == 0 else x)

df_fillna.reset_index(inplace=True)

#patient id 순서대로 정렬, [] 제거
df_fillna["num"] = df_fillna['patientnum'].str.extract('N(\d+)T').astype(int)
df_fillna = df_fillna.sort_values(by="num").drop("num",axis=1).set_index("patientnum")
df_fillna = df_fillna.applymap(lambda x: str(x)[1:-1] if pd.notnull(x) else x)

#string 처리 제거
df_fillna = df_fillna.applymap(lambda x: str(x).strip("'") if pd.notnull(x) else x)

#mapping할 정보 정리
data_map = {'frameshift deletion':"Frame shift",
 'frameshift insertion':"Frame shift",
 'nonframeshift deletion':"In frame Indel",
 'nonsynonymous SNV':"Missense",
 'splicing':"Splice site",
 'stopgain SNV':"Nonsense"}

df_fillna.reset_index(inplace=True)

#환자 id 순서대로 sorting할 경우
final_output = df_fillna.apply(lambda x: x.apply(lambda val: data_map.get(val, val)))

#mutation 순서가 많은 순서
final_output_sort = final_output.sort_values(by="RHOA").sort_values(by="BANP").sort_values(by="ARID1A").sort_values(by="MUC5B").sort_values(by="TP53").sort_values(by="CDH1")
final_output_sort = final_output_sort[["patientnum","CDH1","TP53","MUC5B","ARID1A","BANP","RHOA"]]
final_output_sort = pd.concat([final_output_sort, no_mut],axis=0)

#저장
final_output_sort.fillna(0).to_csv("figure1a_mutationtype.csv",index=False)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#2. patient metadata
#patient metadata 가져오기
patient_meta = pd.read_csv("clinical_info.txt",sep="\t")

#heatmap input과 동일하게 patient id 형식 맞춰주기
patient_meta["Normal"]="N"+patient_meta["Normal"].str[:-1]
patient_meta["Tumor"]="T"+patient_meta["Tumor"].str[:-1]
patient_meta["sampleID"]=patient_meta["Normal"]+patient_meta["Tumor"]

#metadata중 annotation bar을 그릴 컬럼만 선택
patient_meta_selected = patient_meta[["sampleID","EBV","MSI","Gender","Histology (Lauren)"]]
patient_meta_selected = patient_meta_selected.rename(columns={"Histology (Lauren)":"Histology"})

#값 변환
patient_meta_selected.Histology.fillna("Others",inplace=True)
patient_meta_selected.EBV.replace("Negative","EBV-",inplace=True)
patient_meta_selected.EBV.replace("EBV(PIK3CAmut)","EBV+",inplace=True)
patient_meta_selected.EBV.replace("EBV","EBV+",inplace=True)
patient_meta_selected.EBV.replace("EBV(PIK3CAwt)","EBV+",inplace=True)

#저장
patient_meta_selected.set_index('sampleID').loc[final_output_sort['patientnum']].reset_index().to_csv("figure1a_patientmeta.csv",index=False)