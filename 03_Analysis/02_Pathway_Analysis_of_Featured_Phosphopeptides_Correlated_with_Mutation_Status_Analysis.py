#-------------------------------------------------------------------------------------------------
# 02 Correlation between Somatic Mutations and Phosphopeptides
#-------------------------------------------------------------------------------------------------

"""
purpose:
Phosphorylated peptides whose levels were upregulated in GC samples 
with nonsynonymous somatic SNVs in ARID1A, CDH1, and RHOA. 

"""

import pandas as pd
from scipy import stats
import qnorm
from scipy.stats import wilcoxon

#input data 가져오기
df = pd.read_csv("phosphopeptide_change.txt",sep="\t")

#consensuspathdb input 생성을 위해 peptide-symbol 컬럼만 가져오기
peptide_symbol = df[["Peptide","Symbol"]]

#컬럼에 patient id만 남기도록 정리
df = df.set_index("Peptide")
df.drop("Symbol",axis=1,inplace=True)

# NaN이 아닌 값이 50% 이상인 행만 선정
selected_rows = df[df.count(axis=1) / len(df.columns) >= 0.5]

# 새로운 데이터프레임 생성
phospho_fch = df.loc[selected_rows.index]

#quantile normalization
phospho_fch_t = phospho_fch.transpose()
df =qnorm.quantile_normalize(phospho_fch_t).transpose()

mutation = pd.read_excel("SuppleTable_SomaticMutation_80_nonsyn.xls")

#각 유전자에 돌연변이가 존재하는 patient id 가져오기
ARID1A_patient_id = mutation[mutation.Gene == "ARID1A"].Patient.to_list()
RHOA_patient_id = mutation[mutation.Gene == "RHOA"].Patient.to_list()
CDH1_patient_id = mutation[mutation.Gene == "CDH1"].Patient.to_list()

#df와 mutation의 patient id 형식 맞춰주기
arid1a = []
for item in ARID1A_patient_id:
    number_part, letter_part = item[:-1], item[-1]
    arid1a.append(f'N{int(number_part) - 1}{letter_part}{number_part}')

#df와 mutation의 patient id 형식 맞춰주기
RHOA = []
for item in RHOA_patient_id:
    number_part, letter_part = item[:-1], item[-1]
    RHOA.append(f'N{int(number_part) - 1}{letter_part}{number_part}')

#df와 mutation의 patient id 형식 맞춰주기
CDH1 = []
for item in CDH1_patient_id:
    number_part, letter_part = item[:-1], item[-1]
    CDH1.append(f'N{int(number_part) - 1}{letter_part}{number_part}')

#mutation과 wildtype 환자군 분류
arid1a_mutated = df.columns[df.columns.isin(arid1a)]
arid1a_nomutated = df.columns[~df.columns.isin(arid1a)]

CDH1_mutated = df.columns[df.columns.isin(CDH1)]
CDH1_nomutated = df.columns[~df.columns.isin(CDH1)]

RHOA_mutated = df.columns[df.columns.isin(RHOA)]
RHOA_nomutated = df.columns[~df.columns.isin(RHOA)]

arid1a_mt_df = df[arid1a_mutated]
arid1a_wt_df = df[arid1a_nomutated]

CDH1_mt_df = df[CDH1_mutated]
CDH1_wt_df = df[CDH1_nomutated]

RHOA_mt_df = df[RHOA_mutated]
RHOA_wt_df = df[RHOA_nomutated]

#wildtype 대비 mutation이 높게 발현된것들 선별하기 위한 비교 데이터프레임
CDH1_median_df = pd.DataFrame({"wildtype":CDH1_wt_df.median(axis=1),"mutation":CDH1_mt_df.median(axis=1)})
CDH1_median_df['mt>0&wt<0'] = CDH1_median_df.apply(lambda x: "o" if (x['wildtype'] < 0) and (x['mutation'] > 0) else "x", axis=1)

#wildtype 대비 mutation이 높게 발현된것들 선별하기 위한 비교 데이터프레임
RHOA_median_df = pd.DataFrame({"wildtype":RHOA_wt_df.median(axis=1),"mutation":RHOA_mt_df.median(axis=1)})
RHOA_median_df['mt>0&wt<0'] = RHOA_median_df.apply(lambda x: "o" if (x['wildtype'] < 0) and (x['mutation'] > 0) else "x", axis=1)

#wildtype 대비 mutation이 높게 발현된것들 선별하기 위한 비교 데이터프레임
arid1a_median_df = pd.DataFrame({"wildtype":arid1a_wt_df.median(axis=1),"mutation":arid1a_mt_df.median(axis=1)})
arid1a_median_df['mt>0&wt<0'] = arid1a_median_df.apply(lambda x: "o" if (x['wildtype'] < 0) and (x['mutation'] > 0) else "x", axis=1)

from scipy.stats import mannwhitneyu

#mannwhitneyu = wilcoxon sum test 진행
RHOA_wcxs_df = pd.DataFrame(columns=["col","statistic","pvalue"])
CDH1_wcxs_df = pd.DataFrame(columns=["col","statistic","pvalue"])
arid1a_wcxs_df = pd.DataFrame(columns=["col","statistic","pvalue"])

for idx in CDH1_wt_df.index:
    wt_data = CDH1_wt_df.loc[idx].dropna().values
    mt_data = CDH1_mt_df.loc[idx].dropna().values

    if wt_data.size > 0 and mt_data.size > 0:
        statistic, p_value = mannwhitneyu(wt_data, mt_data)
        result = pd.DataFrame([[idx, statistic, p_value]], columns=["col", "statistic", "pvalue"])
        CDH1_wcxs_df = CDH1_wcxs_df.append(result,ignore_index=True)

for idx in arid1a_wt_df.index:
    wt_data = arid1a_wt_df.loc[idx].dropna().values
    mt_data = arid1a_mt_df.loc[idx].dropna().values

    if wt_data.size > 0 and mt_data.size > 0:
        statistic, p_value = mannwhitneyu(wt_data, mt_data)
        result = pd.DataFrame([[idx, statistic, p_value]], columns=["col", "statistic", "pvalue"])
        arid1a_wcxs_df = arid1a_wcxs_df.append(result,ignore_index=True)

for idx in RHOA_wt_df.index:
    wt_data = RHOA_wt_df.loc[idx].dropna().values
    mt_data = RHOA_mt_df.loc[idx].dropna().values

    if wt_data.size > 0 and mt_data.size > 0:
        statistic, p_value = mannwhitneyu(wt_data, mt_data)
        result = pd.DataFrame([[idx, statistic, p_value]], columns=["col", "statistic", "pvalue"])
        RHOA_wcxs_df = RHOA_wcxs_df.append(result,ignore_index=True)

#p값 기준 0.05보다 작은 peptide 선별
RHOA_pval_005 = RHOA_wcxs_df[RHOA_wcxs_df.pvalue < 0.05].col.to_list()
CDH1_pval_005 = CDH1_wcxs_df[CDH1_wcxs_df.pvalue < 0.05].col.to_list()
arid1a_pval_005 = arid1a_wcxs_df[arid1a_wcxs_df.pvalue < 0.05].col.to_list()

#median값이 mutation에서 0보다 크고, wildtype에서 0보다 작은 peptide 선별
RHOA_median_mt_over0_wt_under0 =RHOA_median_df[RHOA_median_df['mt>0&wt<0'] == "o"].index.to_list()
CDH1_median_mt_over0_wt_under0 =CDH1_median_df[CDH1_median_df['mt>0&wt<0'] == "o"].index.to_list()
arid1a_median_mt_over0_wt_under0 =arid1a_median_df[arid1a_median_df['mt>0&wt<0'] == "o"].index.to_list()

#median, p값 교집합 기준 유의미한 peptide 선별
RHOA_sig = list(set(RHOA_median_mt_over0_wt_under0)&set(RHOA_pval_005))
CDH1_sig = list(set(CDH1_median_mt_over0_wt_under0)&set(CDH1_pval_005))
arid1a_sig = list(set(arid1a_median_mt_over0_wt_under0)&set(arid1a_pval_005))

RHOA_final_input = df[df.index.isin(RHOA_sig)]
CDH1_final_input = df[df.index.isin(CDH1_sig)]
arid1a_final_input = df[df.index.isin(arid1a_sig)]

#heatmap에서 색깔로 mutation과 wildtype을 구분하기 위해 mutation과 wildtype끼리 모아서 concat
RHOA_wt = RHOA_final_input[RHOA_final_input.columns[~RHOA_final_input.columns.isin(RHOA)]]
RHOA_mt = RHOA_final_input[RHOA_final_input.columns[RHOA_final_input.columns.isin(RHOA)]]
RHOA_final = pd.concat([RHOA_mt, RHOA_wt],axis=1)

CDH1_wt = CDH1_final_input[CDH1_final_input.columns[~CDH1_final_input.columns.isin(CDH1)]]
CDH1_mt = CDH1_final_input[CDH1_final_input.columns[CDH1_final_input.columns.isin(CDH1)]]
CDH1_final = pd.concat([CDH1_mt, CDH1_wt],axis=1)

arid1a_wt = arid1a_final_input[arid1a_final_input.columns[~arid1a_final_input.columns.isin(arid1a)]]
arid1a_mt = arid1a_final_input[arid1a_final_input.columns[arid1a_final_input.columns.isin(arid1a)]]
arid1a_final = pd.concat([arid1a_mt, arid1a_wt],axis=1)

#heatmap을 그리기 위해 -1 ~ 1 로 스케일링
from sklearn.preprocessing import MaxAbsScaler
scaler = MaxAbsScaler()
RHOA_final_cols = RHOA_final.columns
RHOA_final_scaler = RHOA_final.copy()
RHOA_final_scaler[RHOA_final_cols] = scaler.fit_transform(RHOA_final_scaler[RHOA_final_cols])

from sklearn.preprocessing import MaxAbsScaler
scaler = MaxAbsScaler()
CDH1_final_cols = CDH1_final.columns
CDH1_final_scaler = CDH1_final.copy()
CDH1_final_scaler[CDH1_final_cols] = scaler.fit_transform(CDH1_final_scaler[CDH1_final_cols])

from sklearn.preprocessing import MaxAbsScaler
scaler = MaxAbsScaler()
arid1a_final_cols = arid1a_final.columns
arid1a_final_scaler = arid1a_final.copy()
arid1a_final_scaler[arid1a_final_cols] = scaler.fit_transform(arid1a_final_scaler[arid1a_final_cols])

RHOA_final_scaler.fillna(0).to_csv("./RHOA_heatmapinput.csv",index=False)
CDH1_final_scaler.fillna(0).to_csv("./CDH1_heatmapinput.csv",index=False)
arid1a_final_scaler.fillna(0).to_csv("./arid1a_heatmapinput.csv",index=False)

#consensuspathdb input 생성을 위한 peptide-symbol match
arid1a_peptide = pd.DataFrame(arid1a_final_scaler.reset_index(),columns=["Peptide"])
pd.merge(peptide_symbol, arid1a_peptide, on="Peptide",how="inner").Symbol.to_csv("arid1a_symbol_phospho.csv",index=False)

CDH1_peptide = pd.DataFrame(CDH1_final_scaler.reset_index(),columns=["Peptide"])
pd.merge(peptide_symbol, CDH1_peptide, on="Peptide",how="inner").Symbol.to_csv("CDH1_symbol_phospho.csv",index=False)

RHOA_peptide = pd.DataFrame(RHOA_final_scaler.reset_index(),columns=["Peptide"])
pd.merge(peptide_symbol, RHOA_peptide, on="Peptide",how="inner").Symbol.to_csv("RHOA_symbol_phospho.csv",index=False)


#patient metadata: annotation bar inpt생성
#patient metadata 가져오기
patient_meta = pd.read_csv("clinical_info.txt",sep="\t")

#mutation 표시를 위해 해당 gene이 mutation된 환자들과 wildtype군 구분
RHOA = ["N111T112","N117T118","N213T214","N5357T5358","N63T64"]
CDH1 = ["N115T116","N117T118","N135T136","N137T138","N145T146","N13T14","N179T180","N17T18","N189T190","N215T216","N27T28","N29T30","N31T32","N53T54","N91T92","N95T96"]
arid1a = ["N117T118","N155T156","N231T232","N249T250","N33T34","N39T40","N65T66"]

#heatmap input과 동일하게 patient id 형식 맞춰주기
patient_meta["Normal"]="N"+patient_meta["Normal"].str[:-1]
patient_meta["Tumor"]="T"+patient_meta["Tumor"].str[:-1]
patient_meta["sampleID"]=patient_meta["Normal"]+patient_meta["Tumor"]

#metadata중 annotation bar을 그릴 컬럼만 선택
patient_meta_selected = patient_meta[["sampleID","EBV","MSI","Gender","Histology (Lauren)"]]
patient_meta_selected = patient_meta_selected.rename(columns={"Histology (Lauren)":"Histology"})

patient_meta_selected.Histology.fillna("Others",inplace=True)
patient_meta_selected.EBV.replace("Negative","EBV-",inplace=True)
patient_meta_selected.EBV.replace("EBV(PIK3CAmut)","EBV+",inplace=True)
patient_meta_selected.EBV.replace("EBV","EBV+",inplace=True)
patient_meta_selected.EBV.replace("EBV(PIK3CAwt)","EBV+",inplace=True)

arid1a_patient_meta_selected = patient_meta_selected.copy()
rhoa_patient_meta_selected = patient_meta_selected.copy()
cdh1_patient_meta_selected = patient_meta_selected.copy()

#arid1a에 mutation 여부에 대한 환자 정보 patient_meta_selected에 컬럼 추가하기
arid1a_patient_meta_selected["Mutation"] = arid1a_patient_meta_selected.sampleID.apply(lambda x: "Mut" if x in arid1a else "WT")

#CDH1에 mutation 여부에 대한 환자 정보 patient_meta_selected에 컬럼 추가하기
cdh1_patient_meta_selected["Mutation"] = cdh1_patient_meta_selected.sampleID.apply(lambda x: "Mut" if x in CDH1 else "WT")

#RHOA에 mutation 여부에 대한 환자 정보 patient_meta_selected에 컬럼 추가하기
rhoa_patient_meta_selected["Mutation"] = rhoa_patient_meta_selected.sampleID.apply(lambda x: "Mut" if x in RHOA else "WT")

#mutation과 wildtype을 구분해서 heatmap annotation을 그리기 위해 dataframe sorting
arid1a_patient_mt = arid1a_patient_meta_selected[arid1a_patient_meta_selected.sampleID.isin(arid1a)]
arid1a_patient_wt = arid1a_patient_meta_selected[~arid1a_patient_meta_selected.sampleID.isin(arid1a)]
arid1a_patient_meta_sorted = pd.concat([arid1a_patient_mt,arid1a_patient_wt],axis=0).set_index("sampleID")
arid1a_patient_meta_sorted.to_csv("./arid1a_patient_meta.csv",index=False)

cdh1_patient_mt = cdh1_patient_meta_selected[cdh1_patient_meta_selected.sampleID.isin(CDH1)]
cdh1_patient_wt = cdh1_patient_meta_selected[~cdh1_patient_meta_selected.sampleID.isin(CDH1)]
cdh1_patient_meta_sorted = pd.concat([cdh1_patient_mt,cdh1_patient_wt],axis=0).set_index("sampleID")
cdh1_patient_meta_sorted.to_csv("./CDH1_patient_meta.csv",index=False)

rhoa_patient_mt = rhoa_patient_meta_selected[rhoa_patient_meta_selected.sampleID.isin(RHOA)]
rhoa_patient_wt = rhoa_patient_meta_selected[~rhoa_patient_meta_selected.sampleID.isin(RHOA)]
rhoa_patient_meta_sorted = pd.concat([rhoa_patient_mt,rhoa_patient_wt],axis=0).set_index("sampleID")
rhoa_patient_meta_sorted.to_csv("./RHOA_patient_meta.csv",index=False)