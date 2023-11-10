### 0. Packages installation ###

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

install.packages("circlize")
library(circlize)

install.packages("colorRamp2")
library(colorRamp2)

#-----------------------------------------------------------------------------#
##### Part 5. Integrative clustering #####
#-----------------------------------------------------------------------------#


### 1. Loading dataset to the R and merged indicative matrix 생성###


# Input: 각 molecule별 그룹핑되어 저장해둔 indicative matrix 파일 (csv 파일들)
# Output: Merged indicative matrix (indi_matrix)


# 각 dataset별 indicative matrix 로딩

# Indicative matrix들을 저장해둔 폴더 경로설정
data_dir <- "C:/Users/admin/Desktop/LAIDD/RNAseq/Clustering/MAD filtering data"

# 각 molecule subtype별 indicative matrix 불러오기
indi_rna <- as.matrix(read.csv(file.path(data_dir, "indi_rna.csv"), row.names=1))
indi_ptn <- as.matrix(read.csv(file.path(data_dir, "indi_ptn.csv"), row.names=1))
indi_phos <- as.matrix(read.csv(file.path(data_dir, "indi_phos.csv"), row.names=1))
indi_glyco <- as.matrix(read.csv(file.path(data_dir, "indi_glyco.csv"), row.names=1))

# Part4 과정에서 환자샘플 네이밍을 일치시켰음. 열에 배열된 환자샘플을 순서대로 배열 
# 먼저 이름순으로 정렬하여 벡터에 저장
sorted_cols_rna <- colnames(indi_rna)[order(colnames(indi_rna))]
sorted_cols_ptn <- colnames(indi_ptn)[order(colnames(indi_ptn))]
sorted_cols_phos <- colnames(indi_phos)[order(colnames(indi_phos))]
sorted_cols_glyco <- colnames(indi_glyco)[order(colnames(indi_glyco))]

# 각 벡터의 순서가 일치하는지 확인
identical(sorted_cols_rna, sorted_cols_ptn)
identical(sorted_cols_rna, sorted_cols_phos)
identical(sorted_cols_rna, sorted_cols_glyco)

# 각 indicative matrix의 열 순서를 통일
indi_rna <- indi_rna[, sorted_cols_rna]
indi_ptn <- indi_ptn[, sorted_cols_rna]
indi_phos <- indi_phos[, sorted_cols_rna]
indi_glyco <- indi_glyco[, sorted_cols_rna]

# 모든 matrix의 열은 동일함 -> 행을 합쳐 Merged indicated matrix 생성
indi_matrix <- rbind(indi_rna, indi_ptn, indi_phos, indi_glyco)


#-----------------------------------------------------------------------------#


### 2. Integrative clustering using ConsensusClusterPlus###


# Input: merged indicative matrix (indi_matrix)
# Output: integrative clustering result (results), figures (title 디렉토리)
# Tool: ConsensusClusterPlus package


library(ConsensusClusterPlus)

# 저장될 경로 설정
title=tempdir()

# integrative matrix로 consensus clustering 실행
# k=2~8
# 1000회 반복
# 전체 molecule 특성의 80% 사용
# hierarchical, dissimilarity measure: pearson correlation
results = ConsensusClusterPlus(indi_matrix,maxK=8,reps=1000,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",plot="png")

# Figure들이 저장된 directory 찾아가서 확인
title

# Figure를 통해 최종 3개의 subtype으로 구분하는 것으로 결정


#-----------------------------------------------------------------------------#


### 3. Subtype 정보와 meta data 포함하여 heatmap 그리기###


# Input: 정렬된 indicative matrix (Class_arranged), meta data (clinical_info.csv)
# Output: Integrative clustering heatmap (iheatmap)
# Tool: ComplexHeatmap, circlize, colorRamp2 package


library(ComplexHeatmap)
library(circlize)
library(colorRamp2)

# indi_matrix 정보 마지막 열에 class 정보 추가
Class <- as.matrix(t(results[[3]][["consensusClass"]]))
Class <- rbind(indi_matrix, Class[1,])
rownames(Class)[12] <- "class"
Class <- t(Class)

# Class 정보 추출 
write.csv(Class, "CCP.csv", row.names = TRUE)

# 추출 후 엑셀에서 integrative subtype별로 필터링 후 CCP_arranged.csv로 저장, 다시 불러오기
Class_arranged = as.matrix(read.csv("C:/Users/admin/Desktop/LAIDD/RNAseq/Clustering/CCP_arranged.csv", header=T, row.names=1))
Class_arranged = t(Class_arranged)

# Integrative subtype 정보만 따로 빼놓고, 열이름 class로 지정
Class_anno <- as.data.frame(Class_arranged[12,])
colnames(Class_anno) <- "class" 

# 원래 Class_arrange에서는 subtype 정보는 삭제 (heatmap에 사용될 matrix)
Class_arranged = Class_arranged[-12,]

# meta data 정리: Meta data는 clinical_info.csv에 따로 있어 불러오기
meta <- read.csv("C:/Users/admin/Desktop/LAIDD/RNAseq/Clustering/Subgroup/meta/clinical_info.csv",
           header=T, row.names=1)

# subtype별 환자 순서대로 정렬
class_order <- colnames(Class_arranged)
meta_arranged <- meta[class_order, ]

# 정렬된 meta data csv로 저장 (annotation용)
write.csv(meta_arranged, "meta_final.csv")
# 엑셀에서 1열 삭제

# Color 지정
col_fun = colorRamp2(c(0, 1), c("black", "red"))

# Annotation 지정 > 위에서 meta data 최종 셋팅했던 csv 파일로부터 불러옴
# Meta data (ha1)
anno = read.csv("C:/Users/admin/Desktop/LAIDD/RNAseq/Clustering/Subgroup/meta/meta_final.csv")
ha1=HeatmapAnnotation(df = anno, col = list(
  Histology = c("Diffuse" = "ivory", "Intestinal" = "darkorchid4", "Mixed" = "darkorchid1", "Others" = "gray"),
  MSI = c("MSI-H" = "darkgreen", "MSS/MSI-L" = "ivory"),
  EBV = c("EBV-" = "ivory", "EBV+" = "darkorange2"),
  Gender = c("F" = "pink", "M" = "skyblue"),
  pStage = c("I" = "azure3", "II" = "khaki", "III" = "coral3", "IV" = "firebrick4")),
  show_annotation_name = TRUE, annotation_name_offset = unit(2, "mm"), border=T,
annotation_name_side = "left") 

# Subtype annotation (ha2)
ha2=HeatmapAnnotation(df = Class_anno, which="col", col = list(class=
  c("1" = "dodgerblue3","2" = "seagreen2", "3" = "lightsalmon")),
  show_annotation_name = FALSE, show_legend = FALSE
)

# Heatmap
iheatmap <- Heatmap(Class_arranged,
                    cluster_columns = FALSE,
                    show_column_names=FALSE,
                    clustering_distance_rows = "pearson",
                    top_annotation = ha2,
                    bottom_annotation = ha1,
                    col=col_fun,
                    show_heatmap_legend=FALSE)

plot(iheatmap)
