#-------------------------------------------------------------------------------------------------
# 05_Clustering_and_determining_signature_molecules
#-------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------#
### Packages installation ###
#----------------------------------------------------------------------------#


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

library(dplyr)

BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install("preprocessCore")
library(preprocessCore)

BiocManager::install("CancerSubtypes")
library(CancerSubtypes)
library(sigclust)

# readr package installation
install.packages("readr")
library(readr)

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

library(stats)

install.packages("circlize")
library(circlize)

install.packages("colorRamp2")
library(colorRamp2)

#-----------------------------------------------------------------------------#
##### Part 1. Clustering #####
#-----------------------------------------------------------------------------#


### 1. Loading dataset to R ###


# Input: 각 환자별 T/N glycopeptide (gp) log2 fold change 정보가 있는 csv 파일 (외부 파일)
# Output: 해당 정보가 정리된 matrix (protein_matrix)


# make dataset matrix of proteomics
# 파일 경로 지정
file_path <- "./03_Analysis/00_DB/Proteomics/glycopeptide_change.csv"

# 행렬로 데이터 불러오기
data_frame <- read.csv(file_path, header = TRUE, stringsAsFactors=FALSE)
protein_matrix <- as.matrix(data_frame)


#-----------------------------------------------------------------------------#.


### 2. 50% 이상의 환자에서 결측치가 아닌 gp만 남기고 삭제 / quantile normalization ###


# Input: 각 환자별 T/N gp log2 fold change 정보가 있는 matrix (protein_matrix)  
# Output: 50% 이상의 환자에서 결측치가 아닌 gp 정보 > quantile normalized matrix (norm_p50_matrix)
# Tool: preprocessCore package


BiocManager::install("preprocessCore")
library(preprocessCore)

# 각 행(gp)에 대해 NA가 아닌 값이 50% 이상인 행(gp)만 추출
p50_matrix <- protein_matrix[rowMeans(!is.na(protein_matrix[, 3:82])) >= 0.5, ]

# Quantile normalization은 전체 데이터를 모두 일괄로 처리해야함

# gp 정보 빼고 발현량만 추출 후 numeric으로 변환
p50 <- p50_matrix[, 3:82]
p50 <- apply(p50, 2, as.numeric)

# Quantile normalization 수행
norm_p50 <- normalize.quantiles(p50)

# normalization 전/후 비교
data.checkDistribution(p50)
data.checkDistribution(norm_p50)

# 다시 protein 정보 합치기
colnames(norm_p50) <- colnames(p50_matrix)[3:82]
norm_p50_matrix <- cbind(p50_matrix[, 1:2], norm_p50)


#-----------------------------------------------------------------------------#

### Option: 결측치 처리 - mice package 이용 ###

# mice 패키지 설치
install.packages("mice")
library(mice)

# norm_p50_matrix의 결측치 처리를 위한 데이터 준비
norm_p50_df <- as.data.frame(norm_p50)
# mice를 사용하여 결측치 처리
imputed_data <- mice(norm_p50_df, m = 5, maxit = 50, meth = 'pmm', seed = 123)

# 대체된 데이터 추출
complete_norm_p50_df <- complete(imputed_data)


#-----------------------------------------------------------------------------#


### 3. Median absolute deviation (MAD) 기준으로 상위 10%, 20%, 30%  해당하는 gp 추출 ###


# Input: 50% 이상의 환자에서 결측치가 아닌 gp 정보 > quantile normalized matrix (norm_p50_matrix)
# Output: MAD값 기준 상위 10/20/30%에 해당하는 gp 정보만 담긴 matrix, gp name은 제외 (M10/M20/M30)


# 각 gp의 MAD 계산
p50_mad <- apply(norm_p50, 1, mad, na.rm=TRUE)

# 상위 10%, 20%, 30%에 해당하는 MAD 값을 계산
top_10_percent <- quantile(p50_mad, probs = 0.9, na.rm=TRUE)
top_20_percent <- quantile(p50_mad, probs = 0.8, na.rm=TRUE)
top_30_percent <- quantile(p50_mad, probs = 0.7, na.rm=TRUE)

# 각 기준에 따라 gp 선택 (gp 이름 포함)
selected_ptn_10 <- norm_p50_matrix[p50_mad >= top_10_percent, ]
selected_ptn_20 <- norm_p50_matrix[p50_mad >= top_20_percent, ]
selected_ptn_30 <- norm_p50_matrix[p50_mad >= top_30_percent, ]

# ExecuteCNMF에 사용할 수 있는 데이터 포맷으로 변경: matrix - 각 행에 gp, 각 열에 환자
M10 <- norm_p50[p50_mad >= top_10_percent, ]
M20 <- norm_p50[p50_mad >= top_20_percent, ]
M30 <- norm_p50[p50_mad >= top_30_percent, ]

# matrix에서 NA를 0으로 대체: CancerSubtypes의 CNMF clustering 수행시 결측치가 있으면 분석 불가
M10[is.na(M10)] <- 0
M20[is.na(M20)] <- 0
M30[is.na(M30)] <- 0


#-----------------------------------------------------------------------------#


### 4. Non-negative 치환: 양수 -> 0 / 음수 -> 0  두 matrix로 나눈 후 리스트화 ###


# Input: MAD값 기준 상위 10/20/30%에 해당하는 gp 정보만 담긴 matrix (M10/M20/M30)
# Output: log2(TPM+1) fold change 값의 양수 or 음수만 담고 있는 두 개의 matrix를 포함한 리스트 (List_10/20/30)


# 음수 값을 0으로 바꾸어주는 함수
replace_negative <- function(x) {
  x[x < 0] <- 0
  return(x)
}

# 1. 음수 값을 0으로 바꾼 매트릭스 = 양수값만 남아있는 매트릭스
positive_10 <- replace_negative(M10)
positive_20 <- replace_negative(M20)
positive_30 <- replace_negative(M30)

# 2. 양수 값을 0으로 바꾼 매트릭스 = 음수값만 남아있는 매트릭스
negative_10 <- replace_negative(-M10)
negative_20 <- replace_negative(-M20)
negative_30 <- replace_negative(-M30)

# 두 매트릭스를 리스트로 묶기
List_10 <- list(positive_10, negative_10)
List_20 <- list(positive_20, negative_20)
List_30 <- list(positive_30, negative_30)


#-----------------------------------------------------------------------------#


### 5. CNMF clustering 실행 ###


# Input: log2(TPM+1) fold change 값의 양수 or 음수만 담고 있는 두 개의 matrix를 포함한 리스트 (List_10/20/30)
# Output: Clustering 결과 (group, distance, consensus, cophenetic cor), silhouette width plot (sil_plot)
# Tool: CancerSubtypes, ComplexHeatmap package


library(CancerSubtypes)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# ExecuteCNMF 함수 실행
# k=2~6, nrun = 100 이상
# MAD 10/20/30%, k 개수별로 클러스터링을 수행하고, 결과로부터 사용할 MAD 기준과 cluster 수 결정

# 결과를 저장할 리스트 생성
result10_list <- list()
result20_list <- list()
result30_list <- list()


# k=2부터 6까지 반복 실행 (MAD 10%, 20%, 30% 이상 gp)
for (k in 2:6) {
  result <- ExecuteCNMF(List_10, k, nrun = 100)
  result10_list[[as.character(k)]] <- result
}
names(result10_list) <- paste0("result10_", 2:6)

for (k in 2:6) {
  result <- ExecuteCNMF(List_20, k, nrun = 100)
  result20_list[[as.character(k)]] <- result
}
names(result20_list) <- paste0("result20_", 2:6)

for (k in 2:6) {
  result <- ExecuteCNMF(List_30, k, nrun = 100)
  result30_list[[as.character(k)]] <- result
}
names(result30_list) <- paste0("result30_", 2:6)


# 위 결과를 통해 MAD 20% 이상 protein들에 대해 cluster 2개로 구분하는 것으로 결정

# Group 정보 / Distance 정보 / Silhouette width / Consensus 정보 저장
group20_2 = result20_list[["result20_2"]][["group"]] #그룹핑 정보
distanceMatrix20_2=result20_list[["result20_2"]][["distanceMatrix"]] #Distance 정보
silhouette20_2=silhouette_SimilarityMatrix(group20_2, distanceMatrix20_2) #Silhouette score
Consensus20_2=result20_list[["result20_2"]][["originalResult"]]@consensus #Consensus 정보
plot(silhouette20_2)

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# 간단히 Consensus matrix plotting
Heatmap(Consensus20_2,name="Consensus value", show_column_names = FALSE, show_row_names = FALSE)

# Silhouette width plot
sil_plot <- plot(silhouette20_2, col="darkorange3")

# Cophenetic correlation 계산
d1 <- dist(distanceMatrix20_2)
hclust_obj <- hclust(d1, "ave")
d2 <- cophenetic(hclust_obj)
cor(d1, d2)

# Dendrogram plotting 해보기
plot(hclust_obj)


#-----------------------------------------------------------------------------#
##### Part 2. Determining signature molecules of each subtype #####
#-----------------------------------------------------------------------------#


### 1. Consensus가 잘 되어있는 환자만 filtering: protein expression 데이터 추출 ###


# Input: silhouette positive patients 정보 (positive_sil_matrix), 전체 gp에 대한log2(TPM+1) fold change 정보 (f2c_analysis)
# Output: 각 환자 cluster/non-cluster별 분리된 gp expression 데이터 (f2c_pts_1~2/f2c_pts_not1~2)

# Silhouette score가 양수인 환재 샘플 = 해당 subtype에 consensus가 잘 되어있는 환자만 분석에 사용

# silhouette 정보 matrix 생성
sil_matrix <- as.matrix(silhouette20_2[, 1:3])

# 환자정보 추출
pts <- colnames(M20)

# sil_matrix의 1열에 pts 추가
sil_matrix <- cbind(pts = pts, sil_matrix)

# silhouette width가 양수인 행을 선택
positive_sil_matrix <- sil_matrix[sil_matrix[, "sil_width"] > 0, ]

# 실루엣 양수인 환자들(pts_core)의 clustering 정보 extract  
clusters <- positive_sil_matrix[, "cluster"]
pts_core <- positive_sil_matrix[, "pts"]

# 각 그룹별 pts_core 클러스터링 정보
cluster_1_pts <- positive_sil_matrix[clusters == 1, ]
pts_1 <- cluster_1_pts[, "pts"]
cluster_2_pts <- positive_sil_matrix[clusters == 2, ]
pts_2 <- cluster_2_pts[, "pts"]

# 환자 cluster별 gp expression data 분리
f2c <- norm_p50_matrix #gp name 용
f2c_analysis <- norm_p50
f2c_analysis <- apply(f2c_analysis, 2, as.numeric)
# core 환자 샘플만 따로 expression matrix 분리
f2c_pts_core <- f2c_analysis[, colnames(f2c_analysis) %in% pts_core]

# Group / Non-group별 expression 정보 분리
f2c_pts_1 <- f2c_pts_core[, colnames(f2c_pts_core) %in% pts_1]
f2c_pts_not1 <- f2c_pts_core[, !colnames(f2c_pts_core) %in% pts_1]
f2c_pts_2 <- f2c_pts_core[, colnames(f2c_pts_core) %in% pts_2]
f2c_pts_not2 <- f2c_pts_core[, !colnames(f2c_pts_core) %in% pts_2]


#-----------------------------------------------------------------------------#


### 2. Subtype간 t-test 수행하고 보정 ###


# Input: 각 환자 cluster/non-cluster별 분리된 gp expression 데이터 (f2c_pts_1~2/f2c_pts_not1~2)
# Output: t_test 결과 with adjusted P value (t_test_group1~2)
# Tool: stats package


library(stats)

t_test_group1 <- data.frame(ptn = character(), P_Value = numeric(), stringsAsFactors = FALSE)
t_test_group2 <- data.frame(ptn = character(), P_Value = numeric(), stringsAsFactors = FALSE)


# 각 gp에 대해 group vs non-group t-test 수행 (3987 = 사용할 전체 pp 수)

for (i in 1:3987) {
  ptn <- i
  
  # 그룹 내 모든 환자에 대해 t-test 수행
  n_group <- sum(!is.na(f2c_pts_1[i, ]))
  n_not_group <- sum(!is.na(f2c_pts_not1[i, ]))
  
  # 충분한 데이터가 있는 경우에만 t-test 수행하도록 설정 (결측치가 너무 많으면 수행 불가)
  if (n_group >= 2 & n_not_group >= 2) {
    t_result_group1 <- t.test(f2c_pts_1[i, ], f2c_pts_not1[i, ], na.action = na.omit) 
    
    # 결과 저장
    t_test_group1 <- rbind(t_test_group1, data.frame(ptn = ptn, P_Value = t_result_group1$p.value))
  }
}

for (i in 1:3987) {
  ptn <- i
  
  # 그룹 내 모든 환자에 대해 t-test 수행
  n_group <- sum(!is.na(f2c_pts_2[i, ]))
  n_not_group <- sum(!is.na(f2c_pts_not2[i, ]))
  
  # 충분한 데이터가 있는 경우에만 t-test 수행하도록 설정 (결측치가 너무 많으면 수행 불가)
  if (n_group >= 2 & n_not_group >= 2) {
    t_result_group2 <- t.test(f2c_pts_2[i, ], f2c_pts_not2[i, ], na.action = na.omit)  
    
    # 결과 저장
    t_test_group2 <- rbind(t_test_group2, data.frame(ptn = ptn, P_Value = t_result_group2$p.value))
  }
}

# Bonferroni 보정을 적용하여 유의수준 보정
t_test_group1$Adjusted_P_Value <- p.adjust(t_test_group1$P_Value, method = "bonferroni")
t_test_group2$Adjusted_P_Value <- p.adjust(t_test_group2$P_Value, method = "bonferroni")

# 결과 brief하게 확인: 보정된 t-test 결과가 훨씬 적은 DEP 산출
hist(t_test_group1$P_Value, 100)
hist(t_test_group1$Adjusted_P_Value, 100)
hist(t_test_group2$P_Value, 100)
hist(t_test_group2$Adjusted_P_Value, 100)


#-----------------------------------------------------------------------------#


### 3. Adjusted P value와 median 값 기준으로 각 그룹별 발현되는 pp 결정 ###


# Input: t_test 결과 with adjusted P value (t_test_group1~2)
# Output: adjusted P value와 median 기준으로 필터링된 signature gp 정보 (f2c_DEG_med_combined)
# Tool: dplyr, CancerSubtypes package


library(dplyr)
library(CancerSubtypes)

# Adjusted P Value가 0.05 이하인 gp 추출
sig_ptn_group1 <- t_test_group1[t_test_group1$Adjusted_P_Value <= 0.05, ]
sig_ptn_group2 <- t_test_group2[t_test_group2$Adjusted_P_Value <= 0.05, ]

# significant_gp에 해당되는 gp expression data만 추출
f2c_DEP_name1 <- f2c[sig_ptn_group1$ptn, ]
f2c_DEP_name2 <- f2c[sig_ptn_group2$ptn, ]

# Core pts의 significant_gp에 해당되는 gp expression data만 추출: name 미포함
f2c_DEP_group1 <- f2c_pts_core[sig_ptn_group1$ptn, ]
f2c_DEP_group2 <- f2c_pts_core[sig_ptn_group2$ptn, ]

# Group/non-group별 significant gp expression data 분리
f2c_DEP_pts_1 <- f2c_DEP_group1[ ,colnames(f2c_pts_core) %in% pts_1]
f2c_DEP_pts_not1 <- f2c_DEP_group1[ ,!colnames(f2c_pts_core) %in% pts_1]
f2c_DEP_pts_2 <- f2c_DEP_group2[ ,colnames(f2c_pts_core) %in% pts_2]
f2c_DEP_pts_not2 <- f2c_DEP_group2[ ,!colnames(f2c_pts_core) %in% pts_2]

# 각 그룹에서 gp별 median 값 추출
median_vals_pts_1 <- apply(f2c_DEP_pts_1, 1, median, na.rm=TRUE)
median_vals_pts_not1 <- apply(f2c_DEP_pts_not1, 1, median, na.rm=TRUE)
median_vals_pts_2 <- apply(f2c_DEP_pts_2, 1, median, na.rm=TRUE)
median_vals_pts_not2 <- apply(f2c_DEP_pts_not2, 1, median, na.rm=TRUE)

# 해당 median 값만 비교하기 위해 f2c_median 매트릭스 생성
g1_median <- cbind(ptn = f2c_DEP_name1[ ,1:2], g1_median = median_vals_pts_1, not1_median = median_vals_pts_not1)
g2_median <- cbind(ptn = f2c_DEP_name2[ ,1:2], g2_median = median_vals_pts_2, not2_median = median_vals_pts_not2)

# 데이터프레임 변환
g1_median_df <- as.data.frame(g1_median)
g2_median_df <- as.data.frame(g2_median)

# Group의 median 정보를 numeric으로 변환
g1_median_df$g1_median <- as.numeric(as.character(g1_median_df$g1_median))
g1_median_df$not1_median <- as.numeric(as.character(g1_median_df$not1_median))
g2_median_df$g2_median <- as.numeric(as.character(g2_median_df$g2_median))
g2_median_df$not2_median <- as.numeric(as.character(g2_median_df$not2_median))

# dplyr 패키지 로드
library(dplyr)

# 각 그룹에서 median 값이 0보다 크고, 타 그룹에서 0보다 작으며, 타 그룹보다 큰 gp 선별
filtered_median_g1 <- g1_median_df %>%
  filter(g1_median > 0 & not1_median < 0 & g1_median > not1_median)
filtered_median_g2 <- g2_median_df %>%
  filter(g2_median > 0 & not2_median < 0 & g2_median > not2_median)

# gp name이 #N/A 또는 공백이 아닌 행만 선택
filtered_median_g1 <- filtered_median_g1[filtered_median_g1$Peptide !="#N/A" & filtered_median_g1$Peptide != "", , drop = FALSE]
filtered_median_g2 <- filtered_median_g2[filtered_median_g2$Peptide !="#N/A" & filtered_median_g2$Peptide != "", , drop = FALSE]

# 각 filtered_median_g1~2의 1열의 gp seq 추출
ptn_vector_p1 <- as.character(unique(filtered_median_g1$Peptide))
ptn_vector_p2 <- as.character(unique(filtered_median_g2$Peptide))

# silhouette 양수인 환자에서 DEP 기준을 정했으나, 다시 모든 환자에서 적용해야함.
# 전체 clustering 정보 extract  
clusters <- sil_matrix[, "cluster"]
pts_all <- sil_matrix[, "pts"]

cluster_1_pts <- sil_matrix[clusters == 1, ]
pts_1 <- cluster_1_pts[, "pts"]
cluster_2_pts <- sil_matrix[clusters == 2, ]
pts_2 <- cluster_2_pts[, "pts"]

# f2c에서 Peptide가 ptn_vector와 일치하는 행 선택
ptn_DEG_med_p1 <- f2c[f2c[, 1] %in% ptn_vector_p1, ]
ptn_DEG_med_p2 <- f2c[f2c[, 1] %in% ptn_vector_p2, ]

# Heatmap용 matrix 생성
# 이 때, 환자수가 많은 그룹부터 배열하여 consensus 맵 배열과 일치시키자
# 환자수가 많은 순으로 정렬
ptn_DEG_med_p21 <- rbind(ptn_DEG_med_p2, ptn_DEG_med_p1)

f2c_DEG_med_pts1 <- ptn_DEG_med_p21[,colnames(ptn_DEG_med_p21) %in% pts_1]
f2c_DEG_med_pts2 <- ptn_DEG_med_p21[,colnames(ptn_DEG_med_p21) %in% pts_2]

f2c_DEG_med_combined <- cbind(f2c_DEG_med_pts2, f2c_DEG_med_pts1)
f2c_DEG_med_combined <- apply(f2c_DEG_med_combined, 2, as.numeric)

# NA를 0으로
f2c_DEG_med_combined[is.na(f2c_DEG_med_combined)] <- 0

# Glycopeptide: pts_2 = group 1, pts_1 = group 2

# 결과 확인
data.checkDistribution(f2c_DEG_med_combined)


#-----------------------------------------------------------------------------#
##### Part 3. Visualization by heatmap #####
#-----------------------------------------------------------------------------#


### 1. Heatmap annotation용 데이터 저장 ###


# Input: 필터링된 signature gp 정보 (f2c_DEG_med_combined), 그룹 정보 (sil_matrix$cluster)
# Output: Heatmap 그룹핑에 사용할 annotation 정보 (csv 파일)
# 1행에는 Cluster 기입


# Cluster
write.csv(colnames(f2c_DEG_med_combined), "glyco_subgroup_cluster.csv", row.names = TRUE)
# 파일 추출 후에 group별 환자수만큼 1열에 GLYCO1, GLYCO2 기입

# Signature gp: Signature molecule의 row annotation에 사용
write.csv(f2c_DEG_med_combined[,1], "glyco_subgroup_signature.csv", row.names = TRUE)
# 파일 추출 후에 group별 signature gp 수만큼 Glyco1, Glyco2 기입
# 1행에는 Signature 기입

# Consensus
write.csv(sil_matrix[,1:2], "glyco_subgroup_consensus.csv", row.names = TRUE)
# 파일 추출 후에 group별 환자에 맞게 1열에 GLYCO1, GLYCO2 기입
# 1행에는 Consensus 기입


#-----------------------------------------------------------------------------#


### 2. Heatmap 그리기 ###


# Input: signature gp 정보 (f2c_DEG_med_combined), consensus 정보 (Consensus20_2), annotation 파일들
# Output: Signature gp heatmap (heatmap1), Consensus heatmap (heatmap2) 
# Tool: ComplexHeatmap, circlize, colorRamp2 package


library(ComplexHeatmap)
library(circlize)
library(colorRamp2)

# Signature molecules heatmap

# Color 지정
col_fun = colorRamp2(c(-2, 0, 2), c("green", "black", "red"))

# Annotation 지정 > 위에서 셋팅했던 csv 파일로부터 불러옴
# Clustering annotation (column)
anno_C = read.csv("./03_Analysis/00_DB/RNAseq/Clustering/Subgroup/glyco_subgroup_cluster.csv")
ha_C=HeatmapAnnotation(df = anno_C, which="col", col = list(
  Cluster = c("GLYCO1" = "blueviolet", "GLYCO2" = "coral1")),
  show_annotation_name = FALSE, show_legend = FALSE
)

# Signature annotation (row)
anno_S = read.csv("./03_Analysis/00_DB/RNAseq/Clustering/Subgroup/glyco_subgroup_signature.csv")
ha_S=HeatmapAnnotation(df = anno_S, which="row", col = list(
  Signature = c("Glyco1" = "blueviolet", "Glyco2" = "coral1")),
  show_annotation_name = FALSE, show_legend = FALSE
)


# Heatmap
heatmap1 = Heatmap(f2c_DEG_med_combined,
                   show_row_names = FALSE,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   col=col_fun,
                   show_column_names=FALSE,
                   right_annotation = ha_S,
                   bottom_annotation = ha_C,
                   show_heatmap_legend=FALSE)
# Legend 지정
legend1 = Legend(direction="horizontal", col_fun = col_fun, title = "log2(fold-change)",
                 title_position = "topcenter",legend_width = unit(3, "cm"))
# Heatmap + Legend plotting
plot(heatmap1, annotation_legend_list=legend1, annotation_legend_side="bottom")


# Consensus heatmap

# Color 지정
col_fun2 = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

# Annotation 지정 > 위에서 셋팅했던 csv 파일로부터 불러옴
# Consensus annotation (column)
anno_A = read.csv("./03_Analysis/00_DB/RNAseq/Clustering/Subgroup/glyco_subgroup_consensus.csv")
ha_A=HeatmapAnnotation(df = anno_A, which="col", col = list(
  Consensus = c("GLYCO1" = "blueviolet", "GLYCO2" = "coral1")),
  show_annotation_name = FALSE, show_legend = FALSE
)

# Heatmap
heatmap2 = Heatmap(Consensus20_2,name="Consensus value",
                   show_column_names = FALSE,
                   show_row_names = FALSE,
                   bottom_annotation = ha_A,
                   show_heatmap_legend=FALSE
)
# Legend 지정
legend2 = Legend(direction="horizontal", col_fun = col_fun2, break_dist = 0.5, title = "Consensus value",
                 title_position = "topcenter",legend_width = unit(3, "cm"))
# Heatmap + Legend plotting
plot(heatmap2, annotation_legend_list=legend2, annotation_legend_side = "bottom")


#-----------------------------------------------------------------------------#
##### Part 4. Integrative clustering을 위한 indicative matrix (0 또는 1) 생성 및 저장 #####
#-----------------------------------------------------------------------------#


# Input: 각 환자들의 group 정보 (sil_matrix)
# Output: group에 속한 여부가 1 또는 0으로 표현된 indicative matrix (indi_glyco.csv) 


# sil_matrix transpose
t_sil_matrix <- t(sil_matrix)

# 1행을 rownames로 설정
colnames(t_sil_matrix) <- t_sil_matrix[1,]

# 1,4행 제거
t_sil_matrix <- t_sil_matrix[-c(1,4), ]

# 1행의 내용을 2행에 복사
t_sil_matrix[2, ] <- t_sil_matrix[1, ]

# 1행의 데이터는 1은 1로, 1이 아니면 모두 0으로
t_sil_matrix[1, ] <- ifelse(t_sil_matrix[1, ] == 1, 1, 0)

# 2행의 데이터는 2는 1로, 2가 아니면 모두 0으로
t_sil_matrix[2, ] <- ifelse(t_sil_matrix[2, ] == 2, 1, 0)

# numeric 치환
t_sil_matrix <- apply(t_sil_matrix, 2, as.numeric)

# 저장할 이름 설정 및 행 이름 변경 (subtype 이름별)
indi_glyco <- t_sil_matrix
rownames(indi_glyco)[1:2] <- c("Glyco2", "Glyco1")

# indicator matrix를 csv로 저장
write.csv(indi_glyco, "indi_glyco.csv", row.names = TRUE)



