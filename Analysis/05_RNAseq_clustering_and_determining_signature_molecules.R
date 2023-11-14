### 0. Packages installation ###


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

install.packages("readr")
library(readr)

library(dplyr)

BiocManager::install("DESeq2")
library(DESeq2)

install.packages("httr")
library(httr)

install.packages("jsonlite")
library(jsonlite)

BiocManager::install("biomaRt")
library(biomaRt)

BiocManager::install("preprocessCore")
library(preprocessCore)

BiocManager::install("CancerSubtypes")
library(CancerSubtypes)
library(sigclust)

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


# Input: 각 환자별 normal and tumor RNAseq processed data (외부 파일)
# Output: 해당 정보가 정리된 리스트 (data_list)
# Tool: readr package


install.packages("readr")
library(readr)

# RNA-seq 데이터: 각 환자마다 normal, tumor 파일 따로 있고, TPM 정보 있음
# make dataset tables of 160 pts RNAseq data at once
data_dir <- "C:/Users/admin/Desktop/LAIDD/RNAseq/Processed RNAseq data/RNAseq_RSEM"   # data가 있는 폴더 directory 입력

# 파일 이름마다 "results"가 들어가있어 공통 패턴으로 지정
data_files <- list.files(data_dir, pattern = "results")

# make empty list
data_list <- list()

# 데이터 파일을 읽어와서 리스트의 항목으로 추가
for (file in data_files) {
  file_path <- file.path(data_dir, file)
  data <- read_delim(file_path, delim = "\t")  # tab 구분자 데이터 파일
  
  # 파일 이름에서 필요없는 부분 제거
  data_name <- gsub("_rsem_genes_original_results", "", file)
  
  # 리스트의 항목으로 데이터프레임 추가
  data_list[[data_name]] <- data
}


#-----------------------------------------------------------------------------#


### 2. 50% 이상의 환자 데이터에서 TPM 1 이상인 gene만 남기고 삭제 ###


# Input: 모든 환자와 RNA expression 정보가 담긴 데이터프레임 (combined_data)  
# Output: 50% 이상의 환자에서 TPM 1 이상인 RNA expression 정보가 담긴 데이터프레임 리스트 (data_list)
# Tool: dplyr package


library(dplyr)

# 모든 데이터프레임을 하나로 합치기
combined_data <- bind_rows(data_list, .id = "Dataset")

# gene별로 50% 이상 환자에서 TPM이 1 이상인 행만 필터링
filtered_data <- combined_data %>%
  group_by(gene_id) %>%
  filter(mean(TPM >= 1, na.rm = TRUE) >= 0.5)

# 필터링된 데이터를 다시 각각의 데이터프레임으로 나누기
cleaned_data_list <- split(filtered_data, filtered_data$Dataset)
data_list <- cleaned_data_list


#-----------------------------------------------------------------------------#


### 3. Calculate Log2(TPM+1) ###


# Input: 50% 이상의 환자에서 TPM 1 이상인 RNA expression 정보가 담긴 데이터프레임 리스트 (data_list)
# Output: log2(TPM+1)을 구하고 gene_id와 log2(TPM+1)만 남긴 데이터프레임 리스트 (data_list)
# Tool: DESeq2 package


BiocManager::install("DESeq2")
library(DESeq2)

# log2(TPM+1) 변환 함수 정의
log2_transform <- function(data) {
  return(log2(data + 1))
}


# 모든 TPM 데이터를 log2(TPM+1) 으로 변환
for (i in 1:length(data_list)) {
  data_list[[i]]$TPM <- log2_transform(data_list[[i]]$TPM)
}

# 열 이름 변경
for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[colnames(data_list[[i]]) == "TPM"] <- "log2(TPM+1)"
}


# 각 데이터셋에서 "gene_id"와 "log2(TPM+1)" 열을 제외한 열을 삭제하고 새로운 데이터프레임 생성
for (i in 1:length(data_list)) {
  data_list[[i]] <- data_list[[i]][, c("gene_id", "log2(TPM+1)")]
}


#-----------------------------------------------------------------------------#


### 4. gene_id를 기반으로 gene_name 불러오기 (Ensembl) ###


# Input: gene_id와 log2(TPM+1)만 남긴 데이터프레임 리스트 (data_list)
# Output: gene_id, gene_name, log2(TPM+1) 열로 구성된 데이터프레임 리스트 (data_list_gene_names)
# Tool: httr, jsonlite, BiomaRt package


install.packages("httr")
install.packages("jsonlite")
library(httr)
library(jsonlite)
BiocManager::install("biomaRt")
library(biomaRt)

# gene_id 추출을 위해 하나의 데이터 (111N) 추출
data_111N <- data_list[[1]]  # 111N 데이터 추출 (예시로 첫 번째 데이터를 사용)

# gene_id 열 추출
gene_ids <- data_111N$gene_id

# gene_ids를 문자열 벡터로 변환
gene_ids <- as.character(gene_ids)
head(gene_ids)

# Ensembl human gene 데이터베이스에 연결하기
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# gene_id를 gene name으로 변환
gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", 
                    values = gene_ids, 
                    mart = ensembl)

# gene_names 데이터프레임에 NA로 표시할 열 추가
gene_names$external_gene_name[is.na(gene_names$external_gene_name)] <- "NA"

# gene_name 열 추가된 데이터리스트 생성
data_list_gene_names <- lapply(data_list, function(data) {
  data <- data.frame(data[, 1], gene_name = gene_names$external_gene_name[match(data$gene_id, gene_names$ensembl_gene_id)], data[, 2])
  colnames(data) <- c("gene_id", "gene_name", "log2_TPM+1")
  return(data)
})


#-----------------------------------------------------------------------------#


### 5. Quantile normalization of Log2(TPM+1) ###
# Quantile normalization은 전체 데이터를 모두 일괄로 처리해야함


# Input: 모든 환자의 log2(PTM+1) 정보가 담긴 하나의 matrix (TPM_matrix)
# Output: Quantile normalized 데이터프레임 리스트 (norm_data_list_final)
# Tool: preprocessCore package


BiocManager::install("preprocessCore")
library(preprocessCore)

# 데이터프레임의 이름 유지를 위한 이름 리스트 생성
data_names <- names(data_list_gene_names)

# 모든 환자 데이터프레임에서 log2(TPM+1) 열 추출
TPM_list <- lapply(data_list_gene_names, function(data) data$"log2_TPM+1")

# log2(TPM+1) 열을 행렬로 변환
TPM_matrix <- do.call(cbind, TPM_list)

# Quantile normalization 수행
norm_matrix <- normalize.quantiles(TPM_matrix)

# 다시 데이터프레임 형태로 변환
norm_data_list <- lapply(seq_along(data_list_gene_names), function(i) {
  data <- data_list_gene_names[[i]]
  data$"log2_TPM+1" <- norm_matrix[, i]
  return(data)
})

# 결과 리스트에 데이터프레임의 이름 할당
names(norm_data_list) <- data_names

# 이름 변경
norm_data_list_final <- norm_data_list


#-----------------------------------------------------------------------------#


### 6. Normal 대비 tumor의 log2(TPM+1) fold change 구하기 ###


# Input: Quantile normalized 데이터프레임 리스트 (norm_data_list_final)
# Output: T/N log2(TPM+1) fold change 계산된 데이터프레임 리스트 (gene_fold_change)


# 주어진 데이터에서, Normal sample은 "홀수N"으로, tumor sample은 "짝수T"로 명명되어있는 것을 이용

# 홀수 환자 선택 (최고 큰 홀수 = 5759)
odd_patients <- seq(from = 1, to = 5759, by = 2)

# 결과를 저장할 리스트 초기화
fold_change_list <- list()

# 홀수N, 짝수T 환자 샘플 이름 지정
for (patient in odd_patients) {
  normal_id <- sprintf("%dN", patient)
  tumor_id <- sprintf("%dT", patient + 1)
  
  # 홀수 normal과 짝수 tumor 선택
  normal_sample <- norm_data_list_final[[normal_id]]
  tumor_sample <- norm_data_list_final[[tumor_id]]
  
  # T/N log2 fold change 계산
  fold_change <- tumor_sample$"log2_TPM+1" - normal_sample$"log2_TPM+1"
  
  # 결과가 비어있지 않은 경우에만 리스트에 추가 (환자 ID 및 tumor ID 조합으로 저장)
  if (!all(is.na(fold_change))) {
    result_id <- sprintf("%dN/%dT", patient, patient + 1)
    fold_change_list[[result_id]] <- fold_change
  }
}

# 만들어진 T/N fold change list를 다시 gene_id, gene_name에 대입하기

# data_list_gene_names로부터 대표로 하나의 데이터프레임 추출하기
data_gn_111N <- norm_data_list_final[[1]]

# fold change 값을 포함하는 리스트 (fold_change_list)

# 결과를 저장할 리스트 초기화
gene_fold_change <- list()

# fold_change_list의 이름을 가져와서 처리
for (result_id in names(fold_change_list)) {
  # fold_change 값을 가져옴
  fold_change <- fold_change_list[[result_id]]
  
  # gene_id와 gene_name을 data_gn_111N에서 추출하여 데이터프레임에 추가
  fold_change_df <- data.frame(gene_id = data_gn_111N$gene_id,
                               gene_name = data_gn_111N$gene_name,
                               Fold_change_log2_TPM_plus_1 = fold_change)
  
  # 결과 데이터프레임을 리스트에 추가하고 이름을 result_id로 설정
  gene_fold_change[[result_id]] <- fold_change_df
}

# 결과 리스트 확인
head(gene_fold_change)


#-----------------------------------------------------------------------------#


### 7. gene_fold_change의 모든 환자 데이터프레임 병합 및 matrix화 ###


# Input: T/N log2(TPM+1) fold change 계산된 데이터프레임 리스트 (gene_fold_change)
# Output: 1열-gene_name, 각 열-환자샘플 / 각 행-gene으로 log2(TPM+1) fold change 정보가 들어있는 병합된 데이터프레임 (merged_df)
# Tool: dplyr package


library(dplyr)

# gene_fold_change 리스트에서 Fold_change_log2_TPM_plus_1 값만 추출하여 데이터프레임으로 만들기
fold_change_df_list <- lapply(gene_fold_change, function(df) {
  df_subset <- df[, "Fold_change_log2_TPM_plus_1", drop = FALSE]
  return(df_subset)
})

# 데이터프레임을 합치기 (각 데이터프레임의 이름은 열 이름으로 지정)
merged_df <- bind_cols(fold_change_df_list)
colnames(merged_df) <- names(gene_fold_change)

# gene_name 열을 추가할 위치 정하기
insert_position <- 1

# merged_df에 gene_name 열 추가
merged_df <- cbind(merged_df[, 1:insert_position - 1], gene_fold_change[[1]]$gene_name, merged_df[, insert_position:ncol(merged_df)])

# 열의 이름을 "gene_name"으로 설정
colnames(merged_df)[insert_position] <- "gene_name"


#-----------------------------------------------------------------------------#


### 8. Median absolute deviation (MAD) 기준으로 상위 10%, 20%, 30%  해당하는 gene 추출 ###


# Input: 환자/gene log2(TPM+1) fold change 정보가 있는 병합 데이터프레임 (merged_df = gene_data)
# Output: MAD값 기준 상위 10/20/30%에 해당하는 gene 정보만 담긴 matrix, gene name은 제외 (M10/M20/M30)
# Tool: CancerSubtypes package


library(CancerSubtypes)

# 이름 설정
gene_data <- merged_df

# Distribution 확인 (CancerSubtypes package)
data.checkDistribution(gene_data[, -1])


# 각 gene의 MAD 계산 (1열 gene name 제외)
gene_mad <- apply(gene_data[, -1], 1, mad)

# 상위 10%, 20%, 30%에 해당하는 MAD 값을 계산
top_10_percent <- quantile(gene_mad, probs = 0.9)
top_20_percent <- quantile(gene_mad, probs = 0.8)
top_30_percent <- quantile(gene_mad, probs = 0.7)

# 각 MAD값 기준 이상의 gene 선택
selected_genes_10 <- gene_data[gene_mad >= top_10_percent, ]
selected_genes_20 <- gene_data[gene_mad >= top_20_percent, ]
selected_genes_30 <- gene_data[gene_mad >= top_30_percent, ]

# ExecuteCNMF에 사용할 수 있는 데이터 포맷으로 변경: matrix - 각 행에 gene, 각 열에 환자
M10 <- as.matrix(selected_genes_10[, -1])
M20 <- as.matrix(selected_genes_20[, -1])
M30 <- as.matrix(selected_genes_30[, -1])

# Distribution check
data.checkDistribution(M10)
data.checkDistribution(M20)
data.checkDistribution(M30)

# Gene name 열 추가
M10_with_name <- cbind(gene_name = selected_genes_10$gene_name, M10)
M20_with_name <- cbind(gene_name = selected_genes_20$gene_name, M20)
M30_with_name <- cbind(gene_name = selected_genes_30$gene_name, M30)


#-----------------------------------------------------------------------------#


### 9. Non-negative 치환: 양수 -> 0 / 음수 -> 0  두 matrix로 나눈 후 리스트화 ###


# Input: MAD값 기준 상위 10/20/30%에 해당하는 gene 정보만 담긴 matrix (M10/M20/M30)
# Output: log2(TPM+1) fold change 값의 양수 or 음수만 담고 있는 두 개의 matrix를 포함한 리스트 (List_10/20/30)


# 각 matrix를 음수값이 없는 두 개의 matrix로 나누기 (음/양수값을 0으로 치환)

# 음수 값을 0으로 바꾸어주는 함수 지정
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

### 10. CNMF clustering 실행 ###


# Input: log2(TPM+1) fold change 값의 양수 or 음수만 담고 있는 두 개의 matrix를 포함한 리스트 (List_10/20/30)
# Output: Clustering 결과 (group, distance, consensus, cophenetic cor), silhouette width plot (sil_plot)
# Tool: CancerSubtypes, ComplexHeatmap package


library(CancerSubtypes)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# ExecuteCNMF 함수 실행
# k=2~6, nrun = 100 이상 권장
# MAD 10/20/30%, k 개수별로 클러스터링을 수행하고, 결과로부터 사용할 MAD 기준과 cluster 수 결정

# 결과를 저장할 빈 리스트 생성
result10_list <- list()
result20_list <- list()
result30_list <- list()

# k=2부터 6까지 반복 실행 (MAD 10%, 20%, 30% 이상 genes)
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

# 위 결과를 통해 MAD 20% 상위 gene들에 대해 cluster 2개로 구분하는 것으로 결정

# Group 정보 / Distance 정보 / Silhouette width / Consensus 정보 저장
group20_2 = result20_list[["result20_2"]][["group"]] #그룹핑 정보
distanceMatrix20_2=result20_list[["result20_2"]][["distanceMatrix"]] #Distance 정보
silhouette20_2=silhouette_SimilarityMatrix(group20_2, distanceMatrix20_2) #Silhouette score
Consensus20_2=result20_list[["result20_2"]][["originalResult"]]@consensus #Consensus 정보

# 간단히 Consensus matrix plotting
Heatmap(Consensus20_2, name="Consensus value",
        show_column_names = FALSE,
        show_row_names = FALSE
)
        
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


### 1. Consensus가 잘 되어있는 환자만 filtering: gene expression 데이터 추출 ###


# Input: silhouette positive patients 정보 (positive_sil_matrix), 전체 gene에 대한log2(TPM+1) fold change 정보 (f2c_analysis)
# Output: 각 환자 cluster별 분리된 gene expression 데이터 (f2c_pts_1/_2)


# Silhouette score가 양수인 환재 샘플 = 해당 subtype에 consensus가 잘 되어있는 환자만 분석에 사용

# silhouette 정보 matrix 생성
sil_matrix <- as.matrix(silhouette20_2[, 1:3])

# 환자넘버링 추출
pts <- colnames(M20)

# sil_matrix의 1열에 pts 추가
sil_matrix <- cbind(pts = pts, sil_matrix)

# silhouette score가 양수인 행을 선택
positive_sil_matrix <- sil_matrix[sil_matrix[, "sil_width"] > 0, ]

# 실루엣 양수인 환자들의 clustering 정보 extract  
clusters <- positive_sil_matrix[, "cluster"]
cluster_1_pts <- positive_sil_matrix[clusters == 1, ]
pts_1 <- cluster_1_pts[, "pts"]
cluster_2_pts <- positive_sil_matrix[clusters == 2, ]
pts_2 <- cluster_2_pts[, "pts"]


# 환자 cluster별 gene expression data 분리
f2c <- as.matrix(gene_data) #gene name 포함한 전체 gene exp 데이터
f2c_analysis <- as.matrix(gene_data[,-1]) #분석용 (gene name 없이 수치만)
f2c_analysis <- apply(f2c_analysis, 2, as.numeric) #numeric factor들로 변환
f2c_pts_1 <- f2c_analysis[, colnames(f2c_analysis) %in% pts_1]
f2c_pts_2 <- f2c_analysis[, colnames(f2c_analysis) %in% pts_2]


#-----------------------------------------------------------------------------#


### 2. Subtype간 t-test 수행하고 보정 ###


# Input: 각 환자 cluster별 분리된 gene expression 데이터 (f2c_pts_1/_2)
# Output: t_test 결과 with adjusted P value (t_test)
# Tool: stats package


library(stats)

# t_test 결과 저장할 빈 데이터프레임 생성
t_test <- data.frame(Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)

# 각 gene에 대해 t-test 수행 (16071 = 사용할 전체 gene 수)
for (i in 1:16071) {
  gene <- i
  
  # 두 그룹 간의 t-test 수행
  t_result <- t.test(f2c_pts_1[i, ], f2c_pts_2[i, ])
  
  # 결과 저장
  t_test <- rbind(t_test, data.frame(Gene = gene, P_Value = t_result$p.value))
}

# Bonferroni 보정을 적용하여 유의수준 보정
t_test$Adjusted_P_Value <- p.adjust(t_test$P_Value, method = "bonferroni")

# 결과 brief하게 확인: 보정된 t-test 결과가 훨씬 적은 DEG 산출
hist(t_test$P_Value, 100)
hist(t_test$Adjusted_P_Value, 100)


#-----------------------------------------------------------------------------#


### 3. Adjusted P value와 median 값 기준으로 각 그룹별 발현되는 genes 결정 ###


# Input: t_test 결과 with adjusted P value (t_test)
# Output: adjusted P value와 median 기준으로 필터링된 signature gene 정보 (f2c_DEG_med_combined)
# Tool: dplyr, CancerSubtypes package


library(dplyr)
library(CancerSubtypes)

# Adjusted P Value가 0.05 이하인 gene들만 추출
significant_genes <- t_test[t_test$Adjusted_P_Value <= 0.05, ]

# significant_genes에 해당되는 gene expression data만 추출
f2c_DEG <- f2c[significant_genes$Gene, ] #전체 silhouette score 양수 환자들
# 각 subtype 환자들
f2c_DEG_pts_1 <- f2c_DEG[ ,colnames(f2c) %in% pts_1] 
f2c_DEG_pts_2 <- f2c_DEG[ ,colnames(f2c) %in% pts_2] 
# numeric 치환
f2c_DEG_pts_1 <- apply(f2c_DEG_pts_1, 2, as.numeric) 
f2c_DEG_pts_2 <- apply(f2c_DEG_pts_2, 2, as.numeric)

# 각 그룹에서 gene별 median 값 추출
median_vals_group1 <- apply(f2c_DEG_pts_1, 1, median)
median_vals_group2 <- apply(f2c_DEG_pts_2, 1, median)

# 해당 median 값만 비교하기 위해 f2c_median 매트릭스 생성
f2c_median <- cbind(gene_name = f2c_DEG[ ,1], Group1_Median = median_vals_group1, Group2_Median = median_vals_group2)
# 데이터프레임 변환
f2c_median_df <- as.data.frame(f2c_median)

# Group의 median 정보를 numeric으로 변환
f2c_median_df$Group1_Median <- as.numeric(as.character(f2c_median_df$Group1_Median))
f2c_median_df$Group2_Median <- as.numeric(as.character(f2c_median_df$Group2_Median))

# 각 그룹에서 median 값이 0보다 크고, 타 그룹에선 0보다 작으며, 타 그룹보다 큰 gene들 선별
filtered_median_g1 <- f2c_median_df %>%
  filter(Group1_Median > 0 & Group2_Median < 0 & Group1_Median > Group2_Median)
filtered_median_g2 <- f2c_median_df %>%
  filter(Group2_Median > 0 & Group1_Median < 0 & Group2_Median > Group1_Median)

# gene_name 열이 공백이 아닌 행만 남기기
filtered_median_g1 <- filtered_median_g1 %>% filter(gene_name != "")
filtered_median_g2 <- filtered_median_g2 %>% filter(gene_name != "")

# filtered_median_genes에서 1열의 gene_name 추출
selected_g1 <- filtered_median_g1[, 1]
selected_g2 <- filtered_median_g2[, 1]

# silhouette 양수인 환자에서 DEG 기준을 정했으나, 다시 모든 환자에서 적용해야함.
# 전체 환자에 대한 clustering 정보 extract  
clusters <- sil_matrix[, "cluster"]
pts_all <- sil_matrix[, "pts"]

cluster_1_pts <- sil_matrix[clusters == 1, ]
pts_1 <- cluster_1_pts[, "pts"]
cluster_2_pts <- sil_matrix[clusters == 2, ]
pts_2 <- cluster_2_pts[, "pts"]

# gene_data에서 gene_name이 selected_genes와 일치하는 행 선택
gene_DEG_med_g1 <- gene_data[gene_data$gene_name %in% selected_g1, ]
gene_DEG_med_g2 <- gene_data[gene_data$gene_name %in% selected_g2, ]

# f2c_DEG_med_filtered에 Heatmap용 matrix 생성
# 이 때, 환자수가 많은 그룹부터 배열하여 consensus 맵 배열과 일치시키자
# 해당 순서가 subtype 순서로 결정
f2c_DEG_med_g1 <- as.matrix(gene_DEG_med_g1[,-1])
f2c_DEG_med_g2 <- as.matrix(gene_DEG_med_g2[,-1])
f2c_DEG_med_g1g2 <- rbind(f2c_DEG_med_g1, f2c_DEG_med_g2)
f2c_DEG_med_pts1 <- f2c_DEG_med_g1g2[,colnames(f2c_DEG_med_g1g2) %in% pts_1]
f2c_DEG_med_pts2 <- f2c_DEG_med_g1g2[,colnames(f2c_DEG_med_g1g2) %in% pts_2]
f2c_DEG_med_combined <- cbind(f2c_DEG_med_pts1, f2c_DEG_med_pts2)

# RNAseq: pts_1 = group 1, pts_2 = group 2

# 결과 확인
data.checkDistribution(f2c_DEG_med_combined)


#-----------------------------------------------------------------------------#
##### Part 3. Visualization by heatmap #####
#-----------------------------------------------------------------------------#


### 1. Heatmap annotation용 데이터 저장 ###


# Input: 필터링된 signature gene 정보 (f2c_DEG_med_combined), 그룹 정보 (sil_matrix$cluster)
# Output: Heatmap 그룹핑에 사용할 annotation 정보 (csv 파일)


# Cluster: Signature molecule의 column annotation에 사용
write.csv(colnames(f2c_DEG_med_combined), "rna_subgroup_cluster.csv", row.names = TRUE)
# 파일 추출 후에 group별 환자수만큼 1열에 RNA1, RNA2 기입
# 1행에는 Cluster 기입

# Signature genes: Signature molecule의 row annotation에 사용
write.csv(f2c_DEG_med_combined[,1], "rna_subgroup_signature.csv", row.names = TRUE)
# 파일 추출 후에 group별 signature gene수만큼 1열에 Rna1, Rna2 기입
# 1행에는 Signature 기입

# Consensus
write.csv(sil_matrix[,1:2], "rna_subgroup_consensus.csv", row.names = TRUE)
# 파일 추출 후에 group별 환자에 맞게 1열에 RNA1, RNA2 기입
# 1행에는 Consensus 기입


#-----------------------------------------------------------------------------#


### 2. Heatmap 그리기 ###


# Input: signature gene 정보 (f2c_DEG_med_combined), consensus 정보 (Consensus20_2), annotation 파일들
# Output: Signature gene heatmap (heatmap1), Consensus heatmap (heatmap2) 
# Tool: ComplexHeatmap, circlize, colorRamp2 package


library(ComplexHeatmap)
library(circlize)
library(colorRamp2)

# Signature molecules heatmap

# Color 지정
col_fun = colorRamp2(c(-2, 0, 2), c("green", "black", "red"))

# Annotation 지정 > 위에서 셋팅했던 csv 파일로부터 불러옴
# Clustering annotation (column)
anno_C = read.csv("C:/Users/admin/Desktop/LAIDD/RNAseq/Clustering/Subgroup/rna_subgroup_cluster.csv")
ha_C=HeatmapAnnotation(df = anno_C, which="col", col = list(
  Cluster = c("RNA1" = "firebrick1", "RNA2" = "darkorange3")), show_annotation_name = FALSE, show_legend = FALSE
)
# Signature annotation (row)
anno_S = read.csv("C:/Users/admin/Desktop/LAIDD/RNAseq/Clustering/Subgroup/rna_subgroup_signature.csv")
ha_S=HeatmapAnnotation(df = anno_S, which="row", col = list(
  Signature = c("Rna1" ="firebrick1", "Rna2" = "darkorange3")), show_annotation_name = FALSE, show_legend = FALSE
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
anno_A = read.csv("C:/Users/admin/Desktop/LAIDD/RNAseq/Clustering/Subgroup/rna_subgroup_consensus.csv")
ha_A=HeatmapAnnotation(df = anno_A, which="col", col = list(
  Consensus = c("RNA1" = "firebrick1", "RNA2" = "darkorange3")), show_annotation_name = FALSE, show_legend = FALSE
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
# Output: group에 속한 여부가 1 또는 0으로 표현된 indicative matrix (indi_rna.csv) 


# sil_matrix transpose
t_sil_matrix <- t(sil_matrix)

# 1행을 rownames로 설정
colnames(t_sil_matrix) <- t_sil_matrix[1,]

# 1행 제거
t_sil_matrix <- t_sil_matrix[-c(1,4), ]

# 열이름 변경 (다른 proteomics dataset들과 이름 동일하게)
colnames(t_sil_matrix) <- gsub("/", "", colnames(t_sil_matrix))
colnames(t_sil_matrix) <- sub("(\\d+)N(\\d+)T", "N\\1T\\2", colnames(t_sil_matrix))

# 1행의 데이터가 1이면 1로, 2이면 0으로 치환
t_sil_matrix[1, ] <- ifelse(t_sil_matrix[1, ] == 1, 1, 0)

# 2행의 데이터가 1이면 2행에 1을,  0이면 1을 넣음
t_sil_matrix[2, ] <- ifelse(t_sil_matrix[2, ] == 1, 1, 0)

# numeric 치환
t_sil_matrix <- apply(t_sil_matrix, 2, as.numeric)

# 저장할 이름 설정 및 행 이름 변경 (subtype 이름별)
indi_rna <- t_sil_matrix
rownames(indi_rna)[2] <- "RNA2"
rownames(indi_rna)[1] <- "RNA1"

# indicator matrix를 csv로 저장
write.csv(indi_rna, "indi_rna.csv", row.names = TRUE)









