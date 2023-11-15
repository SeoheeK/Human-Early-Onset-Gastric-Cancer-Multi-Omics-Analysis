#-------------------------------------------------------------------------------------------------
# 04 Correlation between mRNA and Protein Abundance Changes
#-------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------#
### 1. Package Install ###
#----------------------------------------------------------------------------#

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()
a
BiocManager::install(c("CancerSubtypes","GenomicFeatures", "AnnotationDbi"))
a
BiocManager::available()

library(CancerSubtypes)
library(GenomicFeatures)
library(AnnotationDbi)
library(readr)
library(BiocManager)
library(CancerSubtypes)
library(sigclust)

#----------------------------------------------------------------------------#
### 2. mRNA data setting ###
#----------------------------------------------------------------------------#

## 2-1. mRNA 데이터 불러오기
# make dataset tables of 160 patientss RNAseq data at once
data_dir <- "C:/R/Data/231006 RNAseq_RSEM"   # data가 있는 폴더 directory 입력
data_files <- list.files(data_dir, pattern = "results")

# 빈 리스트 생성
data_list <- list()

# 데이터 파일을 읽어와서 리스트의 항목으로 추가
for (file in data_files) {
  file_path <- file.path(data_dir, file)
  data <- read_delim(file_path, delim = "\t")  # tab 구분자 데이터 파일인 경우
  
  # 파일 이름에서 확장자 제거
  data_name <- gsub("_rsem_genes_original_results", "", file)
  
  # 리스트의 항목으로 데이터프레임 추가
  data_list[[data_name]] <- data
}
#----------------------------------------------------------------------------#
## 2-2. 30% 이상의 환자 데이터에서 TPM 1 이상인 gene만 남기고 제거
# dplyr 패키지 로드
library(dplyr)

# 모든 데이터프레임을 하나로 합치기
combined_data <- bind_rows(data_list, .id = "Dataset")

# gene별로 30% 이상 환자에서 TPM이 1 이상인 행만 필터링
filtered_data <- combined_data %>%
  group_by(gene_id) %>%
  filter(mean(TPM >= 1, na.rm = TRUE) >= 0.3)

# 필터링된 데이터를 다시 각각의 데이터프레임으로 나누기
cleaned_data_list <- split(filtered_data, filtered_data$Dataset)
data_list <- cleaned_data_list

#----------------------------------------------------------------------------#
## 2-3. Calculate Log2(TPM+1)
# 필요한 라이브러리 불러오기 (DESeq2 사용)
BiocManager::install("DESeq2")

library(DESeq2)

# log2(TPM+1) 변환 함수 정의
log2_transform <- function(data) {
  return(log2(data + 1))
}


# 모든 TPM+1 데이터에 log2 변환 적용
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

#----------------------------------------------------------------------------#
## 2-4. gene_id to gene_name (Ensembl)
#필요 패키지
install.packages("httr")
install.packages("jsonlite")
library(httr)
library(jsonlite)

# 111N 데이터 추출
data_111N <- data_list[[1]]  # 111N 데이터 추출 (예시로 첫 번째 데이터를 사용)

# gene_id 열 추출
gene_ids <- data_111N$gene_id

# gene_ids를 문자열 벡터로 변환
gene_ids<- as.character(gene_ids)
head(gene_ids)

# BiomaRt 설치
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)

# Ensembl 데이터베이스에 연결하기
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# gene_id를 gene name으로 변환
gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", 
                    values = gene_ids, 
                    mart = ensembl)

# gene_names 데이터프레임에 NA로 표시할 열 추가
gene_names$external_gene_name[is.na(gene_names$external_gene_name)] <- "NA"

# 데이터프레임에 gene_name 열 추가
data_list_gene_names <- lapply(data_list, function(data) {
  data <- data.frame(data[, 1], gene_name = gene_names$external_gene_name[match(data$gene_id, gene_names$ensembl_gene_id)], data[, 2])
  colnames(data) <- c("gene_id", "gene_name", "log2_TPM+1")
  return(data)
})
#----------------------------------------------------------------------------#
## 2-5. Quantile normalization of Log2(TPM+1)
# Quantile normalization은 전체 데이터를 모두 일괄로 처리해야함

# preprocessCore package 설치
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
#----------------------------------------------------------------------------#
## 2-6. Option: Median centering of quantile-normalized data
# Median centering은 각각의 데이터셋별로 수행해야함

# Median centering 수행
for (i in 1:length(norm_data_list)) {
  df <- norm_data_list[[i]]
  df$"log2_TPM+1" <- df$"log2_TPM+1" - median(df$"log2_TPM+1", na.rm = TRUE)
  norm_data_list[[i]] <- df
}

# 결과 확인
head(norm_data_list[[1]])
#----------------------------------------------------------------------------#
## 2-7. 다시 한 번 quantile normalization 수행
# 모든 환자 데이터프레임에서 log2(TPM+1) 열 추출
TPM_list2 <- lapply(norm_data_list, function(data) data$"log2_TPM+1")

# log2(TPM+1) 열을 행렬로 변환
TPM_matrix2 <- do.call(cbind, TPM_list2)

# Quantile normalization 수행
norm_matrix2 <- normalize.quantiles(TPM_matrix2)

# 다시 데이터프레임 형태로 변환
norm_data_list_final <- lapply(seq_along(norm_data_list), function(i) {
  data <- norm_data_list[[i]]
  data$"log2_TPM+1" <- norm_matrix2[, i]
  return(data)
})

# 결과 리스트에 데이터프레임의 이름 할당
names(norm_data_list_final) <- data_names
#----------------------------------------------------------------------------#
## 2-8. Log2(TPM+1) fold change 구하기
# 홀수 환자 선택
odd_patients <- seq(from = 1, to = 5759, by = 2)

# 결과를 저장할 리스트 초기화
fold_change_list <- list()


# 홀수N, 짝수T 환자 샘플 이름 지정
for (patient in odd_patients) {
  normal_id <- sprintf("%dN", patient)
  tumor_id <- sprintf("%dT", patient + 1)
  
  # 홀수 환자와 짝수 환자 선택
  normal_sample <- norm_data_list_final[[normal_id]]
  tumor_sample <- norm_data_list_final[[tumor_id]]
  
  # log2 fold change 계산
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
#----------------------------------------------------------------------------#
## 2-9. gene_fold_change의 모든 환자 데이터프레임 병합 및 matrix화
# dplyr 장착
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
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
### 3. protein data setting ##
#----------------------------------------------------------------------------#

# protein 데이터 불러오기
protein_data <- read.csv("protein_expression_change.txt",header=TRUE,sep="\t")
View(protein_data)
str(protein_data)

# protein 데이터에서 Symbol 결측치 제거하기
protein_data1 <- subset(protein_data, !is.na(Symbol) & Symbol != "#N/A")
View(protein_data1)

# major peptide만 남기기
protein_data2 <- protein_data1[!duplicated(protein_data1$Symbol), ]
View(protein_data2)

# 30% 이상의 환자에서 protein abundance 측정 가능한 protein 데이터 선택
selected_protein_data <- protein_data2[protein_data2$X..of.value >= 80/3,]
View(selected_protein_data)
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
### 4. mRNA-protein mapping by gene symbol
#----------------------------------------------------------------------------#

# protein 데이터 gene symbol 추출
protein_genes <- selected_protein_data[,5]
str(protein_genes) #4675개

# mRNA 데이터 gene symbol 추출
mRNA_genes <- merged_df[,1]
str(mRNA_genes)    #16071개

# protein 데이터 gene symbol 중복 제거
set_A <- unique(protein_genes)

# mRNA 데이터 gene symbol 중복 제거
set_B <- unique(mRNA_genes)

# protein & mRNA gene의 교집합 구하기
intersection <- intersect(set_A,set_B)

# 교집합 gene의 개수 구하기
intersection_count <- length(intersection)
intersection_count  #4286개

# protein gene set 설정
final_protein_data <- selected_protein_data[selected_protein_data$Symbol %in% intersection,]
View(final_protein_data)

# mRNA gene set 설정
final_mRNA_data <- merged_df[merged_df$gene_name %in% intersection,]
View(final_mRNA_data)

# protein data에서 Symbol과 sample별 fold change 값만 남기고 제거
final_protein_data2 <- final_protein_data[,-c(1,2,3,4,6)]
View(final_protein_data2)

# protein data에서 gene symbol을 기준으로 정렬
sorted_protein <- final_protein_data2[order(final_protein_data2$Symbol),]
View(sorted_protein)

# protein data에서 환자 샘플 번호 오름차순으로 정렬 (mRNA 데이터와 환자 샘플 순서 일치시키기 위해)
sorted_protein <- sorted_protein %>% 
  select(order(as.numeric(sub("^N(\\d+).*","\\1",names(sorted_protein))))) %>% 
  select(81,1:80)

# mRNA data에서 gene symbol을 기준으로 행 정렬
sorted_mRNA <- final_mRNA_data[order(final_mRNA_data$gene_name),]
View(sorted_mRNA)

#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
### 5. Spearman Correlation Analysis
#----------------------------------------------------------------------------#
## 5-1. 패키지 설치 및 불러오기
library(ggplot2)
install.packages("reshape2")
library(reshape2)
install.packages("ggdist")
library(ggdist)
#----------------------------------------------------------------------------#
## 5-2. Spearman의 상관 계수 및 p-value 계산 (결측치 제거)

# correlation matrix 초기화 (4286 x 80)
cor_matrix <- matrix(NA, nrow = 4286, ncol = 80)

# p-value matrix 초기화 (4286 x 80)
p_value_matrix <- matrix(NA, nrow = 4286, ncol = 80)

# protein 데이터 fold change 값만 저장하기 (gene symbol 제외)
my_protein <- sorted_protein[,-1]
rownames(my_protein) <- NULL

# mRNA 데이터 fold change 값만 저장하기 (gene symbol 제외)
my_mRNA <- sorted_mRNA[,-1]
rownames(my_mRNA) <- NULL

# Spearman 상관 분석
for (i in 1:4286) {
  for (j in 1:80) {
    cor_test_result <- cor.test(as.numeric(my_mRNA[i, ]), as.numeric(my_protein[i, ]), method = "spearman", use = "complete.obs")
    cor_matrix[i, j] <- cor_test_result$estimate
    p_value_matrix[i, j] <- cor_test_result$p.value
  }
}

View(cor_matrix)
View(p_value_matrix)

# p-value 매트릭스 데이트프레임으로 변환
p_value_df <- as.data.frame(p_value_matrix)

# p-value 데이터프레임 행이름/열이름 설정 
rownames(p_value_df) <- rownames(my_mRNA)
colnames(p_value_df) <- colnames(my_mRNA)
View(p_value_df)

# 상관 계수 데이터 프레임 생성
correlation_df <- data.frame(cor_matrix=cor_matrix)
colnames(correlation_df) <- colnames(my_mRNA)
View(correlation_df)
#----------------------------------------------------------------------------#
## 5-3. 상관계수 평균값 구하기

# 상관계수 데이터프레임 -> 매트릭스로 전환
spearman_corr_matrix <- as.matrix(correlation_df)

# 상관계수 평균값 계산  
overall_mean <- mean(spearman_corr_matrix)
overall_mean  # 0.2926156 (논문: 0.28)
#----------------------------------------------------------------------------#
## 5-4. 양의 상관 계수 % 확인

# 상관계수 데이터프레임 -> 벡터로 전환
spearman_corr_vector <- as.vector(correlation_df)
spearman_corr_vector

# 전체 상관계수 개수 확인
length(rownames(correlation_df))  #4286개

# 양의 상관계수 개수 계산
positive_correlations <- colSums(spearman_corr > 0, na.rm = TRUE)
positive_correlations   # 3968개

# 양의 상관계수 % 계산
percentage_positive_corr <- (3968 / 4286) * 100
percentage_positive_corr   # 92.58% (논문 : 91.4%)
#----------------------------------------------------------------------------#
## 5-5. FDR이 0.01 미만인 (유의미한) positive correlation 계산

# 유전자, 샘플 개수 변수로 저장
total_genes <- 4286
total_samples <- 80

# p-value를 1차원 벡터로 변환하고 NA 값 제외
p_values <- unlist(p_value_df)[!is.na(unlist(p_value_df))]

# Benjamini-Hochberg FDR 계산
bh_fdr <- matrix(p.adjust(p_values, method = "BH"), nrow = total_genes, ncol = total_samples, byrow = TRUE)
View(bh_fdr)
str(bh_fdr)
positive_FDR <- bh_fdr[correlation_df > 0]
range(positive_FDR)

# FDR < 0.01 인 것들의 개수 세기
count1 <- sum(positive_FDR != 0 & positive_FDR < 0.01, na.rm = TRUE) #개수 : 121245개

# FDR < 0.01 인 것들의 비율 계산
count1/(total_genes*total_samples)*100 #비율 : 35.36077%

# FDR < 0.01 인 비율 변수로 저장
percentage_FDR <- 35.36077

#----------------------------------------------------------------------------#
## 5-6. Figure 2-A 그래프 그리기 : Distributions of Spearman's correlation coefficients

# correlation_df의 값들 unlist 하여 벡터로 저장
correlation_values <- unlist(correlation_df)

# 그래프 그리기
ggplot(data = data.frame(correlation_values = correlation_values), # correlation_values 데이터를 데이터프레임으로 만들어 사용
                        aes(x = correlation_values)) + # x축은 correlation_values
  geom_histogram(aes(y = ..density..), binwidth = 0.05, fill = "blue", boundary=0,alpha = 0.6) + # y축은 밀도, 파란색 히스토그램 생성
  geom_density(color = "red", alpha = 0.6, linewidth =1) + # 밀도 그래프는 빨간색으로 생성
  labs(title = "Spearman Correlation Analysis of mRNA and protein", # 그래프 제목 설정
       x = "Spearman Correlation Coefficient", # 그래프 x축 이름 설정
       y = "Probability Density") + # 그래프 y축 이름 설정
  scale_x_continuous(breaks = c(-0.6, -0.3, 0.0, 0.3, 0.6, 0.9), expand = expansion(add = c(0, 0.02))) + # x축 간격 설정
  geom_vline(xintercept = overall_mean, linetype = "dashed", color = "green", linewidth = 1) + # 상관계수 평균값 선 추가
  annotate("text", x = 0.15, y = 2.0, label = paste("Mean =", round(overall_mean, 3)), color = "green", vjust = -0.5) +
  annotate("text", x = -0.28, y = 2.0, label = paste("Positive Correlation :", round(percentage_positive_corr, 3),"%"), vjust = -0.5) +
  annotate("text", x = -0.2, y = 1.9, label = paste("Significant Positive Correlation :", round(percentage_FDR, 3),"%"), vjust = -0.5) +
  theme_minimal()
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
### 6. KEGG Pathway 분석 (DAVID)
#----------------------------------------------------------------------------#

## 6-1. Correlation efficient가 가장 높은 500개, 가장 낮은 500개 유전자 추출

# 상관계수 범위 확인하기
range(cor_values) # -0.4952532~0.8208819

# gene ID(ensembl), symbol 포함한 데이터 저장
sort_a <- normal_sample[normal_sample$gene_name %in% final_proteome_data$Symbol,] %>% arrange(gene_name)

# gene Ensembl ID 추출
id_of_genes <- sort_a$gene_id
id_of_genes

# gene symbol 추출
name_of_genes <- sort_a$gene_name
name_of_genes

# gene Ensembl ID, symbol, cor_values 데이터프레임 만들기
mydata <- data.frame(id_of_genes, name_of_genes,cor_values)
View(mydata)

# cor_values가 가장 큰 top 500 gene 추출
top_500_indices <- order(cor_values, decreasing = TRUE)[1:500]
highest_correlation_genes <- sort_a$gene_id[top_500_indices]

# cor_values가 가장 작은 bottom 500 gene 추출
bottom_500_indices <- order(cor_values, decreasing = FALSE)[1:500]
lowest_correlation_genes <- sort_a$gene_id[bottom_500_indices]

# 엑셀 파일로 저장하기
install.packages("openxlsx")
library(openxlsx)
highest_correlation_genes
data1 <- data.frame(Gene = highest_correlation_genes)
View(data1)
data2 <- data.frame(Gene = lowest_correlation_genes)
write.xlsx(data1, file="highest_correlation_genes.xlsx")
write.xlsx(data2, file="lowest_correlation_genes.xlsx")

#----------------------------------------------------------------------------#
## 6-2.KEGG pathway 분석 (DAVID)

# DAVID를 활용해 top 500, bottom 500에 대한 pathway 3개씩 검색 -> 엑셀 파일 저장


#----------------------------------------------------------------------------#
## 6-3. KEGG pathway 데이터 불러오기

# 첫번째 top pathway 불러오기
top_metabolic_pathway <- read_excel("6_pathways.xlsx",sheet=1) #top1
# 두번째 top pathway 불러오기
top_PPAR_signaling <- read_excel("6_pathways.xlsx",sheet=2)  #top2
# 세번째 top pathway 불러오기
top_Arginine_proline_metabolism <- read_excel("6_pathways.xlsx",sheet=3) #top3
# 첫번째 bottom pathway 불러오기
bottom_Ubiquitin_mediated_proteolysis <- read_excel("6_pathways.xlsx",sheet=4) #bottom1
# 두번째 bottom pathway 불러오기
bottom_Spliceosome <- read_excel("6_pathways.xlsx",sheet=5)  #bottom2
# 세번째 bottom pathway 불러오기
bottom_Ribosome <- read_excel("6_pathways.xlsx",sheet=6) #bottom3

###
sort_e <- sort_d
rownames(sort_e) <- 1:1000
View(sort_e)

sort_s <- sort_bb
View(sort_s)
rownames(sort_s) <- 1:4286
###

# top1의 gene ID 저장
top1 <- top_metabolic_pathway$ID
# top2의 gene ID 저장
top2 <- top_PPAR_signaling$ID
# top3의 gene ID 저장
top3 <- top_Arginine_proline_metabolism$ID
# bottom1의 gene ID 저장
bottom1 <- bottom_Ubiquitin_mediated_proteolysis$ID
# bottom2의 gene ID 저장
bottom2 <- bottom_Spliceosome$ID
# bottom3의 gene ID 저장
bottom3 <- bottom_Ribosome$ID

# sort_b의 상관계수 열을 기준으로 내림차순 정렬
sort_s <- sort_b[order(sort_b$V3, decreasing = TRUE),]
View(sort_s)

# top1 pathway에 해당하는 유전자들 index 저장
top1_index <-rownames(subset(sort_s, sort_s$gene_id %in% top1))
# top2 pathway에 해당하는 유전자들 index 저장
top2_index <-rownames(subset(sort_s, sort_s$gene_id %in% top2))
# top3 pathway에 해당하는 유전자들 index 저장
top3_index <-rownames(subset(sort_s, sort_s$gene_id %in% top3))
# bottom1 pathway에 해당하는 유전자들 index 저장
bottom1_index <-rownames(subset(sort_s, sort_s$gene_id %in% bottom1))
# bottom2 pathway에 해당하는 유전자들 index 저장
bottom2_index <-rownames(subset(sort_s, sort_s$gene_id %in% bottom2))
# bottom3 pathway에 해당하는 유전자들 index 저장
bottom3_index <-rownames(subset(sort_s, sort_s$gene_id %in% bottom3))

top1_index #617
top2_index #27
top3_index #17
bottom1_index #59
bottom2_index #77
bottom3_index #81

# DAVID에서 선택한 KEGG pathway 정보에 대해 데이터 프레임 생성성
pathways <- data.frame(
  pathway = c("top_metabolic_pathway", "top_PPAR_signaling", "top_Arginine_proline_metabolism", "bottom_Ubiquitin_mediated_proteolysis", "bottom_Spliceosome", "bottom_Ribosome"),
  gene_count = c(617,27,17,59,77,81)  # 각 pathway에 속하는 유전자 수
)
#----------------------------------------------------------------------------#
## 6-4. Figure 2-B 그래프 그리기

# 1) over_plot
plot(1:nrow(sort_s), sort_s$coefficient, type = "n",
     ylab = "Spearman's Correlation Coefficient",
     col = "black",
     xlim = c(0, 4286),
     xaxt = "n",
     xlab ="")

# y=0 이상인 부분 파란색으로 채우기
polygon(c(1, which(sort_s$coefficient >= 0), nrow(sort_s)),
        c(0, sort_s$coefficient[sort_s$coefficient >= 0], 0),
        col = "blue")

# y=0 이하인 부분 주황색으로 채우기
polygon(c(1, which(sort_s$coefficient < 0), nrow(sort_s)),
        c(0, sort_s$coefficient[sort_s$coefficient < 0], 0),
        col = "orange")

abline(h = 0, col = "black", lwd=2)

###############################################################################
# 2)under_plot

# 그래프 생성
plot(1:4286, rep(0, 4286), type = "n", 
     xlim = c(1, 4286), ylim = c(0, length(pathways) + 1), frame.plot = FALSE, xaxt = "n", yaxt = "n", xlab="", ylab="")

# 각 pathway에 대해 바코드 그리기
for (i in 1:length(pathways)) {
  barcode_height <- i  # 각 pathway의 높이
  x_values <- pathways[[i]]
  
  # x값에 따라 색상 결정
  color <- ifelse(x_values <= 3968, colors[1], colors[2])
  
  # 바코드 개수를 각 pathway의 값들의 개수에 맞게 표현
  rect(x_values[1:length(x_values)], barcode_height - 0.4,
       x_values[1:length(x_values)] + 2, barcode_height + 0.4,
       col = color, border = NA, lwd = 10)  # lwd 값으로 선의 굵기 설정
}
###############################################################################