#-------------------------------------------------------------------------------------------------
# 02 Correlation between Somatic Mutations and Phosphopeptides
#-------------------------------------------------------------------------------------------------

#load package 
library(ComplexHeatmap)
library(colorRamp2)
library(circlize)

# Set the variable name
name <- "RHOA"

# input data
file_path1=paste0("C:/Users/", name, "_heatmapinput.csv")
file_path2=paste0("C:/Users/", name, "_patient_meta.csv")

data <- read.csv(file_path1)
anno <- read.csv(file_path2)

# annotation bar
ha <- HeatmapAnnotation(df = anno, col = list(
  EBV = c("EBV-" = "white", "EBV+" = "orange"),
  MSI = c("MSI-H" = "green", "MSS/MSI-L" = "white"),
  Gender = c("F" = "pink", "M" = "blue"),
  Histology = c("Diffuse" = "white", "Intestinal" = "darkorchid4", "Mixed" = "darkorchid1", "Others" = "gray"),
  Mutation = c("Mut" = "black", "WT" = "white")
), show_annotation_name = TRUE, annotation_name_offset = unit(2, "mm"), border = TRUE)

# heatmap
heatmap1 <- Heatmap(
  as.matrix(data),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = colorRamp2(c(-1, 0, 1), c("green", "black", "red")),
  show_column_names = FALSE,
  top_annotation = ha,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = "log2foldchange", at = c(-1, 0, 1), title_position = "leftcenter-rot", legend_height = unit(4, "cm"))
)

# save figure
file_path=paste0("C:/Users/", name, ".png")
png(file_path, width = 400, height = 800)
plot(heatmap1)
dev.off()