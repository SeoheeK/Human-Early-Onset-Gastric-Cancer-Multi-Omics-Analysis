#package load
library(ComplexHeatmap)
library(colorRamp2)

#input file: clinical meta data
file_path="patientmeta.csv"
anno = read.csv(file_path)

#annotation bar
ha=HeatmapAnnotation(df = anno[,-1], col = list(
  EBV = c("EBV-" = "white", "EBV+" = "orange"),
  MSI = c("MSI-H" = "green", "MSS/MSI-L" = "white"),
  Gender = c("F" = "pink", "M" = "blue"),
  Histology = c("Diffuse" = "white", "Intestinal" = "darkorchid4", "Mixed" = "darkorchid1", "Others" = "gray")), show_legend=FALSE, show_annotation_name = TRUE, annotation_name_offset = unit(2, "mm"),border=TRUE)
#save png
png("Patient_Meta_Annotationbar.png", width = 500, height = 400)
plot(ha)
dev.off()

#input file: mutation subtype data
file_path="mutationtype.csv"
anno2 = read.csv(file_path)
#annotation bar
ha2=HeatmapAnnotation(df = anno2[,-1], col = list(
  CDH1 = c("Frame shift" = "red", "Missense" = "green", "Nonsense" = "gray", "Splice site" = "purple", "In frame Indel" = "pink","0"="white"),
  TP53 = c("Frame shift" = "red", "Missense" = "green", "Nonsense" = "gray", "Splice site" = "purple", "In frame Indel" = "pink","0"="white"),
  ARID1A = c("Frame shift" = "red", "Missense" = "green", "Nonsense" = "gray", "Splice site" = "purple", "In frame Indel" = "pink","0"="white"),
  BANP = c("Frame shift" = "red", "Missense" = "green", "Nonsense" = "gray", "Splice site" = "purple", "In frame Indel" = "pink","0"="white"),
  RHOA = c("Frame shift" = "red", "Missense" = "green", "Nonsense" = "gray", "Splice site" = "purple", "In frame Indel" = "pink","0"="white"),
  MUC5B = c("Frame shift" = "red", "Missense" = "green", "Nonsense" = "gray", "Splice site" = "purple", "In frame Indel" = "pink","0"="white")
), show_annotation_name = TRUE, annotation_name_offset = unit(2, "mm"),border=TRUE, show_legend=FALSE)
#save png
png("Mutationtype_Annotation.png", width = 500, height = 400)
plot(ha2)
dev.off()
