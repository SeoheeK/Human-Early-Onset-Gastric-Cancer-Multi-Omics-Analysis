#package load
library(ComplexHeatmap)
library(colorRamp2)

#input data: consensuspathdb output
data = read.csv("-log10pvalue_3genes_pathway.csv")
#fill na
data[is.na(data)] <- 0
#set rownames
rownames(df) <- df$pathway.name
#pathway heatmap
heatmap1 = Heatmap(as.matrix(df[,-1]), rect_gp = gpar(col = "black"),border_gp = gpar(col = "black"),col=colorRamp2(c(0,3,6),c("white","yellow","brown")),  show_row_names=TRUE,show_column_names = TRUE,column_names_centered = TRUE,column_names_rot = 45,column_title_gp = gpar(fontsize = 12, fontface = 'bold'), cluster_rows=FALSE,cluster_columns=FALSE, show_heatmap_legend = FALSE,column_names_side = c("top"),width = ncol(df)*unit(5, "mm"),height = nrow(df)*unit(5, "mm"))
legend1 = Legend(col_fun = colorRamp2(c(0, 3, 6), c("white", "yellow", "brown")), title = "-log10(qvalue)", title_position = "topcenter",legend_width = unit(4, "cm"),direction = "horizontal", at = c(0, 3, 6))
#save png
png("C:/Users/smj97/Desktop/pathway_heatmap3.png", width = 800, height = 1000)
plot(heatmap1, annotation_legend_list= legend1,annotation_legend_side="bottom")
dev.off()
