# Title     : TODO
# Objective : TODO
# Created by: simonepuccio
# Created on: 29/01/2019
# package import
# package import
library("pheatmap")
library("RColorBrewer")
library(readr)
library(EnhancedVolcano)
# Read matrix of TPM counts from rnaseq analysis. The matrix is already parsed by cutting only cols of norm counts
matrix_heat <- read.delim("~/Documents/HumanitasProjects/SP010/tables/heatmap_tableqvalue.txt",sep="\t",header = TRUE)
# set gene_id as row.name
row.names(matrix_heat) <- matrix_heat$Id
# delete the duplicate col
matrix_heat$Id <- NULL
# logaritmic trasformation
matrix_heat2 <- matrix_heat+0.000001
Logmatrix_heat <- log2(matrix_heat2)
# set palette, ref to this for change this option https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf
col <- colorRampPalette(c("#8E4E8E","#000000","#FFD700"))(1000)

annotation_col = data.frame(
  Markers = factor(rep(c("CCR8mICOSm","CCR8pICOSp"), c(5,5)))
  #Donors = factor(rep(c("D111","D124","D140","D141","D148","D111","D124","D140","D141","D148"), c(1,1,1,1,1,1,1,1,1,1))
)
row.names(annotation_col) <- colnames(Logmatrix_heat)

ann_colors = list(
  Markers = c(CCR8mICOSm = "blue", CCR8pICOSp = "red")
)
test_label = factor(rep(c("","BATF","BATF3",""), c(300,1,1,2367)))


pheatmap(Logmatrix_heat,
         color = col,
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         cex=1,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 8, # Make fonts smaller
         cluster_cols = T,
         cluster_rows = T,
         cutree_cols = 2,
         annotation_colors = ann_colors,
         annotation_col =  annotation_col,
         show_colnames = F,
         show_rownames = F,
         clustering_method="average",
         main="")