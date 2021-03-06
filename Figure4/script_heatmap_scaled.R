##### Select the work directory (change the path to your files) 
setwd("./")
##### Load the gplots package
library(gplots)
##### Load the txt file (tab delimited)
matrix_r<-read.table("heatmap_template.txt", sep="\t", stringsAsFactors=F, header=T, row.names=1)
##### Define the colors you want to use for the heatmap
heat.colors<-colorRampPalette(c("navy","blue4","blue","skyblue","khaki1","lightgoldenrod1","goldenrod1","orange"))(100)
##### Make a heatmap  ##change parameters if needed (type help(heatmap.2) to get a full description of the options)
scale.data<-as.matrix((t(matrix_r)-apply(t(matrix_r),1,mean))/apply(t(matrix_r),1,sd))
##### Plot heatmap 
heatmap.2(as.matrix(t(scale.data)),
    dendrogram="both", scale="none",  na.color="grey",
      col=heat.colors,trace="none",labRow = rownames(t(scale.data)),key = TRUE,keysize = 1, cexCol=1,
      	 density.info="none",symkey=FALSE, margins=c(7,10), main="iMFI", 
      	 	xlab="Markers", ylab= "CLUSTERS")
