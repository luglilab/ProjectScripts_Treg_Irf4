#####Select the working directory (change the path to your files) 
setwd("~/Desktop/Template")

#####Load the reshape2 and ggplot2 package
library(reshape2)
library(ggplot2)

#####Load the txt file (tab delimited)
balloon<-read.table("Template_table_balloonplot.txt", sep="\t", stringsAsFactors=F, header=T)

#####We need to original order of the clusters on the y-axis, so it matches the order of the heatmap; make sure cell A1 is called Cluster
balloon$Cluster <- factor(balloon$Cluster, levels = rev(unique(balloon$Cluster)))
balloon_melted<-melt(balloon,sort=F)

#####Make a balloonplot  --> change max_size if the balloons are overlapping
p <- ggplot(balloon_melted, aes(x = variable, y = Cluster))

p + geom_point(aes(size=value), colour="grey34") + theme(panel.background=element_blank()) +
  scale_size_area(max_size=12)

#####To save the balloonplot as pdf
pdf("Balloonplot.pdf")

p <- ggplot(balloon_melted, aes(x = variable, y = Cluster))

p + geom_point(aes(size=value), colour="grey34") + theme(panel.background=element_blank()) +
  scale_size_area(max_size=12)
  
dev.off()

###To save the balloonplot as eps for Illustrator
ggsave(file="balloonplot.eps")