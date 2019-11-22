library(DOSE)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)

ImmunoSig_filtered <- read.delim("~/Documents/HumanitasProjects/SP010/ImmunoSig_filtered.txt")

p <- ggplot(ImmunoSig_filtered, aes(x = GeneRatio, y = fct_reorder(Description2, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.20), low="darkred",high = "orange") +
  ylab(NULL) +
  ggtitle("IMMUNO GSEA")

p + facet_grid(.~Cluster) 

ggsave("~/Documents/HumanitasProjects/SP010/merged_GSE2.eps",device="eps", scale = 2)
