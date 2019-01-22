library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
# set working directory
setwd("/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure4/Qvalue005/")
# set file name
irf4 = "Irf4_tumor_intersect.txt"
batf = "Batf_tumor_intersect.txt"
common = "Common.txt"
# read tables
Irf4_genesup <- read.table("~/Documents/HumanitasProjects/SP010/Figure4/Qvalue005/Irf4FCpos.txt", quote="\"", comment.char="")
BATF_genesup <- read.table("~/Documents/HumanitasProjects/SP010/Figure4/Qvalue005/BatfFCpos.txt", quote="\"", comment.char="")
Common_genesup <- read.table("~/Documents/HumanitasProjects/SP010/Figure4/Qvalue005/CommonFCpos.txt", quote="\"", comment.char="")
Irf4_genesdown <- read.table("~/Documents/HumanitasProjects/SP010/Figure4/Qvalue005/Irf4FCneg.txt", quote="\"", comment.char="")
BATF_genesdown <- read.table("~/Documents/HumanitasProjects/SP010/Figure4/Qvalue005/BatfFCneg.txt", quote="\"", comment.char="")
Common_genesdown <- read.table("~/Documents/HumanitasProjects/SP010/Figure4/Qvalue005/CommonFCneg.txt", quote="\"", comment.char="")
#
# row.names(Irf4_genes) <- Irf4_genes$V1
Irf4_genes_entrez_up = bitr(Irf4_genesup$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
BATF_genes_entrez_up = bitr(BATF_genesup$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
Common_genes_entrez_up = bitr(Common_genesup$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
Irf4_genes_entrez_down = bitr(Irf4_genesdown$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
BATF_genes_entrez_down = bitr(BATF_genesdown$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
Common_genes_entrez_down = bitr(Common_genesdown$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# create list of list
gcSample_up = list(Irf4_genes = Irf4_genes_entrez_up$ENTREZID,
                BATF_genes = BATF_genes_entrez_up$ENTREZID,
                Common_genes = Common_genes_entrez_up$ENTREZID)
gcSample_down = list(Irf4_genes = Irf4_genes_entrez_down$ENTREZID,
                   BATF_genes = BATF_genes_entrez_down$ENTREZID,
                   Common_genes = Common_genes_entrez_down$ENTREZID)
#
CompClusEnrichedPath_up <- compareCluster(geneCluster = gcSample_up, fun = "enrichPathway", pvalueCutoff=0.05)
dotplot(CompClusEnrichedPath_up, showCategory = 20)
#
CompClusEnrichedPath_down <- compareCluster(geneCluster = gcSample_down, fun = "enrichPathway", pvalueCutoff=0.05)
dotplot(CompClusEnrichedPath_down, showCategory = 30)

write.table(CompClusEnrichedPath_up,"~/Documents/HumanitasProjects/SP010/Figure4/Qvalue005/path_enrichment_up.txt",sep="\t",quote=F, row.names = F)
write.table(CompClusEnrichedPath_down,"~/Documents/HumanitasProjects/SP010/Figure4/Qvalue005/path_enrichment_down.txt",sep="\t",quote=F, row.names = F)