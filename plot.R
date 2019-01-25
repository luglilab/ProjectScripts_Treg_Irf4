# source("http://bioconductor.org/biocLite.R")
# biocLite("flowCore")
# library(devtools)
# install_github("ParkerICI/premessa")
#devtools::install_github("ParkerICI/grappolo")
# set library
library("flowCore")
library("cytofCore")
library("grappolo")
############################# Premessa 
set.seed(123456)
maindir<-"/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/Input_csv/"
dir <- getwd()
files <- list.files(maindir ,pattern='.csv$', full=F)

for (i in files) {  
  setwd("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/Input_csv/")
  datamatrix<-read.csv(i)
  datamatrix <- datamatrix[,-c(1:4,10,12,15,18,30)] 
  setwd("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/Input/")
  name<-strsplit(i, ".csv")[[1]]
  #name<-unlist(name)
  
  cytofCore.write.FCS(as.matrix(datamatrix), filename=paste( "sp_", name,".fcs", sep=""), what = "numeric")     
  setwd("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/Input_csv/") } 
####################### Grappolo
library("grappolo")
col.names <- c("Comp.APC.A....CXCR3","Comp.APC.H7.A....CD73","Comp.APC.R700.A....CD25","Comp.B.710_50.A....TIGIT","Comp.BUV395.A....CD69","Comp.BUV563.A....CD45RA","Comp.BUV661.A....HLA.DR" ,
               "Comp.BUV737.A....CD95","Comp.BV421.A....PD.1","Comp.BV480.A....Ki.67","Comp.BV570.A....CD27","Comp.BV605.A....CD57","Comp.BV650.A....TIM.3","Comp.BV711.A....CCR7",
               "Comp.BV786.A....CXCR5","Comp.FITC.A....CD98","Comp.YG.586_15.A....IRF4","Comp.YG.610_20.A....Eomes","Comp.YG.670_30.A....CD71","Comp.YG.710_50.A....FoxP3","Comp.YG.780_60.A....T.bet")
# cluster for single sample
setwd("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/ClustersSingle/")
cluster_fcs_files_in_dir("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/Input/", num.cores = 2, num.clusters = 29,col.names = col.names,
                         asinh.cofactor = NULL)
# cluster merged
setwd("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/ClustersSingle/")

files.groups <- list(all.samples = list.files(pattern = "*.fcs$"))
cluster_fcs_files_groups(files.groups, num.cores = 2, num.clusters = 29, asinh.cofactor = NULL,
                         downsample.to = 5000, col.names = col.names,
                         output.dir = "/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/ClusterGrouped/") #col.names mandatory 

setwd("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/Input/")
Tumor.files.list <- list.files(pattern = ".fcs", path = "/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/tumor/", full.names = TRUE)
Tissue.files.list <- list.files(pattern = ".fcs", path = "/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/tissue/", full.names = TRUE)
Blood.files.list <- list.files(pattern = ".fcs$", path = "/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/blood/", full.names = TRUE)

files.groups <- list(
  Tumor.pooled = Tumor.files.list,
  Tissue.pooled = Tissue.files.list,
  Blood.pooled = Blood.files.list
)
cluster_fcs_files_groups(files.groups, num.cores = 2, col.names = col.names, 
                         num.clusters = 29, asinh.cofactor = NULL, downsample.to = 50000, output.dir = "/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/clustered_by_tissue/")
########################### Vite 
library(vite)
setwd("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/clustered_by_tissue/")
input.files <- list.files(path = "/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/clustered_by_tissue", pattern = "*.clustered.txt$", full.names = TRUE)
col.names <- c("Comp.APC.A....CXCR3","Comp.APC.H7.A....CD73","Comp.APC.R700.A....CD25","Comp.B.710_50.A....TIGIT","Comp.BUV395.A....CD69","Comp.BUV563.A....CD45RA","Comp.BUV661.A....HLA.DR" ,
               "Comp.BUV737.A....CD95","Comp.BV421.A....PD.1","Comp.BV480.A....Ki.67","Comp.BV570.A....CD27","Comp.BV605.A....CD57","Comp.BV650.A....TIM.3","Comp.BV711.A....CCR7",
               "Comp.BV786.A....CXCR5","Comp.FITC.A....CD98","Comp.YG.586_15.A....IRF4","Comp.YG.610_20.A....Eomes","Comp.YG.670_30.A....CD71","Comp.YG.710_50.A....FoxP3","Comp.YG.780_60.A....T.bet")

setwd("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/Input/")
filenames <- list.files(dir ,pattern='.fcs$', full=F)
combined_data_transformed_tumor <- cytof_exprsMerge(fcsFiles = Tumor.files.list, comp=FALSE,                                             
                                              transformMethod = "none",
                                              mergeMethod = "all")

cytofCore.write.FCS(combined_data_transformed_tumor, filename="/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/landmarks/combined_data_transformed_tumor.fcs", what = "numeric")

combined_data_transformed_tissue <- cytof_exprsMerge(fcsFiles = Tissue.files.list, comp=FALSE,                                             
                                                    transformMethod = "none",
                                                    mergeMethod = "all")

cytofCore.write.FCS(combined_data_transformed_tissue, filename="/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/landmarks/combined_data_transformed_tissue.fcs", what = "numeric")

combined_data_transformed_blood <- cytof_exprsMerge(fcsFiles = Blood.files.list, comp=FALSE,                                             
                                                    transformMethod = "none",
                                                    mergeMethod = "all")

cytofCore.write.FCS(combined_data_transformed_blood, filename="/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/landmarks/combined_data_transformed_blood.fcs", what = "numeric")
#
landmarks.data <- vite::load_landmarks_from_dir("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/landmarks/", asinh.cofactor = 5, transform.data = T)


#qui

vite::run_scaffold_analysis(input.files, ref.file = input.files[1], 
                            col.names=col.names,landmarks.data,process.clusters.data = F )


G <- vite::get_unsupervised_graph_from_files(input.files, col.names = col.names, filtering.threshold = 15)
setwd("/Users/simonepuccio/Documents/HumanitasProjects/Scaffold_pipeline/")
vite::write_graph(G, "first_unsupervised.graphml")
