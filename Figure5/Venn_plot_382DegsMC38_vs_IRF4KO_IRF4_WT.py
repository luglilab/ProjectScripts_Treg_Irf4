import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn3_circles
import pandas as pd
# set color
# IRF4 = "#ff6666"
# BATF = "#6666ff"
# chipcolor = "#b2b2b2"
# input file
#### Prepare Input list for Venn Plot
Genelist382 = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/382CommonGenes.txt"
# File with all expressed genes
ExpressedMC38vsSpleen = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/EdgeR_Pnas_GSE116347/tables/MC38vsSpleen.complete.anno.txt"
# create dataframe
dfGenelist382 = pd.read_csv(Genelist382, sep="\t",header=None)
dfExpressedMC38vsSpleen = pd.read_csv(ExpressedMC38vsSpleen,sep="\t",header=0)
# Merge
dfMerge = pd.merge(dfGenelist382,dfExpressedMC38vsSpleen,left_on=0,right_on="Id")
# remove duplicates
dfMerge.drop_duplicates(subset =0, inplace = True)
# drop columns
dfMerge = dfMerge.drop([0,'Unnamed: 0'], axis=1)
genenames_MC38 = pd.Series(dfMerge['Id'], index=dfMerge.index)
###### Prepare input genename IRF4Ko and IRFWt
Irf4ko_wt_pval005 = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/IRF4ko_IRF4wt_mouse.Pval.annotated.txt"
Batf4ko_wt_pval005 = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/BATFko_BATFwt_mouse.Pval.annotated.txt"
# create dataframe
dfIrf4ko_wt_pval005 = pd.read_csv(Irf4ko_wt_pval005,sep="\t",header=0)
dfBatf4ko_wt_pval005 = pd.read_csv(Batf4ko_wt_pval005,sep="\t",header=0)
# drop duplicates
dfIrf4ko_wt_pval005.drop_duplicates(subset = "symbol", inplace = True)
dfBatf4ko_wt_pval005.drop_duplicates(subset = "symbol", inplace = True)
# create two duplicates for logFC<0 and logFC>=0
dfIrf4ko = dfIrf4ko_wt_pval005[(dfIrf4ko_wt_pval005['logFC'] >= 0)]
dfIrf4wt = dfIrf4ko_wt_pval005[(dfIrf4ko_wt_pval005['logFC'] < 0)]
#
dfBatfko = dfBatf4ko_wt_pval005[(dfBatf4ko_wt_pval005['logFC'] >= 0)]
dfBatfwt = dfBatf4ko_wt_pval005[(dfBatf4ko_wt_pval005['logFC'] < 0)]
# convert df to series Irf4
genenames_Irf4ko = pd.Series(dfIrf4ko['symbol'], index=dfIrf4ko.index)
genenames_Irf4wt = pd.Series(dfIrf4wt['symbol'], index=dfIrf4wt.index)
# convert df to series Batf
genenames_Batfko = pd.Series(dfBatfko['symbol'], index=dfBatfko.index)
genenames_Batfwt = pd.Series(dfBatfwt['symbol'], index=dfBatfwt.index)
###