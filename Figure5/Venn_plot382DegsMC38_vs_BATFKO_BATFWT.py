import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn3_circles
import pandas as pd
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
# create dataframe
dfIrf4ko_wt_pval005 = pd.read_csv(Irf4ko_wt_pval005,sep="\t",header=0)
# drop duplicates
dfIrf4ko_wt_pval005.drop_duplicates(subset = "symbol", inplace = True)
# create two duplicates for logFC<0 and logFC>=0
dfIrf4ko = dfIrf4ko_wt_pval005[(dfIrf4ko_wt_pval005['logFC'] >= 0)]
dfIrf4wt = dfIrf4ko_wt_pval005[(dfIrf4ko_wt_pval005['logFC'] < 0)]
# convert df to series Irf4
genenames_Irf4ko = pd.Series(dfIrf4ko['symbol'], index=dfIrf4ko.index)
genenames_Irf4wt = pd.Series(dfIrf4wt['symbol'], index=dfIrf4wt.index)
###
c =venn3([set(genenames_Irf4ko),set(genenames_Irf4wt),set(genenames_MC38)], ('Irf4 KO Degs', 'Irf4 WT Degs', 'TITR MC38'))
c.get_patch_by_id('100').set_color("#FFFFFF") # irf4 ko peaks
c.get_patch_by_id('010').set_color("#FFFFFF") # irf4 wt peaks
c.get_patch_by_id('001').set_color("#FFFFFF") # degs low
c.get_patch_by_id('101').set_color("orange") # intersection irf4 degs
c.get_patch_by_id('011').set_color("limegreen") # intersection chip batf degs
v=venn3_circles(subsets = (275, 370, 0, 326, 2,54,0),  linewidth=1, color="black")
#
Irf4ko_MC38_intersect = set(genenames_Irf4ko).intersection(set(genenames_MC38))
Irf4wt_MC38_intersect = set(genenames_Irf4wt).intersection(set(genenames_MC38))
#
plt.savefig("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/VennPlotIRF4KO_IRF4WT_382DEGS.eps",format='eps')

Irf4ko_MC38_intersect_list = list(Irf4ko_MC38_intersect)
Irf4wt_MC38_intersect_list = list(Irf4wt_MC38_intersect)
#
dfIrf4ko_MC38_intersect_list = pd.DataFrame(Irf4ko_MC38_intersect_list, columns=['Gene Name'])
dfIrf4ko_MC38_intersect = pd.merge(dfIrf4ko_MC38_intersect_list, dfExpressedMC38vsSpleen, left_on= "Gene Name",right_on="Id")
#
dfIrf4wt_MC38_intersect_list = pd.DataFrame(Irf4wt_MC38_intersect_list, columns=['Gene Name'])
dfIrf4wt_MC38_intersect = pd.merge(dfIrf4wt_MC38_intersect_list,dfExpressedMC38vsSpleen, left_on= "Gene Name",right_on= "Id")
#
dfIrf4ko_MC38_intersect.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4ko_MC38_intersect_annotated.txt", index=False,sep="\t",header=False)
dfIrf4wt_MC38_intersect_list.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4wt_MC38_intersect_annotated.txt", index=False,sep="\t",header=False)
#
# with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4ko_MC38_intersect.txt",'w') as dfIrf4ko_MC38:
#     for i in range(len(dfIrf4ko_MC38_intersect_list)):
#         dfIrf4ko_MC38.write(str(dfIrf4ko_MC38_intersect_list[i]) + "\n")
# dfIrf4ko_MC38.close()
#
# with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4wt_MC38_intersect.txt",'w') as dfIrf4wt_MC38:
#     for i in range(len(dfIrf4wt_MC38_intersect_list)):
#         dfIrf4wt_MC38.write(str(dfIrf4wt_MC38_intersect_list[i]) + "\n")
# dfIrf4wt_MC38.close()
