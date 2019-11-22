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
Batf4ko_wt_pval005 = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/BATFko_BATFwt_mouse.Pval.annotated.txt"
# create dataframe
dfBatf4ko_wt_pval005 = pd.read_csv(Batf4ko_wt_pval005,sep="\t",header=0)
# drop duplicates
dfBatfko_wt_pval005.drop_duplicates(subset = "symbol", inplace = True)
# create two duplicates for logFC<0 and logFC>=0
dfBatfko = dfBatf4ko_wt_pval005[(dfBatf4ko_wt_pval005['logFC'] >= 0)]
dfBatfwt = dfBatf4ko_wt_pval005[(dfBatf4ko_wt_pval005['logFC'] < 0)]
# convert df to series Batf
genenames_Batfko = pd.Series(dfBatfko['symbol'], index=dfBatfko.index)
genenames_Batfwt = pd.Series(dfBatfwt['symbol'], index=dfBatfwt.index)
###
c =venn3([set(genenames_Batfko),set(genenames_Batfwt),set(genenames_MC38)], ('Batf KO Degs', 'Batf WT Degs', 'TITR MC38'))
c.get_patch_by_id('100').set_color("#FFFFFF") # irf4 ko peaks
c.get_patch_by_id('010').set_color("#FFFFFF") # irf4 wt peaks
c.get_patch_by_id('001').set_color("#FFFFFF") # degs low
c.get_patch_by_id('101').set_color("orange") # intersection irf4 degs
c.get_patch_by_id('011').set_color("limegreen") # intersection chip batf degs
v=venn3_circles(subsets = (693, 744, 0, 298, 20,64,0),  linewidth=1, color="black")
#
Batfko_MC38_intersect = set(genenames_Batfko).intersection(set(genenames_MC38))
Batfwt_MC38_intersect = set(genenames_Batfwt).intersection(set(genenames_MC38))
#
plt.savefig("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/VennPlotBATFKO_BATFWT_382DEGS.eps",format='eps')

Batfko_MC38_intersect_list = list(Batfko_MC38_intersect)
Batfwt_MC38_intersect_list = list(Batfwt_MC38_intersect)
#
Batfko_MC38_intersect_list =  pd.DataFrame(Batfko_MC38_intersect_list, columns=['Gene Name'])
dfBatfko_MC38_intersect_list = pd.merge(Batfko_MC38_intersect_list,dfExpressedMC38vsSpleen, left_on= "Gene Name",right_on="Id")
#
Batfwt_MC38_intersect_list =  pd.DataFrame(Batfwt_MC38_intersect_list, columns=['Gene Name'])
dfBatfwt_MC38_intersect_list = pd.merge(Batfwt_MC38_intersect_list,dfExpressedMC38vsSpleen, left_on= "Gene Name",right_on="Id")
#
dfBatfko_MC38_intersect_list.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Batfko_MC38_intersect_annotated.txt", index=False,sep="\t",header=False)
dfBatfwt_MC38_intersect_list.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Batfwt_MC38_intersect_annotated.txt", index=False,sep="\t",header=False)
#
# with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4ko_MC38_intersect.txt",'w') as dfIrf4ko_MC38:
#     for i in range(len(Irf4ko_MC38_intersect_list)):
#         dfIrf4ko_MC38.write(str(Irf4ko_MC38_intersect_list[i]) + "\n")
# dfIrf4ko_MC38.close()
#
# with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4wt_MC38_intersect.txt",'w') as dfIrf4wt_MC38:
#     for i in range(len(Irf4wt_MC38_intersect_list)):
#         dfIrf4wt_MC38.write(str(Irf4wt_MC38_intersect_list[i]) + "\n")
# dfIrf4wt_MC38.close()
