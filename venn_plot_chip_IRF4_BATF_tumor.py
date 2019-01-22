import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn3_circles
import pandas as pd
# set color
IRF4 = "#ff6666"
BATF = "#6666ff"
chipcolor = "#b2b2b2"
# set path for import file
chipIRF4 = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/peaks_folders/Treg_IRF4_annotated_peaks_annotated_converted.txt"
chipBATF = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/peaks_folders/Th9_BATF_annotated_peaks_01_annotated_converted.txt"
rnaseqTumor = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/peaks_folders/tumorsdegs.txt"
# read files
dfchipIRF4 = pd.read_csv(chipIRF4, sep="\t", header=0)
dfchipBATF = pd.read_csv(chipBATF, sep="\t", header=0)
dfrnaseqTumor = pd.read_csv(rnaseqTumor, sep="\t", header=0)
# filter degs
dfrnaseqTumor = dfrnaseqTumor[dfrnaseqTumor.padj <= 0.05]
#
genenames_chipIRF4 = pd.Series(dfchipIRF4['human_symbol'], index=dfchipIRF4.index)
genenames_chipBATF = pd.Series(dfchipBATF['human_symbol'], index=dfchipBATF.index)
genenames_rnaseqTumor = pd.Series(dfrnaseqTumor['Gene_name'], index=dfrnaseqTumor.index)
# Make the diagram
c =venn3([set(genenames_chipIRF4),set(genenames_chipBATF),set(genenames_rnaseqTumor)], ('Irf4-Bound (ChIP)', 'BATF-Bound (ChIP)', 'ICOS+CCR8+ vs ICOS-CCR8- (RNA-Seq)'))
c.get_patch_by_id('100').set_color("#FFFFFF") # irf4 peaks
c.get_patch_by_id('010').set_color("#FFFFFF") # batf peaks
c.get_patch_by_id('001').set_color("#FFFFFF") # degs low
c.get_patch_by_id('101').set_color("orange") # intersection irf4 degs
c.get_patch_by_id('111').set_color("purple") # common
c.get_patch_by_id('011').set_color("limegreen") # intersection chip batf degs
c.get_patch_by_id('110').set_color("#FFFFFF")# common peaks
v=venn3_circles(subsets = (1535, 2539, 1014, 1740, 408,284,238),  linewidth=1, color="black")
# list generation
Irf4_tumor_intersect = set(genenames_chipIRF4).intersection(set(genenames_rnaseqTumor)).difference(genenames_chipBATF)
Batf_tumor_intersect = set(genenames_chipBATF).intersection(set(genenames_rnaseqTumor)).difference(genenames_chipIRF4)

Common = set(genenames_chipIRF4).intersection(set(genenames_chipBATF)).intersection(set(genenames_rnaseqTumor))
#
plt.savefig("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/venn_plot_Irf4_Batf_Degs_qvalue_filtered.eps",format='eps')
# conversion list
Irf4_tumor_intersect_list = list(Irf4_tumor_intersect)
Batf_tumor_intersect_list = list(Batf_tumor_intersect)
Common_list = list(Common)
# annotation common / intersection
dfIrf4_tumor_intersect_list =  pd.DataFrame(Irf4_tumor_intersect_list, columns=['Gene Name'])
dfIntesrctionIrf4chipTumor = pd.merge(dfIrf4_tumor_intersect_list,dfchipIRF4, left_on= "Gene Name", right_on= "human_symbol")

#
dfBatf_tumor_intersect_list = pd.DataFrame(Batf_tumor_intersect_list, columns=['Gene Name'])
dfIntesrctionBatfchipTumor = pd.merge(dfBatf_tumor_intersect_list,dfchipBATF, left_on= "Gene Name", right_on= "human_symbol")
#
dfCommon_list = pd.DataFrame(Common_list, columns=['Gene Name'])
dfCommon_Irf4 = pd.merge(dfCommon_list,dfchipIRF4, left_on= "Gene Name", right_on= "human_symbol")
dfCommon_Irf4_Batf = pd.merge(dfCommon_Irf4,dfchipBATF, left_on= "Gene Name_x", right_on= "human_symbol")
#
dfIntesrctionIrf4chipTumor.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4_tumor_intersect_annotated.txt", index=False,sep="\t",header=False)
dfIntesrctionBatfchipTumor.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Batf_tumor_intersect_annotated.txt", index=False,sep="\t",header=False)
dfCommon_Irf4_Batf.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Common_annotated.txt", index=False,sep="\t",header=False)
# export files
with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4_tumor_intersect.txt",'w') as Irf4_tumor_file:
    for i in range(len(Irf4_tumor_intersect_list)):
        Irf4_tumor_file.write(str(Irf4_tumor_intersect_list[i]) + "\n")
Irf4_tumor_file.close()

with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Batf_tumor_intersect.txt",'w') as Batf_tumor_file:
    for i in range(len(Batf_tumor_intersect_list)):
        Batf_tumor_file.write(str(Batf_tumor_intersect_list[i]) + "\n")
Batf_tumor_file.close()

with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Common.txt",'w') as Common_file:
    for i in range(len(Common_list)):
        Common_file.write(str(Common_list[i]) + "\n")
Common_file.close()

