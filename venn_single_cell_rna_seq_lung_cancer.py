import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn3_circles
import pandas as pd
# set color
up = "red"
down = "blu"
singlecell = "grey"
# set path for import file
scrnaseq= "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure3/signatures_GSE99254.txt"
bulkrnaseq= "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure3/signatures_CCR8pICOSp_vs_CCR8nICOSn.txt"
# read files
dfbulkrnaseq = pd.read_csv(bulkrnaseq, sep="\t", header=0)
dfscrnaseq = pd.read_csv(scrnaseq, sep="\t", header=0)
# cut and convert columns with gene name to series
rnaseq_icosplus_ccr8 = pd.Series(dfbulkrnaseq ['Genes_UP_CCR8pICOSp'], index=dfbulkrnaseq.index)
rnaseq_icosminus_ccr8 = pd.Series(dfbulkrnaseq['Genes_UP_CCR8mICOSm'], index=dfbulkrnaseq.index)
singlecell_ctl4 = pd.Series(dfscrnaseq['CD4_treg_CD4_C9_CTL4'], index=dfscrnaseq.index)
# drop na
rnaseq_icosplus_ccr8.dropna(inplace=True)
rnaseq_icosminus_ccr8.dropna(inplace=True)
singlecell_ctl4.dropna(inplace=True)
# Make the diagram
c =venn3([set(rnaseq_icosplus_ccr8),set(rnaseq_icosminus_ccr8),set(singlecell_ctl4)], ('ICOS+CCR8+ DEGS', 'ICOS-CCR8- DEGS', 'Cluster 9 CTLA4'))
c.get_patch_by_id('100').set_color("white") # irf4 peaks
c.get_patch_by_id('010').set_color("white") # batf peaks
c.get_patch_by_id('001').set_color("white") # degs low
c.get_patch_by_id('101').set_color("white") # intersection irf4 degs

c.get_patch_by_id('011').set_color("white") # intersection chip batf degs

v=venn3_circles(subsets = (1632, 2070, 0, 197, 279, 29, 0),  linewidth=1, color="black")

#
plt.show()
plt.savefig("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure3/venn_plot_ctla4.eps",format='eps')


# # cut and convert columns with gene name to series
# singlecell_foxp3 = pd.Series(dfscrnaseq['CD4_treg_CD4_C8_FOXP3'], index=dfscrnaseq.index)
# singlecell_foxp3.dropna(inplace=True)
# #
# # Make the diagram
# c =venn3([set(rnaseq_icosplus_ccr8),set(rnaseq_icosminus_ccr8),set(singlecell_foxp3)], ('ICOS+CCR8+ DEGS', 'ICOS-CCR8- DEGS', 'Cluster 8 FOXP3'))
# c.get_patch_by_id('100').set_color("white") # irf4 peaks
# c.get_patch_by_id('010').set_color("white") # batf peaks
# c.get_patch_by_id('001').set_color("white") # degs low
# c.get_patch_by_id('101').set_color("white") # intersection irf4 degs
# c.get_patch_by_id('011').set_color("white") # intersection chip batf degs
# #
# v=venn3_circles(subsets = (1908, 1987, 0, 74, 3, 112, 0),  linewidth=1, color="black")
# intersect = set(rnaseq_icosplus_ccr8).intersection(set(singlecell_foxp3))
# plt.show()
# plt.savefig("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure3/venn_plot_foxp3.eps",format='eps')

