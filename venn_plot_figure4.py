import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import pandas as pd
# set color
up = "#ff6666"
down = "#6666ff"
chipcolor = "#b2b2b2"
# set path for import file
irf4ko_degsup = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/EdgeR_Nature_Imm_GSE49929/tables/IRF4KOvsWT.up.txt"
irf4ko_degsdown = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/EdgeR_Nature_Imm_GSE49929/tables/IRF4KOvsWT.down.txt"
chiptargetIrf4 = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/peaks_folders/Treg_IRF4_annotated_peaks.txt"
# import degs rnaseq Irf4 ko
dfirf4ko_degsup = pd.read_csv(irf4ko_degsup,sep="\t",header=0)
dfirf4ko_degsdown = pd.read_csv(irf4ko_degsdown,sep="\t",header=0)
# import chip peaks annoted
dfchiptargetIrf4 = pd.read_csv(chiptargetIrf4,sep="\t",header=0)
# cut and convert columns with gene name to series
genenames_degsup = pd.Series(dfirf4ko_degsup['Id'], index=dfirf4ko_degsup.index)
genenames_degsdown = pd.Series(dfirf4ko_degsdown['Id'], index=dfirf4ko_degsdown.index)
genenames_chipIrf4 = pd.Series(dfchiptargetIrf4['Gene Name'], index=dfchiptargetIrf4.index)
# # Make the diagram
c =venn3([set(genenames_degsup),set(genenames_degsdown),set(genenames_chipIrf4)], ('Irf4_repressed_genes', 'Irf4_activated_genes', 'ChipTarget'))
c.get_patch_by_id('100').set_color(down)
c.get_patch_by_id('010').set_color(up)
c.get_patch_by_id('001').set_color(chipcolor)
up_chip = set(genenames_degsup).intersection(set(genenames_chipIrf4))
down_chip = set(genenames_degsdown).intersection(set(genenames_chipIrf4))
# plt.show()
plt.savefig("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/venn_plot_irf4.eps",format='eps')
#
up_chip_list = list(up_chip)
down_chip_list = list(down_chip)
with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4_up_chip.txt",'w') as up_chip_file:
    for i in range(len(up_chip_list)):
        up_chip_file.write(str(up_chip_list[i]) + "\n")
up_chip_file.close()

with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4_down_chip.txt",'w') as down_chip_file:
    for i in range(len(down_chip_list)):
        down_chip_file.write(str(down_chip_list[i]) + "\n")
down_chip_file.close()

