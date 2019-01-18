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
batfko_degsup = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/BATFKO_vsWT.up.txt"
batfko_degsdown = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/BATFKO_vsWT.down.txt"
chiptargetBatf = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/peaks_folders/Th9_BATF_annotated_peaks_01_annotated.txt"
# import degs microarray Batf ko
dfbatfko_degsup = pd.read_csv(batfko_degsup,sep="\t",header=0)
dfbatfko_degsdown = pd.read_csv(batfko_degsdown,sep="\t",header=0)
# import chip peaks annoted
dfchiptargetBatf = pd.read_csv(chiptargetBatf,sep="\t",header=0)
# cut and convert columns with gene name to series
genenames_degsupbatf = pd.Series(dfbatfko_degsup['symbol'], index=dfbatfko_degsup.index)
genenames_degsdownbatf = pd.Series(dfbatfko_degsdown['symbol'], index=dfbatfko_degsdown.index)
genenames_chipbatf = pd.Series(dfchiptargetBatf['Gene Name'], index=dfchiptargetBatf.index)
# # Make the diagram
batf =venn3([set(genenames_degsupbatf),set(genenames_degsdownbatf),set(genenames_chipbatf)], ('Batf_repressed_genes', 'Batf_activated_genes', 'ChipTarget'))
batf.get_patch_by_id('100').set_color(down)
batf.get_patch_by_id('010').set_color(up)
batf.get_patch_by_id('001').set_color(chipcolor)
up_chip = set(genenames_degsupbatf).intersection(set(genenames_chipbatf))
down_chip = set(genenames_degsdownbatf).intersection(set(genenames_chipbatf))
# plt.show()
plt.savefig("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/venn_plot_batf.eps",format='eps')
#
up_chip_list = list(up_chip)
down_chip_list = list(down_chip)
with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Batf_up_chip.txt",'w') as up_chip_file:
    for i in range(len(up_chip_list)):
        up_chip_file.write(str(up_chip_list[i]) + "\n")
up_chip_file.close()

with open("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Batf_down_chip.txt",'w') as down_chip_file:
    for i in range(len(down_chip_list)):
        down_chip_file.write(str(down_chip_list[i]) + "\n")
down_chip_file.close()