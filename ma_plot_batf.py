import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

filernaseq = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/edgeRresults/tables/CCR8pICOSpvsCCR8mICOSm.complete.txt"
"""
cut -f 1,9 Batf_tumor_intersect_annotated.txt | awk 'BEGIN{FS=OFS="\t"}{split($2,a,"(");print $1,a[1]}' - | sort -k1,1b | sed -e 's/Intergenic/intergenic/g' - > Batf_tumor_intersect_annotated_maplot.txt
cut -f 1,8 Irf4_tumor_intersect_annotated.txt | awk 'BEGIN{FS=OFS="\t"}{split($2,a,"(");print $1,a[1]}' - | sort -k1,1b | sed -e 's/Intergenic/intergenic/g' - > Irf4_tumor_intersect_annotated_maplot.txt
cut -f 1,8 Common_annotated.txt | awk 'BEGIN{FS=OFS="\t"}{split($2,a,"(");print $1,a[1]}' - | sort -k1,1b | sed -e 's/Intergenic/intergenic/g' - > Common_tumor_intersect_annotated_maplot.txt
"""
batf_chip_annotated = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Batf_tumor_intersect_annotated_maplot.txt"
irf4_chip_annotated = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4_tumor_intersect_annotated_maplot.txt"
common_chip_rnaseq = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Common_tumor_intersect_annotated_maplot.txt"
# load data file
RnaseqGene = pd.read_csv(filernaseq,
                sep="\t",
                header=0)

BatfGene = pd.read_csv(batf_chip_annotated,
                       sep="\t",
                       header=None)

Irf4Gene = pd.read_csv(irf4_chip_annotated,
                       sep="\t",
                       header=None)
CommonGene = pd.read_csv(common_chip_rnaseq,
                         sep="\t",
                         header=None
                         )
# Merge RnaSeq Expression with Venn Genes intersect and common
dfRnaseqBatfmerged = pd.merge(RnaseqGene,BatfGene,left_on="Id", right_on= 0, how = 'outer')
dfRnaseqIrf4Batfmerged = pd.merge(dfRnaseqBatfmerged,Irf4Gene,left_on="Id", right_on= 0, how = 'outer')
dfRnaseqCommonmerge = pd.merge(dfRnaseqIrf4Batfmerged,CommonGene,left_on="Id", right_on= 0, how = 'outer')
# filter out non expressed genes
dfRnaseqCommonmerge = dfRnaseqCommonmerge.loc[(dfRnaseqCommonmerge['CCR8mICOSm'] >= 0) & (dfRnaseqCommonmerge['CCR8pICOSp'] >= 0)]
# drop columns irf4 genes batf anno
dfRnaseqCommonmerge = dfRnaseqCommonmerge.drop(['0_x', '1_y','0_y',0,1], axis=1)
# rename columns
dfRnaseqCommonmerge = dfRnaseqCommonmerge.rename(columns={'1_x': 'BATF_target', 2: 'BATF_common'})
# print (dfRnaseqCommonmerge.head())
# color definition
dfRnaseqCommonmerge['color'] = "lightgrey"   # "'#cccccc' # set grey for all expressed genes
dfRnaseqCommonmerge_intergenic = dfRnaseqCommonmerge.copy()
dfRnaseqCommonmerge_promoter = dfRnaseqCommonmerge.copy()
dfRnaseqCommonmerge_exon = dfRnaseqCommonmerge.copy()
dfRnaseqCommonmerge_intron = dfRnaseqCommonmerge.copy()
# color of intergenic
dfRnaseqCommonmerge_intergenic.loc[dfRnaseqCommonmerge_intergenic['BATF_target'] == "intergenic", 'color'] = "limegreen"      # "#ff1919"
dfRnaseqCommonmerge_intergenic.loc[dfRnaseqCommonmerge_intergenic['BATF_common'] == "intergenic", 'color'] = "purple"          #"#774177"
# color of promoter
dfRnaseqCommonmerge_promoter.loc[dfRnaseqCommonmerge_promoter['BATF_target'] == "promoter-TSS ", 'color'] = "limegreen"      # "#ff1919"
dfRnaseqCommonmerge_promoter.loc[dfRnaseqCommonmerge_promoter['BATF_common'] == "promoter-TSS ", 'color'] = "purple"          #"#774177"
# color of intron
dfRnaseqCommonmerge_intron.loc[dfRnaseqCommonmerge_intron['BATF_target'] == "intron ", 'color'] = "limegreen"      # "#ff1919"
dfRnaseqCommonmerge_intron.loc[dfRnaseqCommonmerge_intron['BATF_common'] == "intron ", 'color'] = "purple"          #"#774177"
# color of exon
dfRnaseqCommonmerge_exon.loc[dfRnaseqCommonmerge_exon['BATF_target'] == "exon ", 'color'] = "limegreen"      # "#ff1919"
dfRnaseqCommonmerge_exon.loc[dfRnaseqCommonmerge_exon['BATF_common'] == "exon ", 'color'] = "purple"          #"#774177"
# Now, data is ready for MA plot
# In MA plot, X-axis is log2 normalized mean of expression counts
dfRnaseqCommonmerge_intergenic['A']=np.log2( (dfRnaseqCommonmerge_intergenic['CCR8mICOSm']+ dfRnaseqCommonmerge_intergenic['CCR8pICOSp'])/2 )
dfRnaseqCommonmerge_promoter['A']=np.log2( (dfRnaseqCommonmerge_promoter['CCR8mICOSm']+ dfRnaseqCommonmerge_promoter['CCR8pICOSp'])/2 )
dfRnaseqCommonmerge_intron['A']=np.log2( (dfRnaseqCommonmerge_intron['CCR8mICOSm']+ dfRnaseqCommonmerge_intron['CCR8pICOSp'])/2 )
dfRnaseqCommonmerge_exon['A']=np.log2( (dfRnaseqCommonmerge_exon['CCR8mICOSm']+ dfRnaseqCommonmerge_exon['CCR8pICOSp'])/2 )
# plot
dfRnaseqCommonmerge_intergenic = dfRnaseqCommonmerge_intergenic.sort_values(by=['color'])
dfRnaseqCommonmerge_promoter = dfRnaseqCommonmerge_promoter.sort_values(by=['color'])
dfRnaseqCommonmerge_intron = dfRnaseqCommonmerge_intron.sort_values(by=['color'])
dfRnaseqCommonmerge_exon = dfRnaseqCommonmerge_exon.sort_values(by=['color'])
#
plt.subplot(2, 2, 1)
plt.subplots_adjust(hspace = 0.1)
plt.scatter(dfRnaseqCommonmerge_promoter['A'], dfRnaseqCommonmerge_promoter['log2FoldChange'], c=dfRnaseqCommonmerge_promoter['color'],alpha=0.7,s=5)
plt.title('Promoter BATF')
plt.axhline(y=0, color='b', linestyle='--')
plt.ylabel('M', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
#
plt.subplot(2, 2, 2)
plt.subplots_adjust(hspace = 0.5)
plt.scatter(dfRnaseqCommonmerge_intergenic['A'], dfRnaseqCommonmerge_intergenic['log2FoldChange'], c=dfRnaseqCommonmerge_intergenic['color'],alpha=0.7,s=5)
plt.axhline(y=0, color='b', linestyle='--')
plt.ylabel('M', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Intergenic BATF')
#
plt.subplot(2, 2, 3)
plt.subplots_adjust(hspace = 0.1)
plt.scatter(dfRnaseqCommonmerge_intron['A'], dfRnaseqCommonmerge_intron['log2FoldChange'], c=dfRnaseqCommonmerge_intron['color'],alpha=0.7,s=5)
plt.title('Intron BATF')
plt.axhline(y=0, color='b', linestyle='--')
plt.xlabel('A',fontsize=10)
plt.ylabel('M', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
#
plt.subplot(2, 2, 4)
plt.subplots_adjust(hspace = 0.5)
plt.scatter(dfRnaseqCommonmerge_exon['A'], dfRnaseqCommonmerge_exon['log2FoldChange'], c=dfRnaseqCommonmerge_exon['color'],alpha=0.7,s=5)
plt.axhline(y=0, color='b', linestyle='--')
plt.xlabel('A',fontsize=10)
plt.ylabel('M', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Exon BATF')
#
plt.savefig("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/MAplot_intergenic_promoter_BATF.eps", format='eps')
#plt.show()
