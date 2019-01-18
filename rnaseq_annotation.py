import os
import pandas as pd
import subprocess
import xlsxwriter

annoMusMusculus = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/EdgeR_Nature_Imm_GSE49929/GRCm38_p6_annotation.txt"
outputEdgeR = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/EdgeR_Nature_Imm_GSE49929/tables/IRF4KOvsWT.complete.txt"
pathoout = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/EdgeR_Nature_Imm_GSE49929/tables/"
dfannoMusMusculus = pd.read_csv(annoMusMusculus,sep="\t",header=0)
dfoutputEdgeR = pd.read_csv(outputEdgeR,sep="\t",header=0)

dfmerged = pd.merge(dfoutputEdgeR,dfannoMusMusculus,left_on=u'Id',right_on=u'Gene stable ID version')
dfmerged = dfmerged.drop([u'Gene stable ID version'], axis=1)
dfmerged['Id'] = dfmerged['Gene name']
dfmerged.to_csv("/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/EdgeR_Nature_Imm_GSE49929/tables/IRF4KOvsWT.complete.anno.txt",
                sep="\t",header=True)
idout = "IRF4KOvsWT"
dfmerged.dropna(inplace=True)
dfMerge = dfmerged
dfPval = dfMerge[dfMerge['padj'] < 0.05]
dfPval = dfPval.sort_values(by='padj', ascending=True)
dfMin = dfPval[dfPval['log2FoldChange'] < 0]
dfMin = dfMin.sort_values(by='padj', ascending=True)
dfMax = dfPval[dfPval['log2FoldChange'] >= 0]
dfMax = dfMax.sort_values(by='padj', ascending=True)
writer = pd.ExcelWriter(pathoout + '/' + idout + '.xlsx', engine='xlsxwriter')
dfMerge.to_excel(writer, sheet_name=idout + '.All', index=False)
dfPval.to_excel(writer, sheet_name=idout + '.padj', index=False)
dfMax.to_excel(writer, sheet_name=idout + '.Max', index=False)
dfMin.to_excel(writer, sheet_name=idout + '.Min', index=False)
writer.save()