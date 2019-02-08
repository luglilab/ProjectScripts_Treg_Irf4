import pandas as pd
import xlsxwriter

# Mus musculus tab file with 3 columns "Gene stable ID version"-"Gene name"-"Gene description"
annoMusMusculus = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/EdgeR_Nature_Imm_GSE49929/GRCm38_p6_annotation.txt"
# set path
pathinput = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/EdgeR_Pnas_GSE116347/tables/"
# Input table
outputEdgeR = {"B16vsSpleen": "B16vsSpleen.complete.txt",
               "MC38vsSpleen": "MC38vsSpleen.complete.txt",
               "MC38vsB16": "MC38vsB16.complete.txt"}
# Read Annotation Table
dfannoMusMusculus = pd.read_csv(annoMusMusculus, sep="\t", header=0)
#
for key, value in outputEdgeR.items():
    dfoutputEdgeR = pd.read_csv("".join([pathinput,value]), sep="\t", header=0)
    dfmerged = pd.merge(dfoutputEdgeR, dfannoMusMusculus, left_on=u'Id', right_on=u'Gene name')
    dfmerged = dfmerged.drop([u'Gene stable ID version'], axis=1)
    dfmerged['Id'] = dfmerged['Gene name']
    dfmerged.to_csv(
        "".join([pathinput, key, ".complete.anno.txt"]), sep="\t", header=True)
    dfmerged.dropna(inplace=True)
    dfMerge = dfmerged
    dfPval = dfMerge[dfMerge['padj'] < 0.05]
    dfPval = dfPval.sort_values(by='padj', ascending=True)
    dfMin = dfPval[dfPval['log2FoldChange'] < 0]
    dfMin = dfMin.sort_values(by='padj', ascending=True)
    dfMax = dfPval[dfPval['log2FoldChange'] >= 0]
    dfMax = dfMax.sort_values(by='padj', ascending=True)
    writer = pd.ExcelWriter(pathinput + '/' + key + '.xlsx', engine='xlsxwriter')
    dfMerge.to_excel(writer, sheet_name=key + '.All', index=False)
    dfPval.to_excel(writer, sheet_name=key + '.Qval', index=False)
    dfMax.to_excel(writer, sheet_name=key + '.Max', index=False)
    dfMin.to_excel(writer, sheet_name=key + '.Min', index=False)
    writer.save()
