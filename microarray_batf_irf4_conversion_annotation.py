import os
import pandas as pd
# set file path
degsBATFmouse = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/BATFko_BATFwt_mouse.txt" #pvalue < 0.05
degsIRF4mouse = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/IRF4ko_IRF4wt_mouse.txt" #pvalue < 0.05
conversiontable = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/human_mouse_hcop_fifteen_column2.txt"
# df creation
dfBATFmouse = pd.read_csv(degsBATFmouse,sep="\t",header=0)
dfIRF4mouse = pd.read_csv(degsIRF4mouse, sep="\t",header=0)
dfconversiontable = pd.read_csv(conversiontable, sep="\t",header=0)
# filter only NCBI converted genes in support columns
dfconversiontable = dfconversiontable[dfconversiontable['support'].str.contains("NCBI")]
# Merging
dfmergeBatf = pd.merge(dfBATFmouse,dfconversiontable,left_on="symbol",right_on="mouse_symbol") # lose 264 genes
dfmergeIrf4 = pd.merge(dfIRF4mouse,dfconversiontable,left_on="symbol",right_on="mouse_symbol") # lose 116 genes
# Subsetting
dfDEGSBatfHuman =  dfmergeBatf[[u'symbol', u'logFC',u'human_symbol']]
dfDEGSIrf4Human =  dfmergeIrf4[[u'symbol', u'logFC',u'human_symbol']]
# file path set
Irf4_tumor_intersect = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4_tumor_intersect.txt"
Batf_tumor_intersect = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Batf_tumor_intersect.txt"
Common = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Common.txt"
# df creation
dfIrf4_tumor_intersect = pd.read_csv(Irf4_tumor_intersect,sep="\t",header=None)
dfBatf_tumor_intersect = pd.read_csv(Batf_tumor_intersect,sep="\t",header=None)
dfCommon = pd.read_csv(Common,sep="\t",header=None)
# intersection batf and irf4
dfmergeIrf4_tumor_intersectwithDEGS = pd.merge(dfIrf4_tumor_intersect,dfDEGSIrf4Human,left_on=0,right_on=u'human_symbol')
dfmergeBatf_tumor_intersectwithDEGS = pd.merge(dfBatf_tumor_intersect,dfDEGSBatfHuman,left_on=0,right_on=u'human_symbol')
# intersection common
dfmergeCommon_intersectwithDEGSIRF4 = pd.merge(dfCommon,dfDEGSIrf4Human,left_on=0,right_on=u'human_symbol')
dfmergeCommon_intersectwithDEGSBATF = pd.merge(dfCommon,dfDEGSBatfHuman,left_on=0,right_on=u'human_symbol')
#
dfmergeIrf4_tumor_intersectwithDEGS.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Irf4_unique_directed_target.txt",sep="\t",header=False,index=False)
dfmergeBatf_tumor_intersectwithDEGS.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Batf_unique_directed_target.txt",sep="\t",header=False,index=False)
dfmergeCommon_intersectwithDEGSIRF4.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Common_and_Irf4DESG_target.txt",sep="\t",header=False,index=False)
dfmergeCommon_intersectwithDEGSBATF.to_csv("/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/Figure4/Common_and_BatfDESG_target.txt",sep="\t",header=False,index=False)


