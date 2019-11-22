import pandas as pd
import os

irf4 = "/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/Irf4_tumor_intersect.txt"
batf = "/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/Batf_tumor_intersect.txt"
common = "/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/Common.txt"
rnaseqexpression = "/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/MC38vsSpleen.complete.anno.txt"
#
dfirf4 = pd.read_csv(irf4,sep="\t",header=None)
dfbatf = pd.read_csv(batf,sep="\t",header=None)
dfcommon = pd.read_csv(common,sep="\t",header=None)
dfrnaseqexpression = pd.read_csv(rnaseqexpression,sep="\t",header=0)
#
print dfirf4.head()
print dfbatf.head()
print dfcommon.head()
print dfrnaseqexpression.head()
#
dfmergeIrf4  = pd.merge(dfirf4,dfrnaseqexpression,left_on=0,right_on="Id")
dfmergeBatf = pd.merge(dfbatf,dfrnaseqexpression,left_on=0,right_on="Id")
dfmergeCommon = pd.merge(dfcommon,dfrnaseqexpression,left_on=0,right_on="Id")
#
Irf4FCpos = dfmergeIrf4[dfmergeIrf4.log2FoldChange >= 0]
Irf4FCneg = dfmergeIrf4[dfmergeIrf4.log2FoldChange < 0]
BatfFCpos = dfmergeBatf[dfmergeBatf.log2FoldChange >= 0]
BatfFCneg = dfmergeBatf[dfmergeBatf.log2FoldChange < 0]
CommonFCpos = dfmergeCommon[dfmergeCommon.log2FoldChange >= 0]
CommonFCneg = dfmergeCommon[dfmergeCommon.log2FoldChange < 0]
#
Irf4FCpos[["Id"]].to_csv("/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/Irf4FCpos.txt",sep="\t",header=False,index=False)
Irf4FCneg[["Id"]].to_csv("/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/Irf4FCneg.txt",sep="\t",header=False,index=False)
BatfFCpos[["Id"]].to_csv("/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/BatfFCpos.txt",sep="\t",header=False,index=False)
BatfFCneg[["Id"]].to_csv("/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/BatfFCneg.txt",sep="\t",header=False,index=False)
CommonFCpos[["Id"]].to_csv("/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/CommonFCpos.txt",sep="\t",header=False,index=False)
CommonFCneg[["Id"]].to_csv("/Users/simonepuccio/Documents/HumanitasProjects/SP010/Figure5/CommonFCneg.txt",sep="\t",header=False,index=False)
