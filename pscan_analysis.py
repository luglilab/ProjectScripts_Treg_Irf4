import os
import sys
import optparse
import pandas as pd
import mygene
import rpy2
from pybiomart import Server
from Bio import SeqIO
import subprocess
#
server = Server(host='http://www.ensembl.org')
dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'])
conversionscript = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/ProjectScripts/utility.py"
# PATH input
EdgeResultsfolder= "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/edgeRresults/tables/"
# PATH output
OutputPAth = "/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/TestPscan/"
# Files name
UPdegs = "CCR8pICOSpvsCCR8mICOSm.up.txt"
DOWNdegs = "CCR8pICOSpvsCCR8mICOSm.down.txt"
expressed = "CCR8pICOSpvsCCR8mICOSm.complete.txt"
fastapromotor = "/home/spuccio/pscan/Human_and_Mouse/Homo_sapiens_950up_50down.fasta"
# Conversion table
annodb = dataset.query(attributes=['external_gene_name','refseq_mrna'])
# remove Gene Symbol with no overlap with Refseq ID
filtered_annodb = annodb[annodb[ u'RefSeq mRNA ID'].notnull()].reset_index(drop=True)
# import to df
dfUP = pd.read_csv("".join([EdgeResultsfolder,UPdegs]),sep="\t",header=0)
dfDOWN = pd.read_csv("".join([EdgeResultsfolder,DOWNdegs]),sep="\t",header=0)
# check duplicates ID
dfUP.drop_duplicates(subset ="Id", inplace = True)
dfDOWN.drop_duplicates(subset ="Id", inplace = True)
# extract gene symbol
dfUP_symb = dfUP[['Id']]
dfDOWN_symb = dfDOWN[['Id']]
#
print "The number of DEGS with positive log2Fold change is %d " % len(dfUP)
print "The number of DEGS with negative log2Fold change is %d " % len(dfDOWN)
# Id conversion UP
dfUP_converted = pd.merge(dfUP_symb,filtered_annodb,left_on=u'Id',right_on=u'Gene name')
dfUP_converted.drop_duplicates(subset ="Id", inplace = True)
dfUP_converted.to_csv("".join([OutputPAth,"Up_regulatedGenesRefseqConverted.txt"]),sep="\t",index=False)
# Id Conversion DOWN
dfDOWN_converted = pd.merge(dfDOWN_symb,filtered_annodb,left_on=u'Id',right_on=u'Gene name')
dfDOWN_converted.drop_duplicates(subset ="Id", inplace = True)
dfDOWN_converted.to_csv("".join([OutputPAth,"Down_regulatedGenesRefseqConverted.txt"]),sep="\t",index=False)
print "The number of DEGS with positive log2Fold change after conversion is %d " % len(dfUP_converted)
print "The number of DEGS with negative log2Fold change after conversion is %d " % len(dfDOWN_converted)
# filter Jaspar Database with RNA-Seq Expressed Genes
# import expressed genes
dfExpressed = pd.read_csv("".join([EdgeResultsfolder,expressed]),sep="\t",header=0)
# filter non expressed genes
dfExpressed = dfExpressed[dfExpressed[ u'pvalue'].notnull()].reset_index(drop=True)
# merge with conversion table
dfExpressedConverted = pd.merge(dfExpressed,filtered_annodb,left_on=u'Id',right_on=u'Gene name')
# drop duplicated
dfExpressedConverted.drop_duplicates(subset ="Id", inplace = True)
# reset index
dfExpressedConverted = dfExpressedConverted.reset_index()
# cut columns
dfExpressedConverted = dfExpressedConverted[['Id','RefSeq mRNA ID']]
# converte fasta to tab
try:
    subprocess(" ".join(["python",conversionscript,"-i", fastapromotor, "-o", "Homo_sapiens_950up_50down", "-p", OutputPAth, "-t", "fasta2tab", "-r 1"]),shell=True)
except:
    print "Complete"
else:
    print "Complete"
# import the tabular file
dfPromoterTab = pd.read_csv("".join([OutputPAth,"Homo_sapiens_950up_50down_1.tab"]), sep="\t",header=None)
# split fasta sequence id
dfPromoterTab['assembly'],dfPromoterTab['consortium'],dfPromoterTab['prefix'], dfPromoterTab['RefseqGenePromotor'] = dfPromoterTab[0].str.split('_', 3).str
# join string to reassembly the refseq id associated
dfPromoterTab['RefseqId'] = dfPromoterTab[['prefix', 'RefseqGenePromotor']].apply(lambda x: '_'.join(x), axis=1)
# delete temp columns
dfPromoterTab = dfPromoterTab.drop(['assembly','consortium','prefix','RefseqGenePromotor'], axis=1)
# create tab filtered for expressed
dfpromoter_expressed = pd.merge(dfExpressedConverted,dfPromoterTab,left_on="RefSeq mRNA ID",right_on="RefseqId")
dfpromoter_expressed = dfpromoter_expressed[[0,1]]
dfpromoter_expressed.to_csv(OutputPAth+"promoter_expressed.tab",sep="\t",header=False,index=False)
subprocess(" ".join(["python",conversionscript,"-i", OutputPAth+"promoter_expressed.tab", "-o", "Homo_sapiens_950up_50down_promoter_expressed", "-p", OutputPAth, "-t", "tab2fasta"]),shell=True)
try:
    subprocess(" ".join(["python",conversionscript,"-i", OutputPAth+"promoter_expressed.tab", "-o", "Homo_sapiens_950up_50down_promoter_expressed", "-p", OutputPAth, "-t", "tab2fasta"]),shell=True)
except:
    print "Complete"
else:
    print "Complete"
print dfPromoterTab.head()