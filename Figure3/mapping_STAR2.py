import sys
import subprocess
import os
import traceback
#
star = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/STAR" # path of STAR alignment
multiqc = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/multiqc" # path of multiqc tool
# set variables
thread = 40 # number of processor
maindir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010" # home of project folder
starindex = "/home/spuccio/AnnotationSTAR/GencodeMusMusculus/STAR/89/" # path of Star genome index
rawdatadir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/raw_data_Nature_Imm_GSE49929/" # raw reads path
mappingdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/star_Nature_Imm_GSE49929" # BAM path
fastqcdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/fastqc_Nature_Imm_GSE49929" # Fastqc path
multiqcdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/multiqc_Nature_Imm_GSE49929" # multiqc output path
STARparameters = " ".join(["--runThreadN",str(thread),"--sjdbOverhang","89","--sjdbGTFtagExonParentGene","gene_id",
                           "--outFilterMultimapNmax","1","--outSAMtype","BAM","SortedByCoordinate",
                          "--outFilterIntronMotifs","RemoveNoncanonical","--quantMode","GeneCounts",
                          "--outFilterScoreMinOverLread","0","--outFilterMatchNminOverLread","0",
                           "--outFilterMatchNmin","0"]) # STAR parameters
#
genomefasta = "/home/spuccio/AnnotationSTAR/GencodeMusMusculus/GRCm38.chr.genome.fa" # Fasta file path
genomegtf = "/home/spuccio/AnnotationSTAR/GencodeMusMusculus/gencode.vM19.annotation.gtf" # annotation file path
# list
fastqreads = [["Wtrep1_r1.fastq","Wtrep1_r2.fastq"],
              ["Wtrep2_r1.fastq","Wtrep2_r2.fastq"],
              ["Korep1_r1.fastq","Korep1_r2.fastq"],
              ["Korep2_r1.fastq","Korep2_r2.fastq"],
              ["Korep3_r1.fastq","Korep3_r2.fastq"]] #paired end modality
# create directory function
def createdir(dirpath):
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
        print(" ".join(["Directory", dirpath.split("/")[-1], "Created"]))
    else:
        print(" ".join(["Directory", dirpath.split("/")[-1], "already exists"]))
# set working directory
os.chdir(maindir)
# create index directoty
createdir(starindex)
# genome index
try:
    subprocess.check_call(" ".join([star,"--runThreadN",str(thread),"--runMode","genomeGenerate","--genomeDir",starindex,"--genomeFastaFiles",
                                    genomefasta,"--sjdbGTFfile",genomegtf,"--sjdbOverhang","89"]),shell=True)
except (subprocess.CalledProcessError, traceback), e:
    print "STAR index of %s fail. Exit" % genomefasta
    sys.exit(1)
else:
    print "STAR index complete."
# create output folder for BAM files
createdir(mappingdir)
# move to mapping dir
os.chdir(mappingdir)
# mapping
for i in range(len(fastqreads)):
    try:
        subprocess.check_call(" ".join([star,"--genomeDir",starindex,"--readFilesIn",rawdatadir+fastqreads[i][0],rawdatadir+fastqreads[i][1],
                                        "--outFileNamePrefix",fastqreads[i][0].split("_")[0],STARparameters,
                                        "--sjdbGTFfile",genomegtf]),shell=True)
    except (subprocess.CalledProcessError, traceback), e:
        print "STAR mapping of %s fail. Exit" % fastqreads[i][0].split("_")[0]
        sys.exit(1)
    print "STAR mapping of %s complete." % fastqreads[i][0].split("_")[0]
# move to main folder
os.chdir(maindir)
# execute multiqc
try:
    subprocess.check_call(" ".join([multiqc,fastqcdir,mappingdir,"-o",multiqcdir]),shell=True)
except (subprocess.CalledProcessError, traceback), e:
    print "MultiQC process fail.Exit"
    sys.exit(1)
else:
    print "MultiQC process complete."
