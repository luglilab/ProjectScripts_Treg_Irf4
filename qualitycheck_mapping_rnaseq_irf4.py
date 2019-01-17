import sys
import subprocess
import os
import traceback
# set path tools used
# change these path to execute the script
fastqdump = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/fastq-dump"
fastqc = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/fastqc"
multiqc = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/multiqc"
#
thread = 40
maindir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010"
rawdatadir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/raw_data_Nature_Imm_GSE49929"
fastqcdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/fastqc_Nature_Imm_GSE49929"
multiqcdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/multiqc_Nature_Imm_GSE49929"
#
sraid = ["SRR953136","SRR953137","SRR953140","SRR953141","SRR953142"]
fastqreads = ["Wtrep1_r1.fastq","Wtrep1_r2.fastq","Wtrep2_r1.fastq","Wtrep2_r2.fastq","Korep1_r1.fastq","Korep1_r2.fastq",
              "Korep2_r1.fastq","Korep2_r2.fastq","Korep3_r1.fastq","Korep3_r2.fastq"]
# set working directory
os.chdir("/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010")
# create directory function
def createdir(dirpath):
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
        print(" ".join(["Directory", dirpath.split("/")[-1], "Created"]))
    else:
        print(" ".join(["Directory", dirpath.split("/")[-1], "already exists"]))
# create directory
createdir(rawdatadir)
# move in raw read folder
os.chdir(rawdatadir)
# download from sra row data
for i in range(len(sraid)):
    try:
        subprocess.check_call(" ".join([fastqdump,"-I","--split-files",sraid[i],"-O",rawdatadir]),shell=True)
    except (subprocess.CalledProcessError,traceback), e:
        print "Fastq-Dump process fail. Exit"
        sys.exit(1)
    else:
        print "Fastq-Dump process complete."
# rename file
os.rename('SRR953136_1.fastq', 'Wtrep1_r1.fastq')
os.rename('SRR953136_2.fastq', 'Wtrep1_r2.fastq')
os.rename('SRR953137_1.fastq', 'Wtrep2_r1.fastq')
os.rename('SRR953137_2.fastq', 'Wtrep2_r2.fastq')
os.rename('SRR953140_1.fastq', 'Korep1_r1.fastq')
os.rename('SRR953140_2.fastq', 'Korep1_r2.fastq')
os.rename('SRR953141_1.fastq', 'Korep2_r1.fastq')
os.rename('SRR953141_2.fastq', 'Korep2_r2.fastq')
os.rename('SRR953142_1.fastq', 'Korep3_r1.fastq')
os.rename('SRR953142_2.fastq', 'Korep3_r2.fastq')
# fastqc generation
# createdir(fastqcdir)
# # compute fastqc
for i in range(len(fastqreads)):
    try:
        subprocess.check_call(" ".join([fastqc,fastqreads[i],"-o",fastqcdir]),shell=True)
    except (subprocess.CalledProcessError, traceback), e:
        print "FastqQC process fail. Exit"
        sys.exit(1)
    else:
        print "FastqQC process complete."
# create directory
createdir(multiqcdir)
# move to project directory
os.chdir(maindir)
# execute multiqc
try:
    subprocess.check_call(" ".join([multiqc,fastqcdir,"-o",multiqcdir]),shell=True)
except (subprocess.CalledProcessError, traceback), e:
    print "MultiQC process fail.Exit"
    sys.exit(1)
else:
    print "MultiQC process complete."
