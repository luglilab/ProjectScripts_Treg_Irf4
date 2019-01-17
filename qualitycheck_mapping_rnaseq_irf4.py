import sys
import subprocess
import os
# set path tools used
paralleldump = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/parallel-fastq-dump"
#
thread = 40
rawdatadir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/raw_data_Nature_Imm_GSE49929"
sraid = ["SRR953136","SRR953137","SRR953140","SRR953141","SRR953142"]
# set working directory
os.chdir("/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010")
# create directory
if not os.path.exists(rawdatadir):
    os.mkdir(rawdatadir)
    print(" ".join(["Directory", rawdatadir.split("/")[-1], "Created"]))
else:
    print(" ".join(["Directory", rawdatadir.split("/")[-1], "already exists"]))
# move in raw read folder
os.chdir(rawdatadir)
# download from sra row data

for i in range(len(sraid)):
    try:
        subprocess.check_call("".join([paralleldump,"--sra-id",sraid[i],"--threads",str(thread),"--outdir",rawdatadir]),
                              shell=True)

