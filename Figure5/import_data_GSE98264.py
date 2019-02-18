import sys
import os
import subprocess
# assign variables
rawdatadir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/raw_data_GSE98264"
fastqdump = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/fastq-dump"
SRR_Acc_list = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/SRR_Acc_list_GSE98264.txt"


def download_fastqdump():
    """
    Execute fastq-dump
    :return:
    """
    # set working directory
    os.chdir(rawdatadir)
    # download raw-data
    try:
        with open(SRR_Acc_list) as f:
            for line in f:
                subprocess.check_call(fastqdump + ' -I --split-files '+line.strip("/n"), shell=True)
    except:
        print('Fastq-Dump process fail. Exit')
        sys.exit(1)
    else:
        print('Fastq-Dump process complete')


if __name__ == "__main__":
    # create folder
    if not os.path.exists(rawdatadir):
        os.mkdir(rawdatadir)
        print(" ".join(["Directory", rawdatadir.split("/")[-1], "Created"]))
        download_fastqdump()
    else:
        print(" ".join(["Directory", rawdatadir.split("/")[-1], "already exists"]))
        download_fastqdump()
