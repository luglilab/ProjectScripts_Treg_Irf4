import sys
import os
import subprocess
#
fastqdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/raw_data_Nature_C_GSE99167"
fastqcout = "/home/spuccio/datadisk2/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/FASTQ_output/"
fastqc = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/fastqc"
os.chdir(fastqdir)

# IRF4
try:
    subprocess.check_call(
        'cat SRR5582849_1.fastq SRR5582850_1.fastq SRR5582851_1.fastq SRR5582852_1.fastq > Th9_IRF4_r1.fastq',
        shell=True)
    subprocess.check_call(
        'cat SRR5582849_2.fastq SRR5582850_2.fastq SRR5582851_2.fastq SRR5582852_2.fastq > Th9_IRF4_r2.fastq',
        shell=True)
except:
    print 'Cat of IRF4 failed.Exit'
    sys.exit(1)
else:
    print 'Cat of IRF4 complete.'
# BATF
try:
    subprocess.check_call(
        'cat SRR5582853_1.fastq SRR5582854_1.fastq SRR5582855_1.fastq SRR5582856_1.fastq > Th9_BATF_r1.fastq',
        shell=True)
    subprocess.check_call(
        'cat SRR5582853_2.fastq SRR5582854_2.fastq SRR5582855_2.fastq SRR5582856_2.fastq > Th9_BATF_r2.fastq',
        shell=True)
except:
    print 'Cat of BATF failed. Exit'
    sys.exit(1)
else:
    print 'Cat of BATF complete.'
# INPUT
try:
    subprocess.check_call(
        'cat SRR5582865_1.fastq SRR5582866_1.fastq SRR5582867_1.fastq SRR5582868_1.fastq > Th9_INPUT_r1.fastq',
        shell=True)
    subprocess.check_call(
        'cat SRR5582865_2.fastq SRR5582866_2.fastq SRR5582867_2.fastq SRR5582868_2.fastq > Th9_INPUT_r2.fastq',
        shell=True)
except:
    print 'Cat of INPUT failed. Exit'
    sys.exit(1)
else:
    print 'Cat of INPUT complete.'

readsamplename = ['Th9_IRF4_r1.fastq','Th9_IRF4_r2.fastq','Th9_BATF_r1.fastq','Th9_BATF_r2.fastq','Th9_INPUT_r1.fastq','Th9_INPUT_r2.fastq']


def fastqcfunc(fastqfiles):
    """
    loop over fastqfiles list and execute fastqc
    :param fastqfiles:
    :return:
    """
    for i in range(len(fastqfiles)):
        try:
            subprocess.check_call(" ".join([fastqc, fastqfiles[i], "-o", fastqcout]),
                                      shell=True)
        except:
            print 'Fastqc of %s failed. Exit' %(fastqfiles[i])
            sys.exit(1)
        else:
            print 'Fastqc of %s complete.' %(fastqfiles[i])


if not os.path.exists(fastqcout):
    os.mkdir(fastqcout)
    print(" ".join(["Directory", fastqcout.split("/")[-1], "Created"]))
    fastqcfunc(readsamplename)
else:
    print(" ".join(["Directory", fastqcout.split("/")[-1], "already exists"]))
    fastqcfunc(readsamplename)
