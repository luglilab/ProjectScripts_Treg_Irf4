import os
import subprocess
import sys
# set path tool
threads = "20"
MAX_FRAG_LEN = "2000"
indexname = "GRCm38"
# the multimapping flag
multimap = "4"
# folder with Mus musculus (house mouse) genome assembly GRCm38 (mm10)
GRCm38indexpath = "/home/spuccio/AnnotationBowtie2/Mus_musculus/GRCm38.p6/"
# bowtie2-build path
bowtiebuild2path = "/home/spuccio/miniconda3/envs/chipseq_env/bin/bowtie2-build"
# bowtie2 path
bowtie2path = "/home/spuccio/miniconda3/envs/chipseq_env/bin/bowtie2"
# path folder
projectdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/"
# path out file sam
mappingout = "".join([projectdir, "bowtie2_mapping_GSE98264"])
# input file path
raw_data_dir = "".join([projectdir, "raw_data_GSE98264/"])
# Path
GRCm38fasta = "".join([GRCm38indexpath, "/GRCm38.primary_assembly.fa"])
raw_fastq = {"Treg_Irf4_r1": ["SRR5483021_1.fastq", "SRR5483021_2.fastq"],
             "Treg_Irf4_r2": ["SRR5483022_1.fastq", "SRR5483022_2.fastq"]}
raw_fastq2 = {"Treg_Irf4_r3": "SRR5483020_1.fastq"}


def createdir(dirpath):
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
        print(" ".join(["Directory", dirpath.split("/")[-1], "Created"]))
    else:
        print(" ".join(["Directory", dirpath.split("/")[-1], "already exists"]))

# check index function


def checkindex(indexpath, fastafile, indexname):
    if os.path.isfile("".join([indexpath, indexname, ".1.bt2"])) == True:
        print("Genome index of Fasta file already exists.")
    else:
        try:
            os.chdir(GRCm38indexpath)
            subprocess.check_call(" ".join([bowtiebuild2path, "--threads", threads, fastafile, indexname]),
                                  shell=True)
        except subprocess.CalledProcessError:
            print("ERROR.Fastqc analysis failed. Stop execution.")
            sys.exit(1)
        else:
            print("Fastq analysis complete.")
    return indexname


def bowtie2mappingpairedend(indexname, fastqname, samname):
    try:
        subprocess.check_call(" ".join([bowtie2path, "-p", threads, "-q", "--local", "-k", multimap,
                                        "-x", indexname, "-X", MAX_FRAG_LEN,
                                        "-1", "".join([raw_data_dir, fastqname[0]]),
                                        "-2", "".join([raw_data_dir, fastqname[1]]),
                                        "-S", "".join([mappingout, samname])]),
                              shell=True)
    except subprocess.CalledProcessError:
        print("ERROR.Mapping of %s with bowtie2 failed. Stop execution." % fastqname)
        sys.exit(1)
    else:
        print("Mapping of %s with bowtie2 complete." % fastqname)


def bowtie2mappingsinglend(indexname, fastqname, samname):
    try:
        subprocess.check_call(" ".join([bowtie2path, "-p", threads, "-q", "--local", "-k", multimap,
                                        "-x", indexname,
                                        "".join([raw_data_dir, fastqname]),
                                        "-S", "".join([mappingout, samname])]),
                              shell=True)
    except subprocess.CalledProcessError:
        print("ERROR.Mapping of %s with bowtie2 failed. Stop execution." % fastqname)
        sys.exit(1)
    else:
        print("Mapping of %s with bowtie2 complete." % fastqname)


if __name__ == "__main__":
    index = checkindex(GRCm38indexpath, GRCm38fasta, indexname)
    createdir(mappingout)
    os.chdir(raw_data_dir)
    for key, value in raw_fastq2.items():
        bowtie2mappingsinglend("".join([GRCm38indexpath, "/", index]), value, key)
    for key, value in raw_fastq.items():
        bowtie2mappingpairedend("".join([GRCm38indexpath, "/", index]), value, key)
