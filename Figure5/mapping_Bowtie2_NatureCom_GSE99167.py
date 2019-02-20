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
GRCm38indexpath = "/home/spuccio/AnnotationBowtie2/Mus_musculus/GRCm38.p6"
# bowtie2-build path
bowtiebuild2path = "/home/spuccio/miniconda3/envs/chipseq_env/bin/bowtie2-build"
# bowtie2 path
bowtie2path = "/home/spuccio/miniconda3/envs/chipseq_env/bin/bowtie2"
# path folder
projectdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/"
# path out file sam
mappingout = "".join([projectdir, "bowtie2_mapping_GSE99167/"])
# input file path
raw_data_dir = "".join([projectdir, "raw_data_Nature_C_GSE99167/"])
# Path
GRCm38fasta = "".join([GRCm38indexpath, "/GRCm38.primary_assembly.fa"])
raw_fastq = {"Th9_BATF_r1": ["SRR5582853_1.fastq", "SRR5582853_2.fastq"],
             "Th9_BATF_r2": ["SRR5582854_1.fastq", "SRR5582854_2.fastq"],
             "Th9_BATF_r3": ["SRR5582855_1.fastq", "SRR5582855_2.fastq"],
             "Th9_BATF_r4": ["SRR5582856_1.fastq", "SRR5582856_2.fastq"],
             "Th9_INPUT_r1": ["SRR5582865_1.fastq", "SRR5582865_2.fastq"],
             "Th9_INPUT_r2": ["SRR5582866_1.fastq", "SRR5582866_2.fastq"],
             "Th9_INPUT_r3": ["SRR5582867_1.fastq", "SRR5582867_2.fastq"],
             "Th9_INPUT_r4": ["SRR5582868_1.fastq", "SRR5582868_2.fastq"]
             }


def createdir(dirpath):
    """
    Make dir function and check if directory is already exists
    :param dirpath: string with path and directory name
    :return:
    """
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
        print(" ".join(["Directory", dirpath.split("/")[-1], "Created"]))
    else:
        print(" ".join(["Directory", dirpath.split("/")[-1], "already exists"]))


def checkindex(indexpath, fastafile, indexname):
    """
    Check if the genome index already exist and in case create a new one
    :param indexpath: Output folder
    :param fastafile: genome fasta file
    :param indexname: index name
    :return:
    """
    for i in range(1, 5):
        if os.path.isfile("".join([indexpath, "/", indexname, ".", str(i), ".bt2"])) == True:
            print("Genome index %s.%d.bt2 already exists." % (indexname, i))
        else:
            try:
                os.chdir(GRCm38indexpath)
                subprocess.check_call(" ".join([bowtiebuild2path, "--threads", threads, fastafile, indexname]),
                                      shell=True)
            except subprocess.CalledProcessError:
                print("ERROR.Fastqc analysis failed. Stop execution.")
                sys.exit(1)
            else:
                print("Index %d OK" % i)
    return indexname


def bowtie2mappingpairedend(indexname, fastqname, pathoutput, samname):
    """
    Bowtie2 mapping Paired-end mode
    :param indexname: path with index name
    :param fastqname: fastq file
    :param samname: output SAM name
    :param pathoutput: path folder output
    :return:
    """
    try:
        subprocess.check_call(" ".join([bowtie2path, "-p", threads, "-q", "--local", "-k", multimap,
                                        "-x", indexname, "-X", MAX_FRAG_LEN,
                                        "-1", "".join([raw_data_dir, fastqname[0]]),
                                        "-2", "".join([raw_data_dir, fastqname[1]]),
                                        "-S", "".join([pathoutput, samname, ".sam"])]),
                              shell=True)
    except subprocess.CalledProcessError:
        print("ERROR.Mapping of %s with bowtie2 failed. Stop execution." % fastqname)
        sys.exit(1)
    else:
        print("Mapping of %s with bowtie2 complete." % fastqname)


if __name__ == "__main__":
    index = checkindex(GRCm38indexpath, GRCm38fasta, "GRCm38")
    createdir(mappingout)
    os.chdir(raw_data_dir)
    for key, value in raw_fastq.items():
        bowtie2mappingpairedend("".join([GRCm38indexpath, "/", index]), value, mappingout, key)
