import sys
import subprocess
import os
import traceback

#
star = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/STAR"  # path of STAR alignment
multiqc = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/multiqc"  # path of multiqc tool
# set variables
thread = 40  # number of processor
maindir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010"  # home of project folder
starindex = "/home/spuccio/AnnotationSTAR/GencodeMusMusculus/STAR/89/"  # path of Star genome index
rawdatadir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/raw_data_GSE116347/"  # raw reads path
mappingdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/star_Pnas_GSE116347"  # BAM path
fastqcdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/fastqc_Pnas_GSE116347"  # Fastqc path
multiqcdir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/multiqc_Pnas_GSE116347"  # multiqc output path
STARparameters = " ".join(
    ["--runThreadN", str(thread), "--sjdbOverhang", "24", "--sjdbGTFtagExonParentGene", "gene_name",
     "--outFilterMultimapNmax", "1", "--outSAMtype", "BAM", "SortedByCoordinate",
     "--outFilterIntronMotifs", "RemoveNoncanonical", "--quantMode", "GeneCounts",
     "--outFilterScoreMinOverLread", "0", "--outFilterMatchNminOverLread", "0",
     "--outFilterMatchNmin", "0"])  # STAR parameters
#
genomefasta = "/home/spuccio/AnnotationSTAR/GencodeMusMusculus/GRCm38.chr.genome.fa" # Fasta file path
genomegtf = "/home/spuccio/AnnotationSTAR/GencodeMusMusculus/gencode.vM19.annotation.gtf" # annotation file path
# list
fastqreads = [["SRR7443376_1.fastq", "SRR7443376_2.fastq"],
              ["SRR7443375_1.fastq", "SRR7443375_2.fastq"],
              ["SRR7443374_1.fastq", "SRR7443374_2.fastq"],
              ["SRR7443370_1.fastq", "SRR7443370_2.fastq"],
              ["SRR7443369_1.fastq", "SRR7443369_2.fastq"],
              ["SRR7443368_1.fastq", "SRR7443368_2.fastq"],
              ["SRR7443371_1.fastq", "SRR7443371_2.fastq"],
              ["SRR7443372_1.fastq", "SRR7443372_2.fastq"],
              ["SRR7443373_1.fastq", "SRR7443373_2.fastq"]]  # paired end modality


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
    subprocess.check_call(" ".join([star, "--runThreadN", str(thread), "--runMode","genomeGenerate", "--genomeDir",
                                    starindex, "--genomeFastaFiles",
                                    genomefasta, "--sjdbGTFfile", genomegtf, "--sjdbOverhang", "24"]), shell=True)
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
        subprocess.check_call(" ".join([star, "--genomeDir", starindex, "--readFilesIn", rawdatadir + fastqreads[i][0],
                                        rawdatadir + fastqreads[i][1],
                                        "--outFileNamePrefix", fastqreads[i][0].split("_")[0], STARparameters,
                                        "--sjdbGTFfile", genomegtf]), shell=True)
    except (subprocess.CalledProcessError, traceback), e:
        print "STAR mapping of %s fail. Exit" % fastqreads[i][0].split("_")[0]
        sys.exit(1)
    print "STAR mapping of %s complete." % fastqreads[i][0].split("_")[0]
# move to main folder
os.chdir(maindir)
# execute multiqc
try:
    subprocess.check_call(" ".join([multiqc, fastqcdir, mappingdir, "-o", multiqcdir]), shell=True)
except (subprocess.CalledProcessError, traceback), e:
    print "MultiQC process fail.Exit"
    sys.exit(1)
else:
    print "MultiQC process complete."
