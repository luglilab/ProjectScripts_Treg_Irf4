import os

# strand shifts at which cross-correlation
command1 = "{Rscript} {phantompeak} -rf -s=-100:5:600 -c={inputfile} -p={thread} -savp -out={output}"
command2 = "{bamcoverage} -b {inputbam} -o {outputbigwig} -p {thread} -of bigwig -bs 10 --effectiveGenomeSize {genomesize} --normalizeUsing RPKM -e 200 "

def phantompeak(RscriptPath,spppath,bamfile,thread,pdffile):
    """
    Phantom peak - strand cross-correlation peak / predominant fragment length
    :param RscriptPath:
    :param spppath:
    :param bamfile:
    :param pdffile:
    :return:
    """
    os.system(command1.format(Rscript=RscriptPath, phantompeak=spppath, inputfile=bamfile,thread=thread, output=pdffile))

def bam2bigwig(bamcoveragePath,bamfile,bigwipout,thread,gensize):
    """
    Conversion BAM to Big WIG normalized
    :param bamcoveragePath:
    :param bamfile:
    :param bigwipout:
    :param gensize:
    :return:
    """
    os.system(command2.format(bamcoverage=bamcoveragePath, inputbam=bamfile,outputbigwig=bigwipout,thread=thread,genomesize=gensize))

def createdir(dirpath):
    """
    Make directory function
    :param dirpath:
    :return:
    """
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
        print(" ".join(["Directory", dirpath.split("/")[-1], "Created"]))
    else:
        print(" ".join(["Directory", dirpath.split("/")[-1], "already exists"]))

if __name__ == "__main__":
    Rscript = "/home/spuccio/miniconda3/envs/chipseq_env/bin/Rscript"
    spp = "/home/spuccio/miniconda3/envs/chipseq_env/bin/run_spp.R"
    BAMfile = {"Treg_Irf4_r1" : "Treg_Irf4_r1_rmdup.bam",
               "Treg_Irf4_r2" : "Treg_Irf4_r2_rmdup.bam",
               "Treg_Irf4_r3": "Treg_Irf4_r3_rmdup.bam"}
    inputpath = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE98264/"
    outputpath = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/phantompeak_bigwig_GSE98264/"
    bamcoverage= "/home/spuccio/miniconda3/bin/bamCoverage"
    genomesize="2308125349"
    thread="30"
    createdir(outputpath)
    for key,value in BAMfile.items():
        phantompeak(RscriptPath=Rscript, spppath=spp, bamfile="".join([inputpath,value]),thread=thread,pdffile="".join([outputpath,key,".pdf"]))
        bam2bigwig(bamcoveragePath=bamcoverage,bamfile="".join([inputpath,value]),bigwipout="".join([outputpath,key,".bw"]),thread=thread,gensize=genomesize)

