import pandas as pd
import os
import subprocess
import datetime

command1 = "samtools view -@ {thread} -S {SAM} | cut -f 1 | sort | uniq | wc -l"
command2 = "samtools view -@ {thread} -S -F 0x04 {SAM} | cut -f 1 | sort | uniq | wc -l"
command3 = "samtools view -@ {thread} -Sh {SAM} | awk \'(substr($1, 1, 1)==\"@\") ||  (( $3 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|M|Y)$/ ) &&  ( ( $7 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|M|Y)$/ ) ||  ($7==\"=\") ||  ($7==\"*\") ))\' - |  samtools view -bhS - > {filteredbam}"
command4 = "samtools view -@ {thread} -h {BAM} | sed \'/chrM/d;/random/d;/chrUn/d\' - | samtools view -Shb - > {NOMITO}"
command5 = "samtools view -@ {thread} -hb -F 1804 {NOMITO} > {UNIQMAP}"
command6 = "samtools view -@ {thread} -hb -q {MAPQ} {UNIQMAP} | samtools sort -o {BAM} - "
command7 = "samtools index -@ {thread} {BAM}"
command8 = "MarkDuplicates INPUT={BAM} OUTPUT={RMDUPBAM} ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE={PICARD} 2>/dev/null"
#
thread = "30"
mapqth= "30"
#
SAMpath = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/"

fileSAM = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE99167/test.sam"
filteredbam = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE99167/test_filtered.bam"
nomitobam = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE99167/test_filtered_nomito.bam"
uniqmapbam = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE99167/test_filtered_nomito_uniq.bam"
outbam = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE99167/test_filtered_nomito_uniq.bam"
rmdupoutbam = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE99167/test_filtered_nomito_uniq_rmdup.bam"
picardmetrix = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE99167/test_filtered_nomito_uniq_rmdup.txt"
logfile = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE99167/test_filtered_nomito_uniq_rmdup.log"


def postalign(thread,bowtiesam,chrfiltered,mitfiltered,uniqmap,mapqth,outputbam,rmdupbam,picardmetrix,logfile):
    """

    :param thread: Number of threads for samtools
    :param bowtiesam: input SAM file generated with Bowtie2
    :param chrfiltered: Path and name of BAM file filtered for contigs chr
    :param mitfiltered: Path and name of BAM file filtered for mitocondrial reads
    :param uniqmap: Path and name of BAM file filtered for multimapping reads
    :param mapqth: Path and name of BAM file filtered for mapping quality
    :param outputbam: Path and name of BAM file filtered with duplicate reads
    :param rmdupbam: Path and name of BAM file filtered for duplicate reads
    :param picardmetrix: Picard metrix file
    :param logfile: Log file
    :return:
    """
    with open(logfile, 'w') as log:
        log.write(datetime.datetime.now().strftime("%y-%m-%d-%H-%M") + "\n")
        # total no of reads - mapped / unmapped
        log.write("Post Alignment log of {INPUTSAM}.\n".format(INPUTSAM=bowtiesam.split("/")[-1]))
        log.write("TotalRawReads: \t {TotalRawReads}".format(
            TotalRawReads=os.popen(command1.format(thread=thread, SAM=bowtiesam)).read()))
        # number of mappable reads (excluding the unmapped reads)
        log.write("NumMappableRead: \t {NumMappableRead}".format(NumMappableRead=os.popen(command2.format(thread=thread, SAM=bowtiesam)).read()))
        # delete the random stuffs (say chr1_... etc.. also) the modified stuff includes chr1 ... chr22 chrX chrY and chrM (adapted for Human and Mouse)
        os.system(command3.format(thread=thread, SAM=bowtiesam, filteredbam=chrfiltered))
        # count the number of reads remaining
        log.write("NumReadNotInChr1-19XYM: \t {READELETECHR}".format(READELETECHR=os.popen(command1.format(thread=thread, SAM=chrfiltered)).read()))
        # delete the mitochondrial reads
        os.system(command4.format(thread=thread, BAM=chrfiltered, NOMITO=mitfiltered))
        # number of reads after mitochondrial read delete
        log.write("NumReadMithocondrial: \t {READELETEMIT}".format(READELETEMIT=os.popen(command1.format(thread=thread, SAM=mitfiltered)).read()))
        # the flag 1804 = read unmapped, mate unmapped, not primary alignment, read quality low, PCR / optical duplicate
        os.system(command5.format(thread=thread, NOMITO=mitfiltered, UNIQMAP=uniqmap))
        # number of uniquely mapped reads
        log.write("UniqMappedRead: \t {UNIQMAPREADS}".format(UNIQMAPREADS=os.popen(command1.format(thread=thread, SAM=uniqmap)).read()))
        # perform quality based thresholding and sorting operation
        os.system(command6.format(thread=thread, MAPQ=mapqth, UNIQMAP=uniqmap, BAM=outputbam))
        # count the number of reads after quality thresholding
        log.write("ReadQualThr: \t {READSQUALTHR}".format(READSQUALTHR=os.popen(command1.format(thread=thread, SAM=outputbam)).read()))
        # index the sorted file
        os.system(command7.format(thread=thread, BAM=outputbam))
        # remove any PCR duplicates using Picard tool
        os.system(command8.format(BAM=outputbam, RMDUPBAM=rmdupbam, PICARD=picardmetrix))
        # count the number of reads after removing duplicate reads
        log.write("ReadAfterRmDup: \t {READSRMDUP}".format(READSRMDUP=os.popen(command1.format(thread=thread, SAM=rmdupbam)).read()))


if __name__ == "__main__":
    samfile = {"Treg_Irf4_r1": "Treg_Irf4_r1.sam",
               "Treg_Irf4_r2": "Treg_Irf4_r2.sam",
               "Treg_Irf4_r3": "Treg_Irf4_r3.sam"}

    postalign(thread = thread,
              bowtiesam = fileSAM,
              chrfiltered = filteredbam,
              mitfiltered = nomitobam,
              uniqmap = uniqmapbam,
              mapqth = mapqth,
              outputbam = outbam,
              rmdupbam = rmdupoutbam,
              picardmetrix = picardmetrix,
              logfile = logfile)




