import os
import shutil

command1 = "{samtools} view -b -F 1548 -q {qualitythr} {inputbam} | {bam2bed} -i -  | gzip -c > {tagalignfile} "
# command2 = r'awk \'BEGIN\{FS=\"\t\";OFS=\"\t\"}{$4=\"N\"; print $0}\' {BEDFILE} | gzip -c > {tagalignfile}'
command3 = "{macs2} callpeak -t {tagalign1} {tagalign2} {tagalign3} -f BED -g mm -n {outputname} -p {pvalthr} " \
           "--outdir {outdir} "
command4 = "{macs2} callpeak -t {tagalign1} {tagalign2} {tagalign3} -f BED -g mm -n {outputname} -p {pvalthr} " \
           "--nomodel --extsize 147 --outdir {outdir} "


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


def createtagalign(samtoolspath, qualthr, inputbam, bamtobed, outputtagfile):
    """

    :param samtoolspath:
    :param qualthr:
    :param inputbam:
    :param bamtobed:
    :param outputtagfile:
    :return:
    """
    os.system(command1.format(samtools=samtoolspath, qualitythr=qualthr, inputbam=inputbam, bam2bed=bamtobed,
                              tagalignfile=outputtagfile))


def peakscalling(macs2path, tagalign1, tagalign2, tagalign3, outputsuffixname, pval, outputdir):
    """

    :param macs2path:
    :param tagalign1:
    :param tagalign2:
    :param tagalign3:
    :param outputsuffixname:
    :param pval:
    :param outputdir:
    :return:
    """
    os.system(command3.format(macs2=macs2path, tagalign1=tagalign1, tagalign2=tagalign2, tagalign3=tagalign3,
                              outputname=outputsuffixname, pvalthr=pval,
                              outdir=outputdir))


def peakscallingnomodel(macs2path, tagalign1, tagalign2, tagalign3, outputsuffixname, pval, outputdir):
    """

    :param macs2path:
    :param tagalign1:
    :param tagalign2:
    :param tagalign3:
    :param outputsuffixname:
    :param pval:
    :param outputdir:
    :return:
    """
    os.system(command4.format(macs2=macs2path, tagalign1=tagalign1, tagalign2=tagalign2, tagalign3=tagalign3,
                              outputname=outputsuffixname, pvalthr=pval,
                              outputdir=outputdir))


if __name__ == "__main__":
    MACS2_P_Val = "1e-3"
    BAMfile = {"Treg_Irf4_r1": "Treg_Irf4_r1_rmdup.bam",
               "Treg_Irf4_r2": "Treg_Irf4_r2_rmdup.bam",
               "Treg_Irf4_r3": "Treg_Irf4_r3_rmdup.bam"}
    inputpath = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/bowtie2_mapping_GSE98264/"
    outputpath = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/macs2_peaks_GSE98264/"
    qualth = "30"
    createdir(outputpath)
    for key, value in BAMfile.items():
        createtagalign(shutil.which('samtools'), qualth, "".join([inputpath, value]), shutil.which('bamToBed'),
        "".join([outputpath, key, ".bed.gz"]))
    listA = []
    for key in BAMfile.items():
        listA.append(key[0])
    peakscalling(shutil.which('macs2'), "".join([outputpath, listA[0], ".bed.gz"]),
                 "".join([outputpath, listA[1], ".bed.gz"]),
                 "".join([outputpath, listA[2], ".bed.gz"]), "Treg_Irf4", MACS2_P_Val, outputpath)
    peakscallingnomodel(shutil.which('macs2'), "".join([outputpath, listA[0], ".bed.gz"]),
                        "".join([outputpath, listA[1], ".bed.gz"]),
                        "".join([outputpath, listA[2], ".bed.gz"]), "Treg_Irf4_nomodel", MACS2_P_Val, outputpath)
