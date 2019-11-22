import os
import subprocess
import sys
import pysam

#
ProjectDir = '/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4'
mappingfolder = '/mnt/datadisk2/spuccio/SP010_RnaSeq_Irf4/STAR_output/'
infer_experiment = '/home/spuccio/miniconda3/envs/rnaseq_env/bin/infer_experiment.py '
bam_stat = '/home/spuccio/miniconda3/envs/rnaseq_env/bin/bam_stat.py '
read_duplication = '/home/spuccio/miniconda3/envs/rnaseq_env/bin/read_duplication.py '
read_distribution = '/home/spuccio/miniconda3/envs/rnaseq_env/bin/read_distribution.py '
geneBody_coverage = '/home/spuccio/miniconda3/envs/rnaseq_env/bin/geneBody_coverage.py '
junction_annotation = '/home/spuccio/miniconda3/envs/rnaseq_env/bin/junction_annotation.py '
bed12 = '/home/spuccio/AnnotationSTAR/gencode.v27.annotation.bed'
houskeeping = '/home/spuccio/AnnotationSTAR/houskeeping_gencode.v27.bed'
rseqcoutputfolder = ProjectDir + '/RSeQC/'
stdeW = open(ProjectDir + '/stderror_mappingqc.txt', 'w')
stdeA = open(ProjectDir + '/stderror_mappingqc.txt', 'a')


#

def directory_creation(path_dir):
    """
    Make folder project
    :param path_dir: string of path that should be checked or created
    :return:
    """
    if not os.path.isdir(path_dir):
        os.makedirs(path_dir)
        print path_dir + ' created.'
    else:
        print path_dir + ' already exist.'


directory_creation(rseqcoutputfolder)

with open(ProjectDir + '/bam.txt', 'w') as bamlist:
    try:
        os.chdir(mappingfolder)
        subprocess.check_call('ls *.bam', stdout=bamlist, shell=True)
    except subprocess.CalledProcessError:
        print 'Error during list of BAM file. Exit'
        sys.exit(1)
    else:
        print 'List of trimmed fastq created.'

lines = [line.rstrip('\n') for line in open(bamlist.name)]
# infer_experiment
for i in range(len(lines)):
    prefix = lines[i].split('.')[0]
    bamfile = mappingfolder + lines[i]
    try:
        with open(rseqcoutputfolder + prefix + '_infer_experiment.txt', 'w') as inferexp_outfile:
            subprocess.check_call(infer_experiment + ' '.join(['-i', bamfile, '-r', bed12]),
                                  stdout=inferexp_outfile, stderr=stdeA, shell=True)
    except subprocess.CalledProcessError:
        print 'Error infer experiment RseQC. Exit'
        sys.exit(1)
    else:
        print 'Infer experiment file for ' + bamfile.split('/')[-1] + ' created.'
# bam_stat
for i in range(len(lines)):
    prefix = lines[i].split('.')[0]
    bamfile = mappingfolder + lines[i]
    try:
        with open(rseqcoutputfolder + prefix + '_bam_stat.txt', 'w') as bam_stat_outfile:
            subprocess.check_call(bam_stat + ' '.join(['-i', bamfile]),
                                  stdout=bam_stat_outfile, stderr=stdeA, shell=True)
    except subprocess.CalledProcessError:
        print 'Error bam stat RseQC. Exit'
        sys.exit(1)
    else:
        print 'Bam stat file experiment for ' + bamfile.split('/')[-1] + ' created.'
# read_duplication
for i in range(len(lines)):
    prefix = lines[i].split('.')[0]
    bamfile = mappingfolder + lines[i]
    try:
        subprocess.check_call(read_duplication + ' '.join(['-i', bamfile, '-o', rseqcoutputfolder + prefix])
                              , stderr=stdeA, shell=True)
    except subprocess.CalledProcessError:
        print 'Error read duplication RseQC. Exit'
        sys.exit(1)
    else:
        print 'Read duplication file experiment for ' + bamfile.split('/')[-1] + ' created.'
# read_distribution
for i in range(len(lines)):
    prefix = lines[i].split('.')[0]
    bamfile = mappingfolder + lines[i]
    try:
        with open(rseqcoutputfolder + prefix + '_read_distribution.txt', 'w') as read_distribution_outfile:
            subprocess.check_call(read_distribution + ' '.join(['-i', bamfile, '-r', bed12]),
                                  stdout=read_distribution_outfile, stderr=stdeA, shell=True)
    except subprocess.CalledProcessError:
        print 'Error read distribution RseQC. Exit'
        sys.exit(1)
    else:
        print 'Read distribution file experiment for ' + bamfile.split('/')[-1] + ' created.'
# # geneBody_coverage
# for i in range(len(lines)):
#     prefix = lines[i].split('.')[0]
#     bamfile = mappingfolder + lines[i]
#     try:
#         with open(rseqcoutputfolder + prefix + '_geneBody_coverage.txt', 'w') as geneBody_coverage_outfile:
#             subprocess.check_call(geneBody_coverage + ' '.join(['-r',houskeeping,'-i', mappingfolder, '-o', geneBody_coverage_outfile.name]),
#                                   stderr=stdeA, shell=True)
#     except subprocess.CalledProcessError:
#         print 'Error geneBody coverage RseQC. Exit'
#         sys.exit(1)
#     else:
#         print 'GeneBody coverage file experiment for ' + bamfile.split('/')[-1] + ' created.'
# junction_annotation
for i in range(len(lines)):
    prefix = lines[i].split('.')[0]
    bamfile = mappingfolder + lines[i]
    try:
        with open(rseqcoutputfolder + prefix + '_junction_annotation.txt', 'w') as junction_annotation_outfile:
            subprocess.check_call(
                geneBody_coverage + ' '.join(['-i', bamfile, '-r', bed12, '-o', junction_annotation_outfile.name]),
                stderr=stdeA, shell=True)
    except subprocess.CalledProcessError:
        print 'Error junction annotation RseQC. Exit'
        sys.exit(1)
    else:
        print 'Junction annotation file experiment for ' + bamfile.split('/')[-1] + ' created.'
