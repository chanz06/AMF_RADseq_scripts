#!/bin/bash

############################################
############################################
##                                        ##
##      Analysis of RADseq data from      ##
##    Rhizophagus irregularis isolates    ##
##                                        ##
############################################
############################################

############################################
############################################
##                                        ##
##          Version October 2018          ##
##                                        ##
##  Author: Fred Masclaux  (August 2017)  ##
##  Edited: Chanz Robbins  (October 2018) ##
##                                        ##
############################################
############################################

############################################
############################################
##                                        ##
##              SLURM OPTIONS             ##
##                                        ##
############################################
############################################

# SLURM option summary:
# ********************
# Options must be prefixed with '#SBATCH'. Lines starting with '##SBATCH' are considered comments.
# --account / -A  : account where the job CPU usage should be charged. If no account is specified, 
#                   the user's default account is charged.
# --job-name / -J : gives a specific name to the job.
# --partition / -p: partition (i.e. queue) where the job should run. Currently 4 partitions are
#                   available: 'normal', 'long', 'ax-normal' and 'ax-long'.
# --output / -o   : write the job’s standard output to the specified file. If this option is 
#                   omitted, SLURM saves the stdout/stderr to a file named "slurm-<jobID>.out".
# --error / -e    : write the job’s standard error to the specified file.
# --nodes / -N    : number of nodes requested for the global allocation. Unless you are running a 
#                   software with MPI support, this value will always be "1".
# --ntasks / -n   : number of tasks requested for the global allocation. In most cases the value 
#                   for this option can be left at "1". If the distribution of tasks across 
#                   nodes is important, use "--ntasks-per-node" instead of "--ntasks".
# --cpus-per-task : number of processors to be allocated to each task. All processors for a task 
#                   are allocated on a same node.
# --mail-user     : email address where notifications should be sent.
# --mail-type     : type of events for which slurm should send an email.
# --mem           : memory limit per node. Use 'M' or 'G' as suffix to specify if the units are 
#                   megabytes or gigabytes, e.g. "12G" for 12 gigabytes. It is also possible to 
#                   use "--mem-per-cpu" instead of "--mem" to link memory request to the number of
#                   processors requested for the job.
# --time          : set limit to the run time of a job (wall time). If the runtime limit is reached,
#                   the job gets killed.
#
# SLURM environment variables:
# ***************************
# $SLURM_JOB_ID       : job ID value.
# $SLURM_JOB_NAME     : job name, i.e. value passed to --job-name.
# $SLURM_CPUS_PER_TASK: number of CPUs allocated to current task.
#
# SLURM array jobs:
# ****************
# To submit an array job, add the "--array" option to the script.
# --array / -a       : submit a job array with indexes “index list”.
#
# The following environment variables can then be used in the array job:
# $SLURM_ARRAY_JOB_ID:  job array master ID value. In array jobs, this must be used instead
#                       of $SLURM_JOB_ID.
# $SLURM_ARRAY_TASK_ID: job array index ID number of the current array replicate.
####################################################################################################


# SLURM options:
# *************
#SBATCH --account isanders_popgen_to_var
#SBATCH --array 1-80
#SBATCH --job-name SNPS_FB
#SBATCH --partition wally
#SBATCH --output %x_%A_%a.out
#SBATCH --error %x_%A_%a.err
#SBATCH --nodes 1 
#SBATCH --mem 10G
#SBATCH --mail-user william.robbins@unil.ch
#SBATCH --mail-type END
#SBATCH --export NONE
#SBATCH --time 9-00:00:00

##SBATCH --reservation=HPC-course
##SBATCH --mem 12G
##SBATCH --time 0-12:00:00

mkdir 06_RAD_SNPS_FREEBAYES

module add Bioinformatics/Software/vital-it
module add UHTS/Analysis/freebayes/1.2.0
module add UHTS/Analysis/samtools/1.8
module add UHTS/Analysis/bamaddrg/2012.05.26

##########################################
##### Before starting you will need: #####
##########################################

#         tagcleaner.pl    -> https://sourceforge.net/projects/tagcleaner/
#             please cite  ->   "Schmieder R, YW Lim, F Rohwer, R Edwards. TagCleaner: Identification and removal of tag sequences 
#                                from genomic and metagenomic datasets. BMC Bioinformatics. (2010) 11:341."
#         prinseq-lite.pl  -> https://sourceforge.net/projects/prinseq/
#             please cite  ->   "Schmieder R & R Edwards. Quality control and preprocessing of metagenomic datasets. Journal 
#                                Bioinformatics. (2011) 27:6."
#         Stacks software  -> http://catchenlab.life.illinois.edu/stacks/
#             please cite  ->   "Catchen J, A Amores, P Hohenlohe, W Cresko, J Postlethwait. Stacks: building and genotyping 
#                                loci de novo from short-read sequences. G3: Genes, Genomes, Genetics, (2011) 1:171-182."
# FastqPairedEndValidator  -> http://www.mcdonaldlab.biology.gatech.edu/bioinformatics/FastqPairedEndValidator.pl
#             please cite  ->   website
#             BWA Aligner  -> http://bio-bwa.sourceforge.net/
#             please cite  ->    "Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler 
#                                 Transform. Bioinformatics, 25:1754-60.

########################################
##### Demultiplexing RAD libraries #####
########################################


# Demultiplex and sort reads by sample
	# STACKS   ->   http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php
	# To properly name files you need a barcode file that looks like this:
		# TCAATATC        sample1
		# ACCGCCTC        sample2
		# ACTTGAAC        sample3
		# etc...

Gen=$(ls 05_RAD_MAPPED/SAM_files/SORTED/*/*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)

Isolate=$(echo $Gen | cut -d'/' -f4)

#########  INDEX GENOME FILES  ##########

samtools faidx  ../ref/$Isolate/$Isolate.fasta


#########  ADD BAM RG HEADERS  ##########

mkdir 06_RAD_SNPS_FREEBAYES/

bam=$(echo $Gen)

sample=$(echo $Gen | cut -d'/' -f4 | cut -d'.' -f1| cut -d'_' -f1,2)

RG=$(echo $Gen | cut -d'/' -f4 | cut -d'.' -f1| cut -d'_' -f3)

bamaddrg -b $bam -s $sample -r $RG > 05_RAD_MAPPED/SAM_files/SORTED/$Isolate/$sample'_'$RG'_RG.bam' 

samtools index -b 05_RAD_MAPPED/SAM_files/SORTED/$Isolate/$sample'_'$RG'_RG.bam'

######## FREEBAYES VARIANT CALLING #########

#### settings
# -J = pooled library
# -K = no scaling of probability based on permutation
# -X = --no-mnps
# -u = --no-complex
# -U = --read-mismatch-limit

####### consider -t = --targets -> bed file from in-silico digest to retrieve best matches

mkdir 06_RAD_SNPS_FREEBAYES/$Isolate

test_alleles=$(echo '4')

keep_alleles=$(echo '3')

freq=$(echo '0.1')

mismatch=$(echo '4')

min_cov=$(echo '20')

F1=$(ls 05_RAD_MAPPED/SAM_files/SORTED/$Isolate/$sample'_'$RG'_RG.bam')

bam=$(echo $F1)

sample=$(echo $F1 | cut -d'/' -f5 | cut -d'.' -f1 | cut -d'_' -f1,2,3)

freebayes -f ../ref/$Isolate/$Isolate*.fasta -p $test_alleles --use-best-n-alleles $keep_alleles -J -X -U $mismatch -F $freq --min-coverage $min_cov -u --vcf 06_RAD_SNPS_FREEBAYES/$Isolate/$sample'.vcf' -b $bam

#gzip $F1
#gzip $F2

exit 0
