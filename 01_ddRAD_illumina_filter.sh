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
#SBATCH --array 1-24
#SBATCH --job-name RADQFILTER
#SBATCH --partition long
#SBATCH --output %x_%A_%a.out
#SBATCH --error %x_%A_%a.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 8G
#SBATCH --mail-user william.robbins@unil.ch
#SBATCH --mail-type END
#SBATCH --export NONE
#SBATCH --time 0-12:00:00

##SBATCH --reservation=HPC-course
##SBATCH --mem 12G
##SBATCH --time 0-12:00:00

mkdir 01_RAD_FILTERED_READS

module add Bioinformatics/Software/vital-it

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


##########################################
###### Download raw data from server #####
##########################################

# Make file with all links

# vim RADseq_files.txt

	# press "i" for insert
	# copy and paste all file links
#	http://uhts-lgtf.vital-it.ch/symlink/RAD2_L7_R1_001_hYRAEsxgNCtr.fastq.gz
#	all
#	file
#	links
#	http://uhts-lgtf.vital-it.ch/symlink/RAD2_L7_R2_029_Ad3bjczSAPpu.fastq.gz
	# press "esc" to escape editing
	# type ":wq" to write the file and exit

# Download, unzip and touch all files:

for file in $(cat RADseq_files.txt); do wget $file; unzip=$(echo $file | cut -d'/' -f5); gunzip $unzip;done

# Check files are present and synchronized

#for file1 in $(ls *_R1_*.fastq); do RAD=$(echo $file1 | cut -d '_' -f1); pair=$(echo $file1 | cut -d'_' -f4); file2=$(ls $RAD*_R2*$pair*.fastq); echo $file1" "$file2 ;done


#######################################
##### Filtering reads for quality #####
#######################################

# CASAVA quality filter, removing :Y:

	# example:
	# @HISEQ:331:CB8W2ANXX:5:1101:1246:1957 2:Y:0:  <- CASAVA quality code :Y:
	# GCATCAACAGGCCACAACCAACCAGAACGTGAAAAAGCGTCCTGCGTGTAGCGAACTGCGATGGGCATACTGTAACCATAAGGCCACGTATTTTGCAAGCT
	# +
	# BBBBBFBFFFFBFBFFFFBF/F7<FBFB/FFB/FFB7FFFFFFFFF/FFFFFB/FFBFFFFFBFF/FFFFFBFFFFFFFFFFBBF/F/<BFFFF7B7FBF#


ZFile=$(ls *.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

gunzip $ZFile

File=$(echo $ZFile | cut -d'.' -f1,2)

RAD=$(echo $File | cut -d '_' -f1)

Read=$(echo $File | cut -d'_' -f3)

cat $File | grep -A 3 '^@.*[^:]*:N:[^:]*:' | grep -v "^--$" > 01_RAD_FILTERED_READS/$RAD'_'$Read'_filtered.fastq'

gzip $File

exit 0