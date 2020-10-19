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
#SBATCH --array 1-12
#SBATCH --job-name RADADTRIM
#SBATCH --partition long
#SBATCH --output %x_%A_%a.out
#SBATCH --error %x_%A_%a.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 10G
#SBATCH --mail-user william.robbins@unil.ch
#SBATCH --mail-type END
#SBATCH --export NONE
#SBATCH --time 2-00:00:00

##SBATCH --reservation=HPC-course
##SBATCH --mem 12G
##SBATCH --time 0-12:00:00

mkdir 02_RAD_TRIMMED_READS

# Load modules needed for your job.
# ********************************
# Reminder: you can search for software with the "vit_soft" command.
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


#######################################
##### Filtering reads for quality #####
#######################################


# Trim adapters from R1 and R2 files 
	# custom adapters consist of 5' and 3' tags
	# mm5 and mm3 are mismatches allowed
	
	# for R1 files
		# tag5 = CGACGCTCTTCCGATCT
		# tag3 = TAAGATCGGAAGAGCGAGA

	# for R2 files
		# tag5 = GTGTGCTCTTCCGATCT
		# tag3 = AGATCGGAAGAGCGTCG
		
#F1=$(ls 01_RAD_FILTERED_READS/*_R1_filtered.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)

F2=$(ls 01_RAD_FILTERED_READS/*_R2_filtered.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)

#RAD1=$(echo $F1 | cut -d'/' -f2 | cut -d'_' -f1)

RAD2=$(echo $F2 | cut -d'/' -f2 | cut -d'_' -f1)

#Read1=$(echo $F1 | cut -d'/' -f2 | cut -d'_' -f2)

Read2=$(echo $F2 | cut -d'/' -f2 | cut -d'_' -f2)

#Tag5_1=$(echo 'CGACGCTCTTCCGATCT')

Tag5_2=$(echo 'GTGTGCTCTTCCGATCT')

#Tag3_1=$(echo 'TAAGATCGGAAGAGCGAGA')

Tag3_2=$(echo 'AGATCGGAAGAGCGTCG')

#/scratch/wally/FAC/FBM/DEE/isanders/popgen_to_var/wrobbins/crobbins/RADseq/scripts/tagcleaner.pl -fastq $F1 -tag5 $Tag5_1 -tag3 $Tag3_1 -matrix exact -trim_within 60 -mm5 2 -mm3 2 -log -verbose -out_format 3 -out 02_RAD_TRIMMED_READS/$RAD1'_'$Read1'_trim' 2> 02_RAD_TRIMMED_READS/'log_'$RAD1'_'$Read1'_trim'

/scratch/wally/FAC/FBM/DEE/isanders/popgen_to_var/wrobbins/crobbins/RADseq/scripts/tagcleaner.pl -fastq $F2 -tag5 $Tag5_2 -tag3 $Tag3_2 -matrix exact -trim_within 60 -mm5 2 -mm3 2 -log -verbose -out_format 3 -out 02_RAD_TRIMMED_READS/$RAD2'_'$Read2'_trim' 2> 02_RAD_TRIMMED_READS/'log_'$RAD2'_'$Read2'_trim'

#gzip $F1
gzip $F2

exit 0
