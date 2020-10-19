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
#SBATCH --job-name RADQS
#SBATCH --partition long
#SBATCH --output %x_%A_%a.out
#SBATCH --error %x_%A_%a.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 15G
#SBATCH --mail-user william.robbins@unil.ch
#SBATCH --mail-type END
#SBATCH --export NONE
#SBATCH --time 2-00:00:00

##SBATCH --reservation=HPC-course
##SBATCH --mem 12G
##SBATCH --time 0-12:00:00

mkdir 03_RAD_QUALSCORE_READS
mkdir 03_RAD_QUALSCORE_READS/bad
mkdir 03_RAD_QUALSCORE_READS/logs

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

# Filter for read quality scores
	# Filter parameters are set for:
		# trim all Ns from left and right
		# trim bases with Phred less than 26 from right and left ###(25)
		# calculate mean Phred of 5 bases and trim by one if less than 30 ###(25)
		# remove all sequences with Ns
		# keep sequences 30 bp and over ###(50)

F1=$(ls 02_RAD_TRIMMED_READS/RAD*_R1_trim.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)

F2=$(ls 02_RAD_TRIMMED_READS/RAD*_R2_trim.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)

RAD=$(echo $F1 | cut -d'/' -f2 | cut -d'_' -f1)

Read=$(echo $F1 | cut -d'/' -f2 | cut -d'_' -f2)

trim_ns_l=$(echo '1')

trim_ns_r=$(echo '1')

trim_qrule=$(echo 'lt')

trim_ql=$(echo '25')

trim_qr=$(echo '25')

trim_qt=$(echo 'mean')

trim_qw=$(echo '5')

trim_qs=$(echo '1')

min_q_mean=$(echo '25')

min_l=$(echo '50')

max_ns=$(echo '0')


/scratch/wally/FAC/FBM/DEE/isanders/popgen_to_var/wrobbins/crobbins/RADseq/scripts/prinseq-lite-0.20.4.pl -fastq $F1 -fastq2 $F2 -trim_ns_left $trim_ns_l -trim_ns_right $trim_ns_r -trim_qual_rule $trim_qrule -trim_qual_left $trim_ql -trim_qual_right $trim_qr -trim_qual_type $trim_qt -trim_qual_window $trim_qw -trim_qual_step $trim_qs -min_qual_mean $min_q_mean -min_len $min_l -ns_max_n $max_ns -log 03_RAD_QUALSCORE_READS/logs/'log_qscore_'$RAD -out_format 3 -out_good 03_RAD_QUALSCORE_READS/$RAD'_good' -out_bad 03_RAD_QUALSCORE_READS/bad/$RAD'_bad'

gzip $F1
gzip $F2
