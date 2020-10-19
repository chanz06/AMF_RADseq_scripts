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
#SBATCH --array 1-2
#SBATCH --job-name RADinsilico
#SBATCH --partition normal
##SBATCH --output %x_%A_%a.out
##SBATCH --error %x_%A_%a.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 10G
#SBATCH --mail-user william.robbins@unil.ch
#SBATCH --mail-type END
#SBATCH --export NONE
#SBATCH --time 0-24:00:00

##SBATCH --reservation=HPC-course
##SBATCH --mem 12G
##SBATCH --time 0-12:00:00

mkdir 00_RAD_in_silico_digest

module add Bioinformatics/Software/vital-it
module add Emboss/EMBOSS/6.6.0

genomes=$(ls 00_RAD_in_silico_digest/*/*_scaffold*_insidi| cut -d'/' -f2 | sed -n ${SLURM_ARRAY_TASK_ID}p)

genome=$(echo $genomes | cut -d'/' -f2)

for file in $(ls 00_RAD_in_silico_digest/$genome/*scaffold*_insidi); do 
	ID=$(echo $file | cut -d'/' -f3 | cut -d'_' -f1,2)
	sed -n '/..Start/,/^#----*/p' $file | awk '{print $1,$1+$9,$9,$4}'| awk -v scaf="$ID " 'BEGIN{FS=OFS=scaf}{print value OFS $0}'| sed -z "s/$ID/Scaffold/" | awk '$5!=""' > 00_RAD_in_silico_digest/$ID'_sub_insidi'
done

awk 'FNR==1 && NR!=1 { while (/^Scaffold*/) getline;} 1 {print}' 00_RAD_in_silico_digest/$genome/*_scaffold0* > 00_RAD_in_silico_digest/$genome/$genome'_insidi'

exit 0