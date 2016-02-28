#!/bin/bash -i
#$ -S /bin/bash
#
# MPI-PKS script for job submission script with 'qsub'.
# Syntax is Bash with special qsub-instructions that begin with '#$'.
# For more detailed documentation, see
#       http://www/closed/getting_started/queuing_system.html
#
# (last change of this file: $Id: f20009313e90b5cdb7168260caf8bf57d96b4f26 $)

# --- Mandatory qsub arguments
# Hardware requirements.
#$ -l h_rss=2000M,h_fsize=10000M,h_cpu=121:00:00,hw=x86_64

# --- Optional qsub arguments
# Change working directory - your job will be run from the directory
# that you call qsub in.  So stdout and stderr will end up there.
#$ -cwd

# --- Job Execution
# For faster disk access copy files to /scratch first.

scratch=/scratch/$USER/$$
mkdir -p $scratch
cd $scratch
cp -r $HOME/Work/Run/code .
cp $HOME/Work/Run/makefile .
cp $HOME/Work/Run/input.in .

# Compile the program 
mkdir -p build
mkdir -p data
make

# change input parameters based on taskID
sed -i.bak "s/\(taskID *= *\)[0-9][0-9]*/\1$SGE_TASK_ID/" input.in

# Execution - running the actual program.
# [Remember: Don't read or write to /home from here.]
echo "Running on $(hostname)"
echo "We are in $(pwd)"
# make run
timeout 120h ./build/a.out $SGE_TASK_ID

# Finish - Copy data back to your home directory, clean up.
result=/data/biophys/wenwen/$JOB_NAME
mkdir -p $result
find . -name r*.dat | xargs -I {} cp {} $result      
echo "========$(date)========" >> $SGE_STDOUT_PATH
echo "====================END=====================" >> $SGE_STDOUT_PATH
cat $SGE_STDOUT_PATH >> $result/${JOB_NAME}.log
cd
rm -rf $scratch
