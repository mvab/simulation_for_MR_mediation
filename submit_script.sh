#!/bin/bash

# input from command line
mediator_row="$1"  # chromosome to process; add ' -F "mediator_row #" ' at the end to submission script, e.g. -F "2"
cores=16
iters=1000
data_path="/newhome/ny19205/simulation_for_MR_mediation/data/"

# request resources:
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00

## not using: #PBS -q himem


# Define working directory
export WORK_DIR=$HOME/simulation_for_MR_mediation

# load latest R
module add languages/R-3.5.1-ATLAS-gcc-6.1

# on compute node, change directory to 'submission directory':
cd $WORK_DIR


echo JOB ID: $PBS_JOBID
echo Working Directory: `pwd`
echo Start Time: `date`

# run it, timing it:
time Rscript process.R  ${cores}  ${iters} ${data_path}  $mediator_row 

echo End Time: `date`
