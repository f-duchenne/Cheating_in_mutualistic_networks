#!/bin/bash

#SBATCH --job-name=solve             # job name
#SBATCH --ntasks=1                        # ask for 1 task
#SBATCH --array=1-500%500
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00                   # time limit in hrs:min:sec
#SBATCH --output=anal_eq_%j.out       # output and error log

SCRIPT="/home/duchenne/cheating/solve_model.r"
singularity exec /home/duchenne/ubuntu2.sif R --vanilla -f $SCRIPT --args $SLURM_ARRAY_TASK_ID
