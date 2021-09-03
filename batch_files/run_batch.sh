#!/bin/bash
# SLURM resource specifications
# (use an extra ‘#’ in front of SBATCH to comment-out any unused options)
#SBATCH --job-name=k235_day2  # shows up in the output of ‘squeue’
#SBATCH --time=1-00:00:00       # specify the requested wall-time
#SBATCH --partition=astro_short  # specify the partition to run on
#SBATCH --nodes=1              # number of nodes allocated for this job
#SBATCH --ntasks-per-node=20    # number of MPI ranks per node
#SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank
#SBATCH --mail-type=ALL
srun dynamics_035_k222_day2.exe