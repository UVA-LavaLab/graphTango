#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --ntasks-per-node=1
#SBATCH --reservation=fas9nw_68
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --nodelist=cortado01
#SBATCH --partition=main,gpu
#SBATCH --job-name="gtango"
#SBATCH --error="err_slurm.log"
#SBATCH --output="out_slurm.log"

module load gcc

../profile.sh