#!/bin/bash
#SBATCH --job-name=SIRsims

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=1000
#SBATCH --time=50:00:00
#SBATCH --output=/home/jbender/Documents/SCI499/slurmout/%j.out

julia --project=/home/jbender/Documents/SCI499  slurmsims2.jl