#!/bin/bash
#SBATCH --job-name=SIRsimsarray
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=1000
#SBATCH --time=50:00:00
#SBATCH --output=/home/jbender/Documents/SCI499/slurmout/%j.out

for mm in 1.0 0.9 0.7 0.5 0.3 0.1
do 
    julia --project=/home/jbender/Documents/SCI499  slurmsimsjulia.jl $mm
done